"""
DIA Data Cleaning and Feature-Filtering Pipeline
==================================================
Supplementary code for: Toward Effect-Based Monitoring of Wastewater Exposure: DNA Adductomics in Gammarus populations Using Complementary Mass Spectrometric Approaches, Section 2.3.2

Applies the SAME six-step filtering pipeline, in the same order, to both:
  - the untargeted 192-window diagnostic-ion screen ("Win_*" features)
  - library-matched target compounds matched against the reference
    spectral library 

Pre-processing (not one of the six steps): calibration-standard replicates
and the automatically generated precursor-mass transition (product m/z equal
to precursor m/z) are removed before Step 1.

The six numbered steps described in the Methods text:

  (1) Blank filter             -- a replicate counts as detected only if its
                                  raw peak area > 3x the pooled blank mean for
                                  that (feature, fragment) combination. This is
                                  a detection threshold: blank areas are NOT
                                  subtracted from sample areas.
  (2) Normalization            -- peak area expressed as a ratio to the
                                  internal dG standard's peak area
                                  (m/z 152.0567 fragment), same injection
  (3) Reproducibility filter   -- per biological group (n = 4-5):
                                  >= 3 replicates passing Step 1 AND
                                  RSD < 35% of the normalized ratio, computed
                                  over the replicates that passed Step 1
  (4) RT consistency filter    -- (a) within-group RT SD <= 0.3 min, AND
                                  (b) for features retained in multiple groups,
                                  each group's mean RT lies within 0.3 min of
                                  the cross-group median RT
  (5) Peak-shape filter        -- mean peak width <= 0.45 min AND
                                  width CV < 35% across the group's replicates
  (6) Correlation-based        -- any two features whose mean RTs differ by
      artifact removal            <= 0.15 min and whose dG-normalized
                                  intensities correlate across all
                                  experimental replicates with Pearson
                                  r >= 0.90 are both excluded

Two additional, OPTIONAL post-processing steps are included below but are
NOT part of the six-step pipeline reported in the Methods text:
  (a) a library dot-product filter, DISABLED by default (see
      APPLY_LIBRARY_DOTP_FILTER below). Evaluated during method development
      and not applied in the reported pipeline. Mean library dot products
      are still written to the output for transparency.
  (b) a single, manually resolved isomer ambiguity (m/z 167.056,
      9-Methyl-guanine vs. N2-Methyl-guanine) -- a one-off structural
      cleanup, not a general filtering rule

INPUT: a Skyline Document Grid export (CSV) containing blank replicates,
       calibration standards, and experimental samples, with columns:
       Molecule List Name, Molecule Name, Precursor Mz, Precursor Charge,
       Replicate Name, Product Mz, Product Charge, Fragment Ion,
       Retention Time, Area, Height, Start Time, End Time, Mass Error PPM,
       sample type, timepoint, concentration, Bioreplicate, and
       (optional) Library Dot Product (blank/#N/A for non-library rows).
       Blank replicates are recognized by sample type "Blank"/"Blank_exp"
       or, if unannotated, by replicate names 1, 2, 3 / exp_1, exp_2, exp_3
       (optionally prefixed "Blank_"). Update INPUT_FILE below.

OUTPUT: FINAL_UNIFIED_CLEAN.csv (written next to this script's working
        directory) -- every (Molecule List Name, Molecule Name, Fragment Ion,
        Group) combination surviving all six steps.

Applying this pipeline to the raw dataset reported in the paper yields a
final, quality-filtered feature set of 177 unique (analyte, fragment)
combinations (337 total rows across groups): 142 unidentified diagnostic-
ion windows and 35 library-matched target compounds.
"""

import csv
import re
import statistics
from collections import defaultdict

# ============================================================
# CONFIGURATION -- update INPUT_FILE to point at your own Skyline export
# ============================================================
INPUT_FILE = 'Export_skyline_DIA_withBlank_with_Library.csv'
OUT_PATH = 'FINAL_UNIFIED_CLEAN.csv'

FOLD_CHANGE_CUTOFF = 3.0          # Step 1
RSD_CUTOFF = 35.0                 # Step 3
MIN_REPLICATES_FOR_RSD = 3        # Step 3
RT_STD_CUTOFF = 0.3               # Step 4a (min)
CROSS_GROUP_RT_TOLERANCE = 0.3    # Step 4b (min), same tolerance as 4a
MAX_PEAK_WIDTH_MIN = 0.45         # Step 5 (min), ~90th percentile of observed widths
WIDTH_CV_THRESHOLD = 35.0         # Step 5 (%)
CORRELATION_RT_TOLERANCE = 0.15   # Step 6 (min)
CORRELATION_MIN = 0.90            # Step 6 (Pearson r)
LIBRARY_DOTP_THRESHOLD = 0.55     # optional post-processing (a) only; not applied

DG_LIST_NAME = 'dG_GlobalStandard'
DG_MOLECULE_NAME = 'dG_std'
DG_PRODUCT_MZ = 152.0567
MZ_TOLERANCE = 0.001

BLANK_SAMPLE_TYPES = {'Blank', 'Blank_exp'}
VALID_GROUPS = ['T0', 'T1_C1', 'T1_C2', 'T1_C3', 'T2_C0', 'T2_C1', 'T2_C2', 'T2_C3']
RT_EPS = 1e-9                     # guards RT comparisons against floating-point rounding


# ============================================================
# HELPERS
# ============================================================
def get_val(row, col):
    v = row.get(col, '')
    return None if v in ('', '#N/A') else float(v)


def is_blank_replicate(row):
    if row.get('sample type', '') in BLANK_SAMPLE_TYPES:
        return True
    name = row['Replicate Name']
    return bool(re.match(r'^(Blank_)?\d+$', name) or re.match(r'^(Blank_)?exp_\d+$', name))


def assign_group(timepoint, concentration):
    if timepoint == 'T0':
        return 'T0'
    return f"{timepoint}_{concentration}"


def feature_key(row):
    return (row['Molecule List Name'], row['Molecule Name'], row['Fragment Ion'])


def pearson(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    num = sum((a - mx) * (b - my) for a, b in zip(x, y))
    denx = sum((a - mx) ** 2 for a in x)
    deny = sum((b - my) ** 2 for b in y)
    if denx == 0 or deny == 0:
        return 0.0
    return num / (denx * deny) ** 0.5


def n_features(rows):
    return len(set((r['MoleculeListName'], r['MoleculeName'], r['Fragment']) for r in rows))


# ============================================================
# LOAD
# ============================================================
with open(INPUT_FILE, newline='') as f:
    raw_rows = list(csv.DictReader(f))
print(f"Loaded {len(raw_rows)} raw rows")


# ============================================================
# PRE-PROCESSING (not one of the six steps)
# ============================================================
rows = [r for r in raw_rows if r.get('sample type') != 'Calibration']
rows = [r for r in rows if abs(float(r['Product Mz']) - float(r['Precursor Mz'])) >= 0.01]
print(f"Pre-processing -- after removing Calibration + precursor-mass transitions: {len(rows)} rows")


# ============================================================
# STEP 1 -- blank filter (fold-change threshold vs pooled blanks)
# ============================================================
blank_values = defaultdict(list)
for row in rows:
    if is_blank_replicate(row):
        blank_values[feature_key(row)].append(get_val(row, 'Area') or 0.0)
if not blank_values:
    raise SystemExit("No blank replicates found -- check 'sample type' or replicate naming "
                     "(see is_blank_replicate). Refusing to run without the blank filter.")
blank_mean = {k: statistics.mean(v) for k, v in blank_values.items()}
n_blank_reps = len(set(r['Replicate Name'] for r in rows if is_blank_replicate(r)))
print(f"Step 1 -- {n_blank_reps} blank replicates; blank means for {len(blank_mean)} "
      f"(list, molecule, fragment) combinations")

for row in rows:
    bmean = blank_mean.get(feature_key(row), 0.0)
    raw = get_val(row, 'Area') or 0.0
    row['_detected'] = raw > FOLD_CHANGE_CUTOFF * bmean if bmean > 0 else raw > 0

exp_rows = [r for r in rows if r.get('sample type') == 'Experimental']
n_det = sum(1 for r in exp_rows if r['_detected'])
print(f"           experimental rows: {len(exp_rows)}, passed >{FOLD_CHANGE_CUTOFF}x blank: "
      f"{n_det} ({100 * n_det / len(exp_rows):.1f}%)")


# ============================================================
# STEP 2 -- normalize to the dG internal standard
# ============================================================
dg_area = {}
for row in raw_rows:
    if row['Molecule List Name'] == DG_LIST_NAME and row['Molecule Name'] == DG_MOLECULE_NAME \
       and abs(float(row['Product Mz']) - DG_PRODUCT_MZ) < MZ_TOLERANCE:
        dg_area[row['Replicate Name']] = get_val(row, 'Area') or 0.0
exp_reps = sorted(set(r['Replicate Name'] for r in raw_rows if r.get('sample type') == 'Experimental'))
missing_dg = [rep for rep in exp_reps if not dg_area.get(rep)]
if missing_dg:
    raise SystemExit(f"dG standard missing or zero in experimental replicates: {missing_dg}")
print(f"Step 2 -- dG standard found for {len(dg_area)} replicates (all {len(exp_reps)} experimental)")

for row in rows:
    dg = dg_area.get(row['Replicate Name'])
    raw = get_val(row, 'Area') or 0.0
    row['_normalized_ratio'] = (raw / dg) if dg else None
    row['_group'] = assign_group(row.get('timepoint', ''), row.get('concentration', '')) \
        if row.get('sample type') == 'Experimental' else None


# ============================================================
# STEP 3 -- reproducibility: >= MIN replicates detected AND RSD < cutoff,
# per (feature, group), over the replicates that passed Step 1
# ============================================================
grouped = defaultdict(list)
for r in exp_rows:
    if r['_detected'] and r['_group'] in VALID_GROUPS and r['_normalized_ratio'] is not None \
       and r['Molecule Name'] != DG_MOLECULE_NAME:
        grouped[feature_key(r) + (r['_group'],)].append(r)

reproducible = []
for (list_name, mol_name, frag, group), members in grouped.items():
    n = len(members)
    if n < MIN_REPLICATES_FOR_RSD:
        continue
    ratios = [m['_normalized_ratio'] for m in members]
    rts = [float(m['Retention Time']) for m in members]
    mean_r, sd_r = statistics.mean(ratios), statistics.stdev(ratios)
    if mean_r == 0:
        continue
    rsd = sd_r / mean_r * 100
    if rsd < RSD_CUTOFF:
        reproducible.append({
            'MoleculeListName': list_name, 'MoleculeName': mol_name, 'Fragment': frag, 'Group': group,
            'N_replicates_detected': n,
            'RT_mean': round(statistics.mean(rts), 3),
            'RT_std': round(statistics.stdev(rts), 3),
            'Mean_Normalized_Ratio': round(mean_r, 6),
            'SD': round(sd_r, 6),
            'RSD_pct': round(rsd, 2),
            '_rt_std_raw': statistics.stdev(rts),
        })
print(f"Step 3 -- reproducibility (n>={MIN_REPLICATES_FOR_RSD}, RSD<{RSD_CUTOFF}%): "
      f"{len(grouped)} feature-group entries -> {len(reproducible)} rows, {n_features(reproducible)} features")


# ============================================================
# STEP 4 -- RT consistency
# (a) within group: RT SD across replicates <= RT_STD_CUTOFF
# (b) across groups: for a feature retained in several groups, each group's
#     mean RT must lie within CROSS_GROUP_RT_TOLERANCE of the median across
#     those groups. Groups where the feature is absent are not considered.
#     Note: with only two groups the median is their midpoint, so either
#     both rows are kept (RTs <= 2x tolerance apart) or both are dropped.
# ============================================================
n_before = len(reproducible)
reproducible = [r for r in reproducible if r['_rt_std_raw'] <= RT_STD_CUTOFF + RT_EPS]
print(f"Step 4a -- within-group RT SD <= {RT_STD_CUTOFF} min: {n_before} -> {len(reproducible)} rows")

feature_group_rt = defaultdict(list)
for r in reproducible:
    feature_group_rt[(r['MoleculeListName'], r['MoleculeName'], r['Fragment'])].append(r)

cross_group_consistent = []
for entries in feature_group_rt.values():
    median_rt = statistics.median(e['RT_mean'] for e in entries)
    cross_group_consistent += [e for e in entries
                               if abs(e['RT_mean'] - median_rt) <= CROSS_GROUP_RT_TOLERANCE + RT_EPS]
print(f"Step 4b -- cross-group RT within {CROSS_GROUP_RT_TOLERANCE} min of median: "
      f"{len(reproducible)} -> {len(cross_group_consistent)} rows")
reproducible = cross_group_consistent


# ============================================================
# STEP 5 -- peak shape: mean width <= MAX_PEAK_WIDTH_MIN AND width CV <
# WIDTH_CV_THRESHOLD. A peak can be consistently wide across replicates and
# still be poorly resolved/co-eluting -- CV alone does not catch that, so
# both criteria are required. Widths = End Time - Start Time over all
# experimental replicates of the group.
# ============================================================
width_lookup = defaultdict(list)
for row in raw_rows:
    if row.get('sample type') != 'Experimental':
        continue
    st, et = get_val(row, 'Start Time'), get_val(row, 'End Time')
    if st is None or et is None:
        continue
    grp = assign_group(row.get('timepoint', ''), row.get('concentration', ''))
    width_lookup[feature_key(row) + (grp,)].append(et - st)

peak_shape_ok, n_dropped_width, n_dropped_cv = [], 0, 0
for r in reproducible:
    widths = width_lookup.get((r['MoleculeListName'], r['MoleculeName'], r['Fragment'], r['Group']), [])
    if not widths:
        peak_shape_ok.append(r)          # no width data -- pass through
        continue
    mean_w = statistics.mean(widths)
    if mean_w > MAX_PEAK_WIDTH_MIN:
        n_dropped_width += 1
        continue
    r['Mean_Peak_Width_min'] = round(mean_w, 3)
    if len(widths) < 2:
        r['Peak_Width_CV_pct'] = ''
        peak_shape_ok.append(r)          # CV undefined with <2 widths -- pass on width alone
        continue
    cv_w = statistics.stdev(widths) / mean_w * 100 if mean_w > 0 else 0.0
    if cv_w < WIDTH_CV_THRESHOLD:
        r['Peak_Width_CV_pct'] = round(cv_w, 1)
        peak_shape_ok.append(r)
    else:
        n_dropped_cv += 1
print(f"Step 5 -- peak shape (mean width <= {MAX_PEAK_WIDTH_MIN} min, width CV < {WIDTH_CV_THRESHOLD}%): "
      f"dropped {n_dropped_width} for width, {n_dropped_cv} for CV "
      f"({len(reproducible)} -> {len(peak_shape_ok)} rows)")
reproducible = peak_shape_ok


# ============================================================
# STEP 6 -- correlation-based artifact removal, across BOTH molecule lists.
# Every pair of features whose mean RTs differ by <= CORRELATION_RT_TOLERANCE
# is tested; if their dG-normalized intensities across all experimental
# replicates correlate with r >= CORRELATION_MIN, both features are excluded.
# ============================================================
feature_rts = defaultdict(list)
for r in reproducible:
    feature_rts[(r['MoleculeListName'], r['MoleculeName'], r['Fragment'])].append(r['RT_mean'])
feature_rt_mean = {k: statistics.mean(v) for k, v in feature_rts.items()}

area_lookup = {}
for row in raw_rows:
    if row.get('sample type') == 'Experimental':
        area_lookup[feature_key(row) + (row['Replicate Name'],)] = get_val(row, 'Area') or 0.0
vectors = {k: [area_lookup.get(k + (rep,), 0.0) / dg_area[rep] for rep in exp_reps]
           for k in feature_rt_mean}

by_rt = sorted(feature_rt_mean, key=feature_rt_mean.get)
to_remove, n_pairs_tested, n_pairs_flagged = set(), 0, 0
for i, a in enumerate(by_rt):
    for b in by_rt[i + 1:]:
        if feature_rt_mean[b] - feature_rt_mean[a] > CORRELATION_RT_TOLERANCE + RT_EPS:
            break
        n_pairs_tested += 1
        if pearson(vectors[a], vectors[b]) >= CORRELATION_MIN:
            n_pairs_flagged += 1
            to_remove.update((a, b))

decorrelated = [r for r in reproducible
                if (r['MoleculeListName'], r['MoleculeName'], r['Fragment']) not in to_remove]
print(f"Step 6 -- {n_pairs_tested} co-eluting pairs tested, {n_pairs_flagged} with r >= {CORRELATION_MIN}: "
      f"removed {len(to_remove)} features ({len(reproducible)} -> {len(decorrelated)} rows, "
      f"{n_features(reproducible)} -> {n_features(decorrelated)} features)")


# ============================================================
# LIBRARY DOT PRODUCTS (informational) and OPTIONAL POST-PROCESSING (a):
# dot-product filter. NOT one of the six steps -- DISABLED by default.
# Mean dot products per (compound, group) are written to the output either way.
# ============================================================
APPLY_LIBRARY_DOTP_FILTER = False

dotp_by_compound_group = defaultdict(list)
for row in raw_rows:
    dotp = get_val(row, 'Library Dot Product')
    if dotp is None or row.get('sample type') != 'Experimental':
        continue
    grp = assign_group(row.get('timepoint', ''), row.get('concentration', ''))
    dotp_by_compound_group[(row['Molecule Name'], grp)].append(dotp)

for r in decorrelated:
    dotps = dotp_by_compound_group.get((r['MoleculeName'], r['Group']), [])
    r['Mean_Library_DotP'] = round(statistics.mean(dotps), 4) if dotps else ''

if APPLY_LIBRARY_DOTP_FILTER:
    qualifying = {r['MoleculeName'] for r in decorrelated
                  if not r['MoleculeName'].startswith('Win_')
                  and r['Mean_Library_DotP'] != '' and r['Mean_Library_DotP'] >= LIBRARY_DOTP_THRESHOLD}
    n_before = len(decorrelated)
    decorrelated = [r for r in decorrelated
                    if r['MoleculeName'].startswith('Win_') or r['MoleculeName'] in qualifying]
    print(f"Optional (a) -- library dotp >= {LIBRARY_DOTP_THRESHOLD}: {n_before} -> {len(decorrelated)} rows")
else:
    print("Optional (a) -- library dotp filter not applied (mean dot products reported in output)")


# ============================================================
# OPTIONAL POST-PROCESSING (b): resolve one known isomer ambiguity at
# m/z 167.056. NOT one of the six steps -- a one-off structural cleanup.
# 9-Methyl-guanine and N2-Methyl-guanine share this fragment mass and report
# identical values (same physical signal); only N2-Methyl-guanine is kept for
# this fragment. 9-Methyl-guanine's unambiguous fragment (149.0455) is untouched.
# ============================================================
n_before = len(decorrelated)
decorrelated = [r for r in decorrelated
                if not (r['MoleculeName'] == '9-Methyl-guanine' and '167.05' in r['Fragment'])]
print(f"Optional (b) -- isomer resolution at m/z 167.056: removed {n_before - len(decorrelated)} rows")


# ============================================================
# OUTPUT
# ============================================================
fieldnames = ['MoleculeListName', 'MoleculeName', 'Fragment', 'Group',
              'N_replicates_detected', 'RT_mean', 'RT_std', 'Mean_Peak_Width_min', 'Peak_Width_CV_pct',
              'Mean_Normalized_Ratio', 'SD', 'RSD_pct', 'Mean_Library_DotP']
with open(OUT_PATH, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
    w.writeheader()
    for r in decorrelated:
        r.setdefault('Mean_Peak_Width_min', '')
        r.setdefault('Peak_Width_CV_pct', '')
        w.writerow(r)

n_window = len({(r['MoleculeName'], r['Fragment']) for r in decorrelated if r['MoleculeName'].startswith('Win_')})
n_target = len({(r['MoleculeName'], r['Fragment']) for r in decorrelated if not r['MoleculeName'].startswith('Win_')})
print(f"\nFinal output: {len(decorrelated)} rows written to {OUT_PATH}")
print(f"  Diagnostic-ion window features: {n_window}")
print(f"  Library-matched target features: {n_target}")
