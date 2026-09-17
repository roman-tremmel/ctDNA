#!/usr/bin/env bash
# Print the count-TSV path for each sample listed in a plain text file -
# for building the --counts argument to build_control_baseline.py (or
# compute_aneuploidy_score.py for a batch of cases) without relying on
# sample names following any particular naming convention (e.g. a
# "control_*" prefix, which real sample IDs like "V91-01" won't have).
#
# Usage: counts_for_samples.sh <results_dir> <suffix> <sample_list_file>
#   results_dir:      the $RUNDIR/results you passed to run_pipeline.sh
#   suffix:            .arm_counts.tsv (mfast-seqs) or .window_counts.tsv (waldo)
#   sample_list_file:  one sample name per line (blank lines and lines
#                      starting with # are skipped) - e.g. a controls.txt
#                      you maintain by hand, listing which of your sample
#                      IDs are healthy controls
#
# Example:
#   printf 'V91-02\nV91-07\nV91-11\n' > "$RUNDIR/controls.txt"
#   python3 build_control_baseline.py \
#     --counts $("$REPO/mfast-seqs/scripts/counts_for_samples.sh" \
#                  "$RUNDIR/results" .arm_counts.tsv "$RUNDIR/controls.txt") \
#     --out "$RUNDIR/results/baseline.tsv"
set -euo pipefail

RESULTS_DIR="$1"
SUFFIX="$2"
SAMPLE_LIST="$3"

missing=0
while IFS= read -r sample; do
  sample="${sample%%$'\r'}"   # strip trailing \r if the list was edited on Windows
  [[ -z "$sample" || "$sample" == \#* ]] && continue
  path="$RESULTS_DIR/$sample/$sample$SUFFIX"
  if [[ ! -f "$path" ]]; then
    echo "WARNING: missing $path (sample '$sample' not yet run, or wrong name?)" >&2
    missing=1
    continue
  fi
  echo "$path"
done < "$SAMPLE_LIST"

if [[ "$missing" == 1 ]]; then
  echo "WARNING: one or more listed samples had no count file (see above) - they were skipped, not included as empty." >&2
fi
