#!/usr/bin/env bash
set -euo pipefail

mkdir -p results

out="results/results.csv"
echo "molecule,E_HF_total,E_MP2_corr,E_MP2_total" > "$out"

for f in data/*.h5; do
  mol=$(basename "$f" .h5)

  hf=$(./mp2_energy "$f" | awk '/E_HF_total:/ {print $2}')
  mp2c=$(./mp2_energy "$f" | awk '/E_MP2_corr:/ {print $2}')
  mp2t=$(./mp2_energy "$f" | awk '/E_MP2_total:/ {print $2}')

  echo "$mol,$hf,$mp2c,$mp2t" >> "$out"
done

echo "Wrote $out"

