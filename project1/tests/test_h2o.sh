#!/usr/bin/env bash
set -euo pipefail

exe="./mp2_energy"
file="data/h2o.h5"

if [[ ! -x "$exe" ]]; then
  echo "ERROR: mp2_energy not built. Run: make"
  exit 1
fi

if [[ ! -f "$file" ]]; then
  echo "ERROR: missing $file (data files are ignored by git)."
  exit 1
fi

out=$($exe "$file")

hf=$(echo "$out"   | awk '/E_HF_total:/  {print $2}')
mp2c=$(echo "$out" | awk '/E_MP2_corr:/  {print $2}')
mp2t=$(echo "$out" | awk '/E_MP2_total:/ {print $2}')

hf_ref="-76.0267987082"
mp2c_ref="-0.2039599741"
mp2t_ref="-76.2307586823"
tol="1e-6"

python3 - << PY
import sys
hf=float("$hf"); mp2c=float("$mp2c"); mp2t=float("$mp2t")
hf_ref=float("$hf_ref"); mp2c_ref=float("$mp2c_ref"); mp2t_ref=float("$mp2t_ref")
tol=float("$tol")

def ok(x, ref): return abs(x-ref) <= tol

errs=[]
if not ok(hf, hf_ref): errs.append(f"HF mismatch: {hf} vs {hf_ref}")
if not ok(mp2c, mp2c_ref): errs.append(f"MP2 corr mismatch: {mp2c} vs {mp2c_ref}")
if not ok(mp2t, mp2t_ref): errs.append(f"MP2 total mismatch: {mp2t} vs {mp2t_ref}")

if errs:
    print("TEST FAILED")
    for e in errs: print(" -", e)
    sys.exit(1)

print("TEST OK (H2O)")
PY
