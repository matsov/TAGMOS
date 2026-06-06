"""
TAGMOS verify_calibration · sanity-check a calibration JSON

Checks:
  1. SHA-256 hash matches the recomputed hash of the calibration content
  2. Engine z quartile cuts are strictly monotonic (Q25 < Q50 < Q75)
  3. Each tertile cut is strictly monotonic (Q33 < Q66)
  4. Each sub-axis has finite mu and strictly positive sigma
  5. The 14 framework axes are all represented
  6. No L0 / Wellmicro proprietary fields are present in the JSON

Usage:
    python verify_calibration.py --calibration my_calibration.json
"""

import argparse, json, hashlib, sys
import math


FRAMEWORK_AXES_14 = ["ENG", "SYN_cof", "SYN_carb",
                     "BA", "TRP", "PROT", "MUC", "HIS", "TMA",
                     "UREM", "IRON", "LPS", "ETU", "HYS"]

FORBIDDEN_KEYS = [
    "wellmicro_italian_rwe", "italian_6508", "italian_5272",
    "L0_calibration", "L0_proprietary", "wmp_internal",
    "frozen_anchor_wellmicro", "internal_calibration_wmp",
]


def fail(msg):
    print(f"  ✗ {msg}")
    return 1


def ok(msg):
    print(f"  ✓ {msg}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--calibration", required=True)
    args = ap.parse_args()

    with open(args.calibration) as f:
        cal = json.load(f)

    errors = 0
    print(f"[TAGMOS verify] verifying: {args.calibration}")
    print(f"[TAGMOS verify] calibration_id: {cal.get('calibration_id', '?')}")
    print()

    # 1. SHA-256 hash check
    stored_hash = cal.pop("calibration_sha256", None)
    if stored_hash is None:
        errors += fail("missing calibration_sha256 field")
    else:
        recomputed = hashlib.sha256(
            json.dumps(cal, sort_keys=True).encode("utf-8")
        ).hexdigest()
        if recomputed == stored_hash:
            ok(f"SHA-256 hash matches ({stored_hash[:16]}...)")
        else:
            errors += fail(f"SHA-256 hash mismatch · stored={stored_hash[:16]}... "
                            f"recomputed={recomputed[:16]}...")
    cal["calibration_sha256"] = stored_hash

    # 2. Engine z quartile cuts
    eng_q = cal.get("engine_quartile_cuts", {})
    q25, q50, q75 = eng_q.get("Q25"), eng_q.get("Q50"), eng_q.get("Q75")
    if None in (q25, q50, q75):
        errors += fail("missing one or more Engine z quartile cuts")
    elif not (q25 < q50 < q75):
        errors += fail(f"Engine z cuts not monotonic: Q25={q25}, Q50={q50}, Q75={q75}")
    else:
        ok(f"Engine z quartile cuts monotonic (Q25 < Q50 < Q75)")

    # 3. tertile cuts monotonicity per axis
    tert = cal.get("tertile_cuts", {})
    n_tert_axes = 0
    for axis, cuts in tert.items():
        q33, q66 = cuts.get("Q33"), cuts.get("Q66")
        if None in (q33, q66) or not (q33 < q66):
            errors += fail(f"{axis} tertile cuts not monotonic: Q33={q33}, Q66={q66}")
        else:
            n_tert_axes += 1
    if n_tert_axes > 0:
        ok(f"tertile cuts monotonic for {n_tert_axes} axes")

    # 4. sub-axis mu/sigma checks
    sub = cal.get("per_sub_axis_mu_sigma", {})
    bad = []
    for sa, vals in sub.items():
        mu, sigma = vals
        if not math.isfinite(mu):
            bad.append(f"{sa}: mu not finite ({mu})")
        if sigma <= 0:
            bad.append(f"{sa}: sigma not positive ({sigma})")
    if bad:
        for b in bad[:5]:
            errors += fail(b)
        if len(bad) > 5:
            errors += fail(f"... {len(bad)-5} more sub-axis mu/sigma issues")
    else:
        ok(f"all {len(sub)} sub-axes have finite mu and positive sigma")

    # 5. framework axes coverage
    axis_mu = cal.get("per_axis_mu_sigma", {})
    missing = [a for a in FRAMEWORK_AXES_14 if a not in axis_mu]
    if missing:
        errors += fail(f"missing framework axes in per_axis_mu_sigma: {missing}")
    else:
        ok(f"all 14 framework axes present in per_axis_mu_sigma")

    # 6. forbidden keys check
    keys = json.dumps(cal).lower()
    found_forbidden = [k for k in FORBIDDEN_KEYS if k.lower() in keys]
    if found_forbidden:
        errors += fail(f"forbidden L0 / Wellmicro keys present: {found_forbidden}")
    else:
        ok(f"no forbidden L0 / Wellmicro keys present")

    # done
    print()
    if errors == 0:
        print(f"[TAGMOS verify] ✅ PASS · {args.calibration} is internally consistent and L0-clean")
        sys.exit(0)
    else:
        print(f"[TAGMOS verify] ✗ {errors} issues detected · review the calibration and re-run")
        sys.exit(1)


if __name__ == "__main__":
    main()
