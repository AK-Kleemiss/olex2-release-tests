# Golden logs

`release_tests.quick.good` and `release_tests.full.good` are the reference
logs, one line per case (sorted) plus `END <pass> <fail> <skip>`. They were
produced on FLOWOFFICE (Windows 11) on 16 Sep 2026 with the NoSpherA2 bundle of
commit 0442ebf1 (XCW_Test), ORCA 6.1, and SALTED Model_V6.

Accept a new golden only after the diff is understood:

    python run_release_tests.py ... --keep-scratch      # look at the case
    python run_release_tests.py ... --update-golden     # then accept

A per-platform file `release_tests.<tier>.<win|linux|mac>.good` overrides the
shared one when present; add one only if a platform drifts beyond 0.0005 in R1.
