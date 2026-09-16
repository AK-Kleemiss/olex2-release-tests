#!/usr/bin/env python3
"""Compare a release-test log with its golden.

A log is what Olex2's run_pipeline.py writes for the "release" group: one
line per case, sorted by name,

    PASS ptb_ecp_malbac  tsc=malbac.tscb cycles=1 R1_sph=0.0312 R1_asp=0.0290 ...
    SKIP salted_v6_sucrose  no .salted file ...
    END 12 0 2

and a final END line with the pass/fail/skip counts. The golden has the same
shape. Cases are matched by name, never by position, so a new or missing case
is reported by name and does not shift everything after it.

Two lines match when they are identical, or when their skeletons (numbers
replaced by <num>, whitespace collapsed) are identical and every number pair
agrees: integer tokens exactly, tokens with a decimal point or exponent within
FLOAT_TOL (0.0005, the R1 tolerance Florian set on 16 Sep 2026). Names like
NA2_0010000.tscb or 12x14x16 are not numbers - a digit run touching a letter
or underscore stays part of the word.

Exit codes: 0 all cases match and none skipped, 2 all match but at least one
SKIP (the banner says which and why), 1 a FAIL line, a mismatch, a case that
only one side has, or a log without END (the run did not complete).
"""
import argparse
import re
import sys

FLOAT_TOL = 0.0005
NUM = re.compile(r"(?<![A-Za-z_.])[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?(?![A-Za-z_.\d])")


def parse(path):
    """{case: (state, detail)}, end tuple or None, in-order case names."""
    cases, order, end = {}, [], None
    with open(path, encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            parts = line.split(None, 2)
            if parts[0] == "END":
                end = tuple(int(x) for x in line.split()[1:4])
                continue
            if parts[0] not in ("PASS", "FAIL", "SKIP") or len(parts) < 2:
                continue
            detail = " ".join(parts[2].split()) if len(parts) > 2 else ""
            cases[parts[1]] = (parts[0], detail)
            order.append(parts[1])
    return cases, end, order


def numbers(text):
    return [m.group(0) for m in NUM.finditer(text)]


def skeleton(text):
    return " ".join(NUM.sub("<num>", text).split())


def is_float(token):
    return "." in token or "e" in token.lower()


def compare_detail(want, got):
    """None when they agree, else a short reason."""
    if want == got:
        return None
    if skeleton(want) != skeleton(got):
        return "text differs"
    for a, b in zip(numbers(want), numbers(got)):
        if is_float(a) or is_float(b):
            if abs(float(a) - float(b)) > FLOAT_TOL:
                return "%s vs %s (|d| %.4g > %.4g)" % (a, b, abs(float(a) - float(b)), FLOAT_TOL)
        elif a != b:
            return "%s vs %s (integers must match)" % (a, b)
    return None


def compare(golden_path, log_path, out=sys.stdout, only=None):
    """Print the report; return the exit code.

    only: names of the cases that were run (a --cases subset); the other golden
    cases are then not reported as missing. A subset run is never a release
    result, so the report says so."""
    want, want_end, want_order = parse(golden_path)
    got, got_end, got_order = parse(log_path)
    problems, skips, matched = [], [], 0
    if only:
        only = set(only) | {"setup"}
        want_order = [n for n in want_order if n in only]

    if got_end is None:
        problems.append("the log has no END line: the run did not complete")

    for name in want_order:
        if name not in got:
            problems.append("%-28s missing from the log" % name)
    for name in got_order:
        if name not in want:
            problems.append("%-28s not in the golden (new case? --update)" % name)

    for name in want_order:
        if name not in got:
            continue
        w_state, w_detail = want[name]
        g_state, g_detail = got[name]
        if g_state == "FAIL":
            problems.append("%-28s FAIL  %s" % (name, g_detail))
        elif g_state == "SKIP":
            skips.append((name, g_detail))
        elif w_state == "SKIP":
            # golden recorded a skip, now it ran: no reference to compare with
            problems.append("%-28s ran (%s) but the golden holds a SKIP - rerun with "
                            "everything installed and --update" % (name, g_detail))
        else:
            why = compare_detail(w_detail, g_detail)
            if why is None:
                matched += 1
            else:
                problems.append("%-28s %s\n%30s golden: %s\n%30s log:    %s"
                                % (name, why, "", w_detail, "", g_detail))

    print("golden %s" % golden_path, file=out)
    print("log    %s" % log_path, file=out)
    print("%d of %d cases match%s" % (matched, len(want_order),
                                       " (subset run, not a release result)" if only else ""), file=out)
    if problems:
        print("", file=out)
        print("PROBLEMS", file=out)
        for p in problems:
            print("  " + p, file=out)
    if skips:
        print("", file=out)
        print("#" * 72, file=out)
        print("#  SKIPPED - these pathways were NOT tested on this machine", file=out)
        print("#" * 72, file=out)
        for name, why in skips:
            print("#  %-28s %s" % (name, why), file=out)
        print("#" * 72, file=out)
    if problems:
        return 1
    return 2 if skips else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("golden")
    ap.add_argument("log")
    ap.add_argument("--update", action="store_true",
                    help="accept the log as the new golden (copied verbatim)")
    ap.add_argument("--only", help="comma-separated case names that were run; the rest of the golden is not expected")
    args = ap.parse_args(argv)
    if args.update:
        with open(args.log, encoding="utf-8") as src, open(args.golden, "w", encoding="utf-8", newline="\n") as dst:
            dst.write(src.read().replace("\r\n", "\n"))
        print("golden %s updated from %s" % (args.golden, args.log))
        return 0
    return compare(args.golden, args.log,
                   only=[c.strip() for c in args.only.split(",") if c.strip()] if args.only else None)


if __name__ == "__main__":
    sys.exit(main())
