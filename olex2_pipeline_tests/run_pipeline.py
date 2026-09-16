"""Entry point, run inside olex2c with py.Run.

  py.Run '<rundir>/util/pyUtil/regression/olex2_pipeline_tests/run_pipeline.py'

See pipeline_tests.py for what the environment variables select. The result
table is printed and also written to OLEX2_TEST_OUT, so a caller can gate a
commit on it without parsing the whole log - olex2c's output is UTF-16 and
carries everything Olex2 says on the way past.
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import time

import olx

# py.Run executes the file without setting __file__, so the location has to
# come from Olex2 rather than from the module
try:
  _here = os.path.dirname(os.path.abspath(__file__))
except NameError:
  _here = os.path.join(olx.BaseDir(), "util", "pyUtil", "regression",
                       "olex2_pipeline_tests")
if _here not in sys.path:
  sys.path.insert(0, _here)

import pipeline_tests
from pipeline_tests import Suite, Result

GROUPS = {
  "api":       "group_api",
  "macros":    "group_macros",
  "model":     "group_model",
  "restraints": "group_restraints",
  "instructions": "group_instructions",
  "validation": "group_validation",
  "software":  "group_software",
  "geometry":  "group_geometry",
  "symmetry":  "group_symmetry",
  "formats":   "group_formats",
  "hkl":       "group_hkl",
  "cif":       "group_cif",
  "solve":     "group_solve",
  "refine":    "group_refine",
  "nosphera2": "group_nosphera2",
  "nsa2_matrix": "group_nsa2_matrix",
  "release":   "group_release",
}
DEFAULT = ("api,macros,model,instructions,restraints,geometry,symmetry,"
           "formats,hkl,cif,solve,refine")


def _write_release_log(suite, path, complete):
  """One line per case sorted by name, then END pass fail skip; and next
  to it <stem>_timings.txt with the seconds that the golden must not hold."""
  if not path:
    return
  lines = sorted(r.log_line() for r in suite.results)
  if complete:
    n = {Result.PASS: 0, Result.FAIL: 0, Result.SKIP: 0}
    for r in suite.results:
      n[r.state] = n.get(r.state, 0) + 1
    lines.append("END %d %d %d" % (n[Result.PASS], n[Result.FAIL], n[Result.SKIP]))
  try:
    with open(path, "w") as f:
      f.write("\n".join(lines) + "\n")
    stem, ext = os.path.splitext(path)
    with open(stem + "_timings.txt", "w") as f:
      for r in sorted(suite.results, key=lambda r: r.name):
        f.write("%s %.1f\n" % (r.name, r.seconds))
    if complete:
      print("release log written to %s" % path)
  except Exception as e:
    print("could not write %s: %s" % (path, e))


def main():
  wanted = os.environ.get("OLEX2_TEST_GROUPS", DEFAULT)
  wanted = [w.strip() for w in wanted.split(",") if w.strip()]
  out = os.environ.get("OLEX2_TEST_OUT",
                       os.path.join(os.path.expanduser("~"), "olex2_pipeline_tests.txt"))

  print("")
  print("=" * 78)
  print("Olex2 pipeline tests: %s" % ", ".join(wanted))
  print("=" * 78)
  sys.stdout.flush()

  suite = Suite()
  t0 = time.time()
  try:
    for name in wanted:
      mod_name = GROUPS.get(name)
      if mod_name is None:
        print("no such group: %s (have %s)" % (name, ", ".join(sorted(GROUPS))))
        continue
      try:
        mod = __import__(mod_name)
      except Exception as e:
        # a group that will not even import is a failure of that group, not of
        # the run - the others still have something to say
        suite.results.append(Result(name, "import", Result.FAIL, str(e)))
        print("FAIL  %-10s import  %s" % (name, e))
        continue
      mod.register(suite)
  finally:
    # the release gate reads this file; a crash half-way must still leave
    # one, without the END line that says the run completed
    _write_release_log(suite, os.environ.get("OLEX2_TEST_RELEASE_LOG"),
                       complete=False)

  n = suite.summary()
  _write_release_log(suite, os.environ.get("OLEX2_TEST_RELEASE_LOG"), complete=True)
  elapsed = time.time() - t0
  lines = [str(r) for r in suite.results]
  lines.append("")
  lines.append("%d passed, %d failed, %d skipped in %.0fs"
               % (n[Result.PASS], n[Result.FAIL], n[Result.SKIP], elapsed))

  print("")
  print("-" * 78)
  for line in lines[-1:]:
    print(line)
  if n[Result.SKIP]:
    print("")
    print("skipped, and why:")
    for r in suite.results:
      if r.state == Result.SKIP:
        print("  %-34s %s" % (r.name, r.detail))
  if n[Result.FAIL]:
    print("")
    print("failed:")
    for r in suite.results:
      if r.state == Result.FAIL:
        print("  %-34s %s" % (r.name, r.detail))
  print("-" * 78)

  try:
    with open(out, "w") as f:
      f.write("\n".join(lines) + "\n")
    print("table written to %s" % out)
  except Exception as e:
    print("could not write %s: %s" % (out, e))
  sys.stdout.flush()
  return 1 if n[Result.FAIL] else 0


main()
