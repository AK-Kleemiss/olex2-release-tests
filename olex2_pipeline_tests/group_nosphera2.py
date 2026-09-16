"""NoSpherA2, across the backends that produce form factors.

Each backend gets the same shape of test: refine spherically for a baseline,
then run **one** aspherical refinement and check that R1 moved and did not get
worse. A tsc that is produced but never reaches the refinement looks identical
to one that works, and the R1 is the only thing that tells them apart.

The aspherical run is a plain `spy.refine.do_refine`, not a hand-built
sequence. With `use_aspherical` set and olex2.refine selected,
RunRefinementPrg hands the whole thing to `aaff.deal_with_AAFF`, which is the
production HAR workflow: it sets the hydrogen treatment, calls
`NoSpherA2_instance.launch()`, requires a table to appear, resolves the table's
columns against the model, validates the table's hash and origin, refines, and
loops on the shifts until it converges or hits Max_HAR_Cycles.

Driving that instead of orchestrating the steps here is the point of these
tests. It covers the parts of the pipeline a hand-built sequence skips
entirely - the hash and origin checks, the scatterer resolution, the
convergence loop - and it means disorder is handled the way it is in
production, by one call that hands every part's wavefunction to `-mtc`, rather
than by this file deciding how parts should be combined.

`deal_with_AAFF` keeps its own record in `<name>.NoSpherA2`, a line per cycle
with R1 and wR2, so that file existing with a cycle in it is the evidence the
workflow ran rather than being skipped.

The quantum backends are minutes to hours, so they are behind OLEX2_TEST_FULL.
discamb is a databank and SALTED a trained model, so both are quick.

A backend that is not offered here is SKIPped by name, with the list that was
offered. It is never dropped silently: this suite exists to say what was
covered, and "no ORCA here" is a different statement from "ORCA works".
"""
from __future__ import absolute_import, division, print_function

import os
import sys

import olx
import olex
from olexFunctions import OV

from variableFunctions import nsa2_get_param, nsa2_set_param

from pipeline_tests import (macro, SkipTest, load, model_in, has_hkl,
                            r1_of_last_refinement)

# Neutral and closed shell, which charge 0 and multiplicity 1 below assert. The
# sample called water is a Mn complex - its SCF diverged within five cycles
# under those settings, because they are the wrong ones for it.
SAMPLE = "sucrose"

# carries PART, so NoSpherA2 has to take the multi-part route
MULTIPART_SAMPLE = "THPP"

# label -> is it a quantum backend, and so slow enough to need OLEX2_TEST_FULL.
# Whether a backend can run is not asked here: _source_name settles it from
# what NoSpherA2 is offering. Testing for an executable got OCC wrong - it needs
# none, the job writes a toml and NoSpherA2.exe reads it.
BACKENDS = [
  ("DiscaMB", False),
  ("SALTED",  False),
  ("xTB",     False),   # semi-empirical, so seconds rather than minutes
  ("pTB",     False),
  ("ORCA",    True),
  ("OCC",     True),
]

# For the disorder case, which needs a wavefunction backend but not an
# expensive one: ORCA on a three-part structure is three full SCFs and ran
# past ten minutes twice. Tried in order.
MULTIPART_BACKENDS = ["xTB", "pTB", "ORCA"]


def _offered():
  """The backends NoSpherA2 is offering, as the settings dropdown sees them."""
  listed = olex.f("spy.NoSpherA2.getwfn_softwares()") or ""
  return [s.strip() for s in str(listed).split(";") if s.strip()]


def _source_name(label):
  """The string NoSpherA2 compares its source against, resolved not guessed.

  None of these names is a fixed constant. discamb is compared against the
  *value* of user.NoSpherA2.discamb_exe, and ORCA carries its version -
  "ORCA 6.1" here, "ORCA 5.0" or bare "ORCA" elsewhere, from what orca -v
  reported at startup. Both of the obvious labels fell past every branch into
  Wfn_Job, which refused them as wavefunction programs.

  So take the name from the list the dropdown is built from: an entry that
  starts with the label is this machine's spelling of it, and no entry means
  the backend is not offered here, which is a skip rather than a failure.
  """
  if label == "DiscaMB":
    label = str(OV.GetParam('user.NoSpherA2.discamb_exe', '') or 'discambMATTS2tsc')
  for entry in _offered():
    # "Get ORCA" is the download prompt, not an installed backend
    if entry.startswith("Get "):
      continue
    if entry == label or entry.startswith(label):
      return entry
  raise SkipTest("%s is not offered by NoSpherA2 here (offering: %s)"
                 % (label, ", ".join(_offered()) or "nothing"))


def register(suite):
  full = os.environ.get("OLEX2_TEST_FULL", "") == "1"
  # A quantum backend is minutes; re-running the others to reach the one being
  # worked on is not free. OLEX2_TEST_BACKENDS names the ones to run.
  wanted = [s.strip().lower()
            for s in os.environ.get("OLEX2_TEST_BACKENDS", "").split(",")
            if s.strip()]
  suite.run("nsa2", "NoSpherA2 is present", t_present)
  # The disorder case runs a wavefunction backend too, and on a three-part
  # structure, so it obeys the same selection - asking for one backend and
  # getting three ORCA jobs as well is a surprise, and a slow one.
  if not wanted or "multipart" in wanted:
    suite.run("nsa2", "disorder goes through -mtc, not a merge",
              t_multipart, suite, full)
  for name, slow in BACKENDS:
    if wanted and name.lower() not in wanted:
      continue
    suite.run("nsa2", "tsc and refine: %s" % name,
              t_backend, suite, name, slow, full)


def t_present():
  exe = os.path.join(OV.BaseDir(),
                     "NoSpherA2" + (".exe" if sys.platform.startswith("win") else ""))
  if not os.path.exists(exe):
    raise SkipTest("NoSpherA2 executable is not in the run directory")
  return os.path.basename(exe)


def _available(name):
  """Resolve the backend and select what it needs, or skip saying why."""
  source = _source_name(name)
  if name == "SALTED":
    _select_salted_model()
  elif name == "OCC":
    # Offered unconditionally by setup_occ_executables, but OCC lives inside
    # NoSpherA2.exe and needs its data beside it. Without it the job is built,
    # launched, and dies with "OCC_DATA_PATH not set or invalid" - and
    # utilities.py sets that variable whether the directory is there or not,
    # so the variable being present says nothing. This is the same test
    # NoSpherA2 makes at startup.
    basis = os.path.join(OV.BaseDir(), "occ", "basis")
    if not os.path.isdir(basis):
      raise SkipTest("OCC data is not installed (%s)" % basis)
  return source


def _select_salted_model():
  """Pick a model, the way the settings page would.

  A trained model is a separate download, so no model is a skip rather than a
  failure - but a model that is installed and simply not chosen is the state a
  fresh run directory is in, and the GUI is the only thing that ever chooses
  one. The list is user.NoSpherA2.salted_models_list, the selection is
  snum.NoSpherA2.selected_salted_model.
  """
  if str(nsa2_get_param('selected_salted_model') or '').strip():
    return
  listed = str(OV.GetParam('user.NoSpherA2.salted_models_list', '') or '')
  for cand in listed.replace(",", ";").split(";"):
    cand = cand.strip()
    if cand and os.path.isdir(cand):
      nsa2_set_param('selected_salted_model', cand)
      return
  raise SkipTest("no SALTED model is installed (user.NoSpherA2.salted_models_list)")


def t_backend(suite, name, slow, full):
  source = _available(name)
  if slow and not full:
    raise SkipTest("quantum backend, set OLEX2_TEST_FULL=1 to include it")

  folder = suite.sample(SAMPLE)
  if not has_hkl(folder):
    raise SkipTest("no hkl with the %s sample" % SAMPLE)
  model = model_in(folder)
  load(model)
  macro("user '%s'" % folder.replace("\\", "/"))

  # spherical first, as the thing the aspherical result has to beat
  _refine_spherically()
  r_spherical = r1_of_last_refinement()

  _set_aspherical(source)
  said = _refine_capturing()
  r_aspherical = r1_of_last_refinement()

  cycles = _har_log_cycles(folder)
  if not cycles:
    raise AssertionError("%s: deal_with_AAFF left no HAR log - the aspherical "
                         "path was not taken" % name)
  tsc = _tsc_beside(folder)
  if tsc is None:
    # say why, not just that. pTB runs to completion and hands back a
    # wavefunction with zero MO occupations; "produced no tsc" hid that.
    raise AssertionError("%s produced no tsc or tscb%s"
                         % (name, _reason(said)))

  # An unchanged R1 is the failure this test exists to catch: the tsc is
  # written, the refinement runs, and the form factors are never used. It
  # happened here - the label resync read the wrong header line, reported it,
  # and the refinement continued spherically to exactly the same R1.
  if abs(r_aspherical - r_spherical) < 1e-5:
    raise AssertionError(
      "%s: R1 %.4f unchanged by the tsc - the table did not reach the "
      "refinement" % (name, r_aspherical))
  if r_aspherical > r_spherical + 0.01:
    raise AssertionError(
      "%s: aspherical R1 %.4f is worse than spherical %.4f"
      % (name, r_aspherical, r_spherical))
  return "%s, %d HAR cycle(s), R1 %.4f -> %.4f" % (
    os.path.basename(tsc), cycles, r_spherical, r_aspherical)


def _refine_capturing():
  """Run the refinement and keep what it printed, still showing it."""
  import io
  import sys
  buf = io.StringIO()
  real_out, real_err = sys.stdout, sys.stderr

  class tee(object):
    def __init__(self, a, b):
      self.a, self.b = a, b

    def write(self, s):
      try:
        self.a.write(s)
      except Exception:
        pass
      self.b.write(s)

    def flush(self):
      try:
        self.a.flush()
      except Exception:
        pass

    def __getattr__(self, n):
      return getattr(self.a, n)

  sys.stdout, sys.stderr = tee(real_out, buf), tee(real_err, buf)
  try:
    macro("spy.refine.do_refine")
  finally:
    sys.stdout, sys.stderr = real_out, real_err
  return buf.getvalue()


def _reason(text):
  """The backend's own complaint, if it made one."""
  for line in str(text).splitlines():
    l = line.strip()
    if not l or len(l) > 160:
      continue
    low = l.lower()
    if ("error" in low or "failed" in low or "empty" in low
        or "aborted" in low) and "traceback" not in low:
      return " - %s" % l
  return ""


def _set_aspherical(source):
  """The settings deal_with_AAFF needs, which only the settings page ever sets.

  charge and multiplicity have no usable defaults - multiplicity ships as 0,
  which is even, and launch refuses an even multiplicity with an even electron
  count.
  """
  nsa2_set_param('charge', '0')
  nsa2_set_param('multiplicity', '1')
  nsa2_set_param('source', source)
  nsa2_set_param('full_HAR', False)     # one cycle; the loop is exercised, not iterated
  nsa2_set_param('Max_HAR_Cycles', '1')
  nsa2_set_param('run_refine', True)
  nsa2_set_param('make_fcf_only', False)
  nsa2_set_param('use_aspherical', True)


def _har_log_cycles(folder):
  """Cycles recorded in deal_with_AAFF's own log, which it writes per run."""
  n = 0
  for root, dirs, files in os.walk(folder):
    for f in files:
      if not f.endswith(".NoSpherA2"):
        continue
      for line in open(os.path.join(root, f), "r", errors="ignore"):
        # 'run' is written left-padded to 3, so a cycle line starts with digits
        if line[:3].strip().isdigit() and int(line[:3]) > 0:
          n += 1
  return n


def t_multipart(suite, full):
  """A disordered structure gets one table, from one run.

  NoSpherA2 passes every part's wavefunction to a single cuqct_tsc call, which
  becomes `-mtc`, and gets back one table covering all of them. The other route
  - a tsc per part and `-merge` afterwards - is the one not to take. So the
  assertion is that a table came back and that no per-part tables were left
  behind to be merged.

  This needs a wavefunction backend. `-mtc` takes wavefunctions, and discamb
  emits finished tables instead, so a discamb run on a disordered structure
  has nothing to hand `-mtc` and merges by construction: THPP came back as
  thpp_part_1.tsc and thpp_part_2.tsc.
  """
  if not full:
    raise SkipTest("needs a wavefunction backend, set OLEX2_TEST_FULL=1")
  source, tried = None, []
  for label in MULTIPART_BACKENDS:
    try:
      source = _available(label)
      break
    except SkipTest as e:
      tried.append("%s (%s)" % (label, e))
  if source is None:
    raise SkipTest("no wavefunction backend available: %s" % "; ".join(tried))
  folder = suite.sample(MULTIPART_SAMPLE)
  if not has_hkl(folder):
    raise SkipTest("no hkl with the %s sample" % MULTIPART_SAMPLE)
  model = model_in(folder)
  parts = set()
  for line in open(model, "r", errors="ignore"):
    if line[:5].strip().upper() == "PART":
      parts.add(line.split()[1])
  if len(parts) < 2:
    raise SkipTest("%s no longer carries more than one part" % MULTIPART_SAMPLE)
  load(model)
  macro("user '%s'" % folder.replace("\\", "/"))
  # refine first, and not only for a baseline: launch does
  # os.mkdir('olex2\\Wfn_job'), a relative two-level path, so the olex2 folder
  # has to exist already and it is the refinement that creates it
  _refine_spherically()
  r_spherical = r1_of_last_refinement()

  _set_aspherical(source)
  macro("spy.refine.do_refine")

  if not _har_log_cycles(folder):
    raise AssertionError("deal_with_AAFF left no HAR log for a %d-part "
                         "structure" % len(parts))
  # The flag itself is never printed, so this asserts its signature: one
  # wavefunction per part, one table covering them all, and nothing left over
  # for a -merge. The merge route would leave a table per part instead.
  tsc = _tsc_beside(folder)
  if tsc is None:
    raise AssertionError("no tsc for a %d-part structure" % len(parts))
  leftovers = _per_part_tables(folder)
  if leftovers:
    raise AssertionError(
      "per-part tables were left to be merged: %s"
      % ", ".join(sorted(os.path.basename(p) for p in leftovers)))
  wfns = _wavefunctions(folder)
  if len(wfns) < 2:
    raise AssertionError(
      "a %d-part structure produced %d wavefunction(s) - the parts were not "
      "computed separately" % (len(parts), len(wfns)))
  tables = _all_tables(folder)
  if len(tables) != 1:
    raise AssertionError("expected one table, found %d: %s"
                         % (len(tables),
                            ", ".join(sorted(os.path.basename(t) for t in tables))))
  r_aspherical = r1_of_last_refinement()
  if abs(r_aspherical - r_spherical) < 1e-5:
    raise AssertionError("R1 %.4f unchanged - the table did not reach the "
                         "refinement" % r_aspherical)
  return "parts %s -> one %s, R1 %.4f -> %.4f" % (
    ",".join(sorted(parts)), os.path.basename(tsc), r_spherical, r_aspherical)


def _refine_spherically():
  nsa2_set_param('use_aspherical', False)
  OV.SetParam('snum.refinement.program', 'olex2.refine')
  # set explicitly: a SHELX method left over from an earlier case is not one
  # olex2.refine accepts, and the set is refused rather than corrected
  OV.SetParam('snum.refinement.method', 'Gauss-Newton')
  OV.SetParam('snum.refinement.max_cycles', 3)
  macro("spy.refine.do_refine")


def _wavefunctions(folder):
  """One per part is what the multi-part route produces before the table."""
  found = set()
  for root, dirs, files in os.walk(folder):
    for f in files:
      # .xtb is what pTB writes (thpp.xtb per part), .molden what xTB writes
      if f.lower().endswith((".molden", ".gbw", ".wfn", ".wfx", ".fchk",
                             ".xtb")):
        found.add(os.path.join(root, f))
  return found


def _all_tables(folder):
  found = set()
  for root, dirs, files in os.walk(folder):
    for f in files:
      if f.lower().endswith((".tsc", ".tscb")):
        found.add(os.path.join(root, f))
  return found


def _per_part_tables(folder):
  found = set()
  for root, dirs, files in os.walk(folder):
    for f in files:
      if "_part_" in f.lower() and f.lower().endswith((".tsc", ".tscb")):
        found.add(os.path.join(root, f))
  return found


def _tsc_beside(folder):
  newest, newest_t = None, -1
  for root, dirs, files in os.walk(folder):
    for f in files:
      if f.lower().endswith((".tsc", ".tscb")):
        p = os.path.join(root, f)
        t = os.path.getmtime(p)
        if t > newest_t:
          newest, newest_t = p, t
  return newest
