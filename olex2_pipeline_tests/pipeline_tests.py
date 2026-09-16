"""End-to-end tests that run inside Olex2, with nothing stubbed.

The suite in the parent directory uses dummy olx/olex modules, so it can only
reach code that does not touch a structure. Everything that matters about
Olex2 - loading a model, running a macro, solving, refining, generating form
factors - needs the real thing behind it. This runs inside olex2c, which is
Olex2 without a window, so `olx`, `spy` and the macro layer are genuine.

  <olex2c> -b <rundir> pipeline.olx

where pipeline.olx is a single line:

  py.Run '<this directory>/pipeline_tests.py'

Selecting what to run, by environment variable, since py.Run takes no
arguments of its own:

  OLEX2_TEST_GROUPS   comma separated: macros, solve, refine, nosphera2
                      default: macros,solve,refine
  OLEX2_TEST_FULL     1 to include the slow quantum backends
  OLEX2_TEST_OUT      where to write the result table

Every case reports PASS, FAIL or SKIP with a reason. A backend that is not
installed is a SKIP and says so - it is never silently dropped, because a
suite that quietly tests nothing looks exactly like a suite that passes.
"""
from __future__ import absolute_import, division, print_function

import os
import shutil
import sys
import tempfile
import time
import traceback

import olx
import olex
from olexFunctions import OV


# ---------------------------------------------------------------------------
# the harness
# ---------------------------------------------------------------------------

class Result(object):
  PASS, FAIL, SKIP = "PASS", "FAIL", "SKIP"

  def __init__(self, group, name, state, detail="", seconds=0.0):
    self.group, self.name = group, name
    self.state, self.detail, self.seconds = state, detail, seconds

  def __str__(self):
    return "%-5s %-10s %-34s %6.1fs  %s" % (
      self.state, self.group, self.name, self.seconds, self.detail)

  def log_line(self):
    """The golden form: no group, no timing, single-spaced detail."""
    return "%s %s  %s" % (self.state, self.name, " ".join(str(self.detail).split()))


class Suite(object):
  def __init__(self):
    self.results = []
    self._copy_n = 0
    self.scratch = os.environ.get("OLEX2_TEST_SCRATCH") or os.path.join(
      tempfile.gettempdir(), "olex2_pipeline_tests")
    # rmtree is best-effort: a QM job that is still finishing, or a file left
    # open by a previous run, keeps its directory. Tolerate what survives
    # rather than failing the whole suite before the first case - sample()
    # takes a fresh copy per call anyway.
    if os.path.isdir(self.scratch):
      shutil.rmtree(self.scratch, ignore_errors=True)
    if not os.path.isdir(self.scratch):
      os.makedirs(self.scratch)

  def run(self, group, name, fn, *args, **kwds):
    """Run one case. A raised SkipTest is a skip, any other exception a fail."""
    t0 = time.time()
    try:
      detail = fn(*args, **kwds) or ""
      state = Result.PASS
    except SkipTest as e:
      state, detail = Result.SKIP, str(e)
    except Exception as e:
      state = Result.FAIL
      detail = "%s: %s" % (type(e).__name__, e)
      # the innermost frame: a FAIL line is never a golden, so it may say
      # where the exception came from
      tb = traceback.extract_tb(sys.exc_info()[2])
      if tb:
        detail += " [%s:%d]" % (os.path.basename(tb[-1][0]), tb[-1][1])
      if OV.GetParam('user.debug', False):
        traceback.print_exc()
    r = Result(group, name, state, detail, time.time() - t0)
    self.results.append(r)
    print(r)
    sys.stdout.flush()
    return r

  def sample(self, name):
    """A fresh copy of a sample structure, per call.

    Copied rather than used in place: a refinement overwrites its own res, so
    running against the samples directly compares each run with the previous
    run's output instead of with the structure as committed.

    Re-copied rather than reused, which is the same trap one level up. Keeping
    the first copy let one test inherit the next one's starting point: the
    NoSpherA2 group ran a spherical refinement that picked up the previous
    backend's tsc and its refined res, so the baseline the aspherical result
    was compared against had already had form factors applied.
    """
    src = os.path.join(OV.BaseDir(), "samples", name)
    if not os.path.isdir(src):
      raise SkipTest("no sample %r" % name)
    dst = os.path.join(self.scratch, name)
    if os.path.isdir(dst):
      shutil.rmtree(dst, ignore_errors=True)
    if os.path.isdir(dst):
      # something still holds a file - a distinct name beats a stale structure
      self._copy_n += 1
      dst = "%s_%d" % (dst, self._copy_n)
    shutil.copytree(src, dst)
    return dst

  def summary(self):
    n = {Result.PASS: 0, Result.FAIL: 0, Result.SKIP: 0}
    for r in self.results:
      n[r.state] += 1
    return n


class SkipTest(Exception):
  """Raised by a case that cannot run here - a missing backend, say."""


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def macro(cmd, must_return=False):
  """Run a macro.

  olex.m returns 0 for a macro that does not exist - but also for AddIns and
  sel, which return 0 having worked - so the return value is not a usable
  success signal on its own and is only checked when a caller asks.

  What a test has to do instead is assert on the effect: that the coordinates
  changed, that an R1 was recorded, that the atom count moved. Trusting the
  call to have worked is how spy.solve.do_solve - which does not exist - let
  three solution tests report the model they had loaded a moment earlier as a
  solution, in 0.0 seconds each.
  """
  rv = olex.m(cmd)
  if must_return and (not rv or str(rv).strip() in ("0", "false", "False")):
    raise AssertionError("Olex2 rejected the macro: %s" % cmd)
  return rv


def atom_count():
  return int(olx.xf.au.GetAtomCount())


def _is_true(v):
  return str(v).lower() in ("true", "1")


def deleted(i):
  """Whether atom i has been removed from the model.

  GetAtomCount keeps counting a deleted atom, so a test that counts by index
  sees no change after `kill`: sucrose reported 22 hydrogens before and after
  `kill $H`, and 44 after HAdd put 22 fresh ones beside the dead ones.
  """
  return _is_true(olx.xf.au.IsAtomDeleted(i))


def real_atom_count():
  """Atoms in the model: neither deleted nor Q peaks.

  A Q peak is a residual density maximum, not part of the model, and writing
  a file drops them - so the count that has to survive a round trip is this
  one, not GetAtomCount.
  """
  n = 0
  for i in range(atom_count()):
    if not deleted(i) and not _is_true(olx.xf.au.IsPeak(i)):
      n += 1
  return n


def space_group():
  # the function form, as the GUI asks for it; olx.SG() answers a status
  return olex.f("sg(%h)")


def load(path):
  """Load a structure and return how many atoms arrived."""
  macro("@reap '%s'" % path.replace("\\", "/"))
  n = atom_count()
  if n == 0:
    raise AssertionError("loaded %s and got no atoms" % os.path.basename(path))
  return n


def model_in(folder):
  """The largest res or ins in a folder - the structure, not a fragment."""
  best, best_size = None, -1
  for root, dirs, files in os.walk(folder):
    for f in files:
      if os.path.splitext(f)[1].lower() in (".res", ".ins"):
        p = os.path.join(root, f)
        if os.path.getsize(p) > best_size:
          best, best_size = p, os.path.getsize(p)
  if best is None:
    raise SkipTest("no res or ins in %s" % os.path.basename(folder))
  return best


def has_hkl(folder):
  for root, dirs, files in os.walk(folder):
    for f in files:
      if f.lower().endswith(".hkl"):
        return True
  return False


def _r1_raw():
  # The parameter first: olex2.refine stores it after every cycle and never
  # reloads the file, so olx.Ins("R1") - the REM R1 line of the res that was
  # *loaded* - still says what the sample shipped with. sucrose ships with
  # 0.0280, which is also what three cycles of Gauss-Newton produce, so the
  # stale value looked right until an aspherical run reached 0.0266 and the
  # test still read 0.0280. SHELXL reloads its res, so Ins is right for it,
  # and doHistoryCreation stores the parameter for it too.
  for key in ("snum.refinement.last_R1", "snum.refinement.R1"):
    v = OV.GetParam(key, None)
    if v not in (None, "", "n/a"):
      try:
        return float(v)
      except (TypeError, ValueError):
        continue
  if str(OV.GetParam('snum.refinement.program', '')).lower() == "olex2.refine":
    # nothing stored means it did not run; the file's line would be a lie
    return None
  v = olx.Ins("R1")
  try:
    return float(v)
  except (TypeError, ValueError):
    return None


def clear_r1():
  """Forget the recorded R1, so the next read cannot return a stale one.

  Called before every refinement. A refinement that does not run leaves the
  previous one's R1 in place, and a test that reads it afterwards reports the
  earlier result as though it were this one - which is how a SHELXL case
  passed in 0.1s with the R1 that olex2.refine had just produced.
  """
  try:
    OV.SetParam('snum.refinement.last_R1', '')
  except Exception:
    pass


def r1_of_last_refinement():
  """R1 as the refinement left it, from Olex2 rather than by parsing a log."""
  v = _r1_raw()
  if v is None:
    raise AssertionError("no R1 recorded - the refinement did not run")
  return v


def program_available(exe):
  """Whether an external program is where Olex2 would look for it."""
  if shutil.which(exe):
    return True
  p = os.path.join(OV.BaseDir(), exe + (".exe" if sys.platform.startswith("win") else ""))
  return os.path.exists(p)
