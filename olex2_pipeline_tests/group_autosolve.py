"""Auto-Solve on the sample structures: the multi-trial charge-flipping
pipeline with its space-group shortlist and element assignment.

Each case loads the deposited model as the reference, runs
`olex2.solve` / `Auto-Solve` through spy.RunSolutionPrg() exactly as the GUI
does, lets the deferred tidy-up (compaq, four cycles, ADP prune, re-typing)
run, and then compares what came back with what was deposited:

  sg_rank      where the deposited space group sits in the shortlist (0 = absent)
  cc           best correlation over the trials (seeds are fixed, so stable)
  peaks        peaks the solution posted, and how many of them the geometry
               classifier moved away from the density call
  allowed      the formula restriction the assignment ran under
  pruned       peaks the ADP prune removed
  retyped      labels the post-cleanup re-typing changed
  types/typed  element histogram of the final model, and the fraction of the
               deposited non-H atoms it accounts for

Auto-Solve is opt-in (user.solution.auto_solve) and needs a cctbx with
smtbx.ab_initio; without that the group skips. The geometry classifier needs
etc/geometry_aid_model.npz beside NoSpherA2.exe, and a run that falls back to
density only is a failure, not a variant.

OLEX2_TEST_AUTOSOLVE_SAMPLES  comma-separated sample names, default SAMPLES
"""
from __future__ import absolute_import, division, print_function

import os
import re

import olx
from olexFunctions import OV

from pipeline_tests import (macro, SkipTest, atom_count, space_group, has_hkl)
from group_nsa2_matrix import sample_copy, _model_file, _load_model
from group_nosphera2 import _refine_capturing

GROUP = "autosolve"
SAMPLES = ("sucrose", "epoxide", "water", "malbac")
NPZ = os.path.join("etc", "geometry_aid_model.npz")


def register(suite):
  wanted = os.environ.get("OLEX2_TEST_AUTOSOLVE_SAMPLES", "").strip()
  samples = [s.strip() for s in wanted.split(",") if s.strip()] or SAMPLES
  for s in samples:
    suite.run(GROUP, "autosolve %s" % s, t_autosolve, suite, s)


def _method():
  """The Auto-Solve method, registered for this session."""
  try:
    import smtbx.ab_initio  # noqa: F401
  except ImportError:
    raise SkipTest("this cctbx has no smtbx.ab_initio")
  import ExternalPrgParameters as EPP
  OV.SetParam('user.solution.auto_solve', True)
  EPP.SPD, EPP.RPD = EPP.defineExternalPrograms()
  prg = EPP.SPD.programs.get('olex2.solve')
  method = prg.methods.get('Auto-Solve') if prg else None
  if method is None:
    raise AssertionError("Auto-Solve is not registered with olex2.solve")
  return method


def _types():
  """element -> atoms of the model, H and Q peaks left out."""
  out = {}
  for i in range(atom_count()):
    t = str(olx.xf.au.GetAtomType(i))
    if t in ("H", "D", "Q"):
      continue
    out[t] = out.get(t, 0) + 1
  return out


def _histogram(types):
  return ",".join("%s:%d" % (k, types[k]) for k in sorted(types))


def _grab(pattern, text, default=None, cast=int):
  m = re.search(pattern, text)
  return cast(m.group(1)) if m else default


def t_autosolve(suite, sample):
  method = _method()
  if not os.path.isfile(os.path.join(OV.BaseDir(), NPZ)):
    raise AssertionError("%s is not in the run directory" % NPZ)
  folder = sample_copy(suite, sample)
  if not has_hkl(folder):
    raise SkipTest("no hkl with the %s sample" % sample)
  _load_model(folder, _model_file(folder))
  ref_types, ref_sg = _types(), space_group()
  n_ref = sum(ref_types.values())

  OV.SetParam('snum.solution.program', 'olex2.solve')
  OV.SetParam('snum.solution.method', 'Auto-Solve')
  OV.SetParam('snum.solution.retype_after_tidy', True)
  # The GUI asks "solve this again?" once a model is loaded; that would block
  # a headless run in the message box.
  ask = OV.GetParam('user.alert_solve_anyway')
  OV.SetParam('user.alert_solve_anyway', 'N')
  try:
    text = _refine_capturing("spy.RunSolutionPrg()")
  finally:
    OV.SetParam('user.alert_solve_anyway', ask)

  if "No solution found" in text or atom_count() == 0:
    raise AssertionError("Auto-Solve found no solution for %s" % sample)
  trials = _grab(r"Best of (\d+) trial", text)
  if not trials or trials < 2:
    raise AssertionError("expected several trials, log says %r" % trials)
  cc = _grab(r"Best of \d+ trial\(s\): correlation ([0-9.]+)", text, cast=float)
  if "using density only" in text or "Geometry step unavailable" in text:
    raise AssertionError("element assignment fell back to density only")
  peaks = _grab(r"Element assignment: (\d+) peaks", text)
  if peaks is None:
    raise AssertionError("no element assignment ran")
  changed = _grab(r"Element assignment: \d+ peaks, (\d+) where geometry", text, 0)
  allowed = _grab(r"Element assignment restricted to: ([^\n]+)", text, "any",
                  cast=lambda s: ",".join(s.replace(",", " ").split()))
  pruned = _grab(r"Pruned (\d+) peak", text, 0)
  if "Re-typed" in text:
    retyped = _grab(r"(\d+) label\(s\) changed", text, 0)
  elif "Re-typing skipped" in text:
    retyped = "skipped"
  else:
    raise AssertionError("the post-cleanup re-typing did not run")

  solver = method.cctbx_solver
  sugg = getattr(solver, 'solution_suggestions', None)
  entries = list(getattr(sugg, 'suggestions', []) or [])
  if not entries:
    raise AssertionError("no space-group suggestions")
  ref_no = solver.solution_f_obs.space_group().type().number()
  numbers = [e.space_group_info.type().number() for e in entries]
  sg_rank = numbers.index(ref_no) + 1 if ref_no in numbers else 0
  if not sg_rank:
    raise AssertionError("deposited %s (No. %d) is not among the suggestions %s"
                         % (ref_sg, ref_no, numbers))

  got = _types()
  typed = sum(min(got.get(k, 0), ref_types[k]) for k in ref_types)/float(n_ref)
  if len(ref_types) > 1 and list(got) == ["C"]:
    raise AssertionError("every peak stayed carbon; the formula is %s"
                         % _histogram(ref_types))
  if typed < 0.5:
    raise AssertionError("final model %s accounts for %.2f of the deposited %s"
                         % (_histogram(got), typed, _histogram(ref_types)))
  return " ".join("%s=%s" % kv for kv in [
    ("sg", ref_sg), ("sg_rank", sg_rank), ("n_suggest", len(entries)),
    ("trials", trials), ("cc", "%.3f" % cc), ("peaks", peaks),
    ("geometry_changed", changed), ("allowed", allowed), ("pruned", pruned),
    ("retyped", retyped), ("atoms", sum(got.values())),
    ("types", _histogram(got)), ("typed", "%.2f" % typed)])
