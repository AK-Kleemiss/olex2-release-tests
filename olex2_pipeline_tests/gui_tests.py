"""The GUI layer of the release gate: drive the panels the way a user does and
keep a picture of each.

Two things come out of every case here. The golden line holds what a script
can decide - the controls a panel must render, the number of sources the
NoSpherA2 combo offers, and that Olex2 logged no error while the panel was
built. The rest is for eyes: each state is written to OLEX2_TEST_GUI_DIR as a
screenshot of the Olex2 window (<state>.png) next to the html tree the panel
was rendered from (<state>.html, html.Dump), and the driver lists them in
index.html. A clipped label, a control drawn over another or a panel that
comes up empty is visible there and nowhere else.

Without a GUI (olex2c) every case skips.
"""
from __future__ import absolute_import, division, print_function

import os
import re
import struct
import subprocess
import sys
import zlib

import olx
from olexFunctions import OV

from pipeline_tests import SkipTest

# what the NoSpherA2 block of the refine panel must render once it is on with
# pTB (always in the bundle) as the source: the extras block is built per
# source by hybrid_GUI.make_*_GUI, and these are the names it emits
REFINE_SOURCE = "pTB"
REFINE_CONTROLS = ["USE_ASPHERICAL", "SNUM_REFINEMENT_NSFF_SOURCE", "reset-wfn",
                   "NoSpherA2_cpus", "NosPherA2_Memory", "NoSpherA2_purification",
                   "NoSpherA2_Charge", "NoSpherA2_Multiplicity"]

# lines Olex2 writes while a panel is built that count as a defect; the
# second pattern lists what merely contains the word
ERROR_RE = re.compile(r"(?i)\b(traceback|error|exception|unknown macro|unknown function|"
                      r"undefined|wrong number of arguments|failed|not found|could not locate)\b")
BENIGN_RE = re.compile(r"(?i)DIIS Error|Errors?\s*:|error_list|no errors|error model|"
                       r"error contribution|using built-in defaults")

# the main tabs and the tool-bar panels, each an html.ItemState from the GUI's
# own blocks (gui/blocks/tab-off.htm, gui/blocks/cbtn-on.htm)
TABS = ["home", "work", "view", "tools", "info"]
PANELS = ["solve", "refine", "draw", "report"]


def gui_dir():
  d = os.environ.get("OLEX2_TEST_GUI_DIR")
  if not d:
    log = os.environ.get("OLEX2_TEST_RELEASE_LOG")
    d = os.path.join(os.path.dirname(log) if log else OV.DataDir(), "gui")
  if not os.path.isdir(d):
    os.makedirs(d)
  return d


def _need_gui(suite=None, sample="epoxide"):
  """Skip without a GUI; with one, a structure on screen (the tabs only
  exist for a loaded file) and the log flushed line by line so a case's
  window of it is what Olex2 wrote during that case."""
  if not OV.HasGUI():
    raise SkipTest("no GUI (console build)")
  olx.app.AutoFlushLog("true")
  print("gui tests: log flushed")  # pushes the buffered start-up lines out now
  if suite is not None and not os.path.isfile(olx.FileFull() or ""):
    from group_nsa2_matrix import sample_copy, _model_file, _load_model
    folder = sample_copy(suite, sample)
    _load_model(folder, _model_file(folder))


# ---------------------------------------------------------------------------
# what Olex2 logged meanwhile
# ---------------------------------------------------------------------------

def _log_file():
  return olx.app.GetLogName()


def log_errors(start=0):
  """Error lines the Olex2 log gained since byte offset `start`."""
  path = _log_file()
  if not path or not os.path.isfile(path):
    return []
  with open(path, "rb") as f:
    f.seek(start)
    text = f.read().decode("utf-8", "replace")
  return [l.strip() for l in text.splitlines()
          if ERROR_RE.search(l) and not BENIGN_RE.search(l)]


def log_size():
  path = _log_file()
  return os.path.getsize(path) if path and os.path.isfile(path) else 0


def _assert_clean(errors, what):
  if errors:
    raise AssertionError("%s logged %d error line(s): %s" % (
      what, len(errors), " | ".join(errors[:4])))


# ---------------------------------------------------------------------------
# a picture and the html tree of the current state
# ---------------------------------------------------------------------------

def _write_png(path, w, h, bgra):
  def chunk(tag, data):
    return (struct.pack(">I", len(data)) + tag + data +
            struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff))
  raw = bytearray()
  row = bytearray(w * 3)
  for y in range(h):
    src = bgra[y * w * 4:(y + 1) * w * 4]
    row[0::3], row[1::3], row[2::3] = src[2::4], src[1::4], src[0::4]
    raw.append(0)
    raw += row
  with open(path, "wb") as f:
    f.write(b"\x89PNG\r\n\x1a\n")
    f.write(chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)))
    f.write(chunk(b"IDAT", zlib.compress(bytes(raw), 6)))
    f.write(chunk(b"IEND", b""))


def _shot_windows(path):
  import ctypes
  from ctypes import wintypes
  user32, gdi32 = ctypes.windll.user32, ctypes.windll.gdi32
  pid, best = os.getpid(), [None, 0]

  @ctypes.WINFUNCTYPE(ctypes.c_bool, wintypes.HWND, wintypes.LPARAM)
  def visit(hwnd, lp):
    owner = wintypes.DWORD()
    user32.GetWindowThreadProcessId(hwnd, ctypes.byref(owner))
    if owner.value == pid and user32.IsWindowVisible(hwnd):
      r = wintypes.RECT()
      user32.GetWindowRect(hwnd, ctypes.byref(r))
      area = (r.right - r.left) * (r.bottom - r.top)
      if area > best[1]:
        best[:] = [hwnd, area]
    return True
  user32.EnumWindows(visit, 0)
  hwnd = best[0]
  if not hwnd:
    return False
  r = wintypes.RECT()
  user32.GetWindowRect(hwnd, ctypes.byref(r))
  w, h = r.right - r.left, r.bottom - r.top

  class BITMAPINFOHEADER(ctypes.Structure):
    _fields_ = [("biSize", wintypes.DWORD), ("biWidth", ctypes.c_long),
                ("biHeight", ctypes.c_long), ("biPlanes", wintypes.WORD),
                ("biBitCount", wintypes.WORD), ("biCompression", wintypes.DWORD),
                ("biSizeImage", wintypes.DWORD), ("biXPelsPerMeter", ctypes.c_long),
                ("biYPelsPerMeter", ctypes.c_long), ("biClrUsed", wintypes.DWORD),
                ("biClrImportant", wintypes.DWORD)]
  hdc = user32.GetWindowDC(hwnd)
  mdc = gdi32.CreateCompatibleDC(hdc)
  bmp = gdi32.CreateCompatibleBitmap(hdc, w, h)
  gdi32.SelectObject(mdc, bmp)
  # PW_RENDERFULLCONTENT: the GL canvas and the html window paint too
  ok = user32.PrintWindow(hwnd, mdc, 2)
  bi = BITMAPINFOHEADER(biSize=ctypes.sizeof(BITMAPINFOHEADER), biWidth=w,
                        biHeight=-h, biPlanes=1, biBitCount=32)
  buf = ctypes.create_string_buffer(w * h * 4)
  gdi32.GetDIBits(mdc, bmp, 0, h, buf, ctypes.byref(bi), 0)
  gdi32.DeleteObject(bmp)
  gdi32.DeleteDC(mdc)
  user32.ReleaseDC(hwnd, hdc)
  if not ok:
    return False
  _write_png(path, w, h, buf.raw)
  return True


def _shot_posix(path):
  # the whole screen: under xvfb-run or on a test desktop that is Olex2
  for argv in (["screencapture", "-x", path], ["import", "-window", "root", path]):
    try:
      if subprocess.call(argv, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
        return True
    except OSError:
      pass
  return False


def snapshot(state):
  """Screenshot + html tree of the current GUI state into the gui dir."""
  d = gui_dir()
  # the page is rebuilt synchronously but painted by the event loop, and
  # PrintWindow hands back what the compositor last saw: pump, then a frame
  import time
  for _ in range(3):
    olx.Refresh()
    time.sleep(0.05)
  olx.html.Dump(os.path.join(d, state + ".html"))
  png = os.path.join(d, state + ".png")
  try:
    ok = _shot_windows(png) if sys.platform == "win32" else _shot_posix(png)
  except Exception as e:
    print("screenshot %s: %s" % (state, e))
    ok = False
  if not ok:
    print("no screenshot for %s" % state)
  with open(os.path.join(d, state + ".html"), "rb") as f:
    return f.read().decode("utf-8", "replace")


def _controls_in(dump, names):
  missing = [n for n in names
             if not re.search(r"name\s*=\s*['\"]?%s\b" % re.escape(n), dump)]
  return missing


# ---------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------

def show_panel(panel):
  olx.html.ItemState("cbtn*", "1", "cbtn-" + panel, "2", "*settings", "0",
                     panel + "-settings", "1")


def show_tab(tab):
  olx.html.ItemState("*", "0", "tab*", "2", "tab-" + tab, "1", "logo1", "1",
                     "index-" + tab + "*", "1", "info-title", "1")


def c_gui_tabs_and_panels(suite):
  """Every main tab and every tool-bar panel renders without an error line;
  a picture of each."""
  _need_gui(suite)
  seen = []
  for tab in TABS:
    at = log_size()
    show_tab(tab)
    OV.UpdateHtml()
    snapshot("tab_" + tab)
    _assert_clean(log_errors(at), "tab " + tab)
    seen.append(tab)
  show_tab("work")
  for panel in PANELS:
    at = log_size()
    show_panel(panel)
    OV.UpdateHtml()
    snapshot("panel_" + panel)
    _assert_clean(log_errors(at), "panel " + panel)
    seen.append(panel)
  return "rendered=%s errors=0" % ",".join(seen)


def c_gui_refine_nosphera2(suite):
  """The refine panel with NoSpherA2 on and its three blocks (options,
  properties, XCW) expanded holds the controls a user needs."""
  _need_gui(suite)
  from variableFunctions import nsa2_set_param
  from NoSpherA2.NoSpherA2 import change_tsc_generator
  at = log_size()
  nsa2_set_param('use_aspherical', True)
  OV.SetParam('user.NoSpherA2.show_XCW', True)
  show_tab("work")
  show_panel("refine")
  change_tsc_generator(REFINE_SOURCE)  # what the combo's onchange does
  olx.html.ItemState("h3-NoSpherA2-extras", "1", "h3-NoSpherA2-Properties", "1",
                     "h3-NoSpherA2-XCW", "1")
  OV.UpdateHtml()
  dump = snapshot("refine_nosphera2")
  missing = _controls_in(dump, REFINE_CONTROLS)
  if missing:
    raise AssertionError("refine panel lacks %s" % ", ".join(missing))
  for block in ("h3-NoSpherA2-extras", "h3-NoSpherA2-Properties", "h3-NoSpherA2-XCW"):
    if block not in dump:
      raise AssertionError("block %s is not in the rendered panel" % block)
  _assert_clean(log_errors(at), "refine panel")
  return "controls=%d/%d blocks=3 errors=0" % (
    len(REFINE_CONTROLS) - len(missing), len(REFINE_CONTROLS))


def _assert_choices_offered(src):
  """The method and basis the panel shows must be entries of the new
  source's lists, or their combos render blank."""
  from NoSpherA2.NoSpherA2 import NoSpherA2_instance as nsp2, get_functional_list
  from NoSpherA2.utilities import source_is_tsc
  from variableFunctions import nsa2_get_param
  if source_is_tsc() or src == "Thakkar IAM":
    return
  for param, offered in (('method', get_functional_list(src)),
                         ('basis_name', nsp2.getBasisListStr())):
    have = nsa2_get_param(param)
    if have not in [c.strip() for c in (offered or "").split(";")]:
      raise AssertionError("%s shows %s=%r, not in its list %s" % (src, param, have, offered))


def c_gui_sources(suite):
  """Choosing each wavefunction source in the combo rebuilds the panel, as
  the onchange of SNUM_REFINEMENT_NSFF_SOURCE does, without an error line;
  a picture per source."""
  _need_gui(suite)
  from NoSpherA2.NoSpherA2 import NoSpherA2_instance as nsp2, change_tsc_generator
  from variableFunctions import nsa2_get_param, nsa2_set_param
  show_tab("work")
  show_panel("refine")
  olx.html.ItemState("h3-NoSpherA2-extras", "1")
  entries = [s for s in nsp2.getwfn_softwares().split(";") if s.strip()]
  sources = [s for s in entries
             if not s.strip().startswith("Get ") and " -- " not in s]
  if not sources:
    raise AssertionError("the source combo offers nothing (%s)" % entries)
  old = nsa2_get_param('source')
  done, dumps = [], set()
  try:
    for src in sources:
      at = log_size()
      change_tsc_generator(src)
      # its "html.itemstate h3-NoSpherA2-extras 2 1" toggles the block, so
      # open it again for the picture
      olx.html.ItemState("h3-NoSpherA2-extras", "1")
      OV.UpdateHtml()
      slug = re.sub(r"[^A-Za-z0-9]+", "_", src.strip()).strip("_")
      dump = snapshot("source_" + slug)
      if nsa2_get_param('source').strip() != src.strip():
        raise AssertionError("choosing %r left the source at %r"
                             % (src.strip(), nsa2_get_param('source')))
      _assert_clean(log_errors(at), "source " + src.strip())
      _assert_choices_offered(src.strip())
      dumps.add(hash(dump))
      done.append(src.strip())
  finally:
    nsa2_set_param('source', old)
    OV.UpdateHtml()
  return "sources=%s distinct_panels=%d errors=0" % (",".join(done), len(dumps))


def c_gui_dispradial(suite):
  """The DispRadial tool panel, opened from the Tools tab. disp_radial is
  imported at start-up only when user.refinement.dispradial is set, so the
  case cannot switch it on for itself."""
  _need_gui(suite)
  if not OV.GetParam('user.refinement.dispradial', False):
    raise SkipTest("user.refinement.dispradial is off")
  at = log_size()
  show_tab("tools")
  olx.html.ItemState("h2-tools*", "2", "h2-tools-DispRadial-DispRadial", "1")
  OV.UpdateHtml()
  dump = snapshot("tool_dispradial")
  if "DispRadial.htm" not in dump:
    raise AssertionError("the Tools tab has no DispRadial panel")
  _assert_clean(log_errors(at), "DispRadial panel")
  return "panel=DispRadial errors=0"


def c_gui_solve_flint(suite):
  """The solve panel with FLINT chosen, as the method combo's onchange does,
  and its extra settings open."""
  _need_gui(suite)
  if not OV.GetParam('user.solution.flint', False):
    raise SkipTest("user.solution.flint is off")
  old = (OV.GetParam('snum.solution.program'), OV.GetParam('snum.solution.method'))
  at = log_size()
  try:
    show_tab("work")
    show_panel("solve")
    OV.set_solution_program("olex2.solve", "FLINT")
    olx.html.ItemState("solution-settings-extra", "1")
    OV.UpdateHtml()
    snapshot("solve_flint")
    if OV.GetParam('snum.solution.method') != "FLINT":
      raise AssertionError("choosing FLINT left the method at %r"
                           % OV.GetParam('snum.solution.method'))
    _assert_clean(log_errors(at), "solve panel with FLINT")
  finally:
    OV.set_solution_program(*old)
    OV.UpdateHtml()
  return "method=FLINT errors=0"


def c_gui_log_clean(suite):
  """What Olex2 logged since it started - the GUI build-up and the
  refinements included - that reads as an error, as distinct kinds (numbers
  collapsed). The golden holds the known ones, so a new kind is a diff;
  the goal is log_errors=0."""
  _need_gui()
  kinds = sorted(set(re.sub(r"\d+", "<n>", e)[:70] for e in log_errors(0)))
  return "log_errors=%d %s" % (len(kinds), " | ".join(kinds))


CASES = [
  ("gui_tabs_and_panels",   c_gui_tabs_and_panels),
  ("gui_refine_nosphera2",  c_gui_refine_nosphera2),
  ("gui_sources",           c_gui_sources),
  ("gui_dispradial",        c_gui_dispradial),
  ("gui_solve_flint",       c_gui_solve_flint),
  ("gui_log_clean",         c_gui_log_clean),
]
