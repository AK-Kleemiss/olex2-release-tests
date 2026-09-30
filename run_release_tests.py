#!/usr/bin/env python3
"""Run the Olex2 release tests inside a shipped Olex2 and compare the log
with the golden.

    python run_release_tests.py --olex2-dir /path/to/olex2 [--full]
                                [--salted-model /path/to/Model_V6]

The bundle under test is the unpacked Olex2 (the directory holding olex2.tag,
the NoSpherA2 executable and util/pyUtil). Nothing in it is modified beyond a
runonce.release_tests.txm macro file that Olex2 executes once on start-up and
that is removed again when the run is over. All state - data directory,
refinement scratch, a copy of the sample data, the log - goes to --out-dir
(default: release_out next to this script).

Olex2 starts with its GUI (on a Linux box without a display wrap the call in
xvfb-run), runs olex2_pipeline_tests/run_pipeline.py from this repository
(the same files live in the Olex2 SVN under util/pyUtil/regression, but they
are not shipped to users) for the "release" group and exits by itself. Before
anything starts the repository is fetched from GitHub and fast-forwarded when
it is behind, so the goldens are the current ones (--no-fetch skips that).
The log Olex2 writes is compared with expected/release_tests.<tier>.good by
compare_golden.py:

    exit 0   every case matches the golden
    exit 2   every case that ran matches, but some were SKIPPED (an external
             program or the SALTED model is missing): the banner says which
    exit 1   a case failed, a number drifted, a case is missing, the run did
             not finish, or Olex2 did not start

--update-golden accepts the fresh log as the new golden instead of comparing.

Windows developers can pass --olex2c <olex2c.exe> to run the same thing
through the console build (no window, needs PYTHONHOME of the Python the
bundle was built for, or --pythonhome).
"""
import argparse
import os
import platform
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare_golden  # noqa: E402

RUNONCE = "runonce.release_tests.txm"
TESTS_DIR = os.path.join(HERE, "olex2_pipeline_tests")


def platform_tag():
    s = sys.platform
    if s.startswith("win"):
        return "win"
    if s == "darwin":
        return "mac"
    return "linux"


def die(msg, code=1):
    print("ERROR: " + msg, file=sys.stderr)
    sys.exit(code)


def banner(title, lines):
    print()
    print("#" * 72)
    print("#  " + title)
    print("#" * 72)
    for l in lines:
        print("#  " + l)
    print("#" * 72)
    print()


def find_exe(olex2_dir, explicit):
    if explicit:
        if not os.path.isfile(explicit):
            die("--olex2-exe %s does not exist" % explicit)
        return explicit
    # olex2.exe / the olex2 launcher script replace themselves with the real
    # program (execv), so waiting on them returns at once: run the real one
    if platform_tag() == "win":
        names = ["olex2.dll"]
    else:
        names = ["olex2_exe"]
    for n in names:
        p = os.path.join(olex2_dir, n)
        if os.path.isfile(p):
            return p
    die("no Olex2 program (%s) in %s; pass --olex2-exe" % ("/".join(names), olex2_dir))


def check_bundle(olex2_dir):
    problems = []
    if not os.path.isfile(os.path.join(olex2_dir, "olex2.tag")):
        problems.append("no olex2.tag - is this an unpacked Olex2?")
    nsa2 = [n for n in ("NoSpherA2.exe", "NoSpherA2") if os.path.isfile(os.path.join(olex2_dir, n))]
    if not nsa2:
        problems.append("no NoSpherA2 executable next to olex2.tag")
    if os.path.exists(os.path.join(olex2_dir, RUNONCE)):
        problems.append("%s already in the bundle: another run is in progress, or a "
                        "previous one was killed; remove it by hand" % RUNONCE)
    if problems:
        die("\n  ".join(["the bundle in %s cannot be tested:" % olex2_dir] + problems))
    tag = open(os.path.join(olex2_dir, "olex2.tag"), errors="ignore").read().strip()
    return tag, nsa2[0]


def check_tests_dir(tests_dir):
    for n in ("run_pipeline.py", "pipeline_tests.py", "group_nosphera2.py",
              "group_nsa2_matrix.py", "group_release.py"):
        if not os.path.isfile(os.path.join(tests_dir, n)):
            die("%s missing from --tests-dir %s" % (n, tests_dir))


def fetch_updates():
    """The goldens are the contract: make sure this clone is at the remote's
    head before a single case runs. A clone with local edits is left alone
    but said so, loudly."""
    def git(*a):
        return subprocess.run(["git", "-C", HERE] + list(a), capture_output=True, text=True)
    if git("rev-parse", "--is-inside-work-tree").returncode != 0:
        banner("NOT A GIT CLONE", ["%s is not a clone of the release-tests repository," % HERE,
                                   "the goldens may be stale - clone it to get update checks"])
        return
    if git("fetch", "origin").returncode != 0:
        banner("NO GITHUB", ["git fetch origin failed, cannot check the goldens for updates"])
        return
    behind = git("rev-list", "--count", "HEAD..@{u}").stdout.strip()
    if behind in ("", "0"):
        print("release-tests repository is up to date with GitHub")
        return
    if git("status", "--porcelain", "--untracked-files=no").stdout.strip():
        banner("GOLDENS OUT OF DATE", ["GitHub is %s commit(s) ahead but this clone has local changes;" % behind,
                                       "not pulling - commit or discard them, or pass --no-fetch"])
        return
    r = git("pull", "--ff-only")
    if r.returncode != 0:
        banner("GOLDENS OUT OF DATE", ["GitHub is %s commit(s) ahead and the fast-forward failed:" % behind]
               + r.stderr.strip().splitlines())
        return
    banner("UPDATED", ["pulled %s commit(s) from GitHub; the run uses the new goldens" % behind,
                       "restart if run_release_tests.py itself changed"])


def fresh_dir(path):
    if os.path.isdir(path):
        shutil.rmtree(path)
    os.makedirs(path)


def write_text(path, text):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


def prepare_out(args):
    out = os.path.abspath(args.out_dir)
    data = os.path.join(out, "data")
    scratch = os.path.join(out, "scratch")
    samples = os.path.join(out, "samples")
    for d in (data, scratch, samples):
        fresh_dir(d)
    # leftovers of an earlier run must not pass for this run's results
    for name in os.listdir(out):
        path = os.path.join(out, name)
        if os.path.isfile(path) and (name.startswith("entry.") or name.startswith("release_tests.")
                                     or name.startswith("olex2")):
            os.remove(path)
    # nothing is refined inside the repository: the cases copy from here.
    # A flint-only --cases run stages just its own samples: the whole
    # set is 1 GB per out-dir, and a sweep makes hundreds of out-dirs
    only = [c.strip() for c in (args.cases or "").split(",") if c.strip()]
    if only and all(c.startswith("flint_") for c in only):
        only = set(c[len("flint_"):] for c in only)
    else:
        only = None
    for name in sorted(os.listdir(args.data_dir)):
        src = os.path.join(args.data_dir, name)
        if os.path.isdir(src) and (only is None or name.lower() in only):
            shutil.copytree(src, os.path.join(samples, name))
    # a fresh data directory would ask things on start-up and on exit
    write_text(os.path.join(data, ".options"), "confirm_on_close=false\n")
    write_text(os.path.join(data, "usettings.dat"),
               "proxy=\nrepository=\nupdate=Never\nlastupdate=%d\nexceptions=\n"
               "dest_repository=\nsrc_for_dest=\nskip=\nolex-port=\nask_update=false\n"
               % int(time.time()))
    return out, data, scratch, samples


ENTRY = r'''# written by run_release_tests.py; Olex2 runs it once through py.Run
import os, sys, traceback
out = {out!r}
open(os.path.join(out, "entry.started"), "w").write("started\n")
try:
  script = {script!r}
  g = {{"__name__": "olex2_release_tests_entry", "__file__": script}}
  exec(compile(open(script).read(), script, "exec"), g)
except SystemExit:
  pass
except Exception:
  with open(os.path.join(out, "olex2_python_error.txt"), "w") as f:
    traceback.print_exc(file=f)
  traceback.print_exc()
finally:
  open(os.path.join(out, "entry.finished"), "w").write("finished\n")
'''


def write_entry(out, tests_dir):
    path = os.path.join(out, "data", "release_tests_entry.py")
    write_text(path, ENTRY.format(out=out, script=os.path.join(tests_dir, "run_pipeline.py")))
    return path


def build_env(args, olex2_dir, out, data, scratch, samples, log_path, table_path):
    env = dict(os.environ)
    tag = platform_tag()
    env["OLEX2_DIR"] = olex2_dir
    # a 64-bit Debug build of the GUI reads OLEX2_DEBUG_DIR instead
    # (olex/xglapp.cpp, TGlXApp::OnInit); harmless for release builds
    env["OLEX2_DEBUG_DIR"] = olex2_dir
    env["OLEX2_DATADIR"] = data
    env["OLEX2_DATADIR_STATIC"] = "TRUE"
    env.pop("OLEX2_CCTBX_DIR", None)
    env["PATH"] = olex2_dir + os.pathsep + env.get("PATH", "")
    if args.pythonhome:
        env["PYTHONHOME"] = args.pythonhome
    elif tag == "win":
        bundled = os.path.join(olex2_dir, "Python")
        if os.path.isdir(bundled):
            env["PYTHONHOME"] = bundled
    else:
        # what the olex2 launcher script sets before it execs olex2_exe
        env["PYTHONHOME"] = olex2_dir
        env.setdefault("BOOST_ADAPTBX_FPE_DEFAULT", "1")
        env.setdefault("BOOST_ADAPTBX_SIGNALS_DEFAULT", "1")
        env.setdefault("OLEX2_GL_STEREO", "FALSE")
        ld = "DYLD_LIBRARY_PATH" if tag == "mac" else "LD_LIBRARY_PATH"
        env[ld] = os.pathsep.join([os.path.join(olex2_dir, "lib"), os.path.join(olex2_dir, "ilib"),
                                   os.path.join(olex2_dir, "cctbx", "cctbx_build", "lib")]
                                  + ([env[ld]] if env.get(ld) else []))
    if args.olex2c:
        # the console build gets no start-up python of its own
        bld = os.path.join(olex2_dir, "cctbx", "cctbx_build")
        env["LIBTBX_BUILD"] = bld
        pp = [os.path.join(olex2_dir, "cctbx", "cctbx_sources"),
              os.path.join(olex2_dir, "cctbx", "cctbx_sources", "boost_adaptbx"),
              os.path.join(bld, "lib")]
        if env.get("PYTHONHOME"):
            pp.append(os.path.join(env["PYTHONHOME"], "Lib", "site-packages"))
        env["PYTHONPATH"] = os.pathsep.join(pp)
    env["OLEX2_TEST_GROUPS"] = "release"
    env["OLEX2_TEST_FULL"] = "1" if args.full else ""
    env["OLEX2_TEST_CASES"] = args.cases or ""
    env["OLEX2_TEST_NCPUS"] = str(args.ncpus)
    env["OLEX2_TEST_SALTED_MODEL"] = os.path.abspath(args.salted_model) if args.salted_model else ""
    env["OLEX2_TEST_SAMPLE_DIR"] = samples
    env["OLEX2_TEST_SCRATCH"] = scratch
    env["OLEX2_TEST_OUT"] = table_path
    env["OLEX2_TEST_RELEASE_LOG"] = log_path
    env["OLEX2_TEST_GUI_DIR"] = os.path.join(out, "gui")
    return env


def kill_tree(proc):
    if proc.poll() is not None:
        return
    try:
        if platform_tag() == "win":
            subprocess.call(["taskkill", "/T", "/F", "/PID", str(proc.pid)],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            import signal
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
    except Exception as e:
        print("could not kill the Olex2 tree: %s" % e)
    try:
        proc.wait(30)
    except Exception:
        pass


def wait_for(proc, out, timeout, start_timeout, gui):
    """True when Olex2 left on its own. Kills it when the start-up or the
    whole run overruns, or when the hook is done and exit did not work."""
    t0 = time.time()
    started = os.path.join(out, "entry.started")
    finished = os.path.join(out, "entry.finished")
    finished_at = None
    while True:
        rc = proc.poll()
        if rc is not None:
            print("Olex2 exited with %s after %.0f s" % (rc, time.time() - t0))
            return True
        now = time.time()
        if finished_at is None and os.path.exists(finished):
            finished_at = now
        if gui and not os.path.exists(started) and now - t0 > start_timeout:
            print("Olex2 did not reach the test hook within %d s (a dialog on start-up? "
                  "see %s)" % (start_timeout, os.path.join(out, "olex2_stdout.txt")))
            kill_tree(proc)
            return False
        if finished_at is not None and now - finished_at > 60:
            print("the tests finished but Olex2 did not exit within 60 s; killing it")
            kill_tree(proc)
            return True
        if now - t0 > timeout:
            print("timeout after %d s; killing Olex2" % timeout)
            kill_tree(proc)
            return False
        time.sleep(1)


def launch(args, olex2_dir, exe, env, out, entry):
    stdout = open(os.path.join(out, "olex2_stdout.txt"), "wb")
    kw = dict(cwd=olex2_dir, env=env, stdout=stdout, stderr=subprocess.STDOUT)
    if platform_tag() != "win":
        kw["start_new_session"] = True
    if args.olex2c:
        cmds = os.path.join(out, "olex2c_cmds.txt")
        write_text(cmds, "@py \"exec(open(r'%s').read())\"\n" % entry)
        argv = [args.olex2c, "script", cmds]
    else:
        argv = [exe]
    print("running: " + " ".join(argv))
    return subprocess.Popen(argv, **kw)


def write_gui_index(out):
    """gui/index.html: every screenshot the GUI cases took, for a look."""
    d = os.path.join(out, "gui")
    if not os.path.isdir(d):
        return None
    shots = sorted(f for f in os.listdir(d) if f.endswith(".png"))
    rows = []
    for s in shots:
        rows.append("<h2>%s</h2><img src='%s' style='max-width:100%%;border:1px solid #888'>"
                    % (s[:-4], s))
        if os.path.isfile(os.path.join(d, s[:-4] + ".html")):  # the render cases dump no panel
            rows.append("<p><a href='%s.html'>html tree</a></p>" % s[:-4])
    write_text(os.path.join(d, "index.html"),
               "<html><body style='font-family:sans-serif'>%s</body></html>\n" % "\n".join(rows))
    return d, len(shots)


def print_timings(log_path):
    stem = os.path.splitext(log_path)[0] + "_timings.txt"
    if not os.path.isfile(stem):
        return
    rows = []
    for line in open(stem):
        parts = line.split()
        if len(parts) == 2:
            try:
                rows.append((float(parts[1]), parts[0]))
            except ValueError:
                pass
    rows.sort(reverse=True)
    total = sum(r[0] for r in rows)
    print("timings: %.0f s in all; slowest: %s" % (
        total, ", ".join("%s %.0f s" % (n, s) for s, n in rows[:5])))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("--olex2-dir", required=True, help="unpacked Olex2 to test")
    ap.add_argument("--olex2-exe", help="the program to start (default: olex2.dll on Windows, olex2_exe elsewhere)")
    ap.add_argument("--olex2c", help="Windows: run through this olex2c.exe console build instead of the GUI")
    ap.add_argument("--pythonhome", help="PYTHONHOME for the run (default: <olex2-dir>/Python when present)")
    ap.add_argument("--data-dir", default=os.path.join(HERE, "sample_data"))
    ap.add_argument("--tests-dir", default=TESTS_DIR,
                    help="the olex2_pipeline_tests modules to run (default: the copy in this repository)")
    ap.add_argument("--no-fetch", action="store_true", help="do not check GitHub for newer goldens first")
    ap.add_argument("--salted-model", help="directory holding the SALTED model (.salted); without it the SALTED cases are skipped")
    tier = ap.add_mutually_exclusive_group()
    tier.add_argument("--full", action="store_true", help="also the ORCA cases (ORCA must be installed)")
    tier.add_argument("--quick", action="store_true", help="the default tier")
    ap.add_argument("--out-dir", default=os.path.join(HERE, "release_out"))
    ap.add_argument("--golden", help="golden log (default: expected/release_tests.<tier>[.<platform>].good)")
    ap.add_argument("--update-golden", action="store_true", help="accept the log as the golden")
    ap.add_argument("--ncpus", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=3600, help="seconds for the whole run")
    ap.add_argument("--start-timeout", type=int, default=180, help="seconds for Olex2 to reach the hook")
    ap.add_argument("--cases", help="comma separated case names, default all of the tier")
    ap.add_argument("--keep-scratch", action="store_true", help="leave the refinement scratch in the out dir")
    args = ap.parse_args(argv)

    olex2_dir = os.path.abspath(args.olex2_dir)
    if not os.path.isdir(olex2_dir):
        die("--olex2-dir %s is not a directory" % olex2_dir)
    if not os.path.isdir(args.data_dir):
        die("--data-dir %s is not a directory" % args.data_dir)
    if args.olex2c and not os.path.isfile(args.olex2c):
        die("--olex2c %s does not exist" % args.olex2c)
    tier_name = "full" if args.full else "quick"
    golden = args.golden
    if not golden:
        per_platform = os.path.join(HERE, "expected", "release_tests.%s.%s.good" % (tier_name, platform_tag()))
        golden = per_platform if os.path.isfile(per_platform) else \
            os.path.join(HERE, "expected", "release_tests.%s.good" % tier_name)

    if not args.no_fetch:
        fetch_updates()
    tests_dir = os.path.abspath(args.tests_dir)
    check_tests_dir(tests_dir)
    tag, nsa2 = check_bundle(olex2_dir)
    exe = None if args.olex2c else find_exe(olex2_dir, args.olex2_exe)
    out, data, scratch, samples = prepare_out(args)
    log_path = os.path.join(out, "release_tests.%s.log" % tier_name)
    table_path = os.path.join(out, "release_tests.%s.table.txt" % tier_name)
    entry = write_entry(out, tests_dir)
    env = build_env(args, olex2_dir, out, data, scratch, samples, log_path, table_path)

    print("Olex2 %s in %s (%s), tier %s, %s" % (tag, olex2_dir, nsa2, tier_name,
                                               platform.platform()))
    print("out dir %s" % out)
    if not args.salted_model:
        print("no --salted-model: the SALTED cases will be skipped")

    runonce = os.path.join(olex2_dir, RUNONCE)
    left_on_its_own = False
    try:
        if not args.olex2c:
            write_text(runonce, "py.Run \"%s\"\nexit\n" % entry.replace("\\", "/"))
        proc = launch(args, olex2_dir, exe, env, out, entry)
        left_on_its_own = wait_for(proc, out, args.timeout, args.start_timeout, gui=not args.olex2c)
    finally:
        if os.path.exists(runonce):
            os.remove(runonce)

    err = os.path.join(out, "olex2_python_error.txt")
    if os.path.isfile(err):
        print("the test hook raised:")
        print(open(err, errors="replace").read())
    if not os.path.isfile(log_path):
        banner("NO LOG", ["Olex2 wrote no %s" % log_path,
                          "see %s" % os.path.join(out, "olex2_stdout.txt")])
        return 1
    print()
    print(open(log_path, errors="replace").read().rstrip())
    print()
    print_timings(log_path)
    gui = write_gui_index(out)
    if gui:
        print("gui    %s (%d screenshots, open index.html)" % gui)
    if not args.keep_scratch:
        shutil.rmtree(scratch, ignore_errors=True)

    only = [c.strip() for c in args.cases.split(",") if c.strip()] if args.cases else None
    if args.update_golden:
        if only:
            banner("NOT UPDATED", ["--update-golden needs the whole tier, not a --cases subset"])
            return 1
        cases, end, order = compare_golden.parse(log_path)
        if end is None:
            banner("NOT UPDATED", ["the log has no END line, the run did not complete"])
            return 1
        skips = [n for n in order if cases[n][0] == "SKIP"]
        fails = [n for n in order if cases[n][0] == "FAIL"]
        if fails:
            banner("NOT UPDATED", ["a golden may not hold failures:"] + fails)
            return 1
        os.makedirs(os.path.dirname(golden), exist_ok=True)
        compare_golden.main([golden, log_path, "--update"])
        if skips:
            banner("GOLDEN HOLDS SKIPS", ["these cases were not run and have no reference:"]
                   + skips + ["rerun with everything installed and --update-golden"])
            return 2
        return 0

    if not os.path.isfile(golden):
        banner("NO GOLDEN", ["%s does not exist" % golden,
                             "run once with --update-golden to create it"])
        return 1
    rc = compare_golden.compare(golden, log_path, only=only)
    if not left_on_its_own and rc == 0:
        rc = 1
        print("Olex2 had to be killed, treating the run as failed")
    print("result: %s" % {0: "PASS", 1: "FAIL", 2: "PASS with skips"}[rc])
    return rc


if __name__ == "__main__":
    sys.exit(main())
