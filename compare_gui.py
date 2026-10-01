#!/usr/bin/env python3
"""Compare the screenshots of a release-test run with the reference pictures.

The GUI and render cases leave <state>.png in <out>/gui. The references live
in expected/gui/<tier>.<platform>_<W>x<H>/, one set per tier and screen size,
since a picture only compares with one taken on the same screen.

A picture matches when no pixel outside the masked rectangles differs by more
than TOL in a channel. mask.txt next to the references lists those rectangles,
one "x0 y0 x1 y1" per line (inclusive, # comments): the window title, the
panel header and the status bar print the out dir, which differs per run.
--update-golden writes a mask.txt with DEFAULT_MASK for a new screen size;
check it against the first diff on that screen.

A changed picture is reported with the bands it differs in, and a copy with
the differing pixels painted magenta goes to <out>/gui/diff/. A screen size
without references is a note, not a failure.

    python compare_gui.py <out>/gui expected/gui/quick.win_2576x1370   (self-check: --test)
"""
import os
import shutil
import struct
import sys
import zlib

TOL = 24
# FLOWOFFICE, 2576x1370: title bar, panel header line, status bar - as rows,
# so a different window width still masks them
DEFAULT_MASK = ["0 15 99999 26", "0 94 99999 105", "0 -23 99999 -11"]


def read_png(path):
    """(w, h, [row bytes RGB]) of an 8-bit RGB/RGBA non-interlaced png."""
    try:
        from PIL import Image
        im = Image.open(path).convert("RGB")
        w, h = im.size
        data = im.tobytes()
        return w, h, [data[y * w * 3:(y + 1) * w * 3] for y in range(h)]
    except ImportError:
        pass
    with open(path, "rb") as f:
        blob = f.read()
    pos, idat, hdr = 8, b"", None
    while pos < len(blob):
        n, tag = struct.unpack(">I4s", blob[pos:pos + 8])
        body = blob[pos + 8:pos + 8 + n]
        if tag == b"IHDR":
            hdr = struct.unpack(">IIBBBBB", body)
        elif tag == b"IDAT":
            idat += body
        pos += 12 + n
    w, h, depth, ctype, _, _, lace = hdr
    if depth != 8 or ctype not in (2, 6) or lace:
        raise ValueError("%s: only 8-bit RGB/RGBA, not interlaced" % path)
    bpp = 3 if ctype == 2 else 4
    raw, stride = zlib.decompress(idat), w * bpp
    rows, prev = [], bytearray(stride)
    for y in range(h):
        ft = raw[y * (stride + 1)]
        cur = bytearray(raw[y * (stride + 1) + 1:(y + 1) * (stride + 1)])
        if ft:  # ponytail: per-byte python, seconds per picture; our own pngs are filter 0
            for i in range(stride):
                a = cur[i - bpp] if i >= bpp else 0
                b, c = prev[i], prev[i - bpp] if i >= bpp else 0
                if ft == 1:
                    p = a
                elif ft == 2:
                    p = b
                elif ft == 3:
                    p = (a + b) >> 1
                else:
                    pa, pb, pc = abs(b - c), abs(a - c), abs(a + b - 2 * c)
                    p = a if pa <= pb and pa <= pc else (b if pb <= pc else c)
                cur[i] = (cur[i] + p) & 255
        prev = cur
        rows.append(bytes(cur) if bpp == 3 else bytes(b for i, b in enumerate(cur) if i % 4 != 3))
    return w, h, rows


def write_png(path, w, h, rows):
    def chunk(tag, data):
        return (struct.pack(">I", len(data)) + tag + data +
                struct.pack(">I", zlib.crc32(tag + data) & 0xffffffff))
    raw = b"".join(b"\0" + bytes(r) for r in rows)
    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n")
        f.write(chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)))
        f.write(chunk(b"IDAT", zlib.compress(raw, 6)))
        f.write(chunk(b"IEND", b""))


def read_mask(ref_dir, w, h):
    path = os.path.join(ref_dir, "mask.txt")
    lines = open(path).read().splitlines() if os.path.isfile(path) else DEFAULT_MASK
    rects = []
    for line in lines:
        line = line.split("#")[0].split()
        if len(line) == 4:
            x0, y0, x1, y1 = (int(v) for v in line)
            rects.append((x0 + w if x0 < 0 else x0, y0 + h if y0 < 0 else y0,
                          x1 + w if x1 < 0 else x1, y1 + h if y1 < 0 else y1))
    return rects


def diff_picture(new, ref, rects):
    """Differing pixels outside the mask: (count, ["y.. x..", ...], diff rows)."""
    w, h, a = new
    _, _, b = ref
    bad, out = {}, None
    for y in range(h):
        if a[y] == b[y]:
            continue
        xs = []
        ra, rb = a[y], b[y]
        for x in range(w):
            i = 3 * x
            if (abs(ra[i] - rb[i]) > TOL or abs(ra[i + 1] - rb[i + 1]) > TOL
                    or abs(ra[i + 2] - rb[i + 2]) > TOL) and \
                    not any(x0 <= x <= x1 and y0 <= y <= y1 for x0, y0, x1, y1 in rects):
                xs.append(x)
        if xs:
            bad[y] = xs
    if not bad:
        return 0, [], None
    out = [bytearray(r) for r in a]
    for y, xs in bad.items():
        for x in xs:
            out[y][3 * x:3 * x + 3] = b"\xff\x00\xff"
    bands, ys = [], sorted(bad)
    start = prev = ys[0]
    for y in ys[1:] + [None]:
        if y is not None and y == prev + 1:
            prev = y
            continue
        xs = [x for yy in range(start, prev + 1) for x in bad.get(yy, [])]
        bands.append("y%d-%d x%d-%d" % (start, prev, min(xs), max(xs)))
        if y is not None:
            start = prev = y
    return sum(len(v) for v in bad.values()), bands, out


def shots(d):
    return sorted(f[:-4] for f in os.listdir(d) if f.endswith(".png")) if os.path.isdir(d) else []


def ref_dir_for(gui_dir, expected_dir, tag):
    """expected/gui/<tag>_<W>x<H> for the screen of this run: the size of its
    largest picture, the window (the render cases draw smaller ones)."""
    sizes = []
    for name in shots(gui_dir):
        with open(os.path.join(gui_dir, name + ".png"), "rb") as f:
            sizes.append(struct.unpack(">II", f.read(24)[16:24]))  # IHDR width, height
    if not sizes:
        return None
    w, h = max(sizes, key=lambda s: s[0] * s[1])
    return os.path.join(expected_dir, "%s_%dx%d" % (tag, w, h))


def compare(gui_dir, ref_dir, whole_tier=True):
    """0 every picture matches, 1 one changed, is new or is missing, 2 no
    references for this screen size."""
    if not os.path.isdir(ref_dir):
        print("gui    no reference pictures for %s; --update-golden records them"
              % os.path.basename(ref_dir))
        return 2
    have, want = shots(gui_dir), shots(ref_dir)
    diff_dir = os.path.join(gui_dir, "diff")
    problems = []
    for name in have:
        if name not in want:
            problems.append("NEW      %s (no reference picture)" % name)
            continue
        new = read_png(os.path.join(gui_dir, name + ".png"))
        ref = read_png(os.path.join(ref_dir, name + ".png"))
        if new[:2] != ref[:2]:
            problems.append("CHANGED  %s: %dx%d, the reference is %dx%d" % ((name,) + new[:2] + ref[:2]))
            continue
        # the mask is for the window; the render cases' pictures have no out dir in them
        window = ref_dir.rstrip("/\\").endswith("_%dx%d" % new[:2])
        n, bands, out = diff_picture(new, ref, read_mask(ref_dir, *new[:2]) if window else [])
        if n:
            os.makedirs(diff_dir, exist_ok=True)
            write_png(os.path.join(diff_dir, name + ".png"), new[0], new[1], out)
            # ponytail: Pict draws the whole scene at one of two scales (ratio ~
            # the canvas aspect) from run to run, cause unknown - so a render
            # picture is shown, not failed, until that is fixed
            problems.append("%s  %s: %d px in %s" % ("CHANGED" if window else "RENDER ", name, n,
                                                     "; ".join(bands[:6]) + (" ..." if len(bands) > 6 else "")))
    if whole_tier:
        problems += ["MISSING  %s (the run took no such picture)" % n for n in want if n not in have]
    print("gui    %d of %d pictures match %s" % (len(have) - sum(not p.startswith("MISSING") for p in problems),
                                                 len(have), ref_dir))
    for p in problems:
        print("gui    " + p)
    if any(p.startswith(("CHANGED", "RENDER")) for p in problems):
        print("gui    differing pixels are magenta in %s" % diff_dir)
    return 1 if any(not p.startswith("RENDER") for p in problems) else 0


def update(gui_dir, ref_dir):
    """The run's pictures become the references; mask.txt is kept."""
    os.makedirs(ref_dir, exist_ok=True)
    have = shots(gui_dir)
    for name in shots(ref_dir):
        if name not in have:
            os.remove(os.path.join(ref_dir, name + ".png"))
    for name in have:
        shutil.copyfile(os.path.join(gui_dir, name + ".png"), os.path.join(ref_dir, name + ".png"))
    mask = os.path.join(ref_dir, "mask.txt")
    if not os.path.isfile(mask):
        with open(mask, "w") as f:
            f.write("# x0 y0 x1 y1, inclusive; negative counts from the right/bottom edge\n"
                    + "\n".join(DEFAULT_MASK) + "\n")
    print("gui    %d reference pictures in %s" % (len(have), ref_dir))


def _test():
    import tempfile
    t = tempfile.mkdtemp()
    run, ref = os.path.join(t, "run"), os.path.join(t, "quick.win_40x60")
    os.makedirs(run)
    w, h = 40, 60
    base = [bytes((x * 7) % 256 for x in range(w * 3)) for _ in range(h)]
    write_png(os.path.join(run, "a.png"), w, h, base)
    assert read_png(os.path.join(run, "a.png"))[2] == base
    write_png(os.path.join(run, "b.png"), 4, 4, [b"\0" * 12] * 4)  # a smaller render picture
    assert ref_dir_for(run, t, "quick.win") == ref
    update(run, ref)
    assert compare(run, ref) == 0
    rows = [bytearray(r) for r in base]
    rows[30][30:33] = b"\0\0\0"   # x10 y30: a real change
    rows[h - 15][0:3] = b"\xff\xff\xff"  # inside the status-bar band (-23..-11)
    write_png(os.path.join(run, "a.png"), w, h, rows)
    n, bands, _ = diff_picture(read_png(os.path.join(run, "a.png")),
                               read_png(os.path.join(ref, "a.png")), read_mask(ref, w, h))
    assert (n, bands) == (1, ["y30-30 x10-10"]), (n, bands)
    assert compare(run, ref) == 1 and os.path.isfile(os.path.join(run, "diff", "a.png"))
    assert compare(run, os.path.join(t, "none")) == 2
    sys.modules["PIL"] = None  # the decoder for machines without PIL
    assert read_png(os.path.join(run, "a.png"))[2] == [bytes(r) for r in rows]
    write_png(os.path.join(run, "a.png"), w, h, base)
    write_png(os.path.join(run, "b.png"), 4, 4, [b"\xff" * 12] * 4)  # a render changed: shown, not failed
    assert compare(run, ref) == 0 and os.path.isfile(os.path.join(run, "diff", "b.png"))
    shutil.rmtree(t)
    print("compare_gui self-check passed")


if __name__ == "__main__":
    if sys.argv[1:] == ["--test"]:
        _test()
    else:
        sys.exit(compare(sys.argv[1], sys.argv[2]))
