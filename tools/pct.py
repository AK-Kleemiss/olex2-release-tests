# pct.py <pass> [pass...]: percentage summary from <pass>.match.txt (or .interim.match.txt)
import sys, re, os, collections
for name in sys.argv[1:]:
  f = max([x for x in (name + '.match.txt', name + '.interim.match.txt') if os.path.exists(x)], key=os.path.getmtime)
  cases, cur = [], None
  for l in open(f, errors='replace'):
    if l.startswith('## '):
      cur = dict(id=l.split()[1], n=0, m=0, wrong=0, noise=0); cases.append(cur)
    elif cur is None: continue
    elif l.startswith('deposited '):
      v = l.split(); cur['n'], cur['m'], cur['right'] = float(v[1]), float(v[3]), float(v[5])
    elif l.rstrip().endswith('<-- NOISE'): cur['noise'] += 1
    elif re.search(r'<-- [A-Z][a-z]?$', l.rstrip()):
      cur['wrong'] += 1
      pair = (l.split()[1], l.split()[-1])
      cur.setdefault('pairs', []).append(pair)
  solved = [c for c in cases if c['n']]
  N = sum(c['n'] for c in solved); M = sum(c['m'] for c in solved)
  W = sum(c['m'] - c['right'] for c in solved); Z = sum(c['noise'] for c in solved)
  full = sum(1 for c in solved if c['right'] >= c['n'] - 0.5 and not c['noise'])
  typed = sum(1 for c in solved if c['right'] >= c['n'] - 0.5)
  ok_all = sum(1 for c in solved if c['m'] == c['n'])
  print("%s: %d cases, %d solved and matched (%.1f%%)" % (name, len(cases), len(solved), 100.0*len(solved)/max(1, len(cases))))
  print("  atoms: %.0f/%.0f found (%.1f%%), %.0f/%.0f correctly typed (%.1f%%), %.0f found but mistyped (%.1f%%), %d noise atoms (%.1f%% of model atoms)"
        % (M, N, 100.0*M/N, M - W, N, 100.0*(M - W)/N, W, 100.0*W/N, Z, 100.0*Z/max(1, sum(c['m'] + c['noise'] + c['wrong'] for c in solved))))
  print("  structures (of solved): every atom found %d (%.1f%%), all found and typed right %d (%.1f%%), also noise-free %d (%.1f%%)"
        % (ok_all, 100.0*ok_all/len(solved), typed, 100.0*typed/len(solved), full, 100.0*full/len(solved)))
  print("  structures (of all %d): all found and typed right %.1f%%, also noise-free %.1f%%" % (len(cases), 100.0*typed/len(cases), 100.0*full/len(cases)))
  conf = collections.Counter(p for c in solved for p in c.get('pairs', []))
  print("  confusions (ours->deposited): " + ", ".join("%s->%s %d" % (a, b, n) for (a, b), n in conf.most_common(8)))
