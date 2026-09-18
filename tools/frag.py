import olex, olx, os
def sites():
  return [(str(olx.xf.au.GetAtomName(i)), tuple(map(float, str(olx.xf.au.GetAtomCrd(i)).split(','))))
          for i in range(int(olx.xf.au.GetAtomCount())) if olx.xf.au.IsAtomDeleted(i) != 'true']
for res in os.environ['DIAG_RES'].split(';'):
  olex.m("reap '%s'" % res)
  a = dict(sites()); olex.m("compaq -a"); b = dict(sites())
  moved = [n for n in a if n in b and max(abs(x - y) for x, y in zip(a[n], b[n])) > 0.01]
  print("###", os.path.basename(os.path.dirname(res)), "moved", len(moved), "of", len(a))
