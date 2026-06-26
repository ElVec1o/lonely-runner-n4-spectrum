#!/usr/bin/env python3
"""
Diagnostic on the saturated bands where the active wall at the global maximizer is NOT
primitive. Question: does the (looser) Lean `PrimitiveActiveWall` witness still exist?

The Lean def needs, for saturated M with mD<1/4:  SOME τ with mgap τ > 1/4 (a *lonely point*,
not necessarily the maximizer) and SOME pair (i,j,ε) (not necessarily active) with integer wall
functional and gcd(A,B)=1.

For each flagged band we report, AT EVERY GLOBAL MAXIMIZER VERTEX:
  - the active runner set,
  - ALL pairs (i,j,ε) whose wall functional is an integer there (active OR not), with gcd(A,B),
  - whether a PRIMITIVE one exists at that vertex.
Then a verdict: is a primitive integer wall available at some maximizer (so the claim holds with
τ = a maximizer), or must one leave the maximizer set entirely (or is it a true refutation)?
"""
from math import gcd
from itertools import combinations
from fractions import Fraction as F


def jk(a, b, c, d, e, f):
    return [(a, d), (a, -d), (b, e), (c, f)]


def minor(rows, p, q):
    return rows[p][0] * rows[q][1] - rows[p][1] * rows[q][0]


def sat_gcd(rows):
    g = 0
    for p, q in combinations(range(4), 2):
        g = gcd(g, minor(rows, p, q))
    return g


def nf(x):
    fr = x - (x.numerator // x.denominator)
    return min(fr, 1 - fr)


def all_maximizers(rows):
    walls = []
    for i, j in combinations(range(4), 2):
        for eps in (1, -1):
            w = (rows[i][0] - eps * rows[j][0], rows[i][1] - eps * rows[j][1])
            if w == (0, 0):
                continue
            lo = min(0, w[0]) + min(0, w[1])
            hi = max(0, w[0]) + max(0, w[1])
            walls.append((w, lo, hi))
    best = F(0)
    verts = set()
    for a in range(len(walls)):
        wa, loa, hia = walls[a]
        for b in range(a + 1, len(walls)):
            wb, lob, hib = walls[b]
            det = wa[0] * wb[1] - wa[1] * wb[0]
            if det == 0:
                continue
            for ma in range(loa, hia + 1):
                for mb in range(lob, hib + 1):
                    t0 = F(ma * wb[1] - wa[1] * mb, det)
                    t1 = F(wa[0] * mb - ma * wb[0], det)
                    if 0 <= t0 < 1 and 0 <= t1 < 1:
                        g = min(nf(F(r[0]) * t0 + F(r[1]) * t1) for r in rows)
                        if g > best:
                            best = g
                            verts = {(t0, t1)}
                        elif g == best:
                            verts.add((t0, t1))
    return best, sorted(verts)


def integer_walls_at(rows, tau):
    """All (i,j,ε,(A,B),gcd) whose wall functional is an integer at tau (active or not)."""
    t0, t1 = tau
    out = []
    for i, j in combinations(range(4), 2):
        for eps in (1, -1):
            A = rows[i][0] - eps * rows[j][0]
            B = rows[i][1] - eps * rows[j][1]
            if (A, B) == (0, 0):
                continue
            val = F(A) * t0 + F(B) * t1
            if val.denominator == 1:
                out.append((i, j, eps, (A, B), gcd(A, B)))
    return out


def active_set(rows, tau, mu):
    t0, t1 = tau
    return [k for k in range(4) if nf(F(rows[k][0]) * t0 + F(rows[k][1]) * t1) == mu]


CASES = [(1, 2, 3, 2, 6, 9), (1, 3, 6, 1, 2, 4), (1, 2, 5, 3, 4, 6),
         (1, 3, 4, 3, 6, 8), (1, 2, 3, 3, 4, 6)]


def main():
    n_max_prim = 0     # primitive integer wall available at SOME maximizer vertex
    n_only_off = 0     # no primitive at any maximizer (would need a non-maximizer lonely point)
    for params in CASES:
        rows = jk(*params)
        sg = sat_gcd(rows)
        mu, verts = all_maximizers(rows)
        print(f"jk{params}: minor-gcd={sg}, mgap={mu} (>1/4: {mu > F(1,4)}), "
              f"{len(verts)} maximizer vertex(es)")
        found_prim = False
        for tau in verts:
            act = active_set(rows, tau, mu)
            walls = integer_walls_at(rows, tau)
            prim = [w for w in walls if w[4] == 1]
            tag = "  PRIMITIVE here" if prim else "  (none primitive)"
            print(f"   τ={tuple(str(x) for x in tau)} active={act}"
                  f" int-walls={[(i,j,e,g) for (i,j,e,_,g) in walls]}{tag}")
            if prim:
                print(f"       -> primitive witnesses: {[(i,j,e,AB) for (i,j,e,AB,g) in prim]}")
                found_prim = True
        if found_prim:
            n_max_prim += 1
            print("   VERDICT: primitive integer wall EXISTS at a maximizer (claim holds, τ=maximizer)")
        else:
            n_only_off += 1
            print("   VERDICT: NO primitive integer wall at ANY maximizer "
                  "(needs a non-maximizer lonely point, or refutes)")
        print()
    print(f"summary over {len(CASES)} flagged cases: "
          f"primitive-at-maximizer={n_max_prim}, none-at-maximizer={n_only_off}")


if __name__ == "__main__":
    main()
