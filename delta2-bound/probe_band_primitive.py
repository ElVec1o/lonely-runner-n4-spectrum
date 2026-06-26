#!/usr/bin/env python3
"""
Refined PrimitiveActiveWall probe, RESTRICTED to the relevant band.

The δ₂(4) ≤ 3/14 statement is  S₂(4) ∩ (3/14, 1/4) = ∅: no saturated band has
mD ∈ (3/14, 1/4), equivalently mgap ∈ (1/4, 2/7).  Bands with mgap ≥ 2/7 (mD ≤ 3/14)
ALREADY satisfy the bound — the wall argument is only NEEDED for mgap ∈ (1/4, 2/7).

The current Lean `PrimitiveActiveWall` hypothesises only `mD < 1/4` (mgap > 1/4), which
also covers the deep-lonely regime mgap ≥ 2/7 where (this probe shows) the maximizer's
active wall can be NON-primitive — e.g. jk(1,2,5,3,4,6), mgap=5/12.  So the natural
"active wall at the maximizer is primitive" claim is FALSE as stated.

HYPOTHESIS UNDER TEST:  restricted to the relevant band  1/4 < mgap < 2/7,  every saturated
band has a PRIMITIVE active wall at SOME global maximizer (active-active pair).  If true, the
clean residual is the BAND-restricted statement, and the over-broad `mD < 1/4` hypothesis is
the bug.

Run:  python3 probe_band_primitive.py [N]      (default N = 16)
"""
import sys
from math import gcd
from itertools import combinations
from fractions import Fraction as F

LO, HI = F(1, 4), F(2, 7)   # relevant band for mgap


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
                            best, verts = g, {(t0, t1)}
                        elif g == best:
                            verts.add((t0, t1))
    return best, sorted(verts)


def active_walls_at(rows, tau, mu):
    """active-active integer walls (the intended deep-hole structure)."""
    t0, t1 = tau
    act = [k for k in range(4) if nf(F(rows[k][0]) * t0 + F(rows[k][1]) * t1) == mu]
    out = []
    for i, j in combinations(act, 2):
        for eps in (1, -1):
            A = rows[i][0] - eps * rows[j][0]
            B = rows[i][1] - eps * rows[j][1]
            if (A, B) == (0, 0):
                continue
            if (F(A) * t0 + F(B) * t1).denominator == 1:
                out.append(gcd(A, B))
    return act, out


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    tested = 0
    refute = []          # saturated band in (1/4,2/7) with NO primitive active wall at any maximizer
    nactive_hist = {}    # active-set size at the maximizer that carries the primitive wall
    for a in range(1, N + 1):
        for b in range(a + 1, N + 1):
            for c in range(b + 1, N + 1):
                for d in range(1, N + 1):
                    for e in range(d + 1, N + 1):
                        for f in range(e + 1, N + 1):
                            rows = jk(a, b, c, d, e, f)
                            if sat_gcd(rows) != 1:
                                continue
                            mu, verts = all_maximizers(rows)
                            if not (LO < mu < HI):     # relevant band only
                                continue
                            tested += 1
                            best_sz = None
                            for tau in verts:
                                act, gs = active_walls_at(rows, tau, mu)
                                if 1 in gs:
                                    best_sz = len(act)
                                    break
                            if best_sz is None:
                                refute.append((a, b, c, d, e, f, str(mu)))
                            else:
                                nactive_hist[best_sz] = nactive_hist.get(best_sz, 0) + 1
    print(f"=== relevant-band probe (1/4 < mgap < 2/7), saturated Case-A, entries ≤ {N} ===")
    print(f"saturated bands with mgap in (1/4, 2/7) tested        : {tested}")
    print(f"  with a PRIMITIVE active wall at a maximizer          : {tested - len(refute)}")
    print(f"  with NONE (refutes band-restricted claim)            : {len(refute)}")
    if refute:
        print(f"    refuters: {refute[:12]}")
    print(f"  active-set size at the primitive-carrying maximizer  : {dict(sorted(nactive_hist.items()))}")
    if not refute:
        print("  ✓ on the RELEVANT band, every saturated form has a primitive active wall at a maximizer.")
        print("    => the clean residual is the BAND-restricted statement; `mD<1/4` over-reaches.")


if __name__ == "__main__":
    main()
