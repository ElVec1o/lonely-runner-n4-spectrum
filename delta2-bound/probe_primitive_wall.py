#!/usr/bin/env python3
"""
Independent probe of the PrimitiveActiveWall residual for  δ₂(4) ≤ 3/14.

The Lean reduction (`delta2sat_le_of_primitiveActiveWall`) shows the bound follows from:

    PrimitiveActiveWall :  every SATURATED rank-2 band M (minor-gcd 1) with mD M < 1/4
        has a maximizer τ (mgap M τ > 1/4) and an active wall (i,j,ε)
        whose integer vector (A,B) = (M_i 0 − ε M_j 0,  M_i 1 − ε M_j 1) is PRIMITIVE: gcd(A,B)=1.

Everything else in the chain is already PROVED in Lean (maximizer attainment, the active wall with
integer functional value, the whole reduction). The ONLY open piece is primitivity, and the
`Saturated2` guard is essential: the un-guarded version is FALSE (2·Lrz: all entries even).

This script does two things, from scratch, in exact integer/rational arithmetic:
  (1) RE-DERIVES the claim independently: over saturated bands with mgap>1/4, it finds the
      maximizer τ, the active runners, the active walls, and checks that at least one is primitive.
  (2) MINES THE MECHANISM: for every active wall it records gcd(A,B), and it asks the structural
      question a proof must answer — *is the primitive active wall forced, and by what?* It tabulates
      (a) how often the FIRST active wall is already primitive, (b) whether a non-primitive active
      wall can occur at all on a saturated band (if so, primitivity needs a CHOICE among walls),
      (c) the relationship gcd(active walls) vs the minor-gcd (=1 for saturated).

Run:  python3 probe_primitive_wall.py [N]      (default N = 14)
"""
import sys
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
    """distance from x to nearest integer, exact Fraction in [0,1/2]."""
    fr = x - (x.numerator // x.denominator)
    return min(fr, 1 - fr)


def maximizers(rows):
    """Exact: enumerate wall-line intersection vertices in [0,1)^2; return
    (best_mgap, list of (τ0,τ1) attaining it)."""
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
    verts = []
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
                            verts = [(t0, t1)]
                        elif g == best:
                            verts.append((t0, t1))
    return best, verts


def active_walls(rows, tau, mu):
    """At maximizer tau with min-distance mu, return the list of active walls
    (i,j,eps,(A,B)) for active runner pairs whose wall functional is an integer."""
    t0, t1 = tau
    active = [k for k in range(4) if nf(F(rows[k][0]) * t0 + F(rows[k][1]) * t1) == mu]
    out = []
    for i, j in combinations(active, 2):
        for eps in (1, -1):
            A = rows[i][0] - eps * rows[j][0]
            B = rows[i][1] - eps * rows[j][1]
            if (A, B) == (0, 0):
                continue
            val = F(A) * t0 + F(B) * t1
            if val.denominator == 1:  # integer functional => genuine active wall
                out.append((i, j, eps, (A, B)))
    return active, out


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 14

    tested = 0
    no_active_wall = 0
    no_primitive = []          # saturated bands where NO active wall is primitive (would refute)
    first_primitive = 0        # at least one active wall, and the FIRST listed one is primitive
    needs_choice = 0           # primitive exists but some active wall is non-primitive (choice matters)
    gcd_hist = {}              # histogram of gcd(A,B) over ALL active walls
    tau_denoms = {}            # histogram of lcm(den τ0, den τ1) at maximizers

    # Case A generic saturated bands; mgap>1/4 is the relevant regime (mD<1/4).
    for a in range(1, N + 1):
        for b in range(a + 1, N + 1):
            for c in range(b + 1, N + 1):
                for d in range(1, N + 1):
                    for e in range(d + 1, N + 1):
                        for f in range(e + 1, N + 1):
                            rows = jk(a, b, c, d, e, f)
                            if sat_gcd(rows) != 1:
                                continue
                            mu, verts = maximizers(rows)
                            if mu <= F(1, 4):       # mD = 1/2 - mu >= 1/4: not in the regime
                                continue
                            tested += 1
                            # use the first maximizer vertex
                            tau = verts[0]
                            dd = (tau[0].denominator * tau[1].denominator
                                  // gcd(tau[0].denominator, tau[1].denominator))
                            tau_denoms[dd] = tau_denoms.get(dd, 0) + 1
                            act, walls = active_walls(rows, tau, mu)
                            if not walls:
                                no_active_wall += 1
                                continue
                            gs = [gcd(A, B) for (_, _, _, (A, B)) in walls]
                            for g in gs:
                                gcd_hist[g] = gcd_hist.get(g, 0) + 1
                            if 1 not in gs:
                                no_primitive.append((a, b, c, d, e, f, gs))
                            else:
                                if gs[0] == 1:
                                    first_primitive += 1
                                if any(g != 1 for g in gs):
                                    needs_choice += 1

    print(f"=== PrimitiveActiveWall probe, saturated Case-A bands, entries ≤ {N} ===")
    print(f"saturated bands with mgap>1/4 (mD<1/4) tested : {tested}")
    print(f"  with NO integer active wall at maximizer     : {no_active_wall}")
    print(f"  with a PRIMITIVE active wall (claim holds)    : {tested - no_active_wall - len(no_primitive)}")
    print(f"  with NO primitive active wall (REFUTES!)      : {len(no_primitive)}")
    if no_primitive:
        print(f"    counterexamples: {no_primitive[:8]}")
    print()
    print("--- mechanism mining ---")
    print(f"first listed active wall already primitive     : {first_primitive}")
    print(f"primitive exists but a non-primitive wall too  : {needs_choice}")
    print(f"gcd(A,B) histogram over ALL active walls       : {dict(sorted(gcd_hist.items()))}")
    print(f"maximizer denominator lcm(den τ0,τ1) histogram : {dict(sorted(tau_denoms.items()))}")
    print()
    if not no_primitive:
        print("✓ every saturated band in the regime has a primitive active wall (claim re-derived).")
        if needs_choice == 0:
            print("  STRUCTURE: EVERY active wall is primitive — primitivity is not a choice,")
            print("  it is forced for all active pairs. A proof need not select the wall.")
        else:
            print("  STRUCTURE: some bands have a non-primitive active wall too, so a proof must")
            print("  CHOOSE the primitive one — primitivity is a property of the right pair, not all.")


if __name__ == "__main__":
    main()
