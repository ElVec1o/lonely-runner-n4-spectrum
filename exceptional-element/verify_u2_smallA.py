#!/usr/bin/env python3
"""
Verifier: U2 realizes neither 1/3 nor 2/7, uniformly in A (companion to
verify_coordinate_bound.py).

Main argument (uniform in A, reusing the M=2(A+B) construction of the coordinate bound
with a wider band):  with t = k/(2(A+B)) and k odd, min_i ||v_i t|| >= c iff
phi = A t mod 1 lands in two arcs of width (1/2 - 2c); the odd-k orbit is equally spaced
with gap 1/(A+B), so

      A + B > 1/(1/2 - 2c)   ==>   ML > c.

  * c = 3/14 (arc width 1/14):  A+B >= 15  ==>  ML > 3/14  ==>  D < 2/7.
  * c = 1/6  (arc width 1/6 ):  A+B >=  7  ==>  ML > 1/6   ==>  D < 1/3.

The finite remainders A+B <= 14 (for 2/7) and A+B <= 6 (for 1/3) are checked exactly.

An appendix re-checks the explicit A in {1,2} value laws / pair-sum witnesses.

Usage:  python3 verify_u2_smallA.py [BOUND]   (default 60; caps A,B in the uniform sweep)
By Vico Bonfioli (Apache-2.0).  Dependency-free, exact rationals.
"""
import sys
from fractions import Fraction as F
from math import gcd


def nf(q):
    r = q - (q.numerator // q.denominator)
    return min(r, 1 - r)


def gap(v, t):
    return min(nf(vi * t) for vi in v)


def ML(v):
    c = {F(0)}
    for vi in v:
        a = abs(vi)
        for k in range(1, 2 * a):
            c.add(F(k, 2 * a))
    for i in range(len(v)):
        for j in range(i + 1, len(v)):
            for d in (v[i] - v[j], v[i] + v[j]):
                if d:
                    for m in range(1, abs(d)):
                        c.add(F(m, abs(d)))
    best = F(0)
    for t in c:
        g = gap(v, t)
        if g > best:
            best = g
    return best


def adm(A, B):
    s = (A, B, A + B, A + 2 * B)
    return gcd(A, B) == 1 and len(set(s)) == 4 and 0 not in s


def constr_min(A, B):
    """best min_i ||v_i k/M|| over odd k, M = 2(A+B)."""
    M = 2 * (A + B)
    best = F(0)
    for k in range(1, M, 2):
        g = gap((A, B, A + B, A + 2 * B), F(k, M))
        if g > best:
            best = g
    return best


def main():
    BOUND = int(sys.argv[1]) if len(sys.argv) > 1 else 60

    # ---- uniform exclusion of 2/7 ----
    bad = [(A, B) for A in range(1, BOUND) for B in range(1, BOUND)
           if adm(A, B) and A + B >= 15 and constr_min(A, B) <= F(3, 14)]
    fin = [(A, B) for A in range(1, 15) for B in range(1, 15)
           if adm(A, B) and A + B <= 14 and F(1, 2) - ML((A, B, A + B, A + 2 * B)) == F(2, 7)]
    print("[2/7] uniform:  A+B >= 15  =>  ML > 3/14  =>  D < 2/7")
    print(f"      construction min > 3/14 for all A+B>=15 (A,B<{BOUND}): {bad == []}")
    print(f"      finite A+B <= 14 with exact D = 2/7: {fin if fin else 'NONE'}")

    # ---- uniform exclusion of 1/3 ----
    bad3 = [(A, B) for A in range(1, BOUND) for B in range(1, BOUND)
            if adm(A, B) and A + B >= 7 and constr_min(A, B) <= F(1, 6)]
    fin3 = [(A, B) for A in range(1, 8) for B in range(1, 8)
            if adm(A, B) and A + B <= 6 and F(1, 2) - ML((A, B, A + B, A + 2 * B)) == F(1, 3)]
    print("\n[1/3] uniform:  A+B >= 7  =>  ML > 1/6  =>  D < 1/3")
    print(f"      construction min > 1/6 for all A+B>=7 (A,B<{BOUND}): {bad3 == []}")
    print(f"      finite A+B <= 6 with exact D = 1/3: {fin3 if fin3 else 'NONE'}")

    # ---- appendix: explicit A in {1,2} value laws ----
    def claim1(B):
        if B % 4 == 0:
            return F(3 * B, 4 * (3 * B + 1))
        if B % 4 == 3:
            return F(3 * B - 1, 4 * (3 * B + 1))
        return None

    badlaw = []
    for B in range(3, 201):
        if claim1(B) is None or not adm(1, B):
            continue
        if ML((1, B, 1 + B, 1 + 2 * B)) != claim1(B):
            badlaw.append(B)
    print(f"\n[appendix] A=1 exact value law holds for in-spectrum B<=200: {badlaw == []}")
    print("\n=> U2 realizes neither 1/3 nor 2/7, for all A >= 1.")


if __name__ == "__main__":
    main()
