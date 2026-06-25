#!/usr/bin/env python3
"""
Verifier for the U1-family characterization  (companion to verify_coordinate_bound.py).

U1 speeds: {A, 2A, 3A, B}, coprime (A,B).  Checks, with exact rationals:

  1. A >= 2  =>  D = 1/4  (the {A,2A,3A} deep hole always survives).
  2. A = 1   =>  D > 1/4 exactly when 4 | B, and then D = 1/4 + 1/(16 j + 4), B = 4 j.
  3. The realized k-set (D = 1/4 + 1/k) is exactly { k = 16 j + 4 : j >= 1 },
     i.e. k = 4 (mod 16);  in particular 12 (=1/3) and 28 (=2/7) are absent.

Usage:  python3 verify_u1_family.py [BOUND]      (default 120; caps B)
By Vico Bonfioli (Apache-2.0).
"""
import sys
from fractions import Fraction as F
from math import gcd


def nfrac(q):
    r = q - (q.numerator // q.denominator)
    return min(r, 1 - r)


def ML(v):
    cand = {F(0)}
    for vi in v:
        a = abs(vi)
        for k in range(1, 2 * a):
            cand.add(F(k, 2 * a))
    for i in range(len(v)):
        for j in range(i + 1, len(v)):
            for d in (v[i] - v[j], v[i] + v[j]):
                if d:
                    ad = abs(d)
                    for m in range(1, ad):
                        cand.add(F(m, ad))
    best = F(0)
    for t in cand:
        g = F(1, 2)
        for vi in v:
            dd = nfrac(vi * t)
            if dd < g:
                g = dd
            if g <= best:
                break
        if g > best:
            best = g
    return best


def admissible(A, B):
    s = (B, A, 2 * A, 3 * A)
    return gcd(A, B) == 1 and len(set(s)) == 4 and 0 not in s


def Dval(A, B):
    return F(1, 2) - ML((B, A, 2 * A, 3 * A))


def main():
    BOUND = int(sys.argv[1]) if len(sys.argv) > 1 else 120

    # ---- 1: A >= 2 => D = 1/4 ----
    bad = []
    for A in range(2, BOUND // 4 + 2):
        for B in range(1, 14 * A + 1):
            if admissible(A, B) and Dval(A, B) != F(1, 4):
                bad.append((A, B, str(Dval(A, B))))
    print("[1] A>=2 => D = 1/4 (deep hole survives):")
    print(f"    counterexamples: {bad if bad else 'NONE'}")

    # ---- 2 & 3: A = 1 family ----
    law_ok, viol = True, []
    ks = set()
    for B in range(1, BOUND + 1):
        if not admissible(1, B):
            continue
        d = Dval(1, B)
        if B % 4 == 0:
            j = B // 4
            pred = F(1, 4) + F(1, 16 * j + 4)
            if d != pred:
                law_ok = False
                viol.append((B, str(d), str(pred)))
            ks.add(16 * j + 4)
        else:
            if d > F(1, 4):
                viol.append((B, "D>1/4 but 4 does not divide B", str(d)))
                law_ok = False
    print("\n[2] A=1: D>1/4 iff 4|B, and D = 1/4 + 1/(16j+4):")
    print(f"    law holds: {law_ok}   violations: {viol if viol else 'none'}")

    print("\n[3] realized k-set at A=1 (D=1/4+1/k):")
    kl = sorted(ks)
    allmod4 = all((k - 4) % 16 == 0 for k in kl)
    print(f"    k = {kl[:14]}{' ...' if len(kl) > 14 else ''}")
    print(f"    all k = 4 (mod 16): {allmod4}")
    print(f"    12 (=1/3) realized? {12 in ks}     28 (=2/7) realized? {28 in ks}")

    # ---- cross: U2 small-A families also avoid 12, 28 ----
    def D2(A, B):
        return F(1, 2) - ML((A, B, A + B, A + 2 * B))
    bad2 = []
    for A in (1, 2):
        for B in range(1, BOUND + 1):
            s = (A, B, A + B, A + 2 * B)
            if gcd(A, B) != 1 or len(set(s)) < 4 or 0 in s:
                continue
            x = D2(A, B) - F(1, 4)
            if x > 0 and x.numerator == 1 and x.denominator in (12, 28):
                bad2.append((A, B, x.denominator))
    print("\n[cross] U2 families A in {1,2}: any realization of k=12 or k=28?")
    print(f"    {bad2 if bad2 else 'NONE (U2 also excludes 1/3 and 2/7 at A<=2)'}")


if __name__ == "__main__":
    main()
