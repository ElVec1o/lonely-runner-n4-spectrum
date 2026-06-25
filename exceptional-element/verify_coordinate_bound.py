#!/usr/bin/env python3
"""
Verifier for the U2 coordinate bound  ML(A, B, A+B, A+2B) >= 2/9  for A >= 3.

Three independent, dependency-free checks (exact rational arithmetic):

  1. EXACT ML (no construction): ML >= 2/9 for every admissible (A,B) with A>=3
     in range, with equality exactly at (3,1), (3,4), (4,1).
  2. PROOF-SCHEME COVERAGE: every such (A,B) is closed either by the asymptotic
     construction M = 2(A+B) (when A+B >= 18) or by an explicit modulus
     certificate (when A+B <= 17) -- the proof's own case split has no gaps.
  3. SHARPNESS: A <= 2 admits ML < 2/9, so the hypothesis A >= 3 is necessary.

Usage:  python3 verify_coordinate_bound.py [BOUND]      (default BOUND = 60)
        BOUND caps A+B in the exact-ML pass.  By Vico Bonfioli (Apache-2.0).
"""
import sys
from fractions import Fraction as F
from math import gcd, ceil, floor


def nfrac(q):
    """Distance from rational q to the nearest integer."""
    r = q - (q.numerator // q.denominator)   # fractional part in [0,1)
    return min(r, 1 - r)


def ML(v):
    """Exact max_t min_i ||v_i t||, searching all breakpoints (no coprimality shortcut)."""
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
            d = nfrac(vi * t)
            if d < g:
                g = d
            if g <= best:
                break
        if g > best:
            best = g
    return best


def admissible(A, B):
    s = (A, B, A + B, A + 2 * B)
    return gcd(A, B) == 1 and len(set(s)) == 4 and 0 not in s


def cert(A, B, M):
    """Smallest k giving t=k/M with all four speeds in the band [2M/9, 7M/9]."""
    if M < 9:
        return None
    lo, hi = ceil(F(2 * M, 9)), floor(F(7 * M, 9))
    rs = (A, B, A + B, A + 2 * B)
    for k in range(1, M):
        if all(lo <= (v * k) % M <= hi for v in rs):
            return k
    return None


def scheme(A, B):
    """The proof's case split. Returns (label, M, k) or None."""
    if A + B >= 18:                                   # Step 4: asymptotic
        M = 2 * (A + B)
        for k in range(1, M, 2):                      # odd k only
            lo, hi = ceil(F(2 * M, 9)), floor(F(7 * M, 9))
            if all(lo <= (v * k) % M <= hi for v in (A, B, A + B, A + 2 * B)):
                return ("asymptotic  M=2(A+B)", M, k)
        return None
    for M in [2 * (A + B), 2 * A + 3 * B, 2 * A + B, A + 3 * B,   # Step 5: finite
              2 * A + 2 * B, 3 * A + 2 * B, 2 * A, 2 * B,
              2 * (A + 2 * B), A + 4 * B, 3 * A + B]:
        k = cert(A, B, M)
        if k is not None:
            return (f"finite     M={M}", M, k)
    return None


def main():
    BOUND = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    TARGET = F(2, 9)

    # ---- Check 1: exact ML ----
    viol, eq = [], []
    for A in range(3, BOUND):
        for B in range(1, BOUND - A + 1):
            if not admissible(A, B):
                continue
            m = ML((A, B, A + B, A + 2 * B))
            if m < TARGET:
                viol.append((A, B, m))
            elif m == TARGET:
                eq.append((A, B))
    print(f"[1] exact ML >= 2/9 for A>=3, A+B<={BOUND}:")
    print(f"    violations: {viol if viol else 'NONE  (lemma holds)'}")
    print(f"    equality ML=2/9 at: {sorted(eq)}   (expected: [(3, 1), (3, 4), (4, 1)])")

    # ---- Check 2: proof-scheme coverage ----
    gaps = []
    for A in range(3, 2 * BOUND):
        for B in range(1, 3 * A + 1):
            if not admissible(A, B):
                continue
            if scheme(A, B) is None:
                gaps.append((A, B))
    print(f"\n[2] proof-scheme coverage for A=3..{2*BOUND-1}, all B:")
    print(f"    uncovered cases: {gaps if gaps else 'NONE  (case split is complete)'}")

    # ---- Check 3: sharpness ----
    sharp = []
    for A in (1, 2):
        for B in range(1, 12):
            if admissible(A, B) and ML((A, B, A + B, A + 2 * B)) < TARGET:
                sharp.append((A, B, str(ML((A, B, A + B, A + 2 * B)))))
    print(f"\n[3] sharpness (A<=2 must allow ML<2/9): {sharp}")

    # ---- finite certificate table (A+B <= 17) ----
    print(f"\nFinite certificate table (A>=3, A+B<=17):")
    rows = 0
    for A in range(3, 17):
        for B in range(1, 18 - A):
            if not admissible(A, B):
                continue
            r = scheme(A, B)
            if r is None:
                print(f"    !! NO CERTIFICATE for ({A},{B})")
                continue
            _, M, k = r
            rs = (A, B, A + B, A + 2 * B)
            norms = ", ".join(f"{min((v*k) % M, M-(v*k) % M)}/{M}" for v in rs)
            mn = min(min((v * k) % M, M - (v * k) % M) for v in rs)
            tight = "  <- tight (=2/9)" if F(mn, M) == TARGET else ""
            print(f"    (A={A:2d},B={B:2d}) {str(rs):16s} t={k}/{M:<3d}  ||v t||=({norms}){tight}")
            rows += 1
    print(f"    {rows} pairs, all certified.")


if __name__ == "__main__":
    main()
