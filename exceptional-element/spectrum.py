#!/usr/bin/env python3
"""
Reproducible exact computation: D = 2/7 appears to be an exceptional element of
the relative Lonely Runner spectrum  S_1(4) cap (1/4, 1/2].

Background (Jain-Kravitz, "Relative Lonely Runner spectra", arXiv:2411.12684):
Theorem 1.3 states that  S_1(4) cap (1/4, 1/2]  has finite symmetric difference
with the set  1/4 + (1/4) Prog(2,3), and that  S_1(4) cap (1/4, 1/2] = S_1(U1) ∪ S_1(U2),
where, up to symmetry,
    U1 = < (0,1,2,3), (1,0,0,0) >   and   U2 = < (1,0,1,1), (1,1,0,2) >
are the only 2-dimensional subtori of (R/Z)^4 with D-value 1/4.  The paper notes the
exact symmetric difference is an uncomputed finite calculation, with numerical
experiments suggesting no exceptional elements.

The 1-dimensional subtori of U1, U2 are parameterized by coprime integers (A,B):
    U1:  v(A,B) = (B, A, 2A, 3A)
    U2:  v(A,B) = (A+B, B, A, A+2B)
and we compute the maximal-loneliness deficit
    D(v) = 1/2 - max_t min_i || v_i t ||         ( ||x|| = distance to nearest integer ).

Empirically, over |A|,|B| <= N, the realized values  D = 1/4 + 1/k  with k in [20, 84]
are exactly  { k congruent 4 mod 8 }  EXCEPT  k = 28.  That is, every nearby value
1/4 + 1/20, 1/4 + 1/36, 1/4 + 1/44, ... occurs, but  1/4 + 1/28 = 2/7  does NOT.

This script computes D exactly (rational arithmetic, no floating point, no
dependencies) and prints the realized k for both subtori.  Run:  python3 spectrum.py [N]

IMPORTANT correctness note: the candidate t = p/q maximizing min_i ||v_i t|| must be
searched over ALL k/d for d in the breakpoint set { |v_i|, |v_i +- v_j|, 2|v_i| } and
ALL 0 <= k < d -- NOT only gcd(k,d)=1 -- because an optimal reduced t = p/q can have
q not in that set while q | d for some d in it (e.g. t = 1/4 = 5/20 with 20 = 7+13).
A coprimality shortcut here silently overestimates D.
"""

from fractions import Fraction as F
from math import gcd
import sys


def nearest_int_dist(x: F) -> F:
    """|| x ||: distance from rational x to the nearest integer, as a Fraction."""
    frac = x - (x.numerator // x.denominator)   # in [0,1)
    if frac < 0:
        frac += 1
    return min(frac, 1 - frac)


def D_value(v):
    """Exact D-value of an integer speed-vector v = (v_0, ..., v_{n-1}):
       D(v) = 1/2 - max_t min_{i : v_i != 0} || v_i t ||."""
    speeds = sorted({abs(s) for s in v if s != 0})
    if not speeds:
        return F(0)
    # Candidate denominators for the optimal t.
    denoms = {1}
    for i in range(len(speeds)):
        denoms.add(2 * speeds[i])
        for j in range(i + 1, len(speeds)):
            denoms.add(speeds[i] + speeds[j])
            denoms.add(abs(speeds[i] - speeds[j]))
    best = F(0)
    for d in denoms:
        for k in range(d + 1):          # ALL k -- see correctness note above
            t = F(k, d)
            m = min(nearest_int_dist(c * t) for c in speeds)
            if m > best:
                best = m
    return F(1, 2) - best


def realized_k(name, direction, N):
    """For coprime (A,B) with |A|,|B| <= N, collect k such that D = 1/4 + 1/k."""
    ks = {}
    k28 = None
    for A in range(-N, N + 1):
        for B in range(-N, N + 1):
            if A == 0 and B == 0:
                continue
            if gcd(abs(A), abs(B)) != 1:
                continue
            D = D_value(direction(A, B))
            if D > F(1, 4):
                diff = D - F(1, 4)
                k = 1 / diff
                if k.denominator == 1:
                    k = int(k)
                    ks[k] = ks.get(k, 0) + 1
                    if k == 28 and k28 is None:
                        k28 = (A, B)
    realized = sorted(k for k in ks if 20 <= k <= 84)
    print(f"{name}:")
    print(f"    realized k in [20,84]: {realized}")
    print(f"    k = 28 (D = 2/7): {'PRESENT at (A,B)=' + str(k28) if k28 else 'ABSENT'}")
    return ks


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    print(f"Exact relative spectra S_1(U1), S_1(U2) cap (1/4, 1/2],  |A|,|B| <= {N}\n")
    realized_k("S_1(U1)   v(A,B) = (B, A, 2A, 3A)", lambda A, B: (B, A, 2 * A, 3 * A), N)
    realized_k("S_1(U2)   v(A,B) = (A+B, B, A, A+2B)", lambda A, B: (A + B, B, A, A + 2 * B), N)
    print("\nNote: S_1(4) cap (1/4,1/2] = S_1(U1) ∪ S_1(U2).  Over this range the union is")
    print("{ k congruent 4 mod 8, k >= 20 } with the single omission k = 28, i.e. D = 2/7.")


if __name__ == "__main__":
    main()
