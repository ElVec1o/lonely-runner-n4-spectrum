# A candidate exceptional element D = 2/7 in S₁(4) ∩ (1/4, 1/2]

An exact, dependency-free computation suggesting that **D = 2/7** is an exceptional
element of the relative Lonely Runner spectrum `S₁(4) ∩ (1/4, 1/2]` — i.e. a member of
the finite symmetric difference in **Theorem 1.3 of Jain–Kravitz**, *Relative Lonely
Runner spectra* ([arXiv:2411.12684](https://arxiv.org/abs/2411.12684)).

By Vico Bonfioli. Computation assisted by Anthropic's Claude. Apache-2.0.

## Setup

For an integer speed-vector `v`, the maximal-loneliness deficit is
`D(v) = 1/2 − max_t min_i ‖v_i t‖`, where `‖x‖` is the distance to the nearest integer.
Jain–Kravitz Theorem 1.3 gives `S₁(4) ∩ (1/4, 1/2] = S₁(U¹) ∪ S₁(U²)` up to a finite
symmetric difference with `1/4 + (1/4)·Prog(2,3)`, where (up to symmetry)

```
U¹ = ⟨(0,1,2,3), (1,0,0,0)⟩,   U² = ⟨(1,0,1,1), (1,1,0,2)⟩
```

are the only 2-dimensional subtori of `(ℝ/ℤ)⁴` with `D = 1/4`. Their 1-dimensional
subtori are parameterized by coprime `(A,B)`:

```
U¹:  v(A,B) = (B, A, 2A, 3A)
U²:  v(A,B) = (A+B, B, A, A+2B)
```

## The observation

Computing `D(v(A,B))` exactly over `|A|, |B| ≤ 600`, the realized values of the form
`D = 1/4 + 1/k` with `k ∈ [20, 84]` are

```
S₁(U¹):  k ∈ {20, 36, 52, 68, 84}
S₁(U²):  k ∈ {20, 36, 44, 52, 60, 68, 76, 84}
union :  k ∈ {20, 36, 44, 52, 60, 68, 76, 84}
```

i.e. exactly `{ k ≡ 4 (mod 8), k ≥ 20 }` **with the single omission `k = 28`**. Every
neighbour `1/4 + 1/20, 1/4 + 1/36, 1/4 + 1/44, …` occurs (the smallest realizations
appear already by `|A|,|B| ≈ 40`), but `1/4 + 1/28 = 2/7` does **not** occur anywhere in
the searched range, for either subtorus.

This is consistent with `D = 2/7` being an exceptional element of the finite symmetric
difference in Theorem 1.3. It is an **observation**, not a proof: ruling out a
realization at `|A|, |B| > 600` requires the structural description of the realizable
set, not an enlarged search.

## Reproduce

```
python3 spectrum.py 60      # exact rationals, no dependencies, ~seconds
```

(Use a larger argument for a wider range; runtime grows with it.)

## Correctness note (a real pitfall)

The candidate `t = p/q` maximizing `min_i ‖v_i t‖` must be searched over **all** `k/d`
for `d` in the breakpoint set `{|v_i|, |v_i ± v_j|, 2|v_i|}` and **all** `0 ≤ k < d` —
*not* only `gcd(k,d) = 1`. An optimal reduced `t = p/q` can have `q` outside that set yet
`q | d` for some `d` in it (e.g. `t = 1/4 = 5/20`, with `20 = 7 + 13`). A coprimality
shortcut silently **overestimates** `D` and manufactures spurious large-`D` values
(including a false `k = 28`). The verifier here avoids that shortcut; the result was
cross-checked against an independent fine-grid evaluation.

## Acknowledgements

Computation developed by Vico Bonfioli with substantial assistance from Anthropic's
Claude. The arithmetic is elementary and fully reproducible from `spectrum.py`;
nothing here relies on trusting that assistance.
