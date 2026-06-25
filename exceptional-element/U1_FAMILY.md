# The U¹ family is completely characterized: `D > 1/4 ⟺ A = 1 and 4 ∣ B`

Companion to [`COORDINATE_BOUND.md`](COORDINATE_BOUND.md) (which handles `U²`). Together
they show that neither `1/3` nor `2/7` is realized by any `1`-dimensional subtorus of
`U¹` or `U²`.

By Vico Bonfioli. Proof developed with assistance from Anthropic's Claude. Apache-2.0.
Mechanically re-checked by [`verify_u1_family.py`](verify_u1_family.py).

## Setup

`U¹` is parameterized by coprime `(A, B)` as `v(A,B) = (B, A, 2A, 3A)`; as a speed
multiset this is `{A, 2A, 3A, B}`. As before `ML(v) = max_t min_i ‖vᵢ t‖` and
`D = 1/2 − ML`, and we care about `D ∈ (1/4, 1/2]`, i.e. `ML < 1/4`.

The three speeds `A, 2A, 3A` are the `{1,2,3}` system scaled by `A`. Its loneliness
`max_s min(‖s‖, ‖2s‖, ‖3s‖) = 1/4` is attained **only** at `s ≡ 1/4, 3/4 (mod 1)` — i.e.
at `t = m/(4A)` with `m` odd. So `ML({A,2A,3A,B}) ≤ 1/4` always, with equality iff some
odd-`m` deep hole `t = m/(4A)` also satisfies `‖B·t‖ ≥ 1/4`.

## Theorem

> Let `A, B` be coprime with `{A, 2A, 3A, B}` distinct and nonzero. Then
> ```
>   D(v(A,B)) > 1/4   ⟺   A = 1  and  4 ∣ B.
> ```
> Writing `B = 4j` (`j ≥ 1`) in that case,
> ```
>   ML({1,2,3,4j}) = j/(4j+1),     D = 1/4 + 1/(16j+4).
> ```
> Consequently `S(U¹) ∩ (1/4, 1/2] = { 1/4 + 1/(16j+4) : j ≥ 1 }`, i.e. exactly the
> values `1/4 + 1/k` with `k ≡ 4 (mod 16)`.

Since `28 ≡ 12 (mod 16)` and `12 ≡ 12 (mod 16)`, neither `k = 28` (`D = 2/7`) nor
`k = 12` (`D = 1/3`) is `≡ 4 (mod 16)`: **`U¹` realizes neither `2/7` nor `1/3`.**

## Proof

### Part 1 — `A ≥ 2 ⟹ D = 1/4` (a deep hole always survives)

`ML ≤ 1/4` is automatic. For `ML ≥ 1/4` it suffices to find an **odd** `m ∈ [1, 4A)`
with `‖B·m/(4A)‖ ≥ 1/4`, equivalently `Bm mod 4A ∈ [A, 3A]`; then `t = m/(4A)` gives
`‖At‖ = ‖m/4‖ = 1/4`, `‖2At‖ = ‖m/2‖ = 1/2`, `‖3At‖ = ‖3m/4‖ = 1/4` (all since `m` is
odd), and `‖Bt‖ ≥ 1/4`, so `min_i ‖vᵢ t‖ = 1/4`.

Let `r = B mod 4A`. Because `gcd(A,B) = 1` and `A ≥ 2`, we have `A ∤ B`, so `4A ∤ B` and
`r ≠ 0`.

* If `r ∈ [A, 3A]`, take `m = 1`.
* If `r ∈ [1, A−1]`, the interval `[A/r, 3A/r]` has length `2A/r > 2`, hence contains an
  odd integer `m`; then `rm ∈ [A, 3A]` with no wraparound (`rm ≤ 3A < 4A`), so
  `Bm mod 4A = rm ∈ [A, 3A]`.
* If `r ∈ [3A+1, 4A−1]`, put `r' = 4A − r ∈ [1, A−1]`; since `Bm ≡ −r'm (mod 4A)` and
  `‖·‖` is even, the previous case applied to `r'` gives the required `m`.

So such an `m` always exists, `ML = 1/4`, and `D = 1/4 ∉ (1/4, 1/2]`. The only way to
block every deep hole is `r = 0`, i.e. `4A ∣ B`, which forces `A ∣ B` and hence `A = 1`.

### Part 2 — `A = 1`: `D > 1/4 ⟺ 4 ∣ B`, and the exact value

At `A = 1` the deep holes are `t = 1/4, 3/4`, where `‖Bt‖ = ‖B/4‖`. Both are blocked
(`‖B/4‖ < 1/4` and `‖3B/4‖ < 1/4`) iff `B ≡ 0 (mod 4)`; otherwise one survives and
`D = 1/4`. So `D > 1/4 ⟺ 4 ∣ B`.

For `B = 4j` (`j ≥ 1`), the time `t = j/(4j+1)` gives
```
‖1·t‖ = j/(4j+1),   ‖2·t‖ = 2j/(4j+1),   ‖3·t‖ = (j+1)/(4j+1),   ‖4j·t‖ = j/(4j+1),
```
the last because `4j·j/(4j+1) = j − j/(4j+1)`. The minimum is `j/(4j+1)`, so
`ML ≥ μ` with `μ := j/(4j+1)`.

For the matching **upper bound** `ML ≤ μ`, fix any `t` and write `g(t) = min` over the
four runners. If `min(‖t‖, ‖2t‖, ‖3t‖) ≤ μ` then `g(t) ≤ μ`. Otherwise `t` lies where the
scaled `{1,2,3}` system exceeds `μ`. Since `μ < 1/4` and `min(‖t‖,‖2t‖,‖3t‖)` attains its
maximum `1/4` only at `t = 1/4, 3/4`, that region is exactly the two intervals
`(μ, (1−μ)/3)` and its reflection `(1 − (1−μ)/3, 1 − μ)` (the binding constraints near
`1/4` are `‖t‖ > μ`, i.e. `t > μ`, and `‖3t‖ > μ`, i.e. `t < (1−μ)/3`). On the closure
`[μ, (1−μ)/3]` the value `4j·t` runs over `[j − μ, j + j/(12j+3)]`, and since
`j/(12j+3) < j/(4j+1) = μ`, this lies in `[j − μ, j + μ]`; hence `‖4j·t‖ = |4j·t − j| ≤ μ`,
with equality only at the left endpoint `t = μ`. So `g(t) ≤ ‖4j·t‖ ≤ μ` on this interval
too (reflection identical). Therefore `ML = μ`, attained at `t = μ = j/(4j+1)`, and
`D = 1/2 − μ = 1/4 + 1/(16j+4)`. ∎

So `U¹` realizes **exactly** the deficits `1/4 + 1/(16j+4)`, `j ≥ 1`. The value
`D = 1/4 + 1/(16j+4)` is also confirmed independently for `1 ≤ j ≤ 12`, and the
upper-bound region/`‖4j·t‖` claims for `1 ≤ j ≤ 8`, in the verifier.

## Verification

`verify_u1_family.py` checks, with exact rational arithmetic and no dependencies:

1. `A ≥ 2 ⟹ D = 1/4` (the deep hole survives) for all admissible `(A,B)` in range;
2. at `A = 1`, `D > 1/4` exactly when `4 ∣ B`, and then `D = 1/4 + 1/(16j+4)`;
3. the realized `k`-set is exactly `{k ≡ 4 (mod 16), k ≥ 20}`, so `12` and `28` are
   absent.

```
python3 verify_u1_family.py 120
```

## Lean formalization

The constructive half — `ML(1,2,3,4j) ≥ j/(4j+1)`, hence `D ≤ 1/4 + 1/(16j+4)` — is
formalized sorry-free in Lean 4 / Mathlib as `ML_u1_family_ge` and `D_u1_family_le`
(`LRSpectrumFull/U1FamilyBound.lean` in the n=3 Lean project). The fourth runner is
discharged exactly as in the paper: `4j·t = -t + j`, so
`nearestIntDist (4j·t) = nearestIntDist (-t) = nearestIntDist t` via the library lemmas
`nearestIntDist_add_int` and `nearestIntDist_neg`. Axiom footprint:
`[propext, Classical.choice, Quot.sound]` — the three standard axioms, no `sorryAx`, no
`native_decide`.

## Combined consequence (U¹ and U²)

* `U¹` (here): realizes exactly `k ≡ 4 (mod 16)`; excludes `12` and `28`.
* `U²`, `A ≥ 3` ([`COORDINATE_BOUND.md`](COORDINATE_BOUND.md)): `D ≤ 5/18 < 2/7 < 1/3`.
* `U²`, all `A` ([`U2_SMALL_A.md`](U2_SMALL_A.md)): the `M = 2(A+B)` construction at band
  `3/14` resp. `1/6` gives `A+B ≥ 15 ⟹ D < 2/7` and `A+B ≥ 7 ⟹ D < 1/3`, with the finite
  remainders checked exactly. This covers `A = 1, 2` (and reproves `A ≥ 3` more weakly).

Therefore **`D = 1/3` and `D = 2/7` are realized by no `1`-dimensional subtorus of
`U¹ ∪ U²`.** Combined with Jain–Kravitz Theorem 1.3 (that `S₁(4) ∩ (1/4, 1/2]` is, up to
a finite symmetric difference, carried by these two subtori), this identifies `1/3` and
`2/7` as exceptional elements. Establishing that they are the *only* exceptional elements
additionally requires the forward (realization) direction for every other progression
value, which is documented numerically in [`README.md`](README.md) and is not proved
here.
