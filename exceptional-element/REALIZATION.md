# The forward (realization) direction: status

The exclusion direction ([`COORDINATE_BOUND.md`](COORDINATE_BOUND.md),
[`U2_SMALL_A.md`](U2_SMALL_A.md), [`U1_FAMILY.md`](U1_FAMILY.md)) proves that `1/3` and
`2/7` are realized by **no** subtorus of `U¹ ∪ U²`. The *forward* direction is the
complementary claim needed for "the symmetric difference is **exactly** `{1/3, 2/7}`":

> every `k ≡ 4 (mod 8)`, `k ≥ 20`, `k ≠ 28` is realized, i.e. `1/4 + 1/k ∈ S(U¹) ∪ S(U²)`.

By Vico Bonfioli. Developed with assistance from Anthropic's Claude. Apache-2.0.
Checked by [`verify_realization.py`](verify_realization.py).

## What splits the problem: residues mod 16

The target `k ≡ 4 (mod 8)` splits into `k ≡ 4 (mod 16)` and `k ≡ 12 (mod 16)`. Each is
realized by a single clean coprime family.

### `k ≡ 4 (mod 16)` — closed, with full proof

`U¹` realizes **exactly** these. By [`U1_FAMILY.md`](U1_FAMILY.md) (both bounds proved),
`D(1,2,3,4j) = 1/4 + 1/(16j+4)`, so `(A,B) = (1, 4j)` realizes `k = 16j+4` for every
`j ≥ 1`:

```
k = 20, 36, 52, 68, 84, 100, …   (k ≡ 4 mod 16, k ≥ 20)
```

This half is a single explicit infinite family with a fully proved exact value (the `U¹`
upper bound uses the rigidity of the scaled `{1,2,3}` system). It is also the half
formalized (lower bound) in Lean.

### `k ≡ 12 (mod 16)` — one clean coprime family; exact value pending one lemma

A single `U²` family realizes all of these, with **no coprimality gaps**:

```
k = 44:            (A, B) = (1, 7)        speeds {1, 7, 8, 9}
k = 16m + 60:      (A, B) = (4m+3, 8)     speeds {4m+3, 8, 4m+11, 4m+19},  m ≥ 0
```

`k = 16m+60` for `m = 0,1,2,…` enumerates `60, 76, 92, 108, …` = exactly
`{k ≡ 12 (mod 16), k ≥ 60}`, and `gcd(4m+3, 8) = 1` **unconditionally** because `4m+3` is
odd and `8 = 2³`. Putting the free parameter in the odd coordinate and fixing `B` to a
power of two is what removes the coprimality gaps that fragmented the naive families
(`A=3, B=8i` failed when `3 ∣ i`; `A=7, B=8i` when `7 ∣ i`). The formerly-degenerate cases
are now clean members: `k = 156 = (27, 8)`, `k = 300 = (63, 8)`.

The exact value is
```
ML(4m+3, 8, 4m+11, 4m+19) = (2m+7)/(8m+30),   D = 1/4 + 1/(16m+60).
```

* **Lower bound `ML ≥ (2m+7)/(8m+30)` (so `D ≤ 1/4 + 1/(16m+60)`): provable.** Take
  `t = p/(8m+30)` with `8m+30 = 2A+3B`, where `p ≡ (2m+7)(4m+19)^{-1} (mod 8m+30)` (the
  inverse exists since `gcd(4m+19, 8m+30) = gcd(4m+19, 8) = 1`). At this `t` the four
  distances are `(A+8)t, (A+16)t → (2m+7)/(8m+30)` (the tied active pair),
  `‖8t‖ = 2(2m+7)/(8m+30)`, `‖At‖ = (2m+9)/(8m+30)`, all `≥ (2m+7)/(8m+30)`. This is a
  clean per-runner evaluation, parametric in `m`.

* **Upper bound `ML ≤ (2m+7)/(8m+30)` (the exact realization): proved** in
  [`U2_REALIZATION_PROOF.md`](U2_REALIZATION_PROOF.md). The `1`-D maximizer sits at `t = p/q`
  with `q` one of twelve tie-moduli, so the bound is a finite per-modulus check. The single
  hard (tight, irreducibly-4-runner) case is the dominant modulus `q = 8m+30 = 4w+2`,
  `w = 2m+7`: the substitution `n = (A+8)p` sends the four runners to the **`m`-independent**
  multipliers `(3, −2, 1, −1)`, collapsing the covering to the `m`-free lemma
  `min(|n|, |2n|, |3n|)_q ≤ w` (proved by a two-interval covering, gap-free for `w ≥ 3`, and
  `w = 2m+7 ≥ 7`). The other eleven moduli are dominated with explicit positive margins. So
  `ML(v_m) = (2m+7)/(8m+30)` exactly, for every `m`.

## Status summary

| direction | statement | status |
|---|---|---|
| exclusion | `1/3, 2/7 ∉ S(U¹)∪S(U²)` | **proved** (all three exclusion notes) |
| forward, `k ≡ 4 (mod 16)` | `U¹` `(1,4j)`, `D = 1/4+1/(16j+4)` | **proved** (exact value, both bounds) |
| forward, `k ≡ 12 (mod 16)` | `U²` `(1,7)` and `(4m+3, 8)`, gap-free coprime family | **proved** (exact value, both bounds; [`U2_REALIZATION_PROOF.md`](U2_REALIZATION_PROOF.md)) |

Both directions are now proved. `2/7` and `1/3` are exceptional, and every other
`k ≡ 4 (mod 8)`, `k ≥ 20`, `k ≠ 28` is realized — by `U¹` `(1, 4j)` for `k ≡ 4 (mod 16)` and
by `U²` `(4m+3, 8)` (plus `(1,7)`) for `k ≡ 12 (mod 16)`, each with a proved exact value.
**The finite symmetric difference of Jain–Kravitz Theorem 1.3 is exactly `{1/3, 2/7}`.**
