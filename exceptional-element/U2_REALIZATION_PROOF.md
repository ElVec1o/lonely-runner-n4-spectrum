# The `U²` realization exact value: `ML(4m+3, 8, 4m+11, 4m+19) = (2m+7)/(8m+30)`

This proves the exact value for the family that realizes `k ≡ 12 (mod 16)` (see
[`REALIZATION.md`](REALIZATION.md)), closing the forward direction. With `A = 4m+3` and
`v_m = (A, 8, A+8, A+16)`,

```
ML(v_m) = (2m+7)/(8m+30),   so   D = 1/2 − ML = 1/4 + 1/(16m+60),   m ≥ 0.
```

The lower bound `ML ≥ (2m+7)/(8m+30)` is the explicit witness of [`REALIZATION.md`](REALIZATION.md).
This note proves the matching **upper bound** `ML ≤ (2m+7)/(8m+30)`, which is the part that
makes the realization exact (not merely `D ≤ …`). It is checked at every step by
[`verify_u2_realization.py`](verify_u2_realization.py).

By Vico Bonfioli. Proof developed with assistance from Anthropic's Claude. Apache-2.0.

## Reduction to twelve tie-moduli

`t ↦ min_i ‖v_i t‖` is continuous, `1`-periodic and piecewise linear, with breakpoints only
where two runners tie (`‖v_i t‖ = ‖v_j t‖`) or one vanishes (`‖v_i t‖ = 0`); both happen only
at `t = p/q` with `q` in the tie-modulus set `{2|v_i|} ∪ {|v_i ± v_j|}`. So
`ML(v_m) = max_q R(q,m)` over

```
q ∈ {8, 16, A, A+8, A+16, |A−8|, A+24, 2A, 2A+8, 2A+16, 2A+24, 2A+32},
R(q,m) := max_{1 ≤ p < q} min_i ‖v_i p/q‖.
```

It suffices to show `R(q,m) ≤ T(m) := (2m+7)/(8m+30)` for each `q`, with equality only at the
**dominant** modulus `q = 2A+24 = 8m+30`.

## The covering lemma

The engine is a single `m`-free statement. Write `|x|_q` for the distance from `x` to the
nearest multiple of `q`, and `M(q) := max_n min(|n|_q, |2n|_q, |3n|_q)`.

> **Covering Lemma.** For every `q ≥ 6`, `M(q) = ⌊q/4⌋`.

**Proof.** `min(|n|,|2n|,|3n|)_q ≤ |n|_q`, and `|n|_q = ⌊q/4⌋` is attained, so `M(q) ≥ ⌊q/4⌋`;
we show `M(q) ≤ ⌊q/4⌋`. Put `D = ⌊q/4⌋` and take `r = |n|_q ∈ [0, ⌊q/2⌋]`. If `r ≤ D`, runner
`1` covers. Otherwise `D < r ≤ ⌊q/2⌋`:
* Runner `2`: `2r ∈ (q/2, q]`, so `|2n|_q = q − 2r ≤ D ⟺ r ≥ (q−D)/2`.
* Runner `3`: `|3n|_q ≤ D ⟺ 3r ∈ [q−D, q+D] ⟺ r ∈ [(q−D)/3, (q+D)/3]`, and `(q−D)/3 ≤ D+1`
  (since `q ≤ 4D+3`), so runner `3` covers `[D+1, ⌊(q+D)/3⌋]`.

The two cover `(D, ⌊q/2⌋]` with no gap because `⌈(q−D)/2⌉ ≤ ⌊(q+D)/3⌋ + 1`: using
`D ≥ (q−3)/4`, one has `(q−D)/2 ≤ (3q+3)/8` and `(q+D)/3 ≥ (5q−3)/12`, and
`(3q+3)/8 ≤ (5q−3)/12 ⟺ 9q+9 ≤ 10q−6 ⟺ q ≥ 15`. So the bands cover for all `q ≥ 15`; the
finitely many `6 ≤ q ≤ 14` are checked directly. Hence `min(|n|,|2n|,|3n|)_q ≤ ⌊q/4⌋` for all
`n`. ∎

## The dominant and the other non-degenerate moduli

A modulus is **degenerate** if some runner equals `q/2`. The non-degenerate ones are
`q ∈ {|A−8|, A+24, 2A+8, 2A+24}`. For each, `gcd(A+8, q) = 1`, and the substitution
`n = (A+8)p` sends the four runner residues to a multiset that is `{1,2,3}` up to sign and
repetition — for the dominant `q = 2A+24 = 8m+30`, using `2(A+8) = q−8 ≡ −8`, the multipliers
are exactly `(3, −2, 1, −1)`; for `|A−8|, A+24, 2A+8` they are `(1,1,2,3)`. Since `|−x|_q = |x|_q`,

```
R(q,m) = M(q)/q = ⌊q/4⌋ / q.
```

The dominant modulus is the largest of the four, and `⌊q/4⌋/q` there equals
`⌊(8m+30)/4⌋/(8m+30) = (2m+7)/(8m+30) = T(m)` — and `T(m)` is the **largest multiple of `1/q`
below `1/4`** (the next, `(2m+8)/q`, exceeds `1/4` since `8m+32 > q`). For the three smaller
non-degenerate moduli, `⌊q/4⌋/q < T(m)` by direct comparison (margins `4m+25, 4m+9, 8`
respectively). So among the non-degenerate moduli the maximum is exactly `T(m)`, attained only
at `q = 8m+30`.

## The degenerate and trivial moduli

* `q ∈ {8, A, A+8, A+16}` is a runner speed, so that runner is `≡ 0`: `R = 0`.
* `q = 16`: `R = 3/16 < T(m)` (margin `8m+22`).
* `q ∈ {2A, 2A+32} = {2A, 2(A+16)}`: the runner `A` resp. `A+16` is `q/2` (`∈ {0, ½}`); the other
  three reduce to `g(x) = min(‖x‖, ‖x+½‖, ‖2x+½‖)`, whose maximum `1/6` is attained on the
  `(1/24)ℤ` grid (all kinks and pairwise crossings lie there). So `R ≤ 1/6 < T(m)` (margin `4m+12`).
* `q = 2A+16 = 2(A+8)`: the runner `A+8` is `q/2`; the remaining three give `R = (2m+5)/(8m+22)
  = ⌊q/4⌋/q < T(m)` (margin `4`).

The small base cases `m = 0, 1, 2` (where `|A−8| ∈ {5,1,3}` is irregular) are checked directly:
`ML = T` holds. All closed forms and margins are confirmed exactly in the script.

## Conclusion

`ML(v_m) = max_q R(q,m) = T(m) = (2m+7)/(8m+30)` for every `m ≥ 0`. Hence
`D(U²(4m+3, 8)) = 1/4 + 1/(16m+60)` **exactly**, so `U²` realizes every `k ≡ 12 (mod 16)`,
`k ≥ 60`; together with `k = 44` at `(1,7)` and the `U¹` family for `k ≡ 4 (mod 16)`, every
`k ≡ 4 (mod 8)`, `k ≥ 20`, `k ≠ 28` is realized. Combined with the exclusion notes, the finite
symmetric difference of Jain–Kravitz Theorem 1.3 is **exactly `{1/3, 2/7}`**.
