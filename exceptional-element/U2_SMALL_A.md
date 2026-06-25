# `U²` realizes neither `1/3` nor `2/7` — uniformly in `A`

[`COORDINATE_BOUND.md`](COORDINATE_BOUND.md) proves the sharp bound `ML ≥ 2/9`
(`D ≤ 5/18`) for `A ≥ 3`, which already excludes `2/7` and `1/3` there. This note closes
the remaining cases `A ∈ {1, 2}` — and in fact gives a single uniform argument valid for
**all** `A ≥ 1`, by reusing the *same* construction with a wider target band.

By Vico Bonfioli. Proof developed with assistance from Anthropic's Claude. Apache-2.0.
Mechanically re-checked by [`verify_u2_smallA.py`](verify_u2_smallA.py).

## The construction, at a general band

Recall from [`COORDINATE_BOUND.md`](COORDINATE_BOUND.md) the modulus `M = 2(A+B)` with
`t = k/M`, `k` odd: runner `A+B` sits at `1/2`, runner `A+2B ≡ −A (mod M)` ties runner
`A`, and the four-runner minimum collapses to `min(‖A·t‖, ‖B·t‖)`. Writing
`φ = A·t mod 1`, for odd `k` one has `‖B·t‖ = |φ − 1/2|`, so for any threshold `c < 1/4`,

```
min(‖A·t‖, ‖B·t‖) ≥ c   ⟺   φ ∈ [c, 1/2 − c] ∪ [1/2 + c, 1 − c],
```

two arcs **each of width `1/2 − 2c`**. The odd-`k` orbit `{φ}` is equally spaced with gap
`1/(A+B)` (proved, both parities, in `COORDINATE_BOUND.md`). So a gap strictly smaller
than the arc width forces a point in the *open* arc, giving a strict inequality:

```
A + B > 1/(1/2 − 2c)   ⟹   ML(A,B,A+B,A+2B) > c.
```

(With `c = 2/9` this is the `A+B ≥ 18` bound of the coordinate note; here we use larger
`c`-gaps, hence smaller thresholds.)

## Excluding `2/7`

`D = 2/7 ⟺ ML = 3/14`. Take `c = 3/14`: the arc width is `1/2 − 3/7 = 1/14`, so

```
A + B ≥ 15   ⟹   ML > 3/14   ⟹   D < 2/7.
```

For the finite remainder `A + B ≤ 14` (finitely many coprime pairs), an exact computation
finds **no** `(A,B)` with `D = 2/7`. Hence `D = 2/7` is realized by `U²` for no `A, B`.

## Excluding `1/3`

`D = 1/3 ⟺ ML = 1/6`. Take `c = 1/6`: the arc width is `1/2 − 1/3 = 1/6`, so

```
A + B ≥ 7   ⟹   ML > 1/6   ⟹   D < 1/3.
```

The finite remainder `A + B ≤ 6` contains no `(A,B)` with `D = 1/3`. Hence `D = 1/3` is
realized by `U²` for no `A, B`.

## Conclusion (with `U¹`)

* `U²`: neither `1/3` nor `2/7` (this note; uniform in `A`).
* `U¹`: realizes exactly `D = 1/4 + 1/(16j+4)` (`k ≡ 4 mod 16`), containing neither
  ([`U1_FAMILY.md`](U1_FAMILY.md)).

**`D = 1/3` and `D = 2/7` are realized by no `1`-dimensional subtorus of `U¹ ∪ U²`.** This
is the exclusion (backward) direction of "the finite symmetric difference of Jain–Kravitz
Theorem 1.3 is exactly `{1/3, 2/7}`". The forward direction (every other `k ≡ 4 (mod 8)`,
`k ≥ 20`, `k ≠ 28` *is* realized) is the numerical observation in [`README.md`](README.md).

## Appendix: the exact `U²` value laws (for reference)

The in-spectrum `U²` deficits obey explicit linear laws (verified to `B ≤ 200` by exact
`ML`; certified below by pair-sum witnesses):

```
A = 1:  B ≡ 1,2 (mod 4): D ≤ 1/4;   B ≡ 0 (4): D = 1/4 + 1/(12B+4);   B ≡ 3 (4): D = 1/4 + 1/(6B+2)
A = 2:  B ≡ 1 (mod 4):   D = 1/4 + 1/(12B+8);   B ≡ 3 (mod 4): D = 1/4 + 1/(12B+16)
```

These are not needed for the exclusion above, but they characterize *which* values the
small-`A` families realize. The `A=1, B ≡ 0 (mod 4)` law has a clean closed-form witness:
with `B = 4c`, `t = (3c+1)/(12c+1)` gives runner distances
`(3c+1, 3c, 6c, 3c)/(12c+1)`, so `ML = 3c/(12c+1) = 3B/(4(3B+1))` exactly.
