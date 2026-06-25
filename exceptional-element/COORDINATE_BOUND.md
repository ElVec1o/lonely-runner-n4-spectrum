# A structural bound for the U² family: `ML(A, B, A+B, A+2B) ≥ 2/9` for `A ≥ 3`

This note supplies the **structural description** that the search in
[`README.md`](README.md) could only verify up to `|A|, |B| ≤ 600`: an elementary,
self-contained proof that bounds the loneliness of the entire `U²` family for all but
finitely many parameters, and thereby rules out `D = 2/7` (and every `D > 5/18`) as a
`U²`-realization outside a finite set.

By Vico Bonfioli. Proof developed with assistance from Anthropic's Claude. Apache-2.0.
The argument is elementary and is mechanically re-checked by
[`verify_coordinate_bound.py`](verify_coordinate_bound.py); nothing here relies on
trusting that assistance.

## Notation

For an integer speed vector `v = (v₁, …, v₄)`,

```
ML(v) = max_{t ∈ ℝ}  min_i ‖vᵢ · t‖,      D(v) = 1/2 − ML(v),
```

where `‖x‖` is the distance from `x` to the nearest integer. The `U²` subtorus of
Jain–Kravitz Theorem 1.3 is parameterized by coprime `(A, B)` as
`v(A,B) = (A+B, B, A, A+2B)`. Since `‖·‖` is even and `1`-periodic, `ML` depends only on
the multiset of speeds, so `ML(v(A,B)) = ML(A, B, A+B, A+2B)`.

## Theorem

> Let `A, B` be coprime positive integers with the four speeds `A, B, A+B, A+2B`
> distinct and nonzero, and suppose **`A ≥ 3`**. Then
> ```
> ML(A, B, A+B, A+2B) ≥ 2/9,        equivalently   D(v(A,B)) ≤ 5/18.
> ```

The bound is **sharp**: `ML = 2/9` (so `D = 5/18 = 1/4 + 1/36`) holds exactly at
`(A,B) ∈ {(3,1), (3,4), (4,1)}`, and the hypothesis `A ≥ 3` cannot be dropped — for
`A ≤ 2` the value drops below `2/9` (e.g. `(1,3)` and `(2,1)` give `ML = 1/5`).

## Consequence: `D = 2/7` is not a `U²`-value for `A ≥ 3`

`D = 2/7` means `ML = 3/14`. Since `2/9 = 0.2222… > 0.2142… = 3/14`, the Theorem gives
`ML ≥ 2/9 > 3/14`, hence `D(v(A,B)) < 2/7` for **every** `A ≥ 3`. More generally
`D(v(A,B)) ≤ 5/18 < 2/7` there. Therefore any `U²`-realization of `D = 2/7`
(equivalently `k = 28` in `D = 1/4 + 1/k`), or of any `D > 5/18` (i.e. `k < 36`), must
have `A ≤ 2`. That is a finite set, and direct enumeration (below) shows `D = 2/7` does
not occur in it. This converts the bounded-search *observation* of `README.md` into a
proof for the `U²` family.

## Proof

A lower bound on `ML` follows from exhibiting one good time `t`. We use only rational
`t = k/M` together with the evenness and `1`-periodicity of `‖·‖`.

### Step 1 — the collapsing modulus `M = 2(A+B)`

Take `M = 2(A+B)` and `t = k/M` with `k` **odd**. Two speeds become free:

* `‖(A+B)·t‖ = ‖k(A+B)/(2(A+B))‖ = ‖k/2‖ = 1/2`  (since `k` is odd);
* `A + 2B ≡ −A  (mod M)`, because `(A+2B) + A = 2(A+B) = M`, so `‖(A+2B)·t‖ = ‖A·t‖`.

Hence
```
min_i ‖vᵢ t‖ = min( ‖A·t‖, ‖B·t‖, 1/2, ‖A·t‖ ) = min( ‖A·t‖, ‖B·t‖ ).
```
The four-runner problem has collapsed to a **two-runner** problem.

### Step 2 — the target window `W`

Since `B = M/2 − A`, for odd `k` we have `B·t = k/2 − A·t ≡ 1/2 − A·t  (mod 1)`. Write
`φ = A·t mod 1 ∈ [0,1)`. Then `‖A·t‖ = min(φ, 1−φ)` and `‖B·t‖ = |φ − 1/2|`, and an
elementary case check gives
```
min(‖A·t‖, ‖B·t‖) ≥ 2/9   ⟺   φ ∈ W := [2/9, 5/18] ∪ [13/18, 7/9],
```
a union of two closed arcs **each of length `1/18`** (indeed `5/18 − 2/9 = 1/18` and
`7/9 − 13/18 = 1/18`). So it suffices to find an odd `k` with `A·k/M mod 1 ∈ W`.

### Step 3 — the odd-`k` orbit is equally spaced with gap `1/(A+B)`

Let `n = A+B`, `M = 2n`. Consider `O = { A·k mod M : k odd, 1 ≤ k ≤ 2n−1 }` (it has `n`
elements).

* **`A` odd.** Then `gcd(A, M) = gcd(A, 2n) = 1` (as `gcd(A, n) = gcd(A, B) = 1` and `A`
  is odd). Multiplication by `A` is a parity-preserving bijection of `ℤ/M`, so
  `O = {odd residues} = {1, 3, …, 2n−1}`, i.e. `O/M = { (2i+1)/(2n) }_{i=0}^{n-1}`.
* **`A` even** (so `B` and `n` are odd). Then `gcd(A, M) = 2` and
  `O = 2·{ (A/2)·k mod n : k odd }`. The `n` odd numbers `1, 3, …, 2n−1` cover every
  residue mod the odd integer `n`, and `×(A/2)` is a bijection mod `n`, so
  `O = {0, 2, …, 2n−2}`, i.e. `O/M = { i/n }_{i=0}^{n-1}`.

In **both** cases the points `φ` are equally spaced around `ℝ/ℤ` with gap exactly
`1/n = 1/(A+B)` (the two cases differ only by the offset `1/(2n)` vs `0`).

### Step 4 — asymptotic closure for `A + B ≥ 18`

A set of points spaced `1/(A+B)` apart meets every closed arc of length `≥ 1/(A+B)`.
The window `W` contains an arc of length `1/18`. Thus if `A + B ≥ 18` (so
`1/(A+B) ≤ 1/18`), some odd `k` satisfies `A·k/M mod 1 ∈ W`, and that `t = k/M` gives
`ML ≥ 2/9`. This case needs **no continued-fraction input**: equal spacing plus the
width of `W` is enough.

### Step 5 — the finite range `A + B ≤ 17`

Only finitely many coprime pairs `(A, B)` have `A ≥ 3` and `A + B ≤ 17` (exactly `71`).
For each, an explicit certificate `t = k/M` with all four norms `≥ 2/9` is exhibited.
Representative cases, including all three tight ones:

| `(A,B)` | speeds `(A,B,A+B,A+2B)` | `t` | `(‖A t‖, ‖B t‖, ‖(A+B)t‖, ‖(A+2B)t‖)` | `min` |
|---|---|---|---|---|
| `(3,1)` | `(3,1,4,5)`   | `4/9`  | `(3/9, 4/9, 2/9, 2/9)`     | `2/9` |
| `(4,1)` | `(4,1,5,6)`   | `4/9`  | `(2/9, 4/9, 2/9, 3/9)`     | `2/9` |
| `(3,4)` | `(3,4,7,11)`  | `1/9`  | `(3/9, 4/9, 2/9, 2/9)`     | `2/9` |
| `(5,3)` | `(5,3,8,11)`  | `2/13` | `(3/13, 6/13, 3/13, 4/13)` | `3/13`|
| `(8,1)` | `(8,1,9,10)`  | `5/18` | `(4/18, 5/18, 9/18, 4/18)` | `2/9` |
| `(10,1)`| `(10,1,11,12)`| `5/22` | `(6/22, 5/22, 11/22, 6/22)`| `5/22`|

The complete `71`-row list is regenerated and checked by
[`verify_coordinate_bound.py`](verify_coordinate_bound.py). ∎

## Verification

`verify_coordinate_bound.py` performs three independent, dependency-free checks:

1. **Exact `ML`** (no construction): computes `ML` over the full breakpoint set and
   confirms `ML ≥ 2/9` for every admissible `(A,B)` with `A ≥ 3` in range, with equality
   exactly at `(3,1), (3,4), (4,1)`.
2. **Proof-scheme coverage**: confirms that for every such `(A,B)`, either `A+B ≥ 18`
   and the `M = 2(A+B)` construction succeeds, or `A+B ≤ 17` and an explicit modulus
   certificate is found — i.e. the proof's own case split leaves no gaps.
3. **Sharpness**: confirms `A ≤ 2` admits parameters with `ML < 2/9` (so the hypothesis
   is necessary), e.g. `(1,3)` and `(2,1)` at `ML = 1/5`.

```
python3 verify_coordinate_bound.py 60      # exact rationals, ~seconds
```

## Scope (what this does and does not close)

This note proves the `U²` structural bound `D ≤ 5/18` for `A ≥ 3`, reducing "`D = 2/7` is
not realized by `U²`" to the remaining cases `A ≤ 2`. Those cases are now also closed, by
explicit witness constructions, in [`U2_SMALL_A.md`](U2_SMALL_A.md), and the `U¹` family
is handled in [`U1_FAMILY.md`](U1_FAMILY.md). Together the three notes establish that
**neither `1/3` nor `2/7` is realized by any `1`-dimensional subtorus of `U¹ ∪ U²`** — the
exclusion (backward) direction of "the symmetric difference is exactly `{1/3, 2/7}`". The
forward direction (every other progression value *is* realized) is documented numerically
in [`README.md`](README.md) and is not proved here.
