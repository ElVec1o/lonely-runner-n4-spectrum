# δ₂(4) = 3/14 via a Hermite-normal-form coordinate bound

Author: Vico Bonfioli, with substantial assistance from Anthropic's Claude for derivation and
computational verification. The argument is elementary and the load-bearing finite check is verified
with exact rational arithmetic (`rust/`).

## 1. Statement

`D(T) = 1/2 − ML(T)`, `ML(T) = max_t min_{i : v_i ≠ 0} ‖v_i t‖` (nearest-integer norm, nonzero-speed
convention). For a 2-dimensional subtorus `U ⊆ (ℝ/ℤ)⁴`,
`D(U) = 1/2 − sup_{(s,t)} min_{i : L_i ≠ 0} ‖L_i(s,t)‖`.

> **Theorem.** No 2-dimensional **saturated** subtorus `U ⊆ (ℝ/ℤ)⁴` has `D(U) ∈ (3/14, 1/4)`.

"Saturated" means `⟨u,v⟩_ℝ ∩ ℤ⁴ = ⟨u,v⟩_ℤ`, equivalently the integer `2×4` generator matrix is
primitive (its six `2×2` minors have gcd 1). Saturation is essential: `D(U)` depends only on the
rational 2-plane, so a non-primitive presentation (e.g. `15·(1,1,0,2)`) is a large-coordinate basis
for a small subtorus and must be reduced first. (Note `3/14 ∈ S₁(4)` already, from `{1,2,6}`; so only
the open-gap phrasing is meaningful, and "δ₂(4)=3/14" means this gap contains no new, i.e.
saturated 2-dimensional, value.)

## 2. Reduction to Hermite normal form

Every saturated rank-2 sublattice of `ℤ⁴` has a **unique** row-Hermite-normal-form basis with pivots
in some coordinate pair `(i,j)`. Since `D` is invariant under coordinate permutation and negation,
and every rank-2 plane projects isomorphically onto *some* coordinate pair (some `2×2` minor ≠ 0),
**WLOG the pivots are `(1,2)`**:

> `U = ⟨(1,0,x,y), (0,1,z,w)⟩`,  `x,y,z,w ∈ ℤ`,  coordinate forms `L = (s, t, xs+zt, ys+wt)`.

This parameterization is exact and non-redundant (no scaling copies); it captures
`U¹ = ⟨(0,1,2,3),(1,0,0,0)⟩` as `(x,y,z,w) = (0,0,2,3)` with `D = 1/4`.

## 3. The structural facts (exact)

All values below are from the **exact** sweep of §5 (no floating point).

**(a) The only saturated `D`-values `≥ 3/14` are `3/14`, `1/4`, `1/2`.** Exhaustive exact over all
`(1,2)`-pivot planes with `|entries| ≤ 70` (3.95·10⁸ planes): every exact `D ≥ 3/14` equals one of
these three. Both the open gap `(3/14,1/4)` **and** the band `(1/4,1/2)` are empty.

**(b) `D = 1/2` is the only unbounded value.** `D(U)=1/2 ⟺ ML(U)=0`, which (in the verifier's
min-over-all-four-forms convention) is exactly `x=z=0` or `y=w=0`, i.e. a coordinate form vanishes
identically — a slice (§4a). These occur at all sizes but have `D ≥ 1/4`, outside the gap. By
contrast `{D=1/4}` (**56 planes, `|entries| ≤ 3`**) and `{D=3/14}` (**144 planes, `|entries| ≤ 6`**)
are **closed finite sets** — counts and max-coords are identical from `B=16` through `B=70`, hence
complete. The `≤6` bound matches `max L₃ = 6` (the `(1,2,6)` triple), exactly as the §4 cascade
predicts.

## 4. The coordinate bound (the engine) — 3-runner doubling cascade

> **Lemma (coordinate bound).** If `U = ⟨(1,0,x,y),(0,1,z,w)⟩` is saturated with `3/14 ≤ D(U) < 1/2`,
> then `|x|,|y|,|z|,|w| ≤ 9`.

(The exact sweep shows the true bound is `6 = max L₃`; `9` is the clean bound the cascade gives in the
degenerate corners, and is all the finite check needs.)

*Mechanism.* Consider the two **doubled-coordinate** integer directions in `U`:

- `d₊ = (1,1, x+z, y+w)` — speeds `(1,1,x+z,y+w)`, i.e. the 3-runner shadow `{1, |x+z|, |y+w|}`;
- `d₋ = (1,−1, x−z, y−w)` — shadow `{1, |x−z|, |y−w|}`.

If `d₊` is **zero-free** (`x+z ≠ 0`, `y+w ≠ 0`) then `D(U) ≤ D(⟨d₊⟩) = D(\{1,|x+z|,|y+w|\})`, a pure
**3-runner** value. By the 3-runner classification, a 3-runner with `D ≥ 3/14` has primitive triple
in `L₃ = {(1,2,3),(1,2,6),(1,3,4),(1,5,6),(2,3,5)}` (independently due to Y-G Chen; a Lean
formalization of this input is at `github.com/ElVec1o/kravitz-lonely-runner-n3`). Since
`{1,|x+z|,|y+w|}` contains a `1`, it is one of the four L₃ triples containing `1`, all with entries
`≤ 6`. Hence `|x+z|,|y+w| ≤ 6`; symmetrically (if `d₋` zero-free) `|x−z|,|y−w| ≤ 6`. From
`x = ((x+z)+(x−z))/2` etc., `|x|,|y|,|z|,|w| ≤ 6`.

The cascade travels **only through doubled (3-runner) directions**, where the classification applies.
This is essential: a cascade through generic 4-runner directions would not work, because `S₁(4)`
genuinely meets `(3/14,1/4)` (e.g. `D(1,2,4,9) = 5/22`) — the 1-dimensional gap is a 3-runner fact
that does not extend to 4 runners. Doubling sidesteps this entirely.

### 4a. Slices (a coordinate form is identically zero)

If `coord-3 form ≡ 0` (i.e. `x = z = 0`) then `U = ⟨(1,0,0,y),(0,1,0,w)⟩` is a 2-torus of the
3-coordinate subspace `{x₃=0} ≅ (ℝ/ℤ)³`, namely `⟨(1,0,y),(0,1,w)⟩`. The lattice direction
`(1,1,y+w)` (pivots doubled, `A=B=1`) has shadow `{1,|y+w|}` — a **2-runner**, so
`D(U) ≤ D(⟨(1,1,y+w)⟩) = 1/(2(1+|y+w|)) ≤ 1/6` (use `(1,−1,2y)` if `y+w=0`, and `D=0` if also
`y=0`). Hence `D(U) ≤ 1/6 < 3/14`: **slices never enter the gap.** Same for `coord-4 form ≡ 0`
(`y=w=0`). These are the planes the verifier reports as `D=1/2` under its min-over-all-four-forms
convention; the δ₂-relevant value is the `≤ 1/6` above. For all other `U`, both forms `xs+zt` and
`ys+wt` are `≢ 0`, so `U` has a zero-free direction and `D(U) < 1/2`.

### 4b. Degenerate cases of the cascade (non-slice; all closed analytically)

Write `P=x+z, Q=y+w, R=x−z, S=y−w` (so `x=(P+R)/2` etc.). The generic case is `P,Q,R,S ≠ 0`
(both `d₊,d₋` zero-free). The coordinate permutations/negations that fix the HNF form act on
`(P,Q,R,S)` as the Klein four-group `⟨(P R)(Q S),\,(P Q)(R S)⟩`, which is **transitive** on
`{P,Q,R,S}` and identifies the pair-orbits `{P,R}~{Q,S}`, `{P,Q}~{R,S}`, `{P,S}~{Q,R}`. Non-slice
means `(x,z)≠(0,0)` and `(y,w)≠(0,0)`. A vanishing of three or more of `P,Q,R,S`, or of the
"diagonal" pair `{P,R}` (`=x=z=0`) or `{Q,S}` (`=y=w=0`), forces a slice (§4a). So up to symmetry
the only non-slice degenerate cases are:

- **(I) `P=0`** (`z=−x`, `x≠0`, and `Q,R,S≠0`). `d₋=(1,−1,2x,y−w)` is zero-free →
  `{1,2|x|,|y−w|}∈L₃` → `|x|≤3` and `|y−w|≤6`. For `Q=y+w`: the `c₃=c₄` direction is
  `(w+x,\,x−y,\,xQ,\,xQ)` (using `wx−yz = x(w+y) = xQ`). If zero-free, `{|w+x|,|x−y|,|xQ|}∈L₃`,
  so `|xQ|≤6`, and as `|x|≥1` we get `|Q|≤6`, hence `|y|,|w|≤6`. The only way it is *not* zero-free
  is `w+x=0` or `x−y=0` (as `xQ≠0`): then `|w|=|x|≤3` (resp. `|y|=|x|≤3`), and with `|y−w|≤6` the
  remaining one is `≤9`. So **`|entries| ≤ 9`**.
- **(II) `P=Q=0`** (`z=−x, w=−y`, `x,y≠0` for non-slice). `d₋=(1,−1,2x,2y)` is zero-free with shadow
  `{1,2|x|,2|y|}` — all-even apart from the `1`, so the only L₃ triple it can be is `(1,2,6)` →
  `{|x|,|y|}={1,3}` → **`|entries| ≤ 3`**.
- **(III) `P=S=0`** (`z=−x, w=y`, `x,y≠0`). Both `d₊,d₋` degenerate; use `c₃=c₄`,
  `(x+y,\,x−y,\,2xy,\,2xy)` (since `wx−yz = yx+yx = 2xy`). If zero-free (`|x|≠|y|`),
  `{|x+y|,|x−y|,2|xy|}∈L₃` → `|xy|≤3` → `|x|,|y|≤3`. If `|x|=|y|`, the `c₁=c₃` direction
  `(−x,1−x,−x,x−2x²)` has shadow `{|x|,|x−1|,|x||2x−1|}`; its largest entry `|x|(2|x|−1)` exceeds
  `6` once `|x|≥3`, forcing `D<3/14` — so `|x|≤2`. Either way **`|entries| ≤ 3`**.

So in **every** non-slice case `3/14 ≤ D(U) < 1/2` gives `|x|,|y|,|z|,|w| ≤ 9`. Every step is a
3-runner cascade resting only on the L₃ classification — no floating point, no continuum.

## 5. Exact exhaustive verification and conclusion

All computation is **exact rational** (`rust/`, `d_2dim` = tie-line vertex enumeration; no floating
point). The sound `B1+` cheap upper bound proves `D < 3/14` for a plane when some small zero-free
combination `a·u+b·v` has 1-dim `D < 3/14`.

> **Sweep result (`delta2_hnf_sweep`).** Over **all** `(1,2)`-pivot saturated planes with
> `|x|,|y|,|z|,|w| ≤ B`:
>
> | `B` | planes | `D ∈ (3/14,1/4)` | `D=1/4` (count, maxcoord) | `D=3/14` (count, maxcoord) | `D=1/2` |
> |----|--------|----|----|----|----|
> | 16 | 1.19·10⁶ | **0** | 56, ≤3 | 144, ≤6 | unbounded |
> | 40 | 4.30·10⁷ | **0** | 56, ≤3 | 144, ≤6 | unbounded |
> | 70 | 3.95·10⁸ | **0** | 56, ≤3 | 144, ≤6 | unbounded |
>
> The only exact `D`-values `≥ 3/14` are `3/14`, `1/4`, `1/2`. Nothing in `(3/14,1/4)` **or**
> `(1/4,1/2)`. The `D=1/4` and `D=3/14` sets are closed and bounded; only `D=1/2` grows with `B`.

*Conclusion.* By §2, WLOG `U = ⟨(1,0,x,y),(0,1,z,w)⟩`. Suppose `D(U) ∈ (3/14,1/4)`. By §4a, `U` is
not a slice (slices have `D ≤ 1/6`); so `D(U) < 1/2`, and `3/14 ≤ D(U) < 1/2`. By the Lemma (§4
generic + §4b degenerate), `|x|,|y|,|z|,|w| ≤ 9` — a finite set of `19⁴ = 130321` planes, each
checked **exactly**: every value is `≤ 3/14` or `1/4` (the verifier's `D=1/2` entries are the §4a
slices, with δ₂-value `≤ 1/6`), **none in `(3/14,1/4)`**. Contradiction. Hence no saturated 2-dim
`U ⊆ (ℝ/ℤ)⁴` has `D(U) ∈ (3/14,1/4)`; equivalently `δ₂(4) = 3/14`. ∎

**Status of rigor.**
- §2 reduction (HNF + coordinate symmetry): rigorous.
- §4a slices: rigorous (`D ≤ 1/6` via a doubled-pivot 2-runner direction).
- §4 generic case + §4b degenerate cases (I,II,III): rigorous, each a 3-runner cascade resting only
  on the L₃ classification; the symmetry reduction shows these exhaust all non-slice degeneracies.
- §5 finite check (`|entries| ≤ 9`, 130321 planes) and the corroborating sweep to `B=70`
  (3.95·10⁸ planes): exact, no floating point.

This is an elementary completion of Jain–Kravitz §3.3 (which classifies `D(U)=1/4`) down to the
threshold `3/14` that defines `δ₂(4)`, with HNF and a 3-runner doubling cascade making the bound
concrete (no continuum / `Prog` machinery), and the bound `6 = max L₃` falling out naturally.

## 6. Remarks: three pitfalls the argument must avoid

These are recorded because each is a natural-looking approach that fails, and seeing why clarifies
where the difficulty lies.

1. **The cascade is valid only for zero-free directions.** `D(U) ≤ D(⟨w⟩)` can fail when `w` has a
   zero coordinate, because the nonzero-speed convention drops that coordinate and breaks the
   monotonicity (e.g. `u=(1,1,3,2), v=(1,-1,1,0)`). Every cascade above uses a zero-free `w`.
2. **The 1-dimensional spectral gap is a 3-runner fact, not a 4-runner one.** It is tempting to argue
   `D(U) ∈ closure(S₁(4))` and invoke a gap, but `S₁(4) ∩ (3/14,1/4) ≠ ∅` (752 four-runners lie
   there, e.g. `D(1,2,4,9)=5/22`, `D(1,2,6,17)=2/9`). The 2-dimensional infimum avoids the gap, but
   that is not a 1-dimensional spectral statement; this is exactly why the cascade stays in 3-runner
   territory.
3. **Saturation is required.** Without it, the value `1/4` and the band `(1/4,1/2)` appear at
   arbitrarily large coordinates — but only as non-primitive scalings of small subtori. The HNF of §2
   removes this redundancy and is what makes the coordinate bound finite.
