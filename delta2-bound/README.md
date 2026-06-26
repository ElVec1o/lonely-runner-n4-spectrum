# δ₂(4) = 3/14 :  S₂(4) ∩ (3/14, 1/4) = ∅

No saturated 2-dimensional subtorus of `(ℝ/ℤ)⁴` has D-value strictly between `3/14` and
`1/4`.

- **`PROOF.md`** — the proof. Hermite-normal-form reduction of saturated 2-planes to
  `U = ⟨(1,0,x,y),(0,1,z,w)⟩`; slices give `D ≤ 1/6`; a 3-runner doubling cascade bounds
  the non-slice case to `|entries| ≤ 9`; an exact finite check over `|entries| ≤ 9` closes
  the gap. The cascade rests only on the 3-runner classification `L₃` (independently due to
  Y-G Chen; see also the Lean formalization at
  `github.com/ElVec1o/kravitz-lonely-runner-n3`). The note closes with three remarks on
  natural-looking approaches that fail, and why.

- **`rust/`** — exact verification (no floating point). `d_2dim` computes `D(U)` exactly by
  tie-line vertex enumeration; `delta2_hnf_sweep` runs it over all `(1,2)`-pivot saturated
  planes with `|entries| ≤ B`, screening with the sound `B1+` upper bound.

## Run

```
cd rust
cargo run --release --bin delta2_hnf_sweep 9     # the proof's finite check (|entries| ≤ 9)
cargo run --release --bin delta2_hnf_sweep 40    # corroboration to |entries| ≤ 40
```

Expected (`B = 9`): `EXACT D in open gap (3/14, 1/4): 0`, and the only exact D-values
`≥ 3/14` are `1/4` (56 planes, coords ≤ 3), `3/14` (144 planes, coords ≤ 6), and `1/2`
(the coordinate-degenerate slices). The `1/4` and `3/14` sets are closed and bounded; only
`1/2` grows with `B`.

`d_2dim` was cross-checked against an independent fine-grid evaluation (0 disagreements on
the gap-relevant planes); the earlier float grid is used nowhere in the proof.

## Independent `WallCore` verifier (`verify_wallcore.py`)

A second, from-scratch implementation (dependency-free exact integer arithmetic, no shared
code with `rust/`) that re-derives the *combinatorial core* the formalization reduces to —
not the maximizer, but the **good-wall** property. The Lean reduction
(`delta2sat_le_of_wallCoreComb`) shows `δ₂(4) ≤ 3/14` follows from: every saturated band has a
wall `(i,j,ε)` whose merged minor-triple is full-support and not `(1,2,3)`-shaped; and
`jkForm_bounded` bounds the normal form to `c,f ≤ 30`. So the bound reduces to a **pure
integer check** — no real-valued `D` needed for the generic family.

```
python3 verify_wallcore.py 30      # the proven bound
```

Result: over **10,045,542** saturated Case-A normal forms (`a<b<c`, `d<e<f ≤ 30`), **every one
has a good wall** (`0` exceptions) ⟹ all are `≤ 3/14`, nothing in the open gap. The script's
guarded `mD` spot-check shows *why* this is sound where a pure-minor abstraction is not: the
boundary form `(1,2,3,1,2,3)` (`mD = 1/4`) is excluded as **non-saturated** (minor-gcd `2`, the
`2·Lrz` scaling). The degenerate special cases (coincident/zero parameters) are bounded
separately (Lean `NormalForm`: caseB/C/D). This independently confirms the WallCore residual
across the full bounded family.

## Active-wall probes (`probe_primitive_wall.py`, `diag_primitive_wall.py`, `probe_band_primitive.py`)

These mine the structure of the *other* clean residual the Lean library exposes,
`PrimitiveActiveWall`: every saturated band with `mD < 1/4` has a maximizer `τ` and a
**primitive** active wall `gcd(A,B) = 1`, `(A,B) = M_i − ε·M_j`. Each probe computes the exact
loneliness maximizer (rational vertex of the wall-line arrangement, exact `Fraction` arithmetic)
and inspects the active walls there.

```
python3 probe_primitive_wall.py 12     # re-derive + mine the gcd / maximizer structure
python3 diag_primitive_wall.py         # diagnose the flagged bands
python3 probe_band_primitive.py 14     # restrict to the relevant band mgap ∈ (1/4, 2/7)
```

Finding: the maximizer active wall is primitive on **most** bands but **not all** — it fails on
the *deep-lonely* regime `mgap > 2/7` (`mD < 3/14`). Concrete refuter: `jk(1,2,5,3,4,6)` is
saturated (six minors `[-6,-2,-9,10,21,-8]`, gcd `1`), has `mgap = 5/12`, and carries **no**
primitive integer wall at any maximizer. So the unrestricted "primitive active wall at the deep
hole" claim over-reaches: those deep bands are exactly the ones where the `δ₂(4)` bound is already
trivial (`mD = 1/2 − mgap ≤ 1/2 − 2/7 = 3/14`). The wall argument is only *needed* on the gap band
`mgap ∈ (1/4, 2/7)` ⟺ `mD ∈ (3/14, 1/4)` — and that band is **empty** of saturated forms
(`0` over entries `≤ 14`, matching the Rust sweep), which is exactly `S₂(4) ∩ (3/14, 1/4) = ∅`.
The honest minimal residual is therefore the band-restricted gap-emptiness, not a primitivity
claim over all `mD < 1/4` (formalized as `BandGapEmpty` in the Lean library).
