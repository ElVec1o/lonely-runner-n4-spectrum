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
