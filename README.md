# The n = 4 Lonely Runner spectrum near the accumulation point 1/4

[![DOI](https://zenodo.org/badge/1258849475.svg)](https://doi.org/10.5281/zenodo.20534835)

Two elementary finite-calculation results about the n = 4 view-obstruction / Lonely
Runner D-spectrum, on either side of its first accumulation point `1/4`, in the framework
of Jain–Kravitz, *Relative Lonely Runner spectra*
([arXiv:2411.12684](https://arxiv.org/abs/2411.12684)).

By Vico Bonfioli. Developed with substantial assistance from Anthropic's Claude.
Apache-2.0. Everything here is elementary and fully reproducible from the included code.

Notation: for an integer speed-vector `v`, `D(v) = 1/2 − max_t min_i ‖v_i t‖` where `‖x‖`
is the distance to the nearest integer; `S₁(n)` is the 1-dimensional D-spectrum,
`S₂(n)` the set of D-values of 2-dimensional subtori.

## Result A — below 1/4: `δ₂(4) = 3/14`  (`delta2-bound/`)

No saturated 2-dimensional subtorus `U ⊆ (ℝ/ℤ)⁴` has `D(U) ∈ (3/14, 1/4)`. Equivalently
the open gap below the accumulation point `1/4` contains no new (2-dimensional) value.

This is an elementary completion of Jain–Kravitz §3.3 (which classifies the two subtori
with `D(U) = 1/4`) down to the threshold `3/14`. The proof (`delta2-bound/PROOF.md`) is a
Hermite-normal-form reduction plus a 3-runner "doubling cascade" resting only on the
3-runner near-equality classification `L₃ = {(1,2,3),(1,2,6),(1,3,4),(1,5,6),(2,3,5)}`;
it uses no continuum/`Prog` machinery. The load-bearing finite check is verified with
exact rational arithmetic (`delta2-bound/rust/`, `cargo run --release --bin
delta2_hnf_sweep 9`): over all `(1,2)`-pivot saturated planes the only exact D-values
`≥ 3/14` are `3/14`, `1/4`, `1/2`, with nothing in the open gap (corroborated exactly out
to `|entries| ≤ 70`, ~4·10⁸ planes).

## Result B — above 1/4: `D = 2/7` is a candidate exceptional element  (`exceptional-element/`)

In `S₁(4) ∩ (1/4, 1/2]` (Jain–Kravitz Theorem 1.3, whose finite symmetric difference is
left uncomputed), the value `D = 2/7` appears to be **omitted**: over `|A|,|B| ≤ 600` the
realized `D = 1/4 + 1/k` with `k ∈ [20,84]` are exactly `{k ≡ 4 mod 8, k ≥ 20}` except
`k = 28`. Run `python3 exceptional-element/spectrum.py 60` (exact, no dependencies).

This is an *observation* (a question for the structure theory), not a proof.

## The two results meet at 3/14

`D = 2/7 ⟺ ML = 3/14`, the same `3/14` that bounds Result A from below. The gap just
below `1/4` is empty (A), and the value just above `1/4` whose loneliness equals `3/14`
is the one omitted from the relative spectrum (B).

## A correctness pitfall (worth flagging)

When searching for the optimal `t = k/d` in `max_t min_i ‖v_i t‖`, one must check **all**
`k`, not only `gcd(k,d) = 1`: an optimal reduced `t = p/q` can have `q` outside the
breakpoint set `{|v_i|, |v_i ± v_j|, 2|v_i|}` yet `q | d` for some `d` in it (e.g.
`1/4 = 5/20`, `20 = 7+13`). The coprime shortcut silently overestimates `D`. The verifiers
here avoid it, and results were cross-checked against independent floating-point grids.

## Citation

If you refer to this work, please cite it via its Zenodo DOI:

> Vico Bonfioli. *The n=4 Lonely Runner spectrum near 1/4: δ₂(4)=3/14 and a candidate
> exceptional element D=2/7*. Zenodo. <https://doi.org/10.5281/zenodo.20534835>

`10.5281/zenodo.20534835` is the concept DOI (always resolves to the latest version); the
`v1.0.0` release is archived at `10.5281/zenodo.20534834`.

## Acknowledgements

Developed by Vico Bonfioli with substantial assistance from Anthropic's Claude for
derivation and computational verification. The arithmetic is elementary and reproducible;
nothing rests on trusting that assistance.
