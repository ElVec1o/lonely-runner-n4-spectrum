//! Exhaustive EXACT sweep over saturated 2-dim subtori of (R/Z)^4 in
//! Hermite normal form with pivots in coords (1,2):
//!
//!     U = < (1,0,x,y), (0,1,z,w) >_R,   |x|,|y|,|z|,|w| <= B.
//!
//! Every saturated rank-2 sublattice of Z^4 is, up to a coordinate
//! permutation/sign (which preserve D), of this form for some pivot
//! pair; D is symmetric under those, so this family covers all
//! saturated 2-tori up to symmetry (those whose projection to coords
//! (1,2) is an isomorphism — i.e. some 2x2 minor is +-1, which holds
//! for a primitive lattice after the right coordinate choice).
//!
//! For each plane we either prove D(U) < 3/14 by the sound B1+ cheap
//! upper bound (a small zero-free integer direction a*u+b*v with
//! D(<a u + b v>) < 3/14, valid since that 1-torus sits inside U), or
//! compute D(U) EXACTLY via d_2dim (tie-line vertex enumeration; no
//! floating point, so no grid artifacts).
//!
//! Reports: any plane with EXACT D in the open gap (3/14, 1/4) (must be
//! none), and the full distribution of exact D-values >= 3/14 with the
//! max coordinate at which each occurs (to expose which values are
//! coordinate-bounded and which escape to large coordinates).
//!
//! Usage: delta2_hnf_sweep <B>     (default 16)

use rayon::prelude::*;
use std::collections::HashMap;
use std::env;

use delta2_verify::d_1dim::b1_plus_below_threshold;
use delta2_verify::d_2dim::d_2dim;

fn main() {
    let b: i64 = env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(16);
    let (gn, gd) = (3i64, 14i64); // 3/14
    let (qn, qd) = (1i64, 4i64); // 1/4

    // Parallel over the outer coordinate x.
    let per_x: Vec<(Vec<(i64, i64, i64, i64)>, HashMap<(i64, i64), (u64, i64)>, u64)> = (-b..=b)
        .into_par_iter()
        .map(|x| {
            let mut gap: Vec<(i64, i64, i64, i64)> = Vec::new();
            let mut dist: HashMap<(i64, i64), (u64, i64)> = HashMap::new();
            let mut filtered: u64 = 0;
            for y in -b..=b {
                for z in -b..=b {
                    for w in -b..=b {
                        let u = [1, 0, x, y];
                        let v = [0, 1, z, w];
                        if b1_plus_below_threshold(&u, &v, gn, gd) {
                            filtered += 1;
                            continue;
                        }
                        let (n, d) = d_2dim(&u, &v);
                        // record only D >= 3/14 (n/d >= 3/14  <=>  n*14 >= 3*d)
                        if (n as i128) * (gd as i128) >= (gn as i128) * (d as i128) {
                            let mc = x.abs().max(y.abs()).max(z.abs()).max(w.abs());
                            let e = dist.entry((n, d)).or_insert((0, 0));
                            e.0 += 1;
                            if mc > e.1 {
                                e.1 = mc;
                            }
                            // open gap: 3/14 < n/d < 1/4
                            let gt = (n as i128) * (gd as i128) > (gn as i128) * (d as i128);
                            let lt = (n as i128) * (qd as i128) < (qn as i128) * (d as i128);
                            if gt && lt {
                                gap.push((x, y, z, w));
                            }
                        } else {
                            filtered += 1;
                        }
                    }
                }
            }
            (gap, dist, filtered)
        })
        .collect();

    // Reduce.
    let mut gap: Vec<(i64, i64, i64, i64)> = Vec::new();
    let mut dist: HashMap<(i64, i64), (u64, i64)> = HashMap::new();
    let mut filtered: u64 = 0;
    for (g, d, f) in per_x {
        gap.extend(g);
        filtered += f;
        for (k, (c, mc)) in d {
            let e = dist.entry(k).or_insert((0, 0));
            e.0 += c;
            if mc > e.1 {
                e.1 = mc;
            }
        }
    }

    let total = (2 * b + 1).pow(4);
    println!("B = {b},  planes = {total},  proven < 3/14 by B1+ = {filtered}");
    println!("EXACT D in open gap (3/14, 1/4): {}", gap.len());
    for g in gap.iter().take(20) {
        println!("    GAP  (x,y,z,w) = {:?}", g);
    }
    let mut keys: Vec<((i64, i64), (u64, i64))> = dist.into_iter().collect();
    // sort by value descending
    keys.sort_by(|a, b| {
        let lhs = (b.0 .0 as i128) * (a.0 .1 as i128);
        let rhs = (a.0 .0 as i128) * (b.0 .1 as i128);
        lhs.cmp(&rhs)
    });
    println!("distinct exact D >= 3/14   [value : count, maxcoord]:");
    for ((n, d), (c, mc)) in keys {
        println!(
            "    {}/{} = {:.5}   count = {:<8}  maxcoord = {}",
            n,
            d,
            n as f64 / d as f64,
            c,
            mc
        );
    }
}
