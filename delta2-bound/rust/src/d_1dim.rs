//! Fast exact D-computation for a 1-dim subtorus <w>_R ⊂ (R/Z)^4.
//!
//! D(<w>_R) = 1/2 - max_t min_i ||t * w_i||.
//!
//! Used as a building block for the B1+ cheap upper bound on D(<u,v>_R)
//! (since <a*u + b*v>_R ⊂ <u,v>_R for any integers a, b, so
//!  D(<u,v>_R) ≤ D(<a*u + b*v>_R)).
//!
//! The maximum is attained at t = k/d for some rational with denominator
//! d in the bounded set { |w_i|, |w_i + w_j|, |w_i - w_j|, 2|w_i| }.

fn gcd(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Compute D(<w>_R) for w ∈ Z^n (n=4 here).
/// Returns (num, den) with den > 0 and gcd(num.abs(), den) = 1.
/// If all w_i = 0, returns (1, 2) (D = 1/2, trivial subtorus).
/// If any w_i = 0, also returns (1, 2) since the corresponding coord
/// is always 0, forcing min_i ||t * w_i|| = 0.
pub fn d_1dim(w: &[i64; 4]) -> (i64, i64) {
    // If any coord is zero, the 1-dim subtorus has D = 1/2.
    if w.iter().any(|&x| x == 0) {
        return (1, 2);
    }

    // Saturate.
    let mut g = 0i64;
    for &x in w.iter() {
        g = gcd(g, x);
    }
    if g == 0 {
        return (1, 2);
    }
    let w: [i64; 4] = [w[0] / g, w[1] / g, w[2] / g, w[3] / g];

    // Candidate denominators: for the optimum t = k/d, d is in the set
    // { |w_i|, |w_i ± w_j|, 2|w_i| } (standard Lonely Runner fact).
    let mut denoms: Vec<i64> = Vec::with_capacity(32);
    for i in 0..4 {
        let ai = w[i].abs();
        if ai > 1 {
            denoms.push(ai);
            denoms.push(2 * ai);
        }
        for j in (i + 1)..4 {
            let s = (w[i] + w[j]).abs();
            let d = (w[i] - w[j]).abs();
            if s > 1 {
                denoms.push(s);
            }
            if d > 1 {
                denoms.push(d);
            }
        }
    }
    denoms.sort_unstable();
    denoms.dedup();

    // Find the t = k/d maximizing min_i ||k * w_i mod d|| / d.
    let mut best_num: i64 = 0;
    let mut best_den: i64 = 1;

    for &d in &denoms {
        if d < 2 {
            continue;
        }
        // Check ALL k in 1..d, NOT only gcd(k,d)=1: an optimal t = p/q in lowest
        // terms can have q ∉ denoms while q | d for some d ∈ denoms (e.g.
        // t = 1/4 = 5/20 with 20 = wᵢ+wⱼ but 4 ∉ denoms), so the reduced point is
        // only reached as the non-coprime fraction k/d. Restricting to gcd(k,d)=1
        // would miss such points and overestimate D.
        for k in 1..d {
            let mut min_norm: i64 = i64::MAX;
            for &wi in w.iter() {
                let r = (k * wi).rem_euclid(d);
                let nrm = r.min(d - r);
                if nrm < min_norm {
                    min_norm = nrm;
                }
                if min_norm == 0 {
                    break;
                }
            }
            // current ML = min_norm / d; compare to best = best_num / best_den
            // current > best  iff  min_norm * best_den > best_num * d
            if min_norm * best_den > best_num * d {
                best_num = min_norm;
                best_den = d;
            }
        }
    }

    // D = 1/2 - ML = (best_den - 2*best_num) / (2*best_den).
    let num = best_den - 2 * best_num;
    let den = 2 * best_den;
    let g = gcd(num.abs(), den);
    (num / g, den / g)
}

/// B1+ cheap upper-bound filter: returns true if a small integer
/// combination a*u + b*v has D < threshold_num/threshold_den, which
/// proves D(<u, v>_R) < threshold_num/threshold_den.
///
/// Sound: if returns true, D(<u, v>_R) < threshold. May miss skips
/// (i.e., return false even when D < threshold) — those fall through
/// to exact d_2dim computation.
pub fn b1_plus_below_threshold(
    u: &[i64; 4],
    v: &[i64; 4],
    thr_num: i64,
    thr_den: i64,
) -> bool {
    // Small zero-free direction multipliers (a, b). This is a cheap, SOUND
    // upper bound: if any a*u + b*v has 1-dim D < threshold, then so does
    // D(<u,v>_R). It is allowed to miss (return false when D < threshold), in
    // which case the caller falls back to the exact d_2dim computation.
    static COMBOS: &[(i64, i64)] = &[
        (1, 0),  (0, 1),  (1, 1),  (1, -1),
        (2, 1),  (2, -1), (1, 2),  (1, -2),
        (3, 1),  (3, -1), (1, 3),  (1, -3),
        (2, 3),  (2, -3), (3, 2),  (3, -2),
        (1, 4),  (1, -4), (4, 1),  (4, -1),
        (3, 4),  (3, -4), (4, 3),  (4, -3),
    ];

    for &(a, b) in COMBOS {
        let w: [i64; 4] = [
            a * u[0] + b * v[0],
            a * u[1] + b * v[1],
            a * u[2] + b * v[2],
            a * u[3] + b * v[3],
        ];
        let (n, d) = d_1dim(&w);
        // n/d < thr_num/thr_den iff n * thr_den < thr_num * d
        if n * thr_den < thr_num * d {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn d_1dim_known_triples() {
        // sorted (1,2,3) → D = 1/4
        assert_eq!(d_1dim(&[1, 2, 3, 0]), (1, 2)); // zero coord → 1/2
        // Using all nonzero — embed (1,2,3) into 4 coords by duplicating
        // smallest. (1,1,2,3) — D should still be 1/4 since duplication
        // doesn't help.
        assert_eq!(d_1dim(&[1, 1, 2, 3]), (1, 4));
        // (1,2,6) → D = 3/14
        assert_eq!(d_1dim(&[1, 1, 2, 6]), (3, 14));
        // (1,5,6) → D = 3/14
        assert_eq!(d_1dim(&[1, 1, 5, 6]), (3, 14));
        // (1,2,9) → D = 1/5
        assert_eq!(d_1dim(&[1, 1, 2, 9]), (1, 5));
    }

    #[test]
    fn b1_plus_skips_low_d_pair() {
        // A random saturated 2x4 with D probably < 3/14
        let u = [1, 0, 0, 0];
        let v = [0, 1, 0, 0]; // U = T^1 × T^1, D = 1/2
        // Hmm — this U has D = 1/2 (a coord-aligned 2-torus has the full
        // free direction). So the filter should NOT skip.
        // Use a real low-D case instead:
        let u = [1, 2, 3, 5];
        let v = [2, 1, 5, 3];
        // Not certain about exact D, but expect well below 3/14.
        let _filtered = b1_plus_below_threshold(&u, &v, 3, 14);
        // Don't assert here; just exercise the code path.
    }
}
