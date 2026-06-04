//! Exact computation of D(U) for a 2-dim saturated subtorus U = <u, v>_R
//! of (R/Z)^4, where u, v ∈ Z^4 are integer generators (assumed
//! saturated; if not, the answer is for the R-span which agrees with
//! the closure of the integer span).
//!
//! Computes
//!   D(U) = 1/2 - max_{(α, β) ∈ [0,1)^2} min_i ||α u_i + β v_i||,
//! and the inner max is attained at vertices where two coordinate norms
//! simultaneously hit the min ("tie lines"). All candidates are
//! enumerated as intersections of distinct tie-line pairs.
//!
//! Arithmetic is i64-only in inner loops; results in i128 for safety
//! when comparing fractions across line pairs.
//!
//! Pure-rust, no extra deps.

use std::collections::HashSet;

/// Reduce (num, denom) to lowest terms with denom > 0.
fn reduce(mut num: i128, mut denom: i128) -> (i128, i128) {
    if denom < 0 {
        num = -num;
        denom = -denom;
    }
    let g = gcd_i128(num.abs(), denom);
    if g == 0 {
        return (0, 1);
    }
    (num / g, denom / g)
}

fn gcd_i128(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// A line in (R/Z)^2: N0·α + N1·β ≡ c_half/2 (mod 1), with
/// c_half ∈ {0, 1}. Normalized so (N0, N1) has its first nonzero
/// component positive, to canonicalize for dedup.
type Line = (i64, i64, u8);

fn normalize_line(mut n0: i64, mut n1: i64, c_half: u8) -> Option<Line> {
    if n0 == 0 && n1 == 0 {
        return None;
    }
    if n0 < 0 || (n0 == 0 && n1 < 0) {
        n0 = -n0;
        n1 = -n1;
    }
    Some((n0, n1, c_half))
}

/// Compute min_i ||(p · u_i + q · v_i) / d|| · d, returning an i64 in
/// [0, d/2]. Here `d > 0`, and the result is the integer ml such that
/// the actual ml-value is `ml / d`.
fn min_norm_scaled(p: i64, q: i64, d: i64, u: &[i64; 4], v: &[i64; 4]) -> i64 {
    let mut best = d / 2 + 1; // upper bound on possible result
    for i in 0..4 {
        // x = (p * u_i + q * v_i) mod d, in [0, d)
        // u_i, v_i are bounded by ~coord_max (small); p, q in [0, d).
        // Product fits in i64 for d up to ~2^31 / coord_max.
        let prod = (p as i128) * (u[i] as i128) + (q as i128) * (v[i] as i128);
        let x = prod.rem_euclid(d as i128) as i64;
        let nx = if x > d / 2 { d - x } else { x };
        if nx < best {
            best = nx;
            if best == 0 {
                return 0;
            }
        }
    }
    best
}

/// Compute D(U) for U = <u, v>_R. Returns (num, denom) in lowest terms.
pub fn d_2dim(u: &[i64; 4], v: &[i64; 4]) -> (i64, i64) {
    // 1. Collect tie lines.
    let mut lines: HashSet<Line> = HashSet::new();

    // p_i ≡ ±p_j tie lines, plus single-coord zero / half lines.
    for i in 0..4 {
        for j in (i + 1)..4 {
            // p_i ≡ p_j: (u_i - u_j) α + (v_i - v_j) β ≡ 0
            if let Some(L) = normalize_line(u[i] - u[j], v[i] - v[j], 0) {
                lines.insert(L);
            }
            // p_i ≡ -p_j: (u_i + u_j) α + (v_i + v_j) β ≡ 0
            if let Some(L) = normalize_line(u[i] + u[j], v[i] + v[j], 0) {
                lines.insert(L);
            }
        }
        // Single-coord breakpoints: p_i ≡ 0 and p_i ≡ 1/2
        if let Some(L0) = normalize_line(u[i], v[i], 0) {
            lines.insert(L0);
        }
        if let Some(L1) = normalize_line(u[i], v[i], 1) {
            lines.insert(L1);
        }
    }

    let line_list: Vec<Line> = lines.into_iter().collect();

    // 2. Track the best ml as a fraction (best_num, best_den), best_num/best_den.
    //    The "result" D = 1/2 - best_ml; we want best_ml as large as possible.
    let mut best_num: i128 = 0;
    let mut best_den: i128 = 1;

    // Helper: candidate ml = candidate_num / candidate_den; update best if larger.
    macro_rules! update_best {
        ($cnum:expr, $cden:expr) => {{
            let cn = $cnum as i128;
            let cd = $cden as i128;
            // Compare cn/cd > best_num/best_den ⟺ cn * best_den > best_num * cd
            // (with both denoms positive).
            if cn * best_den > best_num * cd {
                let (r_n, r_d) = reduce(cn, cd);
                best_num = r_n;
                best_den = r_d;
            }
        }};
    }

    // Always include (0, 0) and (1/2, 1/2) as candidates.
    update_best!(min_norm_scaled(0, 0, 2, u, v), 2);
    update_best!(min_norm_scaled(1, 1, 2, u, v), 2);

    // 3. Enumerate intersections of distinct line pairs.
    for idx_a in 0..line_list.len() {
        let (a0, a1, ca) = line_list[idx_a];
        for idx_b in (idx_a + 1)..line_list.len() {
            let (b0, b1, cb) = line_list[idx_b];
            // det = a0 * b1 - b0 * a1
            let det_i128 = (a0 as i128) * (b1 as i128) - (b0 as i128) * (a1 as i128);
            if det_i128 == 0 {
                continue;
            }
            let det = det_i128 as i64; // line entries are small; det fits in i64
            let d_abs = det.unsigned_abs() as i64;

            // Common denominator: D = 2 * |det| (because c ∈ {0, 1/2}).
            let big_d: i64 = 2 * d_abs;

            // We enumerate (k_A, k_B) ∈ [0, |det|)^2, and compute
            //   α_num = b1 * (ca + 2 k_A) - a1 * (cb + 2 k_B)
            //   β_num = a0 * (cb + 2 k_B) - b0 * (ca + 2 k_A)
            // both modulo big_d. Then the candidate (α, β) =
            //   (α_num / det, β_num / det) where α_num/β_num are taken
            //   mod (2 |det|) signed-then-reduced — equivalently we
            //   directly compute the integer-scaled min_norm with d=big_d.
            //
            // Cleaner restatement: let s_A = ca + 2 k_A and s_B = cb + 2 k_B.
            //   α = (b1·s_A - a1·s_B) / (2·det)
            //   β = (a0·s_B - b0·s_A) / (2·det)
            // So if we set
            //   p ≡ (b1·s_A - a1·s_B) / sign(det) · sign-correction ... — easier to
            //   compute α_num_modD = (b1·s_A - a1·s_B) · sign(det), reduced mod big_d,
            //   then α = α_num_modD / big_d.
            //
            // We can fold sign(det) into the modular reduction: since
            // α = num/det = (sign(det)*num)/(|det|) (numerically same as
            // num/det), and 2·det vs 2·|det| differ only by sign which gets
            // absorbed by reducing mod big_d.

            let det_sign: i64 = if det > 0 { 1 } else { -1 };

            for k_a in 0..d_abs {
                let s_a: i64 = (ca as i64) + 2 * k_a;
                for k_b in 0..d_abs {
                    let s_b: i64 = (cb as i64) + 2 * k_b;
                    // α_num = b1 * s_A - a1 * s_B  (this is over 2·det)
                    // → p = (α_num) * det_sign mod big_d, where big_d = 2·|det|
                    let alpha_num_i128 =
                        (b1 as i128) * (s_a as i128) - (a1 as i128) * (s_b as i128);
                    let beta_num_i128 =
                        (a0 as i128) * (s_b as i128) - (b0 as i128) * (s_a as i128);
                    // Multiply by det_sign to convert "over 2·det" to "over 2·|det|"
                    let p = (alpha_num_i128 * det_sign as i128).rem_euclid(big_d as i128) as i64;
                    let q = (beta_num_i128 * det_sign as i128).rem_euclid(big_d as i128) as i64;

                    let ml_int = min_norm_scaled(p, q, big_d, u, v);
                    if ml_int > 0 {
                        update_best!(ml_int, big_d);
                    }
                }
            }
        }
    }

    // 4. D = 1/2 - best_ml = (best_den - 2 * best_num) / (2 * best_den).
    let d_num = best_den - 2 * best_num;
    let d_den = 2 * best_den;
    let (n, d) = reduce(d_num, d_den);
    (n as i64, d as i64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check(u: [i64; 4], v: [i64; 4], expected_num: i64, expected_den: i64) {
        let (n, d) = d_2dim(&u, &v);
        assert_eq!(
            (n, d),
            (expected_num, expected_den),
            "u={:?}, v={:?}: got {}/{}, expected {}/{}",
            u,
            v,
            n,
            d,
            expected_num,
            expected_den
        );
    }

    #[test]
    fn u1_quarter() {
        check([0, 1, 2, 3], [1, 0, 0, 0], 1, 4);
    }

    #[test]
    fn u2_quarter() {
        check([1, 0, 1, 1], [1, 1, 0, 2], 1, 4);
    }

    #[test]
    fn delta2_fifth() {
        check([1, 1, -1, 0], [2, -2, -3, -1], 1, 5);
    }

    #[test]
    fn case_1_sixth() {
        check([1, 1, -3, -2], [1, -1, 0, -1], 1, 6);
    }

    #[test]
    fn former_5_24_sixth() {
        check([1, 1, -1, 0], [2, -2, -1, -3], 1, 6);
    }
}
