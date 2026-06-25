// Independent large-scale cross-check of: S1(4) cap (1/4,1/2] differs from the JK24
// progression by exactly {1/3, 2/7}.  Recomputes EXACT ML from scratch (no trust in the
// closed forms of the proof notes) and confirms:
//   FORWARD : every k = 4 (mod 8), 20 <= k <= KMAX, k != 28 is realized at D = 1/4 + 1/k
//             by the designated coprime family (U1 (1,4j) for k=4 mod16; U2 (1,7) and
//             (4m+3,8) for k=12 mod16).
//   EXCLUDE : no coprime (A,B) with A,B <= BND, via U1 {B,A,2A,3A} or U2 {A,B,A+B,A+2B},
//             realizes 1/3 (k=12) or 2/7 (k=28).
//
// Parallel (all cores), early-pruned exact i128 arithmetic, live progress + ETA.
// NOTE: total cost is ~O(KMAX^2) (each k's config grows), so big KMAX is genuinely slow;
//       KMAX=50000 (a few seconds on multicore) is already far past any plausible
//       counterexample, and the theorem is proved uniformly anyway. Bigger is just for show.
//
// Build & run:   rustc -O verify_spectrum.rs && ./verify_spectrum [KMAX] [BND]
//   e.g.  ./verify_spectrum 50000 300
// By Vico Bonfioli (Apache-2.0).

use std::env;
use std::io::Write;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

#[inline]
fn nid(num: i128, den: i128) -> i128 {
    let r = num.rem_euclid(den);
    if r > den - r { den - r } else { r }
}

// Exact ML(v) = (num, den) with ML = num/den; early-prunes candidates against running best.
fn ml(v: &[i128]) -> (i128, i128) {
    let n = v.len();
    let mut dens: Vec<i128> = Vec::with_capacity(16);
    for i in 0..n {
        let d2 = 2 * v[i].abs();
        if d2 > 0 { dens.push(d2); }
        for j in (i + 1)..n {
            let a = (v[i] - v[j]).abs();
            if a > 0 { dens.push(a); }
            let b = v[i] + v[j];
            if b > 0 { dens.push(b); }
        }
    }
    dens.sort_unstable();
    dens.dedup();
    let (mut bn, mut bd) = (0i128, 1i128);
    for &d in &dens {
        for a in 0..d {
            let mut g = d;
            let mut beat = true;
            for &vi in v {
                let r = nid(vi * a, d);
                if r < g { g = r; }
                if g * bd <= bn * d { beat = false; break; } // can't beat best; g only drops
            }
            if beat { bn = g; bd = d; }
        }
    }
    (bn, bd)
}

fn gcd(a: i128, b: i128) -> i128 { if b == 0 { a.abs() } else { gcd(b, a % b) } }

fn distinct4(v: &[i128]) -> bool {
    for i in 0..4 {
        if v[i] == 0 { return false; }
        for j in (i + 1)..4 { if v[i] == v[j] { return false; } }
    }
    true
}

fn realizes(v: &[i128], k: i128) -> bool {
    if !distinct4(v) { return false; }
    let (n, d) = ml(v);
    n * 4 * k == (k - 4) * d
}

fn config_for(k: i128) -> Vec<i128> {
    if k % 16 == 4 {
        vec![1, 2, 3, 4 * ((k - 4) / 16)] // U1 (1, 4j)
    } else if k == 44 {
        vec![1, 7, 8, 15] // U2 (A,B)=(1,7)
    } else {
        let m = (k - 60) / 16; // U2 (4m+3, 8)
        vec![4 * m + 3, 8, 4 * m + 11, 4 * m + 19]
    }
}

fn bar(label: &str, frac: f64, el: f64, eta: f64, tail: &str) {
    let w = 30usize;
    let f = (frac.clamp(0.0, 1.0) * w as f64) as usize;
    let bar: String = "#".repeat(f) + &"-".repeat(w - f);
    print!("\r{} [{}] {:5.1}%  elapsed {:6.1}s  ETA {:>6}  {}   ",
           label, bar, frac * 100.0,
           el, if eta.is_finite() { format!("{:.1}s", eta) } else { "  ...".into() }, tail);
    let _ = std::io::stdout().flush();
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let kmax: i128 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(20000);
    let bnd: i128 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(300);
    let nt = thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
    eprintln!("threads={}  KMAX={}  BND={}", nt, kmax, bnd);

    // ===== FORWARD (parallel, work-stealing; ETA weighted by k since cost ~ k) =====
    let mut targets: Vec<i128> = Vec::new();
    let mut k = 20i128;
    while k <= kmax { if k != 28 { targets.push(k); } k += 8; }
    let ntar = targets.len();
    let total_w: u64 = targets.iter().map(|&x| x as u64).sum::<u64>().max(1);
    let targets = Arc::new(targets);
    let idx = Arc::new(AtomicUsize::new(0));
    let donew = Arc::new(AtomicU64::new(0));
    let fails = Arc::new(Mutex::new(Vec::<i128>::new()));
    let t0 = Instant::now();
    let mut hs = Vec::new();
    for _ in 0..nt {
        let (tg, ix, dw, fl) = (targets.clone(), idx.clone(), donew.clone(), fails.clone());
        hs.push(thread::spawn(move || loop {
            let i = ix.fetch_add(1, Ordering::Relaxed);
            if i >= tg.len() { break; }
            let kk = tg[i];
            if !realizes(&config_for(kk), kk) { fl.lock().unwrap().push(kk); }
            dw.fetch_add(kk as u64, Ordering::Relaxed);
        }));
    }
    loop {
        let dw = donew.load(Ordering::Relaxed);
        let el = t0.elapsed().as_secs_f64();
        let frac = dw as f64 / total_w as f64;
        let eta = if frac > 0.01 { el * (1.0 - frac) / frac } else { f64::NAN };
        bar("[forward]", frac, el, eta, &format!("{} targets, k<= {}", ntar, kmax));
        if dw >= total_w { break; }
        thread::sleep(Duration::from_millis(500));
    }
    for h in hs { h.join().unwrap(); }
    let mut fwd = fails.lock().unwrap().clone();
    fwd.sort();
    println!("\r[forward] DONE  {} targets realized at D=1/4+1/k; failures: {}{}",
             ntar, if fwd.is_empty() { "NONE".to_string() } else { format!("{:?}", fwd) },
             " ".repeat(30));

    // ===== EXCLUSION (parallel over A) =====
    let aix = Arc::new(AtomicUsize::new(1));
    let adone = Arc::new(AtomicU64::new(0));
    let ebad = Arc::new(Mutex::new(Vec::<(i128, i128, i128)>::new()));
    let te = Instant::now();
    let mut eh = Vec::new();
    for _ in 0..nt {
        let (ax, ad, eb) = (aix.clone(), adone.clone(), ebad.clone());
        eh.push(thread::spawn(move || loop {
            let a = ax.fetch_add(1, Ordering::Relaxed) as i128;
            if a > bnd { break; }
            for b in 1..=bnd {
                if gcd(a, b) != 1 { continue; }
                let u1 = [b, a, 2 * a, 3 * a];
                let u2 = [a, b, a + b, a + 2 * b];
                for kk in [12i128, 28] {
                    if realizes(&u1, kk) { eb.lock().unwrap().push((a, b, -kk)); }
                    if realizes(&u2, kk) { eb.lock().unwrap().push((a, b, kk)); }
                }
            }
            ad.fetch_add(1, Ordering::Relaxed);
        }));
    }
    loop {
        let ad = adone.load(Ordering::Relaxed);
        let el = te.elapsed().as_secs_f64();
        let frac = ad as f64 / bnd as f64;
        let eta = if frac > 0.02 { el * (1.0 - frac) / frac } else { f64::NAN };
        bar("[exclude]", frac, el, eta, &format!("A,B <= {}", bnd));
        if ad >= bnd as u64 { break; }
        thread::sleep(Duration::from_millis(500));
    }
    for h in eh { h.join().unwrap(); }
    let eb = ebad.lock().unwrap().clone();
    println!("\r[exclude] DONE  A,B <= {}; realizations of 12 or 28: {}{}",
             bnd, if eb.is_empty() { "NONE".to_string() } else { format!("{:?}", eb) },
             " ".repeat(30));

    println!("\n=> exactly {{1/3, 2/7}} confirmed: {}", fwd.is_empty() && eb.is_empty());
}
