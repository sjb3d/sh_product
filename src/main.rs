use std::f64::consts::*;

trait Factorial {
    fn factorial(self) -> Self;
}

impl Factorial for i32 {
    fn factorial(self) -> Self {
        (1..=self).product()
    }
}

fn legendre_impl(l: i32, m: i32, cos_theta: f64, sin_theta: f64) -> f64 {
    if l == m {
        if l == 0 {
            1.0
        } else {
            ((2 * m - 1) as f64) * sin_theta * legendre_impl(l - 1, m - 1, cos_theta, sin_theta)
        }
    } else if l == m + 1 {
        ((2 * m + 1) as f64) * cos_theta * legendre_impl(l - 1, m, cos_theta, sin_theta)
    } else {
        (((2 * l - 1) as f64) * cos_theta * legendre_impl(l - 1, m, cos_theta, sin_theta)
            + ((1 - l - m) as f64) * legendre_impl(l - 2, m, cos_theta, sin_theta))
            / ((l - m) as f64)
    }
}

fn legendre(l: i32, m: i32, cos_theta: f64, sin_theta: f64) -> f64 {
    assert!(m.abs() <= l);
    legendre_impl(l, m.abs(), cos_theta, sin_theta)
}

fn lm_from_index(index: usize) -> (i32, i32)
{
    let mut l = 0;
    loop {
        let start = l*l;
        let end = start + 2*l + 1;
        if index < end {
            let l = l as i32;
            let m = ((index - start) as i32) - l;
            return (l, m);
        }
        l += 1;
    }
}

fn spherical_harmonic(l: i32, m: i32, phi: f64, cos_theta: f64, sin_theta: f64) -> f64 {
    let m_abs = m.abs();
    let const_sq = ((2 * l + 1) * (l - m_abs).factorial()) as f64
        / ((4 * (l + m_abs).factorial()) as f64 * PI);
    let legendre_part = legendre(l, m, cos_theta, sin_theta);
    let exp_part = match m {
        m if m > 0 => SQRT_2 * (m as f64 * phi).cos(),
        m if m < 0 => SQRT_2 * (-m as f64 * phi).sin(),
        _ => 1.0,
    };
    const_sq.sqrt() * legendre_part * exp_part
}

fn integrate<F>(f: F) -> f64
where
    F: Fn(f64, f64, f64) -> f64,
{
    const PHI_STEPS: usize = 256;
    const COS_THETA_STEPS: usize = PHI_STEPS / 2;
    let mut sum = 0.0;
    for cos_theta_index in 0..COS_THETA_STEPS {
        let cos_theta = 1.0 - 2.0 * ((cos_theta_index as f64) + 0.5) / (COS_THETA_STEPS as f64);
        let sin_theta = (1.0 - cos_theta * cos_theta).max(0.0).sqrt();
        for phi_index in 0..PHI_STEPS {
            let phi = 2.0 * PI * ((phi_index as f64) + 0.5) / (PHI_STEPS as f64);
            sum += f(phi, cos_theta, sin_theta);
        }
    }
    sum * (4.0 * PI) / ((PHI_STEPS * COS_THETA_STEPS) as f64)
}

const MAX_ORDER: usize = 2;
const COEFF_COUNT: usize = (MAX_ORDER + 1)*(MAX_ORDER + 1);

const TOLERANCE: f64 = 0.001;

fn main() {
    for i1 in 0..COEFF_COUNT {
        let (l1, m1) = lm_from_index(i1);
        for i2 in i1..COEFF_COUNT {
            let (l2, m2) = lm_from_index(i2);
            let integral = integrate(|phi, cos_theta, sin_theta| {
                spherical_harmonic(l1, m1, phi, cos_theta, sin_theta)
                    * spherical_harmonic(l2, m2, phi, cos_theta, sin_theta)
            });
            let expected = if i1 == i2 { 1.0 } else { 0.0 };
            if (integral - expected).abs() > TOLERANCE {
                println!(
                    "integral y_{}^{} * y_{}^{} failed: got {}, expected {}",
                    l1, m1, l2, m2, integral, expected
                );
            }
        }
    }

    for i1 in 0..COEFF_COUNT {
        let (l1, m1) = lm_from_index(i1);
        for i2 in i1..COEFF_COUNT {
            let (l2, m2) = lm_from_index(i2);
            for i3 in i2..COEFF_COUNT {
                let (l3, m3) = lm_from_index(i3);
                let integral = integrate(|phi, cos_theta, sin_theta| {
                    spherical_harmonic(l1, m1, phi, cos_theta, sin_theta)
                        * spherical_harmonic(l2, m2, phi, cos_theta, sin_theta)
                        * spherical_harmonic(l3, m3, phi, cos_theta, sin_theta)
                });
                if integral.abs() > TOLERANCE {
                    println!(
                        "integral y_{}^{} * y_{}^{} * y_{}^{} = {}",
                        l1, m1, l2, m2, l3, m3, integral
                    );
                }
            }
        }
        println!();
    }
}
