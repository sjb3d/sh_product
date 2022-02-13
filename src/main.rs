use std::{f64::consts::*, mem};

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

fn lm_from_index(index: usize) -> (i32, i32) {
    let mut l = 0;
    loop {
        let start = l * l;
        let end = start + 2 * l + 1;
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

#[derive(Debug, Clone, Copy)]
struct ConstantReference {
    index: usize,
    is_negated: bool,
}

#[derive(Debug)]
struct ConstantStore {
    unique_abs_values: Vec<f64>,
}

impl ConstantStore {
    fn new() -> Self {
        Self {
            unique_abs_values: Vec::new(),
        }
    }

    fn add(&mut self, value: f64) -> ConstantReference {
        let abs_value = value.abs();
        for (index, check_abs_value) in self.unique_abs_values.iter().enumerate() {
            if (abs_value - check_abs_value).abs() < TOLERANCE {
                return ConstantReference {
                    index,
                    is_negated: (value < 0.0),
                };
            }
        }
        let index = self.unique_abs_values.len();
        self.unique_abs_values.push(abs_value);
        ConstantReference {
            index,
            is_negated: (value < 0.0),
        }
    }
}

#[derive(Debug)]
struct Element {
    triple: (usize, usize, usize),
    cref: ConstantReference,
}

impl Element {
    fn new(mut i1: usize, mut i2: usize, mut i3: usize, cref: ConstantReference) -> Self {
        let sort_pair = |a: &mut usize, b: &mut usize| {
            if *a > *b {
                mem::swap(a, b)
            }
        };
        sort_pair(&mut i1, &mut i2);
        sort_pair(&mut i1, &mut i3);
        sort_pair(&mut i2, &mut i3);
        assert!(i1 <= i2 && i2 <= i3);
        Self {
            triple: (i1, i2, i3),
            cref,
        }
    }

    fn complete_pair(&self, i1: usize, i2: usize) -> Option<usize> {
        match self.triple {
            (p1, p2, k) if p1 == i1 && p2 == i2 => Some(k),
            (k, p1, p2) if p1 == i1 && p2 == i2 => Some(k),
            (p2, k, p1) if p1 == i1 && p2 == i2 => Some(k),
            (k, p2, p1) if p1 == i1 && p2 == i2 => Some(k),
            (p2, p1, k) if p1 == i1 && p2 == i2 => Some(k),
            (p1, k, p2) if p1 == i1 && p2 == i2 => Some(k),
            _ => None,
        }
    }
}

#[derive(Debug)]
struct PairTracker {
    counters: Vec<usize>,
}

impl PairTracker {
    fn new() -> Self {
        Self {
            counters: vec![0usize; COEFF_COUNT * (COEFF_COUNT + 1) / 2],
        }
    }

    fn index(mut i1: usize, mut i2: usize) -> usize {
        // packed index for symmetric storage
        if i2 < i1 {
            mem::swap(&mut i1, &mut i2);
        }
        i2 * (i2 + 1) / 2 + i1
    }

    fn add_pair(&mut self, i1: usize, i2: usize) {
        self.counters[Self::index(i1, i2)] += 1;
    }

    fn remove_pair(&mut self, i1: usize, i2: usize) {
        self.counters[Self::index(i1, i2)] -= 1;
    }

    fn foreach_triple<F>(triple: (usize, usize, usize), mut f: F)
    where
        F: FnMut(usize, usize),
    {
        match triple {
            (i1, i2, _) if i1 == i2 => {
                f(i1, i1);
            }
            (_, i2, i3) if i2 == i3 => {
                f(i2, i2);
            }
            (i1, _, i3) if i3 == i1 => {
                f(i3, i3);
            }
            (i1, i2, i3) => {
                f(i1, i2);
                f(i2, i3);
                f(i3, i1);
            }
        }
    }

    fn add_triple(&mut self, triple: (usize, usize, usize)) {
        Self::foreach_triple(triple, |i1, i2| self.add_pair(i1, i2));
    }

    fn remove_triple(&mut self, triple: (usize, usize, usize)) {
        Self::foreach_triple(triple, |i1, i2| self.remove_pair(i1, i2));
    }

    fn find_best(&self) -> (usize, usize) {
        let mut best_counter = 0;
        let mut best_pair = (0, 0);
        let mut index = 0;
        for i2 in 0..COEFF_COUNT {
            for i1 in 0..=i2 {
                assert_eq!(Self::index(i1, i2), index);
                let counter = self.counters[index];
                if counter > best_counter {
                    best_counter = counter;
                    best_pair = (i1, i2);
                }
                index += 1;
            }
        }
        assert!(best_counter > 0);
        best_pair
    }
}

struct Term {
    i1: usize,
    i2: usize,
    kd: Vec<(usize, ConstantReference)>,
}

impl Term {
    fn new(i1: usize, i2: usize) -> Self {
        Self {
            i1,
            i2,
            kd: Vec::new(),
        }
    }
}

trait ForEachInterspersed {
    type Item;
    fn for_each_interspersed<F, G>(self, f: F, g: G)
    where
        F: FnMut(Self::Item),
        G: FnMut();
}

impl<I> ForEachInterspersed for I
where
    I: Iterator,
{
    type Item = I::Item;
    fn for_each_interspersed<F, G>(self, mut f: F, mut g: G)
    where
        F: FnMut(Self::Item),
        G: FnMut(),
    {
        let mut is_first = true;
        for item in self {
            if is_first {
                is_first = false;
            } else {
                g();
            }
            f(item);
        }
    }
}

const MAX_ORDER: usize = 2;
const COEFF_COUNT: usize = (MAX_ORDER + 1) * (MAX_ORDER + 1);

const TOLERANCE: f64 = 0.001;

fn main() {
    // numerically integrate all unique pairs, check results
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

    // numerically integrate all unique triples, gather non-zero ones
    let mut constants = ConstantStore::new();
    let mut elements = Vec::<Element>::new();
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
                    elements.push(Element::new(i1, i2, i3, constants.add(integral)));
                }
            }
        }
    }

    // greedy pair algorithm
    // ref: "Code Generation and Factoring for Fast Evaluation of Low-order Spherical Harmonic Products and Squares" by John Snyder
    let mut pairs = PairTracker::new();
    for elem in elements.iter() {
        pairs.add_triple(elem.triple)
    }
    let mut terms = Vec::new();
    while !elements.is_empty() {
        // find the most used pair
        let (i1, i2) = pairs.find_best();
        let mut term = Term::new(i1, i2);

        // handle all elements for this pair
        let mut elem_index = 0;
        while elem_index < elements.len() {
            let elem = &elements[elem_index];
            if let Some(k) = elem.complete_pair(i1, i2) {
                term.kd.push((k, elem.cref));
                pairs.remove_triple(elem.triple);
                elements.swap_remove(elem_index);
            } else {
                elem_index += 1;
            }
        }
        terms.push(term);
    }

    // sort and emit the terms
    // ref: "Code Generation and Factoring for Fast Evaluation of Low-order Spherical Harmonic Products and Squares" by John Snyder
    terms.sort_by(|a: &Term, b: &Term| {
        if a.i2 != b.i2 {
            return a.i2.cmp(&b.i2);
        };
        a.i1.cmp(&b.i1)
    });
    for (index, value) in constants.unique_abs_values.iter().enumerate() {
        println!("const float C{index} = {value}f;");
    }
    println!();
    let print_cref = |d: &ConstantReference| {
        let index = d.index;
        if d.is_negated {
            print!("(-C{index})")
        } else {
            print!("C{index}")
        }
    };
    for term in terms.iter() {
        let Term { i1, i2, .. } = term;
        if i1 == i2 {
            let i = i1;

            if term.kd.iter().any(|(k, _)| k != i) {
                print!("ta = ");
                term.kd
                    .iter()
                    .filter(|(k, _)| k != i)
                    .for_each_interspersed(
                        |(k, d)| {
                            print_cref(d);
                            print!("*a[{k}]");
                        },
                        || print!(" + "),
                    );
                println!(";");

                print!("tb = ");
                term.kd
                    .iter()
                    .filter(|(k, _)| k != i)
                    .for_each_interspersed(
                        |(k, d)| {
                            print_cref(d);
                            print!("*b[{k}]");
                        },
                        || print!(" + "),
                    );
                println!(";");

                println!("c[{i}] += ta*b[{i}] + tb*a[{i}];");
            }
            println!("t = a[{i}]*b[{i}];");
            term.kd.iter().for_each(|(k, d)| {
                print!("c[{k}] += ");
                print_cref(d);
                println!("*t;");
            });
        } else {
            print!("ta = ");
            term.kd.iter().for_each_interspersed(
                |(k, d)| {
                    print_cref(d);
                    print!("*a[{k}]");
                },
                || print!(" + "),
            );
            println!(";");

            print!("tb = ");
            term.kd.iter().for_each_interspersed(
                |(k, d)| {
                    print_cref(d);
                    print!("*b[{k}]");
                },
                || print!(" + "),
            );
            println!(";");

            println!("c[{i1}] += ta*b[{i2}] + tb*a[{i2}];");
            println!("c[{i2}] += ta*b[{i1}] + tb*a[{i1}];");
            
            println!("t = a[{i1}]*b[{i2}] + a[{i2}]*b[{i1}];");
            term.kd.iter().for_each(|(k, d)| {
                print!("c[{k}] += ");
                print_cref(d);
                println!("*t;");
            });
        }
        println!();
    }
}
