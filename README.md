# sh_product

An implementation of [Code Generation and Factoring for Fast Evaluation of Low-order Spherical Harmonic Products and Squares](https://www.microsoft.com/en-us/research/publication/code-generation-and-factoring-for-fast-evaluation-of-low-order-spherical-harmonic-products-and-squares/) (paper by John Snyder).

The program does the following:

- Numerically integrates spherical harmonic (dual) products for all unique combinations of basis functions
  - Checks that the results are close enough to 1 or 0 as expected
  - If not you may need to increase the step count (see `PHI_STEPS` in code)
- Numerically integrates spherical harmonic triple products for all unique combinations of basis functions
  - This is parallelised using [rayon](https://github.com/rayon-rs/rayon) but may take a while for higher orders!
- Runs the algorithm from the paper to greedily pick basis function pairs that are most common
    - Emits C-like code in blocks for each picked pair, until all coefficients are done

Due to factoring out shared parts of these basis function pairs, the generated code has fewer multiplies compared with code that just sums all the non-zero parts of the triple product.

The output is order 2 (9 terms) spherical harmonics by default, but can emit higher order by changing the `MAX_ORDER` constant. (I have tried up to order 6, which generates code with 5351 multiplies in agreement with the paper, but I have not tried running the output!)

Here is the output for an order 2 multiply pasted into a C function:

```C
// compute c = a * b for functions represented as 9-term SH
void sh_product(float *c, const float *a, const float *b)
{
    float ta, tb, t;

    const float C0 = 0.2820948f;
    const float C1 = 0.12615182f;
    const float C2 = 0.21850969f;
    const float C3 = 0.25228918f;
    const float C4 = 0.18022375f;
    const float C5 = 0.15607835f;
    const float C6 = 0.09013593f;

    t = a[0]*b[0];
    c[0] = C0*t;

    ta = C0*a[0] + (-C1)*a[6] + (-C2)*a[8];
    tb = C0*b[0] + (-C1)*b[6] + (-C2)*b[8];
    c[1] = ta*b[1] + tb*a[1];
    t = a[1]*b[1];
    c[0] += C0*t;
    c[6] = (-C1)*t;
    c[8] = (-C2)*t;

    ta = C2*a[5];
    tb = C2*b[5];
    c[1] += ta*b[2] + tb*a[2];
    c[2] = ta*b[1] + tb*a[1];
    t = a[1]*b[2] + a[2]*b[1];
    c[5] = C2*t;

    ta = C0*a[0] + C3*a[6];
    tb = C0*b[0] + C3*b[6];
    c[2] += ta*b[2] + tb*a[2];
    t = a[2]*b[2];
    c[0] += C0*t;
    c[6] += C3*t;

    ta = C2*a[4];
    tb = C2*b[4];
    c[1] += ta*b[3] + tb*a[3];
    c[3] = ta*b[1] + tb*a[1];
    t = a[1]*b[3] + a[3]*b[1];
    c[4] = C2*t;

    ta = C2*a[7];
    tb = C2*b[7];
    c[2] += ta*b[3] + tb*a[3];
    c[3] += ta*b[2] + tb*a[2];
    t = a[2]*b[3] + a[3]*b[2];
    c[7] = C2*t;

    ta = C0*a[0] + (-C1)*a[6] + C2*a[8];
    tb = C0*b[0] + (-C1)*b[6] + C2*b[8];
    c[3] += ta*b[3] + tb*a[3];
    t = a[3]*b[3];
    c[0] += C0*t;
    c[6] += (-C1)*t;
    c[8] += C2*t;

    ta = (-C4)*a[6] + C0*a[0];
    tb = (-C4)*b[6] + C0*b[0];
    c[4] += ta*b[4] + tb*a[4];
    t = a[4]*b[4];
    c[6] += (-C4)*t;
    c[0] += C0*t;

    ta = C5*a[7];
    tb = C5*b[7];
    c[4] += ta*b[5] + tb*a[5];
    c[5] += ta*b[4] + tb*a[4];
    t = a[4]*b[5] + a[5]*b[4];
    c[7] += C5*t;

    ta = C0*a[0] + (-C5)*a[8] + C6*a[6];
    tb = C0*b[0] + (-C5)*b[8] + C6*b[6];
    c[5] += ta*b[5] + tb*a[5];
    t = a[5]*b[5];
    c[0] += C0*t;
    c[8] += (-C5)*t;
    c[6] += C6*t;

    ta = C0*a[0];
    tb = C0*b[0];
    c[6] += ta*b[6] + tb*a[6];
    t = a[6]*b[6];
    c[6] += C4*t;
    c[0] += C0*t;

    ta = C5*a[8] + C0*a[0] + C6*a[6];
    tb = C5*b[8] + C0*b[0] + C6*b[6];
    c[7] += ta*b[7] + tb*a[7];
    t = a[7]*b[7];
    c[8] += C5*t;
    c[0] += C0*t;
    c[6] += C6*t;

    ta = C0*a[0] + (-C4)*a[6];
    tb = C0*b[0] + (-C4)*b[6];
    c[8] += ta*b[8] + tb*a[8];
    t = a[8]*b[8];
    c[0] += C0*t;
    c[6] += (-C4)*t;

    // multiply count = 120
    // addition count = 74
}
```
