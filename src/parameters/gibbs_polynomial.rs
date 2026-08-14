use majordome_numerical::prelude::*;

// Here defined as:
// G = a + bT + cT ln T + dT^2 + eT^3 + fT^4 + g/T

pub fn cp_gibbs_polynomial<T: Numeric>(a: &[T], t: T) -> T {
    // {c: -1, T**3*f: -12, T**2*e: -6, T*d: -2, g/T**2: -2}
    let c2 = T::from_f64(2.0);
    let c6 = T::from_f64(6.0);
    let c12 = T::from_f64(12.0);

    let mut val: T = T::from_f64(0.0);

    val = val - a[2]; // -c
    val = val - c2 * a[3] * t; // -2*d*T
    val = val - c6 * a[4] * (t * t); // -6*e*T^2
    val = val - c12 * a[5] * (t * t * t); // -12*f*T^3
    val = val - c2 * a[6] / (t * t); // -2*g/T^2

    val
}

pub fn enthalpy_gibbs_polynomial<T: Numeric>(a: &[T], t: T) -> T {
    // {a: 1, T*c: -1, T**2*d: -1, T**3*e: -2, T**4*f: -3, g/T: 2}
    let c2 = T::from_f64(2.0);
    let c3 = T::from_f64(3.0);

    let mut val: T = T::from_f64(0.0);

    val = val + a[0]; // a
    val = val - a[2] * t; // -c*T
    val = val - a[3] * (t * t); // -d*T^2
    val = val - c2 * a[4] * (t * t * t); // -2*e*T^3
    val = val - c3 * a[5] * (t * t * t * t); // -3*f*T^4
    val = val + c2 * a[6] / t; // 2*g/T

    val
}

pub fn entropy_gibbs_polynomial<T: Numeric>(a: &[T], t: T) -> T {
    // {b: -1, c: -1, g/T**2: 1, c*log(T): -1, T**3*f: -4, T**2*e: -3, T*d: -2}
    let c2 = T::from_f64(2.0);
    let c3 = T::from_f64(3.0);
    let c4 = T::from_f64(4.0);

    let mut val: T = T::from_f64(0.0);

    val = val - a[1]; // -b
    val = val - a[2] * (t.ln() + T::from_f64(1.0)); // -c(ln T + 1)
    val = val - c2 * a[3] * t; // -2*d*T
    val = val - c3 * a[4] * (t * t); // -3*e*T^2
    val = val - c4 * a[5] * (t * t * t); // -4*f*T^3
    val = val + a[6] / (t * t); // g/T^2

    val
}
