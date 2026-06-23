use majordome_numerical::prelude::*;

pub fn cp_gibbs_polynomial<T: Numeric>(a: &[T], t: T) -> T {
    // G = a + bT + cT ln T + dT^2 + eT^3 + fT^4 + g/T
    // S = -dG/dT = -b - c(ln T + 1) - 2dT - 3eT^2 - 4fT^3 + g/T^2
    // H = G + TS = a - cT + dT^2 + 2eT^3 + 3fT^4 + 2g/T
    // Cp = dH/dT = -c + 2dT + 6eT^2 + 12fT^3 - 2g/T^2
    let c2 = T::from_f64(2.0);
    let c6 = T::from_f64(6.0);
    let c12 = T::from_f64(12.0);
    (-a[2]) + c2 * a[3] * t + c6 * a[4] * (t * t) + c12 * a[5] * (t * t * t) - c2 * a[6] / (t * t)
}

pub fn enthalpy_gibbs_polynomial<T: Numeric>(a: &[T], t: T) -> T {
    // H = a - cT + dT^2 + 2eT^3 + 3fT^4 + 2g/T
    let c2 = T::from_f64(2.0);
    let c3 = T::from_f64(3.0);
    a[0] - a[2] * t
        + a[3] * (t * t)
        + c2 * a[4] * (t * t * t)
        + c3 * a[5] * (t * t * t * t)
        + c2 * a[6] / t
}

pub fn entropy_gibbs_polynomial<T: Numeric>(a: &[T], t: T) -> T {
    // S = -b - c(ln T + 1) - 2dT - 3eT^2 - 4fT^3 + g/T^2
    let c2 = T::from_f64(2.0);
    let c3 = T::from_f64(3.0);
    let c4 = T::from_f64(4.0);
    (-a[1])
        - a[2] * (t.ln() + T::from_f64(1.0))
        - c2 * a[3] * t
        - c3 * a[4] * (t * t)
        - c4 * a[5] * (t * t * t)
        + a[6] / (t * t)
}
