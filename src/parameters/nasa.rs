use crate::core::R_GAS;
use majordome_numerical::prelude::*;

pub fn cp_nasa7<T: Numeric>(a: &[T], t: T) -> T {
    let poly = a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4])));
    T::from_f64(self::R_GAS) * poly
}

pub fn enthalpy_nasa7<T: Numeric>(a: &[T], t: T) -> T {
    let c2 = T::from_f64(1.0 / 2.0);
    let c3 = T::from_f64(1.0 / 3.0);
    let c4 = T::from_f64(1.0 / 4.0);
    let c5 = T::from_f64(1.0 / 5.0);
    let poly =
        a[5] / t + a[0] + t * (c2 * a[1] + t * (c3 * a[2] + t * (c4 * a[3] + t * (c5 * a[4]))));
    T::from_f64(self::R_GAS) * poly * t
}

pub fn entropy_nasa7<T: Numeric>(a: &[T], t: T) -> T {
    let c2 = T::from_f64(1.0 / 2.0);
    let c3 = T::from_f64(1.0 / 3.0);
    let c4 = T::from_f64(1.0 / 4.0);
    let poly =
        a[6] + a[0] * t.ln() + t * (a[1] + t * (c2 * a[2] + t * (c3 * a[3] + t * (c4 * a[4]))));
    T::from_f64(self::R_GAS) * poly
}

// --------------------------------------------------------------------------

pub fn cp_nasa9<T: Numeric>(a: &[T], t: T) -> T {
    let poly = a[0] / (t * t) + a[1] / t + a[2] + t * (a[3] + t * (a[4] + t * (a[5] + t * a[6])));
    T::from_f64(self::R_GAS) * poly
}

pub fn enthalpy_nasa9<T: Numeric>(a: &[T], t: T) -> T {
    let c2 = T::from_f64(1.0 / 2.0);
    let c3 = T::from_f64(1.0 / 3.0);
    let c4 = T::from_f64(1.0 / 4.0);
    let c5 = T::from_f64(1.0 / 5.0);
    let poly = -a[0] / (t * t)
        + a[1] * t.ln() / t
        + a[2]
        + a[3] * t * c2
        + a[4] * (t * t) * c3
        + a[5] * (t * t * t) * c4
        + a[6] * (t * t * t * t) * c5
        + a[7] / t;
    T::from_f64(self::R_GAS) * poly * t
}

pub fn entropy_nasa9<T: Numeric>(a: &[T], t: T) -> T {
    let c2 = T::from_f64(1.0 / 2.0);
    let c3 = T::from_f64(1.0 / 3.0);
    let c4 = T::from_f64(1.0 / 4.0);
    let poly = -a[0] / (c2 * t * t) - a[1] / t
        + a[2] * t.ln()
        + a[3] * t
        + a[4] * (t * t) * c2
        + a[5] * (t * t * t) * c3
        + a[6] * (t * t * t * t) * c4
        + a[8];
    T::from_f64(self::R_GAS) * poly
}
