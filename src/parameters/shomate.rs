use majordome_numerical::prelude::*;

pub fn cp_shomate<T: Numeric>(a: &[T], t: T) -> T {
    let tt = t / T::from_f64(1000.0);
    let poly = a[0] + tt * (a[1] + tt * (a[2] + tt * a[3]));
    poly + a[4] / (tt * tt)
}

pub fn enthalpy_shomate<T: Numeric>(a: &[T], t: T) -> T {
    let tt = t / T::from_f64(1000.0);
    let c2 = T::from_f64(1.0 / 2.0);
    let c3 = T::from_f64(1.0 / 3.0);
    let c4 = T::from_f64(1.0 / 4.0);
    let poly =
        a[5] - a[4] / tt + tt * (a[0] + tt * (c2 * a[1] + tt * (c3 * a[2] + tt * (c4 * a[3]))));
    T::from_f64(1000.0) * poly
}

pub fn entropy_shomate<T: Numeric>(a: &[T], t: T) -> T {
    let tt = t / T::from_f64(1000.0);
    let c2 = T::from_f64(1.0 / 2.0);
    let c3 = T::from_f64(1.0 / 3.0);
    a[6] + a[0] * tt.ln() - a[4] / (T::from_f64(2.0) * tt * tt)
        + tt * (a[1] + tt * (c2 * a[2] + tt * (c3 * a[3])))
}
