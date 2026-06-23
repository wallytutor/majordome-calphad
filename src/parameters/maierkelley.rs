use majordome_numerical::prelude::*;

pub fn cp_maierkelley<T: Numeric>(a: &[T], t: T) -> T {
    a[0] + a[1] * t + a[2] / (t * t)
}

pub fn enthalpy_maierkelley<T: Numeric>(a: &[T], t: T, t_ref: T, h_ref: T) -> T {
    let half = T::from_f64(0.5);
    let delta_h = a[0] * (t - t_ref) + half * a[1] * (t * t - t_ref * t_ref)
        - a[2] * (T::from_f64(1.0) / t - T::from_f64(1.0) / t_ref);
    h_ref + delta_h
}

pub fn entropy_maierkelley<T: Numeric>(a: &[T], t: T, t_ref: T, s_ref: T) -> T {
    let half = T::from_f64(0.5);
    let delta_s = a[0] * (t / t_ref).ln() + a[1] * (t - t_ref)
        - half * a[2] * (T::from_f64(1.0) / (t * t) - T::from_f64(1.0) / (t_ref * t_ref));
    s_ref + delta_s
}
