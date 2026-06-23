use super::constant::T_REF;
use super::substance::Substance;
use crate::parameters::*;
use majordome_numerical::prelude::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AggregationType {
    Solid,
    Liquid,
    Gas,
}

#[derive(Debug, Clone)]
pub struct TemperatureRange {
    pub t_min: f64,
    pub t_max: f64,
    pub model: Parameterization,
}

#[derive(Debug, Clone)]
pub enum Parameterization {
    MaierKelley {
        a: f64,
        b: f64,
        c: f64,
        h_ref: f64,
        s_ref: f64,
    },

    NASA7 {
        a1: f64,
        a2: f64,
        a3: f64,
        a4: f64,
        a5: f64,
        a6: f64,
        a7: f64,
    },

    NASA9 {
        a1: f64,
        a2: f64,
        a3: f64,
        a4: f64,
        a5: f64,
        a6: f64,
        a7: f64,
        a8: f64,
        a9: f64,
    },

    Shomate {
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        e: f64,
        f: f64,
        g: f64,
        h: f64,
    },

    GibbsPolynomial {
        a: f64,
        b: f64,
        c: f64,
        d: f64,
        e: f64,
        f: f64,
        g: f64,
    },

    Compound {
        components: Vec<(Box<Substance>, f64)>,
        deviation: Box<Parameterization>,
    },
}

impl Parameterization {
    pub fn cp<T: Numeric>(&self, t: T) -> T {
        match self {
            Parameterization::MaierKelley { a, b, c, .. } => {
                let coeffs = [T::from_f64(*a), T::from_f64(*b), T::from_f64(*c)];
                cp_maierkelley(&coeffs, t)
            }
            Parameterization::NASA7 {
                a1,
                a2,
                a3,
                a4,
                a5,
                a6,
                a7,
            } => {
                let coeffs = [
                    T::from_f64(*a1),
                    T::from_f64(*a2),
                    T::from_f64(*a3),
                    T::from_f64(*a4),
                    T::from_f64(*a5),
                    T::from_f64(*a6),
                    T::from_f64(*a7),
                ];
                cp_nasa7(&coeffs, t)
            }
            Parameterization::NASA9 {
                a1,
                a2,
                a3,
                a4,
                a5,
                a6,
                a7,
                a8,
                a9,
            } => {
                let coeffs = [
                    T::from_f64(*a1),
                    T::from_f64(*a2),
                    T::from_f64(*a3),
                    T::from_f64(*a4),
                    T::from_f64(*a5),
                    T::from_f64(*a6),
                    T::from_f64(*a7),
                    T::from_f64(*a8),
                    T::from_f64(*a9),
                ];
                cp_nasa9(&coeffs, t)
            }
            Parameterization::Shomate {
                a,
                b,
                c,
                d,
                e,
                f,
                g,
                h,
            } => {
                let coeffs = [
                    T::from_f64(*a),
                    T::from_f64(*b),
                    T::from_f64(*c),
                    T::from_f64(*d),
                    T::from_f64(*e),
                    T::from_f64(*f),
                    T::from_f64(*g),
                    T::from_f64(*h),
                ];
                cp_shomate(&coeffs, t)
            }
            Parameterization::GibbsPolynomial {
                a,
                b,
                c,
                d,
                e,
                f,
                g,
            } => {
                let coeffs = [
                    T::from_f64(*a),
                    T::from_f64(*b),
                    T::from_f64(*c),
                    T::from_f64(*d),
                    T::from_f64(*e),
                    T::from_f64(*f),
                    T::from_f64(*g),
                ];
                cp_gibbs_polynomial(&coeffs, t)
            }
            Parameterization::Compound {
                components,
                deviation,
            } => {
                let mut cp = deviation.cp(t);
                for (substance, coeff) in components {
                    cp = cp + T::from_f64(*coeff) * substance.cp(t);
                }
                cp
            }
        }
    }

    pub fn enthalpy<T: Numeric>(&self, t: T) -> T {
        match self {
            Parameterization::MaierKelley { a, b, c, h_ref, .. } => {
                let coeffs = [T::from_f64(*a), T::from_f64(*b), T::from_f64(*c)];
                enthalpy_maierkelley(&coeffs, t, T::from_f64(T_REF), T::from_f64(*h_ref))
            }
            Parameterization::NASA7 {
                a1,
                a2,
                a3,
                a4,
                a5,
                a6,
                a7,
            } => {
                let coeffs = [
                    T::from_f64(*a1),
                    T::from_f64(*a2),
                    T::from_f64(*a3),
                    T::from_f64(*a4),
                    T::from_f64(*a5),
                    T::from_f64(*a6),
                    T::from_f64(*a7),
                ];
                enthalpy_nasa7(&coeffs, t)
            }
            Parameterization::NASA9 {
                a1,
                a2,
                a3,
                a4,
                a5,
                a6,
                a7,
                a8,
                a9,
            } => {
                let coeffs = [
                    T::from_f64(*a1),
                    T::from_f64(*a2),
                    T::from_f64(*a3),
                    T::from_f64(*a4),
                    T::from_f64(*a5),
                    T::from_f64(*a6),
                    T::from_f64(*a7),
                    T::from_f64(*a8),
                    T::from_f64(*a9),
                ];
                enthalpy_nasa9(&coeffs, t)
            }
            Parameterization::Shomate {
                a,
                b,
                c,
                d,
                e,
                f,
                g,
                h,
            } => {
                let coeffs = [
                    T::from_f64(*a),
                    T::from_f64(*b),
                    T::from_f64(*c),
                    T::from_f64(*d),
                    T::from_f64(*e),
                    T::from_f64(*f),
                    T::from_f64(*g),
                    T::from_f64(*h),
                ];
                enthalpy_shomate(&coeffs, t)
            }
            Parameterization::GibbsPolynomial {
                a,
                b,
                c,
                d,
                e,
                f,
                g,
            } => {
                let coeffs = [
                    T::from_f64(*a),
                    T::from_f64(*b),
                    T::from_f64(*c),
                    T::from_f64(*d),
                    T::from_f64(*e),
                    T::from_f64(*f),
                    T::from_f64(*g),
                ];
                enthalpy_gibbs_polynomial(&coeffs, t)
            }
            Parameterization::Compound {
                components,
                deviation,
            } => {
                let mut h = deviation.enthalpy(t);
                for (substance, coeff) in components {
                    h = h + T::from_f64(*coeff) * substance.enthalpy(t);
                }
                h
            }
        }
    }

    pub fn entropy<T: Numeric>(&self, t: T) -> T {
        match self {
            Parameterization::MaierKelley { a, b, c, s_ref, .. } => {
                let coeffs = [T::from_f64(*a), T::from_f64(*b), T::from_f64(*c)];
                entropy_maierkelley(&coeffs, t, T::from_f64(T_REF), T::from_f64(*s_ref))
            }
            Parameterization::NASA7 {
                a1,
                a2,
                a3,
                a4,
                a5,
                a6,
                a7,
            } => {
                let coeffs = [
                    T::from_f64(*a1),
                    T::from_f64(*a2),
                    T::from_f64(*a3),
                    T::from_f64(*a4),
                    T::from_f64(*a5),
                    T::from_f64(*a6),
                    T::from_f64(*a7),
                ];
                entropy_nasa7(&coeffs, t)
            }
            Parameterization::NASA9 {
                a1,
                a2,
                a3,
                a4,
                a5,
                a6,
                a7,
                a8,
                a9,
            } => {
                let coeffs = [
                    T::from_f64(*a1),
                    T::from_f64(*a2),
                    T::from_f64(*a3),
                    T::from_f64(*a4),
                    T::from_f64(*a5),
                    T::from_f64(*a6),
                    T::from_f64(*a7),
                    T::from_f64(*a8),
                    T::from_f64(*a9),
                ];
                entropy_nasa9(&coeffs, t)
            }
            Parameterization::Shomate {
                a,
                b,
                c,
                d,
                e,
                f,
                g,
                h,
            } => {
                let coeffs = [
                    T::from_f64(*a),
                    T::from_f64(*b),
                    T::from_f64(*c),
                    T::from_f64(*d),
                    T::from_f64(*e),
                    T::from_f64(*f),
                    T::from_f64(*g),
                    T::from_f64(*h),
                ];
                entropy_shomate(&coeffs, t)
            }
            Parameterization::GibbsPolynomial {
                a,
                b,
                c,
                d,
                e,
                f,
                g,
            } => {
                let coeffs = [
                    T::from_f64(*a),
                    T::from_f64(*b),
                    T::from_f64(*c),
                    T::from_f64(*d),
                    T::from_f64(*e),
                    T::from_f64(*f),
                    T::from_f64(*g),
                ];
                entropy_gibbs_polynomial(&coeffs, t)
            }
            Parameterization::Compound {
                components,
                deviation,
            } => {
                let mut s = deviation.entropy(t);
                for (substance, coeff) in components {
                    s = s + T::from_f64(*coeff) * substance.entropy(t);
                }
                s
            }
        }
    }
}
