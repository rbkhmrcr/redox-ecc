//! This is documentation for the `curve` module.
//!
//! The curve module is meant to be used for bar.

extern crate num_bigint;
use num_bigint::{BigInt, BigUint, ToBigInt};

use num_traits::identities::Zero;

use std::str::FromStr;

use crate::field::{FpElt, PrimeField};
use crate::montgomery::point::{Point, ProyCoordinates};
use crate::montgomery::scalar::Scalar;
use crate::montgomery::CurveID;
use crate::EllipticCurve;
use crate::Field;
use crate::{do_if_eq, FromFactory};

/// This is an elliptic curve defined in Montgomery from and defined by the equation:
/// by^2=x^3+ax^2+x.
///
#[derive(Clone, PartialEq)]
pub struct Curve {
    f: PrimeField,
    pub(super) a: FpElt,
    pub(super) b: FpElt,
    pub(super) s: FpElt,
    pub(super) r: BigUint,
    gx: FpElt,
    gy: FpElt,
    h: BigUint,
}

impl EllipticCurve for Curve {
    type F = PrimeField;
    type P = Point;
    type Coordinates = ProyCoordinates;
    type S = Scalar;
    fn new_point(&self, c: Self::Coordinates) -> Self::P {
        let e = self.clone();
        let pt = Point { e, c };
        do_if_eq!(self.is_on_curve(&pt), true, pt, ERR_ECC_NEW)
    }
    fn new_scalar(&self, k: BigInt) -> Self::S {
        Scalar::new(k, &self.r)
    }
    fn identity(&self) -> Self::P {
        let f = &self.f;
        self.new_point(ProyCoordinates {
            x: f.zero(),
            y: f.one(),
            z: f.zero(),
        })
    }
    fn is_on_curve(&self, p: &Self::P) -> bool {
        let p = &p.c;
        let l = &self.b * &(&p.y ^ 2u32) * &p.z;
        let r = &p.x * &((&p.x ^ 2u32) + &self.a * &p.x * &p.z + &(&p.z ^ 2u32));
        let e = l - r;
        e.is_zero()
    }
    fn get_order(&self) -> BigUint {
        self.r.clone()
    }
    fn get_field(&self) -> Self::F {
        self.f.clone()
    }
    fn get_cofactor(&self) -> BigInt {
        self.h.to_bigint().unwrap()
    }
    fn get_generator(&self) -> Self::P {
        self.new_point(ProyCoordinates {
            x: self.gx.clone(),
            y: self.gy.clone(),
            z: self.f.one(),
        })
    }
}

impl std::convert::From<&CurveID> for Curve {
    fn from(id: &CurveID) -> Curve {
        let params = id.0;
        let f = PrimeField::create(BigUint::from_str(params.p).unwrap());
        let a = f.from(params.a);
        let b = f.from(params.b);
        let s = f.from(params.s);
        let gx = f.from(params.gx);
        let gy = f.from(params.gy);
        let r = BigUint::from_str(params.r).unwrap();
        let h = BigUint::from_str(params.h).unwrap();
        Curve {
            f,
            a,
            b,
            s,
            r,
            h,
            gx,
            gy,
        }
    }
}

impl std::fmt::Display for Curve {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Montgomery Curve by^2=x^3+ax^2+x\na: {}\nb: {}",
            self.a, self.b,
        )
    }
}

const ERR_ECC_NEW: &str = "not valid point";
