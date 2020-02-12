//! This is documentation for the `field` module.
//!
//! The field module is meant to be used for bar.

extern crate num_bigint;
extern crate num_integer;

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_bigint::ToBigInt;
use num_integer::Integer;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, std::cmp::PartialEq)]
pub struct PrimeField {
    modulus: BigInt,
}

impl std::fmt::Display for PrimeField {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "GF({})", self.modulus)
    }
}

impl PrimeField {
    pub fn new(modulus: BigUint) -> PrimeField {
        let modulus = modulus.to_bigint().unwrap();
        PrimeField { modulus }
    }
    pub fn elt(&self, n: i64) -> PrimeFieldElement {
        self.from(BigInt::from(n))
    }
    #[inline]
    pub fn from_str(&self, s: &str) -> PrimeFieldElement {
        use std::str::FromStr;
        PrimeFieldElement {
            f: self.clone(),
            n: BigInt::from_str(s).unwrap(),
        }
    }
    pub fn zero(&self) -> PrimeFieldElement {
        self.elt(0)
    }
    pub fn one(&self) -> PrimeFieldElement {
        self.elt(1)
    }
    #[inline]
    fn neg_mod(&self, x: &PrimeFieldElement) -> PrimeFieldElement {
        self.from(-&x.n)
    }
    #[inline]
    fn add_mod(&self, x: &PrimeFieldElement, y: &PrimeFieldElement) -> PrimeFieldElement {
        self.from(&x.n + &y.n)
    }
    #[inline]
    fn sub_mod(&self, x: &PrimeFieldElement, y: &PrimeFieldElement) -> PrimeFieldElement {
        self.from(&x.n - &y.n)
    }
    #[inline]
    fn mul_mod(&self, x: &PrimeFieldElement, y: &PrimeFieldElement) -> PrimeFieldElement {
        self.from(&x.n * &y.n)
    }
    #[inline]
    fn from(&self, n: BigInt) -> PrimeFieldElement {
        let n = n.mod_floor(&self.modulus);
        PrimeFieldElement { f: self.clone(), n }
    }
}

#[derive(Clone, std::cmp::PartialEq)]
pub struct PrimeFieldElement {
    f: PrimeField,
    n: BigInt,
}

impl std::fmt::Display for PrimeFieldElement {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.n)
    }
}

const ERR_BIN_OP: &'static str = "elements of different fields";
const ERR_INV_OP: &'static str = "numerator must be 1u32";

macro_rules! impl_bin_op {
    ($trait:ident, $name:ident, $method:ident) => {
        impl<'a, 'b, 'c> $trait<&'b PrimeFieldElement> for &'a PrimeFieldElement {
            type Output = PrimeFieldElement;
            #[inline]
            fn $name(self, other: &PrimeFieldElement) -> Self::Output {
                do_if_eq!(self.f, other.f, self.f.$method(&self, &other), ERR_BIN_OP)
            }
        }
        impl<'a, 'c> $trait<&'a PrimeFieldElement> for PrimeFieldElement {
            type Output = PrimeFieldElement;
            #[inline]
            fn $name(self, other: &Self) -> Self::Output {
                do_if_eq!(self.f, other.f, self.f.$method(&self, &other), ERR_BIN_OP)
            }
        }
        impl<'c> $trait for PrimeFieldElement {
            type Output = PrimeFieldElement;
            #[inline]
            fn $name(self, other: Self) -> Self::Output {
                do_if_eq!(self.f, other.f, self.f.$method(&self, &other), ERR_BIN_OP)
            }
        }
    };
}

impl_bin_op!(Add, add, add_mod);
impl_bin_op!(Sub, sub, sub_mod);
impl_bin_op!(Mul, mul, mul_mod);

impl<'a> Neg for &'a PrimeFieldElement {
    type Output = PrimeFieldElement;
    fn neg(self) -> Self::Output {
        self.f.neg_mod(&self)
    }
}

impl<'a> Div<&'a PrimeFieldElement> for u32 {
    type Output = PrimeFieldElement;
    #[inline]
    fn div(self, other: &PrimeFieldElement) -> Self::Output {
        do_if_eq!(
            self,
            1u32,
            {
                let p = &other.f.modulus;
                let p_minus_2 = p.sub(2i32);
                PrimeFieldElement {
                    f: other.f.clone(),
                    n: other.n.modpow(&p_minus_2, &p),
                }
            },
            ERR_INV_OP
        )
    }
}

impl num_traits::identities::Zero for PrimeFieldElement {
    fn zero() -> Self {
        unimplemented!()
    }
    fn is_zero(&self) -> bool {
        self.n.is_zero()
    }
    fn set_zero(&mut self) {
        self.n.set_zero();
    }
}

impl num_traits::identities::One for PrimeFieldElement {
    fn one() -> Self {
        unimplemented!()
    }
    fn is_one(&self) -> bool {
        self.n.is_one()
    }
    fn set_one(&mut self) {
        self.n.set_one();
    }
}
