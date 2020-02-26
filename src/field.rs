//! This is documentation for the `field` module.
//!
//! The field module is meant to be used for bar.

use num_bigint::BigInt;

use crate::ops::{AddRef, DivRef, MulRef, NegRef, SubRef};

pub trait FromFactory<T: Sized>: Field {
    fn from(&self, _: T) -> <Self as Field>::Elt;
}
pub trait IntoFactory<T: Field>: Sized {
    fn lift(&self, _: T) -> <T as Field>::Elt;
}

/// Sqrt trait adds square-root calculation and quadratic-residue testing.
pub trait Sqrt {
    /// Determines whether an element is a quadratic residue.
    fn is_square(&self) -> bool;
    /// Returns one square-root if the element is a quadratic residue. Otherwise, output is arbitrary.
    /// There is no notion of principal square-root in the ouput.
    fn sqrt(&self) -> Self;
}

#[derive(Copy, Clone)]
pub enum Sgn0Endianness {
    BigEndian,
    LittleEndian,
}

pub trait Sgn0 {
    fn sgn0_be(&self) -> i32;
    fn sgn0_le(&self) -> i32;
    #[inline]
    fn sgn0(&self, s: Sgn0Endianness) -> i32 {
        match s {
            Sgn0Endianness::BigEndian => self.sgn0_be(),
            Sgn0Endianness::LittleEndian => self.sgn0_le(),
        }
    }
}

pub trait CMov: Clone {
    fn cmov(x: &Self, y: &Self, b: bool) -> Self {
        match b {
            true => y.clone(),
            false => x.clone(),
        }
    }
}

use std::ops::{Add, Div, Mul, Sub};

/// FieldElement is an element of a finite field.
pub trait FieldElement {}
impl<T> FieldElement for T where T: Add<Output = T> + Sub<Output = T> + Mul<Output = T> {}
// pub trait FieldElement: AddRef + SubRef + MulRef + DivRef + NegRef + Sqrt + CMov + Sgn0 {}

/// Field is a fabric to instante a finite field.
pub trait Field {
    /// `Elt` determines the type of field elements.
    type Elt: FieldElement;
    fn elt(&self, _: BigInt) -> Self::Elt;
    fn zero(&self) -> Self::Elt;
    fn one(&self) -> Self::Elt;
}
