use num_traits::identities::Zero;

use crate::ellipticcurve::{EllipticCurve, MapToCurve};
use crate::field::{CMov, Field, Sgn0, Sqrt};
use crate::ops::FromFactory;
use crate::primefield::FpElt;
use crate::weierstrass::Curve;

#[derive(Clone)]
pub struct SwiftECA0 {
    e: Curve,
    m3: FpElt,
    t: FpElt,
}

impl SwiftECA0 {
    pub fn new(e: Curve) -> SwiftECA0 {
        if !SwiftECA0::verify(&e) {
            panic!("wrong input parameters")
        } else {
            let f = e.get_field();
            let m3 = f.from("4602937940656409685400179041082242364498080236264115595900560044423621507154");
            let t = f.from("12345");
            SwiftECA0 { e, m3, t }
        }
    }
    fn verify(e: &Curve) -> bool {
        let precond1 = e.a.is_zero();               // A = 0
        let precond2 = !e.b.is_zero();              // B != 0
        precond1 && precond2
    }
}

impl MapToCurve for SwiftECA0 {
    type E = Curve;
    fn map(
        &self,
        u: &<<Self::E as EllipticCurve>::F as Field>::Elt,
    ) -> <Self::E as EllipticCurve>::Point {

        // define g(u) = u^3 + b & h(u) = 3u^2

        let mut g = u ^ 2u32;               //  1. g = u^2
        g = g * u;                          //  2. g = u^3
        g = g + &self.e.b;                  //  3. g = u^3 + b

        let mut h = u ^ 2u32;               // 4. h = u^2
        h = &h + &h + &h;                   // 5. h = 3u^2

        // compute new point (X(u,t), Y(u,t)) on the conic
        // X(u,t) = ( u^3 + b - t^2 ) / 2*t
        // Y(u,t) = ( X(u,t) + t ) / u*sqrt(-3)

        let mut x1 = u ^ 3u32;
        let mut y1 = &self.t ^ 2u32;
        x1 = x1 + &self.e.b;
        x1 = x1 - y1;
        y1 = y1 + y1;
        y1 = y1 + x1;
        let mut z1 = u * &self.m3;
        x1 = x1 * z1;
        z1 = z1 * &self.t;
        z1 = z1 * z1;

        // Compute affine point in V via ψ : S -> V
        //  x1 = x/2y - u/2, x2 = -x/2y - u/2, x3 = u + 4y^2

        // Compute projective point in surface S
        // y = (2Y1)^2, v = X1*Z1 - u*Y1*Z1, w = 2*Y1*Z1
        let mut y = &y1 + &y1;
        y = &y ^ 2u32;
        let mut v = &y1 * u;
        v = &x1 - &v;
        v = &v * &z1;
        let mut w = &y1 * &z1;
        w = &w + &w;

        // Compute affine point in V
        //  x1 = v/w, x2 = -u-v/w, x3 = u + y^2/w^2
        w  = 1u32 / &w;
        x1 = v * &w;
        let x2 = -(u + &x1);
        let mut x3 = y * w;
        x3 = x3 ^ 2u32;

        // TODO : come back to this so it's not so repetitive?
        let mut gx1 = &x1 ^ 2u32;           // gx1 = x1^2
        gx1 = gx1 + &self.e.a;              // gx1 = gx1 + a
        gx1 = gx1 * &x1;                    // gx1 = gx1 * x1
        gx1 = gx1 + &self.e.b;              // gx1 = x1^3 + ax1 + b

        let mut gx2 = &x2 ^ 2u32;           // gx2 = x2^2
        gx2 = gx2 + &self.e.a;              // gx2 = gx2 + a
        gx2 = gx2 * &x2;                    // gx2 = gx2 * x2
        gx2 = gx2 + &self.e.b;              // gx2 = x2^3 + ax2 + b

        /*
        let mut gx3 = &x3 ^ 2u32;           // gx3 = x3^2
        gx3 = gx3 + &self.e.a;              // gx3 = gx3 + a
        gx3 = gx3 * &x3;                    // gx3 = gx3 * x3
        gx3 = gx3 + &self.e.b;              // gx3 = x3^3 + ax3 + b
        */
        let e2 = gx1.is_square();
        let x = cmov(&x2, &x1, e2);
        let y2 = cmov(&gx2, &gx1, e2);
        let mut y = y2.sqrt();
        let e3 = u.sgn0() == y.sgn0();
        y = cmov(&(-&y), &y, e3);
        self.e.new_point(x, y)
    }
}
