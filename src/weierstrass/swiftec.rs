use num_traits::identities::Zero;

use crate::ellipticcurve::{EllipticCurve, MapToCurve};
use crate::field::{CMov, Field, Sgn0, Sqrt};
use crate::ops::FromFactory;
use crate::primefield::FpElt;
use crate::weierstrass::Curve;

#[derive(Clone)]
pub struct SwiftEC {
    e: Curve,
    az: FpElt,
    bz: FpElt,
    cz: FpElt,
    dz: FpElt,
    ez: FpElt,
    t:  FpElt,
}

impl SwiftEC {
    pub fn new(e: Curve) -> SwiftEC {
        if !SwiftEC::verify(&e) {
            panic!("wrong input parameters")
        } else {
            let f = e.get_field();
            let (mut az, mut bz, mut cz, mut dz, mut ez, mut t) = (f.one(), f.one(), f.one(), f.one(), f.one(), f.one());

            if e.a == f.from("115792089210356248762697446949407573530086143415290314195533631308867097853951") {

            // p256 case : h reducible, g irreducible. A.2
            // let s = 2 * (e.a / -3).sqrt();
            // let ub1 = (-1 / (e.b - e.a * s/3)).sqrt();
            // let ub2 = (-1 / (e.b + e.a * s/3)).sqrt();
            // let sqrtdet = (-16 * (4 * e.a * e.a + 27 * e.b * e.b)).sqrt();

            // let z = 4 * e.a * s * ub1 * ub2 * sqrtdet;
            // let az = -3 * s * (ub1 + ub2) * sqrtdet / z;
            // let bz = -2 * e.a * (ub1 - ub2) * sqrtdet / z;
            // let cz = -2 * e.a * s * (ub1 + ub2) * sqrtdet / z;
            // let dz = 16 * e.a * e.a * (ub1 - ub2) + 36 * e.b * s * (ub1 + ub2) / z;
            // let ez = 24 * e.a * e.b * (ub1 - ub2) - 8 * e.a * e.a * s * (ub1 + ub2) / z;
            az = f.from("79590650173511991483349173699886758896967575878758521100255198758162002409904");
            bz = f.from("87259034199983582448131960955643677496937372051237793188553049556356375626956");
            cz = f.from("72402878073688514558696546499041629266237135073063586190556865101410190888094");
            dz = f.from("102868106456433133884483737613882731040771060075621085360470254629499268085043");
            ez = f.from("26895936680464261060020821327532362641796128279102205464958515732833961445805");
            t = f.from("12345");

            } else if e.a == f.from("57896044618658097711785492504343953926634992332820282019728792003956564819949") {

            az = f.from("40450099932775124723841776485528502281840213886211106156430698406316951711711");
            bz = f.from("47686436168683827353029616919495895934469740993643549276545205584580349974358");
            cz = f.from("43773275001684433843520557773841752863220665495070228024025088987991857099671");
            dz = f.from("16031328859210480743063058373158317041785230125655642687073861157940246142853");
            ez = f.from("44217469855591667576766889010674012963334650726642028109726614458823388315468");
            t = f.from("12345");
            }

            SwiftEC { e, az, bz, cz, dz, ez, t }
        }
    }
    fn verify(e: &Curve) -> bool {
        let precond1 = !e.a.is_zero();              // A != 0
        let precond2 = !e.b.is_zero();              // B != 0
        precond1 && precond2
    }
}

impl MapToCurve for SwiftEC {
    type E = Curve;
    fn map(
        &self,
        u: &<<Self::E as EllipticCurve>::F as Field>::Elt,
    ) -> <Self::E as EllipticCurve>::Point {
        let cmov = FpElt::cmov;

        // finding a fixed point (X0(u), Y0(u)) on the conic
        // X0 = ( Au^2 + Bu + C ) / Z
        // Y0 = ( Du + E ) / Z
        let mut x0 = u ^ 2u32;              //  0. x0 = u^2
        x0 = &self.az * &x0;                //  1. x0 = A/Z * u^2
        let mut t1 = &self.bz * u;          //  2. t1 = B/Z * u
        x0 = &x0 + &t1 + &self.cz;          //  3. x0 = x0 + t1 + C/Z
        let mut y0 = &self.dz * u;          //  4. y0 = D/Z * u
        y0 = &y0 + &self.ez;                //  5. y0 = y0 + E/Z

        // define g(u) = u^3 + au + b & h(u) = 3u^2 + 4a (pg 12)

        let mut g = u ^ 2u32;               //  6. g = u^2
        g = g + &self.e.a;                  //  7. g = u^2 + a
        g = g * u;                          //  8. g = u^3 + au
        g = g + &self.e.b;                  //  9. g = u^3 + au + b

        let mut h = u ^ 2u32;               // 10. h = u^2
        h = &h + &h + &h;                   // 11. h = 3u^2
        t1 = &self.e.a + &self.e.a + &self.e.a + &self.e.a; // 12. t1 = 4a
        h = h + &t1;                        // 13. h = 3u^2 + 4a

        // compute new point (X(u,t), Y(u,t)) on the conic
        // X(u,t) = ( g + (h * (y0 - t * x0)^2) ) / x0 * (1 + t^2 * h)
        // Y(u,t) = y0 + t * ( X - x0 )

        let mut z1 = &self.t * &x0;         // 14. z1 = t * x0
        let mut y1 = y0 - &z1;              // 15. y1 = y0 - z1
        let mut x1 = &y1 ^ 2u32;            // 16. x1 = (y0 - t * x0)^2
        x1 = x1 * &h;                       // 17. x1 = h * (y0 - t * x0)^2
        x1 = x1 + &g;                       // 18. x1 = x1 + g
        z1 = z1 * &self.t;                  // 19. z1 = (t * x0) * t
        z1 = z1 * &h;                       // 20. z1 = z1 * h
        z1 = z1 + &x0;                      // 21. z1 = z1 + x0
        y1 = y1 * &z1;                      // 22. y1 = (y0 - t*x0) * z1
        t1 = &self.t * &x1;                 // 23. t1 = t * x1
        y1 = &y1 + &t1;                     // 24. y1 = y0 + t * ( X - x0 )

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
