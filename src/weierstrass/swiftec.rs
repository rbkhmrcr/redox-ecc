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
}

impl SwiftEC {
    pub fn new(e: Curve) -> SwiftEC {
        if !SwiftEC::verify(&e) {
            panic!("wrong input parameters")
        } else {
            // p256 case : h reducible, g irreducible. A.2
            // let s = 2 * (e.a / -3).sqrt();
            // let ub1 = (-1 / (e.b - e.a * s/3)).sqrt();
		        // let ub2 = (-1 / (e.b + e.a * s/3)).sqrt();
            // let sqrtdet = (-16 * (4 * e.a * e.a + 27 * e.b * e.b)).sqrt();

            // let z = 4 * e.a * s * ub1 * ub2 * sqrtdet;
            // let az = -3 * s * (ub1 + ub2) * sqrtdet / z;
            let f = e.get_field();
            let az = f.from("79590650173511991483349173699886758896967575878758521100255198758162002409904");
            // let bz = -2 * e.a * (ub1 - ub2) * sqrtdet / z;
            let bz = f.from("87259034199983582448131960955643677496937372051237793188553049556356375626956");
            // let cz = -2 * e.a * s * (ub1 + ub2) * sqrtdet / z;
            let cz = f.from("72402878073688514558696546499041629266237135073063586190556865101410190888094");
            // let dz = 16 * e.a * e.a * (ub1 - ub2) + 36 * e.b * s * (ub1 + ub2) / z;
            let dz = f.from("102868106456433133884483737613882731040771060075621085360470254629499268085043");
            // let ez = 24 * e.a * e.b * (ub1 - ub2) - 8 * e.a * e.a * s * (ub1 + ub2) / z;
            let ez = f.from("26895936680464261060020821327532362641796128279102205464958515732833961445805");
            SwiftEC { e, az, bz, cz, dz, ez }
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
        let f = self.e.get_field();
        let cmov = FpElt::cmov;

        let mut t1 = u ^ 2u32; //         0.   t1 = u^2
        t1 = &self.z * &t1; //            1.   t1 = Z * u^2
        let mut t2 = &t1 ^ 2u32; //       2.   t2 = t1^2
        let mut x1 = &t1 + &t2; //        3.   x1 = t1 + t2
        x1 = 1u32 / &x1; //               4.   x1 = inv0(x1)
        let e1 = x1.is_zero(); //         5.   e1 = x1 == 0
        x1 = x1 + f.one(); //             6.   x1 = x1 + 1
        x1 = cmov(&x1, &self.c2, e1); //  7.   x1 = CMOV(x1, c2, e1)
        x1 = x1 * &self.c1; //            8.   x1 = x1 * c1
        let mut gx1 = &x1 ^ 2u32; //      9.  gx1 = x1^2
        gx1 = gx1 + &self.e.a; //         10. gx1 = gx1 + A
        gx1 = gx1 * &x1; //               11. gx1 = gx1 * x1
        gx1 = gx1 + &self.e.b; //         12. gx1 = gx1 + B
        let x2 = &t1 * &x1; //            13.  x2 = t1 * x1
        t2 = t1 * t2; //                  14.  t2 = t1 * t2
        let gx2 = &gx1 * &t2; //          15. gx2 = gx1 * t2
        let e2 = gx1.is_square(); //      16.  e2 = is_square(gx1)
        let x = cmov(&x2, &x1, e2); //    17.   x = CMOV(x2, x1, e2)
        let y2 = cmov(&gx2, &gx1, e2); // 18.  y2 = CMOV(gx2, gx1, e2)
        let mut y = y2.sqrt(); //         19.   y = sqrt(y2)
        let e3 = u.sgn0() == y.sgn0(); // 20.  e3 = sgn0(u) == sgn0(y)
        y = cmov(&(-&y), &y, e3); //      21.   y = CMOV(-y, y, e3)
                                  //
//        self.e.new_point(x, y)
//       #Evalaute initial point in conic X(u),Y(u)
//    X0 = fp_mul(X[2], u)
//    X0 = fp_add(X0, X[1])
//    X0 = fp_mul(X0, u)
//    X0 = fp_add(X0, X[0])
//    Y0 = fp_mul(Y[1], u)
//    Y0 = fp_add(Y0, Y[0])
//
//    #print("Test X0,Y0,Z0",X0**2+(3*u**2+4*a)*Y0**2+(u**3+a*u+b)==0)
//
//    # Evaluate f(u)=3u^2+4a and g(u)=u^3+au+b
//    f = fp_sqr(u)
//    g = fp_add(f, a)
//    g = fp_mul(g, u)
//    g = fp_add(g, b)
//    Z1 = fp_add(f, f)
//    f = fp_add(Z1, f)
//    f = fp_add(f, ax4)
//
//    #print("Test f,g",X0**2+f*Y0**2+g==0)
//
//    # Compute new point in conic
//    #   X1 = f*(Y0-t*X0)^2 + g
//    #   Z1 = X0(1 + f*t^2)
//    #   Y1 = Z1*Y0 + t*(X - Z*X0)
//    Z1 = fp_mul(t, X0)
//    Y1 = fp_sub(Y0, Z1)
//    X1 = fp_sqr(Y1)
//    X1 = fp_mul(X1, f)
//    X1 = fp_add(X1, g)
//    Z1 = fp_mul(Z1, t)
//    Z1 = fp_mul(Z1, f)
//    Z1 = fp_add(Z1, X0)
//    Y1 = fp_mul(Y1, Z1)
//    tX = fp_mul(t, X1)
//    Y1 = fp_add(Y1, tX)
//    #assert(not Z1 == 0)
//
//    #print("Test X1,Y1,Z1",X1**2+(3*u**2+4*a)*Y1**2+(u**3+a*u+b)*Z1**2==0)
//
//    # Compute projective point in surface S
//    #   y = (2Y1)^2
//    #   v = X1*Z1 - u*Y1*Z1
//    #   w = 2*Y1*Z1
//    y = fp_add(Y1, Y1)
//    y = fp_sqr(y)
//    v = fp_mul(Y1, u)
//    v = fp_sub(X1, v)
//    v = fp_mul(v, Z1)
//    w = fp_mul(Y1, Z1)
//    w = fp_add(w, w)
//
//    #print("Tests v,y,w", y**2*(w**2*u**2 + w*v*u + v**2 + w**2*a) == -w**4*(u**3 + a*u + b))
//
//
//
//    # Compute affine point in V
//    #   x1 = v/w
//    #   x2 = -u-v/w
//    #   x3 = u + y^2/w^2
//    try:
//        w = fp_inv(w)
//    except ZeroDivisionError:
//        raise PointAtInfinity
//    x1 = fp_mul(v, w)
//    x2 = fp_add(u, x1)
//    x2 = fp_neg(x2)
//    x3 = fp_mul(y, w)
//    x3 = fp_sqr(x3)
//    x3 = fp_add(u, x3)
//
//    # Compute g(x_i)
//    y21 = fp_sqr(x1)
//    y21 = fp_add(y21, a)
//    y21 = fp_mul(y21, x1)
//    y21 = fp_add(y21, b)
//
//    y22 = fp_sqr(x2)
//    y22 = fp_add(y22, a)
//    y22 = fp_mul(y22, x2)
//    y22 = fp_add(y22, b)
//
//    y23 = fp_sqr(x3)
//    y23 = fp_add(y23, a)
//    y23 = fp_mul(y23, x3)
//    y23 = fp_add(y23, b)
//
//    # Find the square
//    c2 = fp_jacobi(y22)
//    c3 = fp_jacobi(y23)
//    x1, x2 = fp_cswap(c2, x1, x2)
//    y21, y22 = fp_cswap(c2, y21, y22)
//    x1, x3 = fp_cswap(c3, x1, x3)
//    y21, y23 = fp_cswap(c3, y21, y23)
//
//    # Find the square-root and choose sign
//    y21 = fp_sqrt(y21)
//    y22 = fp_neg(y21)
//    c1 = ((int(y21) % 2) ^ (int(s) % 2))
//    y21, y22 = fp_cswap(c1, y21, y22)
//    return x1, y21

    }
}
