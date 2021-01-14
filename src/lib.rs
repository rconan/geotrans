mod quaternion;
mod vector;

use quaternion::Quaternion;
use vector::Vector;

const BETA: f64 = 13.601685f64;
const L: f64 = 8.71;

pub enum Frame<T: Into<Quaternion>> {
    OSS(T),
    M1S1(T),
    M1S2(T),
    M1S3(T),
    M1S4(T),
    M1S5(T),
    M1S6(T),
    M1S7(T),
}
impl<T: Into<Quaternion>> Frame<T> {
    pub fn to(self, frame: Frame<&mut Vector>) {
        use Frame::*;
        fn to_local<U: Into<Quaternion>>(u: U, v: &mut Vector, o: f64) {
            let p: Quaternion = u.into();
            let (s, c) = (90. - o).to_radians().sin_cos();
            let t = Vector::from([L * c, L * s, 0.]);
            let q = Quaternion::unit(o.to_radians(), Vector::k())
                * Quaternion::unit(BETA.to_radians(), Vector::i());
            *v = Vector::from((&q * p * q.complex_conjugate() + t.into()).vector_as_slice())
        }
        {
            if let OSS(u) = self {
                match frame {
                    M1S1(v) => to_local(u, v, 0.),
                    M1S2(v) => to_local(u, v, 60.),
                    M1S3(v) => to_local(u, v, 120.),
                    M1S4(v) => to_local(u, v, 180.),
                    M1S5(v) => to_local(u, v, 240.),
                    M1S6(v) => to_local(u, v, 300.),
                    M1S7(v) | OSS(v) => {
                        *v = Vector::from(u.into().vector_as_slice());
                    }
                }
                Some(())
            } else {
                None
            }
        }
        .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use Frame::*;

    macro_rules! to_frame {
         ($($name:ident: $value:expr,)*) => {
    $(
        #[test]
        fn $name() {
            for oo in [Vector::null(), Vector::i(), Vector::j(), Vector::k()].iter() {
                println!("{:*<5}", "");
                println!("OSS  : {}", oo);
                let mut oo1 = Vector::null();
                match $value {
                    (1,) =>                 OSS(oo).to(M1S1(&mut oo1)),
                    (2,) =>                 OSS(oo).to(M1S2(&mut oo1)),
                    (3,) =>                 OSS(oo).to(M1S3(&mut oo1)),
                    (4,) =>                 OSS(oo).to(M1S4(&mut oo1)),
                    (5,) =>                 OSS(oo).to(M1S5(&mut oo1)),
                    (6,) =>                 OSS(oo).to(M1S6(&mut oo1)),
                    (7,) =>                 OSS(oo).to(M1S7(&mut oo1)),
                    (0,) =>                 OSS(oo).to(OSS(&mut oo1)),
                    _ => ()

                }
                println!("M1SX : {}", oo1);
                println!("{:*<5}", "");
            }
        }
    )*
    }
    }

    to_frame! {
        m1s1: (1,),
        m1s2: (2,),
        m1s3: (3,),
        m1s4: (4,),
        m1s5: (5,),
        m1s6: (6,),
        m1s7: (7,),
        oss: (7,),
    }
}
