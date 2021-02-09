mod quaternion;
mod vector;

pub use quaternion::Quaternion;
pub use vector::Vector;

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
        fn operators(o: f64) -> (Vector, Quaternion) {
            let (s, c) = (90. + o).to_radians().sin_cos();
            let t = Vector::from([L * c, L * s, 0.]);
            let q = Quaternion::unit(o.to_radians(), Vector::k())
                * Quaternion::unit(BETA.to_radians(), Vector::i());
            (t, q)
        }
        fn direct<U: Into<Quaternion>>(u: U, v: &mut Vector, o: f64) {
            let (t, q) = operators(o);
            let p: Quaternion = u.into();
            *v = Vector::from((&q * p * q.complex_conjugate() + t.into()).vector_as_slice())
        }
        fn inverse<U: Into<Quaternion>>(u: U, v: &mut Vector, o: f64) {
            let (t, q) = operators(o);
            let p: Quaternion = u.into();
            *v = Vector::from((q.complex_conjugate() * (p - t.into()) * &q).vector_as_slice())
        }
        {
            match self {
                OSS(u) => match frame {
                    M1S1(v) => direct(u, v, 0.),
                    M1S2(v) => direct(u, v, -60.),
                    M1S3(v) => direct(u, v, -120.),
                    M1S4(v) => direct(u, v, -180.),
                    M1S5(v) => direct(u, v, -240.),
                    M1S6(v) => direct(u, v, -300.),
                    M1S7(v) | OSS(v) => {
                        *v = Vector::from(u.into().vector_as_slice());
                    }
                },
                M1S1(u) => match frame {
                    OSS(v) => inverse(u, v, 0.),
                    _ => unimplemented!(),
                },
                M1S2(u) => match frame {
                    OSS(v) => inverse(u, v, -60.),
                    _ => unimplemented!(),
                },
                M1S3(u) => match frame {
                    OSS(v) => inverse(u, v, -120.),
                    _ => unimplemented!(),
                },
                M1S4(u) => match frame {
                    OSS(v) => inverse(u, v, -180.),
                    _ => unimplemented!(),
                },
                M1S5(u) => match frame {
                    OSS(v) => inverse(u, v, -240.),
                    _ => unimplemented!(),
                },
                M1S6(u) => match frame {
                    OSS(v) => inverse(u, v, -300.),
                    _ => unimplemented!(),
                },
                M1S7(u) => match frame {
                    OSS(v) => {
                        *v = Vector::from(u.into().vector_as_slice());
                    }
                    _ => unimplemented!(),
                },
            }
        }
    }
}
pub fn m1_any_to_oss(sid: usize, u: Vector) -> Vector {
    use Frame::*;
    let mut v = Vector::null();
    match sid {
        1 =>                 M1S1(u).to(OSS(&mut v)),
        2 =>                 M1S2(u).to(OSS(&mut v)),
        3 =>                 M1S3(u).to(OSS(&mut v)),
        4 =>                 M1S4(u).to(OSS(&mut v)),
        5 =>                 M1S5(u).to(OSS(&mut v)),
        6 =>                 M1S6(u).to(OSS(&mut v)),
        7 =>                 M1S7(u).to(OSS(&mut v)),
        _ => ()
    }
    v
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
        oss2m1s1: (1,),
        oss2m1s2: (2,),
        oss2m1s3: (3,),
        oss2m1s4: (4,),
        oss2m1s5: (5,),
        oss2m1s6: (6,),
        oss2m1s7: (7,),
        oss2oss: (7,),
    }

    #[test]
    fn oss_m1sx_oss() {
        for oo in [Vector::null(), Vector::i(), Vector::j(), Vector::k()].iter() {
            println!("{:*<5}", "");
            println!("OSS  : {}", oo);
            let mut oo1 = Vector::null();
            OSS(oo.clone()).to(M1S4(&mut oo1));
            M1S4(oo1.clone()).to(OSS(&mut oo1));
            println!("M1SX : {}", oo1);
            println!("{:*<5}", "");
        }
    }
}
