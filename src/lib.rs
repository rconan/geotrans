//! Geometric transformation for the GMT segmented mirrors

use std::marker::PhantomData;

mod quaternion;
mod vector;

pub use quaternion::Quaternion;
pub use vector::Vector;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("The segment id {0} is not in the range [1,7]")]
    SegmentId(i32),
}

const BETA: f64 = 13.601685f64;
const L: f64 = 8.71;

/// Conic surface
#[derive(Debug, Clone)]
pub struct Conic {
    /// Radius of curvature
    radius: f64,
    /// Conic constant
    constant: f64,
}
impl Conic {
    /// GMT M1 conic surface
    pub fn m1() -> Self {
        Self {
            radius: 36f64,
            constant: -0.9982857,
        }
    }
    /// GMT M2 conic surface
    pub fn m2() -> Self {
        Self {
            radius: -4.1639009,
            constant: -0.71692784,
        }
    }
    /// Conic surface height
    pub fn height(&self, rho: f64) -> f64 {
        let c = self.radius.abs();
        let rho2 = rho * rho;
        self.radius.signum() * rho2 / (c + (c * c - (self.constant + 1f64) * rho2).sqrt())
    }
}

/// Type representing the GMT primary mirror
#[derive(Default, Debug, Clone)]
pub struct M1;
/// Type representing the GMT secondary mirror
#[derive(Default, Debug, Clone)]
pub struct M2;

/// GMT segmented mirror
#[derive(Debug, Clone)]
pub struct Segment<M> {
    /// Segment # id
    id: i32,
    /// Mirror height
    height: f64,
    /// Radial inclination angle
    beta: Option<f64>,
    /// Distance to the origin
    distance: Option<f64>,
    /// Segment clocking angle \[degree\]
    cloking: Option<i32>,
    /// Conic surface
    conic: Conic,
    mirror: PhantomData<M>,
}

/// Segment specialization traits
pub trait SegmentTrait {
    fn new(o: i32) -> Result<Self, Error>
    where
        Self: Sized;
    fn rotation(&self) -> Option<Quaternion>;
}
impl<M> Segment<M> {
    /// Returns a [`Vector`] with the segment origin coordinates in the OSS
    pub fn translation(&self) -> Vector {
        if self.id < 7 {
            let o = self.cloking.unwrap() as f64;
            let d = self.distance.unwrap();
            let z = self.conic.height(d);
            let (s, c) = (90. + o).to_radians().sin_cos();
            Vector::from([d * c, d * s, self.height + z])
        } else {
            Vector::from([0., 0., self.height])
        }
    }
}
impl SegmentTrait for Segment<M1> {
    /// Returns [`M1`] [`Segment`] `id`
    fn new(id: i32) -> Result<Self, Error> {
        match id {
            7 => Ok(Self {
                id: 7,
                height: 3.9,
                beta: None,
                distance: None,
                cloking: None,
                conic: Conic::m1(),
                mirror: PhantomData,
            }),
            id @ 1..=6 => Ok(Self {
                id,
                height: 3.9,
                beta: Some(13.601685f64),
                distance: Some(8.71),
                cloking: Some(-60i32 * (id - 1)),
                conic: Conic::m1(),
                mirror: PhantomData,
            }),
            _ => Err(Error::SegmentId(id)),
        }
    }
    /// Returns a [`Quaternion`] representing the 3D rotation of a [`M1`] [`Segment`] frame in the OSS
    fn rotation(&self) -> Option<Quaternion> {
        if self.id < 7 {
            let o = self.cloking.unwrap() as f64;
            Some(
                Quaternion::unit(o.to_radians(), Vector::k())
                    * Quaternion::unit(self.beta.unwrap().to_radians(), Vector::i()),
            )
        } else {
            None
        }
    }
}
impl SegmentTrait for Segment<M2> {
    /// Returns [`M2`] [`Segment`] `id`
    fn new(id: i32) -> Result<Self, Error> {
        match id {
            7 => Ok(Self {
                id: 7,
                height: 3.9 + 20.26247614,
                beta: None,
                distance: None,
                cloking: None,
                conic: Conic::m2(),
                mirror: PhantomData,
            }),
            id @ 1..=6 => Ok(Self {
                id,
                height: 3.9 + 20.26247614,
                beta: Some(14.777498),
                distance: Some(1.08774),
                cloking: Some(180i32 - 60i32 * (id - 1)),
                conic: Conic::m2(),
                mirror: PhantomData,
            }),
            _ => Err(Error::SegmentId(id)),
        }
    }
    /// Returns a [`Quaternion`] representing the 3D rotation of a [`M2`] [`Segment`] frame in the OSS
    fn rotation(&self) -> Option<Quaternion> {
        if self.id < 7 {
            let o = self.cloking.unwrap() as f64;
            Some(
                Quaternion::unit(o.to_radians(), Vector::k())
                    * Quaternion::unit(std::f64::consts::PI, Vector::j())
                    * Quaternion::unit(self.beta.unwrap().to_radians(), Vector::i()),
            )
        } else {
            None
        }
    }
}

/// Geometric transformation with respect to the OSS coordinate system
pub trait Transform {
    /// Transforms the coordinates given in the OSS into a segment
    fn fro<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<Self, Error>
    where
        Self: Into<Vector>,
        Vector: Into<Self>,
        Segment<M>: SegmentTrait,
    {
        maybe_segment.map(|segment| {
            let u: Vector = self.into();
            let t = segment.translation();
            if let Some(q) = segment.rotation() {
                let p: Quaternion = From::<Vector>::from(u);
                let v = Vector::from(
                    (q.complex_conjugate() * (p - From::<Vector>::from(t)) * &q).vector_as_slice(),
                );
                v.into()
            } else {
                (u - t).into()
            }
        })
    }
    /// Transforms a the vector given in the OSS into a segment
    fn vfrov<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<Self, Error>
    where
        Self: Into<Vector>,
        Vector: Into<Self>,
        Segment<M>: SegmentTrait,
    {
        maybe_segment.map(|segment| {
            if let Some(q) = segment.rotation() {
                let u: Vector = self.into();
                let p: Quaternion = From::<Vector>::from(u);
                let v = Vector::from((q.complex_conjugate() * p * &q).vector_as_slice());
                v.into()
            } else {
                self
            }
        })
    }
    /// Transforms the coordinates of a segment into the OSS
    fn to<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<Self, Error>
    where
        Self: Into<Vector>,
        Vector: Into<Self>,
        Segment<M>: SegmentTrait,
    {
        maybe_segment.map(|segment| {
            let u: Vector = self.into();
            let t = segment.translation();
            if let Some(q) = segment.rotation() {
                let p: Quaternion = From::<Vector>::from(u);
                let v = Vector::from(
                    (&q * p * q.complex_conjugate() + From::<Vector>::from(t)).vector_as_slice(),
                );
                v.into()
            } else {
                (u + t).into()
            }
        })
    }
    /// Transforms a vector of a segment into the OSS
    fn vtov<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<Self, Error>
    where
        Self: Into<Vector>,
        Vector: Into<Self>,
        Segment<M>: SegmentTrait,
    {
        maybe_segment.map(|segment| {
            if let Some(q) = segment.rotation() {
                let u: Vector = self.into();
                let p: Quaternion = From::<Vector>::from(u);
                let v = Vector::from((&q * p * q.complex_conjugate()).vector_as_slice());
                v.into()
            } else {
                self
            }
        })
    }
}
impl Transform for [f64; 3] {}
impl Transform for Vec<f64> {}
/// Mutable geometric transformation with respect to the OSS coordinate system
pub trait TransformMut<'a> {
    /// Transforms the coordinates given in the OSS into a segment
    fn fro<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<(), Error>
    where
        Self: Into<Vector>,
        Segment<M>: SegmentTrait;
    /// Transforms a vector coordinates given in the OSS into a segment
    fn vfrov<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<(), Error>
    where
        Self: Into<Vector>,
        Segment<M>: SegmentTrait;
    /// Transforms the coordinates of a segment into the OSS
    fn to<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<(), Error>
    where
        Self: Into<Vector>,
        Segment<M>: SegmentTrait;
    /// Transforms a segment of a segment into the OSS
    fn vtov<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<(), Error>
    where
        Self: Into<Vector>,
        Segment<M>: SegmentTrait;
}
impl<'a> TransformMut<'a> for &'a mut [f64; 3] {
    /// Transforms the coordinates given in the OSS into a segment
    fn fro<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<(), Error>
    where
        Self: Into<Vector>,
        Segment<M>: SegmentTrait,
    {
        self.to_owned().fro(maybe_segment).map(|v| {
            let _ = std::mem::replace(self, v);
        })
    }
    /// Transforms a vector given in the OSS into a segment
    fn vfrov<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<(), Error>
    where
        Self: Into<Vector>,
        Segment<M>: SegmentTrait,
    {
        self.to_owned().vfrov(maybe_segment).map(|v| {
            let _ = std::mem::replace(self, v);
        })
    }
    /// Transforms the coordinates of a segment into the OSS
    fn to<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<(), Error>
    where
        Self: Into<Vector>,
        Segment<M>: SegmentTrait,
    {
        self.to_owned().to(maybe_segment).map(|v| {
            let _ = std::mem::replace(self, v);
        })
    }
    /// Transforms a vector of a segment into the OSS
    fn vtov<M>(self, maybe_segment: Result<Segment<M>, Error>) -> Result<(), Error>
    where
        Self: Into<Vector>,
        Segment<M>: SegmentTrait,
    {
        self.to_owned().vtov(maybe_segment).map(|v| {
            let _ = std::mem::replace(self, v);
        })
    }
}

#[deprecated(since = "0.2.0", note = "Users should instead use the Transform trait")]
fn conic(rho: f64) -> f64 {
    let c = 36f64;
    let kp = -0.998286;
    let rho2 = rho * rho;
    rho2 / (c + (c * c - (kp + 1f64) * rho2).sqrt())
}

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
            let t = Vector::from([L * c, L * s, 3.9 + conic(L)]);
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
#[deprecated(since = "0.2.0", note = "Users should instead use the Transform trait")]
pub fn m1_any_to_oss<T: Into<Vector>>(sid: usize, p: T) -> Vector {
    use Frame::*;
    let u: Vector = p.into();
    let mut v = Vector::null();
    match sid {
        1 => M1S1(u).to(OSS(&mut v)),
        2 => M1S2(u).to(OSS(&mut v)),
        3 => M1S3(u).to(OSS(&mut v)),
        4 => M1S4(u).to(OSS(&mut v)),
        5 => M1S5(u).to(OSS(&mut v)),
        6 => M1S6(u).to(OSS(&mut v)),
        7 => M1S7(u).to(OSS(&mut v)),
        _ => (),
    }
    v
}
#[deprecated(since = "0.2.0", note = "Users should instead use the Transform trait")]
pub fn oss_to_any_m1<T: Into<Vector>>(sid: usize, p: T) -> Vector {
    use Frame::*;
    let u: Vector = p.into();
    let mut v = Vector::null();
    match sid {
        1 => OSS(u).to(M1S1(&mut v)),
        2 => OSS(u).to(M1S2(&mut v)),
        3 => OSS(u).to(M1S3(&mut v)),
        4 => OSS(u).to(M1S4(&mut v)),
        5 => OSS(u).to(M1S5(&mut v)),
        6 => OSS(u).to(M1S6(&mut v)),
        7 => OSS(u).to(M1S7(&mut v)),
        _ => (),
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

    #[test]
    fn transform_imut_m1_to() {
        for sid in 1..=8 {
            let u = [0.1f64, 0.1, 0.];
            let v = u.to(Segment::<M1>::new(sid));
            println!("M1S{} - v: {:#?}", sid, v);
        }
    }
    #[test]
    fn transform_imut_m1_vtov() {
        for sid in 1..=8 {
            let u = [0.1f64, 0.1, 0.];
            let v = u.vtov(Segment::<M1>::new(sid));
            println!("M1S{} - v: {:#?}", sid, v);
        }
    }
    #[test]
    fn transform_mut_m1_to() {
        for sid in 1..=8 {
            let mut u = [0.1f64, 0.1, 0.];
            (&mut u).to(Segment::<M1>::new(sid)).unwrap();
            println!("M1S{} - u: {:#?}", sid, u);
        }
    }
    #[test]
    fn transform_mut_m1_tofro() {
        for sid in 1..=7 {
            let mut u = [0.1f64, 0.1, 0.];
            (&mut u).to(Segment::<M1>::new(sid)).unwrap();
            (&mut u).fro(Segment::<M1>::new(sid)).unwrap();
            println!("M1S{} - v: {:#?}", sid, u);
        }
    }
    #[test]
    fn transform_m1_tofro() {
        for sid in 1..=7 {
            let u = [0.1f64, 0.1, 0.];
            let v = u
                .to(Segment::<M1>::new(sid))
                .unwrap()
                .fro(Segment::<M1>::new(sid));
            println!("M1S{} - v: {:#?}", sid, v);
        }
    }
    #[test]
    fn transform_m2_to() {
        for sid in 1..=8 {
            let u = [0.1f64, 0.1, 0.];
            let v = u.to(Segment::<M2>::new(sid));
            println!("M2S{} - v: {:#?}", sid, v);
        }
    }
    #[test]
    fn transform_m2_tofro() {
        for sid in 1..=7 {
            let u = [0.1f64, 0.1, 0.];
            let v = u
                .to(Segment::<M2>::new(sid))
                .unwrap()
                .fro(Segment::<M2>::new(sid));
            println!("M1S{} - v: {:#?}", sid, v);
        }
    }
}
