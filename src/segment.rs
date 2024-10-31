//! Geometric transformation for the GMT segmented mirrors

use std::marker::PhantomData;

use crate::{Conic, Error, Gmt, Quaternion, Vector, M1, M2};

/// GMT segmented mirror
#[derive(Debug, Clone)]
pub struct Segment<M: Gmt> {
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
impl<M: Gmt> Segment<M> {
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
            Some(
                Quaternion::unit(180f64.to_radians(), Vector::k())
                    * Quaternion::unit(std::f64::consts::PI, Vector::j()),
            )
        }
    }
}
