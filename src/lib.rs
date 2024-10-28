//! Geometric transformation for the GMT segmented mirrors

mod quaternion;
mod segment;
mod transform;
mod vector;

use std::marker::PhantomData;

pub use quaternion::Quaternion;
pub use segment::{Segment, SegmentTrait};
pub use transform::{Transform, TransformMut};
pub use vector::Vector;

#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("The segment id {0} is not in the range [1,7]")]
    SegmentId(i32),
}

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

/// GMT optics interface
pub trait Gmt {}
impl Gmt for M1 {}
impl Gmt for M2 {}

/// GMT mirrors
pub struct Mirror<M: Gmt>(PhantomData<M>);
impl<M: Gmt> Mirror<M>
where
    Segment<M>: SegmentTrait + Clone,
{
    /// Returns the segment rigid body motions for a mirror tip-tilt
    pub fn tiptilt_2_rigidbodymotions((tip, tilt): (f64, f64)) -> Vec<f64> {
        let q_tt = Quaternion::unit(tip, Vector::i()) * Quaternion::unit(tilt, Vector::j());
        let v7 = Vector::null().to(<Segment<M> as SegmentTrait>::new(7).unwrap());
        let mut rbm = vec![];
        for sid in 1..=7 {
            let segment = <Segment<M> as SegmentTrait>::new(sid).unwrap();

            let v = Vector::null().to(segment.clone()) - &v7;
            let vp = &q_tt * v * q_tt.complex_conjugate();
            let vs = (Vector::from(vp.vector_as_slice()) + &v7).fro(segment.clone());
            rbm.extend(vs);

            let (r, p, y) = if let Some(q_s) = segment.rotation() {
                let q = q_s.complex_conjugate() * &q_tt * q_s;
                q.euler_angles()
            } else {
                q_tt.euler_angles()
            };
            rbm.extend([r, p, y]);
        }
        rbm
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn transform_imut_m1_to() {
        for sid in 1..=7 {
            let u = [0.1f64, 0.1, 0.];
            let v = u.to(Segment::<M1>::new(sid).unwrap());
            println!("M1S{} - v: {:#?}", sid, v);
        }
    }
    #[test]
    fn transform_imut_m1_vtov() {
        for sid in 1..=7 {
            let u = [0.1f64, 0.1, 0.];
            let v = u.vtov(Segment::<M1>::new(sid).unwrap());
            println!("M1S{} - v: {:#?}", sid, v);
        }
    }
    #[test]
    fn transform_mut_m1_to() {
        for sid in 1..=7 {
            let mut u = [0.1f64, 0.1, 0.];
            (&mut u).to(Segment::<M1>::new(sid).unwrap());
            println!("M1S{} - u: {:#?}", sid, u);
        }
    }
    #[test]
    fn transform_mut_m1_tofro() {
        for sid in 1..=7 {
            let mut u = [0.1f64, 0.1, 0.];
            (&mut u).to(Segment::<M1>::new(sid).unwrap());
            (&mut u).fro(Segment::<M1>::new(sid).unwrap());
            println!("M1S{} - v: {:#?}", sid, u);
        }
    }
    #[test]
    fn transform_m1_tofro() {
        for sid in 1..=7 {
            let u = [0.1f64, 0.1, 0.];
            let v = u
                .to(Segment::<M1>::new(sid).unwrap())
                .fro(Segment::<M1>::new(sid).unwrap());
            println!("M1S{} - v: {:#?}", sid, v);
        }
    }
    #[test]
    fn transform_m2_to() {
        for sid in 1..=7 {
            let u = [0.1f64, 0.1, 0.];
            let v = u.to(Segment::<M2>::new(sid).unwrap());
            println!("M2S{} - v: {:#?}", sid, v);
        }
    }
    #[test]
    fn transform_m2_tofro() {
        for sid in 1..=7 {
            let u = [0.1f64, 0.1, 0.];
            let v = u
                .to(Segment::<M2>::new(sid).unwrap())
                .fro(Segment::<M2>::new(sid).unwrap());
            println!("M1S{} - v: {:#?}", sid, v);
        }
    }
    #[test]
    fn global_tiptilt_o() {
        let v7 = Vector::null().to(Segment::<M1>::new(7).unwrap());
        let q = Quaternion::unit(10f64.to_radians(), Vector::i());
        for sid in 1..=7 {
            let segment = Segment::<M1>::new(sid).unwrap();
            let v = Vector::null().to(segment.clone()) - &v7;
            let vp = &q * Quaternion::from(v) * q.complex_conjugate();
            let vs = (Vector::from(vp.vector_as_slice()) + &v7).fro(segment);
            println!("M2S{} - v: {:#?}", sid, vs);
        }
    }
    #[test]
    fn global_tiptilt_ijk() {
        let q = Quaternion::unit(10f64.to_radians(), Vector::i());
        for sid in 1..=7 {
            let segment = Segment::<M1>::new(sid).unwrap();
            let v = Vector::k().vtov(segment.clone());
            let vp = &q * Quaternion::from(v) * q.complex_conjugate();
            let vs = Vector::from(vp.vector_as_slice()).vfrov(segment);
            println!("M2S{} - v: {:#?}", sid, vs);
        }
    }
    #[test]
    fn global_tiptilt() {
        let q_tt = Quaternion::unit(10f64.to_radians(), Vector::j());
        for sid in 1..=7 {
            let segment = Segment::<M1>::new(sid).unwrap();
            if let Some(q_s) = segment.rotation() {
                let q = q_s.complex_conjugate() * &q_tt * q_s;
                dbg!(q.norm());
                let (r, p, y) = q.euler_angles();
                println!(
                    "#{sid}: {:+6.2?}deg",
                    (r.to_degrees(), p.to_degrees(), y.to_degrees())
                );
            } else {
                let (r, p, y) = q_tt.euler_angles();
                println!(
                    "#{sid}: {:+6.2?}deg",
                    (r.to_degrees(), p.to_degrees(), y.to_degrees())
                );
            }
        }
    }
    #[test]
    fn rbm_m1() {
        let rbm =
            Mirror::<M1>::tiptilt_2_rigidbodymotions((1f64.to_radians(), -2.5f64.to_radians()));
        rbm.chunks(6).enumerate().for_each(|(i, rbm)| {
            println!(
                "#{}: {:>+6.1?}{:>+6.1?}",
                i + 1,
                rbm[..3].iter().map(|x| x * 1e3).collect::<Vec<_>>(),
                rbm[3..].iter().map(|x| x.to_degrees()).collect::<Vec<_>>()
            )
        });
    }
    #[test]
    fn rbm_m2() {
        let rbm =
            Mirror::<M2>::tiptilt_2_rigidbodymotions((1f64.to_radians(), -2.5f64.to_radians()));
        rbm.chunks(6).enumerate().for_each(|(i, rbm)| {
            println!(
                "#{}: {:>+6.1?}{:>+6.1?}",
                i + 1,
                rbm[..3].iter().map(|x| x * 1e3).collect::<Vec<_>>(),
                rbm[3..].iter().map(|x| x.to_degrees()).collect::<Vec<_>>()
            )
        });
    }
}
