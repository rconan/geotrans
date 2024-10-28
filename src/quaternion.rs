use crate::Vector;
use std::cmp::PartialEq;
use std::f64::consts::FRAC_PI_2;
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

/// Quaternion
#[derive(Clone, Debug)]
pub struct Quaternion {
    scalar: f64,
    vector: Vector,
}
impl Add for Quaternion {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            scalar: self.scalar + other.scalar,
            vector: self.vector + other.vector,
        }
    }
}
impl Sub for Quaternion {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            scalar: self.scalar - other.scalar,
            vector: self.vector - other.vector,
        }
    }
}
impl Quaternion {
    pub fn new<T: AsRef<[f64]>>(scalar: f64, vector: T) -> Self {
        Self {
            scalar,
            vector: vector.as_ref().into(),
        }
    }
    pub fn pure<T: AsRef<[f64]>>(u: T) -> Self {
        Quaternion::new(0f64, u)
    }
    pub fn unit<T: AsRef<[f64]>>(theta: f64, u: T) -> Self {
        let v = Vector::from(u.as_ref());
        Self {
            scalar: (0.5 * theta).cos(),
            vector: (0.5 * theta).sin() * &v / v.norm(),
        }
    }
    pub fn identity() -> Self {
        Self::new(1f64, [0f64; 3])
    }
    pub fn complex_conjugate(&self) -> Self {
        Self {
            vector: -self.vector.clone(),
            ..*self
        }
    }
    pub fn norm_squared(&self) -> f64 {
        self.scalar * self.scalar + self.vector.norm_squared()
    }
    pub fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }
    pub fn inverse(&self) -> Self {
        self.complex_conjugate() / self.norm_squared()
    }
    pub fn vector_as_slice(&self) -> &[f64] {
        self.vector.as_ref()
    }
    pub fn euler_angles(&self) -> (f64, f64, f64) {
        let w = self.scalar;
        let [x, y, z] = self.vector.clone().into();
        // roll (x-axis rotation)
        let sinr_cosp = 2. * (w * x + y * z);
        let cosr_cosp = 1. - 2. * (x * x + y * y);
        let roll = sinr_cosp.atan2(cosr_cosp);

        // pitch (y-axis rotation)
        let sinp = (1. + 2. * (w * y - x * z)).sqrt();
        let cosp = (1. - 2. * (w * y - x * z)).sqrt();
        let pitch = 2. * sinp.atan2(cosp) - FRAC_PI_2;
        // yaw (z-axis rotation)
        let siny_cosp = 2. * (w * z + x * y);
        let cosy_cosp = 1. - 2. * (y * y + z * z);
        let yaw = siny_cosp.atan2(cosy_cosp);
        (roll, pitch, yaw)
    }
}
impl From<Vector> for Quaternion {
    fn from(v: Vector) -> Self {
        Quaternion::pure(v)
    }
}
impl From<&Vector> for Quaternion {
    fn from(v: &Vector) -> Self {
        Quaternion::pure(v.clone())
    }
}
impl From<&mut Vector> for Quaternion {
    fn from(v: &mut Vector) -> Self {
        Quaternion::pure(v.clone())
    }
}
impl From<&[f64]> for Quaternion {
    fn from(v: &[f64]) -> Self {
        Quaternion::pure(Vector::from(v).clone())
    }
}
impl<T: Into<Quaternion>> Mul<T> for Quaternion {
    type Output = Quaternion;
    fn mul(self, rhs: T) -> Quaternion {
        let rhs: Quaternion = rhs.into();
        Quaternion {
            scalar: self.scalar * rhs.scalar - self.vector.dot(&rhs.vector),
            vector: self.scalar * &rhs.vector
                + rhs.scalar * &self.vector
                + self.vector.cross(&rhs.vector),
        }
    }
}
impl<T: Into<Quaternion>> Mul<T> for &Quaternion {
    type Output = Quaternion;
    fn mul(self, rhs: T) -> Quaternion {
        let rhs: Quaternion = rhs.into();
        Quaternion {
            scalar: self.scalar * rhs.scalar - self.vector.dot(&rhs.vector),
            vector: self.scalar * &rhs.vector
                + rhs.scalar * &self.vector
                + self.vector.cross(&rhs.vector),
        }
    }
}
impl Mul for &Quaternion {
    type Output = Quaternion;
    fn mul(self, rhs: &Quaternion) -> Quaternion {
        Quaternion {
            scalar: self.scalar * rhs.scalar - self.vector.dot(&rhs.vector),
            vector: self.scalar * &rhs.vector
                + rhs.scalar * &self.vector
                + self.vector.cross(&rhs.vector),
        }
    }
}
impl Mul<&Quaternion> for Quaternion {
    type Output = Quaternion;
    fn mul(self, rhs: &Quaternion) -> Quaternion {
        Quaternion {
            scalar: self.scalar * rhs.scalar - self.vector.dot(&rhs.vector),
            vector: self.scalar * &rhs.vector
                + rhs.scalar * &self.vector
                + self.vector.cross(&rhs.vector),
        }
    }
}
impl Div<f64> for &Quaternion {
    type Output = Quaternion;
    fn div(self, rhs: f64) -> Quaternion {
        Quaternion {
            scalar: self.scalar / rhs,
            vector: self.vector.clone() / rhs,
        }
    }
}
impl Div<f64> for Quaternion {
    type Output = Quaternion;
    fn div(self, rhs: f64) -> Quaternion {
        Quaternion {
            scalar: self.scalar / rhs,
            vector: self.vector / rhs,
        }
    }
}
impl Mul<&Quaternion> for f64 {
    type Output = Quaternion;
    fn mul(self, rhs: &Quaternion) -> Quaternion {
        Quaternion {
            scalar: self * rhs.scalar,
            vector: self * &rhs.vector,
        }
    }
}
impl PartialEq for Quaternion {
    fn eq(&self, other: &Self) -> bool {
        self.scalar == other.scalar && self.vector == other.vector
    }
}
impl fmt::Display for Quaternion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:.6} + {}", self.scalar, self.vector)
    }
}

#[cfg(test)]
mod tests {
    use std::f64::EPSILON;

    use super::*;

    #[test]
    fn quaternion_new() {
        let a = Quaternion::new(3., [1., -2., 1.]);
        let v = Quaternion::new(3., Vector::from([1., -2., 1.]));
        println!("a: {}", a);
        println!("v: {}", v);
    }

    #[test]
    fn quaternion_addition() {
        let p = Quaternion::new(3., [1., -2., 1.]);
        let q = Quaternion::new(2., [-1., 2., 3.]);
        let s = p + q;
        assert_eq!(s, Quaternion::new(5., [0., 0., 4.]));
    }
    #[test]
    fn quaternion_multiplication() {
        let p = Quaternion::new(3., [1., -2., 1.]);
        let q = Quaternion::new(2., [-1., 2., 3.]);
        let m = &p * &q;
        assert_eq!(m, Quaternion::new(8., [-9., -2., 11.]));
    }
    #[test]
    fn euler_angles_roll() {
        let q = Quaternion::unit(10f64.to_radians(), Vector::i());
        let (r, p, y) = q.euler_angles();
        println!("roll : {}deg", r.to_degrees());
        println!("pitch: {}deg", p.to_degrees());
        println!("yaw  : {}deg", y.to_degrees());
        assert!((r.to_degrees() - 10f64).abs() < 1e2 * EPSILON)
    }
    #[test]
    fn euler_angles_pitch() {
        let q = Quaternion::unit(-20f64.to_radians(), Vector::j());
        let (r, p, y) = q.euler_angles();
        println!("roll : {}deg", r.to_degrees());
        println!("pitch: {}deg", p.to_degrees());
        println!("yaw  : {}deg", y.to_degrees());
        assert!((p.to_degrees() - -20f64).abs() < 1e2 * EPSILON)
    }
    #[test]
    fn euler_angles_yaw() {
        let q = Quaternion::unit(30f64.to_radians(), Vector::k());
        let (r, p, y) = q.euler_angles();
        println!("roll : {}deg", r.to_degrees());
        println!("pitch: {}deg", p.to_degrees());
        println!("yaw  : {}deg", y.to_degrees());
        assert!((y.to_degrees() - 30f64).abs() < 1e2 * EPSILON)
    }
}
