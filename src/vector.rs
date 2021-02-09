use std::cmp::PartialEq;
use std::fmt;
use std::ops::{Add, AddAssign, Deref, Div, Mul, Neg, Sub};
/// Vector
#[derive(Clone, Debug)]
pub struct Vector([f64; 3]);
impl Vector {
    pub fn dot(&self, other: &Vector) -> f64 {
        self.0
            .iter()
            .zip(other.0.iter())
            .fold(0., |a, (x, y)| a + x * y)
    }
    pub fn cross(&self, other: &Vector) -> Vector {
        let [a1, a2, a3] = self.0;
        let [b1, b2, b3] = other.0;
        Vector([a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1])
    }
    pub fn norm_squared(&self) -> f64 {
        self.dot(self)
    }
    pub fn norm(&self) -> f64 {
        self.norm_squared().sqrt()
    }
    pub fn null() -> Self {
        Vector::from(0f64)
    }
    pub fn i() -> Self {
        Vector::from([1, 0, 0])
    }
    pub fn j() -> Self {
        Vector::from([0, 1, 0])
    }
    pub fn k() -> Self {
        Vector::from([0, 0, 1])
    }
}
impl AsRef<[f64]> for Vector {
    fn as_ref(&self) -> &[f64] {
        &self.0
    }
}
impl Deref for Vector {
    type Target = [f64];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
impl From<Vec<f64>> for Vector {
    fn from(v: Vec<f64>) -> Self {
        Vector([v[0], v[1], v[2]])
    }
}
impl From<f64> for Vector {
    fn from(v: f64) -> Self {
        Vector([v; 3])
    }
}
impl From<std::slice::Iter<'_, f64>> for Vector {
    fn from(v: std::slice::Iter<'_, f64>) -> Self {
        Vector::from(v.cloned().collect::<Vec<f64>>())
    }
}
impl From<[f64; 3]> for Vector {
    fn from(v: [f64; 3]) -> Self {
        Vector(v)
    }
}
impl From<[i32; 3]> for Vector {
    fn from(v: [i32; 3]) -> Self {
        Vector::from(v.iter().map(|&x| x as f64).collect::<Vec<f64>>())
    }
}
impl From<&[f64]> for Vector {
    fn from(v: &[f64]) -> Self {
        if v.len() == 2 {
            Vector([v[0], v[1], 0.])
        } else {
            Vector([v[0], v[1], v[2]])
        }
    }
}
impl Add for Vector {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Vector([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
        ])
    }
}
impl AddAssign for Vector {
    fn add_assign(&mut self, other: Self) {
        *self = Self([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
        ])
    }
}
impl Add<&Vector> for Vector {
    type Output = Self;
    fn add(self, other: &Self) -> Self {
        Vector([
            self.0[0] + other.0[0],
            self.0[1] + other.0[1],
            self.0[2] + other.0[2],
        ])
    }
}
impl Sub for Vector {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Vector([
            self.0[0] - other.0[0],
            self.0[1] - other.0[1],
            self.0[2] - other.0[2],
        ])
    }
}
impl Sub<&Vector> for Vector {
    type Output = Self;
    fn sub(self, other: &Self) -> Self {
        Vector([
            self.0[0] - other.0[0],
            self.0[1] - other.0[1],
            self.0[2] - other.0[2],
        ])
    }
}
impl Neg for Vector {
    type Output = Self;
    fn neg(self) -> Self {
        Vector([-self.0[0], -self.0[1], -self.0[2]])
    }
}
impl Mul<f64> for &Vector {
    type Output = Vector;
    fn mul(self, rhs: f64) -> Vector {
        Vector([
            rhs * self.as_ref()[0],
            rhs * self.as_ref()[1],
            rhs * self.as_ref()[2],
        ])
    }
}
impl Div<f64> for &Vector {
    type Output = Vector;
    fn div(self, rhs: f64) -> Vector {
        Vector([
            self.as_ref()[0] / rhs,
            self.as_ref()[1] / rhs,
            self.as_ref()[2] / rhs,
        ])
    }
}
impl Div<f64> for Vector {
    type Output = Vector;
    fn div(self, rhs: f64) -> Vector {
        Vector([
            self.as_ref()[0] / rhs,
            self.as_ref()[1] / rhs,
            self.as_ref()[2] / rhs,
        ])
    }
}
impl Mul<&Vector> for f64 {
    type Output = Vector;
    fn mul(self, rhs: &Vector) -> Vector {
        Vector([
            rhs.as_ref()[0] * self,
            rhs.as_ref()[1] * self,
            rhs.as_ref()[2] * self,
        ])
    }
}
impl Mul<Vector> for f64 {
    type Output = Vector;
    fn mul(self, rhs: Vector) -> Vector {
        Vector([
            rhs.as_ref()[0] * self,
            rhs.as_ref()[1] * self,
            rhs.as_ref()[2] * self,
        ])
    }
}
impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        self.0
            .iter()
            .zip(other.0.iter())
            .fold(true, |a, (x, y)| a && x == y)
    }
}
impl fmt::Display for Vector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:.6}i + {:.6}j + {:.6}k",
            self.0[0], self.0[1], self.0[2]
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vector_dot() {
        let u = Vector::from([1., -2., 1.]);
        let v = Vector::from([-1., 2., 3.]);
        let d = u.dot(&v);
        assert_eq!(d, -2.)
    }

    #[test]
    fn vector_cross() {
        let u = Vector::from([1., -2., 1.]);
        let v = Vector::from([-1., 2., 3.]);
        let c = u.cross(&v);
        assert_eq!(c, Vector::from([-8., -4., 0.]));
    }

    #[test]
    fn vector_scale() {
        let u = Vector::from([1., -2., 1.]);
        let s = 3. * u;
        assert_eq!(s, Vector::from([3., -6., 3.]));
    }
}
