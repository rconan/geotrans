use crate::{Gmt, Quaternion, Segment, SegmentTrait, Vector};

/// Geometric transformation with respect to the OSS coordinate system
pub trait Transform {
    /// Transforms the coordinates given in the OSS into a segment
    fn fro<M>(self, segment: Segment<M>) -> Self
    where
        M: Gmt,
        Self: Into<Vector>,
        Vector: Into<Self>,
        Segment<M>: SegmentTrait,
    {
        let u: Vector = self.into();
        let t = segment.translation();
        if let Some(q) = segment.rotation() {
            // let p: Quaternion = From::<Vector>::from(u);
            let v = Vector::from((q.complex_conjugate() * (u - t) * &q).vector_as_slice());
            v.into()
        } else {
            (u - t).into()
        }
    }
    /// Transforms a the vector given in the OSS into a segment
    fn vfrov<M>(self, segment: Segment<M>) -> Self
    where
        M: Gmt,
        Self: Into<Vector>,
        Vector: Into<Self>,
        Segment<M>: SegmentTrait,
    {
        if let Some(q) = segment.rotation() {
            let u: Vector = self.into();
            let p: Quaternion = From::<Vector>::from(u);
            let v = Vector::from((q.complex_conjugate() * p * &q).vector_as_slice());
            v.into()
        } else {
            self
        }
    }
    /// Transforms the coordinates of a segment into the OSS
    fn to<M>(self, segment: Segment<M>) -> Self
    where
        M: Gmt,
        Self: Into<Vector>,
        Vector: Into<Self>,
        Segment<M>: SegmentTrait,
    {
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
    }
    /// Transforms a vector of a segment into the OSS
    fn vtov<M>(self, segment: Segment<M>) -> Self
    where
        M: Gmt,
        Self: Into<Vector>,
        Vector: Into<Self>,
        Segment<M>: SegmentTrait,
    {
        if let Some(q) = segment.rotation() {
            let u: Vector = self.into();
            let p: Quaternion = From::<Vector>::from(u);
            let v = Vector::from((&q * p * q.complex_conjugate()).vector_as_slice());
            v.into()
        } else {
            self
        }
    }
}
impl Transform for [f64; 3] {}
impl Transform for Vec<f64> {}
impl Transform for Vector {}
/// Mutable geometric transformation with respect to the OSS coordinate system
pub trait TransformMut<'a> {
    /// Transforms the coordinates given in the OSS into a segment
    fn fro<M>(self, segment: Segment<M>)
    where
        M: Gmt,
        Self: Into<Vector>,
        Segment<M>: SegmentTrait;
    /// Transforms a vector coordinates given in the OSS into a segment
    fn vfrov<M>(self, segment: Segment<M>)
    where
        M: Gmt,
        Self: Into<Vector>,
        Segment<M>: SegmentTrait;
    /// Transforms the coordinates of a segment into the OSS
    fn to<M>(self, segment: Segment<M>)
    where
        M: Gmt,
        Self: Into<Vector>,
        Segment<M>: SegmentTrait;
    /// Transforms a segment of a segment into the OSS
    fn vtov<M>(self, segment: Segment<M>)
    where
        M: Gmt,
        Self: Into<Vector>,
        Segment<M>: SegmentTrait;
}
impl<'a> TransformMut<'a> for &'a mut [f64; 3] {
    /// Transforms the coordinates given in the OSS into a segment
    fn fro<M>(self, segment: Segment<M>)
    where
        M: Gmt,
        Self: Into<Vector>,
        Segment<M>: SegmentTrait,
    {
        let v = self.to_owned().fro(segment);
        let _ = std::mem::replace(self, v);
    }
    /// Transforms a vector given in the OSS into a segment
    fn vfrov<M>(self, segment: Segment<M>)
    where
        M: Gmt,
        Self: Into<Vector>,
        Segment<M>: SegmentTrait,
    {
        let v = self.to_owned().vfrov(segment);
        let _ = std::mem::replace(self, v);
    }
    /// Transforms the coordinates of a segment into the OSS
    fn to<M>(self, segment: Segment<M>)
    where
        M: Gmt,
        Self: Into<Vector>,
        Segment<M>: SegmentTrait,
    {
        let v = self.to_owned().to(segment);
        let _ = std::mem::replace(self, v);
    }
    /// Transforms a vector of a segment into the OSS
    fn vtov<M>(self, segment: Segment<M>)
    where
        M: Gmt,
        Self: Into<Vector>,
        Segment<M>: SegmentTrait,
    {
        let v = self.to_owned().vtov(segment);
        let _ = std::mem::replace(self, v);
    }
}
