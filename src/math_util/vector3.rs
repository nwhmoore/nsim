/// Three Cartesian components of a vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vector3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl std::ops::Sub for Vector3 {
    type Output = Vector3;

    fn sub(self, rhs: Vector3) -> Self::Output {
        Vector3 {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl std::ops::Add for Vector3 {
    type Output = Vector3;

    fn add(self, rhs: Vector3) -> Self::Output {
        Vector3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl std::ops::Mul<f64> for Vector3 {
    type Output = Vector3;

    fn mul(self, factor: f64) -> Self::Output {
        Vector3 {
            x: self.x * factor,
            y: self.y * factor,
            z: self.z * factor,
        }
    }
}

impl std::ops::Div<f64> for Vector3 {
    type Output = Vector3;

    fn div(self, factor: f64) -> Self::Output {
        self * (1.0 / factor)
    }
}

impl Vector3 {
    /// Returns the squared Euclidean norm of the vector.
    pub fn square(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Returns the cross product of this vector with `rhs`.
    pub fn cross(&self, rhs: &Vector3) -> Vector3 {
        Vector3 {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }
}

/// Three parallel scalar series representing Cartesian components of a vector.
///
/// The vectors are indexed in lockstep. Depending on the owner, an index may
/// identify a particle or a recorded diagnostic sample.
#[derive(Default, Clone)]
pub struct Vector3Series {
    /// X components.
    x: Vec<f64>,
    /// Y components.
    y: Vec<f64>,
    /// Z components.
    z: Vec<f64>,
}

impl Vector3Series {
    /// Creates a zero-filled series with one vector entry per index.
    pub fn new_zeros(length: usize) -> Self {
        Vector3Series {
            x: vec![0.0; length],
            y: vec![0.0; length],
            z: vec![0.0; length],
        }
    }

    /// Returns the vector stored at `idx`.
    pub fn value_at(&self, idx: usize) -> Vector3 {
        Vector3 {
            x: self.x[idx],
            y: self.y[idx],
            z: self.z[idx],
        }
    }

    /// Stores `value` at the given index.
    pub fn set_value_at(&mut self, idx: usize, value: Vector3) {
        self.x[idx] = value.x;
        self.y[idx] = value.y;
        self.z[idx] = value.z;
    }

    /// Returns the number of stored vectors in the series.
    pub fn len(&self) -> usize {
        debug_assert_eq!(self.x.len(), self.y.len());
        debug_assert_eq!(self.x.len(), self.z.len());

        self.x.len()
    }

    #[allow(unused)]
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    /// Appends one vector value to the end of the series.
    pub fn push(&mut self, vector3: &Vector3) {
        self.x.push(vector3.x);
        self.y.push(vector3.y);
        self.z.push(vector3.z);
    }
}

#[cfg(test)]
mod test {
    use crate::math_util::vector3::Vector3;

    static LHS: Vector3 = Vector3 {
        x: 1.0,
        y: 2.0,
        z: 3.0,
    };
    static RHS: Vector3 = Vector3 {
        x: 4.0,
        y: 5.0,
        z: 6.0,
    };

    #[test]
    fn test_add() {
        assert_eq!(
            LHS + RHS,
            Vector3 {
                x: 5.0,
                y: 7.0,
                z: 9.0,
            }
        )
    }

    #[test]
    fn test_sub() {
        assert_eq!(
            LHS - RHS,
            Vector3 {
                x: -3.0,
                y: -3.0,
                z: -3.0,
            }
        )
    }

    #[test]
    fn test_square() {
        assert_eq!(LHS.square(), 14.0)
    }

    #[test]
    fn test_cross() {
        assert_eq!(
            LHS.cross(&RHS),
            Vector3 {
                x: -3.0,
                y: 6.0,
                z: -3.0,
            }
        )
    }
}
