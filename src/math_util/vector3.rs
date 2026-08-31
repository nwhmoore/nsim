/// Three Cartesian components of a vector.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
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
    #[must_use]
    pub fn square(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    /// returns the magnitude of a vector
    #[must_use]
    pub fn norm(&self) -> f64 {
        self.square().sqrt()
    }

    /// rotate the vector around the z-axis
    #[must_use]
    pub fn rotate_z(self, angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();

        Vector3 {
            x: self.x * cos - self.y * sin,
            y: self.x * sin + self.y * cos,
            z: self.z,
        }
    }

    /// rotate the vector around the x-axis
    #[must_use]
    pub fn rotate_x(self, angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();

        Vector3 {
            x: self.x,
            y: self.y * cos - self.z * sin,
            z: self.y * sin + self.z * cos,
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
    pub x: Vec<f64>,
    /// Y components.
    pub y: Vec<f64>,
    /// Z components.
    pub z: Vec<f64>,
}

impl Vector3Series {
    /// fills each series by cloning [`value`]
    pub fn fill(&mut self, value: f64) {
        self.x.fill(value);
        self.y.fill(value);
        self.z.fill(value);
    }

    /// Creates a zero-filled series with one vector entry per index.
    #[must_use]
    pub fn new_with_zeros(length: usize) -> Self {
        Vector3Series {
            x: vec![0.0; length],
            y: vec![0.0; length],
            z: vec![0.0; length],
        }
    }

    /// creates empty series with capacity
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Vector3Series {
            x: Vec::with_capacity(capacity),
            y: Vec::with_capacity(capacity),
            z: Vec::with_capacity(capacity),
        }
    }
    /// Returns the vector stored at `idx`.
    ///
    /// This creates an entirely new structure, DO NOT USE IN FORCE EVALUATION
    /// OR INTEGRATION. For diagnostic and testing API only.
    #[must_use]
    pub fn vector_at(&self, idx: usize) -> Vector3 {
        Vector3 {
            x: self.x[idx],
            y: self.y[idx],
            z: self.z[idx],
        }
    }

    /// Returns the number of stored vectors in the series.
    #[must_use]
    pub fn len(&self) -> usize {
        debug_assert_eq!(self.x.len(), self.y.len());
        debug_assert_eq!(self.x.len(), self.z.len());

        self.x.len()
    }

    /// checks if empty
    #[must_use]
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
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    const EPSILON: f64 = 1e-12;

    static LHS: Vector3 = Vector3 {
        x: 3.0,
        y: 4.0,
        z: 0.0,
    };
    static RHS: Vector3 = Vector3 {
        x: 4.0,
        y: 5.0,
        z: 6.0,
    };

    fn assert_vector_eq(actual: Vector3, expected: Vector3) {
        assert!((actual.x - expected.x).abs() < EPSILON);
        assert!((actual.y - expected.y).abs() < EPSILON);
        assert!((actual.z - expected.z).abs() < EPSILON);
    }

    #[test]
    fn test_add() {
        assert_vector_eq(
            LHS + RHS,
            Vector3 {
                x: 7.0,
                y: 9.0,
                z: 6.0,
            },
        );
    }

    #[test]
    fn test_sub() {
        assert_vector_eq(
            LHS - RHS,
            Vector3 {
                x: -1.0,
                y: -1.0,
                z: -6.0,
            },
        );
    }

    #[test]
    fn test_square() {
        assert_eq!(LHS.square(), 25.0);
    }

    #[test]
    fn test_norm() {
        assert_eq!(LHS.norm(), 5.0);
    }

    #[test]
    fn rotate_z_90_degrees() {
        let vector = Vector3 {
            x: 1.0,
            y: 0.0,
            z: 5.0,
        };

        let result = vector.rotate_z(FRAC_PI_2);

        assert_vector_eq(
            result,
            Vector3 {
                x: 0.0,
                y: 1.0,
                z: 5.0,
            },
        );
    }

    #[test]
    fn rotate_z_preserves_z() {
        let vector = Vector3 {
            x: 3.0,
            y: 4.0,
            z: 7.0,
        };

        let result = vector.rotate_z(1.234);

        assert!((result.z - vector.z).abs() < EPSILON);
    }

    #[test]
    fn rotate_z_zero_degrees() {
        let vector = Vector3 {
            x: 3.0,
            y: 4.0,
            z: 5.0,
        };

        let result = vector.rotate_z(0.0);

        assert_vector_eq(result, vector);
    }

    #[test]
    fn rotate_x_90_degrees() {
        let vector = Vector3 {
            x: 5.0,
            y: 1.0,
            z: 0.0,
        };

        let result = vector.rotate_x(FRAC_PI_2);

        assert_vector_eq(
            result,
            Vector3 {
                x: 5.0,
                y: 0.0,
                z: 1.0,
            },
        );
    }

    #[test]
    fn rotate_x_preserves_x() {
        let vector = Vector3 {
            x: 7.0,
            y: 3.0,
            z: 4.0,
        };

        let result = vector.rotate_x(1.234);

        assert!((result.x - vector.x).abs() < EPSILON);
    }

    #[test]
    fn rotate_x_zero_degrees() {
        let vector = Vector3 {
            x: 3.0,
            y: 4.0,
            z: 5.0,
        };

        let result = vector.rotate_x(0.0);

        assert_vector_eq(result, vector);
    }
}
