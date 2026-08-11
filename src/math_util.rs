//! Shared structure-of-arrays and numerical-accuracy utilities.

use crate::math_util::vector3::Vector3;

pub mod kahan;
pub mod vector3;

/// Pairwise geometric data.
pub struct Geometry {
    /// Relative displacement vector between two particles.
    pub r_vec: Vector3,
    /// Euclidean separation magnitude.
    pub dist: f64,
    /// Inverse cube of the separation magnitude.
    pub inv_dist_cubed: f64,
}

impl Geometry {
    /// Computes the relative geometry for one particle pair.
    pub fn calculate_geometry(r_vec: Vector3) -> Self {
        let dist_squared = r_vec.square();
        let dist = dist_squared.sqrt();
        let inv_dist_cubed = 1.0 / (dist_squared * dist);

        Geometry {
            r_vec,
            dist,
            inv_dist_cubed,
        }
    }
}