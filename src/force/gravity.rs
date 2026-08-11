//! Gravitational acceleration, potential-energy calculation

use crate::utils::{Geometry, Vector3};
use std::f64::consts::PI;

/// Gravitational constant in AU³ · year⁻² · solar-mass⁻¹.
///
/// The units of this constant currently set the units of the entire simulation.
pub const GRAVITY: f64 = 4.0 * PI * PI;

/// Computes the pairwise gravitational acceleration induced on two massive
/// bodies.
///
/// The returned pair `(first, second)` contains the acceleration vectors
/// applied to the first and second body respectively, with equal-and-opposite
/// forces.
pub(super) fn massive_pair_acceleration(
    first_mass: f64,
    second_mass: f64,
    geometry: &Geometry,
) -> (Vector3, Vector3) {
    let scale_from_first = -GRAVITY * first_mass * geometry.inv_dist_cubed;
    let scale_from_second = -GRAVITY * second_mass * geometry.inv_dist_cubed;

    let first = Vector3 {
        x: geometry.r_vec.x * scale_from_second,
        y: geometry.r_vec.y * scale_from_second,
        z: geometry.r_vec.z * scale_from_second,
    };
    let second = Vector3 {
        x: -geometry.r_vec.x * scale_from_first,
        y: -geometry.r_vec.y * scale_from_first,
        z: -geometry.r_vec.z * scale_from_first,
    };

    (first, second)
}

/// Computes the gravitational acceleration of a massless particle toward one
/// massive attractor.
///
/// The result is the acceleration vector from the attractor's position to the
/// massless particle, scaled by the attractor mass and the pair geometry.
pub(super) fn massless_acceleration(attractor_mass: f64, geometry: &Geometry) -> Vector3 {
    let scale = -GRAVITY * attractor_mass * geometry.inv_dist_cubed;

    Vector3 {
        x: geometry.r_vec.x * scale,
        y: geometry.r_vec.y * scale,
        z: geometry.r_vec.z * scale,
    }
}

/// Computes the Newtonian gravitational potential energy of one massive pair.
///
/// The value is derived from the pair separation stored in `geometry`. The
/// caller must only evaluate each active massive pair once and avoid coincident
/// positions, which would otherwise create a singular potential energy.
pub(super) fn gravitational_potential_energy(
    first_mass: f64,
    second_mass: f64,
    geometry: &Geometry,
) -> f64 {
    -GRAVITY * first_mass * second_mass / geometry.dist
}
