//! Gravitational acceleration, potential-energy calculation

use crate::{
    force::{PairForceContribution, PairwiseForce},
    math_util::{Geometry, vector3::Vector3},
    particle::ParticleState,
};
use std::f64::consts::PI;

/// Gravitational constant in AU³ · year⁻² · solar-mass⁻¹.
///
/// The units of this constant currently set the units of the entire simulation.
pub const GRAVITY: f64 = 4.0 * PI * PI;

pub struct NewtonianGravity;

impl PairwiseForce for NewtonianGravity {
    fn evaluate_pair(
        &self,
        state: &ParticleState,
        first_idx: usize,
        second_idx: usize,
        geometry: &Geometry,
    ) -> PairForceContribution {
        let first_mass = state.masses()[first_idx];
        let second_mass = state.masses()[second_idx];

        // Gravity
        let scale = -GRAVITY * geometry.inv_dist_cubed();

        let first_acceleration = gravity_acceleration(second_mass, geometry.r_vec(), scale);
        let second_acceleration = gravity_acceleration(first_mass, geometry.r_vec(), -scale);

        let potential_energy = gravitational_potential_energy(first_mass, second_mass, geometry);

        PairForceContribution {
            first_acceleration,
            second_acceleration,
            potential_energy,
        }
    }
}

/// Computes the gravitational acceleration of a particle toward a massive
/// attractor.
///
/// The result is the acceleration vector from the attractor's position to the
/// particle, scaled by the attractor mass and the pair geometry.
pub(super) fn gravity_acceleration(attractor_mass: f64, r_vec: &Vector3, scale: f64) -> Vector3 {
    Vector3 {
        x: r_vec.x * scale * attractor_mass,
        y: r_vec.y * scale * attractor_mass,
        z: r_vec.z * scale * attractor_mass,
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
    -GRAVITY * first_mass * second_mass / geometry.dist()
}
