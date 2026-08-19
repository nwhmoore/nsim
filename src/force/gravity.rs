//! Gravitational acceleration, potential-energy calculation

use crate::{
    force::{ForceEvaluation, PairwiseForce},
    particle::ParticleState,
};
use std::f64::consts::PI;

/// Gravitational constant in AU³ · year⁻² · solar-mass⁻¹.
///
/// The units of this constant currently set the units of the entire simulation.
pub const GRAVITY: f64 = 4.0 * PI * PI;

pub struct NewtonianGravity;

impl PairwiseForce for NewtonianGravity {
    fn evaluate(&self, state: &ParticleState, output: &mut ForceEvaluation<'_>) {
        let positions = state.positions();
        let masses = state.masses();
        let n = state.particle_count();

        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions.x[i] - positions.x[j];
                let dy = positions.y[i] - positions.y[j];
                let dz = positions.z[i] - positions.z[j];

                let r2 = dx * dx + dy * dy + dz * dz;
                // TODO: make an explicit collision policy
                assert!(r2 > 0.0, "particles {i} and {j} occupy the same position");
                let inv_r = r2.sqrt().recip();
                let inv_r3 = inv_r * inv_r * inv_r;

                let scale_i = -GRAVITY * masses[j] * inv_r3;
                let scale_j = GRAVITY * masses[i] * inv_r3;

                output.accelerations.x[i] += dx * scale_i;
                output.accelerations.y[i] += dy * scale_i;
                output.accelerations.z[i] += dz * scale_i;

                output.accelerations.x[j] += dx * scale_j;
                output.accelerations.y[j] += dy * scale_j;
                output.accelerations.z[j] += dz * scale_j;

                output
                    .potential_energy
                    .add(-GRAVITY * masses[i] * masses[j] * inv_r)
            }
        }
    }
}
