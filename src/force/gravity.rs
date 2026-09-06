//! Gravitational acceleration and potential-energy calculation.

use crate::{
    force::{Force, ForceEvaluation},
    math_util::KahanAccumulator,
    particle::ParticleState,
};

/// Newtonian gravitational constant. This currently sets the units of the
/// entire simulation.
pub const GRAVITY: f64 = 1.0;

/// Direct Newtonian gravitational force.
#[derive(Clone)]
pub struct NewtonianGravity;

impl Force for NewtonianGravity {
    fn evaluate(&self, state: &ParticleState, output: &mut ForceEvaluation<'_>) {
        let positions = state.positions();
        let masses = state.masses();
        let massive_count = state.massive_count();
        let particle_count = state.particle_count();

        for i in 0..massive_count {
            let pos_mass1_x = positions.x[i];
            let pos_mass1_y = positions.y[i];
            let pos_mass1_z = positions.z[i];

            let mu1 = GRAVITY * masses[i];

            for j in (i + 1)..massive_count {
                let dx = pos_mass1_x - positions.x[j];
                let dy = pos_mass1_y - positions.y[j];
                let dz = pos_mass1_z - positions.z[j];

                let r2 = dx * dx + dy * dy + dz * dz;
                // TODO: make an explicit collision policy
                debug_assert!(r2 > 0.0, "particles {i} and {j} occupy the same position");
                let inv_r = r2.sqrt().recip();
                let inv_r3 = inv_r * inv_r * inv_r;

                let scale_1 = -GRAVITY * masses[j] * inv_r3;
                let scale_2 = mu1 * inv_r3;

                output.accelerations.x[i] += dx * scale_1;
                output.accelerations.y[i] += dy * scale_1;
                output.accelerations.z[i] += dz * scale_1;

                output.accelerations.x[j] += dx * scale_2;
                output.accelerations.y[j] += dy * scale_2;
                output.accelerations.z[j] += dz * scale_2;
            }

            for test_idx in massive_count..particle_count {
                let dx = positions.x[test_idx] - pos_mass1_x;
                let dy = positions.y[test_idx] - pos_mass1_y;
                let dz = positions.z[test_idx] - pos_mass1_z;

                let r2 = dx * dx + dy * dy + dz * dz;
                let inv_r = r2.sqrt().recip();
                let inv_r3 = inv_r * inv_r * inv_r;

                let scale = -mu1 * inv_r3;

                output.accelerations.x[test_idx] += dx * scale;
                output.accelerations.y[test_idx] += dy * scale;
                output.accelerations.z[test_idx] += dz * scale;
            }
        }
    }

    fn calculate_potential_energy(&self, state: &ParticleState) -> Option<f64> {
        let positions = state.positions();
        let masses = state.masses();
        let massive_count = state.massive_count();

        let mut potential_energy = KahanAccumulator::default();

        for i in 0..massive_count {
            for j in (i + 1)..massive_count {
                let dx = positions.x[i] - positions.x[j];
                let dy = positions.y[i] - positions.y[j];
                let dz = positions.z[i] - positions.z[j];

                let r2 = dx * dx + dy * dy + dz * dz;
                // TODO: make an explicit collision policy
                debug_assert!(r2 > 0.0, "particles {i} and {j} occupy the same position");
                let inv_r = r2.sqrt().recip();

                potential_energy.add(-GRAVITY * masses[i] * masses[j] * inv_r);
            }
        }

        Some(potential_energy.total())
    }
}
