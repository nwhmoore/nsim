//! Gravitational acceleration, potential-energy calculation

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
        let massive_indices = state.massive_indices();
        let massless_indices = state.massless_indices();

        for (i, &massive1_idx) in massive_indices.iter().enumerate() {
            let pos_mass1_x = positions.x[massive1_idx];
            let pos_mass1_y = positions.y[massive1_idx];
            let pos_mass1_z = positions.z[massive1_idx];

            let mu1 = GRAVITY * masses[massive1_idx];

            for &massive2_idx in massive_indices.iter().skip(i + 1) {
                let dx = pos_mass1_x - positions.x[massive2_idx];
                let dy = pos_mass1_y - positions.y[massive2_idx];
                let dz = pos_mass1_z - positions.z[massive2_idx];

                let r2 = dx * dx + dy * dy + dz * dz;
                // TODO: make an explicit collision policy
                debug_assert!(
                    r2 > 0.0,
                    "particles {massive1_idx} and {massive2_idx} occupy the same position"
                );
                let inv_r = r2.sqrt().recip();
                let inv_r3 = inv_r * inv_r * inv_r;

                let scale_1 = -GRAVITY * masses[massive2_idx] * inv_r3;
                let scale_2 = mu1 * inv_r3;

                output.accelerations.x[massive1_idx] += dx * scale_1;
                output.accelerations.y[massive1_idx] += dy * scale_1;
                output.accelerations.z[massive1_idx] += dz * scale_1;

                output.accelerations.x[massive2_idx] += dx * scale_2;
                output.accelerations.y[massive2_idx] += dy * scale_2;
                output.accelerations.z[massive2_idx] += dz * scale_2;
            }

            for &test_idx in massless_indices {
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
        let massive_indices = state.massive_indices();

        let mut potential_energy = KahanAccumulator::default();

        for (i, &massive1_idx) in massive_indices.iter().enumerate() {
            for &massive2_idx in massive_indices.iter().skip(i + 1) {
                let dx = positions.x[massive1_idx] - positions.x[massive2_idx];
                let dy = positions.y[massive1_idx] - positions.y[massive2_idx];
                let dz = positions.z[massive1_idx] - positions.z[massive2_idx];

                let r2 = dx * dx + dy * dy + dz * dz;
                // TODO: make an explicit collision policy
                debug_assert!(
                    r2 > 0.0,
                    "particles {massive1_idx} and {massive2_idx} occupy the same position"
                );
                let inv_r = r2.sqrt().recip();

                potential_energy
                    .add(-GRAVITY * masses[massive1_idx] * masses[massive2_idx] * inv_r);
            }
        }

        Some(potential_energy.total())
    }
}
