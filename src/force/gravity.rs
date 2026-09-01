//! Gravitational acceleration, potential-energy calculation

use crate::{
    force::{Force, ForceEvaluation},
    math_util::{Kahan3Series, KahanAccumulator},
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
                let inv_r3 = inv_r * inv_r * inv_r;

                let scale_i = -GRAVITY * masses[massive2_idx] * inv_r3;
                let scale_j = GRAVITY * masses[massive1_idx] * inv_r3;

                output.accelerations.x[massive1_idx] += dx * scale_i;
                output.accelerations.y[massive1_idx] += dy * scale_i;
                output.accelerations.z[massive1_idx] += dz * scale_i;

                output.accelerations.x[massive2_idx] += dx * scale_j;
                output.accelerations.y[massive2_idx] += dy * scale_j;
                output.accelerations.z[massive2_idx] += dz * scale_j;
            }
        }

        // preliminary benches suggest to NOT fold the test particle loop into
        // the above `part1` loop
        let massless_indices = state.massless_indices();
        for &test_idx in massless_indices {
            for &massive_idx in massive_indices {
                let dx = positions.x[test_idx] - positions.x[massive_idx];
                let dy = positions.y[test_idx] - positions.y[massive_idx];
                let dz = positions.z[test_idx] - positions.z[massive_idx];

                let r2 = dx * dx + dy * dy + dz * dz;
                // TODO: make an explicit collision policy
                debug_assert!(
                    r2 > 0.0,
                    "particles {test_idx} and {massive_idx} occupy the same position"
                );
                let inv_r = r2.sqrt().recip();
                let inv_r3 = inv_r * inv_r * inv_r;

                let scale = -GRAVITY * masses[massive_idx] * inv_r3;

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

/// Compensated accumulated Newtonian gravity.
#[derive(Clone)]
pub struct CompensatedNewtonianGravity;

impl Force for CompensatedNewtonianGravity {
    fn evaluate(&self, state: &ParticleState, output: &mut ForceEvaluation<'_>) {
        let positions = state.positions();
        let masses = state.masses();
        let massive_indices = state.massive_indices();

        let mut accumulator = Kahan3Series::with_len(state.particle_count());

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
                let inv_r3 = inv_r * inv_r * inv_r;

                let scale_i = -GRAVITY * masses[massive2_idx] * inv_r3;
                let scale_j = GRAVITY * masses[massive1_idx] * inv_r3;

                accumulator.x[massive1_idx].add(dx * scale_i);
                accumulator.y[massive1_idx].add(dy * scale_i);
                accumulator.z[massive1_idx].add(dz * scale_i);

                accumulator.x[massive2_idx].add(dx * scale_j);
                accumulator.y[massive2_idx].add(dy * scale_j);
                accumulator.z[massive2_idx].add(dz * scale_j);
            }
            output.accelerations.x[massive1_idx] = accumulator.x[massive1_idx].total();
            output.accelerations.y[massive1_idx] = accumulator.y[massive1_idx].total();
            output.accelerations.z[massive1_idx] = accumulator.z[massive1_idx].total();
        }

        let massless_indices = state.massless_indices();
        for &test_idx in massless_indices {
            for &source_idx in massive_indices {
                let dx = positions.x[test_idx] - positions.x[source_idx];
                let dy = positions.y[test_idx] - positions.y[source_idx];
                let dz = positions.z[test_idx] - positions.z[source_idx];

                let r2 = dx * dx + dy * dy + dz * dz;
                // TODO: make an explicit collision policy
                debug_assert!(
                    r2 > 0.0,
                    "particles {test_idx} and {source_idx} occupy the same position"
                );
                let inv_r = r2.sqrt().recip();
                let inv_r3 = inv_r * inv_r * inv_r;

                let scale_i = -GRAVITY * masses[source_idx] * inv_r3;

                accumulator.x[test_idx].add(dx * scale_i);
                accumulator.y[test_idx].add(dy * scale_i);
                accumulator.z[test_idx].add(dz * scale_i);
            }
            output.accelerations.x[test_idx] = accumulator.x[test_idx].total();
            output.accelerations.y[test_idx] = accumulator.y[test_idx].total();
            output.accelerations.z[test_idx] = accumulator.z[test_idx].total();
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
