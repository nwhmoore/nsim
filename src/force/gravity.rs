//! Gravitational acceleration, potential-energy calculation

use crate::{
    force::{Force, ForceEvaluation},
    math_util::kahan::{Kahan3Series, KahanAccumulator},
    particle::ParticleState,
};

pub const GRAVITY: f64 = 1.0;

#[derive(Clone)]
pub struct NewtonianGravity;

impl Force for NewtonianGravity {
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
                debug_assert!(r2 > 0.0, "particles {i} and {j} occupy the same position");
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
            }
        }
    }

    fn calculate_potential_energy(&self, state: &ParticleState) -> Option<f64> {
        let positions = state.positions();
        let masses = state.masses();
        let n = state.particle_count();

        let mut potential_energy = KahanAccumulator::default();

        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions.x[i] - positions.x[j];
                let dy = positions.y[i] - positions.y[j];
                let dz = positions.z[i] - positions.z[j];

                let r2 = dx * dx + dy * dy + dz * dz;
                // TODO: make an explicit collision policy
                debug_assert!(r2 > 0.0, "particles {i} and {j} occupy the same position");
                let inv_r = r2.sqrt().recip();

                potential_energy.add(-GRAVITY * masses[i] * masses[j] * inv_r)
            }
        }

        Some(potential_energy.total())
    }
}

#[derive(Clone)]
pub struct CompensatedNewtonianGravity;

impl Force for CompensatedNewtonianGravity {
    fn evaluate(&self, state: &ParticleState, output: &mut ForceEvaluation<'_>) {
        let positions = state.positions();
        let masses = state.masses();
        let n = state.particle_count();

        let mut accumulator = Kahan3Series::with_len(n);

        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions.x[i] - positions.x[j];
                let dy = positions.y[i] - positions.y[j];
                let dz = positions.z[i] - positions.z[j];

                let r2 = dx * dx + dy * dy + dz * dz;
                // TODO: make an explicit collision policy
                debug_assert!(r2 > 0.0, "particles {i} and {j} occupy the same position");
                let inv_r = r2.sqrt().recip();
                let inv_r3 = inv_r * inv_r * inv_r;

                let scale_i = -GRAVITY * masses[j] * inv_r3;
                let scale_j = GRAVITY * masses[i] * inv_r3;

                accumulator.x[i].add(dx * scale_i);
                accumulator.y[i].add(dy * scale_i);
                accumulator.z[i].add(dz * scale_i);

                accumulator.x[j].add(dx * scale_j);
                accumulator.y[j].add(dy * scale_j);
                accumulator.z[j].add(dz * scale_j);
            }
            output.accelerations.x[i] = accumulator.x[i].total();
            output.accelerations.y[i] = accumulator.y[i].total();
            output.accelerations.z[i] = accumulator.z[i].total();
        }
    }

    fn calculate_potential_energy(&self, state: &ParticleState) -> Option<f64> {
        let positions = state.positions();
        let masses = state.masses();
        let n = state.particle_count();

        let mut potential_energy = KahanAccumulator::default();

        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions.x[i] - positions.x[j];
                let dy = positions.y[i] - positions.y[j];
                let dz = positions.z[i] - positions.z[j];

                let r2 = dx * dx + dy * dy + dz * dz;
                // TODO: make an explicit collision policy
                debug_assert!(r2 > 0.0, "particles {i} and {j} occupy the same position");
                let inv_r = r2.sqrt().recip();

                potential_energy.add(-GRAVITY * masses[i] * masses[j] * inv_r)
            }
        }

        Some(potential_energy.total())
    }
}
