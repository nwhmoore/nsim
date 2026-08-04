use crate::{GRAVITY, particle::ParticleState};

pub struct ForceBuffer {
    pub ax: Vec<f64>,
    pub ay: Vec<f64>,
    pub az: Vec<f64>,
}

impl ForceBuffer {
    pub fn update_accelerations(&mut self, state: &ParticleState) {
        for target_index in 0..state.mass.len() {
            let mut ax = 0.0;
            let mut ay = 0.0;
            let mut az = 0.0;
            for source_index in 0..state.mass.len() {
                if target_index == source_index {
                    continue;
                }

                let Some(source_mass) = state.mass[source_index] else {
                    continue;
                };

                let dx = state.x[target_index] - state.x[source_index];
                let dy = state.y[target_index] - state.y[source_index];
                let dz = state.z[target_index] - state.z[source_index];

                let dist_squared = dx * dx + dy * dy + dz * dz;
                ax += gravity_acceleration(state.x[target_index], dist_squared, source_mass);
                ay += gravity_acceleration(state.y[target_index], dist_squared, source_mass);
                az += gravity_acceleration(state.z[target_index], dist_squared, source_mass);
            }
            self.ax[target_index] = ax;
            self.ay[target_index] = ay;
            self.az[target_index] = az;
        }
    }
}

pub fn gravity_acceleration(dimension_dist: f64, dist_squared: f64, attractor_mass: f64) -> f64 {
    // fvec = m1 avec = g m1 m2 / rmag^3 rvec

    -GRAVITY * attractor_mass * dimension_dist / (dist_squared * dist_squared.sqrt())
}
