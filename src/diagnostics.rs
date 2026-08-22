//! Whole-system self recorded as structure-of-arrays time series.
//!
//! Global quantities include only active particles with a `Some` mass value.
//! Massless test particles are intentionally excluded from total mass,
//! kinetic energy, momentum, angular momentum, and center-of-mass quantities.

use crate::{
    math_util::{
        kahan::{Kahan3, KahanAccumulator},
        {Vector3, Vector3Series},
    },
    particle::ParticleState,
};

/// Time series of global quantities derived from simulation states.
///
/// Every field has one entry per call to [`Diagnostics::record`]. The entry at
/// a given index therefore refers to the time at the same index in [`Self::time`].
#[derive(Default, Clone)]
pub struct Diagnostics {
    /// Simulation time associated with each diagnostic sample.
    time: Vec<f64>,

    /// Total mass of active massive bodies.
    total_mass: Vec<f64>,

    /// Total kinetic energy of active massive bodies.
    kinetic_energy: Vec<f64>,
    /// Pairwise Newtonian gravitational potential energy.
    grav_potential_energy: Vec<f64>,
    /// Sum of kinetic and gravitational potential energy.
    total_energy: Vec<f64>,

    /// Total linear momentum, stored as parallel component series.
    linear_momentum: Vector3Series,
    /// Total angular momentum about the simulation origin, stored as parallel
    /// component series.
    angular_momentum: Vector3Series,

    /// Center-of-mass position of the active massive bodies.
    center_of_mass_position: Vector3Series,
    /// Center-of-mass velocity of the active massive bodies.
    center_of_mass_velocity: Vector3Series,
}

impl Diagnostics {
    pub fn number_samples(&self) -> usize {
        self.time.len()
    }

    pub fn linear_momentum(&self) -> &Vector3Series {
        &self.linear_momentum
    }

    pub fn angular_momentum(&self) -> &Vector3Series {
        &self.angular_momentum
    }

    pub fn total_energy(&self) -> &[f64] {
        &self.total_energy
    }

    /// Records one diagnostic sample for a simulation state.
    ///
    /// `potential_energy` must have been evaluated for the same positions in
    /// `state`. It is supplied by the force model so that self do not
    /// duplicate the force-law calculation. If `state` contains no active
    /// massive bodies, the total mass is zero and the center-of-mass values are
    /// computed as division-by-zero results (typically `NaN` or `Inf`).
    pub fn record(&mut self, time: f64, particle_state: &ParticleState, potential_energy: f64) {
        self.time.push(time);

        let mut total_mass = KahanAccumulator::default();

        let mut kinetic_energy = KahanAccumulator::default();

        let mut momentum = Kahan3::default();

        let mut angular_momentum = Kahan3::default();

        let mut mass_position = Kahan3::default();

        for particle_index in 0..particle_state.particle_count() {
            if particle_state.alive_statuses()[particle_index] {
                let position = Vector3 {
                    x: particle_state.positions().x[particle_index],
                    y: particle_state.positions().y[particle_index],
                    z: particle_state.positions().z[particle_index],
                };
                let velocity = Vector3 {
                    x: particle_state.velocities().x[particle_index],
                    y: particle_state.velocities().y[particle_index],
                    z: particle_state.velocities().z[particle_index],
                };
                let mass = particle_state.masses()[particle_index];

                total_mass.add(mass);
                kinetic_energy.add(0.5 * mass * (velocity.square()));

                momentum.add(&(velocity * mass));

                angular_momentum.add(&(position.cross(&velocity) * mass));

                mass_position.add(&(position * mass));
            }
        }

        let all_mass = total_mass.total();
        self.total_mass.push(all_mass);

        let kinetic_energy = kinetic_energy.total();
        self.kinetic_energy.push(kinetic_energy);
        self.grav_potential_energy.push(potential_energy);
        self.total_energy.push(kinetic_energy + potential_energy);

        let total_linear_momentum = momentum.total();
        self.linear_momentum.push(&momentum.total());

        self.angular_momentum.push(&angular_momentum.total());

        self.center_of_mass_position
            .push(&(mass_position.total() / all_mass));

        // v_cm = (sum_i m_i v_i) / (sum_i m_i) = P / M.
        self.center_of_mass_velocity
            .push(&(total_linear_momentum / all_mass));

        debug_assert_eq!(self.number_samples(), self.total_mass.len());
        debug_assert_eq!(self.number_samples(), self.kinetic_energy.len());
        debug_assert_eq!(self.number_samples(), self.grav_potential_energy.len());
        debug_assert_eq!(self.number_samples(), self.total_energy.len());
        debug_assert_eq!(self.number_samples(), self.total_mass.len());
        debug_assert_eq!(self.number_samples(), self.linear_momentum.len());
        debug_assert_eq!(self.number_samples(), self.angular_momentum.len());
        debug_assert_eq!(self.number_samples(), self.center_of_mass_position.len());
        debug_assert_eq!(self.number_samples(), self.center_of_mass_velocity.len());
    }

    pub fn validate_diagnostics(&self) -> std::io::Result<usize> {
        let sample_count = self.number_samples();
        if sample_count == 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "empty self record",
            ));
        }

        for (name, series_len) in [
            ("total_mass", self.total_mass.len()),
            ("kinetic_energy", self.kinetic_energy.len()),
            ("grav_potential_energy", self.grav_potential_energy.len()),
            ("total_energy", self.total_energy.len()),
            ("linear_momentum", self.linear_momentum.len()),
            ("angular_momentum", self.angular_momentum.len()),
            (
                "center_of_mass_position",
                self.center_of_mass_position.len(),
            ),
            (
                "center_of_mass_velocity",
                self.center_of_mass_velocity.len(),
            ),
        ] {
            if series_len != sample_count {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!("self series {name} has {series_len} samples; expected {sample_count}"),
                ));
            }
        }

        Ok(sample_count)
    }

    pub fn diagnostics_values_at(&self, sample_index: usize) -> [f64; 17] {
        [
            self.time[sample_index],
            self.total_mass[sample_index],
            self.kinetic_energy[sample_index],
            self.grav_potential_energy[sample_index],
            self.total_energy[sample_index],
            self.linear_momentum.x[sample_index],
            self.linear_momentum.y[sample_index],
            self.linear_momentum.z[sample_index],
            self.angular_momentum.x[sample_index],
            self.angular_momentum.y[sample_index],
            self.angular_momentum.z[sample_index],
            self.center_of_mass_position.x[sample_index],
            self.center_of_mass_position.y[sample_index],
            self.center_of_mass_position.z[sample_index],
            self.center_of_mass_velocity.x[sample_index],
            self.center_of_mass_velocity.y[sample_index],
            self.center_of_mass_velocity.z[sample_index],
        ]
    }
}
