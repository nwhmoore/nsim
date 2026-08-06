//! Whole-system diagnostics recorded as structure-of-arrays time series.
//!
//! Global quantities include only active particles with a `Some` mass value.
//! Massless test particles are intentionally excluded from total mass,
//! kinetic energy, momentum, angular momentum, and center-of-mass quantities.

use crate::{
    particle::ParticleState,
    utils::{KahanAccumulator, VectorSeries},
};

/// Time series of global quantities derived from simulation states.
///
/// Every field has one entry per call to [`Diagnostics::record`]. The entry at
/// a given index therefore refers to the time at the same index in [`Self::time`].
#[derive(Default)]
pub struct Diagnostics {
    /// Simulation time associated with each diagnostic sample.
    pub time: Vec<f64>,

    /// Total mass of active massive bodies.
    pub total_mass: Vec<f64>,

    /// Total kinetic energy of active massive bodies.
    pub kinetic_energy: Vec<f64>,
    /// Pairwise Newtonian gravitational potential energy.
    pub grav_potential_energy: Vec<f64>,
    /// Sum of kinetic and gravitational potential energy.
    pub total_energy: Vec<f64>,

    /// Total linear momentum, stored as parallel component series.
    pub linear_momentum: VectorSeries,
    /// Total angular momentum about the simulation origin, stored as parallel
    /// component series.
    pub angular_momentum: VectorSeries,

    /// Center-of-mass position of the active massive bodies.
    pub center_of_mass_position: VectorSeries,
    /// Center-of-mass velocity of the active massive bodies.
    pub center_of_mass_velocity: VectorSeries,
}

impl Diagnostics {
    /// Records one diagnostic sample for a simulation state.
    ///
    /// `potential_energy` must have been evaluated for the same positions in
    /// `state`. It is supplied by the force model so that diagnostics do not
    /// duplicate the force-law calculation. The center-of-mass values are
    /// undefined if `state` contains no active massive bodies.
    pub fn record(&mut self, time: f64, state: &ParticleState, potential_energy: f64) {
        self.time.push(time);

        let mut total_mass = KahanAccumulator::default();

        let mut kinetic_energy = KahanAccumulator::default();

        let mut momentum_x = KahanAccumulator::default();
        let mut momentum_y = KahanAccumulator::default();
        let mut momentum_z = KahanAccumulator::default();

        let mut angular_momentum_x = KahanAccumulator::default();
        let mut angular_momentum_y = KahanAccumulator::default();
        let mut angular_momentum_z = KahanAccumulator::default();

        let mut mass_position_x = KahanAccumulator::default();
        let mut mass_position_y = KahanAccumulator::default();
        let mut mass_position_z = KahanAccumulator::default();

        for particle_index in 0..state.mass.len() {
            if let Some(mass) = state.mass[particle_index]
                && state.alive[particle_index]
            {
                let x = state.position.x[particle_index];
                let y = state.position.y[particle_index];
                let z = state.position.z[particle_index];

                let vx = state.velocity.x[particle_index];
                let vy = state.velocity.y[particle_index];
                let vz = state.velocity.z[particle_index];

                total_mass.add(mass);
                kinetic_energy.add(0.5 * mass * (vx * vx + vy * vy + vz * vz));

                momentum_x.add(mass * vx);
                momentum_y.add(mass * vy);
                momentum_z.add(mass * vz);

                angular_momentum_x.add(mass * (y * vz - z * vy));
                angular_momentum_y.add(mass * (z * vx - x * vz));
                angular_momentum_z.add(mass * (x * vy - y * vx));

                mass_position_x.add(mass * x);
                mass_position_y.add(mass * y);
                mass_position_z.add(mass * z);
            }
        }

        let mass = total_mass.total();
        self.total_mass.push(mass);

        let kinetic_energy = kinetic_energy.total();
        self.kinetic_energy.push(kinetic_energy);
        self.grav_potential_energy.push(potential_energy);
        self.total_energy.push(kinetic_energy + potential_energy);

        let total_linear_momentum_x = momentum_x.total();
        let total_linear_momentum_y = momentum_y.total();
        let total_linear_momentum_z = momentum_z.total();

        self.linear_momentum.x.push(total_linear_momentum_x);
        self.linear_momentum.y.push(total_linear_momentum_y);
        self.linear_momentum.z.push(total_linear_momentum_z);

        self.angular_momentum.x.push(angular_momentum_x.total());
        self.angular_momentum.y.push(angular_momentum_y.total());
        self.angular_momentum.z.push(angular_momentum_z.total());

        self.center_of_mass_position
            .x
            .push(mass_position_x.total() / mass);
        self.center_of_mass_position
            .y
            .push(mass_position_y.total() / mass);
        self.center_of_mass_position
            .z
            .push(mass_position_z.total() / mass);

        // v_cm = (sum_i m_i v_i) / (sum_i m_i) = P / M.
        self.center_of_mass_velocity
            .x
            .push(total_linear_momentum_x / mass);
        self.center_of_mass_velocity
            .y
            .push(total_linear_momentum_y / mass);
        self.center_of_mass_velocity
            .z
            .push(total_linear_momentum_z / mass);
    }
}
