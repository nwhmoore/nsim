//! Particles in the simulation.

use crate::math_util::vector3::{Vector3, Vector3Series};

/// Collection of particles stored as catalog metadata and simulation state.
///
/// The catalog and state use a structure-of-arrays layout. Their vectors must
/// remain aligned: a particle at index `i` has its name and metadata in the
/// catalog's index `i` and its position, velocity, and mass in the state's
/// index `i`.
#[derive(Default)]
pub struct ParticleSystem {
    /// Stable catalog metadata for every particle.
    catalog: ParticleCatalog,
    /// Mutable numerical state used by the integrator.
    state: ParticleState,
    /// Particle ID number to be assigned next.
    next_particle_id: usize,
}

impl ParticleSystem {
    /// Creates an empty particle system.
    pub fn new_system() -> Self {
        ParticleSystem::default()
    }

    pub fn catalog(&self) -> &ParticleCatalog {
        &self.catalog
    }

    pub fn state(&self) -> &ParticleState {
        &self.state
    }

    pub fn state_mut(&mut self) -> &mut ParticleState {
        &mut self.state
    }

    pub fn particle_count(&self) -> usize {
        self.catalog.id.len()
    }

    /// Adds a particle and assigns it the next available catalog ID.
    ///
    /// The particle's metadata and state are appended at the same index in
    /// their respective arrays.
    pub fn new_particle(&mut self, particle: Particle) {
        self.catalog.id.push(self.next_particle_id);
        self.next_particle_id += 1;

        self.catalog.name.push(particle.name);
        self.catalog.radius.push(particle.radius);

        self.state.masses.push(particle.mass);

        self.state.positions.push(&particle.position);

        self.state.velocities.push(&particle.velocity);

        self.state.alive_statuses.push(true);

        debug_assert_eq!(self.particle_count(), self.catalog.name.len());
        debug_assert_eq!(self.particle_count(), self.catalog.radius.len());
        debug_assert_eq!(self.particle_count(), self.state.masses.len());
        debug_assert_eq!(self.particle_count(), self.state.positions.len());
        debug_assert_eq!(self.particle_count(), self.state.velocities.len());
        debug_assert_eq!(self.particle_count(), self.state.alive_statuses.len());
    }
}

/// Persistent metadata associated with each particle.
///
/// Every vector is indexed by the same particle index.
#[derive(Default)]
pub struct ParticleCatalog {
    /// Stable numeric ID assigned when the particle is added.
    id: Vec<usize>,
    /// Particle names, also used as output filename stems.
    name: Vec<String>,
    /// Particle radii.
    radius: Vec<f64>,
}

impl ParticleCatalog {
    pub fn get_particle_name(&self, particle_index: usize) -> std::io::Result<&str> {
        self.name
            .get(particle_index)
            .map(String::as_str)
            .ok_or_else(|| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidInput,
                    format!("particle index {particle_index} is out of bounds"),
                )
            })
    }
}

/// Time-varying numerical state stored for all particles.
///
/// A mass of `None` marks a massless test particle.
#[derive(Default)]
pub struct ParticleState {
    /// Particle masses; `None` denotes a massless test particle.
    masses: Vec<Option<f64>>,

    /// Cartesian positions.
    positions: Vector3Series,

    /// Cartesian velocities.
    velocities: Vector3Series,

    /// Whether each particle is active in the system.
    alive_statuses: Vec<bool>,
}

impl ParticleState {
    pub fn particle_count(&self) -> usize {
        self.alive_statuses.len()
    }

    pub fn masses(&self) -> &[Option<f64>] {
        &self.masses
    }

    pub fn positions(&self) -> &Vector3Series {
        &self.positions
    }

    pub fn positions_mut(&mut self) -> &mut Vector3Series {
        &mut self.positions
    }

    pub fn velocities(&self) -> &Vector3Series {
        &self.velocities
    }

    pub fn velocities_mut(&mut self) -> &mut Vector3Series {
        &mut self.velocities
    }

    pub fn alive_statuses(&self) -> &[bool] {
        &self.alive_statuses
    }
}

/// Initial metadata and state used to add one particle to a [`ParticleSystem`].
pub struct Particle {
    /// Name of the particle, also used as the output filename stem.
    pub name: String,
    /// Radius of the particle.
    pub radius: f64,
    /// Initial position `(x, y, z)`.
    pub position: Vector3,
    /// Initial velocity `(u, v, w)`.
    pub velocity: Vector3,
    /// Mass, or `None` for a massless test particle.
    pub mass: Option<f64>,
}
