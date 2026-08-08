//! Particles in the simulation.

use crate::utils::VectorSeries;

/// Collection of particles stored as catalog metadata and simulation state.
///
/// The catalog and state use a structure-of-arrays layout. Their vectors must
/// remain aligned: a particle at index `i` has its name and metadata in the
/// catalog's index `i` and its position, velocity, and mass in the state's
/// index `i`.
pub struct ParticleSystem {
    /// Stable catalog metadata for every particle.
    pub catalog: ParticleCatalog,
    /// Mutable numerical state used by the integrator.
    pub state: ParticleState,
    next_particle_id: usize,
}

impl ParticleSystem {
    /// Creates an empty particle system.
    pub fn new_system() -> Self {
        ParticleSystem {
            catalog: ParticleCatalog {
                id: Vec::new(),
                name: Vec::new(),
                radius: Vec::new(),
            },
            state: ParticleState {
                mass: Vec::new(),
                position: VectorSeries {
                    x: Vec::new(),
                    y: Vec::new(),
                    z: Vec::new(),
                },
                velocity: VectorSeries {
                    x: Vec::new(),
                    y: Vec::new(),
                    z: Vec::new(),
                },
                alive: Vec::new(),
            },
            next_particle_id: 0,
        }
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

        self.state.mass.push(particle.mass);

        self.state.position.x.push(particle.pos.0);
        self.state.position.y.push(particle.pos.1);
        self.state.position.z.push(particle.pos.2);

        self.state.velocity.x.push(particle.vel.0);
        self.state.velocity.y.push(particle.vel.1);
        self.state.velocity.z.push(particle.vel.2);

        self.state.alive.push(true);
    }
}

/// Persistent metadata associated with each particle.
///
/// Every vector is indexed by the same particle index.
pub struct ParticleCatalog {
    /// Stable numeric ID assigned when the particle is added.
    pub id: Vec<usize>,
    /// Particle names, also used as output filename stems.
    pub name: Vec<String>,
    /// Particle radii.
    pub radius: Vec<f64>,
}

/// Time-varying numerical state stored for all particles.
///
/// A mass of `None` marks a massless test particle.
pub struct ParticleState {
    /// Particle masses; `None` denotes a massless test particle.
    pub mass: Vec<Option<f64>>,

    /// Cartesian positions.
    pub position: VectorSeries,

    /// Cartesian velocities.
    pub velocity: VectorSeries,

    /// Whether each particle is active in the system.
    pub alive: Vec<bool>,
}

impl ParticleState {
    pub fn alive_particle_groups(&self) -> (Vec<(usize, f64)>, Vec<usize>) {
        let mut massive = Vec::with_capacity(self.alive.len());
        let mut massless = Vec::with_capacity(self.alive.len());

        for particle_idx in 0..self.alive.len() {
            if !self.alive[particle_idx] {
                continue;
            }

            if let Some(mass) = self.mass[particle_idx] {
                massive.push((particle_idx, mass));
            } else {
                massless.push(particle_idx);
            }
        }
        (massive, massless)
    }
}

/// Initial metadata and state used to add one particle to a [`ParticleSystem`].
pub struct Particle {
    /// Name of the particle, also used as the output filename stem.
    pub name: String,
    /// Radius of the particle.
    pub radius: f64,
    /// Initial position `(x, y, z)`.
    pub pos: (f64, f64, f64),
    /// Initial velocity `(u, v, w)`.
    pub vel: (f64, f64, f64),
    /// Mass, or `None` for a massless test particle.
    pub mass: Option<f64>,
}
