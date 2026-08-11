//! Particles in the simulation.

use crate::utils::{Vector3, Vector3Series};

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
                position: Vector3Series {
                    x: Vec::new(),
                    y: Vec::new(),
                    z: Vec::new(),
                },
                velocity: Vector3Series {
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

        self.state.position.x.push(particle.pos.x);
        self.state.position.y.push(particle.pos.y);
        self.state.position.z.push(particle.pos.z);

        self.state.velocity.x.push(particle.vel.x);
        self.state.velocity.y.push(particle.vel.y);
        self.state.velocity.z.push(particle.vel.z);

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
    pub position: Vector3Series,

    /// Cartesian velocities.
    pub velocity: Vector3Series,

    /// Whether each particle is active in the system.
    pub alive: Vec<bool>,
}

/// Initial metadata and state used to add one particle to a [`ParticleSystem`].
pub struct Particle {
    /// Name of the particle, also used as the output filename stem.
    pub name: String,
    /// Radius of the particle.
    pub radius: f64,
    /// Initial position `(x, y, z)`.
    pub pos: Vector3,
    /// Initial velocity `(u, v, w)`.
    pub vel: Vector3,
    /// Mass, or `None` for a massless test particle.
    pub mass: Option<f64>,
}
