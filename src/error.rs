use std::{error::Error, fmt};

#[derive(Debug)]
#[non_exhaustive]
pub enum SimError {
    /// Simulation cannot be built without a valid integrator
    MissingIntegrator,
}

impl fmt::Display for SimError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SimError::MissingIntegrator => {
                f.write_str("Simulation requires an integrator to be built.")
            }
        }
    }
}

impl Error for SimError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            SimError::MissingIntegrator => None,
        }
    }
}
