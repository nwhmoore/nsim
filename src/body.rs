use crate::math::{Acceleration, Position, Velocity};

/// Large body which exerts gravity on other objects
pub struct LargeBody {
    /// Name of object
    pub name: String,
    /// Position
    pub pos: Position,
    /// Velocity
    pub vel: Velocity,
    pub acc: Acceleration,
    /// Mass of `LargeBody` in units of solar mass
    pub mass: f64,
}

pub struct SmallBody {
    pub name: String,
    pub pos: Position,
    pub vel: Velocity,
    pub acc: Acceleration,
}

pub trait Positioned {
    fn position(&self) -> &Position;
}

impl Positioned for SmallBody {
    fn position(&self) -> &Position {
        &self.pos
    }
}

impl Positioned for LargeBody {
    fn position(&self) -> &Position {
        &self.pos
    }
}
