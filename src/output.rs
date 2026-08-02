//! Text-file output for body trajectories.

use std::{
    fs::{self, File, OpenOptions},
    io::{self, Write},
    path::{Path, PathBuf},
};

use crate::{
    body::{LargeBody, SmallBody},
    math::{Position, Velocity, time::Time},
};

const OUTPUT_DIRECTORY: &str = "output";

/// Column headings written to every trajectory file.
const COLUMN_HEADER: &str =
    "time(yr)   x(AU)    y(AU)    z(AU)    u(AU/yr)    v(AU/yr)    w(AU/yr)";

/// The state fields required to write a body's trajectory.
pub trait OutputBody {
    /// Returns the body's name used in the first line and filename.
    fn output_name(&self) -> &str;

    /// Returns the body's current position.
    fn output_position(&self) -> Position;

    /// Returns the body's current velocity.
    fn output_velocity(&self) -> Velocity;
}

impl OutputBody for LargeBody {
    fn output_name(&self) -> &str {
        &self.name
    }

    fn output_position(&self) -> Position {
        self.pos
    }

    fn output_velocity(&self) -> Velocity {
        self.vel
    }
}

impl OutputBody for SmallBody {
    fn output_name(&self) -> &str {
        &self.name
    }

    fn output_position(&self) -> Position {
        self.pos
    }

    fn output_velocity(&self) -> Velocity {
        self.vel
    }
}

fn output_path(name: &str) -> PathBuf {
    Path::new(OUTPUT_DIRECTORY).join(format!("{name}.out"))
}

/// Creates or truncates a body's output file and writes its headers.
///
/// The file is created at `output/<body-name>.out`. The `output` directory is
/// created automatically when needed.
pub fn create_body_file<B>(body: &B) -> io::Result<()>
where
    B: OutputBody,
{
    fs::create_dir_all(OUTPUT_DIRECTORY)?;

    let mut file = File::create(output_path(body.output_name()))?;
    writeln!(file, "{}", body.output_name())?;
    writeln!(file, "{COLUMN_HEADER}")?;

    Ok(())
}

/// Appends one time-step row to a body's output file.
///
/// The row contains time, position, and velocity in the units documented by
/// the column header. The file must have been initialized with
/// [`create_body_file`] first.
pub fn append_body_timestep<B>(body: &B, time: Time) -> io::Result<()>
where
    B: OutputBody,
{
    let mut file = OpenOptions::new()
        .append(true)
        .open(output_path(body.output_name()))?;

    let position = body.output_position();
    let velocity = body.output_velocity();

    writeln!(
        file,
        "{:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}",
        time.as_years(),
        position.x,
        position.y,
        position.z,
        velocity.u,
        velocity.v,
        velocity.w,
    )?;

    Ok(())
}
