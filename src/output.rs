//! Text-file output for particle trajectories.

use std::{
    fs::{self, File, OpenOptions},
    io::{self, Write},
    path::{Path, PathBuf},
};

use crate::particle::ParticleSystem;

const OUTPUT_DIRECTORY: &str = "output";
const COLUMN_HEADER: &str = "time   x    y    z    u    v    w";

fn output_path(name: &str) -> PathBuf {
    Path::new(OUTPUT_DIRECTORY).join(format!("{name}.out"))
}

fn particle_name(system: &ParticleSystem, particle_index: usize) -> io::Result<&str> {
    system
        .catalog
        .name
        .get(particle_index)
        .map(String::as_str)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("particle index {particle_index} is out of bounds"),
            )
        })
}

/// Creates or truncates the output file for one particle.
///
/// The file is written to `output/<particle-name>.out`. The output directory is
/// created automatically when it does not already exist. The first line is the
/// particle name and the second line contains the column headings.
///
/// `particle_index` must refer to an entry in the particle catalog. Existing
/// files with the same particle name are truncated.
///
/// # Errors
///
/// Returns an I/O error if the directory or file cannot be created, or if the
/// particle index is invalid.
pub fn create_particle_file(system: &ParticleSystem, particle_index: usize) -> io::Result<()> {
    let name = particle_name(system, particle_index)?;
    fs::create_dir_all(OUTPUT_DIRECTORY)?;

    let mut file = File::create(output_path(name))?;
    writeln!(file, "{name}")?;
    writeln!(file, "{COLUMN_HEADER}")?;

    Ok(())
}

/// Appends one particle state at the given simulation time.
///
/// Position values are written as `x`, `y`, and `z`; velocity values are
/// written as `u`, `v`, and `w`. Values use scientific notation with enough
/// precision to preserve typical `f64` results when read back.
///
/// `particle_index` must refer to an entry whose position and velocity arrays
/// contain matching entries. The file must have been initialized with
/// [`create_particle_file`] first.
///
/// # Errors
///
/// Returns an I/O error if the particle index or state data is invalid, if the
/// file does not exist, or if the append fails.
pub fn append_particle_timestep(
    system: &ParticleSystem,
    particle_index: usize,
    time: f64,
) -> io::Result<()> {
    let name = particle_name(system, particle_index)?;
    let position = (
        system.state.x.get(particle_index),
        system.state.y.get(particle_index),
        system.state.z.get(particle_index),
    );
    let velocity = (
        system.state.vx.get(particle_index),
        system.state.vy.get(particle_index),
        system.state.vz.get(particle_index),
    );

    let (Some(&x_pos), Some(&y_pos), Some(&z_pos)) = position else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("position data is missing for particle index {particle_index}"),
        ));
    };
    let (Some(&vx), Some(&vy), Some(&vz)) = velocity else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("velocity data is missing for particle index {particle_index}"),
        ));
    };

    let mut file = OpenOptions::new().append(true).open(output_path(name))?;

    writeln!(
        file,
        "{time:.17e} {x_pos:.17e} {y_pos:.17e} {z_pos:.17e} {vx:.17e} {vy:.17e} {vz:.17e}"
    )?;

    Ok(())
}
