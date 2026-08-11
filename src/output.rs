//! Text-file output for particle trajectories and system diagnostics.

use std::{
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
};

use crate::{diagnostics::Diagnostics, particle::ParticleSystem};

const OUTPUT_DIRECTORY: &str = "output";
const PARTICLE_COLUMN_HEADER: &str = "time   x    y    z    u    v    w";
const DIAGNOSTICS_FILE: &str = "diagnostics.out";
const DIAGNOSTICS_COLUMN_HEADER: &str = concat!(
    "time   total_mass   kinetic_energy   grav_potential_energy   total_energy   ",
    "p_x   p_y   p_z   l_x   l_y   l_z   com_x   com_y   com_z   com_u   com_v   com_w"
);

fn output_path(name: &str) -> PathBuf {
    Path::new(OUTPUT_DIRECTORY).join(format!("{name}.out"))
}

fn diagnostics_path() -> PathBuf {
    Path::new(OUTPUT_DIRECTORY).join(DIAGNOSTICS_FILE)
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
pub fn create_particle_file(system: &ParticleSystem, particle_index: usize) -> std::io::Result<()> {
    let name = system.catalog().get_particle_name(particle_index)?;
    fs::create_dir_all(OUTPUT_DIRECTORY)?;

    let mut file = File::create(output_path(name))?;
    writeln!(file, "{name}")?;
    writeln!(file, "{PARTICLE_COLUMN_HEADER}")?;

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
) -> std::io::Result<()> {
    let catalog = system.catalog();
    let state = system.state();
    let name = catalog.get_particle_name(particle_index)?;

    let position = state.positions().value_at(particle_index);
    let velocity = state.velocities().value_at(particle_index);

    let x_pos = position.x;
    let y_pos = position.y;
    let z_pos = position.z;
    let vx = velocity.x;
    let vy = velocity.y;
    let vz = velocity.z;

    if particle_index >= state.positions().len() || particle_index >= state.velocities().len() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("particle index {particle_index} is out of bounds"),
        ));
    }

    let mut file = OpenOptions::new().append(true).open(output_path(name))?;

    writeln!(
        file,
        "{time:.17e} {x_pos:.17e} {y_pos:.17e} {z_pos:.17e} {vx:.17e} {vy:.17e} {vz:.17e}"
    )?;

    Ok(())
}

/// Writes the complete diagnostics report to `output/diagnostics.out`.
///
/// The report starts with initial-state, final-state, and total-change rows,
/// followed by every recorded diagnostic sample. All numeric values use the
/// same scientific-notation precision as the particle trajectory files.
///
/// Every diagnostics series must contain the same non-zero number of samples.
/// Existing diagnostics output is replaced.
///
/// # Errors
///
/// Returns an I/O error if the output directory or file cannot be created, if
/// writing fails, or if the diagnostics series are empty or misaligned.
pub fn write_diagnostics_file(diagnostics: &Diagnostics) -> std::io::Result<()> {
    let sample_count = diagnostics.validate_diagnostics()?;
    let initial = diagnostics.diagnostics_values_at(0);
    let final_values = diagnostics.diagnostics_values_at(sample_count - 1);
    let total_change = std::array::from_fn(|index| final_values[index] - initial[index]);

    fs::create_dir_all(OUTPUT_DIRECTORY)?;
    let mut file = File::create(diagnostics_path())?;

    writeln!(file, "Diagnostics")?;
    writeln!(file, "summary")?;
    writeln!(file, "state   {DIAGNOSTICS_COLUMN_HEADER}")?;
    write_labeled_diagnostics_row(&mut file, "initial_state", &initial)?;
    write_labeled_diagnostics_row(&mut file, "final_state", &final_values)?;
    write_labeled_diagnostics_row(&mut file, "total_change", &total_change)?;

    writeln!(file)?;
    writeln!(file, "record")?;
    writeln!(file, "{DIAGNOSTICS_COLUMN_HEADER}")?;
    for sample_index in 0..sample_count {
        write_diagnostics_values(&mut file, &diagnostics.diagnostics_values_at(sample_index))?;
    }

    Ok(())
}

fn write_labeled_diagnostics_row(
    file: &mut File,
    label: &str,
    values: &[f64; 17],
) -> std::io::Result<()> {
    write!(file, "{label}")?;
    for value in values {
        write!(file, " {value:.17e}")?;
    }
    writeln!(file)
}

fn write_diagnostics_values(file: &mut File, values: &[f64; 17]) -> std::io::Result<()> {
    let (first, rest) = values
        .split_first()
        .expect("diagnostics values are non-empty");
    write!(file, "{first:.17e}")?;
    for value in rest {
        write!(file, " {value:.17e}")?;
    }
    writeln!(file)
}
