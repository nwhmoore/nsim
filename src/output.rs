//! Text-file output for particle trajectories and system diagnostics.

use std::{
    fs::{self, File, OpenOptions},
    io::{self, Write},
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
) -> io::Result<()> {
    let name = particle_name(system, particle_index)?;
    let position = (
        system.state.position.x.get(particle_index),
        system.state.position.y.get(particle_index),
        system.state.position.z.get(particle_index),
    );
    let velocity = (
        system.state.velocity.x.get(particle_index),
        system.state.velocity.y.get(particle_index),
        system.state.velocity.z.get(particle_index),
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
pub fn write_diagnostics_file(diagnostics: &Diagnostics) -> io::Result<()> {
    let sample_count = validate_diagnostics(diagnostics)?;
    let initial = diagnostics_values_at(diagnostics, 0);
    let final_values = diagnostics_values_at(diagnostics, sample_count - 1);
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
        write_diagnostics_values(&mut file, &diagnostics_values_at(diagnostics, sample_index))?;
    }

    Ok(())
}

fn validate_diagnostics(diagnostics: &Diagnostics) -> io::Result<usize> {
    let sample_count = diagnostics.time.len();
    if sample_count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "cannot write an empty diagnostics record",
        ));
    }

    for (name, series_len) in [
        ("total_mass", diagnostics.total_mass.len()),
        ("kinetic_energy", diagnostics.kinetic_energy.len()),
        (
            "grav_potential_energy",
            diagnostics.grav_potential_energy.len(),
        ),
        ("total_energy", diagnostics.total_energy.len()),
        ("linear_momentum.x", diagnostics.linear_momentum.x.len()),
        ("linear_momentum.y", diagnostics.linear_momentum.y.len()),
        ("linear_momentum.z", diagnostics.linear_momentum.z.len()),
        ("angular_momentum.x", diagnostics.angular_momentum.x.len()),
        ("angular_momentum.y", diagnostics.angular_momentum.y.len()),
        ("angular_momentum.z", diagnostics.angular_momentum.z.len()),
        (
            "center_of_mass_position.x",
            diagnostics.center_of_mass_position.x.len(),
        ),
        (
            "center_of_mass_position.y",
            diagnostics.center_of_mass_position.y.len(),
        ),
        (
            "center_of_mass_position.z",
            diagnostics.center_of_mass_position.z.len(),
        ),
        (
            "center_of_mass_velocity.x",
            diagnostics.center_of_mass_velocity.x.len(),
        ),
        (
            "center_of_mass_velocity.y",
            diagnostics.center_of_mass_velocity.y.len(),
        ),
        (
            "center_of_mass_velocity.z",
            diagnostics.center_of_mass_velocity.z.len(),
        ),
    ] {
        if series_len != sample_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "diagnostics series {name} has {series_len} samples; expected {sample_count}"
                ),
            ));
        }
    }

    Ok(sample_count)
}

fn diagnostics_values_at(diagnostics: &Diagnostics, sample_index: usize) -> [f64; 17] {
    [
        diagnostics.time[sample_index],
        diagnostics.total_mass[sample_index],
        diagnostics.kinetic_energy[sample_index],
        diagnostics.grav_potential_energy[sample_index],
        diagnostics.total_energy[sample_index],
        diagnostics.linear_momentum.x[sample_index],
        diagnostics.linear_momentum.y[sample_index],
        diagnostics.linear_momentum.z[sample_index],
        diagnostics.angular_momentum.x[sample_index],
        diagnostics.angular_momentum.y[sample_index],
        diagnostics.angular_momentum.z[sample_index],
        diagnostics.center_of_mass_position.x[sample_index],
        diagnostics.center_of_mass_position.y[sample_index],
        diagnostics.center_of_mass_position.z[sample_index],
        diagnostics.center_of_mass_velocity.x[sample_index],
        diagnostics.center_of_mass_velocity.y[sample_index],
        diagnostics.center_of_mass_velocity.z[sample_index],
    ]
}

fn write_labeled_diagnostics_row(
    file: &mut File,
    label: &str,
    values: &[f64; 17],
) -> io::Result<()> {
    write!(file, "{label}")?;
    for value in values {
        write!(file, " {value:.17e}")?;
    }
    writeln!(file)
}

fn write_diagnostics_values(file: &mut File, values: &[f64; 17]) -> io::Result<()> {
    let (first, rest) = values
        .split_first()
        .expect("diagnostics values are non-empty");
    write!(file, "{first:.17e}")?;
    for value in rest {
        write!(file, " {value:.17e}")?;
    }
    writeln!(file)
}
