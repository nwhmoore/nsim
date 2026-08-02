# nsim

`nsim` is a small Rust N-body simulation prototype using astronomical units:

- distance: AU
- time: years
- mass: solar masses
- velocity: AU/year
- acceleration: AU/year²

The gravitational constant is set to `4π²`, which makes Kepler's third law take
the convenient form

```text
T² = a³ / M
```

where `T` is in years, `a` is in AU, and `M` is in solar masses.

## Usage

Install Rust, then run the simulation from the project directory:

```sh
cargo run
```

The initial conditions and timestep are currently configured in
[`src/main.rs`](src/main.rs). The example places a one-solar-mass `Sol` at the
origin and a test-particle `Jupiter` at 5 AU with the circular-orbit velocity.
The run uses one orbital period and 100 integration steps.

The simulation creates or replaces an `output/` directory and writes one file
per body:

```text
output/Sol.out
output/Jupiter.out
```

Each file contains the body name, a column header, and one row per recorded
state, including the initial state at `t = 0`:

```text
Jupiter
time(yr)   x(AU)    y(AU)    z(AU)    u(AU/yr)    v(AU/yr)    w(AU/yr)
0.00000000000000000e0 5.00000000000000000e0 0.00000000000000000e0 ...
```

Repeated runs overwrite files with the same body names. Body names should be
unique and suitable for use as filenames.

## Integration method

Each timestep uses a kick-drift-kick leapfrog scheme, also commonly described
as velocity Verlet:

```text
aₙ       = acceleration at position xₙ
vₙ₊₁/₂   = vₙ + aₙ · Δt/2
xₙ₊₁     = xₙ + vₙ₊₁/₂ · Δt
aₙ₊₁     = acceleration at position xₙ₊₁
vₙ₊₁     = vₙ₊₁/₂ + aₙ₊₁ · Δt/2
```

The implementation is in [`src/simulation.rs`](src/simulation.rs). Gravity
from a source body is computed as

```text
a = -G · M · (x - x_source) / |x - x_source|³
```

Self-interaction is excluded for `LargeBody` values. A temporary acceleration
buffer is used so all large-body accelerations are calculated from the same
positions before they are written back.

## Accuracy and limitations

The method is second-order accurate: for a fixed smooth problem, the global
integration error is generally `O(Δt²)`. Leapfrog/velocity Verlet usually has
good long-term energy behavior, but its numerical orbital frequency differs
slightly from the continuous solution. Consequently, an orbit can preserve
its radius well while accumulating a phase error and not land exactly on its
initial coordinates after the analytical period.

The current implementation has these modeling limitations:

- `SmallBody` is a massless test particle. It feels gravity from all large
  bodies but does not exert gravity back on them or on other small bodies.
- `LargeBody` values interact with one another, while self-gravity is skipped.
- The timestep is fixed; there is no adaptive stepping or error controller.
- Very close approaches and collisions are not regularized. As the separation
  approaches zero, the inverse-cube force becomes singular.
- The simulation loop can advance slightly past `time_end` when the requested
  interval is not an exact multiple of the timestep.
- The current program is configured by editing `src/main.rs`; it does not yet
  provide a command-line interface or input-file format.

## License

This project is licensed under the Apache License, Version 2.0. See
[`LICENSE`](LICENSE) for the complete license text.
