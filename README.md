# nsim

`nsim` is a small Rust N-body simulation prototype. It currently models massive
bodies and massless test particles using a fixed-timestep
leapfrog/velocity-Verlet integrator.

## Current model

The simulator uses a structure-of-arrays layout:

- `Particle` describes the initial properties of one particle.
- `ParticleSystem` stores particle metadata in a `ParticleCatalog` and numerical
  values in a `ParticleState`.
- `ForceBuffer` stores the Cartesian acceleration of every particle.
- `leapfrog_timestep` advances the state by one fixed timestep.

Particle indices are shared across the catalog, state, force buffer, and output
routines. The vectors for a particle at index `i` must remain aligned.

## Units

By default the gravitational constant is set to `G = 4π²`. This naturally sets
the units of the simulation to astronomical units:

- distance: AU
- time: years
- velocity: AU/year
- acceleration: AU/year²
- mass: solar masses

For other units, change the gravitational constant `GRAVITY` in
[`src/main.rs`](src/main.rs).

## Running the example

Install Rust, then run the binary from the project directory:

```sh
cargo run
```

The example configuration is in [`src/main.rs`](src/main.rs). It creates:

- `Sol` at `(0, 0, 0)` with mass `1.0` solar masses
- `Jupiter` at `(5, 0, 0)` AU with no mass, making it a test particle
- Jupiter's circular-orbit velocity in the positive `y` direction
- one orbital period with approximately 100 integration steps

## Output

At startup, the program creates the `output/` directory if necessary and creates
one file per particle:

```text
output/Sol.out
output/Jupiter.out
```

Existing files with the same particle names are replaced. The initial state at
`t = 0` is recorded before the first integration step, followed by one row for
each timestep. The output format is:

```text
Jupiter
time   x    y    z    u    v    w
0.00000000000000000e0 5.00000000000000000e0 0.00000000000000000e0 ...
```

The values use scientific notation with 17 digits after the decimal point. The
columns are time, position `(x, y, z)`, and velocity `(u, v, w)`, in the units
chosen above.

Global diagnostics are collected in memory by `Diagnostics` while the program
runs.

## Integration method

Each call to `leapfrog_timestep` performs a kick-drift-kick update. The force
buffer must already contain the acceleration at the current positions:

```text
aₙ       = cached acceleration at position xₙ
vₙ₊₁/₂   = vₙ + aₙ · Δt/2
xₙ₊₁     = xₙ + vₙ₊₁/₂ · Δt
aₙ₊₁     = acceleration at position xₙ₊₁
vₙ₊₁     = vₙ₊₁/₂ + aₙ₊₁ · Δt/2
```

The implementation is in [`src/integration.rs`](src/integration.rs). The
force calculation uses

```text
a = -G · M · (x - x_source) / |x - x_source|³
```

Self-interaction is skipped. The initial force evaluation populates the buffer;
each timestep computes the acceleration at the new position and leaves that
value cached for the next timestep. The same force evaluation also returns the
pairwise gravitational potential energy used by the diagnostics.

## Accuracy

Leapfrog/velocity Verlet is second-order accurate for smooth problems, with a
global truncation error that generally scales as `O(Δt²)`. It usually keeps
long-term energy behavior bounded, but it still accumulates phase error. An
orbit can maintain nearly the correct radius while failing to return to its
exact starting coordinates after the analytical period.

The example's 100 steps per orbit are intended for demonstration, not high
precision. Reducing the timestep by a factor of ten generally reduces the
integration error by roughly a factor of 100, at approximately ten times the
runtime cost.
