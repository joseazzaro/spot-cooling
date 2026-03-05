# CFD Simulation Integration - Lattice Boltzmann Method (LBM)

## Overview

The Spot Cooling Designer now includes a 2D Lattice Boltzmann Method (LBM) simulator to visualize jet flow physics. This provides:

- **Velocity field visualization** with streamlines and contours
- **Jet decay analysis** (centerline velocity vs axial distance) 
- **Spreading rate estimation** (jet opening angle)
- **Integration with ASHRAE thermal model** for coupled analysis

## Architecture

### Components

1. **`models/jet_lbm.py`**: Core LBM Solver
   - D2Q9 lattice (9 velocities in 2D)
   - BGK collision operator
   - Inlet/outlet boundary conditions
   - Computes velocity, pressure, and density fields

2. **`controllers/cfd_controller.py`**: Simulation Orchestration
   - Creates simulator from cooling calculation parameters
   - Runs time-stepping simulation
   - Post-processes results (centerline decay, spreading)
   - Bridges ASHRAE analytical model with CFD

3. **`views/cfd_visualizer.py`**: Result Visualization
   - Matplotlib-based visualization widget
   - Velocity field contours + streamlines
   - Statistics panel with key metrics
   - PNG export capability

### Integration with Main Application

When user clicks "Run CFD Simulation":
1. Execute normal cooling calculation → get nozzle conditions (T0, V0, etc.)
2. Create LBM simulator with these parameters
3. Run transient simulation for 500 time steps (~5-10 seconds)
4. Display velocity field and analysis metrics
5. Allow export of results

## Physical Model

### Lattice Boltzmann Method (D2Q9)

**Velocities** (in lattice units):
```
        3
        |
    8-0-1
        |
    7-6-5
  2   4
```

**Key Equations:**
- **Collision**: f_i^new = f_i - (f_i - f_eq_i) / τ
- **Streaming**: f_i(x) ← f_i(x - c_i)
- **Macroscopic**: ρ = Σf_i, u = Σf_i c_i / ρ

### Reynolds Number
Based on jet parameters:
```
Re = (V_inlet × D_jet) / ν
```

In simulation:
- Typical Re = 100-1000 (depending on inlet diameter)
- Controlled via kinematic viscosity parameter
- Scales automatically from cooling calculation

### Mach Number
```
Ma = V_inlet / c_sound ≈ 0.01-0.05 for air jets
```
- Kept << 1 for isothermal incompressible flow assumption
- Clamped to M ≤ 0.3 to maintain solver stability

### Domain Size
- Grid: 80 × 50 cells (x, y)
- Aspect ratio ~1.6 (typical for near-field jet analysis)
- Can be adjusted in `LBM2DJetSimulator.__init__`

## Boundary Conditions

| Boundary | Type | Implementation |
|----------|------|-----------------|
| **Inlet (x=0)** | Dirichlet velocity | Gaussian jet profile |
| **Outlet (x=nx-1)** | Zero-gradient | Extrapolation |
| **Sides (y=0, y=ny-1)** | Periodic | Ring buffer (numpy roll) |

### Zou-He Method

Inlet velocity boundary condition uses Zou-He method for stable velocity fixing:
```
- Known: u_inlet (from cooling calc), ρ_inlet = 1.0
- Compute: f_eq at inlet
- Add: non-equilibrium part from adjacent cell for stability
```

## Running CFD Simulation

### From GUI

1. Set simulation parameters (same as cooling analysis)
2. Click "Run & report" (optional, for verification)
3. Click "Run CFD Simulation (2D LBM)" button
4. Progress dialog shows simulation progress
5. Results window displays:
   - Velocity field (contours + streamlines)
   - Statistics: Re, Ma, max velocity
   - Export button for PNG

### Programmatically

```python
from controllers.cfd_controller import create_cfd_from_cooling_solution

# From cooling result
cfd_ctrl = create_cfd_from_cooling_solution(main_window, cooling_result)

# Run and extract results
results = cfd_ctrl.run_simulation(n_steps=500)

# Analyze
centerline = cfd_ctrl.get_jet_centerline_decay()
spreading = cfd_ctrl.get_jet_spreading_rate()
```

## Output Analysis

### Centerline Velocity Decay

Classical jet behavior: U(x)/U₀ ~ constant / r

Where r = (X - X₀)/D is the normalized axial distance.

Code extracts centerline velocity profile and normalizes by inlet velocity.

### Jet Spreading Rate

Measures half-velocity contour width at progression of x positions.

Typical spreading angle:
- Isothermal circular jet: 12-16°
- With buoyancy: 20-24°

CFD provides actual spreading from simulation.

## Performance

Typical execution time on modern hardware:
- 500 time steps: **5-10 seconds**
- Grid size 80×50 = 4000 nodes
- No JIT compilation cache: ~10s first run
- With Numba JIT: ~8s subsequent runs

Bottleneck: Streaming step (9 field operations per timestep)

## Validation against ASHRAE

### Comparison Points

1. **Maximum velocity at inlet**
   - ASHRAE: V₀ from solve_case()
   - CFD: Peak velocity at x=0
   - Should match (inlet BC is Dirichlet)

2. **Velocity decay law**
   - ASHRAE: V_j/V_0 = 1.464/r (simplified)
   - CFD: Actual field data
   - Compare centerline profiles

3. **Spreading rate**
   - ASHRAE: Implicit in acceptable line slope
   - CFD: Explicit from velocity contours
   - Typical: 12-16° for jets

### Known Limitations

- **Isothermal** (no temperature equation yet)
  - Could add energy equation for future work
  - Density variations from temperature ignored
  
- **2D Geometry**
  - Axially symmetric model assumes cylindrical jet
  - Actual 3D effects not captured
  
- **Short simulation**
  - 500 steps reaches pseudo-steady state in near-field
  - Far-field dynamics may not be fully developed

## Future Enhancements

### Phase 1: Energy Equation
- Add temperature equation to LBM
- Multiple relaxation times for thermal LBM (TLBM)
- Compute temperature field alongside velocity
- Validate against ASHRAE thermal predictions

### Phase 2: 3D Simulation
- Extend D2Q9 to D3Q27 (3D, 27 velocities)
- Full axisymmetric domain
- ~20x more computation
- Realistic 3D spreading

### Phase 3: Turbulence Modeling
- Smagorinsky LES for higher Reynolds numbers
- Turbulent kinetic energy transport
- More accurate far-field prediction

### Phase 4: Coupled Analysis
- Interactive CFD during UI parameter adjustment
- Real-time velocity field update
- Visualization of comfort zone relative to jet field

## Files Modified/Added

**New Files:**
- `models/jet_lbm.py` - LBM solver
- `controllers/cfd_controller.py` - CFD orchestration
- `views/cfd_visualizer.py` - Visualization widget
- `CFD_README.md` - This documentation

**Modified Files:**
- `views/main_window.py` - Added CFD button and callback
- `requirements.txt` - Added `numba` dependency

## References

1. Succi, S. (2001). *The Lattice Boltzmann Equation*
2. Boesch, F. et al. (2015). LBM Method Review
3. ASHRAE RP-884 (Original cooling standard)
4. Barba, L., Forsyth, G. CFDPython Course (Educational LBM)

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Import error for `numba` | Run `pip install numba` |
| Slow first run | Numba compiling JIT functions. Subsequent runs faster |
| NaN in results | Check Re is >0, Mach <0.3, inlet conditions stable |
| Memory error on large grids | Reduce nx, ny in `create_default_simulator()` |

## Performance Tuning

For faster execution:
```python
# Reduce grid resolution
sim = LBM2DJetSimulator(nx=60, ny=40, reynolds_number=100)

# Fewer time steps (converges faster at low Re)
results = cfd_ctrl.run_simulation(n_steps=300)
```

For higher accuracy:
```python
# Finer grid (2x slower)
sim = LBM2DJetSimulator(nx=160, ny=100, reynolds_number=100)

# More steps (better steady-state)
results = cfd_ctrl.run_simulation(n_steps=1000)
```
