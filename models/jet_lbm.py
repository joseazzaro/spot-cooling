# -*- coding: utf-8 -*-
"""
Lattice Boltzmann Method (LBM) 2D Solver for Jet Simulation
Based on D2Q9 velocity lattice model

Physical domain simulation:
- Incompressible isothermal flow (or weakly compressible)
- 2D jet in semi-infinite domain
- Boundary conditions: inlet (jet), outlet (free), symmetry
"""

import numpy as np
from numba import jit
import warnings


class LBM2DJetSimulator:
    """
    2D Lattice Boltzmann Method simulator for air jet.
    
    Uses D2Q9 lattice (9 velocities in 2D):
    Velocities: (0,0), (±1,0), (0,±1), (±1,±1)
    """
    
    def __init__(self, nx=80, ny=50, reynolds_number=100, mach_number=0.1):
        """
        Initialize LBM simulator.
        
        Args:
            nx: Grid size in x direction (axial, downstream)
            ny: Grid size in y direction (radial, perpendicular)
            reynolds_number: Re = V*D/nu (typically 100-1000 for jets)
            mach_number: Ma = V/c_s (Ma < 0.3 for incompressible)
        """
        self.nx = nx
        self.ny = ny
        self.Re = reynolds_number
        self.Ma = mach_number
        
        # D2Q9 lattice vectors (row = velocity index, col = [vx, vy])
        self.c = np.array([
            [0, 0],    # 0: rest
            [1, 0],    # 1: right
            [-1, 0],   # 2: left
            [0, 1],    # 3: up
            [0, -1],   # 4: down
            [1, 1],    # 5: right-up
            [-1, -1],  # 6: left-down
            [1, -1],   # 7: right-down
            [-1, 1],   # 8: left-up
        ], dtype=np.float64)
        
        # Weights for D2Q9
        self.w = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36],
                         dtype=np.float64)
        
        # Speed of sound in lattice units
        self.c_s = 1.0 / np.sqrt(3.0)
        self.c_s2 = self.c_s ** 2
        
        # Physical parameters
        # Jet inlet velocity (lattice units)
        self.u0 = self.Ma * self.c_s
        
        # Kinematic viscosity (lattice units)
        # Re = u0 * D / nu, where D = jet diameter in lattice units
        self.D_lattice = 4.0  # Jet diameter in lattice units
        self.nu = (self.u0 * self.D_lattice) / self.Re
        
        # Relaxation time
        self.tau = 3.0 * self.nu + 0.5
        
        # Distribution functions f[i,y,x] (9 velocities)
        self.f = np.ones((9, self.ny, self.nx), dtype=np.float64)
        self.f_eq = np.ones((9, self.ny, self.nx), dtype=np.float64)
        
        # Macroscopic variables
        self.rho = np.ones((self.ny, self.nx), dtype=np.float64)
        self.ux = np.zeros((self.ny, self.nx), dtype=np.float64)
        self.uy = np.zeros((self.ny, self.nx), dtype=np.float64)
        self.velocity = np.zeros((self.ny, self.nx), dtype=np.float64)
        
        # Initialize with equilibrium
        self._initialize()
    
    def _initialize(self):
        """Initialize distribution functions to equilibrium"""
        self.rho[:] = 1.0
        self.ux[:] = 0.0
        self.uy[:] = 0.0
        self._compute_equilibrium()
        self.f = self.f_eq.copy()
    
    def _compute_equilibrium(self):
        """Compute equilibrium distribution f_eq = w_i * rho * (1 + u·c_i/c_s^2 + (u·c_i)^2/(2*c_s^4) - u^2/(2*c_s^2))"""
        u_sq = self.ux**2 + self.uy**2
        
        for i in range(9):
            cu = self.c[i, 0] * self.ux + self.c[i, 1] * self.uy
            self.f_eq[i] = self.w[i] * self.rho * (
                1.0 + cu / self.c_s2 + 
                0.5 * cu**2 / (self.c_s2**2) - 
                0.5 * u_sq / self.c_s2
            )
    
    def set_inlet_condition(self, jet_center_y, jet_radius, velocity):
        """
        Set inlet boundary condition at x=0 (Dirichlet for velocity).
        
        Args:
            jet_center_y: Center of jet in y direction (lattice units)
            jet_radius: Radius of jet (lattice units)
            velocity: Inlet velocity (lattice units)
        """
        # Parabolic or uniform profile
        y_indices = np.arange(self.ny)
        y_rel = y_indices - jet_center_y
        
        # Gaussian profile (smoother than hard cutoff)
        mask = np.exp(-(y_rel**2) / (2 * jet_radius**2))
        
        self.ux[:, 0] = velocity * mask
        self.uy[:, 0] = 0.0
    
    def collision_step(self):
        """BGK collision operator: f_i^new = f_i - (f_i - f_eq_i) / tau"""
        self._compute_macroscopic()
        self._compute_equilibrium()
        self.f = self.f - (self.f - self.f_eq) / self.tau
    
    def streaming_step(self):
        """Streaming (advection) step: f_i(x) <- f_i(x - c_i dt)"""
        f_new = np.zeros_like(self.f)
        
        for i in range(9):
            cx, cy = self.c[i, 0], self.c[i, 1]
            
            # Periodic boundary in y, outlet in x (extrapolate)
            y_shift = int(cy)
            x_shift = int(cx)
            
            # Shift field
            f_new[i] = np.roll(self.f[i], (-y_shift, -x_shift), axis=(0, 1))
        
        self.f = f_new
    
    def _compute_macroscopic(self):
        """Compute macroscopic variables from distributions: rho = sum(f_i), u = sum(f_i * c_i) / rho"""
        self.rho = np.sum(self.f, axis=0)
        
        # Avoid division by zero
        rho_inv = 1.0 / np.maximum(self.rho, 1e-10)
        
        self.ux = np.sum(self.f * self.c[:, 0:1, np.newaxis, np.newaxis], axis=0) * rho_inv
        self.uy = np.sum(self.f * self.c[:, 1:2, np.newaxis, np.newaxis], axis=0) * rho_inv
        
        self.velocity = np.sqrt(self.ux**2 + self.uy**2)
    
    def apply_boundary_conditions(self):
        """Apply boundary conditions at all borders"""
        # Inlet (x=0): Dirichlet velocity (already set in set_inlet_condition)
        # For inlet nodes, need to reconstruct f from rho and velocity (Zou-He BC)
        self._zou_he_inlet()
        
        # Outlet (x=nx-1): Zero-gradient (extrapolation)
        self.f[:, :, -1] = self.f[:, :, -2]
        
        # Top/Bottom (y): Periodic (numpy roll handles this automatically)
        # or could implement no-slip wall if needed
    
    def _zou_he_inlet(self):
        """Zou-He inlet boundary condition for constant velocity profile"""
        # Only apply to x=0 column
        x = 0
        
        # Macroscopic values at inlet (already set)
        # Need to compute f values from known rho and u at inlet
        
        # For simplicity: use equilibrium + non-equilibrium part from neighbor
        self.rho[:, x] = 1.0  # Assume constant density at inlet
        
        for i in range(9):
            cu = self.c[i, 0] * self.ux[:, x] + self.c[i, 1] * self.uy[:, x]
            u_sq = self.ux[:, x]**2 + self.uy[:, x]**2
            
            f_eq_i = self.w[i] * self.rho[:, x] * (
                1.0 + cu / self.c_s2 + 
                0.5 * cu**2 / (self.c_s2**2) - 
                0.5 * u_sq / self.c_s2
            )
            
            # Use non-eq part from x=1 to ensure smooth inlet
            f_neq_i = self.f[i, :, 1] - self.f_eq[i, :, 1]
            self.f[i, :, x] = f_eq_i + f_neq_i
    
    def step(self):
        """Execute one LBM time step: collision, streaming, boundary conditions"""
        self.collision_step()
        self.streaming_step()
        self.apply_boundary_conditions()
    
    def run_simulation(self, n_steps, jet_center_y=None, jet_radius=None, 
                     jet_velocity=None, callback=None):
        """
        Run simulation for n_steps.
        
        Args:
            n_steps: Number of time steps
            jet_center_y: Center of inlet jet (default: middle of domain)
            jet_radius: Radius of inlet jet (default: 2)
            jet_velocity: Inlet velocity (default: self.u0)
            callback: Function called after each step with progress
        """
        if jet_center_y is None:
            jet_center_y = self.ny / 2.0
        if jet_radius is None:
            jet_radius = 2.0
        if jet_velocity is None:
            jet_velocity = self.u0
        
        self.set_inlet_condition(jet_center_y, jet_radius, jet_velocity)
        
        for step in range(n_steps):
            self.step()
            
            if callback and (step % max(1, n_steps // 10) == 0):
                callback(step, n_steps)
        
        self._compute_macroscopic()  # Final computation
        return self.velocity, self.ux, self.uy, self.rho


def create_default_simulator(reynolds=100, mach=0.1):
    """Create a default LBM simulator with reasonable parameters"""
    return LBM2DJetSimulator(nx=80, ny=50, reynolds_number=reynolds, mach_number=mach)
