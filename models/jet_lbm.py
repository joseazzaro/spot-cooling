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
    
    def set_inlet_condition(self, jet_center_x, jet_radius, velocity):
        """
        Set inlet boundary condition at y=0 (ceiling, Dirichlet for velocity).
        Jet enters from top (y=0) and flows downward (uy < 0).
        
        Args:
            jet_center_x: Center of jet in x direction (lattice units, typically nx/2 for centerline)
            jet_radius: Radius of jet opening (lattice units)
            velocity: Inlet velocity magnitude (lattice units, flows downward as negative uy)
        """
        # Gaussian profile for smooth jet opening
        x_indices = np.arange(self.nx)
        x_rel = x_indices - jet_center_x
        
        # Gaussian radial profile (smooth nozzle)
        mask = np.exp(-(x_rel**2) / (2 * jet_radius**2))
        
        # Jet enters at y=0 with downward velocity (uy = +velocity)
        # No radial velocity at entrance
        self.ux[0, :] = 0.0
        self.uy[0, :] = velocity * mask  # Positive = downward (y increases downward)
    
    def collision_step(self):
        """BGK collision operator: f_i^new = f_i - (f_i - f_eq_i) / tau"""
        self._compute_macroscopic()
        self._compute_equilibrium()
        self.f = self.f - (self.f - self.f_eq) / self.tau
    
    def streaming_step(self):
        """
        Streaming (advection) step: f_i(x,t+1) <- f_i(x - c_i, t)
        - X direction: symmetric (not periodic) - reflect at x=0
        - Y direction: non-periodic (inlet top, outlet bottom)
        """
        f_new = np.zeros_like(self.f)
        
        for i in range(9):
            cx, cy = self.c[i, 0], self.c[i, 1]
            x_shift = int(cx)
            y_shift = int(cy)
            
            # Y dimension: standard non-periodic streaming (top to bottom flow)
            if y_shift > 0:
                f_shifted_y = self.f[i, :-y_shift, :]
                f_temp = np.zeros((self.ny, self.nx))
                f_temp[y_shift:, :] = f_shifted_y
            elif y_shift < 0:
                y_abs = abs(y_shift)
                f_shifted_y = self.f[i, y_abs:, :]
                f_temp = np.zeros((self.ny, self.nx))
                f_temp[:-y_abs, :] = f_shifted_y
            else:
                f_temp = self.f[i, :, :]
            
            # X dimension: symmetric (not periodic) - use mirror reflection
            if x_shift > 0:  # Moving right (+x direction)
                # Shift right: copy from left
                f_new[i, :, x_shift:] = f_temp[:, :-x_shift]
                # Left boundary (x < x_shift): Mirror from right side
                for x_idx in range(x_shift):
                    f_new[i, :, x_idx] = f_temp[:, x_shift - 1 - x_idx]
            elif x_shift < 0:  # Moving left (-x direction)
                x_abs = abs(x_shift)
                # Shift left: copy from right
                f_new[i, :, :-x_abs] = f_temp[:, x_abs:]
                # Right boundary: Mirror from left side
                for x_idx in range(x_abs):
                    f_new[i, :, self.nx - 1 - x_idx] = f_temp[:, self.nx - x_abs + x_idx]
            else:  # x_shift == 0
                f_new[i, :, :] = f_temp[:, :]
        
        self.f = f_new
    
    def _compute_macroscopic(self):
        """Compute macroscopic variables from distributions: rho = sum(f_i), u = sum(f_i * c_i) / rho"""
        self.rho = np.sum(self.f, axis=0)
        
        # Avoid division by zero
        rho_safe = np.maximum(self.rho, 1e-10)
        
        # Reshape velocity components for proper broadcasting
        self.ux = (self.f * self.c[:, 0].reshape(9, 1, 1)).sum(axis=0) / rho_safe
        self.uy = (self.f * self.c[:, 1].reshape(9, 1, 1)).sum(axis=0) / rho_safe
        
        self.velocity = np.sqrt(self.ux**2 + self.uy**2)
    
    def apply_boundary_conditions(self):
        """
        Apply boundary conditions at all borders.
        Called at the START of a timestep, before collision and streaming.
        """
        # Inlet (y=0): Dirichlet velocity (Zou-He scheme) - ceiling jet entrance
        self._zou_he_inlet()
        
    def _apply_outlet_bc(self):
        """
        Apply boundary conditions after streaming:
        1. Outlet (y=ny-1): Zero-gradient extrapolation
        """
        # Outlet (y=ny-1): Zero-gradient (extrapolation) - bottom of domain
        self.f[:, -1, :] = self.f[:, -2, :]
    
    def _zou_he_inlet(self):
        """
        Zou-He inlet boundary condition for constant velocity profile at ceiling (y=0).
        Reconstructs distribution function from known velocity and density.
        """
        y = 0
        
        # Set density at inlet (assume ambient density)
        rho_inlet = 1.0
        self.rho[y, :] = rho_inlet
        
        # Macroscopic velocity is already set by set_inlet_condition()
        # Just compute f_eq from current rho and u
        
        for i in range(9):
            cu = self.c[i, 0] * self.ux[y, :] + self.c[i, 1] * self.uy[y, :]
            u_sq = self.ux[y, :]**2 + self.uy[y, :]**2
            
            # Equilibrium distribution
            self.f_eq[i, y, :] = self.w[i] * rho_inlet * (
                1.0 + cu / self.c_s2 + 
                0.5 * cu**2 / (self.c_s2**2) - 
                0.5 * u_sq / self.c_s2
            )
            
            # For inlet BC, use pure equilibrium (could add non-eq part for stability)
            self.f[i, y, :] = self.f_eq[i, y, :]
    
    def step(self):
        """Execute one LBM time step: collision, streaming, boundary conditions"""
        # Inlet boundary condition (set ceiling jet profile)
        self.apply_boundary_conditions()
        
        # Compute macroscopic from current distribution
        self._compute_macroscopic()
        
        # Collision step
        self._compute_equilibrium()
        self.f = self.f - (self.f - self.f_eq) / self.tau
        
        # Streaming (advection) - allows flow to leave domain
        self.streaming_step()
        
        # Outlet BC: prevent backflow
        self._apply_outlet_bc()
    
    def run_simulation(self, n_steps, jet_center_x=None, jet_radius=None, 
                     jet_velocity=None, callback=None):
        """
        Run simulation for n_steps.
        
        Args:
            n_steps: Number of time steps
            jet_center_x: Center of inlet jet in x-direction at y=0 ceiling (default: center)
            jet_radius: Radius of inlet jet (default: 2)
            jet_velocity: Inlet velocity (default: self.u0)
            callback: Function called after each step with progress
        """
        if jet_center_x is None:
            jet_center_x = self.nx / 2.0
        if jet_radius is None:
            jet_radius = 2.0
        if jet_velocity is None:
            jet_velocity = self.u0
        
        self.set_inlet_condition(jet_center_x, jet_radius, jet_velocity)
        
        for step in range(n_steps):
            self.step()
            
            if callback and (step % max(1, n_steps // 10) == 0):
                callback(step, n_steps)
        
        self._compute_macroscopic()  # Final computation
        return self.velocity, self.ux, self.uy, self.rho


def create_default_simulator(reynolds=100, mach=0.1):
    """Create a default LBM simulator with reasonable parameters"""
    return LBM2DJetSimulator(nx=80, ny=50, reynolds_number=reynolds, mach_number=mach)
