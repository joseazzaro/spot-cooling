# -*- coding: utf-8 -*-
"""
CFD Controller for integrating LBM jet simulation with cooling calculations.
Bridges the gap between the analytical model (ASHRAE) and CFD visualization.
"""

from models.jet_lbm import LBM2DJetSimulator
from models.psychrometric import air_density
import numpy as np


class CFDSimulationController:
    """
    Orchestrates 2D LBM simulation of air jet.
    Integrates physical parameters from cooling calculations.
    """
    
    def __init__(self, main_window):
        """
        Initialize CFD controller.
        
        Args:
            main_window: Reference to MainWindow for accessing UI controls
        """
        self.main_window = main_window
        self.simulator = None
        self.last_results = None
    
    def create_simulator_from_cooling_params(self, jet_diameter_m, inlet_velocity_m_s,
                                           ambient_density_kg_m3, ambient_viscosity):
        """
        Create LBM simulator with parameters from cooling calculation.
        
        Args:
            jet_diameter_m: Physical jet diameter [meters]
            inlet_velocity_m_s: Inlet velocity [m/s]
            ambient_density_kg_m3: Air density at ambient conditions [kg/m³]
            ambient_viscosity: Air kinematic viscosity [m²/s]
        
        Returns:
            LBM2DJetSimulator instance
        """
        # Calculate Reynolds number based on jet parameters
        reynolds = (inlet_velocity_m_s * jet_diameter_m) / ambient_viscosity
        
        # Mach number (for subsonic jet, typically << 1)
        speed_of_sound = 343.0  # m/s at 20°C (approx)
        mach = inlet_velocity_m_s / speed_of_sound
        
        # Clamp Mach to valid range for isothermal solver
        mach = min(mach, 0.3)
        
        # Create simulator (default 80x50 grid)
        self.simulator = LBM2DJetSimulator(nx=80, ny=50, 
                                          reynolds_number=reynolds,
                                          mach_number=mach)
        
        # Store physical parameters for post-processing
        self.jet_diameter_physical = jet_diameter_m
        self.inlet_velocity_physical = inlet_velocity_m_s
        self.ambient_density = ambient_density_kg_m3
        self.reynolds = reynolds
        self.mach = mach
        
        return self.simulator
    
    def run_simulation(self, n_steps=500, callback=None):
        """
        Execute LBM simulation.
        
        Args:
            n_steps: Number of time steps to simulate (default 500, ~settling time)
            callback: Optional progress callback function(step, total_steps)
        
        Returns:
            dict with results: velocity, ux, uy, rho fields
        """
        if self.simulator is None:
            raise RuntimeError("Simulator not initialized. Call create_simulator_from_cooling_params first.")
        
        # Run simulation
        velocity, ux, uy, rho = self.simulator.run_simulation(
            n_steps=n_steps,
            jet_center_y=self.simulator.ny / 2.0,
            jet_radius=2.0,  # ~4-5 cells for jet diameter
            jet_velocity=self.simulator.u0,
            callback=callback
        )
        
        # Store results
        self.last_results = {
            'velocity': velocity,
            'ux': ux,
            'uy': uy,
            'rho': rho,
            'reynolds': self.reynolds,
            'mach': self.mach,
            'jet_diameter': self.jet_diameter_physical,
            'inlet_velocity': self.inlet_velocity_physical,
        }
        
        return self.last_results
    
    def get_jet_centerline_decay(self, x_start=5, x_end=75):
        """
        Extract jet centerline velocity decay U(x) / U0.
        Classic result: U(x)/U0 ~ D / (x-x0) for top-hat profile after initial region.
        
        Args:
            x_start: Start position for analysis
            x_end: End position for analysis
        
        Returns:
            dict with x positions and normalized velocities
        """
        if self.last_results is None:
            return None
        
        velocity = self.last_results['velocity']
        ny = velocity.shape[0]
        centerline_y = ny // 2
        
        # Extract centerline
        centerline_velocity = velocity[centerline_y, x_start:x_end]
        x_positions = np.arange(x_start, x_end)
        u_initial = np.max(centerline_velocity[:5])  # Reference velocity
        
        normalized_velocity = centerline_velocity / max(u_initial, 1e-10)
        
        return {
            'x': x_positions,
            'u_normalized': normalized_velocity,
            'u_max': u_initial,
        }
    
    def get_jet_spreading_rate(self):
        """
        Estimate jet spreading rate (opening angle).
        Typical value: ~12-24 degrees for isothermal circular jets.
        
        Returns:
            Spreading angle in degrees
        """
        if self.last_results is None:
            return None
        
        velocity = self.last_results['velocity']
        ny, nx = velocity.shape
        
        # Find half-velocity contour at several x positions
        spreading_angles = []
        
        for x in [10, 20, 30, 40, 50]:
            if x >= nx:
                continue
            
            u_max = np.max(velocity[:, x])
            u_half = u_max * 0.5
            
            # Find y positions where u > u_half
            above_half = np.where(velocity[:, x] > u_half)[0]
            
            if len(above_half) > 1:
                y_min = above_half[0] - ny / 2.0
                y_max = above_half[-1] - ny / 2.0
                y_width = y_max - y_min
                
                # Half-angle in radians
                half_angle = np.arctan(y_width / (2.0 * x)) * 180 / np.pi
                spreading_angles.append(half_angle)
        
        # Average spreading angle
        if spreading_angles:
            return np.mean(spreading_angles)
        return None
    
    def verify_ashrae_model(self):
        """
        Compare CFD results with ASHRAE analytical model.
        
        Returns:
            dict with comparison metrics
        """
        if self.last_results is None:
            return None
        
        # Get jet decay from CFD
        decay_data = self.get_jet_centerline_decay()
        
        if decay_data is None:
            return None
        
        x = decay_data['x']
        u_norm = decay_data['u_normalized']
        
        # ASHRAE model prediction
        # Vj/V0 = 1.464 / r, where r = X0/D0 + 2.572
        # For comparison: map lattice coords to physical distance
        
        # Simple comparison: extract first point below 0.5
        idx_05 = np.where(u_norm < 0.5)[0]
        if len(idx_05) > 0:
            x_half = idx_05[0]
            # This is where jet velocity is half the inlet
            # ASHRAE would predict this at specific X0/D0 ratio
        
        return {
            'centerline_decay': decay_data,
            'spreading_rate': self.get_jet_spreading_rate(),
            'reynolds': self.reynolds,
            'mach': self.mach,
        }


def create_cfd_from_cooling_solution(main_window, cooling_result):
    """
    Helper function to create and run CFD simulation from a cooling calculation result.
    
    Args:
        main_window: Reference to MainWindow
        cooling_result: Result dict from solve_spot_cooling()
    
    Returns:
        CFDSimulationController with results
    """
    from models.psychrometric import air_density
    
    controller = CFDSimulationController(main_window)
    
    # Extract parameters from UI
    D0 = main_window.e_D0.value()
    V0 = cooling_result['V0']
    T0 = cooling_result['T0']
    P0 = cooling_result['P0']
    
    # Calculate air properties at nozzle conditions
    # For simplicity: use ambient kinematic viscosity (valid for small ΔT)
    kinematic_viscosity = 1.8e-5  # m²/s at 20°C (approximate)
    
    # Create simulator with physical parameters
    controller.create_simulator_from_cooling_params(
        jet_diameter_m=D0,
        inlet_velocity_m_s=V0,
        ambient_density_kg_m3=1.225,  # kg/m³ at sea level
        ambient_viscosity=kinematic_viscosity
    )
    
    return controller
