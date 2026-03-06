# -*- coding: utf-8 -*-
"""
Main controller for Spot Cooling Designer application.
Handles business logic and user interactions.
"""

from models.cooling_calc import solve_spot_cooling, solve_spot_cooling_from_target_point


class MainController:
    """Main application controller"""
    
    def __init__(self, main_window):
        """
        Initialize controller with reference to main window.
        
        Args:
            main_window: Reference to MainWindow instance for accessing UI controls
        """
        self.main_window = main_window
    
    def solve_normal(self):
        """
        Solve spot cooling problem with normal parameters.
        Uses target velocity to find required nozzle conditions.
        
        Returns: dict with solution from solve_spot_cooling
        """
        target_velocity = self.main_window.e_VJ.value()
        result = solve_spot_cooling(
            ambient_temp=self.main_window.e_TA.value(),
            ambient_rh_fraction=self.main_window.e_RHA.value() / 100.0,
            mean_radiant_temp=self.main_window.e_TMR.value(),
            target_velocity=target_velocity,
            metabolic_rate=self.main_window.e_M.value(),
            clothing_icl=self.main_window.e_ICL.value(),
            jet_diameter=self.main_window.e_D0.value(),
            axial_distance=self.main_window.e_X0.value(),
            nozzle_rh_fraction=self.main_window.e_RH0.value() / 100.0,
            p_atm_kpa=self.main_window.e_PATM.value(),
            include_buoyancy=self.main_window.chk_buoy.isChecked()
        )
        result['Vj_used'] = target_velocity
        return result
    
    def solve_from_selected_point(self, target_temp, target_humidity_ratio, selected_velocity=None):
    
        """
        Inverse solve: find nozzle conditions to achieve target point.
        Uses selected point on psychrometric chart as constraint.

        Args:
            target_temp: Target temperature at measurement area [°C]
            target_humidity_ratio: Target humidity ratio at measurement area [kg/kg]
            selected_velocity: Optional Vj from selected chart point [m/s]

        Returns: dict with solution from solve_spot_cooling_from_target_point
        """
        target_velocity = (
            self.main_window.e_VJ.value()
            if selected_velocity is None else float(selected_velocity)
        )

        result = solve_spot_cooling_from_target_point(
            ambient_temp=self.main_window.e_TA.value(),
            ambient_rh_fraction=self.main_window.e_RHA.value() / 100.0,
            mean_radiant_temp=self.main_window.e_TMR.value(),
            target_velocity=target_velocity,
            metabolic_rate=self.main_window.e_M.value(),
            clothing_icl=self.main_window.e_ICL.value(),
            jet_diameter=self.main_window.e_D0.value(),
            axial_distance=self.main_window.e_X0.value(),
            target_temp=target_temp,
            target_humidity_ratio=target_humidity_ratio,
            p_atm_kpa=self.main_window.e_PATM.value(),
            include_buoyancy=self.main_window.chk_buoy.isChecked()
        )
        result['Vj_used'] = target_velocity
        return result
