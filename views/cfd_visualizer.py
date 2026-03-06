# -*- coding: utf-8 -*-
"""
CFD Visualization widget for Matplotlib.
Displays velocity field, streamlines, and contours from LBM simulation.
"""

import numpy as np
from PyQt5 import QtWidgets, QtCore
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class CFDVisualizerWidget(FigureCanvas):
    """Widget to display 2D jet CFD results"""

    @staticmethod
    def compute_focus_xlim(velocity_m_s, x_coords, color_max, threshold_ratio=0.03):
        """
        Estimate symmetric x-limits around the active jet region.

        This keeps the plot readable while preserving physical aspect ratio.
        """
        if velocity_m_s is None or velocity_m_s.size == 0:
            return None

        ny, nx = velocity_m_s.shape
        if nx < 3:
            return None

        threshold = max(float(color_max) * float(threshold_ratio), 1e-6)
        col_max = np.nanmax(velocity_m_s, axis=0)
        active_cols = np.where(col_max >= threshold)[0]
        if active_cols.size < 2:
            return None

        pad = max(2, int(0.08 * nx))
        left = max(0, int(active_cols[0]) - pad)
        right = min(nx - 1, int(active_cols[-1]) + pad)

        half_span = max(abs(float(x_coords[left])), abs(float(x_coords[right])))
        domain_half = max(abs(float(x_coords[0])), abs(float(x_coords[-1])))
        half_span = min(max(half_span, 1e-6), domain_half)

        return (-half_span, half_span)

    @staticmethod
    def compute_visible_velocity_max(velocity_m_s, percentile=99.9, reference_velocity=None):
        """
        Compute a robust maximum for visualization, excluding edge cells and spikes.

        Args:
            velocity_m_s: Velocity magnitude field in m/s
            percentile: High percentile used as visible max
            reference_velocity: Optional reference speed (e.g. V0) in m/s

        Returns:
            Robust visible maximum velocity in m/s
        """
        arr = np.asarray(velocity_m_s, dtype=np.float64)
        if arr.size == 0:
            return 1e-6

        # Drop one-cell border where BC artifacts are most common.
        if arr.shape[0] > 2 and arr.shape[1] > 2:
            arr_core = arr[1:-1, 1:-1]
        else:
            arr_core = arr

        finite_vals = arr_core[np.isfinite(arr_core)]
        if finite_vals.size == 0:
            return 1e-6

        finite_vals = np.clip(finite_vals, 0.0, None)

        # Focus percentile on the active jet region (exclude near-zero ambient values).
        if reference_velocity is not None and reference_velocity > 0.0:
            active_threshold = 0.05 * float(reference_velocity)
            active_vals = finite_vals[finite_vals >= active_threshold]
            if active_vals.size >= 10:
                vals_for_percentile = active_vals
            else:
                vals_for_percentile = finite_vals[finite_vals > 0.0]
        else:
            vals_for_percentile = finite_vals[finite_vals > 0.0]

        if vals_for_percentile.size == 0:
            vals_for_percentile = finite_vals

        robust_max = float(np.percentile(vals_for_percentile, percentile))
        if not np.isfinite(robust_max) or robust_max <= 0.0:
            robust_max = float(np.max(finite_vals))

        # Keep physical readability when a reference velocity exists.
        if reference_velocity is not None and reference_velocity > 0.0:
            robust_max = max(robust_max, float(reference_velocity))

        return max(robust_max, 1e-6)
    
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(13, 8), dpi=100)
        self.ax_vel = None
        self.ax_temp = None
        super().__init__(self.fig)
        self.setParent(parent)
        
        self.parent_window = parent
        self.velocity_data = None
        self.ux_data = None
        self.uy_data = None
        self.rho_data = None
        self.lattice_spacing = 1.0
        self.color_max_m_s = None
        self.x0_distance_m = None
        self.x_zero_index = None
        self.x_domain_limits = None
        self.ambient_temp_c = None
        self.nozzle_temp_c = None
        self.inlet_velocity_m_s = None
        self.temperature_data = None
        self.jet_diameter_m = None
        self.include_buoyancy = False
        
    def set_simulation_data(
        self,
        velocity,
        ux,
        uy,
        rho,
        lattice_spacing=1.0,
        velocity_scale=1.0,
        color_max_m_s=None,
        x0_distance_m=None,
        x_zero_index=None,
        x_domain_limits=None,
        ambient_temp_c=None,
        nozzle_temp_c=None,
        inlet_velocity_m_s=None,
        jet_diameter_m=None,
        include_buoyancy=False,
    ):
        """
        Set simulation results for visualization.
        
        Data is used directly without transposition:
        - x (horizontal) = radial direction
        - y (vertical) = axial direction (downward when inverted)
        
        Args:
            velocity: Velocity magnitude field [ny, nx] (lattice units)
            ux: X-velocity component (radial) [ny, nx] (lattice units)
            uy: Y-velocity component (axial) [ny, nx] (lattice units)
            rho: Density field [ny, nx]
            lattice_spacing: Physical spacing between lattice nodes (meters)
            velocity_scale: Conversion factor from lattice to m/s
            color_max_m_s: Optional upper bound for color scale in m/s
        """
        # No transposition - use data as-is
        # Shape is (ny, nx) which matches matplotlib's (y, x) convention
        self.velocity_data = velocity * velocity_scale  # Convert to m/s
        self.ux_data = ux * velocity_scale  # Radial velocity
        self.uy_data = uy * velocity_scale  # Axial velocity (downward)
        self.rho_data = rho
        self.lattice_spacing = lattice_spacing
        self.color_max_m_s = color_max_m_s
        self.x0_distance_m = x0_distance_m
        self.x_zero_index = x_zero_index
        self.x_domain_limits = x_domain_limits
        self.ambient_temp_c = ambient_temp_c
        self.nozzle_temp_c = nozzle_temp_c
        self.inlet_velocity_m_s = inlet_velocity_m_s
        self.jet_diameter_m = jet_diameter_m
        self.include_buoyancy = bool(include_buoyancy)
        self.temperature_data = self._compute_temperature_field()
        self.update_visualization()

    def _compute_temperature_field(self):
        """Build a temperature proxy field for isotherm plotting from velocity mixing."""
        if self.velocity_data is None:
            return None
        if self.ambient_temp_c is None or self.nozzle_temp_c is None:
            return None

        ambient = float(self.ambient_temp_c)
        nozzle = float(self.nozzle_temp_c)

        if self.inlet_velocity_m_s is not None and self.inlet_velocity_m_s > 0.0:
            norm_v = self.velocity_data / float(self.inlet_velocity_m_s)
        else:
            vmax = max(float(np.max(self.velocity_data)), 1e-6)
            norm_v = self.velocity_data / vmax

        norm_v = np.clip(norm_v, 0.0, 1.0)

        # Linear mixing proxy: faster core is closer to nozzle temperature.
        return ambient - (ambient - nozzle) * norm_v
    
    def update_visualization(self):
        """
        Redraw visualization with current data.
        
        Domain layout:
        - X (horizontal): radial direction, centered at 0
        - Y (vertical): physical height, with 0 at floor and max at ceiling
        - Data format: [ny, nx] - rows are y positions, columns are x positions
        """
        self.fig.clear()
        self.ax_vel = self.fig.add_subplot(211)
        self.ax_temp = self.fig.add_subplot(212, sharex=self.ax_vel, sharey=self.ax_vel)

        if self.velocity_data is None:
            self.ax_vel.text(0.5, 0.5, 'No simulation data',
                             ha='center', va='center', transform=self.ax_vel.transAxes)
            self.fig.canvas.draw_idle()
            return
        
        ny, nx = self.velocity_data.shape
        
        # Create meshgrid: x = radial with configurable origin index for x=0.
        x0_idx = float(self.x_zero_index) if self.x_zero_index is not None else (nx / 2.0)
        xx = (np.arange(nx) - x0_idx) * self.lattice_spacing
        # Physical height coordinate: strictly increasing from floor (0) to ceiling (H).
        yy = np.arange(ny) * self.lattice_spacing
        XX, YY = np.meshgrid(xx, yy)
        
        # White background - no fill
        self.ax_vel.set_facecolor('white')
        self.ax_temp.set_facecolor('white')
        self.fig.patch.set_facecolor('white')
        
        # Colored contour LINES only (no fill) with robust max to suppress spikes.
        data_visible_max = self.compute_visible_velocity_max(
            self.velocity_data,
            percentile=99.9,
            reference_velocity=self.color_max_m_s,
        )
        if self.color_max_m_s is not None:
            # Never hide the robust maximum even if an external cap is passed.
            color_max = max(float(self.color_max_m_s), data_visible_max)
        else:
            color_max = data_visible_max
        color_max = max(color_max, 1e-6)
        # Convert storage convention (row 0 at ceiling) to plotting convention (row 0 at floor).
        velocity_for_plot = np.flipud(np.clip(self.velocity_data, 0.0, color_max))

        # Fixed velocity isoline spacing: every 0.2 m/s.
        velocity_step = 0.2
        if color_max >= velocity_step:
            levels = np.arange(velocity_step, color_max + 0.5 * velocity_step, velocity_step)
        else:
            # Fallback for low-speed cases where 0.2 m/s spacing cannot form contours.
            level_min = max(color_max / 200.0, 1e-6)
            levels = np.linspace(level_min, color_max, 6)
        contour = self.ax_vel.contour(XX, YY, velocity_for_plot, levels=levels,
                          cmap='plasma', linewidths=2.0, alpha=1.0)

        # Label only a subset of levels and skip the top strip to avoid clutter at ceiling.
        label_levels = levels[1::2]
        label_field = velocity_for_plot.copy()
        top_rows_to_skip = max(1, int(0.08 * ny))
        label_field[-top_rows_to_skip:, :] = np.nan
        label_contour = self.ax_vel.contour(
            XX,
            YY,
            label_field,
            levels=label_levels,
            colors='black',
            linewidths=0.0,
            alpha=0.0,
        )
        self.ax_vel.clabel(
            label_contour,
            label_levels,
            inline=True,
            inline_spacing=2,
            fontsize=10,
            fmt='%.2f',
            colors='black',
        )
        
        # Add colorbar to show velocity scale
        self.fig.colorbar(contour, ax=self.ax_vel, label='Velocity Magnitude (m/s)')
        
        # Streamlines: white, subtle, just to show flow direction
        stride = max(1, ny // 20)
        ux_for_plot = np.flipud(self.ux_data)
        # y-axis is upward-positive in plot coordinates, so downward CFD velocity is negative.
        uy_for_plot = -np.flipud(self.uy_data)
        ux_sample = ux_for_plot[::stride, ::stride]
        uy_sample = uy_for_plot[::stride, ::stride]
        XX_sample = XX[::stride, ::stride]
        YY_sample = YY[::stride, ::stride]
        
        self.ax_vel.streamplot(XX_sample, YY_sample, ux_sample, uy_sample,
                       color='lightgray', linewidth=0.7, density=1.5,
                       arrowsize=1.5, arrowstyle='->')
        
        # Formatting
        self.ax_vel.set_ylabel('Axial Distance (m)', fontsize=11, weight='bold')
        self.ax_vel.set_title('2D Jet Velocity Field from Ceiling Diffuser (LBM Simulation)',
                      fontsize=12, weight='bold')
        self.ax_vel.set_aspect('equal', adjustable='box')

        # Respect explicit user-defined domain limits when available.
        if self.x_domain_limits is not None:
            x_left, x_right = self.x_domain_limits
            if x_left < x_right:
                self.ax_vel.set_xlim(float(x_left), float(x_right))
        else:
            # Auto-focus x-range around active jet so equal aspect remains readable.
            focus_xlim = self.compute_focus_xlim(velocity_for_plot, xx, color_max)
            if focus_xlim is not None:
                self.ax_vel.set_xlim(*focus_xlim)

        self.ax_vel.grid(True, alpha=0.2, linestyle=':', color='gray')
        
        # Temperature panel (isotherms)
        if self.temperature_data is not None:
            temperature_for_plot = np.flipud(self.temperature_data)
            t_min = float(np.min(temperature_for_plot))
            t_max = float(np.max(temperature_for_plot))
            temp_step = 2.0
            if t_max - t_min < 1e-6:
                t_levels = np.linspace(t_min - 0.1, t_max + 0.1, 8)
            else:
                t_start = np.ceil(t_min / temp_step) * temp_step
                t_levels = np.arange(t_start, t_max + 0.5 * temp_step, temp_step)
                if t_levels.size < 2:
                    # Fallback to preserve visible contours in narrow temperature bands.
                    t_levels = np.linspace(t_min, t_max, 6)

            contour_t = self.ax_temp.contour(
                XX,
                YY,
                temperature_for_plot,
                levels=t_levels,
                cmap='coolwarm',
                linewidths=1.8,
                alpha=0.95,
            )
            self.ax_temp.clabel(
                contour_t,
                contour_t.levels[::2],
                inline=True,
                inline_spacing=2,
                fontsize=9,
                fmt='%.1f',
                colors='black',
            )
            self.fig.colorbar(contour_t, ax=self.ax_temp, label='Temperature (degC)')
        else:
            self.ax_temp.text(
                0.5,
                0.5,
                'Isotherms unavailable (missing TA/T0)',
                ha='center',
                va='center',
                transform=self.ax_temp.transAxes,
                fontsize=10,
                color='gray',
            )

        self.ax_temp.set_xlabel('Radial Distance (m)', fontsize=11, weight='bold')
        self.ax_temp.set_ylabel('Axial Distance (m)', fontsize=11, weight='bold')
        self.ax_temp.set_title('Isotherm Field (Temperature Proxy)', fontsize=11, weight='bold')
        self.ax_temp.set_aspect('equal', adjustable='box')
        self.ax_temp.grid(True, alpha=0.2, linestyle=':', color='gray')

        if self.x_domain_limits is not None:
            x_left, x_right = self.x_domain_limits
            if x_left < x_right:
                self.ax_temp.set_xlim(float(x_left), float(x_right))

        self.ax_vel.tick_params(labelbottom=False)

        # Mark X0 distance from ceiling if provided and inside visible y-range.
        if self.x0_distance_m is not None:
            y_min, y_max = self.ax_vel.get_ylim()
            y_low = min(y_min, y_max)
            y_high = max(y_min, y_max)
            domain_height = (ny - 1) * self.lattice_spacing
            x0_from_ceiling = float(self.x0_distance_m)
            x0 = domain_height - x0_from_ceiling
            if y_low <= x0 <= y_high:
                self.ax_vel.axhline(
                    y=x0,
                    color='dimgray',
                    linestyle='--',
                    linewidth=1.5,
                    alpha=0.9,
                    zorder=6,
                )
                x_left, x_right = self.ax_vel.get_xlim()
                self.ax_vel.text(
                    x_left + 0.01 * (x_right - x_left),
                    x0,
                    f'X0 = {x0_from_ceiling:.2f} m',
                    fontsize=10,
                    color='dimgray',
                    ha='left',
                    va='bottom',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.65, pad=1.5),
                    zorder=7,
                )

                self.ax_temp.axhline(
                    y=x0,
                    color='dimgray',
                    linestyle='--',
                    linewidth=1.2,
                    alpha=0.8,
                    zorder=6,
                )
                x_left_t, x_right_t = self.ax_temp.get_xlim()
                self.ax_temp.text(
                    x_left_t + 0.01 * (x_right_t - x_left_t),
                    x0,
                    f'X0 = {x0_from_ceiling:.2f} m',
                    fontsize=10,
                    color='dimgray',
                    ha='left',
                    va='bottom',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.65, pad=1.5),
                    zorder=7,
                )
        
        self.fig.tight_layout()
        self.fig.canvas.draw_idle()


class CFDResultsDialog(QtWidgets.QDialog):
    """Dialog to display and interact with CFD simulation results"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('CFD Simulation Results - 2D Jet')
        self.resize(1350, 900)
        
        layout = QtWidgets.QVBoxLayout(self)
        
        # Visualizer widget
        self.visualizer = CFDVisualizerWidget(self)
        layout.addWidget(self.visualizer)
        
        # Statistics panel
        stats_layout = QtWidgets.QHBoxLayout()
        
        self.lbl_max_velocity = QtWidgets.QLabel('Max Velocity: —')
        self.lbl_reynolds = QtWidgets.QLabel('Reynolds: —')
        self.lbl_jet_diameter = QtWidgets.QLabel('Jet Diameter: —')
        self.lbl_mach = QtWidgets.QLabel('Mach: —')
        
        stats_layout.addWidget(self.lbl_max_velocity)
        stats_layout.addWidget(self.lbl_reynolds)
        stats_layout.addWidget(self.lbl_jet_diameter)
        stats_layout.addWidget(self.lbl_mach)
        
        layout.addLayout(stats_layout)
        
        # Export button
        self.btn_export = QtWidgets.QPushButton('Export Results as PNG')
        self.btn_export.clicked.connect(self.export_image)
        layout.addWidget(self.btn_export)
        
        # Close button
        self.btn_close = QtWidgets.QPushButton('Close')
        self.btn_close.clicked.connect(self.close)
        layout.addWidget(self.btn_close)
    
    def set_results(
        self,
        velocity,
        ux,
        uy,
        rho,
        reynolds,
        mach,
        jet_diameter_physical,
        inlet_velocity=None,
        axial_distance_x0=None,
        x_zero_index=None,
        x_domain_limits=None,
        ambient_temp_c=None,
        nozzle_temp_c=None,
        include_buoyancy=False,
        lattice_spacing=1.0,
        velocity_scale=1.0,
    ):
        """
        Set and display simulation results.
        
        Args:
            velocity_scale: Conversion factor from lattice units to m/s
        """
        velocity_m_s = velocity * velocity_scale
        raw_max = float(np.max(velocity_m_s))
        visible_max = self.visualizer.compute_visible_velocity_max(
            velocity_m_s,
            percentile=99.9,
            reference_velocity=inlet_velocity,
        )

        self.visualizer.set_simulation_data(
            velocity,
            ux,
            uy,
            rho,
            lattice_spacing,
            velocity_scale,
            color_max_m_s=visible_max,
            x0_distance_m=axial_distance_x0,
            x_zero_index=x_zero_index,
            x_domain_limits=x_domain_limits,
            ambient_temp_c=ambient_temp_c,
            nozzle_temp_c=nozzle_temp_c,
            inlet_velocity_m_s=inlet_velocity,
            jet_diameter_m=jet_diameter_physical,
            include_buoyancy=include_buoyancy,
        )
        
        # Update statistics
        max_vel = raw_max
        if inlet_velocity is not None:
            self.lbl_max_velocity.setText(
                f'Max Velocity (visible, P99.9): {visible_max:.4f} m/s | '
                f'Raw: {max_vel:.4f} m/s | V0: {inlet_velocity:.4f} m/s'
            )
        else:
            self.lbl_max_velocity.setText(
                f'Max Velocity (visible, P99.9): {visible_max:.4f} m/s | Raw: {max_vel:.4f} m/s'
            )
        self.lbl_reynolds.setText(f'Reynolds: {reynolds:.1f}')
        self.lbl_mach.setText(f'Mach: {mach:.4f}')
        self.lbl_jet_diameter.setText(f'Jet Diameter: {jet_diameter_physical:.4f} m')
    
    def export_image(self):
        """Export visualization as PNG"""
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Export Results', 'jet_simulation.png', 'PNG Images (*.png)'
        )
        if path:
            self.visualizer.fig.savefig(path, dpi=150, bbox_inches='tight')
            QtWidgets.QMessageBox.information(self, 'Success', 
                                             f'Results exported to:\n{path}')
