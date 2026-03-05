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
    
    def __init__(self, parent=None):
        self.fig = Figure(figsize=(10, 6), dpi=90)
        self.ax = self.fig.add_subplot(111)
        super().__init__(self.fig)
        self.setParent(parent)
        
        self.parent_window = parent
        self.velocity_data = None
        self.ux_data = None
        self.uy_data = None
        self.rho_data = None
        self.lattice_spacing = 1.0
        
    def set_simulation_data(self, velocity, ux, uy, rho, lattice_spacing=1.0):
        """
        Set simulation results for visualization.
        
        Args:
            velocity: Velocity magnitude field [ny, nx]
            ux: X-velocity component [ny, nx]
            uy: Y-velocity component [ny, nx]
            rho: Density field [ny, nx]
            lattice_spacing: Physical spacing between lattice nodes (meters)
        """
        self.velocity_data = velocity
        self.ux_data = ux
        self.uy_data = uy
        self.rho_data = rho
        self.lattice_spacing = lattice_spacing
        self.update_visualization()
    
    def update_visualization(self):
        """Redraw visualization with current data"""
        self.ax.clear()
        
        if self.velocity_data is None:
            self.ax.text(0.5, 0.5, 'No simulation data', 
                        ha='center', va='center', transform=self.ax.transAxes)
            self.fig.canvas.draw_idle()
            return
        
        ny, nx = self.velocity_data.shape
        xx = np.arange(nx) * self.lattice_spacing
        yy = np.arange(ny) * self.lattice_spacing - (ny / 2.0) * self.lattice_spacing
        XX, YY = np.meshgrid(xx, yy)
        
        # Plot velocity magnitude as contours + colormap
        levels = np.linspace(0, np.max(self.velocity_data) * 0.9, 20)
        contourf = self.ax.contourf(XX, YY, self.velocity_data, levels=levels, cmap='hot')
        contour = self.ax.contour(XX, YY, self.velocity_data, levels=levels[::2], 
                                 colors='black', alpha=0.3, linewidths=0.5)
        
        cbar = self.fig.colorbar(contourf, ax=self.ax, label='Velocity Magnitude')
        
        # Streamlines
        # Sample velocity field for better visualization
        stride = 2
        ux_sample = self.ux_data[::stride, ::stride]
        uy_sample = self.uy_data[::stride, ::stride]
        XX_sample = XX[::stride, ::stride]
        YY_sample = YY[::stride, ::stride]
        
        self.ax.streamplot(XX_sample, YY_sample, ux_sample, uy_sample, 
                          color='white', linewidth=0.8, density=1.5, arrowsize=1.5)
        
        # Formatting
        self.ax.set_xlabel('Axial Distance (m)', fontsize=11, weight='bold')
        self.ax.set_ylabel('Radial Distance (m)', fontsize=11, weight='bold')
        self.ax.set_title('2D Jet Velocity Field (LBM Simulation)', 
                         fontsize=12, weight='bold')
        self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3, linestyle=':')
        
        self.fig.tight_layout()
        self.fig.canvas.draw_idle()


class CFDResultsDialog(QtWidgets.QDialog):
    """Dialog to display and interact with CFD simulation results"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('CFD Simulation Results - 2D Jet')
        self.resize(1000, 700)
        
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
    
    def set_results(self, velocity, ux, uy, rho, reynolds, mach, 
                   jet_diameter_physical, lattice_spacing=1.0):
        """Set and display simulation results"""
        self.visualizer.set_simulation_data(velocity, ux, uy, rho, lattice_spacing)
        
        # Update statistics
        max_vel = np.max(velocity)
        self.lbl_max_velocity.setText(f'Max Velocity: {max_vel:.4f} m/s')
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
