# -*- coding: utf-8 -*-
"""
Main application window for Spot Cooling Designer.
PyQt5-based GUI for ASHRAE spot cooling calculations.
"""

import math
import numpy as np
from PyQt5 import QtWidgets, QtCore, QtGui, QtPrintSupport

from views.psychro_chart import PsychroChart


class MainWindow(QtWidgets.QMainWindow):
    """Main application window"""
    
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Spot Cooling - ASHRAE')
        self.resize(1600, 900)
        self.last_report_html = None
        
        # Central widget layout
        cw = QtWidgets.QWidget()
        self.setCentralWidget(cw)
        layout = QtWidgets.QHBoxLayout(cw)
        
        # Left panel: Form controls
        form = QtWidgets.QFormLayout()
        form_widget = QtWidgets.QWidget()
        form_widget.setLayout(form)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(form_widget)
        left_panel = QtWidgets.QWidget()
        left_panel.setLayout(QtWidgets.QVBoxLayout())
        left_panel.layout().addWidget(scroll)
        
        # Right panel: Psychrometric chart
        self.chart = PsychroChart(parent=self)
        
        # Splitter for resizable panels
        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        splitter.addWidget(left_panel)
        splitter.addWidget(self.chart)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 2)
        splitter.setSizes([350, 1000])
        layout.addWidget(splitter)
        
        # Input controls
        self._create_controls()
        
        # Add controls to form
        self._populate_form(form)
        
        # Connect signals
        self._connect_signals()
        
        # Chart update timer (debounce)
        # Prevents rapid redraws when user adjusts spinbox values
        # Only chart-affecting parameters trigger redraws (220ms after last change)
        self._chart_update_timer = QtCore.QTimer(self)
        self._chart_update_timer.setSingleShot(True)
        self._chart_update_timer.setInterval(220)
        self._chart_update_timer.timeout.connect(self._apply_chart_update)
        
        # Track previous parameters to detect changes that invalidate solutions
        # If Icl, M, Tmr, or Vj change, the operating solution becomes invalid
        # and must be cleared to prevent displaying incorrect results
        self._prev_icl = self.e_ICL.value()
        self._prev_m = self.e_M.value()
        self._prev_tmr = self.e_TMR.value()
        self._prev_vj = self.e_VJ.value()
        
        # Disable keyboard tracking for spinboxes
        spinboxes = [
            self.e_TA, self.e_RHA, self.e_TMR, self.e_VJ, self.e_M, self.e_ICL,
            self.e_D0, self.e_X0, self.e_RT, self.e_ang, self.e_RH0,
            self.e_PATM
        ]
        for spin in spinboxes:
            spin.setKeyboardTracking(False)
    
    def _create_controls(self):
        """Create all input controls"""
        # Ambient conditions
        self.e_TA = QtWidgets.QDoubleSpinBox()
        self.e_TA.setRange(-20, 80)
        self.e_TA.setValue(40.0)
        self.e_TA.setSuffix(' C')
        
        self.e_RHA = QtWidgets.QDoubleSpinBox()
        self.e_RHA.setRange(0, 100)
        self.e_RHA.setValue(50.0)
        self.e_RHA.setSuffix(' %')
        
        self.e_TMR = QtWidgets.QDoubleSpinBox()
        self.e_TMR.setRange(-20, 120)
        self.e_TMR.setValue(45.0)
        self.e_TMR.setSuffix(' C')
        
        self.e_VJ = QtWidgets.QDoubleSpinBox()
        self.e_VJ.setRange(0.2, 3.0)
        self.e_VJ.setSingleStep(0.1)
        self.e_VJ.setValue(2.0)
        self.e_VJ.setSuffix(' m/s')
        
        # Work & clothing
        self.e_M = QtWidgets.QDoubleSpinBox()
        self.e_M.setRange(40, 400)
        self.e_M.setValue(87.0)
        self.e_M.setSuffix(' W/m2')
        
        self.e_ICL = QtWidgets.QDoubleSpinBox()
        self.e_ICL.setRange(0.0, 2.0)
        self.e_ICL.setSingleStep(0.05)
        self.e_ICL.setValue(0.6)
        self.e_ICL.setSuffix(' clo')
        
        # Jet geometry
        self.e_D0 = QtWidgets.QDoubleSpinBox()
        self.e_D0.setRange(0.03, 0.5)
        self.e_D0.setSingleStep(0.005)
        self.e_D0.setValue(0.3048)
        self.e_D0.setSuffix(' m')
        
        self.e_X0 = QtWidgets.QDoubleSpinBox()
        self.e_X0.setRange(0.3, 4.0)
        self.e_X0.setSingleStep(0.01)
        self.e_X0.setValue(3.048)
        self.e_X0.setSuffix(' m')
        
        self.e_RT = QtWidgets.QDoubleSpinBox()
        self.e_RT.setRange(0.1, 1.0)
        self.e_RT.setSingleStep(0.01)
        self.e_RT.setValue(0.3048)
        self.e_RT.setSuffix(' m')
        
        self.e_ang = QtWidgets.QDoubleSpinBox()
        self.e_ang.setRange(5, 45)
        self.e_ang.setValue(22.0)
        self.e_ang.setSuffix(' deg')
        
        self.cb_geom_mode = QtWidgets.QComboBox()
        self.cb_geom_mode.addItems(['Calculate X0 from Rt', 'Calculate Rt from X0'])
        
        self.btn_calcX0 = QtWidgets.QPushButton('Calculate geometry')
        
        # Coil/nozzle
        self.e_RH0 = QtWidgets.QDoubleSpinBox()
        self.e_RH0.setRange(80, 100)
        self.e_RH0.setValue(95.0)
        self.e_RH0.setSuffix(' %')
        
        # Inverse results (read-only)
        self.e_T0_inv = QtWidgets.QDoubleSpinBox()
        self.e_T0_inv.setRange(-50, 80)
        self.e_T0_inv.setDecimals(2)
        self.e_T0_inv.setReadOnly(True)
        self.e_T0_inv.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.e_T0_inv.setSuffix(' C')
        
        self.e_P0_inv = QtWidgets.QDoubleSpinBox()
        self.e_P0_inv.setRange(0, 100)
        self.e_P0_inv.setDecimals(3)
        self.e_P0_inv.setReadOnly(True)
        self.e_P0_inv.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.e_P0_inv.setSuffix(' mmHg')
        
        self.e_V0_inv = QtWidgets.QDoubleSpinBox()
        self.e_V0_inv.setRange(0, 100)
        self.e_V0_inv.setDecimals(3)
        self.e_V0_inv.setReadOnly(True)
        self.e_V0_inv.setButtonSymbols(QtWidgets.QAbstractSpinBox.NoButtons)
        self.e_V0_inv.setSuffix(' m/s')
        
        # Other
        self.e_PATM = QtWidgets.QDoubleSpinBox()
        self.e_PATM.setRange(70, 110)
        self.e_PATM.setValue(101.325)
        self.e_PATM.setDecimals(3)
        self.e_PATM.setSuffix(' kPa')
        
        self.chk_buoy = QtWidgets.QCheckBox('Include buoyancy (Ar)')
        
        # Buttons
        self.btn_run = QtWidgets.QPushButton('Run & report')
        self.btn_use_selected = QtWidgets.QPushButton('Use selected point')
        self.btn_use_selected.setEnabled(False)
        self.btn_pdf = QtWidgets.QPushButton('Export PDF...')
        self.btn_cfd = QtWidgets.QPushButton('Run CFD Simulation (2D LBM)')
        
        # Labels
        self.lbl_selected = QtWidgets.QLabel('Selected from chart: None')
    
    def _populate_form(self, form):
        """Add controls to form layout"""
        # Ambient
        form.addRow(QtWidgets.QLabel('<b>Ambient</b>'))
        form.addRow('T_A:', self.e_TA)
        form.addRow('RH_A:', self.e_RHA)
        form.addRow('T_mr:', self.e_TMR)
        form.addRow('V_j:', self.e_VJ)
        
        # Work & clothing
        form.addRow(QtWidgets.QLabel('<b>Work & clothing</b>'))
        form.addRow('M:', self.e_M)
        form.addRow('I_cl:', self.e_ICL)
        
        # Jet geometry
        form.addRow(QtWidgets.QLabel('<b>Jet geometry</b>'))
        form.addRow('D_0:', self.e_D0)
        form.addRow('X_0:', self.e_X0)
        form.addRow('R_t:', self.e_RT)
        form.addRow('Angle:', self.e_ang)
        form.addRow('Geometry mode:', self.cb_geom_mode)
        form.addRow(self.btn_calcX0)
        
        # Coil/nozzle
        form.addRow(QtWidgets.QLabel('<b>Coil/nozzle</b>'))
        form.addRow('rh_0:', self.e_RH0)
        
        # Inverse
        form.addRow(QtWidgets.QLabel('<b>From selected point (inverse)</b>'))
        form.addRow('T_0*:', self.e_T0_inv)
        form.addRow('P_0*:', self.e_P0_inv)
        form.addRow('V_0*:', self.e_V0_inv)
        
        # Other
        form.addRow(QtWidgets.QLabel('<b>Other</b>'))
        form.addRow('P_atm:', self.e_PATM)
        form.addRow(self.chk_buoy)
        form.addRow(self.btn_run)
        form.addRow(self.btn_use_selected)
        form.addRow(self.btn_pdf)
        form.addRow(self.btn_cfd)
        
        # Psychrometric
        form.addRow(QtWidgets.QLabel('<b>Psychrometric</b>'))
        form.addRow(self.lbl_selected)
    
    def _connect_signals(self):
        """Connect button and value change signals"""
        self.btn_calcX0.clicked.connect(self.calc_geometry)
        self.cb_geom_mode.currentIndexChanged.connect(self.on_geometry_mode_changed)
        self.btn_run.clicked.connect(self.run_calc)
        self.btn_use_selected.clicked.connect(self.use_selected_point)
        self.btn_pdf.clicked.connect(self.export_pdf)
        self.btn_cfd.clicked.connect(self.run_cfd_simulation)
        
        # Chart update on parameter changes
        self.e_TA.valueChanged.connect(self.on_chart_params_changed)
        self.e_RHA.valueChanged.connect(self.on_chart_params_changed)
        self.e_VJ.valueChanged.connect(self.on_chart_params_changed)
        self.e_RH0.valueChanged.connect(self.on_chart_params_changed)
        self.chk_buoy.toggled.connect(self.on_chart_params_changed)
        self.e_TMR.valueChanged.connect(self.on_chart_params_changed)
        self.e_M.valueChanged.connect(self.on_chart_params_changed)
        self.e_ICL.valueChanged.connect(self.on_chart_params_changed)
        self.e_PATM.valueChanged.connect(self.on_chart_params_changed)
        
        self.on_geometry_mode_changed()
    
    def on_geometry_mode_changed(self):
        """Switch geometry calculation mode"""
        mode = self.cb_geom_mode.currentIndex() if hasattr(self, 'cb_geom_mode') else 0
        calc_x0 = (mode == 0)
        self.e_X0.setReadOnly(calc_x0)
        self.e_RT.setReadOnly(not calc_x0)
        self.e_X0.setButtonSymbols(
            QtWidgets.QAbstractSpinBox.NoButtons if calc_x0
            else QtWidgets.QAbstractSpinBox.UpDownArrows
        )
        self.e_RT.setButtonSymbols(
            QtWidgets.QAbstractSpinBox.NoButtons if not calc_x0
            else QtWidgets.QAbstractSpinBox.UpDownArrows
        )
    
    def calc_geometry(self):
        """Calculate jet geometry (X0 or Rt)"""
        Rt = self.e_RT.value()
        D0 = self.e_D0.value()
        ang_deg = self.e_ang.value()
        tan_half = math.tan(math.radians(ang_deg / 2.0))
        
        if tan_half <= 1e-12:
            QtWidgets.QMessageBox.warning(
                self, 'Geometry',
                'Angle must be greater than 0° for geometry calculation.'
            )
            return
        
        if self.cb_geom_mode.currentIndex() == 0:
            # Calculate X0 from Rt
            num = Rt - 0.5 * D0
            if num <= 0.0:
                QtWidgets.QMessageBox.warning(
                    self, 'Geometry',
                    'Rt must be greater than D0/2 to compute X0.'
                )
                return
            X0 = num / tan_half
            self.e_X0.setValue(X0)
        else:
            # Calculate Rt from X0
            X0 = self.e_X0.value()
            Rt_calc = X0 * tan_half + 0.5 * D0
            if Rt_calc <= 0:
                QtWidgets.QMessageBox.warning(
                    self, 'Geometry',
                    'Calculated Rt is not physically valid.'
                )
                return
            self.e_RT.setValue(Rt_calc)
    
    def on_chart_params_changed(self):
        """Debounce chart updates when input parameters change"""
        if hasattr(self, '_chart_update_timer'):
            self._chart_update_timer.start()
    
    def _apply_chart_update(self):
        """Apply deferred chart update when chart parameters change.
        
        Only clears operating solution if parameters that affect the solution changed
        (Icl, M, Tmr, Vj). Other parameters (TA, RH, PATM) don't invalidate the solution.
        This prevents the graph from disappearing when user makes minor adjustments.
        """
        if hasattr(self, 'chart'):
            # Check if parameters changed that would invalidate current solution
            icl_changed = abs(self.e_ICL.value() - self._prev_icl) > 1e-6
            m_changed = abs(self.e_M.value() - self._prev_m) > 1e-6
            tmr_changed = abs(self.e_TMR.value() - self._prev_tmr) > 1e-6
            vj_changed = abs(self.e_VJ.value() - self._prev_vj) > 1e-6
            
            # Update stored values
            self._prev_icl = self.e_ICL.value()
            self._prev_m = self.e_M.value()
            self._prev_tmr = self.e_TMR.value()
            self._prev_vj = self.e_VJ.value()
            
            # If any critical parameter changed, solution is invalid
            # (Icl, M, Tmr, Vj all affect the acceptable line or solution)
            if icl_changed or m_changed or tmr_changed or vj_changed:
                self.chart.clear_operating_solution()
            
            # Update chart parameters and redraw
            self.chart.Icl = self.e_ICL.value()
            self.chart.M = self.e_M.value()
            self.chart.Tmr = self.e_TMR.value()
            self.chart.p_atm_kpa = self.e_PATM.value()
            self.chart.update_chart()
    
    def on_chart_selected(self, T_j, W_j):
        """Handle selection from psychrometric chart"""
        self.chart.selected_Tj_point = (T_j, W_j)
        self.lbl_selected.setText(f'Selected: T_j={T_j:.1f}°C, W_j={W_j*1000:.2f} g/kg')
        self.btn_use_selected.setEnabled(True)
    
    def on_invalid_chart_selection(self, msg):
        """Handle invalid point selection from chart"""
        self.chart.clear_selected_overlay()
        self.lbl_selected.setText(f'Selected from chart: None ({msg})')
        self.btn_use_selected.setEnabled(False)
        QtWidgets.QMessageBox.warning(self, 'Invalid selection', msg)
    
    def _is_physical_solution(self, result, TA=None):
        """Validate solution is physically consistent"""
        vals = [
            result.get('T0'), result.get('Tj'), result.get('P0'),
            result.get('Pj'), result.get('RH_0'), result.get('RH_j')
        ]
        if any((v is None or not np.isfinite(v)) for v in vals):
            return False
        if result['P0'] <= 0.0 or result['Pj'] <= 0.0:
            return False
        if not (0.0 <= result['RH_0'] <= 1.0 and 0.0 <= result['RH_j'] <= 1.0):
            return False
        
        if TA is not None and np.isfinite(TA):
            t_min = min(TA, result['T0']) - 1e-6
            t_max = max(TA, result['T0']) + 1e-6
            if not (t_min <= result['Tj'] <= t_max):
                return False
        elif result['Tj'] < result['T0'] - 1e-6:
            return False
        
        return True
    
    def use_selected_point(self):
        """Calculate from selected point on chart"""
        if not hasattr(self.chart, 'selected_Tj_point') or self.chart.selected_Tj_point is None:
            QtWidgets.QMessageBox.warning(self, 'Error', 'No point selected on chart')
            return
        
        T_j_selected, W_j_selected = self.chart.selected_Tj_point
        
        try:
            from controllers.main_controller import MainController
            controller = MainController(self)
            result = controller.solve_from_selected_point(T_j_selected, W_j_selected)
            
            if not self._is_physical_solution(result, TA=self.e_TA.value()):
                self.chart.clear_operating_solution()
                self.chart.update_chart()
                QtWidgets.QMessageBox.warning(
                    self, 'Invalid selection',
                    'Selected point does not produce a physically feasible solution.'
                )
                return
            
            # Update UI with results
            self.e_RH0.setValue(max(self.e_RH0.minimum(),
                                   min(self.e_RH0.maximum(), result['RH_0'] * 100.0)))
            self.e_T0_inv.setValue(result['T0'])
            self.e_P0_inv.setValue(result['P0'])
            self.e_V0_inv.setValue(result['V0'])
            
            self.chart.set_operating_solution(result, self.e_VJ.value(), source='selected')
            self.chart.update_chart()
            
            html = self.build_html(result)
            self.last_report_html = html
            self.show_report_dialog(html)
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, 'Error', str(e))
    
    def run_calc(self):
        """Run calculation with current parameters"""
        try:
            self.chart.clear_selected_overlay()
            self.lbl_selected.setText('Selected from chart: None')
            self.btn_use_selected.setEnabled(False)
            
            from controllers.main_controller import MainController
            controller = MainController(self)
            result = controller.solve_normal()
            
            if not self._is_physical_solution(result, TA=self.e_TA.value()):
                self.chart.clear_operating_solution()
                self.chart.update_chart()
                QtWidgets.QMessageBox.warning(
                    self, 'Invalid solution',
                    'Run result is not physically feasible.'
                )
                return
            
            self.chart.set_operating_solution(result, self.e_VJ.value(), source='run')
            self.chart.update_chart()
            
            html = self.build_html(result)
            self.last_report_html = html
            self.show_report_dialog(html)
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, 'Error', str(e))
    
    def build_html(self, result):
        """Build HTML report from results"""
        def fmt(x, d=3):
            try:
                return f"{x:.{d}f}"
            except Exception:
                return str(x)
        
        html = """
        <html><head><meta charset='utf-8'>
        <style>body{font-family:Arial;margin:20px}table{border-collapse:collapse;width:100%}
        th,td{border:1px solid #ccc;padding:6px;text-align:right}th{background:#f2f2f2}
        .left{text-align:left}</style>
        </head><body>
        <h1>Spot Cooling Report</h1>
        <h2>Results</h2>
        <table>
        <tr><th class='left'>Variable</th><th>Value</th><th>Units</th></tr>
        """
        
        rows = [
            ('V_j (used)', fmt(self.e_VJ.value(), 3), 'm/s'),
            ('m (calculated)', fmt(result['m'], 4), 'mmHg/C'),
            ('C (calculated)', fmt(result['C'], 3), 'mmHg'),
            ('T_a(0.5)', fmt(result['Ta50'], 2), 'C'),
            ('T_0', fmt(result['T0'], 3), 'C'),
            ('P_0', fmt(result['P0'], 3), 'mmHg'),
            ('RH_0', fmt(result['RH_0'] * 100, 2), '%'),
            ('T_j', fmt(result['Tj'], 3), 'C'),
            ('P_j', fmt(result['Pj'], 3), 'mmHg'),
            ('RH_j', fmt(result['RH_j'] * 100, 2), '%'),
            ('V_0', fmt(result['V0'], 3), 'm/s'),
            ('V_j/V_0', fmt(result['Vratio'], 4), '-'),
            ('(T_A-T_j)/(T_A-T_0)', fmt(result['Tratio'], 4), '-'),
            ('Q_0', fmt(result['Q0'], 4), 'm3/s'),
            ('Q_j', fmt(result['Qj'], 4), 'm3/s'),
            ('Entrainment', fmt(result['Qe'], 4), 'm3/s'),
            ('m_dot0', fmt(result['m_dot0'], 4), 'kg/s'),
            ('Q_total', fmt(result['Q_total'], 3), 'kW'),
            ('Q_sensible', fmt(result['Q_sens'], 3), 'kW'),
        ]
        
        for k, v, u in rows:
            html += f"<tr><td class='left'>{k}</td><td>{v}</td><td>{u}</td></tr>"
        
        html += "</table></body></html>"
        return html
    
    def show_report_dialog(self, html):
        """Show report in popup dialog"""
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle('Spot Cooling Report')
        dlg.resize(700, 600)
        layout = QtWidgets.QVBoxLayout(dlg)
        view = QtWidgets.QTextBrowser()
        view.setHtml(html)
        layout.addWidget(view)
        btn_export = QtWidgets.QPushButton('Export as PDF')
        btn_export.clicked.connect(lambda: self.export_pdf_from_dialog(html))
        layout.addWidget(btn_export)
        dlg.exec_()
    
    def export_pdf_from_dialog(self, html):
        """Export report dialog to PDF"""
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Export to PDF', 'spot_cooling_report.pdf', 'PDF (*.pdf)'
        )
        if not path:
            return
        
        printer = QtPrintSupport.QPrinter(QtPrintSupport.QPrinter.HighResolution)
        printer.setOutputFormat(QtPrintSupport.QPrinter.PdfFormat)
        printer.setOutputFileName(path)
        printer.setPageMargins(12, 12, 12, 12, QtPrintSupport.QPrinter.Millimeter)
        doc = QtGui.QTextDocument()
        doc.setHtml(html)
        doc.print_(printer)
        QtWidgets.QMessageBox.information(self, 'Success', f'PDF exported to:\n{path}')
    
    def export_pdf(self):
        """Legacy export PDF button"""
        if self.last_report_html is None:
            QtWidgets.QMessageBox.warning(
                self, 'Warning',
                'No report generated yet. Run a calculation first.'
            )
            return
        self.export_pdf_from_dialog(self.last_report_html)
    
    def run_cfd_simulation(self):
        """Execute 2D LBM CFD simulation of the jet"""
        try:
            # First, run normal calculation to get the solution
            from controllers.main_controller import MainController
            controller = MainController(self)
            result = controller.solve_normal()
            
            if not self._is_physical_solution(result, TA=self.e_TA.value()):
                QtWidgets.QMessageBox.warning(
                    self, 'Invalid solution',
                    'Current parameters do not produce a physical solution. '
                    'Please adjust parameters and try again.'
                )
                return
            
            # Now create CFD simulation from the result
            from controllers.cfd_controller import create_cfd_from_cooling_solution
            from views.cfd_visualizer import CFDResultsDialog
            
            # Show progress dialog while running simulation
            progress = QtWidgets.QProgressDialog(
                'Running 2D Lattice Boltzmann Simulation...', 
                'Cancel', 0, 100, self
            )
            progress.setWindowModality(QtCore.Qt.WindowModal)
            progress.setAutoClose(True)
            progress.setAutoReset(True)
            
            def update_progress(step, total):
                progress.setValue(int(100 * step / total))
                QtWidgets.QApplication.processEvents()
            
            cfd_controller = create_cfd_from_cooling_solution(self, result)
            
            # Run simulation
            results = cfd_controller.run_simulation(n_steps=500, callback=update_progress)
            
            # Show results dialog
            dlg = CFDResultsDialog(self)
            dlg.set_results(
                velocity=results['velocity'],
                ux=results['ux'],
                uy=results['uy'],
                rho=results['rho'],
                reynolds=results['reynolds'],
                mach=results['mach'],
                jet_diameter_physical=results['jet_diameter'],
                lattice_spacing=0.001  # 1mm lattice spacing
            )
            
            progress.close()
            dlg.exec_()
            
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, 'CFD Error', f'Simulation failed:\n{str(e)}')
