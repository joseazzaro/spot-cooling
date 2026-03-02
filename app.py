# -*- coding: utf-8 -*-
# Spot Cooling Designer - single file
import sys, math
from PyQt5 import QtWidgets, QtCore, QtGui, QtPrintSupport

# psychro utils
import math

def psat_mmhg(T):
    return math.exp(18.6686 - 4030.183/(T + 235.0))

def pv_from_rh_T(rh, T):
    return rh * psat_mmhg(T)

def rh_from_pv_T(Pv, T):
    return Pv / psat_mmhg(T)

MMHG_TO_KPA = 0.133322368

def humidity_ratio_from_pv(Pv_mmhg, p_atm_kpa=101.325):
    Pv_kpa = Pv_mmhg * MMHG_TO_KPA
    return 0.62198 * Pv_kpa / max(1e-6, (p_atm_kpa - Pv_kpa))

def moist_air_enthalpy_kJkg(T, w):
    return 1.006*T + w*(2501.0 + 1.86*T)

def air_density_approx(T, w, p_atm_kpa=101.325):
    T_k = T + 273.15
    P_pa = p_atm_kpa*1000.0
    R_d = 287.055
    return P_pa/(R_d*T_k*(1.0 + 1.607*w))

# regression table
ICL_LEVELS = [0.3, 0.6, 0.9]
M_LEVELS   = [87.0, 174.0, 232.0]
V_LEVELS   = [0.5, 1.0, 1.5, 2.0]

TABLE = {
    0.3:{
        87.0:  {0.5:(-0.293,48.1), 1.0:(-0.185,46.2), 1.5:(-0.143,45.7), 2.0:(-0.118,45.2)},
        174.0: {0.5:(-0.263,35.9), 1.0:(-0.211,39.4), 1.5:(-0.136,38.4), 2.0:(-0.118,39.2)},
        232.0: {0.5:(-0.385,32.6), 1.0:(-0.193,32.8), 1.5:(-0.118,32.9), 2.0:(-0.118,35.0)}
    },
    0.6:{
        87.0:  {0.5:(-0.277,45.9), 1.0:(-0.158,43.3), 1.5:(-0.129,43.3), 2.0:(-0.118,43.3)},
        174.0: {0.5:(-0.291,32.0), 1.0:(-0.139,31.5), 1.5:(-0.118,32.8), 2.0:(-0.118,34.9)},
        232.0: {0.5:(-0.505,30.0), 1.0:(-0.248,28.7), 1.5:(-0.185,29.8), 2.0:(-0.118,28.5)}
    },
    0.9:{
        87.0:  {0.5:(-0.248,42.5), 1.0:(-0.191,43.3), 1.5:(-0.118,40.9), 2.0:(-0.100,41.1)},
        174.0: {0.5:(-0.467,34.0), 1.0:(-0.227,31.2), 1.5:(-0.139,29.9), 2.0:(-0.118,30.2)},
        232.0: {1.0:(-0.361,25.9), 1.5:(-0.229,24.4), 2.0:(-0.216,26.4)}
    }
}

def _bracket(value, grid):
    if value <= grid[0]: return grid[0], grid[0], 0.0
    if value >= grid[-1]: return grid[-1], grid[-1], 0.0
    for i in range(len(grid)-1):
        lo, hi = grid[i], grid[i+1]
        if lo <= value <= hi:
            t = 0.0 if hi==lo else (value-lo)/(hi-lo)
            return lo, hi, t
    return grid[-2], grid[-1], 1.0

def regress_Ta50_coeffs(Icl, M, V):
    ic_lo, ic_hi, t_ic = _bracket(Icl, ICL_LEVELS)
    m_lo, m_hi, t_m  = _bracket(M,   M_LEVELS)
    v_lo, v_hi, t_v  = _bracket(V,   V_LEVELS)
    def ab_at(ic, mm, vv):
        tab_ic = TABLE.get(ic, {})
        tab_m  = tab_ic.get(mm, {})
        if vv in tab_m: return tab_m[vv]
        if not tab_m: raise ValueError('Missing table data')
        vv2, ab = min(tab_m.items(), key=lambda kv: abs(kv[0]-vv))
        return ab
    def lerp(a,b,t): return a+(b-a)*t
    a_ll,b_ll = ab_at(ic_lo,m_lo,v_lo); a_lh,b_lh = ab_at(ic_lo,m_lo,v_hi)
    a_hl,b_hl = ab_at(ic_lo,m_hi,v_lo); a_hh,b_hh = ab_at(ic_lo,m_hi,v_hi)
    a_lm,b_lm = lerp(a_ll,a_lh,t_v), lerp(b_ll,b_lh,t_v)
    a_hm,b_hm = lerp(a_hl,a_hh,t_v), lerp(b_hl,b_hh,t_v)
    a_im,b_im = lerp(a_lm,a_hm,t_m), lerp(b_lm,b_hm,t_m)
    a_ll,b_ll = ab_at(ic_hi,m_lo,v_lo); a_lh,b_lh = ab_at(ic_hi,m_lo,v_hi)
    a_hl,b_hl = ab_at(ic_hi,m_hi,v_lo); a_hh,b_hh = ab_at(ic_hi,m_hi,v_hi)
    a_lm2,b_lm2 = lerp(a_ll,a_lh,t_v), lerp(b_ll,b_lh,t_v)
    a_hm2,b_hm2 = lerp(a_hl,a_hh,t_v), lerp(b_hl,b_hh,t_v)
    a_im2,b_im2 = lerp(a_lm2,a_hm2,t_m), lerp(b_lm2,b_hm2,t_m)
    a = lerp(a_im,a_im2,t_ic); b = lerp(b_im,b_im2,t_ic)
    return a,b

# core

def compute_factors(Vj, Icl, Tmr):
    hc = 8.3*(max(1e-6,Vj)**0.6)
    hr = 3.87 + 0.031*Tmr
    h  = hc + hr
    fcl = 1.0 + 0.2*Icl
    Fcl  = 1.0/(1.0 + 0.155*fcl*h*Icl)
    Fpcl = 1.0/(1.0 + 0.143*hc*Icl)
    return hc,hr,h,fcl,Fcl,Fpcl

def ta50_from_regression(Tmr,Icl,M,Vj):
    a,b = regress_Ta50_coeffs(Icl,M,Vj)
    return a*Tmr + b

def acceptable_line(Icl, M, Vj, Tmr):
    W=0.5
    hc,hr,h,fcl,Fcl,Fpcl = compute_factors(Vj,Icl,Tmr)
    m = (h*fcl*Fcl)/(2.2*hc*Fpcl*W)
    Ta50 = ta50_from_regression(Tmr,Icl,M,Vj)
    Ps50 = psat_mmhg(Ta50)
    C = m*Ta50 + 0.5*Ps50
    return m,C,Ta50

def jet_ratios(X0,D0, include_buoyancy, TA, TO, V0_guess=10.0):
    r = (X0/D0) + 2.572
    if not include_buoyancy:
        return 1.464/r, 4.539/r
    beta=1.0/273.15; g=9.81; dT=max(0.0,TA-TO)
    Ar=g*beta*dT*D0/max(1e-6,V0_guess**2)
    Vratio=(1.464/r)*((1.0+0.21*Ar*r*r)**(1.0/3.0))
    Tratio=(4.539/r)*((1.0+0.21*Ar*r*r)**(1.0/3.0))
    return Vratio,Tratio

def solve_case(TA,RH_A,Tmr,Vj,M,Icl,D0,X0, rh0=0.95, p_atm_kpa=101.325, include_buoyancy=False):
    PA = pv_from_rh_T(RH_A,TA)
    m,C,Ta50 = acceptable_line(Icl,M,Vj,Tmr)
    def residual_T0(T0):
        P0 = rh0*psat_mmhg(T0)
        Vratio,Tratio = jet_ratios(X0,D0,False,TA,T0)
        Ti = TA - Tratio*(TA - T0)
        Pi = PA - Tratio*(PA - P0)
        return Pi - (m*Ti + C)
    lo,hi=-5.0,40.0
    def f(x): return residual_T0(x)
    for (a,b) in [(-5.0,40.0),(-20.0,40.0),(-10.0,50.0),(0.0,60.0)]:
        if f(a)*f(b)<=0: lo,hi=a,b; break
    for _ in range(100):
        mid=0.5*(lo+hi); fm=f(mid)
        if abs(fm)<1e-4 or (hi-lo)<1e-4: T0=mid; break
        if f(lo)*fm<=0: hi=mid
        else: lo=mid
    else:
        T0=mid
    P0 = rh0*psat_mmhg(T0)
    Vratio,Tratio = jet_ratios(X0,D0,include_buoyancy,TA,T0)
    Ti = TA - Tratio*(TA - T0)
    Pi = PA - Tratio*(PA - P0)
    RH_i = rh_from_pv_T(Pi,Ti)
    V0 = Vj/max(1e-9,Vratio)
    area0 = math.pi*(D0**2)/4.0
    Q0 = V0*area0
    Qj = Q0/max(1e-9,Tratio)
    Qe = max(0.0, Qj-Q0)
    wA = humidity_ratio_from_pv(PA, p_atm_kpa)
    w0 = humidity_ratio_from_pv(P0, p_atm_kpa)
    hA = moist_air_enthalpy_kJkg(TA,wA)
    h0 = moist_air_enthalpy_kJkg(T0,w0)
    rho0 = air_density_approx(T0,w0,p_atm_kpa)
    m_dot0 = rho0*Q0
    Q_total = m_dot0*(hA - h0)
    Q_sens  = m_dot0*1.006*(TA - T0)
    return dict(m=m, C=C, Ta50=Ta50, PA=PA, Pi=Pi, P0=P0, Ti=Ti, T0=T0,
                RH_i=RH_i, RH_0=rh0, Vratio=Vratio, Tratio=Tratio,
                V0=V0, Q0=Q0, Qj=Qj, Qe=Qe, m_dot0=m_dot0, Q_total=Q_total, Q_sens=Q_sens)

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Spot Cooling - ASHRAE (Azer)')
        self.resize(1100, 720)
        cw=QtWidgets.QWidget(); self.setCentralWidget(cw)
        layout=QtWidgets.QHBoxLayout(cw)
        form=QtWidgets.QFormLayout(); left=QtWidgets.QWidget(); left.setLayout(form)
        scroll=QtWidgets.QScrollArea(); scroll.setWidgetResizable(True); scroll.setWidget(left)
        self.view=QtWidgets.QTextBrowser()
        layout.addWidget(scroll,0); layout.addWidget(self.view,1)
        # inputs
        self.e_TA=QtWidgets.QDoubleSpinBox(); self.e_TA.setRange(-20,80); self.e_TA.setValue(40.0); self.e_TA.setSuffix(' C')
        self.e_RHA=QtWidgets.QDoubleSpinBox(); self.e_RHA.setRange(0,100); self.e_RHA.setValue(50.0); self.e_RHA.setSuffix(' %')
        self.e_TMR=QtWidgets.QDoubleSpinBox(); self.e_TMR.setRange(-20,120); self.e_TMR.setValue(45.0); self.e_TMR.setSuffix(' C')
        self.e_VJ=QtWidgets.QDoubleSpinBox(); self.e_VJ.setRange(0.2,3.0); self.e_VJ.setSingleStep(0.1); self.e_VJ.setValue(2.0); self.e_VJ.setSuffix(' m/s')
        self.e_M=QtWidgets.QDoubleSpinBox(); self.e_M.setRange(40,400); self.e_M.setValue(87.0); self.e_M.setSuffix(' W/m2')
        self.e_ICL=QtWidgets.QDoubleSpinBox(); self.e_ICL.setRange(0.0,2.0); self.e_ICL.setSingleStep(0.05); self.e_ICL.setValue(0.6); self.e_ICL.setSuffix(' clo')
        self.e_D0=QtWidgets.QDoubleSpinBox(); self.e_D0.setRange(0.03,0.5); self.e_D0.setSingleStep(0.005); self.e_D0.setValue(0.127); self.e_D0.setSuffix(' m')
        self.e_X0=QtWidgets.QDoubleSpinBox(); self.e_X0.setRange(0.3,4.0); self.e_X0.setSingleStep(0.01); self.e_X0.setValue(1.259); self.e_X0.setSuffix(' m')
        self.e_RT=QtWidgets.QDoubleSpinBox(); self.e_RT.setRange(0.1,1.0); self.e_RT.setSingleStep(0.01); self.e_RT.setValue(0.3048); self.e_RT.setSuffix(' m')
        self.e_ang=QtWidgets.QDoubleSpinBox(); self.e_ang.setRange(5,45); self.e_ang.setValue(22.0); self.e_ang.setSuffix(' deg')
        self.btn_calcX0=QtWidgets.QPushButton('Calc X0 from Rt & angle')
        self.e_RH0=QtWidgets.QDoubleSpinBox(); self.e_RH0.setRange(80,100); self.e_RH0.setValue(95.0); self.e_RH0.setSuffix(' %')
        self.e_PATM=QtWidgets.QDoubleSpinBox(); self.e_PATM.setRange(70,110); self.e_PATM.setValue(101.325); self.e_PATM.setDecimals(3); self.e_PATM.setSuffix(' kPa')
        self.chk_buoy=QtWidgets.QCheckBox('Include buoyancy (Ar)')
        self.btn_run=QtWidgets.QPushButton('Run & report')
        self.btn_pdf=QtWidgets.QPushButton('Export PDF...')
        # layout
        form.addRow(QtWidgets.QLabel('<b>Ambient</b>'))
        form.addRow('T_A:', self.e_TA); form.addRow('RH_A:', self.e_RHA); form.addRow('T_mr:', self.e_TMR); form.addRow('V_j:', self.e_VJ)
        form.addRow(QtWidgets.QLabel('<b>Work & clothing</b>'))
        form.addRow('M:', self.e_M); form.addRow('I_cl:', self.e_ICL)
        form.addRow(QtWidgets.QLabel('<b>Jet geometry</b>'))
        form.addRow('D_0:', self.e_D0); form.addRow('X_0:', self.e_X0); form.addRow('R_t:', self.e_RT); form.addRow('Angle:', self.e_ang); form.addRow(self.btn_calcX0)
        form.addRow(QtWidgets.QLabel('<b>Coil/nozzle</b>'))
        form.addRow('rh_0:', self.e_RH0)
        form.addRow(QtWidgets.QLabel('<b>Other</b>'))
        form.addRow('P_atm:', self.e_PATM); form.addRow(self.chk_buoy); form.addRow(self.btn_run); form.addRow(self.btn_pdf)
        # signals
        self.btn_calcX0.clicked.connect(self.calc_x0_from_rt)
        self.btn_run.clicked.connect(self.run_calc)
        self.btn_pdf.clicked.connect(self.export_pdf)
        self.view.setHtml('<h3>Fill inputs and press Run & report.</h3>')

    def calc_x0_from_rt(self):
        Rt=self.e_RT.value(); D0=self.e_D0.value(); ang=self.e_ang.value()*math.pi/180.0
        a=0.5/D0; b=math.tan(ang); c=-Rt
        disc=b*b-4*a*c
        if disc<=0:
            QtWidgets.QMessageBox.warning(self,'Geometry','No real solution for X0.')
            return
        x1=(-b+disc**0.5)/(2*a); x2=(-b-disc**0.5)/(2*a)
        X0=max(x1,x2); self.e_X0.setValue(X0)

    def run_calc(self):
        try:
            res = solve_case(self.e_TA.value(), self.e_RHA.value()/100.0, self.e_TMR.value(), self.e_VJ.value(), self.e_M.value(), self.e_ICL.value(), self.e_D0.value(), self.e_X0.value(), rh0=self.e_RH0.value()/100.0, p_atm_kpa=self.e_PATM.value(), include_buoyancy=self.chk_buoy.isChecked())
            self.view.setHtml(self.build_html(res))
        except Exception as e:
            QtWidgets.QMessageBox.critical(self,'Error',str(e))

    def build_html(self, r):
        def fmt(x,d=3):
            try: return f"{x:.{d}f}"
            except: return str(x)
        html = """
        <html><head><meta charset='utf-8'>
        <style>body{font-family:Arial;margin:20px}table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccc;padding:6px;text-align:right}th{background:#f2f2f2}.left{text-align:left}</style>
        </head><body>
        <h1>Spot Cooling Report</h1>
        <h2>Results</h2>
        <table>
        <tr><th class='left'>Var</th><th>Value</th><th>Units</th></tr>
        """
        rows = [
            ('T_0', fmt(r['T0'],2), 'C'),
            ('P_0', fmt(r['P0'],3), 'mmHg'),
            ('RH_0', fmt(r['RH_0']*100,2), '%'),
            ('T_i', fmt(r['Ti'],2), 'C'),
            ('P_i', fmt(r['Pi'],3), 'mmHg'),
            ('RH_i', fmt(r['RH_i']*100,2), '%'),
            ('V_0', fmt(r['V0'],3), 'm/s'),
            ('V_i/V_0', fmt(r['Vratio'],4), '-'),
            ('(T_A-T_i)/(T_A-T_0)', fmt(r['Tratio'],4), '-'),
            ('Q_0', fmt(r['Q0'],4), 'm3/s'),
            ('Q_j', fmt(r['Qj'],4), 'm3/s'),
            ('Entrainment', fmt(r['Qe'],4), 'm3/s'),
            ('m_dot0', fmt(r['m_dot0'],4), 'kg/s'),
            ('Q_total', fmt(r['Q_total'],3), 'kW'),
            ('Q_sensible', fmt(r['Q_sens'],3), 'kW'),
        ]
        for k,v,u in rows:
            html += f"<tr><td class='left'>{k}</td><td>{v}</td><td>{u}</td></tr>"
        html += "</table>"
        html += f"<p>m={fmt(r['m'],4)} mmHg/C; C={fmt(r['C'],3)} mmHg; T_a(0.5)={fmt(r['Ta50'],2)} C.</p>"
        html += "</body></html>"
        return html

    def export_pdf(self):
        path,_=QtWidgets.QFileDialog.getSaveFileName(self,'Export to PDF','spot_cooling_report.pdf','PDF (*.pdf)')
        if not path: return
        printer=QtPrintSupport.QPrinter(QtPrintSupport.QPrinter.HighResolution)
        printer.setOutputFormat(QtPrintSupport.QPrinter.PdfFormat)
        printer.setOutputFileName(path)
        printer.setPageMargins(12,12,12,12, QtPrintSupport.QPrinter.Millimeter)
        doc=QtGui.QTextDocument(); doc.setHtml(self.view.toHtml()); doc.print_(printer)

if __name__=='__main__':
    app=QtWidgets.QApplication(sys.argv)
    w=MainWindow(); w.show()
    sys.exit(app.exec_())
