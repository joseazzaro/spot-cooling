from __future__ import annotations
from typing import Optional, Union, Tuple
import numpy as np

ArrayLike = Union[float, np.ndarray]

# --- Optional deps ---
try:
    import CoolProp.CoolProp as CP  # type: ignore
    _HAS_COOLPROP = True
except Exception:
    _HAS_COOLPROP = False

try:
    from iapws import IAPWS97  # type: ignore
    _HAS_IAPWS = True
except Exception:
    _HAS_IAPWS = False

# Try to import IAPWS-95 for subfreezing liquid water (< triple point)
try:
    from iapws import IAPWS95  # type: ignore
    _HAS_IAPWS95 = True
except Exception:
    _HAS_IAPWS95 = False

# Policies
DEFAULT_PRIMARY = "IAPWS"
DEFAULT_COOLPROP_THRESHOLD = 10000
_TRIPLE_K = 273.16  # triple point temperature (K)


def _to_np(x) -> np.ndarray:
    return np.atleast_1d(np.asarray(x, dtype=float))


def _choose_backend(n_points: Optional[int], prefer: Optional[str] = None,
                    threshold: int = DEFAULT_COOLPROP_THRESHOLD) -> str:
    """Choose computational backend for IF97-valid region.

    Preference order:
    - explicit 'prefer' if available and installed
    - if data size >= threshold and CoolProp is available -> 'COOLPROP'
    - otherwise fall back to IAPWS if available
    - finally CoolProp if available
    """
    if prefer:
        p = prefer.upper()
        if p.startswith("COOL") and _HAS_COOLPROP:
            return "COOLPROP"
        if p.startswith("IAP") and _HAS_IAPWS:
            return "IAPWS"
    n = 1 if n_points is None else int(n_points)
    if n >= threshold and _HAS_COOLPROP:
        return "COOLPROP"
    if _HAS_IAPWS:
        return "IAPWS"
    if _HAS_COOLPROP:
        return "COOLPROP"
    raise RuntimeError("No se encontró backend: instale 'iapws' o 'CoolProp'.")


def _needs_subfreezing_liquid(T_K_arr: np.ndarray) -> bool:
    """True si alguna T está por debajo del triple point (zona no cubierta por IF97)."""
    return bool(np.any(np.asarray(T_K_arr, dtype=float) < _TRIPLE_K))


# ---------- Sublimación sobre hielo (Murphy & Koop 2005) ----------
# ln(p/Pa) = 9.550426 - 5723.265/T + 3.53068*ln(T) - 0.00728332*T
# Válida aprox. 110 K – 273.16 K

def _pws_ice_MK(T_K: np.ndarray) -> np.ndarray:
    T = np.asarray(T_K, dtype=float)
    return np.exp(9.550426 - 5723.265/T + 3.53068*np.log(T) - 0.00728332*T)


# ------------------------- Public API ----------------------------

def pws_T(T_C: ArrayLike, prefer: Optional[str] = None,
          threshold: int = DEFAULT_COOLPROP_THRESHOLD) -> ArrayLike:
    """Presión de saturación (Pa) vs T (°C).
    T<0°C: sublimación (Murphy & Koop 2005). T>=0°C: saturación LV (IF97).
    Devuelve escalar si la entrada era escalar, o ndarray en caso contrario.
    """
    T_arr = _to_np(T_C)
    T_K = T_arr + 273.15
    backend = _choose_backend(T_arr.size, prefer, threshold)
    out = np.empty_like(T_arr)
    mask_liq = T_arr >= 0.0
    mask_ice = ~mask_liq
    if backend == "IAPWS":
        if np.any(mask_liq):
            Tk = T_K[mask_liq]
            out[mask_liq] = np.array([IAPWS97(T=float(t), x=0).P for t in Tk], dtype=float) * 1e6
        if np.any(mask_ice):
            out[mask_ice] = _pws_ice_MK(T_K[mask_ice])
    else:
        if np.any(mask_liq):
            out[mask_liq] = CP.PropsSI("P", "T", T_K[mask_liq], "Q", 0, "IF97::Water")
        if np.any(mask_ice):
            out[mask_ice] = _pws_ice_MK(T_K[mask_ice])
    return float(out[0]) if np.ndim(T_C) == 0 else out


def tws_P(p_Pa: ArrayLike, prefer: Optional[str] = None,
          threshold: int = DEFAULT_COOLPROP_THRESHOLD) -> ArrayLike:
    """Temperatura de saturación (K) vs p (Pa) (líquido-vapor IF97).
    *No* invierte la curva de sublimación.
    """
    p_arr = _to_np(p_Pa)
    backend = _choose_backend(p_arr.size, prefer, threshold)
    if backend == "IAPWS":
        T_K = np.array([IAPWS97(P=float(p)/1e6, x=0).T for p in p_arr], dtype=float)
    else:
        T_K = np.array(CP.PropsSI("T", "P", p_arr, "Q", 0, "IF97::Water"), dtype=float)
    return float(T_K[0]) if np.ndim(p_Pa) == 0 else T_K


def water_specific_volume(T_C: ArrayLike, P_kPa: ArrayLike,
                          prefer: Optional[str] = None,
                          threshold: int = DEFAULT_COOLPROP_THRESHOLD) -> ArrayLike:
    """Específico v (m³/kg) del agua a T (°C) y P (kPa)."""
    T_arr, P_arr = _to_np(T_C), _to_np(P_kPa)
    T_K = T_arr + 273.15
    P_Pa = P_arr * 1000.0

    # --- helper: detectar si estamos (muy) cerca de ps(T) ---
    # usaremos el mismo pws_T del backend para comparar
    def _near_sat(Tk_vec, Ppa_vec, rtol=1e-6, atol=5.0):
        ps = pws_T(Tk_vec - 273.15)  # pws en Pa
        return np.isfinite(ps) & (np.abs(Ppa_vec - ps) <= np.maximum(atol, rtol * np.maximum(Ppa_vec, ps)))

    near_sat = _near_sat(T_K, P_Pa)

    # --- Subcero: preferir HEOS, luego IAPWS95 ---
    if _needs_subfreezing_liquid(T_K):
        if _HAS_COOLPROP:
            try:
                # si estamos en saturación: usar Q=0
                if np.any(near_sat):
                    v = np.empty_like(T_arr)
                    idx = near_sat
                    v[idx] = np.array(CP.PropsSI("V", "T", T_K[idx], "Q", 0, "Water"), dtype=float)
                    idx = ~near_sat
                    if np.any(idx):
                        v[idx] = np.array(CP.PropsSI("V", "T", T_K[idx], "P", P_Pa[idx], "Water"), dtype=float)
                else:
                    v = np.array(CP.PropsSI("V", "T", T_K, "P", P_Pa, "Water"), dtype=float)
                return float(v[0]) if np.ndim(T_C) == 0 and np.ndim(P_kPa) == 0 else v
            except Exception:
                # Fallback a IAPWS95
                if _HAS_IAPWS95:
                    v = np.empty_like(T_arr)
                    for i, (t, p) in enumerate(zip(T_K, P_Pa)):
                        st = IAPWS95(T=float(t), P=float(p)/1e6)  # P en MPa
                        v[i] = st.v
                    return float(v[0]) if np.ndim(T_C) == 0 and np.ndim(P_kPa) == 0 else v
                raise
        elif _HAS_IAPWS95:
            v = np.empty_like(T_arr)
            for i, (t, p) in enumerate(zip(T_K, P_Pa)):
                st = IAPWS95(T=float(t), P=float(p)/1e6)
                v[i] = st.v
            return float(v[0]) if np.ndim(T_C) == 0 and np.ndim(P_kPa) == 0 else v
        else:
            raise RuntimeError("Se requieren CoolProp o IAPWS95 para T < 0°C (líquido subenfriado).")

    # --- Zona IF97 válida ---
    backend = _choose_backend(max(T_arr.size, P_arr.size), prefer, threshold)
    if backend == "IAPWS":
        # IAPWS97 maneja saturación con x=0 naturalmente
        try:
            v = np.empty_like(T_arr)
            for i, (t, p) in enumerate(zip(T_K, P_Pa)):
                if _near_sat(np.array([t]), np.array([p]))[0]:
                    st = IAPWS97(T=float(t), x=0)   # líquido saturado
                else:
                    st = IAPWS97(T=float(t), P=float(p)/1e6)
                v[i] = st.v
        except Exception:
            # Fallback a CoolProp
            if _HAS_COOLPROP:
                pass  # seguimos con la rama COOLPROP de abajo
            else:
                raise
        else:
            return float(v[0]) if np.ndim(T_C) == 0 and np.ndim(P_kPa) == 0 else v

    # --- COOLPROP (IF97::Water), con manejo de saturación ---
    if not _HAS_COOLPROP:
        raise RuntimeError("CoolProp no disponible y no fue posible usar IAPWS97.")

    try:
        if np.any(near_sat):
            v = np.empty_like(T_arr)
            idx = near_sat
            # pedir con calidad (líquido saturado)
            try:
                v[idx] = np.array(CP.PropsSI("V", "T", T_K[idx], "Q", 0, "IF97::Water"), dtype=float)
            except Exception:
                # si IF97 no deja, usar HEOS
                v[idx] = np.array(CP.PropsSI("V", "T", T_K[idx], "Q", 0, "Water"), dtype=float)
            # el resto fuera de saturación: (T,P)
            idx = ~near_sat
            if np.any(idx):
                try:
                    v[idx] = np.array(CP.PropsSI("V", "T", T_K[idx], "P", P_Pa[idx], "IF97::Water"), dtype=float)
                except Exception:
                    v[idx] = np.array(CP.PropsSI("V", "T", T_K[idx], "P", P_Pa[idx], "Water"), dtype=float)
        else:
            # todo fuera de saturación
            try:
                v = np.array(CP.PropsSI("V", "T", T_K, "P", P_Pa, "IF97::Water"), dtype=float)
            except Exception:
                v = np.array(CP.PropsSI("V", "T", T_K, "P", P_Pa, "Water"), dtype=float)
    except Exception:
        # último recurso: pequeña perturbación de P para alejar ~1 Pa de la curva
        dP = np.maximum(P_Pa * 1e-8, 1.0)
        v = np.array(CP.PropsSI("V", "T", T_K, "P", P_Pa + dP, "Water"), dtype=float)

    return float(v[0]) if np.ndim(T_C) == 0 and np.ndim(P_kPa) == 0 else v


def water_enthalpy_liq(T_C: ArrayLike, P_kPa: ArrayLike | None = None,
                        prefer: Optional[str] = None,
                        threshold: int = DEFAULT_COOLPROP_THRESHOLD) -> ArrayLike:
    """Entalpía del agua líquida (kJ/kg) a T (°C); si se provee P (kPa) usa estado (T,P),
    de lo contrario asume saturación líquida a esa T para consistencia con psicrometría.
    """
    T_arr = _to_np(T_C)
    if P_kPa is None:
        P_Pa = pws_T(T_arr, prefer=prefer, threshold=threshold)
        P_kPa_arr = P_Pa / 1000.0
    else:
        P_kPa_arr = _to_np(P_kPa)
    T_K = T_arr + 273.15
    P_Pa = P_kPa_arr * 1000.0

    # Subcero: preferir HEOS, luego IAPWS95
    if _needs_subfreezing_liquid(T_K):
        if _HAS_COOLPROP:
            h = np.array(CP.PropsSI("H", "T", T_K, "P", P_Pa, "Water"), dtype=float) / 1000.0
            return float(h[0]) if np.ndim(T_C) == 0 else h
        elif _HAS_IAPWS95:
            h = np.empty_like(T_arr)
            for i, (t, p) in enumerate(zip(T_K, P_Pa)):
                st = IAPWS95(T=float(t), P=float(p)/1e6)
                h[i] = st.h  # kJ/kg
            return float(h[0]) if np.ndim(T_C) == 0 else h
        else:
            raise RuntimeError(
                "Se requieren CoolProp o 'iapws' con IAPWS95 para T < 0°C (líquido subenfriado)."
            )

    # Zona IF97 válida
    backend = _choose_backend(max(T_arr.size, np.atleast_1d(P_kPa_arr).size), prefer, threshold)
    if backend == "IAPWS":
        h = np.empty_like(T_arr)
        for i, (t, p) in enumerate(zip(T_K, P_Pa)):
            st = IAPWS97(T=float(t), P=float(p)/1e6)
            h[i] = st.h  # kJ/kg
        return float(h[0]) if np.ndim(T_C) == 0 else h
    else:
        h = np.array(CP.PropsSI("H", "T", T_K, "P", P_Pa, "IF97::Water"), dtype=float) / 1000.0
        return float(h[0]) if np.ndim(T_C) == 0 else h


def isothermal_compressibility(
    T_C: ArrayLike, P_kPa: ArrayLike,
    prefer: Optional[str] = None,
    threshold: int = DEFAULT_COOLPROP_THRESHOLD
) -> ArrayLike:
    """
    Compresibilidad isotérmica k_T (1/Pa) del agua a T (°C), P (kPa).
    - Subcero: HEOS (Water) directo; o derivada con IAPWS95 si no hay CoolProp.
    - IF97 válido: intentar primero IF97 directo; si falla, probar HEOS; si aún falla, derivada numérica con HEOS.
    - Siempre existe fallback a derivada numérica vía V para máxima robustez.
    """
    T_arr = _to_np(T_C)
    P_arr = _to_np(P_kPa)
    T_K   = T_arr + 273.15
    P_Pa  = P_arr * 1000.0

    # --- Subcero: preferir HEOS (Water), luego derivada con IAPWS95 ---
    if _needs_subfreezing_liquid(T_K):
        if _HAS_COOLPROP:
            try:
                kT = np.array(
                    CP.PropsSI("isothermal_compressibility", "T", T_K, "P", P_Pa, "Water"),
                    dtype=float
                )
                return float(kT[0]) if np.ndim(T_C) == 0 and np.ndim(P_kPa) == 0 else kT
            except Exception:
                # Fallback numérico con HEOS
                v0 = np.array(CP.PropsSI("V", "T", T_K, "P", P_Pa, "Water"), dtype=float)
                dP = np.maximum(P_Pa * 1e-6, 5.0)  # Pa
                v_plus = np.array(CP.PropsSI("V", "T", T_K, "P", P_Pa + dP, "Water"), dtype=float)
                dv_dP = (v_plus - v0) / dP
                kT = - (1.0 / v0) * dv_dP
                return float(kT[0]) if np.ndim(T_C) == 0 and np.ndim(P_kPa) == 0 else kT
        elif _HAS_IAPWS95:
            # Derivada numérica con IAPWS95
            v0 = water_specific_volume(T_arr, P_arr, prefer=prefer, threshold=threshold)
            dP = np.maximum(P_Pa * 1e-6, 5.0)  # Pa
            v_plus = water_specific_volume(T_arr, (P_Pa + dP)/1000.0, prefer=prefer, threshold=threshold)
            dv_dP = (v_plus - v0) / dP
            kT = - (1.0 / v0) * dv_dP
            return float(kT[0]) if np.ndim(T_C) == 0 and np.ndim(P_kPa) == 0 else kT
        else:
            raise RuntimeError("Se requieren CoolProp o IAPWS95 para k_T a T < 0°C (líquido subenfriado).")

    # --- Zona IF97 válida ---
    backend = _choose_backend(max(T_arr.size, P_arr.size), prefer, threshold)
    if backend == "IAPWS":
        # Derivada numérica con IAPWS97
        v0 = water_specific_volume(T_arr, P_arr, prefer=prefer, threshold=threshold)
        dP = np.maximum(P_Pa * 1e-6, 5.0)
        v_plus = water_specific_volume(T_arr, (P_Pa + dP)/1000.0, prefer=prefer, threshold=threshold)
        dv_dP = (v_plus - v0) / dP
        kT = - (1.0 / v0) * dv_dP
        return float(kT[0]) if np.ndim(T_C) == 0 and np.ndim(P_kPa) == 0 else kT
    else:
        # CoolProp: primero IF97 directo; si falla, HEOS; si falla, derivada numérica con HEOS
        assert _HAS_COOLPROP, "CoolProp no disponible"
        try:
            kT = np.array(
                CP.PropsSI("isothermal_compressibility", "T", T_K, "P", P_Pa, "IF97::Water"),
                dtype=float
            )
        except Exception:
            try:
                kT = np.array(
                    CP.PropsSI("isothermal_compressibility", "T", T_K, "P", P_Pa, "Water"),
                    dtype=float
                )
            except Exception:
                # Fallback numérico con HEOS
                v0 = np.array(CP.PropsSI("V", "T", T_K, "P", P_Pa, "Water"), dtype=float)
                dP = np.maximum(P_Pa * 1e-6, 5.0)
                v_plus = np.array(CP.PropsSI("V", "T", T_K, "P", P_Pa + dP, "Water"), dtype=float)
                dv_dP = (v_plus - v0) / dP
                kT = - (1.0 / v0) * dv_dP

        return float(kT[0]) if np.ndim(T_C) == 0 and np.ndim(P_kPa) == 0 else kT


def water_kT_and_v(T_C, P_kPa, prefer="COOLPROP", backend="HEOS"):
    """
    Devuelve:
      kT : compresibilidad isotérmica [1/Pa]
      v  : volumen específico [m^3/kg]
    Acepta escalares o arrays en T_C [°C] y P_kPa [kPa].
    """
    # Selección del backend CoolProp
    # - "HEOS::Water" (general)
    # - "IF97::Water" (muy estable en líquido: 0–800°C, 0–1000 bar)
    if backend.upper() in ("HEOS", "IF97"):
        fluid = f"{backend.upper()}::Water"
    else:
        fluid = "Water"  # Default

    T = np.asarray(T_C, dtype=float) + 273.15
    P = np.asarray(P_kPa, dtype=float) * 1000.0  # kPa → Pa

    # Chequeos de sanidad útiles en desarrollo:
    if np.nanmedian(P) > 5e5:  # > 500 kPa
        # Si trabajás con atm de referencia, P_kPa debería ser ~101.3 kPa
        print("[WARN] P_kPa parece estar en Pa (o muy alta). ¿Duplicaste la conversión?")
    if np.any(~np.isfinite(T)) or np.any(~np.isfinite(P)) or np.any(P <= 0.0):
        raise ValueError("T o P contienen valores no finitos o P<=0.")

    # Función escalar segura
    def _one(TK, PP):
        if not np.isfinite(TK) or not np.isfinite(PP) or PP <= 0.0:
            return np.nan, np.nan
        try:
            kT = CP.PropsSI("isothermal_compressibility", "T", float(TK), "P", float(PP), fluid)  # 1/Pa
            rho = CP.PropsSI("D", "T", float(TK), "P", float(PP), fluid)  # kg/m^3
            v = 1.0 / rho
            return kT, v
        except Exception:
            # Guard por región bifásica exacta: probá con un pequeño delta de T
            dT = 1e-3
            try:
                kT = CP.PropsSI("isothermal_compressibility", "T", float(TK) + dT, "P", float(PP), fluid)
                rho = CP.PropsSI("D", "T", float(TK) + dT, "P", float(PP), fluid)
                v = 1.0 / rho
                return kT, v
            except Exception:
                return np.nan, np.nan

    # Aplicar de forma vectorizada (pero evaluando en escalares)
    kT_list, v_list = zip(*[ _one(t, p) for t, p in np.broadcast(T, P) ])
    kT_arr = np.array(kT_list).reshape(T.shape)
    v_arr  = np.array(v_list).reshape(T.shape)
    return kT_arr, v_arr
