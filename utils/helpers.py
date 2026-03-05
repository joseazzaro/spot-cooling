# -*- coding: utf-8 -*-
"""
Helper utilities for Spot Cooling Designer
"""

def bracket_value(value, grid):
    """
    Find bracket indices and interpolation factor for value within grid.
    Returns (lo, hi, t) where lo and hi are grid values and t is interpolation factor [0,1]
    """
    if value <= grid[0]: 
        return grid[0], grid[0], 0.0
    if value >= grid[-1]: 
        return grid[-1], grid[-1], 0.0
    
    for i in range(len(grid)-1):
        lo, hi = grid[i], grid[i+1]
        if lo <= value <= hi:
            t = 0.0 if hi == lo else (value - lo) / (hi - lo)
            return lo, hi, t
    
    return grid[-2], grid[-1], 1.0


def linear_interpolate(a, b, t):
    """Linear interpolation: result = a + (b-a)*t"""
    return a + (b - a) * t


def is_finite_positive(value):
    """Check if value is finite and positive"""
    try:
        return float(value) > 0 and float(value) == float(value)  # NaN check
    except (TypeError, ValueError):
        return False
