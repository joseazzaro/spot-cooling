# Spot Cooling Designer

ASHRAE Standard RP-884 Implementation for Spot Cooling Design

## Project Structure (MVC Architecture)

```
spot_cooling/
├── app.py                 # Application entry point
├── requirements.txt       # Python dependencies
│
├── models/                # Business logic and calculations
│   ├── __init__.py
│   ├── psychrometric.py   # Psychrometric calculations for humid air
│   └── cooling_calc.py    # Spot cooling design calculations
│
├── views/                 # PyQt5 GUI components
│   ├── __init__.py
│   ├── main_window.py     # Main application window with controls
│   └── psychro_chart.py   # Psychrometric chart visualization
│
├── controllers/           # Application logic and user interaction
│   ├── __init__.py
│   └── main_controller.py # Main application controller
│
├── utils/                 # Utility modules
│   ├── __init__.py
│   ├── constants.py       # Constants and regression tables
│   └── helpers.py         # Helper functions
│
├── psychro/               # Legacy psychrometric module (optional)
│   └── ...
│
└── virt/                  # Python virtual environment
    └── ...
```

## Architecture Overview

### Models (`models/`)
Contains the core calculation logic separated from the UI:

- **psychrometric.py**: Psychrometric calculations for humid air
  - Saturation vapor pressure
  - Relative humidity conversions
  - Humidity ratio calculations
  - Enthalpy and air density
  - Heat transfer coefficients
  - Clothing and comfort factors

- **cooling_calc.py**: Spot cooling design calculations
  - Regression coefficients from ASHRAE tables
  - Acceptable line parameters (m, C)
  - Jet velocity and temperature ratios
  - Normal solve (target velocity → nozzle conditions)
  - Inverse solve (target point → nozzle conditions)

### Views (`views/`)
PyQt5 GUI components responsible for presentation:

- **main_window.py**: Main application window
  - Form controls for input parameters
  - Control panel layout with scrolling
  - Report generation and display
  - PDF export functionality
  - Event handling for user input

- **psychro_chart.py**: Interactive psychrometric chart
  - Matplotlib-based visualization
  - Acceptable comfort zone lines
  - Psychrometric background curves
  - Operating line display
  - Point selection for inverse solve

### Controllers (`controllers/`)
Business logic that bridges models and views:

- **main_controller.py**: Main application controller
  - Orchestrates calculations
  - Retrieves parameters from UI
  - Calls model methods
  - Handles normal and inverse solving modes

### Utilities (`utils/`)
Shared constants and helper functions:

- **constants.py**: Global constants
  - ASHRAE regression tables (Icl, M, V)
  - Unit conversion factors
  - Psychrometric coefficients
  - Air property constants

- **helpers.py**: Utility functions
  - Value bracketing for interpolation
  - Linear interpolation
  - Validation helpers

## Running the Application

```bash
# Install dependencies
pip install -r requirements.txt

# Run the application
python app.py
```

## Features

- **Interactive Psychrometric Chart**: Click to select target conditions
- **Dual Calculation Modes**:
  - Normal: Input target velocity → get nozzle conditions
  - Inverse: Select point on chart → get required nozzle conditions
- **Geometry Calculator**: Convert between jet spread distance (Rt) and axial distance (X₀)
- **Buoyancy Effects**: Optional Archimedes number correction
- **PDF Export**: Generate and export calculation reports
- **Real-time Updates**: Chart updates as parameters change

## Key Equations

Based on ASHRAE Standard RP-884:

### Acceptable Comfort Line
$$P_j = -m \cdot T_j + C$$

where:
- m, C determined from comfort regression table
- Depends on: Icl (clothing), M (metabolic rate), Vj (velocity)
- Tmr (mean radiant temperature)

### Jet Decay Model
$$\frac{V_j}{V_0} = \frac{1.464}{r}$$
$$\frac{T_A - T_j}{T_A - T_0} = \frac{4.539}{r}$$

where:
- r = (X₀/D₀) + 2.572
- Optional buoyancy correction using Archimedes number

## Dependencies

- **PyQt5** (≥5.15): GUI framework
- **NumPy** (≥1.20): Numerical computations
- **Matplotlib** (≥3.5): Chart visualization

## Notes

- All temperatures in °C
- All pressures in mmHg (except atmospheric in kPa)
- Humidity ratio in kg_water/kg_dry_air
- Flow rates in m³/s
- Metabolic rate in W/m²
- Clothing in clo units
