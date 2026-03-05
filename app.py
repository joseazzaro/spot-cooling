# -*- coding: utf-8 -*-
"""
Spot Cooling Designer - ASHRAE Standard RP-884 Implementation
Main application entry point

MVC Architecture:
- Models: Psychrometric calculations and cooling design logic
- Views: PyQt5 GUI components (main window, charts)
- Controllers: Business logic and user interaction handling
"""

import sys
from PyQt5 import QtWidgets

from views.main_window import MainWindow


def main():
    """Application entry point"""
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
