import sys
sys.path.insert(0, 'c:/Users/Jose/Soft Projects/spot_cooling')
try:
    from app import MainWindow
    print('✓ App imports successfully - no syntax or import errors')
except Exception as e:
    print(f'✗ Error importing app: {e}')
    import traceback
    traceback.print_exc()
