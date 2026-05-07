#==================================================#
# Install required Python packages for the project
#==================================================#

import subprocess
import sys
import importlib.util
import platform

def install_packages():
    # -------------------------------------------------------
    # List of third-party libraries to install
    # os, warnings, pathlib are built-in libraries, no need to list here
    # -------------------------------------------------------
    # Core packages (required)
    core_packages = [
        "numpy",
        "pandas",
        "geopandas",
        "osmnx"
    ]
    
    # Optional packages (may require additional system dependencies)
    optional_packages = [
        "pandana"  # Requires 'tables' which needs HDF5 libraries
    ]

    # Use Tsinghua mirror for faster downloads (recommended for large geospatial libraries)
    pypi_index_url = "https://pypi.tuna.tsinghua.edu.cn/simple"

    print(f"Current Python path called by Stata: {sys.executable}")
    print(f"Python version: {sys.version}")
    print(f"Platform: {platform.system()} {platform.machine()}")
    print("-" * 60)

    # Install core packages
    print("\n[Core Packages]")
    for package in core_packages:
        _install_package(package, pypi_index_url)

    # Install optional packages
    print("\n[Optional Packages]")
    for package in optional_packages:
        _install_package(package, pypi_index_url, optional=True)

    print("-" * 60)
    _print_summary()

def _install_package(package, pypi_index_url, optional=False):
    """Install a single package."""
    spec = importlib.util.find_spec(package)
    if spec is not None:
        print(f"[Installed] {package} detected, skipping.")
        return True
    
    print(f"[Installing] {package} ... This may take a few minutes...")
    try:
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", 
            package, 
            "-i", pypi_index_url
        ], stderr=subprocess.STDOUT)
        print(f"[Success] {package} installation complete.")
        return True
    except subprocess.CalledProcessError as e:
        if optional:
            print(f"[Skipped] {package} installation failed (optional package).")
            _print_optional_package_help(package)
        else:
            print(f"[Failed] {package} installation error.")
            print("Common causes: missing system dependencies or network timeout.")
        return False

def _print_optional_package_help(package):
    """Print help for optional package installation failure."""
    if package == "pandana":
        print("         pandana requires 'tables' which needs HDF5 libraries.")
        print("         To install manually:")
        print("         1. Install HDF5: brew install hdf5 c-blosc2")
        print("         2. Install tables: HDF5_DIR=/opt/homebrew pip3 install tables")
        print("         3. Install pandana: pip3 install pandana")
        print("         Or use conda: conda install -c conda-forge pandana")

def _print_summary():
    """Print installation summary."""
    print("\nInstallation Summary:")
    print("-" * 60)
    
    packages_status = {
        "numpy": importlib.util.find_spec("numpy") is not None,
        "pandas": importlib.util.find_spec("pandas") is not None,
        "geopandas": importlib.util.find_spec("geopandas") is not None,
        "osmnx": importlib.util.find_spec("osmnx") is not None,
        "pandana": importlib.util.find_spec("pandana") is not None,
    }
    
    for pkg, installed in packages_status.items():
        status = "✓ Installed" if installed else "✗ Not installed"
        print(f"  {pkg}: {status}")
    
    all_core_installed = all(packages_status[p] for p in ["numpy", "pandas", "geopandas", "osmnx"])
    
    if all_core_installed:
        print("\n[OK] All core packages are installed successfully!")
    else:
        print("\n[Warning] Some core packages failed to install.")
        print("Please ensure your system has basic geospatial libraries installed.")
        print("On macOS, recommend using: brew install gdal")

if __name__ == "__main__":
    install_packages()