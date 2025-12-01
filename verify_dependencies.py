#!/usr/bin/env python3
"""
Small helper script to check that the main dependencies are installed.
"""

import sys
from importlib import import_module

REQUIRED_PACKAGES = {
    "pandas": "pandas",
    "requests": "requests",
    "beautifulsoup4": "bs4",
    "tenacity": "tenacity",
    "urllib3": "urllib3",
}

OPTIONAL_PACKAGES = {
    "black": "black",
    "flake8": "flake8",
    "flake8-bugbear": "bugbear",
    "flake8-comprehensions": "flake8_comprehensions",
    "mypy": "mypy",
    "isort": "isort",
}


def check_package(label: str, import_name: str) -> bool:
    """Try importing a package and print a short status line."""
    try:
        module = import_module(import_name)
    except ImportError:
        print(f"[missing] {label}")
        return False

    version = getattr(module, "__version__", None)
    if version:
        print(f"[ok]      {label} (version {version})")
    else:
        print(f"[ok]      {label}")
    return True


def check_python_version() -> bool:
    """Require Python >= 3.9."""
    v = sys.version_info
    print(f"Python: {v.major}.{v.minor}.{v.micro}")
    if v < (3, 9):
        print("Python 3.9 or newer is required.")
        return False
    return True


def main() -> int:
    print("Romanian palatalization – dependency check")
    print("-" * 50)

    python_ok = check_python_version()
    print()

    print("Required packages:")
    required_ok = True
    for pkg, import_name in REQUIRED_PACKAGES.items():
        if not check_package(pkg, import_name):
            required_ok = False

    print("\nOptional (dev) packages:")
    optional_count = 0
    for pkg, import_name in OPTIONAL_PACKAGES.items():
        if check_package(pkg, import_name):
            optional_count += 1

    print("\n" + "-" * 50)
    if python_ok and required_ok:
        print("All required dependencies are available.")
        print(f"Optional dev tools: {optional_count}/{len(OPTIONAL_PACKAGES)}")
        print("\nYou can run the pipeline with:")
        print("  ./run_pipeline.sh")
        return 0

    print("Some required dependencies are missing.")
    print("Install them with:")
    print("  pip install .")
    print("Or for editable install with dev tools:")
    print("  pip install -e .[dev]")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
