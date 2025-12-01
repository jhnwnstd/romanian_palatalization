#!/usr/bin/env python3
"""Package setup for Romanian Palatalization Analysis."""

from pathlib import Path

from setuptools import find_packages, setup

ROOT = Path(__file__).parent

# Long description
readme_path = ROOT / "README.md"
if readme_path.exists():
    long_description = readme_path.read_text(encoding="utf-8")
else:
    long_description = ""

# Install requirements
requirements_path = ROOT / "requirements.txt"
install_requires: list[str] = []
if requirements_path.exists():
    for line in requirements_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        install_requires.append(line)

setup(
    name="romanian-palatalization",
    version="1.0.0",
    description="Corpus-based analysis of palatalization patterns in Romanian morphophonology",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="John Winstead",
    author_email="j.hn.w.nst.d@gmail.com",
    url="https://github.com/jhnwnstd/romanian_palatalization",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.9",
    install_requires=install_requires,
    extras_require={
        "dev": [
            "black>=23.0.0",
            "flake8>=6.0.0",
            "mypy>=1.0.0",
            "isort>=5.12.0",
            "pytest>=7.0.0",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Text Processing :: Linguistic",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    keywords=[
        "romanian",
        "palatalization",
        "phonology",
        "morphophonology",
        "linguistics",
        "corpus-analysis",
    ],
    project_urls={
        "Bug Tracker": "https://github.com/jhnwnstd/romanian_palatalization/issues",
        "Source": "https://github.com/jhnwnstd/romanian_palatalization",
    },
)
