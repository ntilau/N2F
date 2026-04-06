"""
Setup configuration for N2F Python package
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="pyn2f",
    version="1.0.0",
    author="Laurent Ntibarikure",
    author_email="",
    description="Near-field to far-field transformations with model order reduction",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ntilau/N2F",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "matplotlib>=3.4.0",
    ],
    project_urls={
        "Bug Tracker": "https://github.com/ntilau/N2F/issues",
        "Source Code": "https://github.com/ntilau/N2F",
    },
)
