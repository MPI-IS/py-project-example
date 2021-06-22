from setuptools import setup


setup(
    name="nlse",
    version="0.1",
    description=(
        "The simulation and the plotting codes for simulating pulse "
        "propagation in systems governed by a Generalized Nonlinear "
        "Schrödinger Equation (GNLSE)"),
    author="Ivan Oreshnikov",
    author_email="oreshnikov.ivan@gmail.com",
    python_requires=">=3.5",
    install_requires=[
        "humanize==3.5.0",
        "matplotlib==3.3.0",
        "numpy==1.19.1",
        "scipy==1.5.2"
    ],
    scripts=[],
    # ^^^ in here one can list the scripts that should be available as
    # stand-alone commands when the package is installed.
    entry_points={
        "console_scripts": [
            "nlse-demo-collision=nlse.demos:soliton_collision"
        ]
    }
)
