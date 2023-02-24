# About

KIAM Astrodynamics Toolbox is an astrodynamical toolbox written in Fortran and Python in Keldysh Institute of Applied Mathematics, Moscow, Russia.
The toolbox contains astronomical constants, physical parameters of celestial bodies, right-hand side of equations of motion in various models, methods for solving ODEs, routines for transformation variables, coordinate systems, and units.
The toolbox also contains a high-level class for trajectory models.

## Installation

`pip install kiam_astro` is the recommended way for installing the toolbox.

One can also download the toolbox using the Releases section on the right.

## Versions

Current version: 1.5.

The next version 2.0 will be released on Feb 27, 2023.

## The toolbox files

The main toolbox files are in the `kiam_astro` directory:

- `kiam.py ` is the main file containing the interface functions for accessing the fortran-compiled modules listed below.

- `FKIAMToolbox.cp39-win_amd64.pyd`, `FKIAMToolbox.cpython-39-darwin.so`, `FKIAMToolbox.cpython-39-x86_64-linux-gnu.so` are python modules containing compiled Fortran routines. It is assumed that users will not import or use these files directly.

- `trajectory.py ` is a high–level class for designing trajectories, it facilitates the propagation of trajectories, transformation of coordinate systems, variables, and units.

- `JPLEPH` is a file with ephemerides of celestial bodies, used by functions in fortran-compiled modules, developed by NASA JPL. This file contains ephemeris D430, which can also be downloaded independently from the site https://ssd.jpl.nasa.gov /.

- `dll files` are auxiliary files for the Fortran–Python interface from Intel(r) Visual Fortran Compiler, which can also be downloaded separately.

- `images` contains images used by the toolbox, e.g., the Earth's and Moon's surfaces.

The `examples` directory contains some typical examples of using the toolbox.

## System requirements

Windows, macOS, Ubuntu.

Python 3.9 only.

See `requirements.txt` for the dependencies.

On macOS and Ubuntu please install gfortran: `sudo apt install gfortran`.

## Information and terms of use

Reference on `kiam.py`: https://shmaxg.github.io/KIAMToolbox/html/kiam.html

Reference on `trajectory.py`: https://shmaxg.github.io/KIAMToolbox/html/Trajectory.html

See also the wiki page of the project on GitHub: https://github.com/shmaxg/KIAMToolbox/wiki

The library is being actively developed and improved, and some parts of it may change significantly in the next year or two.

The library is free for use (MIT License is used), the only request is to mention its use in acknowledgements in scientific papers, for example like this:

"Package KIAM Astrodynamics Toolbox by Maksim Shirobokov, MIT license, see https://github.com/shmaxg/KIAMToolbox for details"

In the future, a link to an article describing the toolbox will appear here, which can also be referenced.

For various questions, you can contact the main author of the library, Maksim Shirobokov: shirobokov@keldysh.ru.