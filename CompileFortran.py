import os

f2py_exec_command = 'python -m numpy.f2py '
module_name = '-m FKIAMToolbox '
f90compiler = '"C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\latest\\windows\\bin\\intel64\\ifort.exe" '
compiler_libs = '-L"C:\\Program Files (x86)\\Intel\\oneAPI\\compiler\\latest\\windows\\compiler\\lib\\intel64" '
working_folder = '-I:"G:\\Мой диск\\Python\\KIAMToolbox" '

file_1 = '".\\Fortran Source Files\\MathAndProgConstants.f90" '
file_2 = '".\\Fortran Source Files\\BaseMeansToolbox.f90" '
file_3 = '".\\Fortran Source Files\\LinearAlgebraLowLevel.f90" '
file_4 = '".\\Fortran Source Files\\LinearAlgebraInterfaces.f90" '
file_5 = '".\\Fortran Source Files\\OdeToolbox.f90" '
file_6 = '".\\Fortran Source Files\\Ephemeris.f90" '
file_7 = '".\\Fortran Source Files\\Translations.f90" '
file_8 = '".\\Fortran Source Files\\GravityMoonCoefficients50.f90" '
file_9 = '".\\Fortran Source Files\\EquationsModule.f90" '
file_10 = '".\\Fortran Source Files\\PropagationModule.f90" '
file_11 = '".\\Fortran Source Files\\ConstantsAndUnits.f90" '

skip_1 = '"norm" "norminf" "zeros" "ones" "eye" "cross" "diff" "concol" "conrow" "colonoperator" "linspace" "double2logical" '
skip_2 = '"dgemm" "xerbla" "lsame" "dgetrf" "dtrsm" "dgetf2" "dlaswp" "dger" "dgetrs" '
skip_3 = '"dscal" "dswap" "dlamch" "dlamc3" "ilaenv" "idamax" "ieeeck" "iparmq" "iparam2stage" '
skip_4 = '"linsolve" "dotmm" "dotvm" "dotmv" "doto" '
skip_5 = '"ode4" "ode8" "ode45" "ode113" "ode87" '
skip_6 = '"fsizer1" "fsizer2" "fsizer3" "pleph" "interp" "split" "state" "const_eph" '
skip_7 = '"getq" "getr" '
skip_8 = '"cmplx_moon" "kgravforce" "ksrpforce" "kintensity" "kgravitymoongradient50" "kearthatm" "kearthj2" '

os.system(f2py_exec_command +
          '-c ' + module_name + '--f90exec=' + f90compiler + compiler_libs + working_folder +
          file_1 + file_2 + file_3 + file_4 + file_5 + file_6 + file_7 + file_8 + file_9 + file_10 + file_11 +
          'skip: ' + skip_1 + skip_2 + skip_3 + skip_4 + skip_5 + skip_6 + skip_7 + skip_8 + ':')
