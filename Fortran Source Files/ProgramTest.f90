program ProgramTest
    
    use MathAndProgConstants
    use BaseMeansToolbox
    use EquationsModule
    
    real(dp) :: dxdt(42), t, x(42)
    
    stm_required = .true.
    
    atm = .true.
    j2 = .true.
    srp = .true.
    sun = .true.
    mercury = .true.
    venus = .true.
    earth = .true.
    mars = .true.
    jupiter = .true.
    saturn = .true.
    uranus = .true.
    neptune = .true.
    moon = .true.
    cmplxmoon = .true.
    
    distunit = moondistunit
    velunit = moonvelunit
    timeunit = moontimeunit
    accunit = moonaccunit
    rsun = moonrsun
    rearth = moonrearth
    rmoon = moonrmoon
    
    JD_Zero = 2459599.5D0
    order = 50.0D0
    area = 1.0D0
    mass = 100.0D0
    
    t = 0.0D0
    x = [1.1D0, 0.0D0, 0.0D0, 0.0D0, 0.9D0, 0.0D0, reshape(eye(6),[36,1])]
    
    dxdt = knbp_rv_Moon(t,x)
    
end program ProgramTest