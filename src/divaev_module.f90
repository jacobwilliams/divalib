!*************************************************************************
!>
!  Holds variables that depend on the environment.

    module divaev_module

    implicit none

    public

    double precision :: eeps10  !! 10 * (machine epsilon).
    double precision :: eeps16  !! 16 * (machine epsilon).
    double precision :: erov10  !! 10 / (largest floating point number).
    double precision :: eeps2   !! 2 * (machine epsilon).
    double precision :: eept75  !! (machine epsilon) ** (.75)
    double precision :: eovep2  !! EEPS2 * (largest floating point number).
    double precision :: ovtm75  !! (largest floating point number) ** (-.75)
    double precision :: ovd10   !! (largest floating point number) / 10.

    double precision :: evc(8) !! Array used for output of variables EEPS2 to EROV10 in
                               !! common block DIVAEV.

    common / divaev / eeps2, eept75, eovep2, ovtm75, ovd10, &
                      eeps10, eeps16, erov10

    equivalence (evc(1), eeps2)

    end module divaev_module
!*************************************************************************
