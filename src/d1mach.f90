    
!````
! d1mach(1) = b**(emin-1), the smallest positive magnitude.
! d1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
! d1mach(3) = b**(-t), the smallest relative spacing.
! d1mach(4) = b**(1-t), the largest relative spacing.
! d1mach(5) = log10(b)
!````

pure function d1mach(i) result(d)

    use iso_fortran_env, only: rp => real64

    implicit none 

    integer,intent(in) :: i
    real(rp) :: d

    real(rp),dimension(5),parameter :: constants = &
        [tiny(1.0_rp), &
         huge(1.0_rp), &
         real(radix(1.0_rp),rp)**(-digits(1.0_rp)), &
         epsilon(1.0_rp), &
         log10(real(radix(1.0_rp),rp))]

    d = constants(i)

end function d1mach
