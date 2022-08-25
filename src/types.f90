module types

    implicit none

    ! Nonrelativistic operator coefficients
    type coeff
        real(kind=8), dimension(:), allocatable :: c
    end type coeff

    type eftheory
        type(coeff) :: isoc(0:1) !c0 = (cp+cn)/2   c1 = (cp-cn)/2
        type(coeff) :: xpnc(0:1) !0 = protons, 1 = neutrons
    end type eftheory

    ! Dark matter particles
    type particle
        real(kind=8) :: j = 0.5d0! darkmatter spin
        real(kind=8) :: mass = 50.0d0
        real(kind=8) :: localdensity = 0.3d0
    end type particle
    
    type onebdenmat
        logical good
        real(kind=8), allocatable :: rho(:,:,:,:)
    end type onebdenmat

    type state
        integer Jx2 ! spin times 2
        integer Tx2 ! isospin times 2
    end type state

    type nucleus
        integer :: Z
        integer :: N
        integer :: Zval, Nval
        integer :: A
        integer :: Mt
        real(kind=8) :: mass
        real(kind=8) :: Nt ! mass density of target nuclei
        type(state) :: groundstate

        type (onebdenmat) :: densitymats
        logical :: evenA
    end type nucleus

end module types
