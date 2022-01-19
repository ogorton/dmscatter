module kinds
    implicit none
    integer,parameter :: dp = kind(1.d0)
end module kinds

module types

    use kinds
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
        real(dp) :: j = 0.5d0! darkmatter spin
        real(dp) :: mass = 50.0d0
        real(dp) :: localdensity = 0.3d0
    end type particle
    
    type onebdenmat
        logical good
        real(dp), allocatable :: rho(:,:,:,:)
    end type onebdenmat

    type state
        integer Jx2 ! spin times 2
        integer Tx2 ! isospin times 2
    end type state

    type nucleus
        integer :: Z
        integer :: N
        integer :: A
        integer :: Mt
        real(dp) :: mass
        real(dp) :: Nt ! mass density of target nuclei
        type(state) :: groundstate

        type (onebdenmat) :: densitymats
        logical :: evenA
    end type nucleus

end module types

module main

    use types
    implicit none

    type(nucleus) :: nuc_target
    type(eftheory) :: eft
    type(particle) :: wimp

end module main    

module settings

    use kinds
    logical :: usemomentum = .false.
    logical :: useenergyfile = .false.
    
    logical :: fillcore = .true.
    logical :: printdens = .false.
    integer :: nlocalstates !,maxlocalstates
    integer :: nsporb   ! # of single-particle orbits

end module settings

module constants

    use kinds
    implicit none
    real(dp) :: pi = datan(1d0)*4d0!3.14159265358979323846264338327950288419716939937510
    real(dp) :: GeV = 1d0
    real(dp) :: kev = 10d0**(-6)
    real(dp) :: femtometer = 5.0677d0 ! /GeV
    real(dp) :: centimeter = 5.0677d0*10d0**13d0
    real(dp) :: kilometerpersecond = 3.33564d0 * 10d0**(-6d0)
    real(dp) :: kilogramday = 7.3634d0*10d0**55d0
    real(dp) :: mn = 0.938272d0
    real(dp) :: mv = 246.2d0

    real(dp) :: vearth  = 232.0d0!km/s
    REAL(dp) :: vescape = 12 * 232.d0! km/s
    REAL(dp) :: vscale  = 220.0d0! km/s 

    real(dp) :: ntscale = 1d0

end module constants
