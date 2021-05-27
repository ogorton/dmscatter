!===============================================================================
module kinds
    implicit none
    integer,parameter :: doublep = kind(1.d0)
end module kinds

!===============================================================================
module parameters

    use kinds

    implicit none

    ! Nonrelativistic operator coefficients
    integer :: num_response_coef = 17
    type coeff
        real(kind=8), dimension(:), allocatable :: c
    end type coeff
 
    type eftheory
        type(coeff) :: isoc(0:1) !c0 = (cp+cn)/2   c1 = (cp-cn)/2
        type(coeff) :: xpnc(0:1) !0 = protons, 1 = neutrons
    end type eftheory
    
    ! Dark matter particles
    type particle
        real(doublep) :: j = 0.5d0! darkmatter spin
        real(doublep) :: mass = 50.0d0
        real(doublep) :: localdensity = 0.3d0
    end type particle

   ! Target nucleus
   real(doublep) :: ntscale = 1d0
   logical :: fillcore = .true.
   logical :: printdens = .false.
   logical :: evenA
   integer nlocalstates !,maxlocalstates
   integer nsporb   ! # of single-particle orbits
   integer jt, tt   ! spin, isospin of transition operator
   integer :: maxeventJt = 5
   integer :: maxoddJt = 6
   logical ::  pndens ! true: proton-neutron format density matrices.
                      ! false: isospin-conserving format density matrices

   type onebdenmat
       logical good
       real(doublep), allocatable :: rho(:,:,:,:)
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
       real(doublep) :: mass
       real(doublep) :: Nt ! mass density of target nuclei
       real(doublep) ::kgday
       type(state) :: groundstate

       type (onebdenmat) :: densitymats
       logical :: evenA
   end type nucleus

end module parameters
!===============================================================================
module main

    use parameters
    implicit none

    type(nucleus) :: nuc_target
    type(eftheory) :: eft
    type(particle) :: wimp

end module

!===============================================================================
module momenta
    use kinds
    implicit none

    logical :: usemomentum = .false.
    logical :: useenergyfile = .false.

end module momenta


!===============================================================================
module constants

    use kinds
    implicit none
    real(doublep) :: pi = 3.14159265358979323846264338327950288419716939937510
    real(doublep) :: GeV = 1d0
    real(doublep) :: kev = 10d0**(-6)
    real(doublep) :: femtometer = 5.0677d0 ! /GeV
    real(doublep) :: centimeter = 5.0677d0*10d0**13d0
    real(doublep) :: kilometerpersecond = 3.33564d0 * 10d0**(-6d0)
    real(doublep) :: kilogramday = 7.3634d0*10d0**55d0
    real(doublep) :: mn = 0.938272d0
    real(doublep) :: mv = 246.2d0

    real(doublep) :: vearth  = 232.0d0!km/s
    REAL(doublep) :: vescape = 12 * 232.d0! km/s
    REAL(doublep) :: vscale  = 220.0d0! km/s 

end module constants

!===============================================================================
module quadrature
    use kinds
    implicit none

    integer :: quadrature_type = 1
    real(doublep) :: quadrature_relerr = 1e-6
    integer :: gaussorder = 12

end module quadrature

!===============================================================================
module keywords
    implicit none

    type keywordpair
        character(len=20) :: key
        real(kind=8) :: val
    end type keywordpair

    integer, parameter :: maxkeywords = 100
    integer :: numkeywords = 0
    type(keywordpair) :: keywordpairs(maxkeywords)

    contains

        subroutine addKeywordpair(key, val)
            implicit none
            character(len=20) :: key
            real(kind=8) :: val
            numkeywords = numkeywords + 1
            keywordpairs(numkeywords)%key = key
            keywordpairs(numkeywords)%val = val
        end subroutine addKeywordpair

end module keywords
