!===============================================================================
module kinds
    implicit none
    integer,parameter :: doublep = kind(1.d0)
    integer :: keylength = 20

end module kinds


!===============================================================================
module response

    implicit none

    integer :: num_response_coef
    
    ! Nonrelativistic operator coefficients
    type coeff
        real(kind=8), dimension(:), allocatable :: c
    end type coeff
 
    ! 0 = protons, 1 = neutrons
    type(coeff) :: pncvec(0:1) 

    ! c0 = (cp+cn)/2   c1 = (cp-cn)/2
    type(coeff) :: cvec(0:1)

end module response

!===============================================================================
module dmparticles

    use kinds
    implicit none

    real(doublep) :: jchi ! darkmatter spin
    real(doublep) :: rhochi ! local dark matter density 

end module dmparticles


!===============================================================================
module targetinfo
   use kinds
   implicit none
   logical :: evenA
   integer nlocalstates !,maxlocalstates
   integer nsporb   ! # of single-particle orbits
   integer jt, tt   ! spin, isospin of transition operator
   integer Jiso, Tiso   ! ang, isospin of ground state times 2
   integer Mtiso ! iso spin projection
   logical ::  pndens
   type onebdenmat
   logical good
   real, allocatable :: rho(:,:,:,:)
   real, allocatable :: rhop(:,:,:),rhon(:,:,:)  ! added in version 9
   end type onebdenmat
   type (onebdenmat) :: densitymats

   real(doublep) :: Nt ! Number of target nuclei

end module targetinfo


!===============================================================================
module spspace
    implicit none
    ! Harmonic oscillator parameter
    real(kind=8) :: bfm
    !------------------SINGLE-PARTICLE STATES---------------------------   

    integer norb(2)                ! # of single-particle j-orbits
                                   ! 1 = p,  2 = n
    integer nsps(2)                ! # of single-particle m-states
                                   ! 1 = p, 2 = n
    integer ncore(2)               ! # of single-particle j-orbits in the core
                                   ! 1 = p,  2 = n
    integer ntotal(2)
    integer,allocatable  :: jorb(:)     ! 2 x J of orbit (note: p,n different)
    integer,allocatable  :: lorb(:)
    integer,allocatable  :: nrorb(:)    ! radial quantum number N of j-orbit
    integer,allocatable  :: torb(:)     ! 2 x Tz of orbit ( p/n = +/-1 )
    integer,allocatable  :: nprincipal(:)
    integer,allocatable  :: nodal(:)
    logical spinless
end module spspace

!===============================================================================
module masses

    implicit none

    ! Mass neutron and mass proton in GeV
    REAL(kind=8) :: mN
    ! 
    REAL(kind=8) :: mV
    ! Mass dark matter particle
    REAL(kind=8) :: mchi
    ! Mass isotope (number of nucleons)
    REAL(kind=8) :: mtarget
    ! Number of protons and neutron in target nucleus
    integer :: num_p, num_n
    ! Reduced mass of the target
    REAL(kind=8) :: muT
    ! Dark matter density
    real(kind=8) :: rhoDM

end module masses

!===============================================================================
module momenta
    use kinds
    implicit none

    real(doublep) :: q ! three-momentum transfer of the dm-nucleus scattering
    real(doublep) :: y ! y^2 is the rescaled q
    logical :: usemomentum = .false.
    logical :: useenergyfile = .false.

end module momenta

!===============================================================================
module velocities

    use kinds
    implicit none
    real(doublep) :: ve !km/s
    REAL(kind=8) :: vesc ! km/s
    REAL(kind=8) :: v0 ! km/s 

end module velocities

!===============================================================================
module constants
    implicit none
    real(kind=8) :: pi = 3.14159265358979323846264338327950288419716939937510
    real(kind=8) :: GeV = 1d0
    real(kind=8) :: kev = 10d0**(-6)
    real(kind=8) :: femtometer
end module constants

!===============================================================================
module quadrature
    use kinds
    implicit none

    integer :: quadrature_type
    integer :: lattice_points
    real(doublep) :: vdist_min
    real(doublep) :: vdist_max

end module quadrature

!===============================================================================
module keywords
    implicit none

    ! Do not change the order of these keywords!!! Their order determines which
    ! parameters are modified when the string is found in the control file.
    character(len=20), dimension(17) :: keyword_array = [character(20) :: &
        "coefnonrel","vearth","dmdens","quadtype","intpoints","gev","femtometer",&
        "wimpmass","vescape","ntarget","weakmscale","vscale","mnucleon","dmspin",&
        "usemomentum","maxwellv0","useenergyfile"]

end module keywords
