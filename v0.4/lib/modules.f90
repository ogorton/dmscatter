!===============================================================================
module response

    implicit none
    
    ! Nonrelativistic operator coefficients
    type coeff
        real(kind=8), dimension(:), allocatable :: c
    end type coeff
 
    ! 0 = protons, 1 = neutrons
    type(coeff) :: cvec(0:1) 

end module response


!===============================================================================
module stateinfo
   implicit none
   logical :: evenA
   integer nlocalstates !,maxlocalstates
   integer nsporb   ! # of single-particle orbits
   integer jt, tt   ! spin, isospin of transition operator
   integer Jiso, Tiso   ! ang, isospin of ground state times 2
   logical ::  pndens
   type onebdenmat
   logical good
   real, allocatable :: rho(:,:,:,:)
   real, allocatable :: rhop(:,:,:),rhon(:,:,:)  ! added in version 9
   end type onebdenmat
   type (onebdenmat) :: densitymats

end module stateinfo


!===============================================================================
module spspace
    implicit none
    !------------------SINGLE-PARTICLE STATES---------------------------   

    integer norb(2)                ! # of single-particle j-orbits
                                   ! 1 = p,  2 = n
    integer nsps(2)                ! # of single-particle m-states
                                   ! 1 = p, 2 = n
    integer ncore(2)               ! # of single-particle j-orbits in the core
                                   ! 1 = p,  2 = n
    integer ntotal(2)
    integer,allocatable  :: jorb(:)        ! 2 x J of orbit (note: p,n different)
    integer,allocatable  :: lorb(:)
    integer,allocatable  :: nrorb(:)       ! radial quantum number N of j-orbit
    integer,allocatable  :: torb(:)        ! 2 x Tz of orbit ( p/n = +/-1 )
    integer,allocatable  :: nprincipal(:)
    integer,allocatable  :: nodal(:)
    logical spinless
end module spspace

!===============================================================================
module masses

    implicit none

    ! Mass neutron and mass proton in GeV
    REAL(kind=8) :: mN = 0.938272
    ! 
    REAL(kind=8) :: mV = 246.2
    ! Mass dark matter particle
    REAL(kind=8) :: mchi = 1.0
    ! Mass isotope (number of nucleons)
    REAL(kind=8) :: Miso = 1.0
    ! Reduced mass of the target
    REAL(kind=8) :: muT = 1.0
    ! Harmonic oscillator parameter
    real(kind=8) :: bfm

end module masses
