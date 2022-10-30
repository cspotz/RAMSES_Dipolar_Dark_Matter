! Stahl 19/08/2020 : modification for dipolar dark matter

module poisson_commons
  use amr_commons
  use poisson_parameters

  real(dp),allocatable,dimension(:)  ::phi,phi_old       ! Potential
  real(dp),allocatable,dimension(:)  ::w_pol_p,w_pol_tograd,w_pol_old   ! Polarization potentials
  real(dp),allocatable,dimension(:)  ::pi_b_mod          ! module of 3-polarization
  real(dp),allocatable,dimension(:)  ::xi_b_mod          ! module of 3-xi_bot
  real(dp),allocatable,dimension(:)  ::rho               ! Density
  real(dp),allocatable,dimension(:)  ::rhod               ! Density only of dipole
  real(dp),allocatable,dimension(:)::divpi               ! to calculate the divergence of pi
  real(dp),allocatable,dimension(:,:)::f, f_old                 ! 3-force
  real(dp),allocatable,dimension(:,:)::f_int             ! 3-force due to dipole
  real(dp),allocatable,dimension(:,:)::ff2               ! grad sourcing xi_b
  real(dp),allocatable,dimension(:,:)::f2               ! to calculate xi . grad phi
  real(dp),allocatable,dimension(:,:)::f3               ! to calculate xi . grad phi
  real(dp),allocatable,dimension(:,:)::f4               ! to calculate xi . grad phi
  real(dp),allocatable,dimension(:,:)::pi_b,pi_b_old     ! 3-polarization
  real(dp),allocatable,dimension(:,:)::pi_dot           ! time variation of pi
  real(dp),allocatable,dimension(:,:)::xi_b              ! 3-xi_bot
  real(dp),allocatable,dimension(:,:)::xi_b_dot          ! derivative 3-xi_bot

  real(dp),allocatable,dimension(:)  ::rho_top   ! Density at last CIC level
  real(dp),allocatable,dimension(:)  ::rhod_top   ! Density at last CIC level

  ! Multigrid lookup table for amr -> mg index mapping
  integer, allocatable, dimension(:) :: lookup_mg   ! Lookup table

  ! Communicator arrays for multigrid levels
  type(communicator), allocatable, dimension(:,:) :: active_mg
  type(communicator), allocatable, dimension(:,:) :: emission_mg

  ! Minimum MG level
  integer :: levelmin_mg
  integer :: isdipol !dipolar DM to swich between newt and dipolar

  ! Multigrid safety switch
  logical, allocatable, dimension(:) :: safe_mode

  ! Multipole coefficients
  real(dp),dimension(1:ndim+1)::multipole

end module poisson_commons
