subroutine interpol_phi(ind_cell,phi_int,ncell,ilevel,icount)
  use amr_commons
  use poisson_commons, only:phi,phi_old
  implicit none
  integer::ncell,ilevel,icount
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is on -
  ! time (linear extrapolation of the change in phi during the last coarse step
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1d0/4d0**ndim
  bb = 3*aa
  cc = 9*aa
  dd = 27*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     tfrac=1d0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0
  end if

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(i,ind)=0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*(phi(ind_cell(i))+(phi(ind_cell(i))-phi_old(ind_cell(i)))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(ind_cell(i))+phi(ind_cell(i)))
           else
              add=coeff*(phi(indice)+(phi(indice)-phi_old(indice))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(indice)+phi(indice))
           endif
           phi_int(i,ind)=phi_int(i,ind)+add
        end do
     end do
  end do

 end subroutine interpol_phi

!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine save_phi_old(ilevel)
  use amr_commons
  use poisson_commons, only:phi,phi_old
  implicit none
  integer ilevel

  !save the old potential for time extrapolation in case of subcycling

  integer::i,ncache,ind,igrid,iskip,istart,ibound
  integer,allocatable,dimension(:)::ind_grid

  do ibound=1,nboundary+ncpu
     if(ibound<=ncpu)then
        ncache=numbl(ibound,ilevel)
        istart=headl(ibound,ilevel)
     else
        ncache=numbb(ibound-ncpu,ilevel)
        istart=headb(ibound-ncpu,ilevel)
     end if
     if(ncache>0)then
        allocate(ind_grid(1:ncache))
        ! Loop over level grids
        igrid=istart
        do i=1,ncache
           ind_grid(i)=igrid
           igrid=next(igrid)
        end do
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           ! save phi
           do i=1,ncache
              phi_old(ind_grid(i)+iskip)=phi(ind_grid(i)+iskip)
           end do
        end do
        deallocate(ind_grid)
     end if
  end do

end subroutine save_phi_old


subroutine interpol_w(ind_cell,phi_int,ncell,ilevel,icount)
  use amr_commons
  use poisson_commons, only:w_pol_p,w_pol_old
  implicit none
  integer::ncell,ilevel,icount
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is on -
  ! time (linear extrapolation of the change in phi during the last coarse step
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1d0/4d0**ndim
  bb = 3*aa
  cc = 9*aa
  dd = 27*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     tfrac=1d0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0
  end if

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(i,ind)=0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*(w_pol_p(ind_cell(i))+(w_pol_p(ind_cell(i))-w_pol_old(ind_cell(i)))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(ind_cell(i))+phi(ind_cell(i)))
           else
              add=coeff*(w_pol_p(indice)+(w_pol_p(indice)-w_pol_old(indice))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(indice)+phi(indice))
           endif
           phi_int(i,ind)=phi_int(i,ind)+add
        end do
     end do
  end do

 end subroutine interpol_w
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine save_w_old(ilevel) !dipolar DM
  use amr_commons
  use poisson_commons, only:w_pol_p,w_pol_old
  implicit none
  integer ilevel

  !save the old potential w for time extrapolation in case of subcycling

  integer::i,ncache,ind,igrid,iskip,istart,ibound
  integer,allocatable,dimension(:)::ind_grid

  do ibound=1,nboundary+ncpu
     if(ibound<=ncpu)then
        ncache=numbl(ibound,ilevel)
        istart=headl(ibound,ilevel)
     else
        ncache=numbb(ibound-ncpu,ilevel)
        istart=headb(ibound-ncpu,ilevel)
     end if
     if(ncache>0)then
        allocate(ind_grid(1:ncache))
        ! Loop over level grids
        igrid=istart
        do i=1,ncache
           ind_grid(i)=igrid
           igrid=next(igrid)
        end do
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           ! save phi
           do i=1,ncache
              w_pol_old(ind_grid(i)+iskip)=w_pol_p(ind_grid(i)+iskip)
           end do
        end do
        deallocate(ind_grid)
     end if
  end do

end subroutine save_w_old

!13/01/21 dipolarDM, pour div pi
subroutine interpol_pi_bx(ind_cell,phi_int,ncell,ilevel,icount)
  use amr_commons
   use poisson_commons, only:pi_b,pi_b_old
  implicit none
  integer::ncell,ilevel,icount
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is on -
  ! time (linear extrapolation of the change in phi during the last coarse step
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1d0/4d0**ndim
  bb = 3*aa
  cc = 9*aa
  dd = 27*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     tfrac=1d0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0
  end if

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(i,ind)=0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*(pi_b(ind_cell(i),1)+(pi_b(ind_cell(i),1)-pi_b_old(ind_cell(i),1))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(ind_cell(i))+phi(ind_cell(i)))
           else
              add=coeff*(pi_b(indice,1)+(pi_b(indice,1)-pi_b_old(indice,1))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(indice)+phi(indice))
           endif
           phi_int(i,ind)=phi_int(i,ind)+add
        end do
     end do
  end do

 end subroutine interpol_pi_bx
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine save_pi_old(ilevel) !dipolar DM
  use amr_commons
  use poisson_commons, only:pi_b,pi_b_old
  implicit none
  integer ilevel

  !save the old potential pi_b for time extrapolation in case of subcycling

  integer::i,ncache,ind,igrid,iskip,istart,ibound,idim
  integer,allocatable,dimension(:)::ind_grid

  do ibound=1,nboundary+ncpu
     if(ibound<=ncpu)then
        ncache=numbl(ibound,ilevel)
        istart=headl(ibound,ilevel)
     else
        ncache=numbb(ibound-ncpu,ilevel)
        istart=headb(ibound-ncpu,ilevel)
     end if
     if(ncache>0)then
        allocate(ind_grid(1:ncache))
        ! Loop over level grids
        igrid=istart
        do i=1,ncache
           ind_grid(i)=igrid
           igrid=next(igrid)
        end do
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           ! Loop over dimensions
           do idim=1,ndim
              ! save phi
              do i=1,ncache
                  pi_b_old(ind_grid(i)+iskip,idim)=pi_b(ind_grid(i)+iskip,idim)
              end do
           end do
        end do
        deallocate(ind_grid)
     end if
  end do

end subroutine save_pi_old



subroutine interpol_pi_bz(ind_cell,phi_int,ncell,ilevel,icount)
  use amr_commons
  use poisson_commons, only:pi_b,pi_b_old
  implicit none
  integer::ncell,ilevel,icount
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is on -
  ! time (linear extrapolation of the change in phi during the last coarse step
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1d0/4d0**ndim
  bb = 3*aa
  cc = 9*aa
  dd = 27*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     tfrac=1d0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0
  end if

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(i,ind)=0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*(pi_b(ind_cell(i),3)+(pi_b(ind_cell(i),3)-pi_b_old(ind_cell(i),3))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(ind_cell(i))+phi(ind_cell(i)))
           else
              add=coeff*(pi_b(indice,3)+(pi_b(indice,3)-pi_b_old(indice,3))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(indice)+phi(indice))
           endif
           phi_int(i,ind)=phi_int(i,ind)+add
        end do
     end do
  end do

 end subroutine interpol_pi_bz

subroutine interpol_pi_by(ind_cell,phi_int,ncell,ilevel,icount)
  use amr_commons
  use poisson_commons, only:pi_b,pi_b_old
  implicit none
  integer::ncell,ilevel,icount
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is on -
  ! time (linear extrapolation of the change in phi during the last coarse step
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1d0/4d0**ndim
  bb = 3*aa
  cc = 9*aa
  dd = 27*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     tfrac=1d0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0
  end if

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(i,ind)=0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*(pi_b(ind_cell(i),2)+(pi_b(ind_cell(i),2)-pi_b_old(ind_cell(i),2))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(ind_cell(i))+phi(ind_cell(i)))
           else
              add=coeff*(pi_b(indice,2)+(pi_b(indice,2)-pi_b_old(indice,2))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(indice)+phi(indice))
           endif
           phi_int(i,ind)=phi_int(i,ind)+add
        end do
     end do
  end do

 end subroutine interpol_pi_by


!25/01/21 dipolarDM, pour xi dot grad F
subroutine interpol_Fx(ind_cell,phi_int,ncell,ilevel,icount)
  use amr_commons
   use poisson_commons
  implicit none
  integer::ncell,ilevel,icount
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is on -
  ! time (linear extrapolation of the change in phi during the last coarse step
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1d0/4d0**ndim
  bb = 3*aa
  cc = 9*aa
  dd = 27*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     tfrac=1d0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0
  end if

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(i,ind)=0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*(f(ind_cell(i),1)+(f(ind_cell(i),1)-f_old(ind_cell(i),1))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(ind_cell(i))+phi(ind_cell(i)))
           else
              add=coeff*(f(indice,1)+(f(indice,1)-f_old(indice,1))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(indice)+phi(indice))
           endif
           phi_int(i,ind)=phi_int(i,ind)+add
        end do
     end do
  end do

 end subroutine interpol_Fx
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine save_f_old(ilevel) !dipolar DM
  use amr_commons
  use poisson_commons, only:f,f_old
  implicit none
  integer ilevel

  !save the old force f for time extrapolation in case of subcycling

  integer::i,ncache,ind,igrid,iskip,istart,ibound,idim
  integer,allocatable,dimension(:)::ind_grid

  do ibound=1,nboundary+ncpu
     if(ibound<=ncpu)then
        ncache=numbl(ibound,ilevel)
        istart=headl(ibound,ilevel)
     else
        ncache=numbb(ibound-ncpu,ilevel)
        istart=headb(ibound-ncpu,ilevel)
     end if
     if(ncache>0)then
        allocate(ind_grid(1:ncache))
        ! Loop over level grids
        igrid=istart
        do i=1,ncache
           ind_grid(i)=igrid
           igrid=next(igrid)
        end do
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           ! Loop over dimensions
           do idim=1,ndim
              ! save phi
              do i=1,ncache
                  f_old(ind_grid(i)+iskip,idim)=f(ind_grid(i)+iskip,idim)
              end do
           end do
        end do
        deallocate(ind_grid)
     end if
  end do

end subroutine save_f_old



subroutine interpol_Fz(ind_cell,phi_int,ncell,ilevel,icount)
  use amr_commons
  use poisson_commons, only:f,f_old
  implicit none
  integer::ncell,ilevel,icount
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is on -
  ! time (linear extrapolation of the change in phi during the last coarse step
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1d0/4d0**ndim
  bb = 3*aa
  cc = 9*aa
  dd = 27*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     tfrac=1d0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0
  end if

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(i,ind)=0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*(f(ind_cell(i),3)+(f(ind_cell(i),3)-f_old(ind_cell(i),3))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(ind_cell(i))+phi(ind_cell(i)))
           else
              add=coeff*(f(indice,3)+(f(indice,3)-f_old(indice,3))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(indice)+phi(indice))
           endif
           phi_int(i,ind)=phi_int(i,ind)+add
        end do
     end do
  end do

 end subroutine interpol_Fz

subroutine interpol_Fy(ind_cell,phi_int,ncell,ilevel,icount)
  use amr_commons
  use poisson_commons, only:f,f_old
  implicit none
  integer::ncell,ilevel,icount
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is on -
  ! time (linear extrapolation of the change in phi during the last coarse step
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1d0/4d0**ndim
  bb = 3*aa
  cc = 9*aa
  dd = 27*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     tfrac=1d0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0
  end if

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(i,ind)=0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*(f(ind_cell(i),2)+(f(ind_cell(i),2)-f_old(ind_cell(i),2))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(ind_cell(i))+phi(ind_cell(i)))
           else
              add=coeff*(f(indice,2)+(f(indice,2)-f_old(indice,2))*tfrac)
              ! add=coeff*(-3d0/8d0*dx**2*boxlen*rho(indice)+phi(indice))
           endif
           phi_int(i,ind)=phi_int(i,ind)+add
        end do
     end do
  end do

 end subroutine interpol_Fy
