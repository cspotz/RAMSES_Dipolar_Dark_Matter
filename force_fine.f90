! Stahl 19/08/2020 : modification for dipolar dark matter
! Here all the relevant quantities for the dipolar DM are calculated
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine force_fine(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  use mpi_mod
  use constants, only : twopi

  use hydro_commons, ONLY: smallr

  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,icount
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::nx_loc,idim
  real(dp)::dx,dx_loc,scale,fact,fourpi
  real(kind=8)::rho_loc,rho_all,epot_loc,epot_all,modef,small_r !Stahl 17/02/21 need du module de f, small_r when rho=0
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx,ff

!small_r=1E-35

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  !-------------------------------------
  ! Compute analytical gravity force
  !-------------------------------------
  if(gravity_type>0)then

     ! Loop over myid grids by vector sweeps
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim

           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Gather cell centre positions
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do
           ! Rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
              end do
           end do

           ! Call analytical gravity routine
           call gravana(xx,ff,dx_loc,ngrid)


           ! Scatter variables
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell(i),idim)=ff(i,idim)
              end do
           end do

        end do
        ! End loop over cells

     end do
     ! End loop over grids

     ! Update boundaries
     do idim=1,ndim
        call make_virtual_fine_dp(f(1,idim),ilevel)
     end do
     if(simple_boundary)call make_boundary_force(ilevel)

  else
  !------------------------------
  ! Compute gradient of potential
  !------------------------------
     ! Update physical boundaries
     call make_boundary_phi(ilevel)

     ! Loop over myid grids by vector sweeps
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do



        ! Compute gradient of potential
        call gradient_phi(ind_grid,ngrid,ilevel,icount)

     end do
     ! End loop over grids
#if NDIM==3
     if (sink)then
        call f_gas_sink(ilevel)
     end if
#endif
     ! Update boundaries
     do idim=1,ndim
        call make_virtual_fine_dp(f(1,idim),ilevel)
     end do
     if(simple_boundary)call make_boundary_force(ilevel)


  
 endif

  !----------------------------------------------
  ! Compute gravity potential and maximum density
  !----------------------------------------------
  rho_loc =0; rho_all =0
  epot_loc=0; epot_all=0
  fourpi=2*twopi
  if(cosmo)fourpi=1.5D0*omega_m*aexp
  fact=-dx_loc**ndim/fourpi/2.0D0


  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Loop over dimensions
        do idim=1,ndim
           do i=1,ngrid
              if(son(ind_cell(i))==0)then
                 epot_loc=epot_loc+fact*f(ind_cell(i),idim)**2
	      end if
           end do
        end do
        ! End loop over dimensions
               do i=1,ngrid
           rho_loc=MAX(rho_loc,dble(abs(rho(ind_cell(i)))))        
	   		
        end do

      
     end do
     ! End loop over cells
  end do
  ! End loop over grids

!something to be done here for for dipolar DM ?
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(epot_loc,epot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(rho_loc ,rho_all ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     epot_loc=epot_all
     rho_loc =rho_all
#endif
     epot_tot=epot_tot+epot_loc
     rho_max(ilevel)=rho_loc

111 format('   Entering force_fine for level ',I2)

end subroutine force_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine dipolar_fine(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  use mpi_mod
  use constants, only : twopi

  use hydro_commons, ONLY: smallr

  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,icount
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::nx_loc,idim
  real(dp)::dx,dx_loc,scale,fact,fourpi
  real(kind=8)::rho_loc,rho_all,epot_loc,epot_all,modef,small_r !Stahl 17/02/21 need du module de f, small_r when rho=0
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx,ff

!small_r=1E-35

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  !-------------------------------------
  ! Compute analytical gravity force
  !-------------------------------------
  if(gravity_type>0)then

   
     if(nstep.ne.0)then

  !-----------------------------------------------
  ! Dipolar : Compute gradient of Fx to update xi_b
  !-----------------------------------------------

 ! Loop over myid grids by vector sweeps
    	 ncache=active(ilevel)%ngrid
     	do igrid=1,ncache,nvector
       	 ngrid=MIN(nvector,ncache-igrid+1)
       	 do i=1,ngrid
           	ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        	end do
        	call grad_Fx(ind_grid,ngrid,ilevel,icount)

     	end do
     ! End loop over grids
     do idim=1,ndim
        call make_virtual_fine_dp(f2(1,idim),ilevel)
     end do


  !------------------------------------------------
  ! Dipolar : Compute gradient of Fy to update xi_b
  !-----------------------------------------------

 ! Loop over myid grids by vector sweeps
     	ncache=active(ilevel)%ngrid
     	do igrid=1,ncache,nvector
        	ngrid=MIN(nvector,ncache-igrid+1)
        	do i=1,ngrid
           	ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        	end do
        	call grad_Fy(ind_grid,ngrid,ilevel,icount) !dipole DM need grad F to update xi_b
     	end do
     ! End loop over grids
    do idim=1,ndim
        call make_virtual_fine_dp(f3(1,idim),ilevel)
     end do



  !------------------------------------------------
  ! Dipolar : Compute gradient of Fz to update xi_b
  !-----------------------------------------------

 ! Loop over myid grids by vector sweeps
     	ncache=active(ilevel)%ngrid
     	do igrid=1,ncache,nvector
        	ngrid=MIN(nvector,ncache-igrid+1)
        	do i=1,ngrid
           	ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        	end do
	       	call grad_Fz(ind_grid,ngrid,ilevel,icount) !dipole DM need grad F to update xi_b

     	end do
     ! End loop over grids
    do idim=1,ndim
        call make_virtual_fine_dp(f4(1,idim),ilevel)
     end do


  !---------------------------------------------
  ! Dipolar : Compute gradient of W_pot_tot_grad 
  !---------------------------------------------

 ! Loop over myid grids by vector sweeps
     	ncache=active(ilevel)%ngrid
     	do igrid=1,ncache,nvector
        	ngrid=MIN(nvector,ncache-igrid+1)
        	do i=1,ngrid
           	ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        	end do
        	call gradient_w(ind_grid,ngrid,ilevel,icount)  

     	end do
     	! End loop over grids

   do idim=1,ndim
        call make_virtual_fine_dp(ff2(1,idim),ilevel)
     end do

     endif


 endif

  !----------------------------------------------
  ! Compute gravity potential and maximum density
  !----------------------------------------------
  fourpi=2*twopi
  if(cosmo)fourpi=1.5D0*omega_m*aexp
  fact=-dx_loc**ndim/fourpi/2.0D0



  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Dipolar DM new Loop over dimensions
        do idim=1,ndim
           do i=1,ngrid
      			if(rhod(ind_cell(i)).ne.0)then
     			xi_b_dot(ind_cell(i),idim)=xi_b_dot(ind_cell(i),idim)+f_int(ind_cell(i),idim)*dtnew(ilevel)-ff2(ind_cell(i),idim)*dtnew(ilevel)/(rhod(ind_cell(i))) !minus because every grad in RAMSES comes with a minus
     			endif
		

!!!!!!Add the non-sperical terms
			if(idim==1)then !not very elegant way to add xi dot grad phi
				xi_b_dot(ind_cell(i),idim)=xi_b_dot(ind_cell(i),idim)+dtnew(ilevel)*(xi_b(ind_cell(i),1)*f2(ind_cell(i),1)+xi_b(ind_cell(i),2)*f2(ind_cell(i),2)+xi_b(ind_cell(i),3)*f2(ind_cell(i),3))
			endif
			if(idim==2)then 
				xi_b_dot(ind_cell(i),idim)=xi_b_dot(ind_cell(i),idim)+dtnew(ilevel)*(xi_b(ind_cell(i),1)*f3(ind_cell(i),1)+xi_b(ind_cell(i),2)*f3(ind_cell(i),2)+xi_b(ind_cell(i),3)*f3(ind_cell(i),3))
			endif
			if(idim==3)then 
				xi_b_dot(ind_cell(i),idim)=xi_b_dot(ind_cell(i),idim)+dtnew(ilevel)*(xi_b(ind_cell(i),1)*f4(ind_cell(i),1)+xi_b(ind_cell(i),2)*f4(ind_cell(i),2)+xi_b(ind_cell(i),3)*f4(ind_cell(i),3))
			endif
!!!!!! Fin Add the non-sperical terms

		xi_b(ind_cell(i),idim)=xi_b(ind_cell(i),idim)+xi_b_dot(ind_cell(i),idim)*dtnew(ilevel) 
		pi_b(ind_cell(i),idim)= pi_b(ind_cell(i),idim)+xi_b_dot(ind_cell(i),idim)*dtnew(ilevel)*(rhod(ind_cell(i))) !dipolar DM, 2)update polarization field
           end do
        end do
        ! End new loop over dimensions


	!dipolar DM update scalar quantities 
        do i=1,ngrid       
          		pi_b_mod(ind_cell(i))=sqrt(pi_b(ind_cell(i),1)*pi_b(ind_cell(i),1)+pi_b(ind_cell(i),2)*pi_b(ind_cell(i),2)+pi_b(ind_cell(i),3)*pi_b(ind_cell(i),3))     !dipolar DM
           		xi_b_mod(ind_cell(i))=sqrt(xi_b(ind_cell(i),1)*xi_b(ind_cell(i),1)+xi_b(ind_cell(i),2)*xi_b(ind_cell(i),2)+xi_b(ind_cell(i),3)*xi_b(ind_cell(i),3))     !dipolar DM


	   		w_pol_p(ind_cell(i))=fourpi*pi_b_mod(ind_cell(i))+fourpi*fourpi*pi_b_mod(ind_cell(i))*pi_b_mod(ind_cell(i))/1.78	!dipolar DM, 3)update DERIVATIVE of the potential (carefull here)

	   		w_pol_tograd(ind_cell(i))=(9.56+twopi*pi_b_mod(ind_cell(i))*pi_b_mod(ind_cell(i))+fourpi*fourpi*pi_b_mod(ind_cell(i))*pi_b_mod(ind_cell(i))*pi_b_mod(ind_cell(i))/(3*1.78))-pi_b_mod(ind_cell(i))*w_pol_p(ind_cell(i)) !Here put Lambda=3*alpha^2*(2pi)^2*a_0^2=242 thus Lambda/8pi=9.56 avec alpha=0.8this changes when units changes
        end do

        ! Dipolar DM other new Loop over dimensions
        do idim=1,ndim
           do i=1,ngrid
		if(xi_b_mod(ind_cell(i))==0)then
			f_int(ind_cell(i),idim)=f(ind_cell(i),idim) !11/10/21 try smth else
		else
               		f_int(ind_cell(i),idim)=w_pol_p(ind_cell(i))*xi_b(ind_cell(i),idim)/xi_b_mod(ind_cell(i)) ! dipolar DM 4)update new force

		end if
           end do
        end do
        ! End new loop over dimensions

     end do
     ! End loop over cells
  end do
  ! End loop over grids

!something to be done here for for dipolar DM !!!!!! MPI

111 format('   Entering dipolar_fine for level ',I2)

end subroutine dipolar_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################


subroutine initdipo_fine(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  use mpi_mod
  use constants, only : twopi

  use hydro_commons, ONLY: smallr

  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,icount
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::nx_loc,idim
  real(dp)::dx,dx_loc,scale,fact,fourpi
  real(kind=8)::rho_loc,rho_all,epot_loc,epot_all,modef,small_r !Stahl 17/02/21 need du module de f, small_r when rho=0
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx,ff

!small_r=1E-35

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do



  !----------------------------------------------
  ! Compute initial condition
  !----------------------------------------------
  fourpi=2*twopi
  if(cosmo)fourpi=1.5D0*omega_m*aexp
  fact=-dx_loc**ndim/fourpi/2.0D0



  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Loop over dimensions
        do idim=1,ndim
           do i=1,ngrid						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
														!0)Fancy way of initializing pi_b
													!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
		modef=sqrt(f(ind_cell(i),1)*f(ind_cell(i),1)+f(ind_cell(i),2)*f(ind_cell(i),2)+f(ind_cell(i),3)*f(ind_cell(i),3))	     
		if(modef.ne.0)then
			pi_b(ind_cell(i),idim)=f(ind_cell(i),idim)*(sqrt(1.78/modef)-1)						
	       endif
           end do
        end do
        ! End loop over dimensions
     end do
     ! End loop over cells
  end do
  ! End loop over grids

!something to be done here for for dipolar DM !!!!!! MPI

111 format('   Entering dipolar_fine for level ',I2)

end subroutine initdipo_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################



subroutine gradient_phi(ind_grid,ngrid,ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons

  use constants, only : twopi

  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right

  real(kind=8)::modef !dipolDM besoin du mod de F 23/11/2021

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_phi(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount)
        call interpol_phi(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=phi(igridn(i,ig1)+ih1)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=phi(igridn(i,ig2)+ih2)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=phi(igridn(i,ig3)+ih3)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=phi(igridn(i,ig4)+ih4)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
           f(ind_cell(i),idim)=a*(phi1(i)-phi2(i)) &
                &             -b*(phi3(i)-phi4(i))
!corrections with dipolar
	   if(nstep_coarse<10)then
	   modef=sqrt(f(ind_cell(i),1)**2+f(ind_cell(i),2)**2+f(ind_cell(i),3)**2)	    
	   if(modef.ne.0)then	  
	  	xi_b(ind_cell(i),idim)=f(ind_cell(i),idim)*sqrt(1.78/modef) ! here F=g_n and xi_b=F_tot
	   endif
	  	pi_b(ind_cell(i),idim)=xi_b(ind_cell(i),idim)/(2*twopi)*(1-sqrt(modef/1.78)) !Rem :sqrt (modef/a_0) because |xi|/a_0=sqrt|g_n/a_0|
	  	f_int(ind_cell(i),idim)=2*twopi*pi_b(ind_cell(i),idim)*(1+sqrt(modef/1.78))
           end if
        end do
     end do
  end do

end subroutine gradient_phi

subroutine gradient_w(ind_grid,ngrid,ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_w(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount) !dipolar DM
        call interpol_w(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=w_pol_tograd(igridn(i,ig1)+ih1)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=w_pol_tograd(igridn(i,ig2)+ih2)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=w_pol_tograd(igridn(i,ig3)+ih3)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=w_pol_tograd(igridn(i,ig4)+ih4)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
           ff2(ind_cell(i),idim)=a*(phi1(i)-phi2(i)) &
                &             -b*(phi3(i)-phi4(i))
        end do
     end do
  end do

end subroutine gradient_w


subroutine grad_pi_bx(ind_grid,ngrid,ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_pi_bx(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount)
        call interpol_pi_bx(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=pi_b(igridn(i,ig1)+ih1,1)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=pi_b(igridn(i,ig2)+ih2,1)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=pi_b(igridn(i,ig3)+ih3,1)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=pi_b(igridn(i,ig4)+ih4,1)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
	   if(idim==1)then
               divpi(ind_cell(i))=a*(phi1(i)-phi2(i)) & !here I put div pix
                  &             -b*(phi3(i)-phi4(i))
!if(divpi(ind_cell(i))>0.1)print*,divpi(ind_cell(i))
           endif
        end do
     end do
  end do

end subroutine grad_pi_bx

subroutine grad_pi_by(ind_grid,ngrid,ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_pi_by(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount) !dipolar DM
        call interpol_pi_by(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=pi_b(igridn(i,ig1)+ih1,2)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=pi_b(igridn(i,ig2)+ih2,2)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=pi_b(igridn(i,ig3)+ih3,2)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=pi_b(igridn(i,ig4)+ih4,2)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
	   if(idim==2)then
               divpi(ind_cell(i))=divpi(ind_cell(i))+a*(phi1(i)-phi2(i)) & !Here I add div piy
                  &             -b*(phi3(i)-phi4(i))
!if(nstep>0)print*,divpi(ind_cell(i))
           endif
        end do
     end do
  end do

end subroutine grad_pi_by

subroutine grad_pi_bz(ind_grid,ngrid,ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_pi_bz(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount) !dipolar DM
        call interpol_pi_bz(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=pi_b(igridn(i,ig1)+ih1,3)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=pi_b(igridn(i,ig2)+ih2,3)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=pi_b(igridn(i,ig3)+ih3,3)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=pi_b(igridn(i,ig4)+ih4,3)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
	   if(idim==3)then
               divpi(ind_cell(i))=divpi(ind_cell(i))+a*(phi1(i)-phi2(i)) & !Here I add div piz
                  &             -b*(phi3(i)-phi4(i))
!if(nstep>0)print*,divpi(ind_cell(i)),rho(ind_cell(i))
           endif
        end do
     end do
  end do

end subroutine grad_pi_bz

! Stahl 25/01/2021 : new routine div pi
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine div_pi(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  use mpi_mod
  use constants, only : twopi
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,icount
  !----------------------------------------------------------
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::nx_loc,idim
  real(dp)::dx,dx_loc,scale,fact,fourpi
  real(kind=8)::rho_loc,rho_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx,ff


  if(numbtot(1,ilevel)==0)return
 ! if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do



     ! Loop over myid grids by vector sweeps
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Compute gradient of potential
        call grad_pi_bx(ind_grid,ngrid,ilevel,icount)  !dipole DM compute div comp by comp everything is stored in divpi
        call grad_pi_by(ind_grid,ngrid,ilevel,icount)  
        call grad_pi_bz(ind_grid,ngrid,ilevel,icount)      
     end do
     ! End loop over grids
!if(nstep>0)print*,divpi  



end subroutine div_pi


subroutine grad_Fx(ind_grid,ngrid,ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_Fx(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount)
        call interpol_Fx(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=f(igridn(i,ig1)+ih1,1)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=f(igridn(i,ig2)+ih2,1)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=f(igridn(i,ig3)+ih3,1)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=f(igridn(i,ig4)+ih4,1)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
               f2(ind_cell(i),idim)=a*(phi1(i)-phi2(i)) & !here stored (dxpix,dypix,dzpix)
                  &             -b*(phi3(i)-phi4(i))
        end do
     end do
  end do

end subroutine grad_Fx

subroutine grad_Fy(ind_grid,ngrid,ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_Fy(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount) !dipolar DM
        call interpol_Fy(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
	      phi1(i)=f(igridn(i,ig1)+ih1,2)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=f(igridn(i,ig2)+ih2,2)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=f(igridn(i,ig3)+ih3,2)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=f(igridn(i,ig4)+ih4,2)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
!Some idea of tests to look for bugs
!if(ind_cell(i)==37458)print*,ind_cell(i),idim,f3(ind_cell(i),idim),'before'
!if(ind_cell(i)==37459)print*,ind_cell(i),idim,f3(ind_cell(i),idim),'before'
!if(ind_cell(i)==10262133)print*,ind_cell(i),idim,f3(ind_cell(i),idim),'before'
!if(ind_cell(i)==10263133)print*,ind_cell(i),idim,f3(ind_cell(i),idim),'before'
               f3(ind_cell(i),idim)=a*(phi1(i)-phi2(i)) &
                  &             -b*(phi3(i)-phi4(i))
!if(ind_cell(i)==37458)print*,ind_cell(i),idim,f3(ind_cell(i),idim),a,phi1(i),-phi2(i),-b,phi3(i),phi4(i),a*(phi1(i)-phi2(i))-b*(phi3(i)-phi4(i)),'after'
!if(ind_cell(i)==37459)print*,ind_cell(i),idim,f3(ind_cell(i),idim),a,phi1(i),-phi2(i),-b,phi3(i),phi4(i),a*(phi1(i)-phi2(i))-b*(phi3(i)-phi4(i)),'after'
!if(ind_cell(i)==10262133)print*,ind_cell(i),idim,f3(ind_cell(i),idim),a,phi1(i),-phi2(i),-b,phi3(i),phi4(i),a*(phi1(i)-phi2(i))-b*(phi3(i)-phi4(i)),'after'
!if(ind_cell(i)==10263133)print*,ind_cell(i),idim,f3(ind_cell(i),idim),a,phi1(i),-phi2(i),-b,phi3(i),phi4(i),a*(phi1(i)-phi2(i))-b*(phi3(i)-phi4(i)),'after'
        end do
     end do
  end do

end subroutine grad_Fy

subroutine grad_Fz(ind_grid,ngrid,ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_Fz(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount) !dipolar DM
        call interpol_Fz(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=f(igridn(i,ig1)+ih1,3)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=f(igridn(i,ig2)+ih2,3)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=f(igridn(i,ig3)+ih3,3)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=f(igridn(i,ig4)+ih4,3)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
               f4(ind_cell(i),idim)=a*(phi1(i)-phi2(i)) &
                  &             -b*(phi3(i)-phi4(i))
        end do
     end do
  end do

end subroutine grad_Fz






