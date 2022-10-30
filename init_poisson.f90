! Stahl 19/08/2020 : modification for dipolar dark matter

subroutine init_poisson
  use pm_commons
  use amr_commons
  use poisson_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer :: info,info2,dummy_io
  integer,parameter::tag=1114
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::ilevel2,numbl2,ilun,ibound,istart
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx,xx2,xx3,xx4
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu

  if(verbose)write(*,*)'Entering init_poisson'

  !------------------------------------------------------
  ! Allocate cell centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(rho (1:ncell))
  allocate(rhod (1:ncell))
  allocate(phi (1:ncell))
  allocate(phi_old (1:ncell))
  allocate(w_pol_p (1:ncell)) !dipolar potential prime
  allocate(w_pol_tograd (1:ncell)) !grad of this gives the contribution to xi_b
  allocate(w_pol_old (1:ncell)) !dipolar potential
  allocate(pi_b_mod (1:ncell)) 
  allocate(xi_b_mod (1:ncell))
  allocate(divpi (1:ncell))  ! divergence (achtung scalar quantity)
  allocate(f   (1:ncell,1:3))
  allocate(f_old   (1:ncell,1:3)) !store old force 
  allocate(xi_b (1:ncell,1:3)) !need the polarization vector
  allocate(xi_b_dot (1:ncell,1:3)) !need the polarization vector
  allocate(pi_b (1:ncell,1:3)) !need the polarization vector
  allocate(pi_dot (1:ncell,1:3)) ! store the time variations of pi
  allocate(pi_b_old (1:ncell,1:3)) !need the polarization vector
  allocate(f_int   (1:ncell,1:3)) !internal force 
  allocate(ff2   (1:ncell,1:3)) !grad sourcing xi_b
  allocate(f2   (1:ncell,1:3)) !div of xi_b
  allocate(f3   (1:ncell,1:3)) !div of xi_b
  allocate(f4   (1:ncell,1:3)) !div of xi_b
!  allocate(my_test1 (1:ncell,1:3)) ! array to check if time evolution has been implemented
!  allocate(my_test2 (1:ncell,1:3)) ! array to check if time evolution has been implemented
!  allocate(my_test3 (1:ncell,1:3)) ! array to check if time evolution has been implemented
  rho=0; rhod=0; phi=0; f=0 ; f_old=0 ; xi_b=0 ; xi_b_dot=0 ; w_pol_p=0 ; w_pol_tograd=0 ; w_pol_old=0 ; f_int=0 ; ff2=0 ; divpi=0 ; f2=0 ; f3=0; f4=0 ; pi_b_mod=0 ; xi_b_mod=0 ; pi_b=0 ; pi_dot=0 ; pi_b_old=0 !initialize variables to zero
  if(cic_levelmax>0)then
     allocate(rho_top(1:ncell))
     allocate(rhod_top(1:ncell))
     rho_top=0
     rhod_top=0
  endif

  !------------------------------------------------------
  ! Allocate multigrid variables
  !------------------------------------------------------
  ! Allocate communicators for coarser multigrid levels
  allocate(active_mg    (1:ncpu,1:nlevelmax-1))
  allocate(emission_mg  (1:ncpu,1:nlevelmax-1))
  do ilevel=1,nlevelmax-1
     do i=1,ncpu
        active_mg   (i,ilevel)%ngrid=0
        active_mg   (i,ilevel)%npart=0
        emission_mg (i,ilevel)%ngrid=0
        emission_mg (i,ilevel)%npart=0
     end do
  end do
  allocate(safe_mode(1:nlevelmax))
  safe_mode = .false.

  !--------------------------------
  ! For a restart, read poisson file
  !--------------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+10
     call title(nrestart,nchar)
     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/grav_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/grav_'//TRIM(nchar)//'.out'
     endif
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)

     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     if(ndim2.ne.ndim+1)then
        if(ndim2.ne.ndim)then
           write(*,*)'File poisson.tmp is not compatible'
           write(*,*)'Found   =',ndim2
           write(*,*)'Expected=',ndim+1
           call clean_stop
        else
           if(myid==1) write(*,*)'Assuming pre commit bce4454 output format'
        endif
     end if
     do ilevel=1,nlevelmax2
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           read(ilun)ilevel2
           read(ilun)numbl2
           if(numbl2.ne.ncache)then
              write(*,*)'File poisson.tmp is not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xx(1:ncache))
              allocate(xx2(1:ncache))
             ! allocate(xx3(1:ncache))
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do
              ! Loop over cells
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 ! Read potential
                 read(ilun)xx
                 do i=1,ncache
                    phi(ind_grid(i)+iskip)=xx(i)
                 end do
                 ! Read force
                 do ivar=1,ndim
                    read(ilun)xx
                    do i=1,ncache
                       f(ind_grid(i)+iskip,ivar)=xx(i)
                    end do
                 end do
		 read(ilun)xx2
                 do i=1,ncache
                    pi_b_mod(ind_grid(i)+iskip)=xx2(i)
                 end do
		! Read new force
                 do ivar=1,ndim
                    read(ilun)xx2
                    do i=1,ncache
                       f_int(ind_grid(i)+iskip,ivar)=xx2(i)
                    end do
                 end do
               ! Read polar vector
                 do ivar=1,ndim
                    read(ilun)xx2
                    do i=1,ncache
                       xi_b(ind_grid(i)+iskip,ivar)=xx2(i)
                    end do
                 end do
               ! Read derivative polar vector
                 do ivar=1,ndim
                    read(ilun)xx2
                    do i=1,ncache
                       xi_b_dot(ind_grid(i)+iskip,ivar)=xx2(i)
                    end do
		! Read pi_b_mod
                 end do
              end do
              deallocate(ind_grid,xx,xx2)
           end if
        end do
     end do
     close(ilun)

     ! Send the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif

#ifndef WITHOUTMPI
     if(debug)write(*,*)'poisson.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'POISSON backup files read completed'
  end if

end subroutine init_poisson


