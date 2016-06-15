      function initrun()
c     ==================================================================
c     initialize the run by setting up
c     the precalculated matrices for the field solver and
c     claculating/normalizing some auxiliary variables
c
c     ------------------------------------------------------------------
c
      include 'genesis.def'
      include 'input.cmn'
      include 'particle.cmn'
      include 'sim.cmn'
      include 'field.cmn'
      include 'time.cmn'
      include 'magnet.cmn'
      include 'mpi.cmn'
c
      real*8 xi,initrun,rslp
      integer ip,i
c
c     seeding of random number generator
c
      initrun=0
      xi=ran1(iseed)          !init ran1
c
c      initiate loop for higher harmonics
c

c
      dedz=eloss*delz*xlamd/eev
      xcuren=curpeak
      npart0=npart                            !save total number or particle 
      if (ione4one.ne.0) then                 !one 4 one simulation. Calculate particle number.
        npart=curpeak*xlamds*zsep/ce
        npart0=npart
      endif
c
c     normalizations
c     ------------------------------------------------------------------
c
      xkw0= twopi/xlamd                       !wiggler wavenumber
      xkper0 = xkw0                           !transverse normalisation
c
c     magnetic field
c     ------------------------------------------------------------------   
c
      call magfield(xkper0,1)                   !magnetic field description
c
      if (inorun.ne.0) then
        ip=PRINTERR(ERRGENWARN,'Termination enforced by user')
        call last
      endif
 
c
c     slipping length
c
      nsep=int(zsep/delz)                      !steps between field slippage
      if (nsep.lt.1) nsep=1
      rslp=0;
      do i=1,nstepz
        rslp=rslp+awslip(i)
      enddo
      rslp=nstepz         ! bypass the implementation with awslip which has errors
      nslp=int(rslp/zsep)                         !total slippage steps
      if (mod(nstepz,nsep).ne.0) nslp=nslp+1   !if not added the effective undulator
c      write(*,*) 'nstepz: ',nstepz, ' nslp: ', nslp
c                                              !would be shorter
c     contruct grid properties (grid spacing, precalculated matrices)
c
      dxy=xkw0*2.d0*dgrid/float(ncar-1)    
      xks=twopi/xlamds

c
c     time dependencies
c    
      if (mpi_id.eq.0) then
c         write(*,*) 'Root loads slippage'
         call loadslpfld(nslp)        !input field for first slice and seeding of ran-function
      endif
c
c     scanning
c
      call scaninit   !initialize scanning
c
c      
c
c     matrix initialization
      call getdiag(delz*xlamd,dxy/xkper0,xks)
c
c     clear space charge field for case that space charge is disabled
c
      do ip=1,npart   !clear space charge term
         ez(ip)=0.d0
      end do       ! ip

c
c
      return
      end
c
