      function output(istepz,islice,xkw0)
c     =============================================
c     calls all output function
c     ---------------------------------------------
c
      include 'genesis.def'
      include 'io.cmn'
      include 'diagnostic.cmn'
c
      integer istepz,islice,output
      real*8  xkw0
      character*11 file_id
c
      output=0
      call diagno(istepz)
      call status(istepz,islice)
c
c     output of particle and field distribution is done in C++
c
c      if (islice.le.firstout) return
c      call outfield(istepz,islice,xkw0)
c      call outpart(istepz,islice)
c
      return
      end

c
      subroutine outpart(istepz,islice)
c     ==================================================================
c     output of global parameter (t-independend): 
c     z, wiggler field
c     ------------------------------------------------------------------
c
      include  'genesis.def'
      include  'io.cmn'
      include  'particle.cmn'
      include  'input.cmn'
      include  'work.cmn'
      include  'sim.cmn'
c
      integer iz,istepz,islice
c
c     ------------------------------------------------------------------ 
c     output of t independent values with first slice
c
      if ((ippart.le.0).or.(ispart.le.0)) return   !no output at all
      if (mod(istepz,ippart).ne.0) return          !output ippartth step
      if (mod(islice,ispart).ne.0) return          !output ispartth slice
c
      if (istepz.eq.0) call rpos(0,xpart,ypart)
      call getpsi(p1)
c
      if (npart.lt.npart0) then     ! check for particle loss
         do iz=npart+1,npart0       ! indicate lost particles with neg. energy
            gamma(iz)=-1.
         enddo
      endif
c
      write(npar,rec=irecpar) (gamma(iz),iz=1,npart0)
      write(npar,rec=irecpar+1) (p1(iz),iz=1,npart0)
c      write(npar,rec=irecpar+1) (theta(iz),iz=1,npart0)
      write(npar,rec=irecpar+2) (xpart(iz)/xkper0,iz=1,npart0)
      write(npar,rec=irecpar+3) (ypart(iz)/xkper0,iz=1,npart0)
      write(npar,rec=irecpar+4) (px(iz),iz=1,npart0)
      write(npar,rec=irecpar+5) (py(iz),iz=1,npart0)
      irecpar=irecpar+6
      return
      end
c
c
c
      subroutine status(istepz,islice)
c     ==================================================================
c     let user know % complete at every 10%.
c     ------------------------------------------------------------------
c
      include 'genesis.def'
      include 'input.cmn'
      include 'magnet.cmn'
      include 'io.cmn'
c
      integer istepz,islice
      real*8 xper,yper
c
      if (itdp.ne.0) return
      xper=100.d0*float(istepz)/float(nstepz)
      yper=100.d0*float(istepz-1)/float(nstepz)
      if (mod(xper,10.0d0).lt.mod(yper,10.0d0)) 
     +       write (nlog,20) islice,int(xper)
   20 format ('Slice ', i5,': Simulation ',i3,'% completed.')
c
      return
      end     !status
c
c
c
      subroutine outfield(istepz,islice,xkper0)
c     ==================================================================
c     dump fieldarray  
c     ------------------------------------------------------------------
c
      include 'genesis.def'
      include 'input.cmn'
      include 'io.cmn'
      include 'field.cmn'
c
      integer i,islice,istepz
      integer ioffset,ih,ifile
      real*8 scltmp,xkper0
c
      if ((ipradi.le.0).or.(isradi.le.0)) return   !no output at all
      if (mod(istepz,ipradi).ne.0) return          !output ipradith step
      if (mod(islice,isradi).ne.0) return          !output isradith slice
c
      scltmp=dxy*eev*xkper0/xks/dsqrt(vacimp)   !
      write(nfld,rec=irecfld) (scltmp* dble(crfield(i)),i=1,ncar*ncar)
      write(nfld,rec=irecfld+1) 
     +        (scltmp*dimag(crfield(i)),i=1,ncar*ncar)
      do ih=2,nhloop
          ioffset=(ih-1)*ncar*ncar
          ifile=nfldh(ih-1)
          write(ifile,rec=irecfld) 
     +	      (scltmp/hloop(ih)* dble(crfield(i)),
     +        i=1+ioffset,ncar*ncar+ioffset)
          write(ifile,rec=irecfld+1) 
     +        (scltmp/hloop(ih)*dimag(crfield(i)),
     +        i=1+ioffset,ncar*ncar+ioffset)
      enddo
c
      irecfld=irecfld+2

      return
      end
c


c
      subroutine closefile(nio)
c     =================================================================
c     closing file
c     ---------------------------------------------------------------
c
      logical isop
c
      if (nio.gt.6) then
        inquire(nio,opened=isop) 
        if (isop) close(nio)                   !close history file
      endif  
      return
      end ! of closefile
c
      function opentextfile(file,status,nio)
c     ==================================================================
c     open ascii file (sequential access)
c     ------------------------------------------------------------------ 
c    
c
      include 'genesis.def'
c
      character*(*) file,status
      integer nio

      opentextfile=nio
      open(nio,file=file,status=status,err=100)
      return
 100  opentextfile=printerr(erropen,file)
      return
      end ! of opentextfile
c
      function openbinfile(root,extension,nio,nsize)
c     ==================================================================
c     open binary file (direct access) as addition output file
c     ------------------------------------------------------------------ 
c    
c
      include 'genesis.def'
c
      character*30 root
      character*4  extension
      character*36 filename
      integer nio,nsize,j,jj

      openbinfile=nio
      j=index(root,' ')
      if (j.eq.0) j=31
      j=j-1
      jj=index(extension,' ')
      if (jj.eq.0) jj=5
      jj=jj-1
      filename=root(1:j)//'.'//extension(1:jj)
      open(nio,file=filename,status='unknown',access='direct',
     +    recl=nsize,err=100)
      return
 100  openbinfile=printerr(erropen,filename)
      return
      end ! of openbinfile
c
c


      subroutine first
c     ============================================
c     initial information for user
c     --------------------------------------------
c
      include 'genesis.def'
      include 'io.cmn'
c 
      write(nlog,100) genver,platf
      return

 100  format('-------------------------------',/,
     c       'Genesis 1.3 has begun execution',/,
     c       '(Version ',f3.1,' ',a,')',/)
      end


      subroutine last()
c     ==================================================================
c     called at end of run.
c     closes all files, which must stay open during the run
c     ------------------------------------------------------------------
c
      include 'genesis.def'
      include 'io.cmn'
      include 'mpi.cmn'
c
      integer ih
c
      write (nlog,100)
      call closefile(nout)   !standard output
      call closefile(nfld)   !field output
      do ih=2,nhmax
         call closefile(nfldh(ih))  !harmonic field output      
         call closefile(ndumph(ih)) !dumped harmonic field     
      enddo
      call closefile(npar)   !particle output 
      call closefile(nfin)   !field input
      call closefile(npin)   !particle input
      call closefile(ndump)  !dumped field
      call closefile(ndmp2)  !dumped particle 
      call closefile(ndis)   !input distribution
      call closetimerec
c      
      write (nlog,200) !genesis has finished
c
      if (nlog.ne.6)  call closefile(nlog)   !log file
c
c      this has to be changed!!!!!!! 
c      pipe out the stop command

c      call MPI_Finalize(mpi_err)
c      stop
c
 100  format('***  closing files')
 200  format(/,'Genesis run has finished',/,
     c         '------------------------') 
      end     !last
c

      function printerr(ierr,text)
c     ========================================================
c     print error messages
c     --------------------------------------------------------
      
      include 'genesis.def'
      include 'io.cmn'

      integer ierr
      character*(*) text

      printerr=ierr
      if (ierr.ge.0) return
      goto (10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
     c      25,26,27,28,29,30,31,32,33)
     c     ,iabs(ierr)
 10   write(nlog,100) text
      return
 11   write(nlog,101) text
      return
 12   write(nlog,102) text
      return
 13   write(nlog,103) text
      return
 14   write(nlog,104) text
      return
 15   write(nlog,105) text
      return
 16   write(nlog,106) text
      return
 17   write(nlog,107) text
      return
 18   write(nlog,108) text
      return
 19   write(nlog,109) text
      return
 20   write(nlog,110) text
      return
 21   write(nlog,111) text
      return
 22   write(nlog,112) text
      return
 23   write(nlog,113) text
      return
 24   write(nlog,114) text
      return
 25   write(nlog,115) text
      return
 26   write(nlog,116) text
      return
 27   write(nlog,117) text
      return
 28   write(nlog,118) text
      return
 29   write(nlog,119) text
      return
 30   write(nlog,120) text
      return
 31   write(nlog,121) text
      return
 32   write(nlog,122) text
      return
 33   write(nlog,123) text
      return
c
c     format statements
c
 100  format('***  File-error: ',a,/,
     c       '***  cannot be opened')
 101  format('***  File-error: ',a,/,
     c       '***  cannot be accessed')
 102  format('***  File-error: ',a,/,
     c       '***  cannot be opened',/,
     c       '***  creating template file: template.in')
 103  format('***  File-error: ',a,/,
     c       '***  error in namelise $newrun')
 104  format('***  Scan-warning: conflict with ITDP',/,
     c       '***  using scan-feature',a)
 105  format('***  Scan-warning: conflict with BEAMFILE',/,
     c       '***  ignoring BEAMFILE: ',a)
 106  format('***  Beamfile-warning: size exceeds NSMAX',/,
     c       '***  ',a)
 107  format('***  Input-error: ',a)
 108  format('***  Input-warning: ',a)
 109  format('***  Input-error: cannot convert to individiual input',/,
     c       '***  ',a)
 110  format('***  Numerical-error: boundary exceeded of',a,/,
     c       '***  ignoring exceeding elements')
 111  format('***  Round-warning: section not multiple of XLAMD',/,
     c       '***  MAGINFILE:',a)    
 112  format('***  Extrapolation-warning: exceeding time window of',/,
     c       '***  BEAMFILE by:',a)
 113  format('***  Scan-error: conflict with MAGINFILE:',a,/,
     c       '***  disabling scan-feature') 
 114  format('***  Scan-error: conflict with MAGOUTFILE:',a,/,
     c       '***  disabling scan-feature')
 115  format('***  Warning: particle loss of ',a,'%')
 116  format('***  Warning: external magnet definition too short for '
     c       ,a)
 117  format('***  Error: invalid filename:',a)
 118  format('***  File-error: cannot read from FIELDFILE:',a)
 119  format('***  Warning: ',a)
 120  format('***  Error: cannot run in background mode.',/,
     c       '***  information needed for ',a)
 121  format('***  Error: CRTIME cannot hold slippage field.',/,
     c       '***  see manual for allocating more memory',a)
 122  format('***  Error: unphysical parameter for loading',/,
     c       '***  ',a)
 123  format('***  ',a)
      end
