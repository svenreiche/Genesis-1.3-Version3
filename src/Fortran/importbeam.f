
c
      subroutine importdispersion 
c     ================================================================= 
c     apply dispersion to imported beam file from readpart
c     subroutine supplied by Atoosa Meseck from Bessy
c     ------------------------------------------------------------------
c
      include 'genesis.def'
      include 'input.cmn'
      include 'particle.cmn'
c
      real*8 iarho,iaphi,mam,imagl_old
      real*8 ypart_old,py_old
      real*8 ma12,ma33,ma34,ma43,ma56
      integer i,ierr
c
      if (ibfield.eq.0.d0) return                     
      if (idril.lt.0.d0) then                      
          ierr=printerr(errinwarn,'IDRIL<0-NO DISP.SECTION')                                     
          return
      endif  
c
      if (igamref.le.0) igamref=gamma0  
c          
      ibfield=abs(ibfield)      
c      
      imagl_old=imagl
      iarho=(igamref*5.11e-04)/(ibfield*0.299793)
      iaphi= asin(imagl/iarho)
      mam=tan(iaphi)/iarho
      imagl=iaphi*iarho 

      ma12=3*idril+4*iarho*sin(iaphi)*cos(iaphi)
     +      +2*idril*cos(iaphi)*cos(iaphi)
           
      ma33= mam*(-10*idril-8*imagl)+
     +      mam*mam*(26*imagl*idril+
     +      15*idril*idril+8*imagl*imagl)+
     +      mam**3*(-12*idril*(imagl**2)-
     +      20*(idril**2)*imagl-7*(idril)**3)
     +      +1+(mam**4)*(4*(idril)**3*imagl+
     +      (2*idril*imagl)**2+idril**4)
           
c ma44=ma33 
           
      ma34=mam*(-28*idril*imagl-
     +      8*(imagl)**2-20*(idril)**2)+
     +      mam*mam*(44*imagl*(idril)**2+
     +      21*(idril)**3+20*idril*imagl**2)+
     +      mam**3*(-16*(idril*imagl)**2-
     +      24*(idril**3)*imagl-8*(idril)**4)
     +      +(mam**4)*(4*(idril)**3*imagl**2+
     +      4*idril**4*imagl+idril**5)+
     +      5*idril+4*imagl

      ma43=tan(iaphi)*(-2*iarho+
     +      2*tan(iaphi)*imagl+tan(iaphi)*idril)*
     +      (2*iarho**2-4*tan(iaphi)*imagl*iarho-
     +      4*tan(iaphi)*idril*iarho+
     +      2*idril*tan(iaphi)*tan(iaphi)*imagl+
     +      (tan(iaphi)*idril)**2)/
     +      iarho**4 


      ma56=8*iarho*sin(iaphi)-
     +      4*iarho*sin(iaphi)*cos(iaphi)+2*idril-
     +      2*idril*cos(iaphi)*cos(iaphi)-4*imagl+
     +      ((5*idril+4*imagl)/igamref**2) 
c
      imagl=imagl_old
c
      do i=1,npart 
           theta(i)=theta(i)+
     +      (ma56*(gamma(i)-igamref)/igamref)*twopi/xlamds/convharm
           xpart(i)=xpart(i)+ma12*px(i)/gamma(i)
           ypart_old=ypart(i)
           py_old= py(i)
           ypart(i)=ma33*ypart_old+ma34*py_old/gamma(i)         
           py(i)=ma43*ypart_old*gamma(i)+ma33*py_old 
      enddo                                       
      return
      end !of import dispersion
c
c      
      subroutine importtransfer 
c     ================================================================= 
c     Transfer matrix calculation supplied by A. Meseck. 
c     ------------------------------------------------------------------
c
      include 'genesis.def'
      include 'io.cmn'
      include 'input.cmn'
      include 'field.cmn'
      include 'sim.cmn'
      include 'particle.cmn'
c
      real*8 ypart_old,py_old,xpart_old,px_old
      real*8 gamma_old,theta_old    
      integer i,ierr
c
      if (trama.eq.0.d0) return
      if (igamref.le.0) igamref=gamma0 

      do i=1,npart 
c  Denormalize      
         px(i)=px(i)/gamma(i)
         py(i)=py(i)/gamma(i)
cc
         xpart_old=xpart(i)
         ypart_old=ypart(i)
         px_old= px(i)
         py_old= py(i)
         theta_old=theta(i)
         gamma_old= gamma(i)

            xpart(i)=itram11*xpart_old+itram12*px_old+
     +      itram13*ypart_old+itram14*py_old+
     +      itram15*theta(i)*xlamds*convharm/twopi+
     +      itram16*(gamma(i)-igamref)/igamref

           px(i)=itram21*xpart_old+itram22*px_old+
     +      itram23*ypart_old+itram24*py_old+
     +      itram25*theta_old*xlamds*convharm/twopi+
     +      itram26*(gamma(i)-igamref)/igamref

         
           ypart(i)=itram31*xpart_old+itram32*px_old+
     +      itram33*ypart_old+itram34*py_old+
     +      itram35*theta_old*xlamds*convharm/twopi+
     +      itram36*(gamma(i)-igamref)/igamref

           py(i)=itram41*xpart_old+itram42*px_old+
     +      itram43*ypart_old+itram44*py_old+
     +      itram45*theta_old*xlamds*convharm/twopi+
     +      itram46*(gamma(i)-igamref)/igamref

          theta(i)=itram55*theta_old+ (itram56*
     +    ((gamma(i)-igamref)/igamref)*twopi/xlamds/convharm)+
     +    (itram51*xpart_old+itram52*px_old+itram53*ypart_old+
     +     itram54*py_old)*twopi/xlamds/convharm

         gamma(i)=(itram61*xpart_old+itram62*px_old+
     +      itram63*ypart_old+itram64*py_old+
     +      itram65*theta_old*xlamds*convharm/twopi)*
     +      igamref + itram66*(gamma(i)-igamref)+igamref

c normalization
           px(i)=px(i)*gamma(i)
           py(i)=py(i)*gamma(i)
cc
      enddo    
                                                                 
      return
      end !of import transfermatrix
c
c            
      function readpart(nread)
c     =================================================================
c     load complete set of particle from file
c     -----------------------------------------------------------------
c
      include 'genesis.def'
      include 'input.cmn'
      include 'particle.cmn'
      include 'io.cmn' 
      include 'sim.cmn'
c
      integer j,i,islice,idel,nread
c

      logical isop
c
      readpart=noerr
      npart0=nread
      npart=nread     !  needed to reset the counter when imported from partfile.

c     apply transfermatrix to particle distribution
      call importtransfer 
c  
c     dispersive section
      call importdispersion    
c
c
c     convert to higher harmonic
c
      if (convharm.gt.1) then
            do i=1,npart
             theta(i)=float(convharm)*theta(i)
            enddo 
      endif
c
c     calculate init. perpendicular velocity (needed in first call of track)
c
      do i=1,npart0      
	btpar(i)=dsqrt(1.d0-(px(i)**2+py(i)**2+1.)/gamma(i)**2)     !parallel velocity
      enddo
c
c     check for particle losses from previous run
c
      idel=0          
      do i=1,npart 
         if (gamma(i).gt.0.) then
            gamma(i-idel)=gamma(i)  
            theta(i-idel)=theta(i)  
            xpart(i-idel)=xpart(i)
            ypart(i-idel)=ypart(i)
            px(i-idel)=px(i)
            py(i-idel)=py(i)
            btpar(i-idel)=btpar(i)
         else    
            idel=idel+1
         endif
      enddo
      npart=npart0-idel
      xcuren=xcuren*float(npart)/float(npart0)

      call chk_loss



c
      return
 100  readpart=printerr(errread,partfile)
      call last
      return
      end
c

