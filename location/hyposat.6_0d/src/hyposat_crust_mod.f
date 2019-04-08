c
c  subroutine version of a SCRIPS program called 'getCN1point'
c  to get the parameters of the model CRUST 1.0 as downloaded from
c
c  http://igppweb.ucsd.edu/~gabi/crust1.html
c
c  Because we do not need any densities to locate an event,
c  this parameter was thrown out from the original data and program
c  files.
c
c     NORSAR, Johannes Schweitzer, November 2017
c
c     Softsediment layers thinner than topography  corrected
c     September 2016
c
c
c layer one and two flipped, after the read statement!
c layer 1: water
c layer 2: ice
c

      subroutine get_mod_c10(itrue,inum,typctl,ierr)

      save

      real*8 alpha1

      integer   typctl,ierr,itrue,inum

c
c     typctl  - verbosity level
c
c     inum    - switch
c               = 1 if model used for crustal corrections
c               = 2 if model used for regional travel times
c
c     itrue   - switch
c               = 0 if model used without topography
c               = 1 if 'true' model is used including topography
c
c     ierr    - ne.0 if any error occured
c

c
      include 'model.h'
c

      parameter(np=9,nlo=360,nla=180)

      real*8    x(4),y(4),u(4),vp(2,maxla),vs(2,maxla),dcolo2,
     *          dz(2,maxla),bilinear,dzw(2),dcolo,dcola,
     *          ele(2),dum,dum2,fac,avp(np,nla,nlo),
     *          avs(np,nla,nlo),athi(np,nla,nlo),bnd(np,nla,nlo),
     *          topo(nla,nlo)

      character file*120

c
      include 'modelc.h'
c

      character*120 file_check

c     grid resolution of crustal model
c     CRUST1.0 has a resolution of 1 degree bwetween 
c     -179.5 and 179.5 long and .89.5 and 89.5 lat
c

      ierr = 0

      if(iread.eq.0) then 

         file = file_check('crust1.vp')
         open(65,file=trim(file),err=999)

         file = file_check('crust1.vs')
         open(66,file=trim(file),err=999)

         file = file_check('crust1.bnds')
         open(67,file=trim(file),err=999)

         if(typctl.gt.4) then
            print*,' ... reading global crustal model CRUST1.0 ...'
         endif

         do 101 j=1,nla
            do 102 i=1,nlo

               read(65,*)(avp(l,j,i),l=1,np)
               read(66,*)(avs(l,j,i),l=1,np)
               read(67,*)(bnd(l,j,i),l=1,np)

               do 103 l = 1,np-1
                  athi(l,j,i) = bnd(l,j,i) - bnd(l+1,j,i)
                  if(dabs(athi(l,j,i)) .lt. 1.d-2) athi(l,j,i) = 0.d0
                  if(dabs(avp(l,j,i))  .lt. 1.d-2) athi(l,j,i) = 0.d0
103            continue

               topo(j,i) = bnd(1,j,i)
               athi(np,j,i) = 20.d0

102         continue
101     continue

       iread = 1
       close (65)
       close (66)
       close (67)
      endif

*-------------------
c     now look up coordinates

      do 80 k = 1,inum

         if(k.eq.1) then
            dcola = 90.d0 - elatc
            dcolo = 180.d0 + alpha1(elonc)
            dcolo2 = alpha1(dcolo+180.d0)
c           print *,k,dcola,dcolo,dcolo2,elatc,elonc
         else if(k.eq.2) then
            dcola = 90.d0 - elat2
            dcolo = 180.d0 + alpha1(elon2)
            dcolo2 = alpha1(dcolo+180.d0)
c           print *,k,dcola,dcolo,dcolo2,elat2,elon2
         endif

         klat = 0

         ilat1 = int (dcola)
         y(1)  = dble(int(dcola))
         ilat2 = ilat1 
         y(2) = y(1)

         if((dcola - ilat1 + 1.d0) .ge. 0.5d0 ) then
            ilat3 = ilat1 + 1
            y(3) = y(1) + 1.d0
         else
            ilat3 = ilat1 - 1
            y(3) = y(1) - 1.d0
         endif

         if(ilat3.eq.0) then
            ilat3 = 1
            klat  = 1
         endif

         if(ilat3.gt.nla) then
            ilat3 = nla
            klat  = 1
         endif

         ilat4 = ilat3
         y(4) = y(3)

         ilon1 = int (dcolo) + 1
         x(1) = dble(int(dcolo))

         if((dcolo - ilon1 + 1.d0) .ge. 0.5d0 ) then
            ilon2 = ilon1 + 1
            x(2) = x(1) + 1.d0
            klon = 0
         else
            ilon2 = ilon1 - 1
            x(2) = x(1) - 1.d0
            klon = 1
         endif

         if(ilon2.gt.nlo) ilon2 = ilon2 - nlo
         if(ilon2.lt.1)   ilon2 = nlo + ilon2 

         if(klat.eq.0) then
            ilon3 = ilon2
            x(3) = x(2)
            ilon4 = ilon1
            x(4) = x(1)
         else
            ilon4 = int (dcolo2)
            x(4) = dble(int(dcolo2))
            if(klon.eq.0) then
               ilon3 = ilon4 - 1
               x(3) = x(4) - 1.d0
            else
               ilon3 = ilon4 + 1
               x(3) = x(4) + 1.d0
            endif

            if(ilon3.gt.nlo) ilon3 = ilon3 - nlo
            if(ilon3.lt.1)   ilon3 = nlo + ilon3 

         endif

c        print *,ilat1,ilat2,ilat3,ilat4
c        print *,ilon1,ilon2,ilon3,ilon4

         u(1) = dble(topo(ilat1,ilon1))
         u(2) = dble(topo(ilat2,ilon2))
         u(3) = dble(topo(ilat3,ilon3))
         u(4) = dble(topo(ilat4,ilon4))

         indbi = 0
         ele(k) = bilinear(x,y,u,dcolo,dcola,indbi)
         indbi = 1

      jmod = 1

      iconr = 0
      do 70 i=1,np

         u(1)  = dble(avp(i,ilat1,ilon1))
         u(2)  = dble(avp(i,ilat2,ilon2))
         u(3)  = dble(avp(i,ilat3,ilon3))
         u(4)  = dble(avp(i,ilat4,ilon4))
         vp(k,i) = bilinear(x,y,u,dcolo,dcola,indbi)

         u(1)  = dble(avs(i,ilat1,ilon1))
         u(2)  = dble(avs(i,ilat2,ilon2))
         u(3)  = dble(avs(i,ilat3,ilon3))
         u(4)  = dble(avs(i,ilat4,ilon4))
         vs(k,i) = bilinear(x,y,u,dcolo,dcola,indbi)

         u(1)  = dble(athi(i,ilat1,ilon1))
         u(2)  = dble(athi(i,ilat2,ilon2))
         u(3)  = dble(athi(i,ilat3,ilon3))
         u(4)  = dble(athi(i,ilat4,ilon4))

         dz(k,i) = bilinear(x,y,u,dcolo,dcola,indbi)

c
c     water layer at top!
c
      if(i.eq.1) then
         if(dz(k,i).ge.1.d-3) then
          dzw(k)=dz(k,i)
          dz(k,i) = 0.d0
         endif
      endif

      if(k.eq.1) elev  = ele(k)
      if(k.eq.2) elev2 = ele(k)

      if(typctl.gt.8) then
         if(i.eq.np) then
           if(k.eq.1) then
             print *,'latitude, longitude, elevation: ',
     +         elatc,elonc,elev,inum,itrue
          else if(k.eq.2) then
             print *,'latitude, longitude, elevation: ',
     +         elat2,elon2,elev2,inum,itrue
          endif
         print 793,'Mantle below Moho: ave. vp, vs:  ',
     +        vp(k,i),vs(k,i)
 793  format(a,2x,2f11.4)
      endif

       if (i.eq.1) then
            print *,' ' 
            print *,' 8-layer crustal model (thickness, vp,vs)'
       endif
         if(i.lt.np) print 794, dz(k,i),vp(k,i),vs(k,i)
794      format (3f10.4)
      endif

      if(k.eq.inum) then

         if(inum.eq.1) then

            if(dz(k,i).ge.1.d-3) then
          
             jmod1 = jmod +1
             if(jmod.eq.1) then
                z(jmod)    = 0.d0
                if(itrue.eq.0) then
                   z(jmod1)   = dz(k,i) - ele(k)
                   if(z(jmod1).le.1.d-3) z(jmod1) = 1.d-3
                else
                   z(jmod1)   = dz(k,i)
                endif
             else
                z(jmod)    = z(jmod-1)
                z(jmod1)   = z(jmod)+dz(k,i)
             endif

             v0(1,jmod) = vp(k,i)
             v0(1,jmod1)= vp(k,i)

             v0(2,jmod) = vs(k,i)
             v0(2,jmod1)= vs(k,i)

             azo(jmod)  = ' '
             azo(jmod1) = ' '

             if(v0(1,jmod).ge.6.3d0 .and.iconr.eq.0) then
                azo(jmod-1) = 'CONR'
                iconr = 1
             endif

             if(i.eq.np) then
                azo(jmod-1) = 'MOHO'
                jmod = jmod1
             else
                jmod = jmod1+1
             endif

            endif

         else if(inum.eq.2) then

            if(dz(1,i).ge.1.d-3 .or. dz(2,i).ge.1.d-3) then

             jmod1 = jmod +1
             if(jmod.eq.1) then
                z(jmod)    = 0.d0
                if(itrue.eq.0) then
                   z(jmod1)   = (dz(1,i)-ele(1)+dz(2,i)-ele(2))/2.d0
                     if(z(jmod1).le.1.d-3) z(jmod1) = 1.d-3
                else
                   z(jmod1)   = (dz(1,i) + dz(2,i) ) / 2.d0
                endif
               else
                z(jmod)    = z(jmod-1)
                z(jmod1)   = z(jmod)+(dz(1,i)+dz(2,i))/2.d0
               endif

             fac  = 0.d0
             dum  = 0.d0
             dum2 = 0.d0

               do 59 i2 =1,2
                if(dz(i2,i).ge.1.d-3) then
                 fac  = fac  + 1.d0
                 dum  = dum  + vp(i2,i)
                 dum2 = dum2 + vs(i2,i)
                endif
59             continue

             v0(1,jmod) = dum  /fac
             v0(1,jmod1)= v0(1,jmod)
             v0(2,jmod) = dum2 /fac
             v0(2,jmod1)= v0(2,jmod)

             azo(jmod)  = ' '
             azo(jmod1) = ' '

             if(v0(1,jmod).ge.6.3d0 .and.iconr.eq.0) then
                azo(jmod-1) = 'CONR'
                iconr = 1
             endif

             if(i.eq.np) then
                azo(jmod-1) = 'MOHO'
                jmod = jmod1
             else
                jmod = jmod1+1
             endif

            endif

         endif

      endif
          
70    continue

80    continue

      zmax = z(jmod)

      if(typctl.gt.8) then
       do 90 i=1,jmod
       print*,i,z(i),v0(1,i),v0(2,i),azo(i)
90       continue
      endif

      return

999   print *,'Something wrong with CRUST 1.0 - input files'
      ierr = 99
      return

      end
c
c  subroutine version of a modified SCRIPS program called 'getCNpoint'
c  to get the parameters of the crust as descriped in:
c  Mooney, Laske and Masters, Crust 5.1: a global
c  crustal model at 5x5 degrees, JGR, January 1998
c
c  Because we do not need any densities to locate an event,
c  this parameter was thrown out from the original data and program
c  files.
c
c     NORSAR, Johannes Schweitzer, March 1999
c
C    here modified to read only the crusts of a standard Earth model
c

      subroutine get_mod_global(typctl,ierr)

      save
c
c     mtype = maximum number allowed different crutsal models of
c             standadrd spherical Earth models

      parameter (mctyp=20)
      integer   typctl,ierr

      real*8 dz, fvel(mctyp,8),fvels(mctyp,8),fthi(mctyp,8)

c
c     typctl  - verbosity level
c
c     ierr    - ne.0 if any error occured
c

      include 'model.h'

      character file*120

      include 'modelc.h'

      character*2 ctype(mctyp)
      character*120 file_check

      file = file_check('std_crusts.dat')
      open(65,file=trim(file))

c...  read in key for crust types and all models
c...............................
      read(65,'(a)')
      read(65,'(a)')
      read(65,'(a)')
      read(65,'(a)')
      if(typctl.gt.4) then
         print*,' ... reading global crustal model file ...'
      endif

      do 101 i=1,mctyp+1
         if(i.gt.mctyp) then
           print *,'Maximum number (',mctyp,') of standadrd spherical',
     +             ' Earth models reached! See ',trim(file)
           stop
         endif
         read(65,'(a)',end=102) ctype(i)
         read(65,*)(fvel(i,l),l=1,8)
         read(65,*)(fvels(i,l),l=1,8)
         read(65,*)(fthi(i,l),l=1,7)
101   continue

102   ityp = i - 1

*-------------------

      do 200 i = 1,ityp

         if (mtyp(1:2).eq.ctype(i)) then
             k = i
             go to 202
          endif

200   continue

      print *,'Could not find standard model type ',mtyp
      go to 999

202   continue

      elev = 0.d0
      iconr = 0
      jmod  = 1

      do 70 i=1,8

c
      dz = fthi(k,i)
      if(i.eq.8) dz = 20.d0

      if(dz.ge.1.d-3) then
          
         jmod1 = jmod + 1
         if(jmod.eq.1) then
            z(jmod)    = 0.d0
            z(jmod1)   = dz
         else
            z(jmod)    = z(jmod-1)
            z(jmod1)   = z(jmod)+dz
         endif

         v0(1,jmod) = fvel(k,i)
         v0(1,jmod1)= fvel(k,i)

         v0(2,jmod) = fvels(k,i)
         v0(2,jmod1)= fvels(k,i)

         azo(jmod)  = ' '
         azo(jmod1) = ' '

         if(v0(1,jmod).ge.6.3d0 .and.iconr.eq.0) then
            azo(jmod-1) = 'CONR'
            iconr = 1
         endif

         if(i.eq.8) then
            azo(jmod-1) = 'MOHO'
            jmod = jmod1
         else
            jmod = jmod1+1
         endif

      endif
          
70    continue

      zmax = z(jmod)

c
      if(typctl.gt.8) then
        print *,'Standard Earth model : ',mtyp
        do 90 i=1,jmod
           print*,i,z(i),v0(1,i),v0(2,i),azo(i)
90      continue
      endif
      ierr = 0
      return

999   ierr = 99
      return

      end
