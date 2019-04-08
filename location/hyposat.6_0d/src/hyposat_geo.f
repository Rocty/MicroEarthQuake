CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C     SUBPROGRAM : DEPI
c
C         PROGRAM COMPUTES :   - EPICENTRAL DISTANCES
C                              - AZIMUTHS AND BACKAZIMUTHS
C
C     From an old main routine of unknwon source (USGS?) with several
C     changes since the 1980s.
C
C     here last changes November 16, 2016 by J. Schweitzer, NORSAR
C
      subroutine depi(rlat,rlon,elat,elon,deld,delk,azid,bazd,d2km)
c
c
c     rlat, rlon = receiver geographical coordinates
c
c     elat, elon = event geographical coordinates
c     
c     deld       = epicentral distance in degree
c
c     delk       = epicentral distance in km
c
c     azid  =  epicenter (elat,elon) to station (rlat,rlon) !!!!!!
c    
c     bazd  =  backazimuth ( from station to event!!)
c
c     d2km  =  scaling factor degrees -> km
c
      IMPLICIT real*8 (A-H,O-Z)

      pi = 4.d0* datan(1.d0)
      pi2= pi/2.d0
      ep=rlat
      el=rlon
      stnp=elat
      stnl=elon
      tol = 1.d-6
C
      DGSPR=180.d0/PI
      RAD=PI/180.d0
C
C     Check, if both coordinates are identical (here
C     assumed a distance less than ~ 1 meter)
C
      IF ( dabs(EP-STNP).lt. tol .and.
     +     dabs(EL-STNL).lt. tol       ) then

        angi  = 0.d0
        range = 0.d0
        az    = 0.d0
        baz   = 0.d0
        d2km  = rad*radloc (EP,1)
        goto 1000

      endif

      R90=90.d0*RAD
C
C **** MAKE ESQ1= 1.d0 FOR SPHERICAL EARTH
C      
c     esq1=1.d0
c     ESQ1=(1.d0-1.d0/298.257d0)**2.d0
c
c     EPIPR=90.d0*RAD-DATAN(ESQ1*dtan(EP*RAD))
c
      EPIPR=r90-convlat(EP,1)*rad
      EPIPC=EPIPR/RAD

C
C   if EL negative, EPILC=360+EL
C
      IF (EL.LT.0.d0) THEN
         EPILC=360.d0+EL
      ELSE
         EPILC=EL
      END IF

cj.s. EPILR=EPILC*RAD
C
C   save pol latitude
C
      E180=(R90-rad*convlat(-90.d0,1))/RAD
C
      if(dabs(dabs(stnp)-90.d0).lt.tol ) stnl=0.d0

      STNPR=R90-rad*convlat(STNP,1)
      STPC=STNPR/RAD

C
C   if STNL negative, STLC=360+STNL
C
      IF (STNL.LT.0.d0) THEN
         STLC=360.d0+STNL
      ELSE
         STLC=STNL
      END IF

      IF (dabs(STPC-E180).lt.tol) THEN
         ANGI=E180-EPIPC
         AZ=180.d0
         BAZ=EL
         GOTO 998
      END IF
      if(dabs(stpc).lt.tol) then
         angi=epipc
         az=180.d0-el
         baz=0.d0
         go to 998
      endif
      if(dabs(epipc).lt.tol) then
         angi=stpc
         baz=0.d0
         az=180.d0-stnl
         go to 998
      endif
      if(dabs(epipc-e180).lt.tol) then
         angi=e180-stpc
         az=stlc
         baz=180.d0
         go to 998
      endif

C
C   DETERMINE POLAR ANGLE
C
      PANG=DABS(STLC-EPILC)

      IF (PANG.GT.180.d0) PANG=360.d0-PANG
      PANG=PANG*RAD

      ANGI=DCOS(EPIPR)*DCOS(STNPR)+DSIN(EPIPR)*DSIN(STNPR)*DCOS(PANG)
      SNANG=DSQRT(DABS(1.d0-q2(ANGI)))
      if(dabs(ANGI).lt.tol) then
        angi = pi2
      else
        ANGI=DATAN(SNANG/ANGI)
      endif

C
C   ANGI negative?
C
      IF (ANGI.LT.0.d0) ANGI=ANGI+PI

      if(dabs(angi).gt.tol) then
         AZ=(DCOS(STNPR)-DCOS(EPIPR)*DCOS(ANGI))/(DSIN(EPIPR)
     *      *DSIN(ANGI))
      else
         AZ= 0.d0
      endif

      SNAZ=DSQRT(DABS(1.d0-q2(AZ)))
      if (dabs(az).lt.tol) then
         az = pi2
      else
         AZ=DATAN(SNAZ/AZ)
      endif

C
C   AZ negative ?
C
      IF (AZ.LT.0.d0) AZ=AZ+PI
      AZ=DGSPR*AZ

      if(dabs(angi).gt.tol) then
         BAZ=(DCOS(EPIPR)-DCOS(STNPR)*DCOS(ANGI))/(DSIN(STNPR)
     *      *DSIN(ANGI))
      else
         BAZ= 0.d0
      endif

      SNBAZ=DSQRT(DABS(1.d0-q2(BAZ)))
      if(dabs(baz).lt.tol) then
         baz= pi2
      else
         BAZ=DATAN(SNBAZ/BAZ)
      endif

C
C   BAZ negative ?
C
      IF (BAZ.LT.0.d0) BAZ=BAZ+PI
      BAZ=BAZ*DGSPR
      ANGI=ANGI*DGSPR

C
C   ADJUST AZIMUTH AND BACKAZIMUTH
C
c      IF (STLC-EPILC) 72,74,73
c 72   IF (STLC-EPILC+180.d0) 80,77,81
c 73   IF (STLC-EPILC-180.d0) 80,77,81
c 74   IF (STPC-EPIPC) 75,78,76

      dum1 = STLC-EPILC
      IF (dum1.le.-tol) goto 72
      IF (dum1.ge.tol) goto 73
      goto 74

 72   dum2 = dum1 + 180.0d0
      IF (dum2.le.-tol) then
         BAZ=360.d0-BAZ
         goto 998
      else if (dum2.ge.tol) then
         AZ=360.d0-AZ
         goto 998
      endif
      goto 77

 73   dum3 = dum1 - 180.0d0
      IF (dum3.le.-tol) then
         BAZ=360.d0-BAZ
         goto 998
      else if (dum3.ge.tol) then
         AZ=360.d0-AZ
         goto 998
      endif
      goto 77

 74   dum4 = STPC-EPILC
      IF (dum4.le.-tol) then
         AZ=0.d0
         BAZ=180.d0
      else if (dum4.ge.tol) then
         AZ=180.d0
         BAZ=0.d0
      else
         AZ=0.d0
         BAZ=0.d0
      endif
      goto 998


c 77   IF (STPC+EPIPC-E180) 78,78,79

 77   dum5 = STPC+EPIPC-E180
      IF (dum5.ge.tol) then
         AZ=180.d0
         BAZ=180.d0
      else
         AZ=0.d0
         BAZ=0.d0
      endif

C
C   CALCULATE DISTANCE IN KILOMETER
C

c
c     First we need a 'mean' radius to convert distance in
c     [deg] on the sperical Earth into [km] with respect
c     to geographical latitude.
c

998   continue

      if (angi.ge.100.d0) then
         drad =  angi / 200.d0
      else if (angi.gt.10.d0 .and. angi.lt.100.d0) then
         drad =  angi / 100.d0
      else if (angi.gt.0.05d0 .and. angi.le.10.d0) then
         drad = angi / 15.d0
      else 
         srad = radloc((ep+stnp)/2.d0,1)
         ndrad = 1
         go to 9991
      endif
      srad = 0.d0

      ndrad = int (0.5d0 + angi / drad) + 1

      ind = 1

      do 999 irad = 1,ndrad
         radl = dble(irad-1) * drad
         call delazd(ep,el,az,radl,ind,elat2,elon2)
         srad = srad + radloc (elat2,1)
999   continue

9991  d2km = rad * srad / dble(ndrad)
      RANGE=ANGI*d2km

1000  deld=angi
      delk=range
      bazd=az
      azid=baz
c
c
      return
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c     Subroutine  delazd.f calculates geographical coordinates from a 
c     source point to another point using azimuth and distance (either 
c     in degree or in km). Input and output are in geographical 
c     coordinates, internally, all calculations are in geocentric 
c     coordinates using the known ellipticity formula.
c
c       date:  July 25, 1996 Johannes Schweitzer
c
c       in original from December 12, 1986
c
c       input:
c
c       elat1,elon1   -   geographical coordinates of starting point
c
c       azi           -   azimuth of new point
c
c       dis           -   distance of new point in degree or km
c
c       ind           -   = 1 dis in degree
c                         = 2 dis in km
c
c       elat2,elon2   -   geographical coordinates of the target
c       
c
      subroutine delazd(elat1,elon1,azi,dis,ind,elat2,elon2)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION P1(2),PX(2),P2(2)
      real*8 radloc
c
      ind2 = 0

      p1(1) = elon1
      p1(2) = elat1
      p2(1) = dis
      p2(2) = azi
c
      PI=4.0d0*DATAN(1.0d0)
      PI2=2.0d0*PI
      RAD=PI/180.d0

      rel = radloc (elat1,1)

10    continue

      if(ind.eq.2) then
         del=p2(1)/(rel*rad)
      else if (ind.eq.1) then
         del=p2(1)
      else
         print *, 'Wrong index for distance'
         return
      endif
c
c     with ellipticity correction esq
c     esq=(1.0d0-1.0d0/298.257d0)**2.0d0
c
c     without ellipticity correction esq=1.0d0 
c
c     A=90.0d0*rad-datan(dtan(P1(2)*RAD)*esq)
c     now done in convlat()

      a= (90.0d0-convlat(P1(2),1))*rad
      C1=DCOS(A)
      C2=DSIN(A)
      AZ=P2(2)*RAD
      bx=del*rad
      C4=DCOS(AZ)
      C5=C4*C2
      IF(P1(2).EQ.-90.0d0) THEN
         CX=PI-BX
         PX(1)=AZ/RAD
         IF(CX.LT.0.0d0) PX(1)=PX(1)+180.0d0
         GO TO 14
      ENDIF
      IF(P1(2).EQ.90.0d0) THEN
         CX=BX
         PX(1)=(PI-AZ)/RAD
         IF(CX.GT.PI) PX(1)=PX(1)+180.0d0
         GO TO 14
      ENDIF
      C6=DCOS(BX)
      C10=C1*C6+C5*DSIN(BX)
      CX=F2(C10,2)
      C11=(C6-C1*C10)/(C2*DSIN(CX))
      C12=DABS(C11)
      IF(C12.GT.1.d0) C11=C11/C12
      BETX=F2(C11,2)
      IF(BX.GT.PI) BETX=PI2-BETX
      IF(AZ.GT.PI2) AZ=AZ-PI2
      IF(AZ.LT.0.d0) AZ=AZ+PI2
      IF(AZ.LE.PI.AND.AZ.GE.0.0d0) PX(1)=P1(1)+BETX/RAD
      IF(AZ.GT.PI.AND.AZ.LE.PI2) THEN
         PX(1)=P1(1)-BETX/RAD
      ENDIF

c14    IF(PX(1).GT.180.0d0) PX(1)=PX(1)-360.0d0
c      IF(PX(1).LT.-180.0d0) PX(1)=360.0d0+PX(1)
c     px(2)=datan(dtan(pi/2.0d0-cx)/esq)/rad

14    px(1) = alpha1(px(1))
      px(2)=convlat((90.d0-cx/rad),2)
      elat2 = px(2)
      elon2 = px(1)

      if(ind.eq.2 .and. ind2.eq.0) then
         ind2 = 1
         rel2 = radloc (elat2,1)
         rel  = (rel + rel2) / 2.d0
         go to 10
      endif

      return
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     FUNCTION F2  limits the ACOS and ASIN for small numerical 
C     instabilities.
C
      FUNCTION F2(A,IND)
      IMPLICIT real*8 (A-H,O-Z)
      B=DABS(A)
      C=A
      IF(B.GT.1.0d0.AND.B.LT.1.00005d0) C=DSIGN(1.d0,A)
      GOTO (1,2),IND
1     F2=DASIN(C)
      RETURN
2     F2=DACOS(C)
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     This function DIRDEL (DIRect DELta) calculates the epicentral 
c     distance for an upgoing phase at a source depth greater than 0.
c
c     We use the AK135 as standard velocity model and assume that the
c     distance differences for such phases in different models is
c     relatively small.
c
c     Model AK135 is here hard coded via DATA statements.
c
c     input:
c     
c               p      -  ray paramter of upgoing phase
c
c               zo     -  source depth
c
c               fa     -  factor how often this phase passes the
c                         structure (here usually fa = 1.d0)
c
c               type   -  character, gives phase type (P or S)
c
c     output 
c               
c               as funtion : distance in degrees
c
c     author:   J.Schweitzer, NORSAR, March 1999
c
c        June 2004: bug in model building of source layer corrected 
c                   subroutine EFAD included
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function dirdel(p0,zo0,fa,type)

      IMPLICIT real*8 (A-H,O-Z)

      PARAMETER (nz=24)

      character   type*1

      dimension v0(2,nz), v(nz), v2(nz), g(nz), z(nz), z0(nz), h(nz)

      data z0/  0. , 20. , 20. , 35. , 35. , 77.5,120. ,165. ,210. ,
     +        210. ,260. ,310. ,360. ,410. ,410. ,460. ,510. ,560. ,
     +        610. ,660. ,660. ,710. ,760. ,809.5/
      data v0/ 5.8   , 5.8   , 6.5   , 6.5   , 8.04  , 8.045 , 8.05  ,
     +         8.175 , 8.3   , 8.3   , 8.4825, 8.665 , 8.8475, 9.03  , 
     +         9.36  , 9.528 , 9.696 , 9.864 ,10.032 ,10.2   ,10.79  ,
     +        10.9229,11.0558,11.1353, 
     +         3.46  , 3.46  , 3.85  , 3.85  , 4.48  , 4.49  , 4.5   ,
     +         4.509 , 4.518 , 4.523 , 4.609 , 4.696 , 4.783 , 4.87  ,
     +         5.08  , 5.186 , 5.292 , 5.398 , 5.504 , 5.61  , 5.96  , 
     +         6.0897, 6.2095, 6.2426/

      PI  = 4.d0*DATAN(1.d0)
      PIM = PI/180.d0
      RE  = 6371.d0
      AA  = PIM*RE

      p = dabs(p0)

      zo = zo0

      dirdel = 0.0d0

      if(type.eq.'P') iph = 1
      if(type.eq.'S') iph = 2

      vmax = 0.d0

      DO 500 I=1,nz

      i1 = i + 1

      call efad(z0(i),v0(iph,i),z(i),v(i))

      if(vmax.lt.v(i)) vmax = v(i)

      if(zo.le.z0(i1)) then
         if(zo.eq.z0(i)) then
           call efad(z0(i1),v0(iph,i1),z(i1),v(i1))
            nl    = i1
            if(vmax.lt.v(i1)) vmax = v(i1)
            go to 550
         else if (zo.gt.z(i)) then
            dz = (zo-Z0(i))/(Z0(i1)-Z0(i))
            vdz = v0(iph,i) + dz * (V0(iph,i1)-v0(iph,i))
           call efad(zo,vdz,z(i1),v(i1))
            if(vmax.lt.v(i1)) vmax = v(i1)
            nl    = i1
            go to 550
         endif
      endif

500   continue

550   do 600 i = 1,nl

      I2=I+1
      H(I)=Z(I2)-Z(I)
      IF(H(I).EQ.0.d0)  H(I)=0.000001d0
      V2(I)=V(I2)
      IF(V2(I).EQ.V(I)) THEN
         V2(I)=1.000001d0*V(I)
         V(I2)=V2(I)
      ENDIF
600   G(I)=(V2(I)-V(I))/H(I)
 
      RVV = p / ((re-zo)*PIM)
      if(rvv*vmax.gt.1.d0) rvv=1.d0/vmax

      R=0.D0

      DO  1100  KK=1,NL-1
      E=V(KK)
      G1=E*RVV
      P=DSQRT(DABS(1.D0-G1*G1))
      O=1.d0/G(KK)
      F=V2(KK)
      G1=F*RVV
      Q=DSQRT(DABS(1.D0-G1*G1))
      R=R+FA*(P-Q)*O
1100  CONTINUE

      dirdel = r/(aa*rvv)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     function radloc
c
c     calculates the 'local' Earth radius using the
c     Earth's geometry and the local (geocentric) latitude.
c
c     input:
c
c     dlat0   =  local latitude
c
c     ind = 1    local latitude given in geographic coordinates
c         = 2    local latitude given in geocentric coordinates
c
      function radloc(dlat0,ind)
      real*8 dlat0,radloc,rada,radb,radian,dpythag,dlat
      real*8 convlat
      integer ind

c
c     Earth geometry after Stacey (1992) 'Physics of the Earth'
c
      rada = 6378.136d0
      radb = 6356.751d0

      radian = datan(1.d0) / 45.d0

      if (ind.eq.1) then 
         dlat = convlat(dlat0,1)
      else if (ind.eq.2) then
         dlat = dlat0
      else
         print *,'ERROR in RADLOC, variable IND wrongly set!!'
         stop
      endif

      radloc = dpythag(rada*dcos(radian*dlat),radb*dsin(radian*dlat))

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     subroutine ellcal
c
c     Here with corrections regarding the local Earth radius at the 
c     source. ELAT = source latitude in geographic coordinates.
c
    
      subroutine ellcal (elat,ax1,ax2,fchi,elmax,elmin,eazi,earea,fmud)
      real*8 elat,ax1(2),ax2(2),aaxes,baxes,earea,
     +       f1,deg2rad,elkm,fchi,elmax,elmin,eazi,pi,fmud,f2
      real*8 radloc, dpythag
      
      pi = 4.d0*datan(1.d0)
      deg2rad = pi/180.d0

      f1 = dcos(deg2rad*elat)

      ax1(2) = ax1(2) * f1
      ax2(2) = ax2(2) * f1

      elkm  = deg2rad*radloc(elat,1)

      f2 = dpythag(ax1(1),ax1(2))
      aaxes = elkm*fchi*dpythag(f2,fmud)
      f2 = dpythag(ax2(1),ax2(2))
      baxes = elkm*fchi*dpythag(f2,fmud)

      earea = aaxes*baxes*pi

      if(aaxes.ge.baxes) then

         elmax = aaxes
         elmin = baxes
         eazi = datan2(ax1(2),ax1(1))/deg2rad

      else

         elmax = baxes
         elmin = aaxes
         eazi = datan2(ax2(2),ax2(1))/deg2rad

      endif

      if(eazi.lt.0.d0)   eazi = eazi + 360.d0
      if(eazi.gt.180.d0) eazi = eazi - 180.d0

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        function convlat(rlat,ind)
c
c       Funtion to convert geographic latitude into geocentric
c       and back.
c
c       ind = 1 geographic to geocentric latitude
c
c           = 2 geocentric to geographic latitude
c
c       Johannes Schweitzer, NORSAR, October 2, 1997
c
        implicit real*8 (a-h,o-z)
        integer ind

        pi      = 4.d0*datan(1.d0)
        deg2rad = pi / 180.d0
        rad2deg = 180.d0 / pi

c
c     Earth figure after Stacey (1992), Physcis of the Earth
c
        rada = 6378.136d0
        radb = 6356.751d0
        eps = q2(radb/rada)

        elat = rlat*deg2rad

        if(ind.eq.1) then

          convlat = rad2deg * datan(eps*dtan(elat))

        else if(ind.eq.2) then

          convlat = rad2deg * datan(dtan(elat)/eps)

        else 

          print *,'ERROR in CONVLAT, variable IND wrongly set!!'
          stop

        endif

        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine efad(zsp,vsp,zfl,vfl)
C
c     input : depth and velocity (zsp, vsp) in a spherical Earth model
c
c     output: depth and velocity (zfl, vfl) in an equivalent flat 
c             Earth model
c
C     author: Johannes Schweitzer, NORSAR
C             June 2004
C
      IMPLICIT real*8 (A-H,O-Z)

      re=6371.d0

      if(zsp.lt.re) then
         f = re/(re-zsp)
         vfl = vsp*f
         zfl = re*dlog(f)
      else
        print *,'Earth-Flattening Approximation is not defined'
        print *,'for depths equal or larger than 6371 km !!! '
        print *,' Check model input!!!'
        stop
      endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     This subroutine calculates the azimuthal gap and the secondary 
c     azimuthal gap of a list of event observations
c
c     taken out from main program hyposat.f
c
      subroutine azigap(epiaz0,dazgap,d1azi,d2azi,dazgap2,d1azi2,d2azi2,
     +                  nobs,indx)

      real*4 dazgap,d1azi,d2azi,dazgap2,d1azi2,d2azi2
      real*4 epiaz0(*)
      integer nobs,indx(*)
      real*4 epiaz(4000)
c
c     internal
c
      integer i,j
      real*4 dmazi


      j = 1

      do 100 i = 1,nobs

         if (epiaz0(indx(i)).lt.0.) go to 100
         
         if (j.eq.1) then
            epiaz(j) = epiaz0(indx(i))
            j = j + 1
         else
            if(epiaz0(indx(i)).gt.epiaz(j-1)) then
              epiaz(j) = epiaz0(indx(i))
              j = j + 1
            else
              go to 100
            endif
         endif

100   continue

      j = j - 1

      dazgap = 0.
      d1azi  = 0.
      d2azi  = 360.

      if(j.gt.1) then

        do 200 i =1,j

           if(i.lt.j) then

             dmazi = epiaz(i+1) - epiaz(i)
             if(dmazi.gt.dazgap) then
                d1azi  = epiaz(i)
                d2azi  = epiaz(i+1)
                dazgap = dmazi
             endif

           else if(i.eq.j) then

             dmazi = epiaz(1) + 360. - epiaz(i)
             if(dmazi.gt.dazgap) then
                 d1azi  = epiaz(i)
                 d2azi  = epiaz(1)
                 dazgap = dmazi
             endif

           endif

200     continue

      else

        dazgap = 360.

      endif

      dazgap2 = 0.
      d1azi2  = 0.
      d2azi2  = 360.

      if(j.gt.2) then

         do 300 i =1,j

           if(i.lt.j-1) then

              dmazi = epiaz(i+2) - epiaz(i)
              if(dmazi.gt.dazgap2) then
                 d1azi2  = epiaz(i)
                 d2azi2  = epiaz(i+2)
                 dazgap2 = dmazi
              endif

           else if(i.eq.j-1) then

              dmazi = epiaz(1) + 360. - epiaz(i)
              if(dmazi.gt.dazgap2) then
                 d1azi2  = epiaz(i)
                 d2azi2  = epiaz(1)
                 dazgap2 = dmazi
              endif

           else if(i.eq.j) then

              dmazi = epiaz(2) + 360. - epiaz(i)
              if(dmazi.gt.dazgap2) then
                 d1azi2  = epiaz(i)
                 d2azi2  = epiaz(2)
                 dazgap2 = dmazi
              endif

           endif

300      continue

      else

         dazgap2 = 360.

      endif

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
