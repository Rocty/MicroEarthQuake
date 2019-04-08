      subroutine tauget_mod(zs,delta,n,phcd,ttc,dtdd,dtdh,dddp,modnam)
c
c     This routine calls IASP91-type tau-spline travel-time tables produced 
c     with the software known as libtau.f and libsun.f .
c
c     If this software is not available, it can be retrieved from 
c     anonymous ftp-server at the USGS, at IRIS, and from the RSES, ANU, 
c     Canberra.
c
c     latest changes 20 March 2017, JS, NORSAR
c
      save 

      include 'ttimes.h'

      real*4 zs, zso, delta

      integer n

      character*20 modnam, modnamo

      logical first
      common /bkin0/first,modnamo,zso

      if(n.lt.0) then
         zso     = -999.
         modnamo = ' '
         n = 0
      endif

      if (modnam.ne.modnamo .or. abs(zs-zso).gt.1.e-3) then

         in = 1 

         if (modnam.ne.modnamo) then
             call tabin(in,modnam)
             first = .true.
         endif

         call depset(zs)

         modnamo = modnam
         zso     = zs

      endif

      call trtm(delta,n,ttc,dtdd,dtdh,dddp,phcd)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tauget_ray(phase,phtyp,rayp0,modnam,depth,delray,ttray,
     +                      rayflag )
c
c     This routine calls IASP91-type tau-spline travel-time tables produced 
c     with the software known as libtau.f and libsun.f .
c
c     If this software is not available, it can be retrieved from 
c     anonymous ftp-server at the USGS, at IRIS, and from the RSES, ANU, 
c     Canberra.
c
c     latest changes 20 March 2017, JS, NORSAR
c
      save

      real*8 rayp0, delray, ttray, depth
      character*8 phase,phtyp*1, phase_type*1

      logical rayflag, rayok

      character*20 modnam,modnamo

      real*4  zs, zso
      logical first
      common /bkin0/first,modnamo,zso

      rayflag = .false.
      rayok      = .false.
      delray     = 0.d0
      ttray      = 0.d0
      tcor       = 0.
      xcor       = 0.
c

      if(phtyp.eq.' ') then
         if(phase.eq.' ' ) then
            if(rayp0.le.20.d0) then
               phtyp = 'P'
               phase = 'PKPdf'
            else
               phtyp = 'S'
               phase = 'S'
            endif
         else
            phtyp = phase_type(phase)
         endif
      endif

      zs   = sngl(depth)

      if (modnam.ne.modnamo .or. abs(zs-zso).le.1.e-3) then

         in = 1

         if(modnam.ne.modnamo) then
            call tabin(in,modnam)
            first = .true.
         endif
         call depset(zs)

         modnamo = modnam
         zso  = zs
      endif

      rayp = sngl(rayp0)

      itest = 0

100   call oneray(phase,rayp,xcor,tcor,rayok)

c     print *,phase,rayp,xcor,tcor,rayok

      if(rayok) then

         delray  = dble(xcor)
         rayflag = .true.

         go to 900

      else

         if(phtyp .eq. 'P') then

            if(rayp0.lt.5.d0 .and. itest.eq.0) then
               phase = 'PKPdf'
               itest = 1
               go to 100
            endif
     
            if(itest.eq.1) then

               if(phase.eq.'PKPdf') then
                 phase = 'PKPdif'
                 go to 100
               endif
     
               if(phase.eq.'PKPdif') then
                  phase = 'PKPbc'
                  go to 100
               endif
     
               if(phase.eq.'PKPbc') then
                  phase = 'PKPab'
                  go to 100
               endif
     
               if(phase.eq.'PKPab') then
                  phase = 'Pdif'
                  go to 100
               endif

               if(phase.eq.'Pdif') then
                  phase = 'PKiKP'
                  go to 100
               endif

            endif
     
            if(rayp0.le.21.d0 .and. itest.ge.0) then
               phase = 'P'
               itest = -1
               go to 100
            endif

            if(itest.lt.0) then

               if(phase.eq.'P') then
                  phase = 'Pn'
                  go to 100
               endif

               if(phase.eq.'Pn') then
                  phase = 'Pb'
                  go to 100
               endif

               if(phase.eq.'Pb') then
                  phase = 'Pg'
                  go to 100
               endif

            endif

         else if(phtyp .eq. 'S') then

            if(rayp0.lt.9.d0 .and. itest.eq.0) then
               phase = 'SKSdf'
               itest = 1
               go to 100
            endif

            if(itest.eq.1) then

               if(phase.eq.'SKSdf') then
                  phase = 'SKSdif'
                  go to 100
               endif

               if(phase.eq.'SKSdif') then
                  phase = 'SKSac'
                  go to 100
               endif

               if(phase.eq.'SKSac') then
                  phase = 'Sdif'
                  go to 100
               endif

            endif

            if(itest.eq.1) then
               phase = 'S'
               itest = -1
               go to 100
            endif

            if(itest.lt.0) then

               if(phase.eq.'S') then
                  phase = 'Sn'
                  go to 100
               endif

               if(phase.eq.'Sn') then
                  phase = 'Sb'
                  go to 100
               endif

               if(phase.eq.'Sb') then
                  phase = 'Sg'
                  go to 100
               endif

            endif

         endif

      endif

900   continue

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ellip(ecolatr,azi,del,zo,phas,p,ecor,ierr)
c
c****6789012345678901234567890123456789012345678901234567890123456789012
c
c     Subroutine ellip calls a routine to calculates the ellipticity 
c     correction for a given source-receiver combination for seismic 
c     phases of the AK135 tables (as far as availablbe). Approximations
c     iare used for several not defined phases.
c
c     After Brian Kennett (pers. communication) does this set of 
c     ellipticity corrections also work fine with IASP91 tables. 
c     Therefore, we will only use this set of tables in HYPOSAT for all 
c     different velocity models.
c
c     Johannes Schweitzer, NORSAR, February 2010
c
c     input:  ecolatr       geocentric colatitude of the event in rad
c
c             azi           azimuth from event to station in deg
c
c             del           distance between event and station in deg
c
c             zo            event depth in km
c
c             phas          phase name
c
c             p             ray paramter of phase in sec/deg
c
c     output: ecor          ellipticity correction of this phase in sec
c
c             ierr          error status
c
c
c version:  25. October 1996,  johannes schweitzer
c

      integer ierr
      real*4  ecolate,azi,p,ecor,p1,delo,zo,del,azi1,ecolatr
      character*8 phas, phas1
      logical errf

c
c     definition of several constants:
c
      deg2rad = atan(1.) / 45.
      ierr    = 0
      ecor    = 0.
      azi1    = azi
      delo    = del
      p1      = p
      ecolate = ecolatr

      phas1 = phas

c
c     Search for multiple core phases observed at distance 0 deg,
c     but not for .PKi... or .SKi... phases.
c
c     In this case, the phases travelled once around the Earth and 
c     ray parameter p1 (=0.!) is not an indication for this!
c
      if(delo.eq.0. .and. 
     *   (phas1(3:3).ne.'i'.and.phas1(4:4).ne.'i') ) then
      	 if (phas1(1:2).eq.'PK'.or.phas1(1:2).eq.'SK') p1=-999.
      	 if (phas1(1:2).eq."P'".or.phas1(1:2).eq."S'") p1=-999.
      	 if (phas1(2:3).eq.'PK'.or.phas1(2:3).eq.'SK') p1=-999.
      	 if (phas1(2:3).eq."P'".or.phas1(2:3).eq."S'") p1=-999.
      endif

      if(p1.lt.0.) then
         delo = 360. - delo
         azi1 = azi1 - 180.
         if(azi1.lt.0.) azi1 = 360. + azi1
      endif

      azi1    = deg2rad*azi1

      call elpcor(phas1,delo,zo,ecolate,azi1,ecor,errf)

      if(errf) then
         ierr = 999
         ecor = 0.
      endif

      return
      end
C=======================================================================
        SUBROUTINE elpcor(phase,edist,edepth,ecolat,azim,tcor,abrt)
c
c       SUBROUTINE ellip()
C                                                                         
C    Ellipticity correction for any given phase using
C    Dziewonski & Gilbert representation
C                                                   
C      The ellipticity corrections are found by linear interpolation       
C    in terms of values calculated for the ak135 model for a wide 
C    range of phases to match the output of the iasp software 
C
Cccj.s.     first call:  ellref(ecolat) 
Cccj.s.                        - to set up source dependent constants
C     2nd call  :  ellcor(phase,edist,depth,ecolat,azim,tcor,abrt) 
C                        - to calculate correction for a station
C
C    Parameters: 
C    character  
C          phase : a  string specifying the PHASE,   -e.g P, ScP etc.  
C                                                        
C    real*4
C          edist  :  epicentral distance to station (in degrees)     
C          edepth :  depth of event         
C          ecolat :  epicentral co-latitude of source (in radians) 
C          azim   :  azimuth from source to station (in radians)
C                                
C          tcor   :  time correction for path to allow for ellipticity
C 
C    logical 
C          abrt   :  a logical variable -usally set to .FALSE.  
C                    which is set to .TRUE. if a phase for      
C                    which no data is available is chosen       
C                                                                         
C=======================================================================
C   B.L.N. Kennett RSES,ANU        May 1995, August 1996                 
C   (based on earlier routine by D.J. Brown)
C   with input from W. Spakman, Utrecht
C=======================================================================
c
c    Slightly changed version:
c           input of data (path via environment variable HYPOSAT_DATA)
c           calling name elpcor
c           no initial call
c
c    October 1996 J. Schweitzer, Bochum
c
c       environment variable name corrected: Jan 27, 1997
c
c
      save sc0,sc1,sc2

      character *(*) phase
      character*8 phcod(57)
      integer phind(57),phspn(57),phnch(57)
      real*4 edist,edepth,ecolat,azim,
     ^       sc0,sc1,sc2,s3,tcor,
     ^       tau0, a0,b0,h0,d0,e0,f0,g0,
     ^       tau1, a1,b1,h1,d1,e1,f1,g1,
     ^       tau2, a2,b2,h2,d2,e2,f2,g2
      real*4 dpth(6),delta(50)
      real*4 t0(50,6),t1(50,6),t2(50,6)
      integer Ne,Nd
      logical abrt
c
c j.s.
c     
      character ic*120, file_check*120
      data phcod/
     & 'Pup   ','P     ','Pdif  ','PKPab ','PKPbc ','PKPdf ',
     & 'PKiKP ','pP    ','pPKPab','pPKPbc','pPKPdf','pPKiKP',
     & 'sP    ','sPKPab','sPKPbc','sPKPdf','sPKiKP','PcP   ',
     & 'ScP   ','SKPab ','SKPbc ','SKPdf ','SKiKP ','PKKPab',
     & 'PKKPbc','PKKPdf','SKKPab','SKKPbc','SKKPdf','PP    ',
     & "P'P'  ",'Sup   ','S     ','Sdif  ','SKSac ','SKSdf ',
     & 'pS    ','pSKSac','pSKSdf','sS    ','sSKSac','sSKSdf',
     & 'ScS   ','PcS   ','PKSab ','PKSbc ','PKSdf ','PKKSab',
     & 'PKKSbc','PKKSdf','SKKSac','SKKSdf','SS    ',"S'S'  ",
     & 'SP    ','PS    ','PnS   '/
      data phind/
     &        1,      14,      91,     136,     165,     178,
     &      235,     364,     433,     462,     475,     532,
     &      661,     742,     771,     784,     841,     970,
     &     1047,    1100,    1113,    1134,    1195,    1316,
     &     1337,    1382,    1507,    1516,    1573,    1702,
     &     1827,    1932,    1945,    2022,    2067,    2132,
     &     2197,    2234,    2295,    2356,    2425,    2490,
     &     2551,    2628,    2681,    2694,    2711,    2772,
     &     2781,    2838,    2967,    3140,    3273,    3398,
     &     3587,    3656,    3697/
      data phspn/
     &        3,      19,      11,       7,       3,      14,
     &       32,      17,       7,       3,      14,      32,
     &       20,       7,       3,      14,      32,      19,
     &       13,       3,       5,      15,      30,       5,
     &       11,      31,       2,      14,      32,      31,
     &       26,       3,      19,      11,      16,      16,
     &        9,      15,      15,      17,      16,      15,
     &       19,      13,       3,       4,      15,       2,
     &       14,      32,      43,      33,      31,      47,
     &       17,      10,       6/ 
      data phnch/
     &        3,       1,       4,       5,       5,       5,
     &        5,       2,       6,       6,       6,       6,
     &        2,       6,       6,       6,       6,       3,
     &        3,       5,       5,       5,       5,       6,
     &        6,       6,       6,       6,       6,       2,
     &        4,       3,       1,       4,       5,       5,
     &        2,       6,       6,       2,       6,       6,
     &        3,       3,       5,       5,       5,       6,
     &        6,       6,       6,       6,       2,       4,
     &        2,       2,       3/ 
      data dpth/ 0.0, 100.0, 200.0, 300.0, 500.0, 700.0 /
c...
c     In addition to the phase names listed above a number of phase 
c     aliases are available in the routine phase_alias, e.g. Pn --> 
c     P etc. The input phase code is first checked against the phcod 
C     array and next against the phase aliases.
c<sc>
c	           initial call to set up source dependent constants
cj.s. entry ellref(ecolat)
c                                            
cj.s. s3 = sqrt(3.0)/2.0
      s3 = sngl(dsqrt(0.75d0))
      sc0 = 0.25*(1.0+3.0*cos(2.0*ecolat))
      sc1 = s3*sin(2.0*ecolat)
      sc2 = s3*sin(ecolat)*sin(ecolat)
c     return
c<sc>
c<ec>                                           phase identification
cj.s. entry ellcor(phase,edist,edepth,ecolat,azim,tcor,abrt)
c      print *, 'phase,edist,edepth,ecolat,azim'
c      print *,  phase,edist,edepth,ecolat,azim
      Nd = 6
      NUMPH = 57
      deldst = 5.0
      abrt = .FALSE.
c                                       check on the length of phase
      l=len(phase)
      if(l.lt.8) then
       stop 
     >    'character variable `phase` should have at least length 8'
      endif

c                                             select phase
      ip = -1
c
c j.s.      nc=min(lnblk(phase),8)
c     to reduce code, use only one function to count characters
c
      nc=min(len_trim(phase),8)
      do 10 i=1,NUMPH
        if(nc.ne.phnch(i)) goto 10
        if (phase(1:nc) .eq. phcod(i)(1:nc)) then
          ip = i
          go to 11
        endif
 10   continue
 11   continue

      if(ip.eq.-1) then
c                                             check phase aliases
        call phase_alias(phase,edist,ip)
      endif
c                                              phase not found
c      print *, 'ip:',ip
      if(ip.lt.0) then
c       print *, phase,'  is not available'
        abrt = .true.
        return
      endif
      Ne = phspn(ip)
c                                          special case of upgoing waves
*
c                                              acquire phase information
c
c     changes for Bochum Oct 23, 1996 J.S.
c
c          open ellipticity correction file saved in varibale 'ic'
c

       ic = file_check('elcordir.tbl')
 
       open(35,file=trim(ic),access='direct',form='formatted',recl=80) 

       nr = phind(ip)
c       print*, 'nrec:',nr
       read(35,61,rec=nr) phcod(ip),np,d1,d2
c       print*, 'phcode,np,d1,d2: ', phcod(ip),np,d1,d2
       nr = nr+1
       if(np.ne.Ne) print*, 'HELP! - index wrong'
       do 15 i=1,np
         read(35,62,rec=nr) delta(i)
         nr = nr+1
         read(35,63,rec=nr) (t0(i,m),m=1,6)
         nr = nr+1
         read(35,63,rec=nr) (t1(i,m),m=1,6)
         nr = nr+1
         read(35,63,rec=nr) (t2(i,m),m=1,6)
         nr = nr+1
 15    continue         
 61    format(a8,i10,2f10.0)
 62    format(f10.0)
 63    format(6f10.4)
c 
       close (35)
c

c                                  distance index
       idist = 1 + int( (edist-d1)/ deldst )
       if(edist.lt.d1) idist =1
       if(edist.gt.d2) idist= np-1
c                                  depth index
       do 25 j = 1,Nd-1
         if ((dpth(j).le.edepth).and.(dpth(j+1).ge.edepth))then
            jdepth = j
            goto 26
         endif
 25    continue
 26    continue
*       print *, 'idist, jdepth;',idist,jdepth
c
*                      need to allow for zero entries (where phase
*                      description strongly depth dependent)
c tau0
         a0 = t0(idist,jdepth)
         b0 = t0(idist,jdepth+1)
         h0 = t0(idist+1,jdepth+1)
         d0 = t0(idist+1,jdepth)
         e0 = a0 + 
     ^       (d0-a0)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         f0 = b0 + 
     ^       (h0-b0)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         g0 = e0 + (f0-e0)*
     ^             (edepth-dpth(jdepth))/(dpth(jdepth+1)-dpth(jdepth))
         tau0 = g0
c tau1
         a1 = t1(idist,jdepth)
         b1 = t1(idist,jdepth+1)
         h1 = t1(idist+1,jdepth+1)
         d1 = t1(idist+1,jdepth)
         e1 = a1 + 
     ^       (d1-a1)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         f1 = b1 + 
     ^       (h1-b1)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         g1 = e1 + (f1-e1)*
     ^             (edepth-dpth(jdepth))/(dpth(jdepth+1)-dpth(jdepth))
         tau1 = g1
c tau2
         a2 = t2(idist,jdepth)
         b2 = t2(idist,jdepth+1)
         h2 = t2(idist+1,jdepth+1)
         d2 = t2(idist+1,jdepth)
         e2 = a2 + 
     ^       (d2-a2)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         f2 = b2 + 
     ^       (h2-b2)*(edist-delta(idist))/(delta(idist+1)-delta(idist))
         g2 = e2 + (f2-e2)*
     ^             (edepth-dpth(jdepth))/(dpth(jdepth+1)-dpth(jdepth))
         tau2 = g2
c
c         print *, 'tau0,tau1,tau2:',tau0,tau1,tau2
c j.s.   caz = cos(azim)
c j.s.   cbz = cos(2.0*azim)
c         print *, 'azim,caz,cbz',azim,caz,cbz    
c
         tcor = sc0*tau0 + sc1*cos(azim)*tau1 + sc2*cos(2.0*azim)*tau2
c
      return
c<ec>
      end
      subroutine phase_alias(phase,delta,ip)

c     check for alternative phase names
c     to get ellipticity corrections
c
c     input phase, delta
c     output ip (index of phcod)

      character*(*) phase

      if(phase(1:3).eq.'Pg ') then
c       phase='P       '
        ip=2
      else if(phase(1:3).eq.'Sg ') then
c       phase='S       '
        ip=33
      else if(phase(1:3).eq.'Lg ') then
c       phase='S       '
        ip=33
      else if(phase(1:4).eq.'pPg ') then
c       phase='pP      '
        ip=8
      else if(phase(1:4).eq.'sPg ') then
c       phase='sP      '
        ip=13
      else if(phase(1:4).eq.'pSg ') then
c       phase='pS      '
        ip=37
      else if(phase(1:4).eq.'sSg ') then
c       phase='sS      '
        ip=40
c
      elseif(phase(1:3).eq.'Pb ') then
c       phase='P       '
        ip=2
      else if(phase(1:3).eq.'Sb ') then
c       phase='S       '
        ip=33
      else if(phase(1:4).eq.'pPb ') then
c       phase='pP      '
        ip=8
      else if(phase(1:4).eq.'sPb ') then
c       phase='sP      '
        ip=13
      else if(phase(1:4).eq.'pSb ') then
c       phase='pS      '
        ip=37
      else if(phase(1:4).eq.'sSb ') then
c       phase='sS      '
c
      elseif(phase(1:3).eq.'Pn ') then
c       phase='P       '
        ip=2
      else if(phase(1:3).eq.'Sn ') then
c       phase='S       '
        ip=33
      else if(phase(1:4).eq.'pPn ') then
c       phase='pP      '
        ip=8
      else if(phase(1:4).eq.'sPn ') then
c       phase='sP      '
        ip=13
      else if(phase(1:4).eq.'pSn ') then
c       phase='pS      '
        ip=37
      else if(phase(1:4).eq.'sSn ') then
c       phase='sS      '
        ip=40
      else if(phase(1:4).eq.'SPn ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'SPb ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'SPg ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'SnP ') then
c       phase='SP      '
        ip=55
      else if(phase(1:4).eq.'PSn ') then
c       phase='PS      '
        ip=56
      else if(phase(1:5).eq.'PnPn ') then
c       phase='PP      '
        ip=30
      else if(phase(1:5).eq.'SnSn ') then
c       phase='SS      '
        ip=53
      elseif(phase(1:4).eq.'PbP ') then
c       phase='P       '
        ip=2
      elseif(phase(1:4).eq.'PmP ') then
c       phase='P       '
        ip=2
      elseif(phase(1:4).eq.'SbS ') then
c       phase='S       '
        ip=33
      elseif(phase(1:4).eq.'SmS ') then
c       phase='S       '
        ip=33
      elseif(phase(1:4).eq.'PbS ') then
c       phase='PS      '
        ip=56
      elseif(phase(1:4).eq.'PmS ') then
c       phase='PS      '
        ip=56
      elseif(phase(1:4).eq.'SbP ') then
c       phase='SP      '
        ip=55
      elseif(phase(1:4).eq.'SmP ') then
c       phase='SP      '
        ip=55
c                                       upgoing P, S
      else if(phase(1:2).eq.'p ') then
c       phase='Pup     '
        ip=1  
      else if(phase(1:2).eq.'s ') then
c       phase='Sup     '
        ip=32 
c                                        
      else if(delta.le.100.0.and.phase.eq.'pPdif   ') then
c       phase='pP      '
        ip=8
      else if(delta.le.100.0.and.phase.eq.'sPdif   ') then
c       phase='sP      '
        ip=13
      else if(delta.le.100.0.and.phase.eq.'pSdif   ') then
c       phase='pS      '
        ip=37
      else if(delta.le.100.0.and.phase.eq.'sSdif   ') then
c       phase='sS      '
        ip=40
      else if(delta.le.165.0.and.phase.eq.'PKPdif  ') then
c       phase='PKPbc '
        ip=5
      else if(delta.le.165.0.and.phase.eq.'pPKPdif ') then
c       phase='pPKPbc '
        ip=10
      else if(delta.le.165.0.and.phase.eq.'sPKPdif ') then
c       phase='sPKPbc '
        ip=15
c                             
      else if(phase(1:6).eq."P'P'P'") then
c       phase="P'P'P'  "
        ip =-1
c                             
      else if(phase(1:4).eq."P'P'") then
c       phase="P'P'    "
        ip =31
      else if(phase(1:6).eq."S'S'S'") then
c       phase="S'S'S'  "
        ip =-1
      else if(phase(1:4).eq."S'S'") then
c       phase="S'S'    "
        ip =54
c                            diffractions (approx)
      else if(delta.gt.100.0.and.phase.eq.'pPdif   ') then
c       phase='Pdif    '
        ip=3
      else if(delta.gt.100.0.and.phase.eq.'sPdif   ') then
c       phase='Pdif    '
        ip=3
      else if(delta.gt.100.0.and.phase.eq.'pSdif   ') then
c       phase='Sdif    '
        ip=34
      else if(delta.gt.100.0.and.phase.eq.'sSdif   ') then
c       phase='Sdif     '
        ip=34
c
      else
        ip=-1
      endif
      return
      end
