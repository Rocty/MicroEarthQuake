c----------------------------------------------------------------------
c  
c          Johannes Schweitzer
c          NORSAR
c          P.O.Box 53
c          N-2027 KJELLER
c          Norway
c
c  e-mail: johannes.schweitzer@norsar.no
c
c----------------------------------------------------------------------
c
c
      program HYPOSAT_6_0d

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      character  version*25
      parameter (version='HYPOSAT Version 6.0d     ')

c
c     last changes: 9 October 2018
c
c----------------------------------------------------------------------
c
c     Short desciption - for more details see HYPOSAT-Manual and
c     HYPOSAT papers in NORSAR Sci. Rep., PAGEOPH and PEPI.
c
c
c     This program locates seismic events by inverting
c     observed travel times, backazimuths, and slowness values.
c
c     Different phases observed at one station can be used to
c     calculate travel-time differences. These differences are
c     then used to invert for the hypocenter.
c
c     A preliminary epicenter will be defined with the backazimuth
c     observations, or with other available information.
c
c     If possible a source time will be estimated by a 
c     Wadati-Approach from the S-P travel-time difference(s) assuming 
c     a constant v(p)/v(s)=sqrt(3.) for each phase type separately.
c
c     The final location is done with a Single-Value-Decomposition
c     algorithm, which results in a least squares fit for
c     all four source parameters using travel-time models from
c     tau-spline-type tables (i.e., IASP91, AK135, PREM, ...),
c     and / or a local/regional model of horizontal layers.
c
c     All travel times are corrected  for the ellipticity of the
c     Earth.
c
c     All available information can be used including standard
c     deviations for observed data and 'a priori' information
c     for the model parameter.
c
c     All calculations are done for the Earth as a sphere and
c     travel times are corrected  for the ellipticity of the
c     Earth.
c
c--------------------------------------------------------------------
c
c               Program History
c
c         see file hyposat_history!
c
c--------------------------------------------------------------------
c
c            Main Program dependencies:
c
c     calls:     delazd, depi, dlsq, fetoh, fhtoe, findrange, get_station,
c                hyposat_cross, ellip, hyposat_gmi, indexx, 
c                tauget_mod, testphase, hyposat_geo, ttloc, get_mod_c10,
c                ellcal, plane, dpythag, isf_out_line, mult_ons,
c                tauget_ray, get_mod_global
c
c     functions: alpha1, alpha2, convlat, phase_type, crustc, 
c                ddmax, dirdel, dmean, q2, radloc
c                file_checkpara, read_event_id, read_origin, 
c                read_origin_head, read_phase_head, read_phase, 
c                write_origin, write_origin_head, lowcas, uppcas,
c                write_isf_error
c
c
c     data exchange by common blocks in include files:
c                
c                include 'gmi.h'
c                include 'lsq.h'
c                include 'ttimes.h'
c                include 'model.h'
c                include 'modelc.h'
c                include 'modelg.h'
c
c     PARAMETER settings for: mstat, mread, mvar, mosci0
c
c****6789012345678901234567890123456789012345678901234567890123456789012
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c     Functions called: variable definitions
c
      real*8           alpha1, alpha2, convlat, crustc, dirdel, q2, 
     +                 radloc, dmean, ddmax, getchi, dpythag
      character        phase_type*1, file_check*120, file_checkpara*120,
     +                 filepara*120, lowcas*30, uppcas*30, chgcas*30

c
c     mstat = maximum number of stations
c
      parameter  (mstat = 2000)

      dimension stala(mstat),stalo(mstat),stael(mstat),del(mstat),
     +          azie(mstat),baz(mstat),stalae(mstat),stavp(mstat),
     +          stavs(mstat),istaph(mstat),stats(mstat),statp(mstat),
     +          delk(mstat), istad(mstat),statr(mstat)

      character*5 sta(mstat),stat,stato,statw,stat1,stationfile*120,
     +          statcorfile*120,outputfile*120,inputfile*120,
     +          magfile*120,magmlfile*120,outputisf*120,inputfilen*120
c
c     mread = maximum number of phases (hyposat-in file)
c
c     mrd2 = maximum number of observations per station
c
      parameter  (mread = 4000, mr2=mread/2, mrd2 = 100)

      parameter  (mloc  = (mread/2 + 3)*mread)

      character phase(mread)*8,phaseu(mread)*8,phid*8,used(mread)*6,
     +          phid2*8,text2(mrd2)*155,phid1*8,phase_t*1,phidr0*8,
     +          string*155,touse(mread)*9,touse0*9,phidr*8,phipl*8,
     +          phai*8,phaj*8,o_string*155,textout*155,text(mread)*155,
     +          arid(mread)*8,statcorstr*80, texth*155

      dimension azi(mread),tt(mread),p(mread),azis(mread),tts(mread),
     +          ps(mread),period(mread),amplit(mread),dinv(mread,4),
     +          tt2(mread),ttu(mread),iev(mread),indx(mread),
     +          indx2(mrd2),idtu(mread),snr(mread),emeran(mread)

ccc  dimension  rays(mread,6)

      real*4    arr(mread),epiaz(mread),epiaz2(mread),dazgap,d1azi,
     +          d2azi,dazgap2,d1azi2,d2azi2

      dimension elo(mloc),ela(mloc),elos(mloc),elas(mloc)
c
c     include file for common block and variable definition GMI 
c
      include 'gmi.h'

      dimension var2(mvar)

c
c     include file for common block and variable definition LSQ 
c
      include 'lsq.h'

c
c     mosci0 is the maximum number of consequent solutions checked for
c     oscillating results.
c
      parameter (mosci0=15)

      dimension dzoos(mosci0),dtos(mosci0),dlaos(mosci0),
     +          dlo1os(mosci0),dloos(mosci0),dlo2os(mosci0),
     +          rzos(mosci0),rtos(mosci0),rlaos(mosci0),rloos(mosci0)

c
c     variables for travel-time calculations
c
      include 'ttimes.h'

      dimension dpdh(mphas), dtdd1(mphas), dtdd2(mphas), ttc1(mphas),
     +          dtdh1(mphas), dddp1(mphas)

      character phcd1(mphas)*8,phcd2(mphas)*8,art*16,modnam*20,
     +          modnam2*20,modn*20,modnam3*20,modnam4*20, mtyp0*3

      real*4 rzo,rdel,ecor,fla1,razi,pa,rzo1,rzo2,
     +       rzoe,rdel1,rmcorr, rdelk

      dimension imodn(4)

      logical first/.true./
      real*4  zso/0./
      character modnamo*20/' '/
      common /bkin0/first,modnamo,zso

c
c     common blocks for local/regional models reside in the following 
c     three include files:
c
      include 'model.h'
      include 'modelc.h'
      include 'modelg.h'

      character imod*1
c
c     variables used for ISF data handling
c
      integer   idum1, yyi, moni, ddi, ievid, itest, isf_null

      parameter (ISF_NULL=9999999)

      character cdum*1, cdum2*20, author*10, onscha*1,cevid*9, cdum3*2,
     +          phisf*10, isf_ref*10, phidd*10, author2*10, dformat*6,
     +          dsformat*5, corid*8, corid2*8, cpick*1

      real*4    rdum, rpa, ramp, rper, rsnr, rlati, rloni, rdepi,
     +          rdmi, rdma, rdum0, relmax, relmin

c
c     Functions called from the ISF software libarary
c
      integer   read_event_id, read_origin, read_origin_head, 
     +          read_phase_head, read_phase, write_origin, 
     +          write_origin_head, write_data_type, write_event_id,
     +          write_comment, write_netmag_head, write_netmag,
     +          write_phase_head, write_isf_error

c
c     other variables
c
      integer   yy,mon,dd,hh,idoy,ierr,typctl,mi,idum, isreg, regnum,
     +          typctlm, y00, mon00, d00, h00

      character mm*4,name*48
      real*8    lat,lon,kmdel,dlati,dloni,ddepi
      real*4    elevs, sec, rlat, rlon, smag

      character title*140, czo*1, region*80, czo1*1, magtypp*3, 
     +          magtyps*6, magtypml*7, statmag*2

c
c     idtmax = number of different travel-time-difference definitions
c              used calculating a initial value for the source time by 
c              using the Wadati-approach
c
      parameter (idtmax = 4)
      dimension dt(idtmax,mr2),idtp(idtmax,mr2),idts(idtmax,mr2),
     +          idt(idtmax),to(idtmax),tos(idtmax),vpvs(idtmax),
     +          vpvss(idtmax),datho(mread),datla(mread),datlo(mread),
     +          ttt(mread)

      logical   zoflag , vlflag  , stcorfl, iloc, surf,
     +          diffflag, dtmflag, epistart, single , output, modout, 
     +          last, magflag, mod2flag, lastfixi, direct, plflag,
     +          conr, rayok, rayokf, mod3flag, mod4flag, 
     +          aziini, kmout, thbaz, thray, tresw, lastfixt, lastfixm,
     +          mttres, isf_in, isf_out, fixinp, lsmu, ref_eve, 
     +          pflag, lgflag, sgflag, wadati, check_mult, aziflag,
     +          sloflag, gapobs, isf_epi, isf_dep, 
     +          new_in, natural, lgsurf, rgsurf, tsurf, lpsurf, isurf,
     +          lrmsisc, firston, firstph, firsts, azionlyf, ldepth0,
     +          old_syntax, emerout, larid

c
c     some constants and initial or default values
c

      pi      = 4.d0*datan(1.d0)
      deg2rad = pi / 180.d0
      rad2deg = 180.d0 / pi
c
c     Earth figure after Stacey (1992), Physcis of the Earth
c
      rada    = 6378.136d0
      radb    = 6356.751d0
      rearth  = 6371.d0

      eps     = q2(radb/rada)
      grad1   = rad2deg/rearth

      dtp0    = 600.d0
      dts0    = dtp0*2.d0
      dtm2    = dts0
      dtmin   = 9999.d0
      dtdt    = 10.d0

      dismin  = pi

      ttray   = 0.d0
      ddel    = 0.d0

      dchang0 = 1.d0

      check   = 9999.d0
      disper  = 0.001d0
      rmso    = 9999.d0
      rmsold  = 9999.d0
      datmax0 = 9999.d0
      nrms1   = 0

      miteras = 0
      iterz   = 0
      iteraz  = 0
      itso    = 0
      ibad0   = 0
      idelw   = 0
      in0sw   = 0
      infind  = 0

      insar   = 0

      nextiter1 = 0
      imaxiter  = 0
      ilastiter = 0
      mosci     = 4
      
      cevid  = '999999999'
      corid  = '_'
      corid2 = ' '
      ievid  = 0
      AUTHOR = 'HYPOSAT'

      vlflag   = .true.
      zoflag   = .false.
      stcorfl  = .false.
      diffflag = .true.
      iloc     = .false.
      epistart = .false.
      single   = .false.
      output   = .false.
      isf_in   = .false.
      isf_out  = .false.
      isf_epi  = .false.
      isf_dep  = .false.
      iellip   = .true.
      modout   = .false.
      magflag  = .false.
      mod2flag = .false.
      mod3flag = .false.
      mod4flag = .false.
      lastfixi = .false.
      lastfixt = .false.
      lastfixm = .false.
      plflag   = .true.
      dtmflag  = .false.
      rayok    = .false.
      rayokf   = .false.
      aziini   = .false.
      kmout    = .false.
      thray    = .false.
      thbaz    = .false.
      tresw    = .false.
      mttres   = .true.
      fixinp   = .false.
      locgeo   = .false.
      locsta   = .false.
      ref_eve  = .false.
      pflag    = .false.
      lgflag   = .false.
      sgflag   = .false.
      lgsurf   = .true.
      rgsurf   = .true.
      lpsurf   = .false.
      tsurf    = .false.
      isurf    = .false.
      check_mult = .false.
      new_in   = .false.
      aziflag  = .true.
      azionlyf = .false.
      sloflag  = .true.
      gapobs  = .false.
      firstph = .false.
      emerout = .false.

      ldepth0 = .false.

      lrmsisc = .false.

      modnam = 'ak135_A'
      modnam2 = modnam
      modnam3 = modnam
      modnam4 = modnam
      imodn(1) = 1
      imodn(2) = 0
      imodn(3) = 0
      imodn(4) = 0
      mtyp0  = 'E6 '
      rmax = 0.d0
      filloc = '_'
      imo     = 0
      stationfile = 'stations.dat'
      outputfile  = 'hyposat-out'
      inputfile   = 'hyposat-in'
      old_syntax = .false.
      statcorfile = ' '
      vpl = 5.8d0
      vsl = 3.46d0
      zo1 =  0.d0
      sdzo1 = 50.d0
      czo = 'F'
      depthmin = 0.d0
      depthmax = 800.d0
      maxiter = 80
      confl = 68.26895d0
      iellipi = 1
      dazim = 45.d0
      dpam  = 5.d0
      typctl = 4
      islow = 1
      setcheck1 = 0.d0
      setcheck  = 1.d0
      thrfixi0 = 0.005d0
      indph0 = 3333
      epilat0 = -999.d0
      epilon0 = -999.d0
      sdlatg  = 10.d0
      sdlat   = 0.d0
      sdlon   = 10.d0
      tome0   = -2840140801.d0
      stome0  = 120.d0
      wadmin  = 0.d0
      wadmax  = 300.d0
      dismaxst = 21001.d0
      disminst = -1.d0

      elmax  = -999.d0
      elmin  = -999.d0
      ieazi  = -999

      dtphmax = 5.d0

      dtmaxazib = 30.d0
      dtmaxazil = 180.d0
      dtmaxslow = 15.d0

      resmaxp = 30.d0
      resmaxs = 30.d0

      sglgdis = -999.d0

      smpu = 0.01d0
      smsu = 0.01d0
      lsmu  =  .false.
      var2(1) = 0.d0
      var2(2) = 0.d0
      var2(3) = 0.d0
      var2(4) = 0.d0

      string   = ' '
      o_string = ' '

      rdmi = 25000.
      rdma = 0.

      sisfi    = 0.1d0
      sisfe    = 1.0d0
      sisfo    = 2.0d0
      sisfaz   = 10.d0
      sisfsl   = 1.d0
      isf_ref  = ' '
      disfmod  = 999.d0

      vrg0    = 2.5d0
      vrg     = vrg0

      vlg0    = 3.5d0
      vlg     = vlg0

      vlr0    = 3.95d0
      vlr     = vlr0

      vlq0    = 4.4d0
      vlq     = vlq0

      vt0     = 1.45d0
      vt      = vt0

      vi0     = 0.33d0
      vi      = vi0

      magtypp  = ' '
      magtyps  = ' '
      magtypml = ' '
      magmlfile = 'MLCORR.TABLE'

      treswf   = 1.d0

c
c     search file 'hyposat-parameter'
c

      filepara = file_checkpara()
      if(filepara.eq.' ') then
         go to 9999
      else
         open (unit=9,file=filepara,status='old')
      endif

c
c     read in steering parameters from parameter file (if found)
c

      do 1 jin = 1,1000

      read (9,'(a)',end=2) string

      if(string(1:1).eq.' ') go to 1
      if(string(1:1).eq.'*') go to 1
      if(string(1:1).eq.'?') go to 1
      if(string(36:37).ne.': ') then
         print *,' Wrong syntax (ignored): ',trim(string)
         go to 1
      endif

      if(string(1:14).eq.'GLOBAL MODEL  ') then
          read (string(38:),'(a)') modnam
          go to 1
      endif

      if(string(1:14).eq.'GLOBAL MODEL 1') then
          read (string(38:),'(a)') modnam
          go to 1
      endif

      if(string(1:14).eq.'GLOBAL MODEL 2') then
          mod2flag = .false.
          read (string(38:),'(a)') modnam2
          if(modnam2 .ne. '_' .and. modnam2.ne.' ') then
            mod2flag = .true.
            imodn(2) = 1
          else
            modnam2 = modnam
          endif
          go to 1
      endif

      if(string(1:14).eq.'GLOBAL MODEL 3') then
          mod3flag = .false.
          read (string(38:),'(a)') modnam3
          if(modnam3 .ne. '_' .and. modnam3.ne.' ') then
            mod3flag = .true.
            imodn(3) = 1
          else
            modnam3 = modnam
          endif
          go to 1
      endif

      if(string(1:14).eq.'GLOBAL MODEL 4') then
          mod4flag = .false.
          read (string(38:),'(a)') modnam4
          if(modnam4 .ne. '_' .and. modnam4.ne.' ') then
            mod4flag = .true.
            imodn(4) = 1
          else
            modnam4 = modnam
          endif
          go to 1
      endif

      if(string(1:25).eq.'GLOBAL CRUSTAL MODEL CODE') then
          read (string(38:),'(a)') mtyp0
          go to 1
      endif

      if(string(1:23).eq.'LOCAL OR REGIONAL MODEL') then
          read (string(38:),'(a)') filloc
          go to 1
      endif

      if(string(1:19).eq.'VERY LOCAL GEOMETRY') then
          intinp = 0
          locgeo = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) locgeo = .true.
          go to 1
      endif

      if(string(1:27).eq.'LOCAL STATION BELOW SURFACE') then
          intinp = 0
          locsta = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) locsta = .true.
          go to 1
      endif

      if(string(1:24).eq.'OUTPUT OF REGIONAL MODEL') then
          intinp = 0
          modout = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) modout = .true.
          go to 1
      endif

      if(string(1:10).eq.'CRUST 5.1 ') then
          print *, 'WARNING: CRUST 5.1 not longer supported,'
          print *, 'we try to use newer CRUST 1.0 instead!'
          string(1:10) = 'CRUST 1.0 '
      endif

      if(string(1:11).eq.'CRUST 1.0  ') then
          read (string(38:),*) imo
          if(imo.lt.0) imo = 0
          go to 1
      endif

      if(string(1:12).eq.'STATION FILE') then
          read (string(38:),'(a)') stationfile
          stationfile = file_check(stationfile)
          go to 1
      endif

      if(string(1:19).eq.'STATION CORRECTIONS') then
          intinp = 0
          vlflag = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) vlflag = .true.
          go to 1
      endif

      if(string(1:23).eq.'STATION CORRECTION FILE') then
          stcorfl = .false.
          read (string(38:),'(a)') statcorfile
          if (trim(statcorfile).ne.'' .and.
     +        trim(statcorfile).ne.'_'   ) then
              statcorfile = file_check(statcorfile)
              stcorfl = .true.
          endif
          go to 1
      endif

      if(string(1:27).eq.'STATION CORR ONLY 1ST PHASE') then
          intinp = 0
          firstph = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) firstph = .true.
          go to 1
      endif

      if(string(1:31).eq.'P-VELOCITY TO CORRECT ELEVATION') then
          read (string(38:),*) vpl
          if(vpl.gt.99.d0 .or. vpl.lt.1.d-3) vpl=5.8d0
          go to 1
      endif

      if(string(1:31).eq.'S-VELOCITY TO CORRECT ELEVATION') then
          read (string(38:),*) vsl
          if(vsl.gt.99.d0) vsl=3.46d0
          if(vsl.lt.1.d-3) then
             if (vsl.gt.-1.d-3 .and. (vpl.ge.2.d-3 .and. vpl.le.99.d0))
     +            then
                vsl = vpl / dsqrt(3.d0)
             else
                vsl=3.46d0
             endif
          endif
          go to 1
      endif

      if(string(1:17).eq.'PLANE WAVE APPROX') then
          intinp = 0
          read (string(38:),*) intinp
          if(intinp.ne.1) plflag = .false.
          go to 1
      endif

      if(string(1:20).eq.'STARTING SOURCE TIME') then
c
c     time formats supported: epochal time 
c          or ASCII formated: yyyy-doy:hh.mi.ss.sss
c                             yyyy-mm-dd:hh.mi.ss.sss
c
c          seconds may be omitted
c
          if(string(38:38) .eq. '_') go to 1
          if(string(38:38) .eq. '0') go to 1
          if(string(38:39) .eq. '0.') go to 1
          
          if(string(42:42).ne.'-') then
             read (string(38:),*) timein
             if (timein.gt.tome0) tome0 = timein
          else

            mm   = ' '
            idum = 0
            idoy = 0
            mon  = 0
            dd   = 0
            hh   = 0
            mi   = 0
            sec  = 0.

            if(string(45:45).ne.'-') then
               if(len_trim(string).lt.51) then
                  print *,'Check format for source time in',
     +                    'hyposat-parameter file'
                  go to 9999
               endif
              read (string(38:51),'(i4,x,i3,2(x,i2))') yy,idoy,hh,mi
              if(len_trim(string).gt.52) read (string(53:),*) sec
            else
              if(len_trim(string).lt.53) then
                  print *,'Check format for source time in',
     +                    'hyposat-parameter file'
                 go to 9999
              endif
              read (string(38:53),'(i4,4(x,i2))') yy,mon,dd,hh,mi
              if(len_trim(string).gt.54) read (string(55:),*) sec
            endif

            call fhtoe(tome0,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

          endif
          go to 1
      endif

      if(string(1:19).eq.'STARTING TIME ERROR') then
          read (string(38:),*) stome0
          go to 1
      endif

      if(string(1:23).eq.'DOUBLE SIDED TIME ERROR') then
          tresw = .false.
          read (string(38:),*) itresw
          if(itresw.eq.1 .or. itresw.eq.2) then
             tresw  = .true.
          endif
          go to 1
      endif

      if(string(1:22).eq.'DBLE SID. ERROR FACTOR') then
          read (string(38:),*) treswf
          if(treswf.le.0.d0) treswf = 1.d0
          go to 1
      endif

      if(string(1:21).eq.'STARTING SOURCE DEPTH') then
          read (string(38:),*) zo1
          go to 1
      endif

      if(string(1:20).eq.'STARTING DEPTH ERROR') then
          read (string(38:),*) sdzo1
          if(sdzo1.eq.0.d0) sdzo1 = 50.d0
          go to 1
      endif

      if(string(1:10).eq.'DEPTH FLAG') then
          read (string(38:),'(a)') czo
          go to 1
      endif

      if(string(1:20).eq.'DEPTH ALLOWED ABOVE 0') then
          intinp = 0
          ldepth0   = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) ldepth0   = .true.
          go to 1
      endif

      if(string(1:13).eq.'MINIMUM DEPTH') then
          read (string(38:),*) depthmin
          go to 1
      endif

      if(string(1:13).eq.'MAXIMUM DEPTH') then
          read (string(38:),*) depthmax
          if (depthmax.gt.700.d0) depthmax = 700.d0
          go to 1
      endif

      if(string(1:24).eq.'STARTING SOURCE LATITUDE') then
          abc = -999.0d0
          read (string(38:),*) abc
          if(abc.ge.-90.d0 .and. abc.le.90.d0) epilat0 = abc
          go to 1
      endif

      if(string(1:23).eq.'STARTING LATITUDE ERROR') then
          read (string(38:),*) sdlatg
          go to 1
      endif

      if(string(1:25).eq.'STARTING SOURCE LONGITUDE') then
          abc = -999.0d0
          read (string(38:),*) abc
          if(abc.ge.-180.d0 .and. abc.le.180.d0) epilon0 = abc
          go to 1
      endif

      if(string(1:24).eq.'STARTING LONGITUDE ERROR') then
          read (string(38:),*) sdlon
          go to 1
      endif

      if(string(1:27).eq.'INCLUDING MODEL UNCERTAINTY') then
          intinp = 0
          lsmu   = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) lsmu   = .true.
          go to 1
      endif

      if(string(1:29).eq.'MEAN P-WAVE MODEL UNCERTAINTY') then
          abc = 0.d0
          read (string(38:),*) abc
          if(abc.gt.0.d0) smpu = abc
          go to 1
      endif

      if(string(1:29).eq.'MEAN S-WAVE MODEL UNCERTAINTY') then
          abc = 0.d0
          read (string(38:),*) abc
          if(abc.gt.0.d0) smsu = abc
          go to 1
      endif

      if(string(1:17).eq.'MIN DT FOR WADATI') then
          read (string(38:),*) wadmin
          if(wadmin.lt.0.d0) wadmin = 0.d0
          go to 1
      endif

      if(string(1:17).eq.'MAX DT FOR WADATI') then
          read (string(38:),*) wadmax
          go to 1
      endif

      if(string(1:20).eq.'MAX EPI DIST OF STAT') then
          read (string(38:),*) dismaxst
          go to 1
      endif

      if(string(1:20).eq.'MIN EPI DIST OF STAT') then
          read (string(38:),*) disminst
          go to 1
      endif

      if(string(1:23).eq.'MAXIMUM # OF ITERATIONS') then
          intinp = 0
          read (string(38:),*) intinp
          if(intinp.gt.0) maxiter = intinp
          go to 1
      endif

      if(string(1:24).eq.'# TO SEARCH OSCILLATIONS') then
          intinp = 0
          read (string(38:),*) intinp
          if(intinp .gt. mosci0) then
              print *, 'MAX # TO SEARCH OSCILLATIONS SET TO:',mosci0
              mosci = mosci0
          endif
          if(intinp .lt. 1) then
             mosci = 1
          else
             mosci = intinp
          endif
          go to 1
      endif

      if(string(1:16).eq.'CONFIDENCE LEVEL') then
          read (string(38:),*) confl
          if(confl.lt.68.26895d0)  confl = 68.26895d0
          if(confl.gt.99.99d0)     confl = 99.99d0
          go to 1
      endif

      if(string(1:18).eq.'CONSTRAIN SOLUTION') then
          intinp = 0
          read (string(38:),*) intinp
          if(intinp.eq.1) then
             lastfixt = .true.
             lastfixi = .false.
             lastfixm = .false.
          else if(intinp.eq.2) then
             lastfixt = .false.
             lastfixi = .true.
             lastfixm = .false.
          else if(intinp.eq.3) then
             lastfixt = .true.
             lastfixi = .true.
             lastfixm = .false.
          else if(intinp.eq.4) then
             lastfixt = .true.
             lastfixi = .false.
             lastfixm = .true.
          else if(intinp.eq.5) then
             lastfixt = .true.
             lastfixi = .true.
             lastfixm = .true.
          else
             lastfixt = .false.
             lastfixi = .false.
             lastfixm = .false.
          endif
          go to 1
      endif

      if(string(1:26).eq.'MAXIMUM ALLOWED P RESIDUUM') then
          abc = -9.d0
          read (string(38:),*) abc
          if (abc.gt.0.d0) then
             resmaxp = abc
          else
             resmaxp = 30.d0
          endif
          go to 1
      endif

      if(string(1:26).eq.'MAXIMUM ALLOWED S RESIDUUM') then
          abc = -9.d0
          read (string(38:),*) abc
          if (abc.gt.0.d0) then
             resmaxs = abc
          else
             resmaxs = 30.d0
          endif
          go to 1
      endif

      if(string(1:20).eq.'MEAN T-T RES. CORREC') then
          intinp = 0
          mttres = .false.
          read (string(38:),*) intinp
          if (intinp.eq.1) mttres = .true.
          go to 1
      endif


      if(string(1:29).eq.'INF. DENSITY MATRIX THRESHOLD') then
          read (string(38:),*) thrfixi0
          if (thrfixi0.gt.1.d0) thrfixi0=1.d0
          if (thrfixi0.le.0.d0) then
             thrfixi0  = -1.d0
          endif
          go to 1
      endif

      if(string(1:23).eq.'EPICENTER ERROR ELLIPSE') then
          iellip = .true.
          read (string(38:),*) iellipi
          if (iellipi.ne.1) iellip = .false.
          go to 1
      endif

      if(string(1:30).eq.'AZIMUTHAL GAP FOR OBSERVATIONS') then
          intinp = 0
          gapobs = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) gapobs = .true.
          go to 1
      endif

      if(string(1:21).eq.'MAXIMUM AZIMUTH ERROR') then
          read (string(38:),*) dazim
          go to 1
      endif

      if(string(1:12).eq.'AZIMUTH ONLY') then
          intinp = 0
          azionlyf = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1)  azionlyf = .true.
          go to 1
      endif

      if(string(1:19).eq.'AZIMUTH AS DEFINING') then
          intinp = 0
          aziflag = .true.
          read (string(38:),*) intinp
          if(intinp.ne.1) aziflag = .false.
          go to 1
      endif

      if(string(1:26).eq.'MAX T RES FOR AZI OF B USE') then
          read (string(38:),*) dtmaxazib 
          go to 1
      endif

      if(string(1:26).eq.'MAX T RES FOR AZI OF L USE') then
          read (string(38:),*) dtmaxazil
          go to 1
      endif

      if(string(1:21).eq.'AZIMUTH ONLY INIT SOL') then
          intinp = 0
          aziini = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) aziini = .true.
          go to 1
      endif

      if(string(1:20).eq.'SLOWNESS AS DEFINING') then
          intinp = 0
          sloflag = .true.
          read (string(38:),*) intinp
          if(intinp.ne.1) sloflag = .false.
          go to 1
      endif

      if(string(1:22).eq.'MAXIMUM SLOWNESS ERROR') then
          read (string(38:),*) dpam
          go to 1
      endif

      if(string(1:26).eq.'MAX T RES FOR SLOWNESS USE') then
          read (string(38:),*) dtmaxslow
          go to 1
      endif

      if(string(1:12).eq.'OUTPUT LEVEL') then
          read (string(38:),*) typctl
          typctlm = typctl
          if (typctl.lt.0)   typctl = 0
          if (typctl.ge.40)  typctl = 4
          if (typctl.gt.10) then
             itypn = mod(typctl,10)
             typctl = 4
             if(itypn.ge.5) typctl = itypn
          endif
          go to 1
      endif

      if(string(1:27).eq.'FLAG EMERGENCE ANGLE OUTPUT') then
          intinp = 0
          emerout = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) emerout = .true.
          go to 1
      endif

      if(string(1:16).eq.'SLOWNESS [S/DEG]') then
          intinp = 0
          read (string(38:),*) intinp
          if(intinp.eq.0 .or. intinp.eq.1) islow = intinp
          if(isf_in) islow = 1
          go to 1
      endif

      if(string(1:17).eq.'LOCATION ACCURACY') then
          read (string(38:),*) setcheck1
          go to 1
      endif

      if(string(1:27).eq.'PHASE INDEX FOR LOCAL MODEL') then
          read (string(38:),*) indph0
          if(indph0.gt.3333) then
            print *,'Wrong input: PHASE INDEX FOR LOCAL MODEL'
            go to 9999
          endif
          go to 1
      endif

      if(string(1:34).eq.'FLAG USING TRAVEL-TIME DIFFERENCES') then
          intinp = 0
          diffflag = .true.
          read (string(38:),*) intinp
          if(intinp.ne.1) diffflag = .false.
          go to 1
      endif

      if(string(1:22).eq.'FLAG USING INPUT FIXED') then
          intinp = 0
          fixinp = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) fixinp = .true.
          go to 1
      endif

      if(string(1:29).eq.'FLAG CHECKING MULTIPLE ONSETS') then
          intinp = 0
          check_mult = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) check_mult = .true.
          go to 1
      endif

      if(string(1:19).eq.'FLAG NEW INPUT FILE') then
          intinp = 0
          new_in = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) new_in = .true.
          go to 1
      endif

      if(string(1:26).eq.'MAX DT FOR MULTIPLE ONSETS') then
          read (string(38:),*) dtphmax
          if(dtphmax.le.0.d0) dtphmax = 5.0d0
          go to 1
      endif

      if(string(1:15).eq.'INPUT FILE NAME') then
          if(string(38:38).ne.' ' .and. string(38:38) .ne.'_')
     +                        read (string(38:),*) inputfile
          go to 1
      endif

      if(string(1:21).eq.'HYPOSAT-IN OLD SYNTAX') then
          intinp = 0
          old_syntax = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) old_syntax = .true.
          go to 1
      endif

      if(string(1:16).eq.'INPUT FORMAT ISF') then
          intinp = 0
          isf_in = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) then
             isf_in = .true.
             islow = 1
          endif
          go to 1
      endif

      if(string(1:22).eq.'ISF REFERENCE LOCATION') then
          read (string(38:),'(a)') isf_ref 
          if (isf_ref.eq.'_') isf_ref = ' '
          go to 1
      endif

      if(string(1:13).eq.'ISF EPICENTER') then
          intinp = 0
          isf_epi = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) isf_epi=.true.
          go to 1
      endif

      if(string(1:9).eq.'ISF DEPTH') then
          intinp = 0
          isf_dep =  .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) isf_dep=.true.
          go to 1
      endif

      if(string(1:5).eq.'ISF_i') then
          read (string(38:),*) sisfi 
          if(sisfi.le.0.d0) sisfi = 0.1d0
          go to 1
      endif

      if(string(1:5).eq.'ISF_e') then
          read (string(38:),*) sisfe
          if(sisfe.le.0.d0) sisfe = 1.0d0
          go to 1
      endif

      if(string(1:5).eq.'ISF_o') then
          read (string(38:),*) sisfo
          if(sisfo.le.0.d0) sisfo = 2.d0
          go to 1
      endif

      if(string(1:6).eq.'ISF_az') then
          read (string(38:),*) sisfaz
          if(sisfaz.le.0.d0) sisfaz = 10.d0
          go to 1
      endif

      if(string(1:7).eq.'ISF_slo') then
          read (string(38:),*) sisfsl
          if(sisfsl.le.0.d0) sisfsl = 1.d0
          go to 1
      endif

      if(string(1:22).eq.'ISF_2ND MODEL DISTANCE') then
          abc = -999.d0
          read (string(38:),*) abc
          if(abc.gt.0.d0 .and. abc.le.180.d0) disfmod = abc
          go to 1
      endif

      if(string(1:12).eq.'ISF EVENT ID') then
          read (string(38:),*) ievid
          if(ievid.ge.999999999 .or. ievid.lt.0) then
             ievid = 0
          else
             write(cevid,'(i9)') ievid
          endif
          go to 1
      endif

      if(string(1:13).eq.'ISF ORIGIN ID') then
          read (string(38:),*) corid
          go to 1
      endif

      if(string(1:17).eq.'OUTPUT FORMAT ISF') then
          intinp = 0
          isf_out = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) isf_out=.true.
          go to 1
      endif

      if(string(1:16).eq.'ISC-TYPE ISF RMS') then
          intinp = 0
          lrmsisc = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) lrmsisc=.true.
          go to 1
      endif

      if(string(1:16).eq.'OUTPUT FILE NAME') then
          if(string(38:38).ne.' ' .and. string(38:38) .ne.'_')
     +                        read (string(38:),*) outputfile
          go to 1
      endif

      if(string(1:13).eq.'OUTPUT SWITCH') then
          intinp = 0
          output = .true.
          read (string(38:),*) intinp
          if(intinp.ne.1) output=.false.
          go to 1
      endif

      if(string(1:12).eq.'OUTPUT IN KM') then
          intinp = 0
          kmout = .false.
          read (string(38:),*) intinp
          if (intinp.eq.1) kmout = .true.
          go to 1
      endif

      if(string(1:21).eq.'OUTPUT OF THEO. BAZ+P') then
          intinp = 0
          thbaz = .false.
          thray = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) thbaz = .true.
          if(intinp.eq.2) thray = .true.
          if(intinp.eq.3) then
             thbaz = .true.
             thray = .true.
          endif
          go to 1
      endif

      if(string(1:18).eq.'AUTHOR OF SOLUTION') then
          read (string(38:),'(a)') author
          go to 1
      endif

      if (string(1:25).eq.'REGIONAL SURFACE WAVES LG') then
          intinp = 0
          lgsurf = .true.
          read (string(38:),*) intinp
          if (intinp.ne.1) lgsurf = .false.
          go to 1
      endif

      if(string(1:17).eq.'LG GROUP-VELOCITY') then
         read (string(38:),*) vlg
         if(vlg.le.0.d0) vlg = vlg0
         go to 1
      endif

      if (string(1:25).eq.'REGIONAL SURFACE WAVES RG') then
          intinp = 0
          rgsurf = .true.
          read (string(38:),*) intinp
          if (intinp.ne.1) rgsurf = .false.
          go to 1
      endif

      if(string(1:17).eq.'RG GROUP-VELOCITY') then
          read (string(38:),*) vrg
          if(vrg.le.0.d0) vrg = vrg0
          go to 1
      endif

      if (string(1:24).eq.'LP SURFACE WAVES (LQ/LR)') then
          intinp = 0
          lpsurf = .false.
          read (string(38:),*) intinp
          if (intinp.eq.1) lpsurf = .true.
          go to 1
      endif

      if(string(1:17).eq.'LR GROUP-VELOCITY') then
          read (string(38:),*) vlr
          if(vlr.le.0.d0) vlr = vlr0
          go to 1
      endif

      if(string(1:17).eq.'LQ GROUP-VELOCITY') then
          read (string(38:),*) vlq
          if(vlq.le.0.d0) vlq = vlq0
          go to 1
      endif

      if (string(1:13).eq.'T-PHASE USAGE') then
          intinp = 0
          tsurf = .false.
          read (string(38:),*) intinp
          if (intinp.eq.1) tsurf = .true.
          go to 1
      endif

      if(string(1:17).eq.'T PHASE GROUP-VEL') then
          read (string(38:),*) vt
          if(vt.le.0.d0) vt = vt0
          go to 1
      endif

      if (string(1:14).eq.'IS-PHASE USAGE') then
          intinp = 0
          isurf = .false.
          read (string(38:),*) intinp
          if (intinp.eq.1) isurf = .true.
          go to 1
      endif

      if(string(1:18).eq.'IS PHASE GROUP-VEL') then
          read (string(38:),*) vi
          if(vi.le.0.d0) vi = vi0
          go to 1
      endif

      if (string(1:18).eq.'P-TYPE ONSETS ONLY') then
          intinp = 0
          pflag = .false.
          read (string(38:),*) intinp
          if (intinp.eq.1) pflag = .true.
          go to 1
      endif

      if (string(1:14).eq.'LG-PHASE TO SG') then
          intinp = 0
          lgflag = .false.
          read (string(38:),*) intinp
          if (intinp.eq.1) lgflag = .true.
          go to 1
      endif

      if (string(1:14).eq.'SG-PHASE TO LG') then
          intinp = 0
          sgflag = .false.
          read (string(38:),*) intinp
          if (intinp.eq.1) sgflag = .true.
          go to 1
      endif

      if (string(1:15).eq.'SG--LG DISTANCE') then
          abc = -999.d0
          read (string(38:),*) abc
          if(abc.gt.0.d0) sglgdis = abc
          go to 1
      endif
      
      if (string(1:21).eq.'MAGNITUDE CALCULATION') then
          intinp = 0
          magflag = .false.
          read (string(38:),*) intinp
          if (intinp.eq.1) magflag = .true.
          go to 1
      endif

      if (string(1:19).eq.'P-ATTENUATION MODEL') then
          read (string(38:),*) magtypp
          go to 1
      endif

      if (string(1:20).eq.'MS-ATTENUATION MODEL') then
          read (string(38:),*) magtyps
          go to 1
      endif

      if (string(1:20).eq.'ML-ATTENUATION MODEL') then
          read (string(38:),*) magtypml
          go to 1
      endif

      if (string(1:18).eq.'ML-CORRECTION FILE') then
          read (string(38:),*) magmlfile
          go to 1
      endif

      if(string(1:15).eq.'REFERENCE EVENT') then
          intinp = 0
          ref_eve = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) ref_eve = .true.
          go to 1
      endif

      if(string(1:26).eq.'REFERENCE SOURCE LONGITUDE') then
          abc = -999.0d0
          read (string(38:),*) abc
          if(abc.ge.-180.d0 .and. abc.le.180.d0) dloni = abc
          go to 1
      endif

      if(string(1:25).eq.'REFERENCE SOURCE LATITUDE') then
          abc = -999.0d0
          read (string(38:),*) abc
          if(abc.ge.-90.d0 .and. abc.le.90.d0) dlati = abc
          go to 1
      endif

      if(string(1:22).eq.'REFERENCE SOURCE DEPTH') then
          read (string(38:),*) ddepi
          go to 1
      endif

1     continue

2     close(9)

      fchi1 = dsqrt(getchi(1.0d0,confl))
      fchi2 = dsqrt(getchi(2.0d0,confl))

      if(zo1.lt.depthmin) zo1 = depthmin
      if(zo1.gt.depthmax) zo1 = depthmax

      if(typctl.gt.4) then
         print *,'modnam   : ',modnam
         if(mod2flag) print *,'modnam 2: ',modnam2
         if(mod3flag) print *,'modnam 3: ',modnam3
         if(mod4flag) print *,'modnam 4: ',modnam4
         print *,'filloc = ',filloc
         print *,'imo = ',imo
         print *,'stationfile = ',stationfile
         print *,'statcorfile = ',statcorfile 
         print *,'inputfile   = ',inputfile 
         print *,'outputfile  = ',outputfile 
         print *,'output switch ',output 
         print *,'output isf ',isf_out 
         print *,'vpl = ',vpl
         print *,'vsl = ',vsl
         print *,'vrg = ',vrg
         print *,'vlg = ',vlg
         print *,'vlq = ',vlq
         print *,'vlr = ',vlr
         print *,'vt  = ',vt 
         print *,'zo1   = ',zo1
         print *,'sdzo1 = ',sdzo1
         print *,'czo   = ',czo
         print *,'epilat0 = ',epilat0
         print *,'sdlatg  = ',sdlatg
         print *,'epilon0 = ',epilon0
         print *,'sdlon   = ',sdlon
         print *,'maxiter = ',maxiter
         print *,'lastfixi = ',lastfixi
         print *,'lastfixi threshold= ',thrfixi0
         print *,'lastfixt = ',lastfixt
         print *,'lastfixm = ',lastfixm
         print *,'confl   = ',confl  
         print *,'dazim = ',dazim 
         print *,'dpam  = ',dpam
         print *,'typctl = ',typctl
         print *,'islow = ',islow
         print *,'setcheck1 = ',setcheck1
         print *,'indph0 = ',indph0 
         print *,'diffflag = ',diffflag
         print *,'Magnitude flags = ', magflag,' ',magtypp,' ',magtyps
     +           ,' ',magtypml
         print *,'Plane wave = ', plflag
      endif

c
c     initializing gobal model
c

      rzo1  = 0.
      rdel1 = 105.
      nphas = -999
      call tauget_mod(rzo1,rdel1,nphas,phcd,ttc,dtdd,dtdh,dddp,modnam)
      pdif  = dble(dtdd(1))

c

      mtyp = mtyp0
      call get_mod_global(typctl,ierr)
      if (ierr.ne.0) then
         rmax = 0.d0
         write(*,'('' Cannot use crustal model corrections '',2a)') 
     +         mtyp,modnam
         if(imo.eq.2) imo = 1
         if(imo.eq.4) imo = 3
         ierr = 0
      else
c
c     save of Global Crustal Model in v0g / zg / ...
c
         jmodg = jmod
         elevg = elev
         zmaxg = zmax
         do 21 i = 1,jmodg
            v0g(1,i) = v0(1,i)
            v0g(2,i) = v0(2,i)
            zg(i)    = z(i)
21       continue
      endif

      if (imo.ge.3) then

         itrue = 1
         inum  = 1
         iread = 0
         elatc = 0.d0
         elonc = 0.d0
         mtyp = 'C10'
         rmax = 1.5d0
         iloc = .true.
         filloc = 'CRUST 1.0'

         call get_mod_c10(itrue,inum,typctl,ierr)

      else if(imo.le.2) then

         mtyp = mtyp0

         if(trim(filloc).ne.' ' .and. indph0.ge.0 .and.
     +      filloc(1:1).ne.'_'                ) then
         
            filloc = file_check(filloc)
            ierr   = 0
            rzo    = 0.
            rdel   = 0.
            czo1   = ' '
            indph  = 0
            elatc  = 0.d0
            elonc  = 0.d0
            sdep   = 0.d0
            rmax = 0.d0
            jmod = 0

            nphas = 0
            call ttloc(rzo,rdel,czo1,nphas,ttc,dtdd,dtdh,dpdh,dddp,
     +                 phcd,rmax,sdep,typctl,ierr,indph,emerout)

            if(ierr.ne.0) then 
               rmax = 0.d0
               write(*,'('' Can only use global model: '',a,
     +              '' ierr= '',i3)') modnam,ierr
               if(imo.eq.1) imo = 0
               ierr = 0
            else

               iloc = .true.

            endif
         endif
      endif

      if(ldepth0 .and. .not.iloc) then

         print *,'Source depth above sea level only possible for ',
     +           'receivers and source within local/regional ',
     +           'model distance!!!'
         go to 9999
      endif
      dazim1 = dazim
      dpam1  = dpam

      if (stcorfl) open(13,file=statcorfile)

      if(setcheck1.gt.0.d0) then
         setcheck  = setcheck1
      else
         setcheck1 = setcheck
      endif

      rminh    = setcheck
      if(rminh.lt.0.1d0) rminh = 0.1d0
      rmint    = rminh/5.d0
      rming    = rminh*grad1

      setcheck2 = 15.d0*setcheck

      chgcas = uppcas(czo)
      czo = chgcas(1:1)
      if(czo.ne.'F' .and. czo.ne.'D' .and. czo.ne.'B' ) czo='F'
      zo = zo1

      if(epilat0.ge.-90.d0 .and. epilon0.ge.-180.0d0 .and.
     +   epilat0.le.90.d0 .and. epilon0.le.180.0d0 ) epistart = .true.

      do 3 i=1,mstat
      sta(i) = ' '
3     continue

c
c     read in all available observed data
c

      open (unit=10,file=inputfile,status='old')

      title = version

      if(output) then
         open (unit=11,file=outputfile)
         write (11,'(a,/)') trim(title)
         write (11,'(''Event solution by '',a,/)') 
     +          trim(author)
      endif

      if (isf_out) then
         if(kmout) then
            isf_out = .false.
            print *, 'ISF output not possible: distance units are',
     +            ' in [km]!'
         endif
      endif

      print *,'PROGRAM ',trim(title)
      print *, ' '

31    read (10,'(a)',end=9999) title

      if(isf_in) then

         if(ievid.eq.0) then
            itest = read_event_id(title, cevid, cdum2) 
            if(itest.eq.20) go to 31

            if(cevid(9:9).eq.' ') then
               read(cevid(1:8),'(i8)') ievid
               cevid(2:9) = cevid(1:8)
               cevid(1:1) = ' '
            else
               read(cevid,'(i9)') ievid
            endif
         endif

         natural = .true.

         if(output) write (11,'(''ISF EVENT ID: '',i9,/)') ievid

32       read (10,'(a)',end=9999) string
         itest = read_origin_head(string)
         if(itest.eq.20) go to 32 

         author2 = ' '

33       o_string = string
         read (10,'(a)',end=9999) string

         if(isf_ref.eq.'PRIME') then

            if (string(1:9).eq.' (#PRIME)' .or.
     +          len_trim(string).le.1) then
               string = o_string
            else
               if(string(1:2).eq.' (') string = o_string
               go to 33
            endif

         endif

         if(string(1:2).eq.' (') then
            string = o_string
            go to 33
         endif

         itest = read_origin(string,yyi,moni,ddi,hh,mi,isec,msec,
     +   cdum,rdum,rdum,rlati,rloni,cdum,rdum,rdum,idum1,rdepi,
     +   cdum,rdum,idum1,idum1,idum1,rdum,rdum,cdum,cdum,
     +   cdum3,author2,corid2)


         if(itest.eq.20) then 
            itest = write_isf_error(1)
            go to 33
         endif

         if(author2.ne.isf_ref .and. isf_ref.ne.' ' .and.
     +      isf_ref.ne.'PRIME') go to 33

         if(cdum3(2:2).ne.'e' .and. cdum3(2:2).ne.'k' .and. 
     +      cdum3(2:2).ne.' ' .and. cdum3(2:2).ne.'n' ) then
            natural = .false.
c           if(output) then
c              write(11,'(a)') string
c              write(11,'(a)') 'Event type wrong: '
c              go to 9999
c           endif
c           go to 33
         endif

         if((corid.eq.' ' .or. corid.eq.'_') .and.
     +      (corid2.ne.' ' .and. corid2.ne.'_') ) corid = corid2

         dlati = dble(rlati)
         dloni = dble(rloni)
         if(rdepi.eq.real(ISF_NULL)) rdepi = 0.
         ddepi = dble(rdepi)

         sec = real(isec)
         if(msec.ne.isf_null) then
            sec = sec + real(msec)/1000.
         endif
         mm = ' '
         idoy = 0
         jdate = 0
         timeoi = 0.d0
         call fhtoe(timeoi,jdate,yyi,moni,mm,ddi,idoy,hh,mi,sec)

         if(isf_epi) then
            epilat0 = dlati
            epilon0 = dloni
            tome0   = timeoi
            epistart = .true.
         endif
         if(isf_dep) then
            zo1     = ddepi
            if(zo1.lt.depthmin) zo1=depthmin
            if(zo1.gt.depthmax) zo1=depthmax
            zo      = zo1
         endif

34       read (10,'(a)',end=9999) string
         itest = read_phase_head(string)
         if(itest.eq.20) go to 34

      else
         if(output) write (11,'(a,/)') trim(title)
         print *,'EVENT ',trim(title)
      endif

      timemin = 9999999999.d0

      stalam = 0.d0
      stalom = 0.d0
      stalom1 = 0.d0
      stalom2 = 0.d0

      terrm = 0.d0
      nobsst = 0

      isnr = 0

      sdpmean = 0.d0
      nsdp    = 0
      sdsmean = 0.d0
      nsds    = 0
      sdmeans = 0.d0

      ii = 0
      string = ' '

      do 12 i=1,mread+200

      read (10,'(a)',end=14) string

      if(string.eq.o_string) go to 12
      if(string(1:4).eq.'STOP') go to 14
      if(string(1:1).eq.'*') go to 12
      if(string(1:1).eq.' ') go to 12
      if(string.eq.' ') go to 12

      ii = ii + 1

      if(ii.gt.mread) then
        print *,'Maximum number of input data reached: ',mread
        go to 9999
      endif

      ierr = 0

4     continue

      azi(ii) = -999.d0
      pin     = -999.d0
      snr(ii) = -999.d0
      amplit(ii) = -999.d0
      period(ii) = -999.d0
      touse0  = 'TASDRM   '
      arid(ii)= ' '
      tt2(ii) = -1.d0

      if(.not.isf_in) then

          lstring = len_trim(string)

          if(lstring.le.34) then
             go to 5
          else if(lstring.eq.35) then
             read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f3.0)',err=5)
     +             stat,phase(ii),yy,mon,dd,hh,mi,sec
          else if(lstring.eq.36) then
             read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f4.1)',err=5)
     +             stat,phase(ii),yy,mon,dd,hh,mi,sec
          else if(lstring.eq.37) then
             read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f5.2)',err=5)
     +             stat,phase(ii),yy,mon,dd,hh,mi,sec
          else if(lstring.eq.38) then
             read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f6.3)',err=5)
     +             stat,phase(ii),yy,mon,dd,hh,mi,sec
          else if(lstring.ge.69) then
             read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f6.3,1x,
     +             f5.3,1x,f6.2,3(1x,f5.2))',err=5) 
     +             stat,phase(ii),yy,mon,dd,hh,mi,sec,
     +             tts(ii),azi(ii),azis(ii),pin,ps(ii)
          else if(lstring.ge.63) then
             read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f6.3,1x,
     +             f5.3,1x,f6.2,2(1x,f5.2))',err=5) 
     +             stat,phase(ii),yy,mon,dd,hh,mi,sec,
     +             tts(ii),azi(ii),azis(ii),pin
          else if(lstring.ge.57) then
             read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f6.3,1x,
     +             f5.3,1x,f6.2,f5.2)',err=5) 
     +             stat,phase(ii),yy,mon,dd,hh,mi,sec,
     +             tts(ii),azi(ii),azis(ii)
          else if(lstring.ge.51) then
             read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f6.3,1x,
     +             f5.3,1x,f6.2)',err=5) 
     +             stat,phase(ii),yy,mon,dd,hh,mi,sec,
     +             tts(ii),azi(ii)
          else if(lstring.ge.44) then
             read (string,'(a5,1x,a8,1x,i4,4(1x,i2),1x,f6.3,1x,
     +             f5.3)',err=5) 
     +             stat,phase(ii),yy,mon,dd,hh,mi,sec,tts(ii)
          else
             go to 5
          endif

          ipos = 71
          ipo2 = ipos + 6
          if (old_syntax) ipo2 = ipos + 5

          if(lstring.ge.ipos) then

             if (old_syntax) then
                read (string(ipos:ipo2),'(a6)',err=5) touse0
                touse0(7:7) = touse0(6:6)
                touse0(6:6) = 'M'
             else
                read (string(ipos:ipo2),'(a7)',err=5) touse0
             endif

             chgcas = uppcas(touse0)
             touse0 = chgcas(1:7) // '  '

          endif

          ipos = ipo2 + 2

          if(lstring.ge.ipos) then
             ipo2 = ipos + 5
             read (string(ipos:ipo2),'(f6.3)',err=5) abc
             if(abc.gt.0.d0) period(ii)= abc
          endif

          ipos = ipo2 + 2

          if(lstring.ge.ipos) then
             ipo2 = ipos + 11
             read (string(ipos:ipo2),'(f12.2)',err=5) abc
             if(abc.gt.0.d0) amplit(ii)= abc
          endif

          ipos = ipo2 + 2
          if(lstring.ge.ipos) then
             ipo2 = ipos + 6
             read (string(ipos:ipo2),'(f7.2)',err=5) abc
             if(abc.gt.0.d0) snr(ii)= abc
          endif

          ipos = ipo2 + 2
          if(lstring.ge.ipos) then
             ipo2 = ipos + 7
             read (string(ipos:ipo2),'(a8)',err=5) arid(ii)
          endif

          ipos = ipo2 + 2
          if(itresw.eq.1 .and. lstring.ge.ipos) then
             ipo2 = ipos + 4
             read (string(ipos:ipo2),'(f5.2)',err=5) abc
             if(abc.gt.0.d0) tt2(ii) = abc
          endif

          if (ii.eq.1) then
             y00   = yy
             mon00 = mon
             d00   = dd
             h00   = hh
          else
             if(string(16:19).eq.'    ') yy = y00
             if(string(21:22).eq.'  ')   mon = mon00
             if(string(24:25).eq.'  ')   dd = d00
             if(string(27:28).eq.'  ')   hh = h00
          endif

          ierr = 0
          go to 6 

5         continue

          if(ierr.le.1) then

             indcs= 200
             indcs = index(string,'#')
             if(indcs.le.120 .and. indcs.gt.0) then
                string(indcs:indcs)=' '
                go to 5
             endif
             ierr = ierr + 1
             go to 4

          else

             ii = ii - 1
             go to 12

          endif
      
      else

          if(string(21:23).eq.' DI') string(21:23) = '_DI'
          itest = read_phase (string,stat,rdum0,rdum,phisf,hh,mi,
     +    isec,msec,rdum,razi,rdum,rpa,rdum,touse0(1:1),touse0(2:2),
     +    touse0(3:3),rsnr,ramp,rper,cpick,cdum,onscha,cdum2,cdum,
     +    smag,arid(ii))

          if(itest.eq.20) then 
c            itest = write_isf_error(0)
             ii = ii - 1
             go to 12
          endif

          if(hh+mi+isec+msec.eq.4*isf_null) then
             ii = ii - 1
             go to 12
          endif

c         print *,stat,hh,mi,isec,msec,phisf,razi,rpa,touse(ii)

          if(dble(rdum0).ge.disfmod) touse0(7:7) = '2'

          chgcas = lowcas(onscha)
          onscha = chgcas(1:1)

          if(onscha.ne.'i' .and. onscha.ne.'e') onscha='_'

          sec = real(isec)
          if(msec.ne.isf_null) then
             sec = sec +real(msec)/1000.
          else
             msec = 0
          endif

          mm = ' '
          idoy = 0
          jdate = 0
          timeop = 0.d0
          call fhtoe(timeop,jdate,yyi,moni,mm,ddi,idoy,hh,mi,sec)

          if( (timeop.lt.timeoi) .and.
     +        (timeop+86400.d0-timeoi.lt.7200.d0) ) then
              timeop = timeop + 86400.d0
          endif
          
          idum = 0
          mm =' '
          call fetoh(timeop,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

          if(phisf(1:1).eq.'(') then
             phidd=phisf(2:)
             ip = len_trim(phidd)
             if(phidd(ip:ip).eq.')') then
                phisf=phidd(1:ip-1)
             else
                phisf=phidd
             endif
          endif

          if(uppcas(phisf(1:1)).eq.'I') then
             phidd=phisf(2:)
             ip = len_trim(phidd)
             if(ip.gt.1) then
                phisf=phidd(1:ip-1)
             else
                phisf='tx'
             endif
             onscha='i'
          endif

          if(uppcas(phisf(1:1)).eq.'E') then
             phidd=phisf(2:)
             ip = len_trim(phidd)
             if(ip.gt.1) then
                phisf=phidd(1:ip-1)
             else
                phisf='tx'
             endif
             onscha='e'
          endif

          if(phisf(1:2).eq.'P ') phisf='P1'
          if(uppcas(phisf(1:3)).eq.'PN ') phisf='P1'
          if(phisf(1:3).eq.'P* ') phisf='Pb'
          if(uppcas(phisf(1:3)).eq.'PG1') phisf='Pg '

          if(phisf(1:2).eq.'S ') phisf='S1'
          if(uppcas(phisf(1:3)).eq.'SN ') phisf='S1'
          if(phisf(1:3).eq.'(S)') phisf='S1'
          if(phisf(1:3).eq.'S* ') phisf='Sb'
          if(uppcas(phisf(1:3)).eq.'SG1') phisf='Sg '
          if(phisf(1:4).eq.'TSG1') phisf='Sg  '
          if(phisf(1:4).eq.'*ESG1') phisf='Sg   '

          if(phisf(1:9).eq.'(S)/(SKS)') phisf='S1'
          if(phisf(1:7).eq.'(S)/SKS') phisf='S1'
          if(phisf(1:7).eq.'S/(SKS)') phisf='S1'
          if(phisf(1:5).eq.'S/SKS') phisf='S1'

          if(phisf(1:3).eq.'SKS') phisf='S1'
          if(phisf(1:4).eq.'SKSa') phisf='S1'
          if(phisf(1:4).eq.'SKS2') phisf='SKSdf'
          if(phisf(1:4).eq.'sSKS') phisf='sSKSac'
          if(phisf(1:4).eq.'pSKS') phisf='pSKSac'

          if(phisf(1:4).eq.'SKKS') phisf='SKKSac'
          if(phisf(1:5).eq.'SKKS2') phisf='SKKSac'
          if(phisf(1:5).eq.'SKKKS') phisf='S3KSac'

          if(phisf(1:4).eq.'SKP ') phisf='SKPab'
          if(phisf(1:5).eq.'SKP2 ') phisf='SKPdf'
          if(phisf(1:5).eq.'SKKP ') phisf='SKKPab'
          if(phisf(1:6).eq.'SKKP2 ') phisf='SKKPdf'

          if(phisf(1:6).eq.'P_DIFF') phisf='Pdif'
          if(phisf(1:5).eq.'P_DIF') phisf='Pdif'

          if(phisf(1:6).eq.'S_DIFF') phisf='Sdif'
          if(phisf(1:5).eq.'S_DIF') phisf='Sdif'

          if(phisf(1:4).eq.'PKP ') phisf='PKPdf'
          if(phisf(1:5).eq.'pPKP ') phisf='pPKPdf'
          if(phisf(1:5).eq.'sPKP ') phisf='sPKPdf'
          if(phisf(1:5).eq.'PKP2 ') phisf='PKPab'
          if(phisf(1:6).eq.'pPKP2 ') phisf='pPKPab'
          if(phisf(1:6).eq.'sPKP2 ') phisf='sPKPab'

          if(phisf(1:4).eq.'PKS ') phisf='PKSdf'
          if(phisf(1:5).eq.'PKS2 ') phisf='PKSab'

          if(phisf(1:6).eq.'PKKP ') phisf='PKKPdf'
          if(phisf(1:6).eq.'PKKP2 ') phisf='PKKPbc'
          if(phisf(1:6).eq.'PKKP3 ') phisf='PKKPdf'

          if(phisf(1:6).eq.'PKKS ') phisf='PKKSdf'
          if(phisf(1:6).eq.'PKKS2 ') phisf='PKKSbc'
          if(phisf(1:6).eq.'PKKS3 ') phisf='PKKSab'

          if(uppcas(phisf(1:1)).ne.'P' .and.  
     +       uppcas(phisf(1:1)).ne.'S' .and. 
     +       uppcas(phisf(1:1)).ne.'L' .and. 
     +       uppcas(phisf(1:1)).ne.'R' .and. 
     +       uppcas(phisf(1:1)).ne.'T'   ) then

             ii = ii - 1
             go to 12

          endif

          phase(ii) = phisf(1:8)

          if(rpa.ne.isf_null) then
              pin = dble(rpa)
              touse0(3:3) = 'S'
              ps(ii) = sisfsl
          else
              pin    = -999.d0
          endif

          if(razi.ne.isf_null) then
              azi(ii)  = dble(razi)
              touse0(2:2) = 'A'
              azis(ii) = sisfaz
          else
              azi(ii)  = -999.d0
          endif

          if(rsnr.ne.isf_null) then
              snr(ii) = dble(rsnr)
          else
              snr(ii) = -999.d0
          endif

          if(ramp.ne.isf_null) then
              amplit(ii) = dble(ramp)
          else
              amplit(ii) = -999.d0
          endif

          if(rper.ne.isf_null) then
              period(ii) = dble(rper)
          else
              period(ii) = -999.d0
          endif

          if(smag.gt.0.) touse0(6:6) = 'M'

          tts(ii) = sisfo
          if(onscha.eq.'i') tts(ii) = sisfi
          if(onscha.eq.'e') tts(ii) = sisfe
          if(mod(msec,100).gt.0 .and. mod(msec+1,100).gt.1) 
     +           tts(ii) = 0.8d0 * tts(ii)

          if(tts(ii).lt.0.1d0) tts(ii)=0.1d0

          if(phase_type(phase(ii)).ne.'P') tts(ii) = tts(ii) * 2.d0

          touse0(1:1) = 'T'
          if(diffflag) touse0(4:4) = 'D'

      endif

6     mm = ' '
      idoy = 0
      jdate = 0
      timeo = 0.d0
c     print*, timeo,jdate,yy,mon,mm,dd,idoy,hh,mi,sec
      call fhtoe(timeo,jdate,yy,mon,mm,dd,idoy,hh,mi,sec)
      tt(ii) = timeo

      touse(ii)=   'TASDRM   '

      if(touse0.ne.'         ') then
         if(touse0(1:1).ne.'T') touse0(1:1)=' '
         if(touse0(2:2).ne.'A') touse0(2:2)=' '
         if(touse0(3:3).ne.'S') touse0(3:3)=' '
         if(touse0(4:4).ne.'D') touse0(4:4)=' '
         if(touse0(5:5).ne.'R') touse0(5:5)=' '
         if(touse0(6:6).ne.'M') touse0(6:6)=' '
         if(touse0(7:7).ne.' ' .and. touse0(7:7).ne.'1' .and. 
     +      touse0(7:7).ne.'2' .and. touse0(7:7).ne.'3' .and.
     +      touse0(7:7).ne.'4') touse0(7:7)=' '

         if(pflag.and.phase_type(phase(ii)).ne.'P') 
     +      touse0(1:6) = '      '
         touse(ii)=touse0

      endif

      if(phase(ii)(1:1).ne.'P' .and. phase(ii)(1:1).ne.'S' .and.
     +   phase(ii)(1:1).ne.'p' .and. phase(ii)(1:1).ne.'s' .and.
     +   phase(ii)(1:1).ne.'R' .and. phase(ii)(1:1).ne.'L' .and.
     +   phase(ii)(1:2).ne.'IS' )   touse(ii)(1:4) = ' A  '

      if(imo.ne.2 .and. imo.ne.4) touse0(5:5) = ' '

      if(.not.aziflag) touse(ii)(2:2) = ' '
      if(.not.sloflag) touse(ii)(3:3) = ' '

      if(touse(ii)(1:1).eq.'T' .or. touse(ii)(4:4).eq.'D') then
         if(tts(ii).le.0.d0) tts(ii) = 2.d0
         terrm  = terrm + tts(ii)
         nobsst = nobsst + 1
      endif

      if(tresw) then
          if(tt2(ii).le.0.d0) tt2(ii) = tts(ii)
          if(itresw.eq.2) then
              tt2(ii) = tts(ii)*treswf
          endif
      else
          tt2(ii) = tts(ii)
      endif

      if(phase(ii).eq.'PKP2')   phase(ii)="P'P'"
      if(phase(ii).eq.'PKPPKP') phase(ii)="P'P'"
      if(phase(ii).eq.'SKS2')   phase(ii)="S'S'"
      if(phase(ii).eq.'SKSSKS') phase(ii)="S'S'"
      if(phase(ii).eq.'PPP')    phase(ii)='P3'
      if(phase(ii).eq.'SSS')    phase(ii)='S3'
      if(phase(ii).eq.'PKhKP')  phase(ii)='PKPpre'

      incap = index(lowcas(phase(ii)),'diff')

      if(incap.ne.0) then 
         phidd = phase(ii)
         phase(ii)(incap:) = 'dif' // phidd(incap+4:)
      endif

62    phid = phase(ii)
      incap =  index(phid,'N')
      if(incap.ne.0) then
         phase(ii)(incap:incap) = 'n'
         go to 62
      endif
      incap =  index(phid,'B')
      if(incap.ne.0) then
         phase(ii)(incap:incap) = 'b'
         go to 62
      endif
      incap =  index(phid,'A')
      if(incap.ne.0) then
         phase(ii)(incap:incap) = 'a'
         go to 62
      endif
      incap =  index(phid,'G')
      if(incap.ne.0) then
         phase(ii)(incap:incap) = 'g'
         go to 62
      endif
      incap =  index(phid,'C')
      if(incap.ne.0) then
         phase(ii)(incap:incap) = 'c'
         go to 62
      endif
      incap =  index(phid,'M')
      if(incap.ne.0) then
         phase(ii)(incap:incap) = 'm'
         go to 62
      endif
      incap =  index(phid,'D')
      if(incap.ne.0) then
         phase(ii)(incap:incap) = 'd'
         go to 62
      endif
      incap =  index(phid,'F')
      if(incap.ne.0) then
         phase(ii)(incap:incap) = 'f'
         go to 62
      endif
      incap =  index(phid,'l')
      if(incap.eq.1) then
         phase(ii)(incap:incap) = 'L'
         go to 62
      endif

      if(lgflag.and.phase(ii).eq.'Lg'.and..not.sgflag) phase(ii)='Sg'
      if(sgflag.and.phase(ii).eq.'Sg'.and..not.lgflag) phase(ii)='Lg'

      if(phid(2:2).eq.' ') touse(ii)(5:5)=' '

      incap =  index(phid,'1')
      if(incap.ne.0 .and. phid(3:3).eq.' ') touse(ii)(5:5)=' '

      if(phid(4:4).eq.' ') then
         if(phid(2:2).eq.'b') touse(ii)(5:5)=' '
         if(phid(2:2).eq.'g') touse(ii)(5:5)=' '
         if(phid(2:2).eq.'n') touse(ii)(5:5)=' '
      endif

      if(phid(2:2).eq.'c') touse(ii)(5:5)=' '
      if(phid(2:2).eq.'d') touse(ii)(5:5)=' '
      if(phid(2:2).eq.'K') touse(ii)(5:5)=' '

      if(touse(ii)(1:1).ne.' ' .or. touse(ii)(4:4).ne.' ') then
         phase_t = phase_type(phid)
         if(phase_t.eq.'P') then
            nsdp    = nsdp + 1
            sdpmean = sdpmean + tts(ii)*tts(ii)
         endif
         if(phase_t.eq.'S') then
            nsds    = nsds + 1
            sdsmean = sdsmean + tts(ii)*tts(ii)
         endif
      endif

      if(string(46:51).eq.'      ') azi(ii) = -999.d0
      if(azi(ii).lt.0.d0) then
         azi(ii)       = -999.d0
         touse(ii)(2:2)= ' '
         azis(ii)      =    0.d0
      endif

c
c     The standard errors for azimuth or ray parameter are yet not 
c     given: set default values (30 or 40 [deg] and 5 [s/deg])!
c
      if(azi(ii).ge.0.d0 .and. azis(ii).le.0.d0) then

          azis(ii)= 30.d0

          if(phase(ii)(1:2).eq.'LR'.or.phase(ii)(1:2).eq.'LQ' )
     +             azis(ii)=40.d0

      endif

      chgcas = uppcas(stat)
      stat = chgcas(1:5)

      do 10 j=1,mstat

      if(stat.eq.sta(j)) then

        iev(ii) = j
        go to 11

      else if (sta(j).eq.' ') then

        call get_station(stationfile,stat,jdate,lat,lon,
     +                   elevs,name,ierr)

        if(ierr.ne.0) then
          print *,'Cannot find station: ',stat,' entry skipped'
          ii = ii - 1
          ierr = 0
          go to 12
        endif

        statpc = 0.d0
        statsc = 0.d0
        statrc = 0.d0

        if(vlflag .and. stcorfl) then
          rewind(13)

7         statcorstr = ' '
          read (13,*,end=81,err=81) statcorstr
          if(statcorstr(1:1).eq.'*') go to 7

          read(statcorstr,*,err=801,end=801) stat1,vpc,vsc,spc,ssc,src
          go to 69
801       src = 0d0
          read(statcorstr,*,err=802,end=802) stat1,vpc,vsc,spc,ssc
          go to 69
802       spc = 0.d0
          ssc = 0.d0
          read(statcorstr,*,err=81,end=81) stat1,vpc,vsc

69        if(stat1.eq.stat) then
            vp = vpc
            vs = vsc
            if(vp.le.0.d0) vp = vpl
            if(vs.le.0.d0) vs = vp / dsqrt(3.d0)
            statpc = spc
            statsc = ssc
            statrc = src
            go to 9
          endif
          go to 7
        endif

81      if(vlflag) then
          vp = vpl
          vs = vsl
        endif

9       continue
c       print *,vlflag,vp,vs

c       if (typctl.gt.8) then
c          print *,j,stat,lat,lon,elevs,name,vp,vs,statpc,statsc
c       endif

        sta(j)   = stat
        stala(j) = lat
        stalae(j)= convlat(lat,1)
        stalo(j) = lon
        stael(j) = dble(elevs)/1000.d0
        stavp(j) = vp
        stavs(j) = vs
        statp(j) = statpc
        stats(j) = statsc
        statr(j) = statrc
        istaph(j) = 0
        istad(j)  = 0
        stalam = stalam + stalae(j)
        p1 = deg2rad*lon
        stalom1 = stalom1 + dcos(p1)
        stalom2 = stalom2 + dsin(p1)
        iev(ii) = j
        istat = j

        go to 11

      endif

10    continue

11    if(timeo.lt.timemin) then
         timemin = timeo
         istatmin = iev(ii)
      endif

      if(string(59:63).eq.'      ') pin = -999.d0
      if(pin.le.0.d0)               pin = -999.d0

      if(islow.eq.0 .and. pin.gt.0.0d0) then
         p(ii)  = radloc(stala(iev(ii)),1)*deg2rad/pin
         ps(ii) = ps(ii)*p(ii)/pin
      else
         p(ii) = pin
      endif

      if(p(ii).le.0.d0) then
         touse(ii)(3:3) = ' '
         ps(ii) =  0.d0
      else
         if(ps(ii).le.0.d0) ps(ii)= 5.d0
      endif

      chgcas = uppcas(phase(ii)(1:1))
      phase_t = chgcas(1:1)

      if(phase_t.eq.'P')  then
         if (istaph(iev(ii)).eq.0 .or. istaph(iev(ii)).eq.2) 
     +       istaph(iev(ii))=istaph(iev(ii)) + 1
      else if(phase_t.eq.'S') then
         if (istaph(iev(ii)).eq.0 .or. istaph(iev(ii)).eq.1)
     +       istaph(iev(ii))=istaph(iev(ii)) + 2
      endif

      if(phase(ii)(1:1).eq. 'P' .and. phase(ii)(3:3).eq. ' ')
     +             touse(ii)(8:8) = 'P'
      if(phase(ii)(1:3).eq. 'PKP' .and. phase(ii)(6:6).eq. ' ') 
     +             touse(ii)(8:8) = 'P'
      if(phase(ii)(1:4).eq. 'Pdif') touse(ii)(7:7) = 'P'

      if(typctl.gt.8) then
         print *,ii,stat,tt(ii),phase(ii),azi(ii),azis(ii),p(ii),ps(ii),
     +           touse(ii)
      endif

12    o_string = string


14    close(10)

      if(stcorfl) close (13)

      if(nobsst.gt.0) then
         terrm = terrm / dble(nobsst)
      else
         terrm = 2.d0
      endif
      terrm = terrm * 10.d0 

      nobs  = ii
      nstat = istat
      stalam  = stalam / dble(nstat)
      stalom1 = stalom1 / dble(nstat)
      stalom2 = stalom2 / dble(nstat)
      stalom  = rad2deg*datan2(stalom2,stalom1)

      if(nsdp.gt.0) then
         sdmeans = sdpmean
         sdpmean = 2.d0*dsqrt(sdpmean / dble(nsdp))
         if(sdpmean.lt.1.d0) sdpmean = 1.d0
      else
         sdpmean = 2.d0
      endif

      if(nsds.gt.0) then
         sdmeans = sdmeans + sdsmean
         sdsmean = 2.d0*dsqrt(sdsmean / dble(nsds))
         if(sdsmean.lt.2.0d0) sdsmean = 2.0d0
      else
         sdsmean = 4.0d0
      endif

      if(nsdp.gt.0 .or. nsds.gt.0 ) then
         sdmeans = dsqrt( sdmeans / dble(nsdp+nsds) )
      else
         sdmeans = 5.d0
      endif

      if(nobs.eq.1) then

         rzo1  = 0.
         nphas = 0

         rdel1 = 3.
         call tauget_mod(rzo1,rdel1,nphas,phcd,ttc,dtdd,
     +                   dtdh,dddp,modnam)
         pmoh  = dble(dtdd(1))

         single = .true.

         if(typctl.gt.4) print *, 'Case: single array observation!'

      endif

      if (check_mult) then
         nobs0 = nobs
         call mult_ons(nobs,nobs0,iev,phase,tt,tts,tt2,azi,azis,p,ps,
     +                 touse,amplit,period,dtphmax,mread,arid)
         if(nobs.ne.nobs0) nobs = nobs0
      endif

c
c     Fix the time system at the earliest onset time for this event.
c

      larid = .false.
      do 15 i = 1,nobs
      tt(i) = tt(i)-timemin
      ttu(i) = tts(i)
c     if(typctl.gt.8) then
c       print*,i,sta(iev(i)),phase(i),tt(i),tts(i),azi(i),azis(i),
c    +       p(i),ps(i),touse(i)
c     endif
      if(arid(i).ne.'        ') larid = .true.
15    continue

      if(epistart) then
      
         elatm  = convlat(epilat0,1)
         elonm  = epilon0
         go to 65

      endif

c
c     At first, let us try to calculate an epicenter from all 
c     available azimuth observations.
c

      sela  = 0.d0
      svla  = 0.d0
      selo1 = 0.d0
      selo2 = 0.d0
      svlo  = 0.d0
      azims = 0.d0
      azimc = 0.d0
      azimr = 0.d0
      iazim = 0
      rpar  = 0.d0
      istater = 0

      jj = 0

      if (nobs.eq.1) then
         azims = dsin(deg2rad*azi(1))
         azimc = dcos(deg2rad*azi(1))
         if(p(1).gt.0.d0) rpar = p(1)
         iazim = 1
         istataz = iev(1)
         go to 51
      endif

      if (azionlyf) typctl = 6

      do 50 i=1,nobs-1

      if(touse(i)(2:2).ne.'A') go to 50
      azi1 = azi(i)
      if(index(phase(i),'pre').gt.0) go to 50
      if(index(phase(i),'KK').gt.0 .and. phase(j)(1:1).ne.'S' .and. 
     +         phase(i)(1:2).ne.'sS')   azi1 = alpha2(azi1-180.d0)
      if(index(phase(i),'P2K').gt.0) azi1 = alpha2(azi1-180.d0)
      if(index(phase(i),'P3K').gt.0) azi1 = alpha2(azi1-180.d0)
      if(index(phase(i),"P'P'").gt.0) azi1 = alpha2(azi1-180.d0)
      if(index(phase(i),"S'S'").gt.0) azi1 = alpha2(azi1-180.d0)

      azims = azims + dsin(deg2rad*azi1)
      azimc = azimc + dcos(deg2rad*azi1)
      iazim = iazim + 1
      istataz = iev(i)

      slat1  = stala(iev(i))
      slat1e = stalae(iev(i))
      slon1  = stalo(iev(i))

      do 20 j=i+1,nobs

      if(touse(j)(2:2).ne.'A') go to 20
      azi2 = azi(j)
      if(index(phase(j),'pre').gt.0) go to 20
      if(index(phase(j),'KK').gt.0 .and. phase(j)(1:1).ne.'S' .and. 
     +         phase(j)(1:2).ne.'sS')   azi1 = alpha2(azi1-180.d0)
      if(index(phase(j),'P2K').gt.0) azi2 = alpha2(azi2-180.d0)
      if(index(phase(j),'P3K').gt.0) azi2 = alpha2(azi2-180.d0)
      if(index(phase(j),"P'P'").gt.0) azi2 = alpha2(azi2-180.d0)
      if(index(phase(j),"S'S'").gt.0) azi2 = alpha2(azi2-180.d0)

      if(i.eq.nobs-1 .and. j.eq.nobs) then
         azims = azims + dsin(deg2rad*azi2)
         azimc = azimc + dcos(deg2rad*azi2)
         iazim = iazim + 1
         istataz = iev(i)
      endif
      if(iev(i).eq.iev(j)) go to 20

      slat2  = stala(iev(j))
      slon2  = stalo(iev(j))

      jj = jj + 1

      if(jj.gt.mloc) then
         print *,'Something wrong with number of locations!'
         go to 9999
      endif

      if(typctl.gt.8 .or. azionlyf) then

         print *,' '
         print *,'    Station BAZ                       dBAZ' 
         print *,'(1) ',sta(iev(i)),azi1,azis(i)
         print *,'(2) ',sta(iev(j)),azi2,azis(j)

      endif

c
c     Calculate distance and angles between the 2 stations
c

      call depi (slat1,slon1,slat2,slon2,del3,dk,ep2,ep1,d2km)

c     if(typctl.gt.7) then
c        print *,'station 1 (lat,lon,distance [km,deg],azimuth: ',
c    +           slat1,slon1,dk,del3,ep1
c        print *,'station 2 (lat,lon,azimuth; ',slat2,slon2,ep2
c     endif

c
c     Now an initial epicenter will be calculated
c

      ierr = 0
      call hyposat_cross(slat1e,slon1,azi1,azis(i),
     +               azi2,azis(j),del3,ep1,ep2,ela(jj),elas(jj),
     +               elo(jj),elos(jj),dismin,typctl,ierr)

      if(ierr.gt.0) then
         jj = jj - 1
         istater = istater + 1
         ierr = ierr + 1
         go to 20
      endif

      sela  = sela + ela(jj)/elas(jj)
      svla  = svla + 1.d0 / elas(jj)

      p1 = deg2rad*elo(jj)
      p2 = 1.d0/(deg2rad*elos(jj))

      selo1  = selo1 + dcos(p1)*p2
      selo2  = selo2 + dsin(p1)*p2
      svlo   = svlo  + p2

20    continue
50    continue

51    nloc = jj

      if(nloc.eq.0) then

        if(rpar .gt. 0.d0) then

           phase_t = ' '
           phidr0 = phase(1)
           if(phidr0.eq.'P1') phidr0 = 'P'
           if(phidr0.eq.'S1') phidr0 = 'S'

           call tauget_ray(phidr0,phase_t,rpar,modnam,zo,ddel,ttray,
     +                     rayokf)

           if(.not.rayokf) then

             ddel =  (14.d0 - rpar)*9.8d0 + 2.0d0
             phidr0 = 'P'

             if(rpar.lt.pdif) then
                ddel = 150.d0
                phidr0 = 'PKPdf'
             endif

             if(rpar.ge.10.0d0) then
                ddel = 23.d0
             endif

             if(rpar.ge.13.0d0) then
                ddel = 10.d0
                phidr0 = 'Pn'
             endif

             if(rpar.ge.14.8d0) then
                ddel = 2.d0
                phidr0 = 'Pb'
             endif

             if(rpar.ge.17.0d0) then
                ddel = 1.d0
                phidr0 = 'Pg'
             endif
           
           endif
         
           if(typctl.gt.0 .and. .not.rayokf) then
                print *,'No distance found. Missing slowness values?'
                print *,'Distance set to ',ddel
           else if(typctl.gt.4) then
              print *,'Initial distance from Station(net): ',ddel
           endif

        endif
  
        if(iazim.le.0 .or. istater.gt.0) then

c
c       Choose a point in the vicinity (1 deg) of the closest
c       station as initial solution:
c

           istatd = istatmin
           azim = 315.d0

           if(typctl.gt.0 .and. .not.rayokf) then
                print *,'No epicenter found. Missing backazimuth ',
     +                'values?'
                print *,'Backazimuth set to ',azim,' degrees'
           endif

        else

           istatd = istataz
           azimr  = datan2(azims,azimc)
           azim   = alpha2(rad2deg*azimr)

        endif

        if (ddel.le.0.d0) ddel = 1.d0

        inddel = 1
        call delazd(stala(istatd),stalo(istatd),azim,ddel,
     +               inddel,elatmg,elonm)
        elatm = convlat(elatmg,1)

        go to 65

      endif

c
c     Now mean source coordinates have to be calculated.
c

      elatm  = sela / svla

      elonm  = rad2deg*datan2(selo2/svlo,selo1/svlo)

      if(nloc.eq.1) then

        sdlat  = elas(1)
        sdlon  = elos(1)

      else if(nloc.gt.1) then

        dla = 0.d0
        dlo = 0.d0
        do 60 i =1,nloc
        dla = dla + q2(ela(i)-elatm) / elas(i)
        p1 = alpha1(elo(i)-elonm)
        dlo = dlo + (p1*p1) / elos(i)
60      continue
        sdlat  = dsqrt(dla / svla)
        sdlon  = dsqrt(dlo / svlo)

      endif

c     The next step is to calculate a first source time. If already 
c     given with input parameter file, this given value will be used.
c
c     We are using the method of Wadati. If we have only one travel
c     time difference S-P we assume Vp/Vs = sqrt(3.). Otherwise
c     we calculate Vp/Vs as a constants for each specific phase type
c     (Pg,Pb,Pn,P,P1)-(Sg,Sb,Sn,S,S1).
c
c     Which travel-time differences between S and P onsets do we have?
c     These travel differences are also used to calculate a source 
c     distance from a station or an array.
c

65    dtkm  = 0.d0
      idtkm = 0

      do 66 i = 1,idtmax
      idt(i)=0
66    continue

      istatd = istatmin
      if(iazim.gt.0 .and. istater.le.0) istatd = istataz

      if (nobs.eq.1) go to 71

      do 702 i = 1,nobs-1

      if(touse(i)(1:1).ne.'T') go to 702

      do 701 j = i+1,nobs

      if(touse(j)(1:1).ne.'T') go to 701

      kmdel = 0.d0

      if(iev(i).eq.iev(j)) then

         phai = phase(i)
         phaj = phase(j)

         dtwad = dabs(tt(j)-tt(i))

        wadati = .true.
        if(dtwad.gt.wadmax) wadati = .false.
        if(dtwad.lt.wadmin) wadati = .false.

        if((phai.eq.'P       ' .and. phaj.eq.'S       ') .or.
     +     (phai.eq.'S       ' .and. phaj.eq.'P       ')) then

          if (wadati) then
             idt(1)        = idt(1) + 1
             dt(1,idt(1))  = dtwad

             if(phai(1:1).eq.'P') then
                idtp(1,idt(1))= i
                idts(1,idt(1))= j
             else
                idtp(1,idt(1))= j
                idts(1,idt(1))= i
             endif
          endif

          if(iev(i).eq.istatd) then
            kmdel = (dtwad/60.d0-1.d0)*1000.d0
            if(kmdel.le.2000.d0) kmdel = dtwad*10.2d0
          endif

        else if((phai.eq.'P1      ' .and. phaj.eq.'S1      ') .or.
     +          (phai.eq.'S1      ' .and. phaj.eq.'P1      ')) then
c
c         For P1 and S1 it is not automatically known which phase it will
c         become. We use AK135 to choose the most
c         presumable phase type for the Wadati Approach.
c
          if (wadati) then
             idtc = 1
             if(dtwad.lt.200.0d0) idtc = 3
             if(dtwad.lt.17.6d0)  idtc = 2

             idt(idtc)          = idt(idtc) + 1
             dt(idtc,idt(idtc)) = dtwad

             if(phai(1:1).eq.'P') then
                idtp(idtc,idt(idtc))= i
                idts(idtc,idt(idtc))= j
             else
                idtp(idtc,idt(idtc))= j
                idts(idtc,idt(idtc))= i
             endif
          endif

          if(iev(i).eq.istatd) then
             if(idtc.eq.1) then
               kmdel = (dtwad/60.d0-1.d0)*1000.d0
               if(kmdel.le.2000.d0) kmdel = dtwad*10.2d0
             else if(idtc.eq.2) then
               kmdel = dtwad*8.58d0
             else if(idtc.eq.3) then
               kmdel = dtwad*10.2d0
             endif
          endif

        else if((phai.eq.'Pg      '                      .and.
     +          (phaj.eq.'Sg      '.or.phaj.eq.'Lg      ')) .or.
     +          (phaj.eq.'Pg      '                      .and.
     +          (phai.eq.'Sg      '.or.phai.eq.'Lg      ')))then

          if (wadati) then
             idt(2)        = idt(2) + 1
             dt(2,idt(2))  = dtwad

             if(phai(1:1).eq.'P') then
                idtp(2,idt(2))= i
                idts(2,idt(2))= j
             else
                idtp(2,idt(2))= j
                idts(2,idt(2))= i
             endif
          endif

          if(iev(i).eq.istatd) kmdel = dtwad*8.58d0

        else if((phai.eq.'Pn      '.and.phaj.eq.'Sn      ') .or.
     +          (phaj.eq.'Pn      '.and.phai.eq.'Sn      '))then

          if (wadati) then
             idt(3)        = idt(3) + 1
             dt(3,idt(3))  = dtwad

             if(phai(1:1).eq.'P') then
                idtp(3,idt(3))= i
                idts(3,idt(3))= j
             else
                idtp(3,idt(3))= j
                idts(3,idt(3))= i
             endif
          endif

          if(iev(i).eq.istatd) kmdel = dtwad*10.2d0

        else if((phai.eq.'Pb      '.and.phaj.eq.'Sb      ') .or.
     +          (phaj.eq.'Pb      '.and.phai.eq.'Sb      '))then

          if (wadati) then
             idt(4)        = idt(4) + 1
             dt(4,idt(4))  = dtwad

             if(phai(1:1).eq.'P') then
                idtp(4,idt(4))= i
                idts(4,idt(4))= j
             else
                idtp(4,idt(4))= j
                idts(4,idt(4))= i
             endif
          endif

          if(iev(i).eq.istatd) kmdel = dtwad*9.47d0

        else if(phai.eq.'Pn      '                          .and.
     +      (phaj.eq.'Sg      '.or. phaj.eq.'Lg      ')) then

          if(iev(i).eq.istatd) kmdel = dtwad*6.02d0

        else if(phaj.eq.'Pn      '                          .and.
     +      (phai.eq.'Sg      '.or. phai.eq.'Lg      ')) then

          if(iev(i).eq.istatd) kmdel = dtwad*6.02d0

        endif

        if(kmdel.gt.0.d0) then
           idtkm = idtkm + 1
           dtkm  = dtkm + kmdel
        endif

      endif

701   continue
702   continue

71    inet = 0

      if(idtkm.ne.0 .and. nloc.eq.0 .and. .not.epistart) then

         dtkm = dtkm / dble(idtkm)

         if(iazim.gt.0 .and. istater.le.0) then

            inddel = 2
            call delazd(stala(istataz),stalo(istataz),azim,dtkm,
     +               inddel,elatmg,elonm)
            elatm = convlat(elatmg,1)
 
            if(typctl.gt.0) then
               if (dtkm.le.0d0) then
                 print *,'Epicenter set to station ', sta(istataz)
               else
                 print *,'Epicenter set from station ',
     +               sta(istataz),': backazimuth',azim,' deg, delta',
     +               dtkm,' km' 
               endif
            endif

            sdlatg = dtkm*grad1
            sdlon  = dtkm*grad1

         else 

            if(plflag) then

               call plane(stala,stalo,iev,tt,nobs,azim,dazi,
     +                    ray,dray,phipl,touse,phase,jref,typctl)

               if(jref.gt.0 .and. dazi.lt.90.d0 .and. dray.lt.4.d0) then

                   inddel = 2
                   call delazd(stala(iev(jref)),stalo(iev(jref)),azim,
     +                         dtkm,inddel,elatmg,elonm)
                   elatm = convlat(elatmg,1)

                   if(typctl.gt.0) then
                      if (dtkm.le.0d0) then
                        print *,'Epicenter set to station ',sta(jref)
                      else
                        print *,'Epicenter set from station ',
     +                            sta(jref),'after plane wave fit: ',
     +                      'backazimuth',azim,' deg, delta',dtkm,' km' 
                      endif
                   endif

                   sdlatg = dazi*dtkm*grad1
                   sdlon  = dazi*dtkm*grad1

                   go to 72

                endif

            endif

              if(dtkm.gt.120.d0 .and. nstat.gt.1) then
            
               elatm = stalam
               elonm = stalom

               inet = 1

               print *,'Epicenter set in center of station net '

            endif

         endif
            
      else if (nloc.eq.0 .and. nstat.gt.1 .and. .not.epistart) then
            
         elatm = stalam
         elonm = stalom
         inet = 1
         print *,'Epicenter set in center of station net '

      endif

72    elatmr  = deg2rad*elatm
      elatmg  = convlat(elatm,2)
      elatmgr = deg2rad*elatmg
      coelatm = 90.d0 - elatm
      coelatmr= deg2rad*coelatm
 
      if(sdlat.ne.0.d0) then
         sdlatg= sdlat / (eps*q2(dcos(elatmr))+q2(dsin(elatmr)))
      else
         sdlat = sdlatg*eps/(q2(dcos(elatmgr))+eps*q2(dsin(elatmgr)))
      endif

      if(dismin.lt.pi) then
         dismin = dismin * rad2deg * 5.d-2
         if(sdlatg.lt.dismin) then
            sdlatg = dismin
            sdlat = sdlatg*eps/(q2(dcos(elatmgr))+
     +                     eps*q2(dsin(elatmgr)))
         endif 
         if(sdlon.lt.dismin) sdlon = dismin
      endif

      if(typctl.gt.0) then
        if(nloc.gt.0) then
           print *,' '
           print*,'Mean epicenter calculated from ',nloc,
     +                  ' observation(s)'
        else
           print*,'Set epicenter '
        endif
        print*,'(Mean) epicenter lat: ',elatmg,' +/- ',sdlatg
        print*,'(Mean) epicenter lon: ',elonm,' +/- ',sdlon
      endif

      if(output) then

        write(11,'(/''Parameters of initial solution ('',
     +              ''+/- 1 standard deviation):''/)') 

        if(nloc.gt.0) then

           write(11,'(''Mean epicenter calculated from'',i5,
     +              '' backazimuth observation pairs'')') nloc
           write(11,'(''Mean epicenter lat:'',f9.3,'' +/- '',f9.3,
     +              '' [deg]'')')        elatmg,sdlatg
           write(11,'(''Mean epicenter lon:'',f9.3,'' +/- '',f9.3,
     +              '' [deg]''/)')        elonm,sdlon

        else if(.not.epistart) then

             write(11,'(''No location from multiple azimuth '',
     +            ''observations found. Missing backazimuth '',
     +            ''values?'')')

           if(inet.ne.0) then
              write(11,'(''Epicenter set in the center of station'',
     +                   '' net '')')
           else if (iazim.gt.0 .and. istater.le.0 .and. dtkm.eq.0.d0)
     +             then
               if (ddel.le.0.d0) then
                  write(11,'(''Epicenter set to station '',a)') 
     +               sta(istatd)
               else
                  write(11,'(''Epicenter set from station '',
     +               a,'' with backazimuth'',f6.1,'' [deg], delta'',
     +               f7.2,'' [deg]'')') sta(istatd),azim,ddel
               endif
           else 
               if (dtkm.le.0.d0) then
                  write(11,'(''Epicenter set to station '',a)') 
     +               sta(istatd)
               else
                  write(11,'(''Epicenter set from station '',
     +               a,'' with backazimuth'',f6.1,'' [deg], delta'',
     +               f7.2,'' [km]'')') sta(istatd),azim,dtkm
               endif
           endif

           write(11,'(''Epicenter lat:'',f9.3,'' [deg]'')')  elatmg
           write(11,'(''Epicenter lon:'',f9.3,'' [deg]''/)') elonm

        else if(epistart) then

           write(11,'(''Initial Epicenter set by input file'')')
           write(11,'(''Epicenter lat:'',f9.3,'' [deg]'')')  elatmg
           write(11,'(''Epicenter lon:'',f9.3,'' [deg]''/)') elonm
        
        endif
      endif

      if(azionlyf) go to 9999
c
      if(typctl.gt.8) then

        do 75 i=1,nstat

        call depi(stala(i),stalo(i),elatmg,elonm,del(i),dk,azie(i),
     +                  baz(i),d2km)

        if(azi(i).ge.0.d0) then
          print *,sta(i),del(i),azie(i),baz(i),azi(i),azi(i)-baz(i)
        else
          print *,sta(i),del(i),azie(i),baz(i)
        endif

75      continue

      endif

      do 80 i=1,idtmax

      vpvs(i) = 0.d0

      if(idt(i).eq.0) go to 80

      if(idt(i).eq.1) then

        vpvs(i) = dsqrt(3.d0)
        f1      = 1.d0/(vpvs(i)-1.d0)

        to(i)   = tt(idtp(i,1)) - dt(i,1)*f1
        tos(i)  = dpythag((1.d0+f1)*tts(idtp(i,1)),f1*tts(idts(i,1)))
        vpvss(i)= 0.5d0

      else if(idt(i).eq.2) then

        f1 = dt(i,2)-dt(i,1)

        if ( dabs(tt(idtp(i,2))-tt(idtp(i,1))).lt.0.01d0) go to 80
        f2 = 1.d0 / (tt(idtp(i,2))-tt(idtp(i,1)))

        am      = f1*f2
        am1     = 1.d0 / am
        vpvs(i) = am + 1.d0
        to(i)   = tt(idtp(i,1)) - dt(i,1)*am1

        f3 = am*f2
        f4 = am1*am1
        vpvss(i)= dsqrt ( q2(( f2+f3)*tts(idtp(i,1))) +
     +                    q2((-f2-f3)*tts(idtp(i,2))) +
     +                    q2(  f2    *tts(idts(i,1))) +
     +                    q2( -f2    *tts(idts(i,2))) )
        
        tos(i) = dsqrt( q2( (1.d0+am1+dt(i,1)*(f2+f3)*f4 )
     +                                         *tts(idtp(i,1)) )   +
     +                  q2( (         dt(i,1)*(-f2-f3)*f4 )
     +                                   *tts(idtp(i,2)) )   +
     +                  q2( (    -am1+dt(i,1)*f2     *f4 )
     +                                   *tts(idtp(i,2)) )   +
     +                  q2( (        -dt(i,1)*f2     *f4 )
     +                                 *tts(idtp(i,2)) )   )

      else 

        do 77 j=1,idt(i)

        f1 = 1.d0
        if(tts(idtp(i,j)).ne.0.d0 .or. tts(idts(i,j)).ne.0.d0)
     +     f1 = 1.d0/dpythag(tts(idtp(i,j)),tts(idts(i,j)))

         f2 = dt(i,j) * f1
         f3 = tt(idtp(i,j)) * f1

        ddl(j) = f2
        ggl(j,1) = f3
        ggl(j,2) = f1

c       if(typctl.gt.6) then
c         print *,j,f1,f3,f2
c       endif

77      continue
        
        im = 2
        in = idt(i)
        call dlsq(in,im)

        if (lsqerr.gt.0) then
          if(typctl.gt.4) print*,'Wadati TYPE ',i,' LSQ-Error: ',lsqerr
          go to 80
        endif

c        if(typctl.gt.6) then
c          print *,aal(1),vvl(1)
c          print *,aal(2),vvl(2)
c          do 78 j=1,idt(i)
c78        print *,ddl(j),rrl(j)
c        endif


        vpvs(i) = aal(1) + 1.d0
        vpvss(i)= vvl(1)
        to(i)   = -aal(2) / aal(1)
        tos(i)  = dpythag(to(i)*vvl(2)/aal(2),to(i)*vvl(1)/aal(1))

      endif

      if(i.eq.1 .and. dabs(to(i)).gt.1300.d0) then
         to(i)=dsign(1300.d0,to(i))
         tos(i)=dabs(to(i))
      endif
      if(i.eq.2 .and. dabs(to(i)).gt.150.d0) then
         to(i)=dsign(150.d0,to(i))
         tos(i)=dabs(to(i))
      endif
      if(i.eq.3 .and. dabs(to(i)).gt.400.d0) then
         to(i)=dsign(400.d0,to(i))
         tos(i)=dabs(to(i))
      endif
      if(i.eq.4 .and. dabs(to(i)).gt.250.d0) then
         to(i)=dsign(250.d0,to(i))
         tos(i)=dabs(to(i))
      endif

      if(typctl.gt.0) then
         print *,'S-P Travel-time difference type ',i
         print *,'Source time from ',idt(i),' observation(s)'
         print *,'   to= ',to(i),  ' +/- ',tos(i)
         print *,'Vp/Vs= ',vpvs(i),' +/- ',vpvss(i)
         if(vpvs(i).lt.1.2d0.or.vpvs(i).gt.2.5d0) print *,'Not used!'
      endif

      if (vpvss(i).gt.9.99d0) vpvss(i)=9.99d0
      if (tos(i).gt.9999.9d0) tos(i)=9999.9d0
      if(output) then
         if(vpvs(i).lt.1.2d0.or.vpvs(i).gt.2.5d0) then
           write(11,'(''S-P Travel-time difference type'',i2,
     +      '' with'',i4,'' observation(s)''/''   to='',f14.1,'' +/-'',
     +      f6.1,'' [s] Vp/Vs= '',f4.2,'' +/- '',f4.2,'' not used!'')')
     +      i,idt(i),to(i)+timemin,tos(i),vpvs(i),vpvss(i)
         else
           write(11,'(''S-P Travel-time difference type'',i2,
     +      '' with'',i4,'' observation(s)''/''   to='',f14.1,'' +/-'',
     +      f6.1,'' [s] Vp/Vs= '',f4.2,'' +/- '',f4.2)') i,idt(i),
     +      to(i)+timemin,tos(i),vpvs(i),vpvss(i)
         endif

      endif

80    continue

c
c     now follows the statistics over all estimated to-values
c
      sto     = 0.d0
      stos    = 0.d0
      svpvs   = 0.d0
      svpvss  = 0.d0
      ito     = 0
 
      do 82 i = 1,idtmax
 
      if(idt(i).eq.0) go to 82
      if(vpvs(i).lt.1.2d0.or.vpvs(i).gt.2.5d0) go to 82

      ito = ito + 1

      sto   = sto  + to(i)/tos(i)
      stos  = stos + 1.d0 /tos(i)
 
      svpvs = svpvs  + vpvs(i)/vpvss(i)
      svpvss= svpvss + 1.d0   /vpvss(i)
 
82    continue
 
      dto   = 0.d0
      dvpvs = 0.d0
 
      if(ito.eq.0) then
         if(nloc.ge.1) then

           call depi(stala(istatmin),stalo(istatmin),elatmg,elonm,
     +               del(istatmin),dk,azie(istatmin),baz(istatmin),d2km)
           rzo1 = 0.
           rdel = real(del(istatmin))
           nphas = 0

           call tauget_mod(rzo1,rdel,nphas,phcd1,ttc1,dtdd1,
     +                         dtdh1,dddp1,modnam)

           tom  = - dble(ttc1(1))

         else

           if(ttray.gt.0d0) then
              tom  = -ttray
           else
              tom  = -dtp0/2.d0
           endif

         endif
         vpvsm = 0.d0
         sdto = dtp0
         go to 85
      else
         tom   = sto / stos
         vpvsm = svpvs / svpvss
      endif
 
      if(ito.eq.1) then
         sdto   = 1.d0 / stos
         sdvpvs = 1.d0 / svpvss
         go to 85
      endif

      do 83 i =1,idtmax

      if(idt(i).eq.0) go to 83
 
      dto   = dto   + q2(to(i)-tom) / tos(i)
      dvpvs = dvpvs + q2(vpvs(i)-vpvsm) / vpvss(i)
 
83    continue
 
      sdto    = dsqrt(dto / svpvs)
      sdvpvs  = dsqrt(dvpvs / svpvss)
 
85    continue

      tome = tom + timemin

      if(ito.gt.0) then
         if(output) then
            if(typctl.gt.0) then
               print*,'Mean source time: ',tome,' +/- ',sdto
               print*,'Mean       vp/vs: ',vpvsm,' +/- ',sdvpvs
            endif
            write(11,'(''Mean source time:'',f15.3,'' +/- '',f7.3,
     +              '' [s]'')')  tome,sdto
            write(11,'(''Mean       vp/vs:'',f15.3,'' +/- '',f7.3)') 
     +               vpvsm,sdvpvs
         endif
      else
         if(output) then
            if(typctl.gt.0) then
              print *,'Source time (set): ',tome,' [s]'
            endif
            write(11,'(''Source time (set):'',f15.3,'' [s]'')') tome
         endif
      endif

c
c     In any case, we use (if set) the initial source time and its standard
c     deviation from hyposat-parameter file or ISF-input file.
c     
      if (tome0 .gt. -2840140800.0d0) then

         tome = tome0

         tom = tome - timemin

         if (stome0 .gt. 0.d0) then
            sdto = stome0
         endif

         if(output) then
            if(typctl.gt.0) then
               print*,'Source time (from input): ',tome,' +/- ',sdto
            endif
            write(11,'(/''Source time (from input):'',f15.3,'' +/- '',
     +                f7.3,'' [s]''/)') tome,sdto
         endif

      endif

c
c     Now new source parameters can be calculated by
c     several iterations using the GMI algorithm. 
c
c     For the first iteration we use as initial solution the read in 
c     source depth, the source time to, and the epicenter coordinates 
c     elatmg and elonm.
c
      iter = 0

      if(czo.eq.'D') sdzo = sdzo1
      if(czo.eq.'F' .or. czo.eq.'B') sdzo = 1.d0

      rs(1) = sdto
      rs(2) = sdlat
      rs(3) = sdlon
      rs(4) = sdzo

      nextiter = 0
      in = 0
      dtp = dtp0
      dts = dts0

      dchang = dchang0

c
c     At first, we build the Jacobi-matrix
c     (loop 300 and 301)
c

100   continue

      last = .false.

      if (ibad0.gt.3) then
         print *,'Could not find a stable solution for these data'
         if(output) then
            write(11,'(/''Could not invert these data!'')')
         endif
         go to 9999
      endif

      iter = iter + 1

      if((check.le.setcheck2 .or. nextiter1.eq.1) .and. 
     +    iteraz.eq.0 ) dtmflag = .true.

      if ((lastfixt .and. check.le.setcheck2) .or. dtmflag) then

         if(zo.le.0.1d0 .and. .not.locgeo) zo=0.d0
         if(zo.lt.0.0d0 .and. locgeo) zo=0.d0

         if(lastfixm) then

            dazim0 = dazim1
            dpam0  = dpam1

         else
         
            f1 = dmax1(0.75d0,rmso)

            dtm    = f1 * 3.d0
            dazim0 = dmin1(15.d0,dazim1)
            dpam0  = dmin1(2.d0,dpam1)

            if (nrms1.gt.5 .and. nrms1.lt.10) then
               dtm  = f1 * 2.d0
               dazim0 = dmin1(10.d0,dazim1)
               dpam0  = dmin1(1.5d0,dpam1)
            else if (nrms1.gt.10) then
               dtm    = f1 * 1.2d0
               dazim0 = dmin1(5.d0,dazim1)
               dpam0  = dmin1(1.d0,dpam1)
            endif
            dtm0 = dtm * 2.d0

         endif

      else
         dtm  = dtm2
         dtm0 = dtm2
         dazim0 = dazim1
         dpam0  = dpam1
      endif

      if(lastfixm .and. dtmflag) then
        dtmp = resmaxp
        dtms = resmaxs
        dtm0 = dtm2
      else
        dtmp = dmin1(dtm,dtp)
        if(dtmp.lt.sdpmean) dtmp = sdpmean
        dtms = dmin1(2.d0*dtm,dts)
        if(dtms.lt.sdsmean) dtms = sdsmean
        dtm0 = dmin1(dtm0,dtm2)
      endif

      ifixaz = 0
      ifixto = 0

      rmsold = rmso
      iremo  = 0

101   continue

      if(dtmp.gt.1000.d0) dtmp = 1000.d0
      if(dtms.gt.1400.d0) dtms = 1400.d0
      if(dtm0.gt.1400.d0) dtm0 = 1400.d0
      if(dpam0.gt.15.d0)  dpam0 = 15.d0
      if(dazim0.gt.90.d0) dazim0 = 90.d0

      stato = ' '
      jj    = 0
      jdt   = 0
      jazi  = 0
      jpa   = 0
      rms1  = 0.d0
      nrms1 = 0
      datmax = 0.d0
      nzo    = 0
      fmud  = 0.0d0
      fmus  = 0.0d0
      fmuh  = 0.0d0
      ismu  = 0

      rzv(1)    = -999.d0
      rzv(2)    = -999.d0

      do 300 i = 1,nobs

      used(i)= ' '

      phaseu(i) = ' '

      dinv(i,1) = 0.d0
      dinv(i,2) = 0.d0
      dinv(i,3) = 0.d0
      dinv(i,4) = 0.d0

      ttt(i) = 0.d0

      if (sta(iev(i)).ne.stato) then

         call depi(stala(iev(i)),stalo(iev(i)),elatmg,elonm,
     +              del(iev(i)),dk,azie(iev(i)),baz(iev(i)),d2km)

         if (kmout) then
            if (dk .gt. dismaxst ) go to 300
            if (dk .lt. disminst ) go to 300
         else
            if (del(iev(i)) .gt. dismaxst ) go to 300
            if (del(iev(i)) .lt. disminst ) go to 300
         endif
         if (dk .lt. 1.d-4 ) del(iev(i))=1.d-4
         if (del(iev(i)) .lt. 1.d-6 ) del(iev(i))=1.d-6

         stato = sta(iev(i))

         if(typctl.gt.8) then
           print *,'STATION EPI: ',i,stato,del(iev(i)),dk
         endif

         imod2 = 0
         if (touse(i)(7:7).eq.'2' .and. mod2flag) then
           imodn(2) = 2
           modn = modnam2
         else if (touse(i)(7:7).eq.'3' .and. mod3flag) then
           imodn(3) = 2
           modn = modnam3
         else if (touse(i)(7:7).eq.'4' .and. mod4flag) then
           imodn(4) = 2
           modn = modnam4
         else
           if (iloc) imod2 = 1
           imodn(1) = 2
           modn = modnam
         endif

         rzo = sngl(zo)
         rdel = sngl(del(iev(i)))

         costalat  = 90.d0-stalae(iev(i))
         costalatr = deg2rad*costalat

         fla1 = sngl(coelatmr)
         razi = sngl(azie(iev(i)))
         
         loctt = 0

         if(imod2.eq.0 .or. del(iev(i)).gt.rmax .or. zo.gt.zmax) then

            nphas0 = 0
            call tauget_mod(rzo,rdel,nphas0,phcd,ttc,dtdd,
     +                      dtdh,dddp,modn)

            if(czo.eq.'D') then

               rzo1 = rzo-1.
c              no event above the surface!
               if(rzo1.lt.0.0) rzo1 = 0.0
c
               rzo2 = rzo+1.
c              no event deeper than 799. km !
               if(rzo2.ge.800.0) rzo2=799.0

               nphas1 = 0
               nphas2 = 0

               if(rzo-rzo1.gt.0.0) then
                  call tauget_mod(rzo1,rdel,nphas1,phcd1,ttc1,dtdd1,
     +                            dtdh1,dddp1,modn)
               endif

               if(rzo2-rzo.gt.0.0) then
                  nphas2 = 0
                  call tauget_mod(rzo2,rdel,nphas2,phcd2,ttc1,dtdd2,
     +                            dtdh1,dddp1,modn)
               endif

               do 280 k=1,nphas0

               dpdh(k) = 0.d0

               if(abs(rzo2-rzo1).le.0.001) go to 280

               if(nphas1.gt.0 .and. nphas2.gt.0) then
                  do 250 ki=1,nphas1
                  do 250 kj=1,nphas2
                  if(phcd1(ki).eq.phcd(k)  .and. 
     +               phcd(k).eq.phcd2(kj)) then
                     dpdh(k) = (dtdd2(kj)-dtdd1(ki))/dble(rzo2-rzo1)
                     go to 280
                  endif
250               continue
               endif

               if(rzo.gt.rzo1 .and. nphas1.gt.0) then
                  do 260 ki=1,nphas1
                  if(phcd(k).eq.phcd1(ki)) then
                     dpdh(k) = (dtdd(k)-dtdd1(ki))/dble(rzo-rzo1)
                     go to 280
                  endif
260               continue
               endif

               if(rzo2.gt.rzo .and. nphas2.gt.0) then
                  do 270 kj=1,nphas2
                  if(phcd(k).eq.phcd2(kj)) then
                     dpdh(k) = (dtdd2(kj)-dtdd(k))/dble(rzo2-rzo)
                     go to 280
                  endif
270               continue
               endif

280            continue

            endif

            nphas = nphas0

         else 

            ierr = 0
            indph = istaph(iev(i))*10000 + indph0
            elatc = elatmg
            elonc = elonm

            elat2 = stala(iev(i))
            elon2 = stalo(iev(i))
            sdep = 0.d0
            if(locsta) sdep  = - stael(iev(i))

            nphas = 0
            call ttloc(rzo,rdel,czo,nphas,ttc,dtdd,dtdh,dpdh,dddp,
     +                 phcd,rmax,sdep,typctl,ierr,indph,emerout)

            if(ierr.ne.0) then
               print *, 'Error in getting travel-time tables for:'
               print *, sta(iev(i)), ' in ',rdel,' deg --> not used!'
               go to 300
            endif
            loctt = 1

         endif

         f1 = dcos(coelatmr)
         f3 = dsin(coelatmr)

         f2 = dcos(costalatr)
         f4 = dsin(costalatr)

         f5 = deg2rad*alpha1(stalo(iev(i))-elonm)
         f5 = dabs(f5)

         f6 = dcos(f5)
         f8 = dsin(f5)

         f7 = dsin(deg2rad*del(iev(i)))

         if(baz(iev(i)).le.180.d0) then
            alpha = deg2rad*baz(iev(i))
         else
            alpha = deg2rad*(360.d0-baz(iev(i)))
         endif

         deldla =  (f4*f1*f6 - f3*f2) / f7
         deldlo =  f4*dsin(alpha)

         f9  = dcos(alpha) * f7*f7
         f10 = dcos(deg2rad*del(iev(i)))

         dazidla =  -f8*(f1*f7+f3*f10*deldla)/f9
         dazidlo =   f3*(f6*f7-f8*f10*deldlo)/f9

          if(baz(iev(i)).gt.180.d0) then
            dazidla = -dazidla
            deldlo  = -deldlo
          endif
          
c         if(typctl.gt.8) then
c             print *,'baz ',baz(iev(i)),deldla,deldlo,dazidla,
c    +              dazidlo,(stalo(iev(i))-elonm),del(iev(i))
c
c         endif
      endif
  
      phid = phase(i)

      dpa = 0.d0
      surf = .false.

      if(uppcas(phid(1:1)).ne.'P' .and. uppcas(phid(1:1)).ne.'S' .and.
     +   phid(1:1).ne.'R' .and. phid(1:1).ne.'L' .and. 
     +   phid(1:1).ne.'T' .and. phid(1:2).ne.'IS'  ) go to 298

      if(sglgdis.gt.0.d0) then
         if(phid.eq.'Sg' .and. dk.ge.sglgdis) phid = 'Lg'
         if(phid.eq.'Lg' .and. dk.lt.sglgdis) phid = 'Sg'
      endif

      if(touse(i)(1:1).eq.'m') go to 298

      if(phid.eq.'Rg') then
         if(.not.rgsurf) go to 298
         if(dk.le.400.d0) then
            surf = .true.
            vsurf = vrg
         else
            phid = 'LR'
         endif
      endif
                
      if(phid.eq.'Lg' ) then
         if(.not.lgsurf) go to 298
         if(dk.le.3000.d0) then
            surf = .true.
            vsurf = vlg
         else 
            phid= 'LR'
         endif
      endif
                
c
c     LR is also assumed for far-distant Lg and Rg
c
      if(phid.eq.'LR') then
         if(.not.lpsurf) go to 298
         surf = .true.
         vsurf = vlr
      endif
                
      if(phid.eq.'LQ') then
         if(.not.lpsurf) go to 298
         surf = .true.
         vsurf = vlq
      endif

      if(phid.eq.'T ') then
         if(.not.tsurf) go to 298
         surf = .true.
         vsurf = vt
      endif

      if(phid.eq.'IS') then
         if(.not.isurf) go to 298
         surf = .true.
         vsurf = vi
      endif

      if (surf) then
         nphas       = nphas + 1
         ttc(nphas)  = dk/vsurf
         dtdd(nphas) = d2km/vsurf
         dtdh(nphas) = 0.d0
         phcd(nphas) = phid
         dddp(nphas) = 0.d0
         dpdh(nphas) = 0.d0
      endif
                
      icha = 0
      
      if((phid(2:3).eq.'n ' .or. phid(2:3).eq.'g ' .or. phid(2:3).eq.
     +   'b ') .and. rdel.gt.30. ) phid(2:3)='  '

295   continue

      if(single .and. insar.le.0 .and. iter.eq.1) then
         phid = phidr0
         phase(1) = phid
      endif
      
      nphass = nphas

      if(phase(i).eq.'P1') then
         nphass = 1
         if (icha.ne.0) then
            if (phid(1:1).eq.'P')   nphass = 2
            if (phid(1:2).eq.'PK' .or. phid(1:3).eq.'Pdi') nphass = 10
         endif
      endif

      do 297 j = 1,nphass

      phid1 = phcd(j)

      firston = .false.

      if (firstph) then

         if (j.eq.1) then
            firston = .true.
            firsts  = .false.
         else
            if (phid1(1:1).eq.'S' .and. .not. firsts) then
                firston = .true.
                firsts  = .true.
            endif
         endif

      endif

c
c     Let's take the first phase from IASPEI-1991-type tables,
c     which fits the phase name of the onset
c
c     print *,'---> phid1, phid ',phid1,phid

      if(phid.eq.'P1' .and. phid1(1:1).eq.'P' .and.
     +   phid1(2:2).ne.'P' .and. phid1(2:2).ne.'S' .and.
     +   (phid1(3:3).eq.' ' .or. phid1(2:4).eq.'dif' .or.
     +    phid1(1:5).eq.'PKPdf')) then

        phid = phid1

        if(p(i).lt.0.9d0*pdif .and. phid1(1:3).eq.'Pdi') then
           if(del(iev(i)).gt.110.d0) then
             phid = 'PKPdf'
             go to 297
           endif
        endif

      endif
 
      if(phid.eq.'S1' .and. phid1(1:1).eq.'S' .and. 
     +   phid1(2:2).ne.'S' .and. phid1(2:2).ne.'P' .and.
     +   (phid1(3:3).eq.' ' .or. phid1(1:5).eq.'SKSac')) then

        phid = phid1

      endif

      if(phid.eq.'PKPdf' .and. nobs.gt.2) then
         
         if(phid1(1:5).eq.'PKiKP' .and. 
     +                       del(iev(i)).ge.90.d0) phid='PKiKP'

      endif

      if(icha.eq.0) then
 
        if(phid1(2:4).eq.'dif') then
          if(phid.eq.'P') phid='Pdif'
          if(phid.eq.'S') phid='Sdif'
        endif

        if(phid(1:5).eq.'PKPab' .and. del(iev(i)).lt.150.d0) 
     +     phid(1:5) = 'PKP  '
 
      endif

      if(phid1.eq.phid) then

         dpa  = dtdd(j)
         pa   = sngl(dpa)
         dpaa = dabs(dpa)

         if(single .and. rayokf) then

            if(dabs(p(i)-dpaa) .ge. 0.01d0) go to 297

         endif

c
c     Phase is identified and now several corrections can be applied
c     (if possible and/or requested).
c

c
c     Any ellipticity correction for this phase?
c

         rzoe = rzo

         ierre = 0
         ecor = 0.

         if (.not.surf .and. .not.locgeo) then
c
c          Ellipticity corrections are yet not available for sources 
c          deeper than 700 km. Therefore we accept a small error
c          and set the source depth to 700. km
c
c          No ellipticity corrections for locations in the case that 
c          the 'local geometry' switch is set.
c
           if(rzoe.gt.700.) rzoe=700.

           call ellip(fla1,razi,rdel,rzoe,phid1,pa,ecor,ierre)

         endif

c
c     Any static station correction for this phase?
c
         phase_t = phase_type(phid1)

         statict = 0.d0
         if(firstph) then
            if (firston) then
               if(phase_t.eq.'P') statict = statp(iev(i))
               if(phase_t.eq.'S') statict = stats(iev(i))
            endif
         else
            if(phase_t.eq.'P') statict = statp(iev(i))
            if(phase_t.eq.'S') statict = stats(iev(i))
            if(surf .and. phid1(2:2).eq.'g') statict = statr(iev(i))
         endif

c
c     Any elevation or Crust 1.0 correction for this phase?
c
         th     = 0.d0
         tcrust = 0.d0
         dph = 0.d0

         if(vlflag .and. .not.surf) then

           hsta = stael(iev(i))

           vloc = 99999.d0
           if(phase_t.eq.'P') vloc = stavp(iev(i))
           if(phase_t.eq.'S') vloc = stavs(iev(i))

           if(vloc .lt. 999.d0) then

              if(imo.eq.2 .or. imo.eq.4 .and. loctt.ne.1) then

                 elatc  = stalae(iev(i))
                 elonc  = stalo(iev(i))
                 indr = 1
                 zoc = 999.d0
                 tcrust = crustc(phase_t,dpaa,zoc,indr,typctl)

                 if(dabs(tcrust).gt.0.d0) hsta = hsta - elev

              endif

              if(locsta .and. hsta.lt.0.d0 .and. loctt.eq.1) hsta = 0.d0

              if(dabs(hsta).gt.0.d0) then

                 radkm = deg2rad*radloc(stala(iev(i)),1)
                 phin = vloc*dpaa/radkm

                 if(phin.lt.1.d0) then

                    dl = hsta / dcos(dasin(phin))

                    ddis2 = q2(dl)-q2(hsta)

                    if(ddis2.gt.0.d0) then

                       ddis = dsqrt(ddis2)/radkm

                       if(dl.gt.0.d0) then
                          th = dl/vloc - dpaa*ddis/radkm
                          dph = ddis * dddp(j)
                       else if(dl.lt.0.d0) then
                          th = dl/vloc + dpaa*ddis/radkm
                          dph = ddis * dddp(j)
                       endif

                    endif

                 endif
              endif
           endif
         endif

c        if(typctl.gt.8) then
c           print *,i,tt(i),tom,ttc(j),ecor,vloc,phase_t,hsta,th
c        endif

c
c        We have eventually to correct this phase for the local 
c        structure at the reflection point at the surface (if 
c        CRUST 1.0 or a local/regional model is available).
c
         trefl  = 0.d0
         trefls = 0.d0
         treflp = 0.d0
         if( .not.surf. and. touse(i)(5:5).eq.'R' .and. loctt.eq.0 .and.
     +      imo.ne.3) then

            if(phid(1:1).ne.'p' .and. phid(1:1).ne.'s') then
               if(phid(1:1).eq.phid(2:2) .and. 
     +            phid(1:1).eq.phid(3:3))  then
                  goto  296
               endif
               if(phid(1:2).eq.phid(3:4) .and. 
     +            phid(1:2).eq.phid(5:6))  then
                  goto  296
               endif
            else
               if(phid(2:2).eq.phid(4:4) .and. 
     +            phid(2:2).eq.phid(5:5))  then
                  goto  296
               endif
               if(phid(2:4).eq.phid(5:6) .and. 
     +            phid(2:4).eq.phid(7:8))  then
                  goto  296
               endif
            endif

            chgcas = uppcas(phid)
            phase_t = chgcas(1:1)

            fmult = 1.d0
            del0  = dirdel(dpaa,zo,fmult,phase_t)
            azi0  = azie(iev(i))

            if(phid(1:1).ne.'p' .and. phid(1:1).ne.'s') then
               if(dpa.ge.0.d0) then
                  del0 = (del(iev(i))+del0)/2.d0
               else
                  del0 = (360.d0 - del(iev(i)) - del0)/2.d0
                  azi0 = alpha2(azie(iev(i))-180.d0)
               endif
            endif
                
            inddel = 1
            call delazd(elatmg,elonm,azi0,del0,inddel,elatc,elonc)
    
c
c     correction for depth phases (e.g.: pP, sP...)
c 
            if(phid(1:1).eq.'p' .or. phid(1:1).eq.'s' .and.
     +         phase_t.ne.' ' .and. del0.le.rmax)  then

              if (zo.ge.zmax) then

                  if ((phid(1:1).eq.'p' .and. phid(2:2).eq.'P') .or.
     +                (phid(1:1).eq.'s' .and. phid(2:2).eq.'S')) then
                      indr = 2
                  else 
                      indr = 3
                  endif

                  zoc    = 9999.d0

              else if(zo.lt.zmax) then

                  zoc = zo

                  if ((phid(1:1).eq.'p' .and. phid(2:2).eq.'P') .or.
     +                (phid(1:1).eq.'s' .and. phid(2:2).eq.'S')) then
                      indr = 4
                  else 
                      indr = 5
                  endif

              endif

              trefls = crustc(phase_t,dpaa,zoc,indr,typctl)

            endif

c
c     correction for surface multiples (e.g.: PnPn,...,PP,SS,P'P')
c
            if( (phid(1:1).eq.phid(2:2) .or. 
     +           phid(1:2).eq.phid(3:4)  ) ) then

                indr = 2
                zoc    = 9999.d0

                treflp = crustc(phase_t,dpaa,zoc,indr,typctl)
   
            endif

c
c      correction for converted surface multiples (e.g.: PnSn,...)
c
            conr = .false.
            if( (phid(1:1).eq.'P' .or. phid(2:2).eq.'P') .and.
     +          (phid(1:1).eq.'S' .or. phid(2:2).eq.'S') .and.
     +          (phid(3:3).eq.' ' .or. phid(3:3).eq.'g' .or.
     +           phid(3:3).eq.'b' .or. phid(3:3).eq.'n')) then
                   conr=.true.
                   phidr = phid(2:)
            endif

            if( (phid(1:1).eq.'P' .or. phid(3:3).eq.'P') .and.
     +          (phid(1:1).eq.'S' .or. phid(3:3).eq.'S') .and.
     +           phid(2:2).ne.'b' .and. phid(2:2).ne.'m' .and.
     +           phid(2:2).ne.'c' .and. phid(2:2).ne.'k' .and.
     +          (phid(2:2).eq.phid(4:4) .or. phid(2:2).eq.'g' .or.
     +           phid(2:2).eq.'n'                      )) then
                   conr=.true.
                   phidr = phid(3:)
            endif

            if(conr) then
                
                zor = 0.d0
                call tauget_ray(phidr,phase_t,dpaa,modn,zor,
     +                        del0,ttray,rayok)

                azi0 = azie(iev(i))
                if(dpa.lt.0.d0) azi0 = alpha2(azie(iev(i))-180.d0)

                inddel = 1
                call delazd(stala(iev(i)),stalo(iev(i)),azi0,del0,
     +                      inddel,elatc,elonc)

                indr = 3
                zoc    = 9999.d0
                treflp = crustc(phase_t,dpaa,zoc,indr,typctl)

            endif

            if (used(i)(5:5).ne.' ' .and.typctl.gt.6) then
                print *,'dirdel: ',phid,' azi ',azi0,' del ',del0
                print *,'  lat ',elatc,' lon ',elonc,' trefl ',trefl
            endif

            trefl = trefls + treflp

            if(dabs(trefl).ge.1.d-3) then
               used(i)(5:5) = 'R'
            else
               trefl = 0.d0
            endif

         endif

296      continue

         ttt(i) = tom + ttc(j) + dble(ecor) 
     +                + th + tcrust + trefl + statict

         dtt     = tt(i) - ttt(i)

         if (typctl.ge.8) then
c            print *,'i, ttt, t0, TT, ECOR, Height, Crust, Refl, '
c                    'Stat, DTT'  
             print *,i,sta(iev(i)),phid,ttt(i),tom,ttc(j),dble(ecor),
     +               th,tcrust,trefl,statict,dtt,touse(i)
         endif

         if(dtt.gt.100.d0 .and. phase(i)(1:2).eq.'P1' .and. 
     +      phid(1:3).eq.'Pdi') then
            phid='PKPdf'
            go to 297
         endif

         dtm = dtm0
         if(phase_t.eq.'P') dtm  = dtmp
         if(phase_t.eq.'S') dtm  = dtms
         dtm = dmin1(dtm,datmax0)

         if((dabs(dtt).le. dtm .or. fixinp ) .and. 
     +        touse(i)(1:1).eq.'T') then

           jj  = jj + 1
           jdt = jdt + 1
           rms1 = rms1 + dtt*dtt
           nrms1= nrms1+ 1

           dat(jj) = dtt

           if (lsmu ) then

              f1 = 0.d0
              if (phase_t.eq.'P') f1 = smpu
              if (phase_t.eq.'S') f1 = smsu

              fmus  = fmus  + q2(f1)
              fmud  = fmud  + q2(f1/dpa)
              fmuh  = fmuh  + q2(f1/dble(dtdh(j)))
              ismu = ismu + 1

           endif
             
c
c          setting uncertainty for this onset time
c
c          if (tresw) was set, tt2(i) may have a different 
c          value than tts(i) !
c

           if(dtt.gt.0.d0) then
               dats(jj) = tts(i)
           else
               dats(jj) = tt2(i)
           endif

c
c          (ierre.ne.0) > No ellipticity correction is available for this
c                         phase. We assume a larger data error!
c
           if(ierre.ne.0) then
              dats(jj)= dats(jj) + 0.5d0
              ierre = 0
           endif
c
c          we have to know the actually used uncertainty for later 
c          statistics
c
           ttu(i) = dats(jj)

           if(dabs(dtt/ttu(i)).gt.datmax) datmax=dabs(dtt/ttu(i))

           a(jj,1) = 1.d0
           a(jj,2) = dpa*deldla
           a(jj,3) = dpa*deldlo

           a(jj,4) = dble(dtdh(j))
           if(dabs(a(jj,4)).lt.1.d-5) then
              a(jj,4)=0.d0
           else
              nzo = nzo + 1
           endif

           used(i)(1:1) = 'T'

           dinv(i,1) = dble(jj)

           datla(i) = a(jj,2)
           datlo(i) = a(jj,3) 
           datho(i) = a(jj,4)

           if(typctl.gt.5) then
                print *,jj,' tt ',tt(i),a(jj,1),a(jj,2),a(jj,3),
     +                     a(jj,4),dat(jj),used(i)
           endif

         else if (dabs(dtt).le.dtm+dtdt .and. touse(i)(4:4).eq.'D') then

c
c          phase can later eventually be used for a travel-time-
c          difference observation.
c

           used(i)(1:1) = 't'

           datla(i) = dpa*deldla
           datlo(i) = dpa*deldlo
           datho(i) = dble(dtdh(j))

         endif

         phaseu(i) = phid

         if(touse(i)(3:3).eq.'S' .and. dabs(dtt).lt.dtmaxslow) then

           ddpa = p(i) - dabs(dpa - dph)
           if(dabs(ddpa).lt.dpam0 .or. fixinp) then
             jj  = jj + 1
             jpa = jpa + 1
             dat(jj)  = ddpa
             dats(jj) = ps(i)
             a(jj,1) = 0.d0
             a(jj,2) = dddp(j)*deldla
             a(jj,3) = dddp(j)*deldlo

             partabl = dddp(j)

             a(jj,4) = dpdh(j)
             if(dabs(a(jj,4)).lt.1.d-5) then 
                a(jj,4)=0.d0
             else
                nzo = nzo + 1
             endif

             used(i)(3:3) = 'S'

             dinv(i,3) = dble(jj)

c            if(typctl.gt.5) then
c                  print *,jj,' p ',p(i),a(jj,1),a(jj,2),a(jj,3),
c    +                    a(jj,4),dat(jj),used(i)
c            endif

         endif

        endif

        go to 298

       endif

297   continue

c
c     Try it with another phase-name from the same phase-type.
c

      if(single) go to 298

      call testphase (phid,icha)

      if(icha.eq.999) go to 298

      go to 295

298   if( ( touse(i)(2:2).eq.'A' .and. .not.aziini ) .and. 
     +    ( dabs(dtt).lt.dtmaxazib .or. 
     +      ( phid(1:1).eq.'L'.and.dabs(dtt).lt.dtmaxazil)) ) then

         if(used(i)(1:1).eq.'T' .and. (dpa.lt.0.d0 
     +                              .or.phase(i)(1:4).eq.'P3KP') ) then
           ddazi = alpha1(azi(i) - alpha2(baz(iev(i))-180.d0))
         else
           ddazi = alpha1(azi(i) - baz(iev(i)))
         endif

         if(dabs(ddazi).lt.dazim0  .or. fixinp) then
           jj = jj + 1
           jazi = jazi + 1
           dat(jj)  = ddazi
           dats(jj) = azis(i)
           a(jj,1) = 0.d0
           a(jj,2) = dazidla
           a(jj,3) = dazidlo
           a(jj,4) = 0.d0
           used(i)(2:2) = 'A'

           dinv(i,2) = dble(jj)

           if(typctl.gt.5) then
                  print *,jj,' baz ',azi(i),a(jj,1),a(jj,2),a(jj,3),
     +                     a(jj,4),dat(jj),used(i)
           endif
         endif

      endif

300   continue

      if(single) go to 302

      if(jazi .gt. jdt .and. jdt.le.3 .and. ifixaz.le.5) then
         ifixaz = ifixaz + 1
         dtmp   = dtmp * 1.5d0
         dtms   = dtms * 1.5d0
         dtm0   = dtmp + dtms
         dazim0 = dazim0 * 1.5d0
         dpam0  = dpam0  * 1.5d0
         go to 101
      endif

      if (nrms1.ge.2 .and. jj.ge.3) then

         rmso = dsqrt(rms1/dble(nrms1))

         if(ilastiter.eq.1 .and. iremo.lt.2 .and. 
     +      nrms1.gt.2 .and. jj.gt.3 .and. lastfixt) then

            rms0 = 1.5d0 * rmso / sdmeans

             if(datmax.ge.rms0) then
               datmax0 = rmso*2.d0
               iremo = iremo+1
               go to 101
            endif

          endif

         iteraz = 0 


      else 
         if(rmso.le.50.d0) then
            rmso = rmso*10.d0
         else 
            rmso = 9999.d0
         endif
         iteraz = iteraz + 1
      endif

      if(.not.diffflag) go to 302

c
c     Add possible travel-time difference(s) as additional 
c     condition(s) to the equation system to be solved.
c
c     Travel-time differences can only be used in the case that we 
c     have more than 2 different phase observations at one station.
c
      
      if(jdt.le.2) go to 302

      ndt = 0

      do 3011 i = 1,nobs-1

      if(touse(i)(4:4).ne.'D') go to 3011
      if((used(i)(1:1).ne.'T' .and. used(i)(1:1).ne.'t') ) go to 3011

      do 301 j = i+1,nobs

         if(sta(iev(i)).ne.sta(iev(j))) go to 301

         if(phaseu(i).eq.phaseu(j)) go to 301
         if (dabs(tt(i)-tt(j)).le.1.d-4)     go to 301
         if (dabs(ttt(i)-ttt(j)).le.1.d-4)   go to 301
         if(touse(j)(4:4).ne.'D' ) go to 301
         if((used(j)(1:1).ne.'T' .and. used(j)(1:1).ne.'t') ) go to 301

         if (((tt(j).gt.tt(i)) .and. (ttt(j).lt.ttt(i))) .or.
     +       ((tt(i).gt.tt(j)) .and. (ttt(i).lt.ttt(j)))   )  go to 301


         dtt = (tt(j) - ttt(j)) - (tt(i) - ttt(i))

         if((dabs(dtt).le.dtm0 .and. dabs(dtt).gt.0.d0).or. fixinp) then

           jj = jj + 1

           ndt = ndt + 1
           idtu(ndt) = j*i + j+i

           dat(jj) = dtt

           dats(jj)= dpythag(ttu(i),ttu(j))

           a(jj,1) = 0.d0
           a(jj,2) = datla(j) - datla(i)
           a(jj,3) = datlo(j) - datlo(i)

           a(jj,4) = datho(j) - datho(i)
           if(dabs(a(jj,4)).lt.1.d-5) then
              a(jj,4)=0.d0
           else
              nzo = nzo + 1
           endif

           dinv(i,4) = dble(jj) + dble(ndt)*1.D-3

c          if(typctl.gt.5) then
c                 print *,jj,' dt ',a(jj,1),a(jj,2),a(jj,3),
c    +                     a(jj,4),dat(jj),used(i),used(j)
c          endif

          endif

301   continue
3011  continue

c
c     Everything is ready for a next or  a 'final' inversion
c
c     hyposat_gmi will do it
c

302   in = jj

      if(czo.eq.'D') then

c        print *,'czo, nobs, jdt ', czo,nobs,jdt
        if(nobs.lt.4 .or. jdt.lt.4) then
           if(zo/(deg2rad*radloc(stala(istatmin),1)).lt.0.25d0 .and. 
     +        iter.gt.2) then
             if( zoflag ) go to 9998
             zo = 0.1d0
             czo = 'B'
             if(typctl.gt.0) then
                print *,'(1) No depth resolution, fixed at',zo,' [km]'
             endif

             if(output) then
                write(11,'(/''No resolution for depth, fixed at''
     +                      ,f7.2)') zo
             endif

             go to 101
           endif
        endif

        im = 4

        if (nzo .le. 0 .or. in.le. 3) then

           if( zoflag ) go to 9998

           czo = 'B'
           im  = 3

           if(typctl.gt.0) then
              print *,'No depth determination possible, fixed at',
     +                zo,' [km]'
           endif
           if(output) then
              write(11,'(/''No depth determination possible, '',
     +                    ''fixed at'',f7.2)') zo
           endif

        endif

      else if (czo.eq.'F' .or. czo.eq.'B') then

        im = 3

      endif

      if(in.gt.mread2 .and. .not.single) then

         print *, 'Inversion matrix: ',in,' (data) ',mread2,
     +            ' (wrong dimension of Jacobian)'
         go to 9999

      endif

      if(in.le.1 .and. .not.single) then

         if(in0sw .ge. 6) then
            print*,'Inversion failed'
            print*,'No data fit within limits of initial model!'
            go to 9999
         endif

         if(plflag .and. in0sw.eq.5) then

           call plane(stala,stalo,iev,tt,nobs,azim,dazi,
     +                ray,dray,phipl,touse,phase,jref,typctl)

           if(jref.gt.0 .and. dazi.lt.90.d0 .and. dray.lt.4.d0) then

             phase_t = ' '
             itray = 1
             rayok = .false.
3021         call tauget_ray(phipl,phase_t,ray,modn,zo,
     +                        ddel,ttray,rayok)

             if(rayok) then

               inddel = 1

c               print *,' ---> [delazd] B',stala(iev(jref)),
c     +              stalo(iev(jref)),azim,ddel,inddel,
c     +              elatmg,elonm
               
               call delazd(stala(iev(jref)),stalo(iev(jref)),azim,
     +                     ddel,inddel,elatmg,elonm)

               
c               print *,' ---> [delazd] A',stala(iev(jref)),
c     +              stalo(iev(jref)),azim,ddel,inddel,
c     +              elatmg,elonm


               elatm = convlat(elatmg,1)

               if(typctl.gt.0) then
                  if(ddel.gt.0d0) then
                    print *,'Epicenter set from station ',
     +               sta(iev(jref)),'after plane wave fit: backazimuth',
     +               azim,' deg, delta',ddel,' deg'
                  else
                    print *,'Epicenter set to station ',sta(iev(jref))
                  endif
               endif

               sdlatg = 45.d0/ray
               sdlon  = 90.d0/ray

               tome = tt(jref) + timemin - ttray

               if(output) then
                  if(ddel.gt.0d0) then
                       write(11,'(''Epicenter set from station '',a8)')
     +                   sta(iev(jref))
                  else
                       write(11,'(''Epicenter set to station '',a8,
     +                   '' deg, delta'',f7.2,'' deg'')') sta(iev(jref))
                  endif
                  write(11,'(''Epicenter lat:'',f9.3,'' [deg]'')')  
     +              elatmg,sdlatg
                  write(11,'(''Epicenter lon:'',f9.3,'' [deg]''/)') 
     +                     elonm,sdlon
                  write(11,'(''Source time set to: '',f15.2)') tome
               endif

               in0sw = in0sw + 1
               go to 101

             else
               
               if(itray.lt.2) then

                  itray = itray + 1
                  ray = ray - dray
                  go to 3021

               endif

             endif

           endif

         endif

         in0sw = in0sw + 1

         if(iter.ge.3) then

           if(iter.gt.mosci) then
               ilas = 1
           else
               ilas = iter-2
           endif

           tom     = dtos(ilas)
           tome    = tom + timemin
           elonm   = dloos(ilas)
           elatm   = dlaos(ilas)
           elatmr  = deg2rad*elatm
           elatmg  = convlat(elatm,2)
           coelatm = 90.d0 - elatm
           coelatmr= deg2rad*coelatm

         endif

         dtmp   = dtmp * 2.d0
         dtms   = dtms * 2.d0
         dtm0   = dtms
         dazim0 = dazim0 * 2.0d0
         dpam0  = dpam0 * 2.0d0

         if(typctl.ge.4) then
            print *, 'Inversion matrix error: in = ',in
            print *, '(too less data to invert!)' 
            print *, 'Time boundaries changed'
         endif

         go to 101

      endif

c
c     special case: only one, single array observation
c
      if(single .and. in.eq.3) then

        dtmflag = .true.

        ddel = dat(2) / partabl

         if((p(1).lt.pdif .or. p(1).gt.9.d0) .and. 
     +                   dabs(ddel).gt.1.d0    )    ddel=ddel/2.d0
         if(p(1).gt.9.d0 .and. ddel.gt.1.d0     )    ddel=1.d0
         if(p(1).gt.9.d0 .and. ddel.lt.-1.d0    )    ddel=-1.d0

        deln = del(1) + ddel 

        inddel = 1
        call delazd(stala(1),stalo(1),azi(1),deln,inddel,elat1,elon1)

        elatm1 = convlat(elat1,1)

        r(1) = dat(1)
        r(2) = elatm1 - elatm
        r(3) = elon1 - elonm
        r(4) = 0.d0

        var(1) = dats(1)

        dvar = dabs(dats(2) / partabl)

        ddel = ddel + dvar
        call delazd(stala(1),stalo(1),azi(1),ddel,inddel,elat1,elon1)
        ddel = ddel - 2.d0*dvar
        call delazd(stala(1),stalo(1),azi(1),ddel,inddel,elat2,elon2)
        vara = (elat2-elat1)/2.d0
        varb = (elon2-elon1)/2.d0

        ddel = ddel + dvar
        aziv1 = azi(1) + dats(3)
        call delazd(stala(1),stalo(1),aziv1,ddel,inddel,elat1,elon1)
        aziv1 = azi(1) - 2.d0*dats(3)
        call delazd(stala(1),stalo(1),aziv1,ddel,inddel,elat2,elon2)

        var(2) = dpythag((elat2-elat1)/2.d0,vara)
        var(3) = dpythag((elon2-elon1)/2.d0,varb)

        var(4) = 0.d0

        res(1) = 0.d0
        res(2) = 0.d0
        res(3) = 0.d0
        res(4) = 0.d0

        go to 304

      else if (single .and. in.lt.3) then

           if(insar.le.15) then

              ddel = deln

              if(insar.le.0) then 

                 phase(1) = 'P'
                 ddel = 50.d0

              else

                if(p(1).le.pdif) then

                   if(insar.ge.7) go to 3028

                   ddel = 148.d0

                   if(phid(1:2).eq.'P ')   phase(1) = 'PKPab'
                   if(phid(1:3).eq.'Pdif') phase(1) = 'PKPab'
                   if(phid.eq.'PKPab')     phase(1) = 'PKPbc'
                   if(phid.eq.'PKPdif')    phase(1) = 'PKPbc'
                   if(phid.eq.'PKPbc')     phase(1) = 'PKPdf'
                   if(phid.eq.'PKPdf')     phase(1) = 'PKiKP'
              
                else if(p(1).gt.9.d0) then

                   if(insar.ge.4) go to 3028

                   if(phid(1:2).eq.'P ' .or. phid(1:2).eq.'P1')  then
                      phase(1) = 'Pn'
                      ddel = 10.d0
                   else if(phid.eq.'Pn') then
                      phase(1) = 'Pb'
                      ddel = 2.d0
                   else if(phid.eq.'Pb') then
                      phase(1) = 'Pg'
                      ddel = 1.d0
                   endif

                endif

              endif

              insar = insar + 1
              inddel = 1
              call delazd(stala(1),stalo(1),azi(1),ddel,inddel,
     +              elatmg,elonm)
              elatm = convlat(elatmg,1)

              phaseu(1) = phase(1)

              go to 100

           endif

3028           print*,'Single phase case, but no inversion possible'
           if(output) then
              write(11,*)'Single phase case, but no inversion ',
     +                         'possible'
           endif

           if(p(1).ge. pmoh) then
              print*, 'due to missing direct',
     +                     ' crustal phase matching theoretical ray',
     +               ' parameter in chosen model.'
              if(output) then
                 write(11,*) 'due to missing direct',
     +                     ' crustal phase matching theoretical ray',
     +               ' parameter in chosen model.'
              endif
           else if(p(1).gt. 10.1d0 .and. p(1).lt. pmoh) then
              print*, 'due to missing direct upper mantle',
     +                     ' crustal phase matching theoretical ray',
     +               ' parameter in chosen model.'
              if(output) then
                 write(11,*) 'due to missing direct upper mantle',
     +                     ' crustal phase matching theoretical ray',
     +               ' parameter in chosen model.'
              endif
           else
              print *, 'due to core phase travel-time-curve ',
     +                     'triplication!'
              if(output) then
                 write(11,*) 'due to core phase travel-time-curve ',
     +                     'triplication!'
              endif
           endif
           go to 9999
           
      endif

      in0sw = 0

      if (iellipi.eq.1 .and. .not.iellip) iellip=.true.

      ilastfixi = 0

c      print *,'--->  [gmi]',in,im,nq, iter
c      do i=1,in,1
c         print *,(a(i,j),j=1,im),dat(i),dats(i)
c      print *,' '
c      enddo
c      print *,' '
      
303   continue

c
      call hyposat_gmi(in,im,nq,ierr,typctlm)

c       do i=1,in,1
c          print *,(a(i,j),j=1,im)
c      print *,' '
c       enddo
c       print *,' '
      
      if(ierr.ne.0) then
         if(output) then
            write(11,*) 'GMI failed: no new solution!'
         endif
         print*,'GMI failed: no new solution!'
         go to 9999
      endif

304   continue
      
      if(typctl.gt.4) then
          print *, 'iter,in,im,nq',iter,in,im,nq
          print *,r(1),var(1)
          print *,r(2),var(2)
          print *,r(3),var(3)
          print *,r(4),var(4)
          do 305 j=1,jj
          f1 = dat(j)-res(j)
305       print *,j,dat(j),res(j),f1
      endif

c
c     estimating the new hypocenter
c
         
c
c     the new source depth
c
      if(im.eq.4) then
        ar4 = dabs(r(4))
        if(ar4.gt.200.d0 ) then
           r40 = r(4)*0.1875d0
        else if(ar4.le.200.d0 .and. ar4.gt.100.d0 ) then
           r40 = r(4)*0.375d0
        else if(ar4.le.100.d0 .and. ar4.gt.40.d0) then
           r40 = r(4)*0.75d0
        else
           r40 = r(4)
        endif

        zo = zo + r40*dchang

        if(var(4).gt.1.d-5) then
           rs(4) = var(4)
        else
           rs(4) = sdzo
        endif

        if(zo.lt.0.1d0) then
           r40 = 0.1d0 - zo
           rs(4) = dpythag(rs(4),r(4)-r40)
           zo = 0.1d0
        endif

        if(zo.ge.800.d0) then
           r40 = 799.d0 - zo
           rs(4) = dpythag(rs(4),r(4)-r40)
           zo = 799.d0
        endif

        if(zo.lt.depthmin) then
           r40 = zo - depthmin
           rs(4) = dpythag(rs(4),r(4)-r40)
           zo = depthmin
        endif

        if(zo.gt.depthmax) then
           r40 = zo - depthmax
           rs(4) = dpythag(rs(4),r(4)-r40)
           zo = depthmax
        endif

        if(rs(4).gt.250.d0) rs(4) = 250.d0

      endif

c
c     the new source time
c
      rsi = 1.d0
      ar1 = dabs(r(1))
      if(r(1).lt.0.d0) rsi = -1.d0
      if(ar1.le.60.d0 .or. single) then
         tome = tome   + r(1) * dchang
      else if(ar1.le.120.d0.and.ar1.gt.60.d0) then
         tome = tome + r(1)/2.d0
      else
         tome = tome + rsi*120.d0
      endif

      tom    = tome - timemin

      if(var(1).gt. 1.d-3) rs(1) = var(1)
      if(rs(1).gt.250.d0) rs(1) = 250.d0

c
c     save the old epicenter solution 
c
      elatmo = elatmg
      elonmo = elonm

c
c     the new source latitude
c
      rsi = 1.d0
      ar2 = dabs(r(2))
      if(r(2).lt.0.d0) rsi = -1.d0
      if(ar2.le.5.d0) then
         elatm = elatm + r(2)*dchang
      else if(ar2.le.10.d0.and.ar2.gt.5.d0) then
         elatm = elatm + r(2)/2.d0
      else
         elatm = elatm + rsi*10.d0
      endif

      ilon = 0
      if(elatm.gt. 90.d0) then
         elatm = 180.d0 - elatm
         ilon  = 1
      else if(elatm.lt.-90.d0) then
         elatm = -(elatm + 180.d0)
         ilon  = 1
      endif

      elatmr  = deg2rad*elatm
      elatmg  = convlat(elatm,2)
      sdlatg  = var(2) /( eps*q2(dcos(elatmr))
     +                            +q2(dsin(elatmr)) )

      coelatm = 90.d0 - elatm
      coelatmr= deg2rad*coelatm

      if(var(2).gt. 1.d-5) rs(2) = var(2)
      if(rs(2) .gt.180.d0)  rs(2) = 180.d0

c
c     the new source longitude
c
      rsi = 1.d0
      ar3 = dabs(r(3))
      if(r(3).lt.0.d0) rsi = -1.d0
      if(ar3.lt.5.d0) then
         elonm   = elonm  + r(3)*dchang
      else if(ar3.le.10.d0.and.ar3.gt.5.d0) then
         elonm   = elonm  + r(3)/2.d0
      else
         elonm   = elonm  + rsi*10.d0
      endif

      elonm = alpha1(elonm)

      if(ilon.eq.1) elonm = alpha1(elonm+180.d0)

      if(var(3).gt. 1.d-5) rs(3) = var(3)
      if(rs(3).gt.180.d0) rs(3) = 180.d0

      last = .false.

      if(iter.eq.maxiter) then
         dtmflag = .true.
         if(zo.le.0.1d0) zo = 0.d0
      else if(iter.gt.maxiter) then
         print *,'Location stopped, more than ',maxiter,
     +           ' iterations'
         if(output) then
            write(11,'(/''Location stopped, more than '',i4,
     +              '' iterations'')') maxiter
         endif

         mosci2 = 1
         go to 399
      endif

      if(iter.le.mosci) then

        dtos(iter)  = tom
        dlaos(iter) = elatm
        dloos(iter) = elonm
        dlo1os(iter)= dcos(deg2rad*elonm)
        dlo2os(iter)= dsin(deg2rad*elonm)
        dzoos(iter) = zo
        rtos(iter)  = dabs(r(1)) + var(1)
        rlaos(iter) = dabs(r(2)) + var(2)
        rloos(iter) = dabs(r(3)) + var(3)
        rzos(iter)  = dabs(r(4)) + var(4)

      endif
 
c
c     We will check if the new solution is close to the former solution.
c     If CHECK [km] is smaller than SETCHECK [km] we will stop.
c

c     The change in the horizontal plane (in the epicenter)

      dk = 0.d0
      call depi (elatmg,elonm,elatmo,elonmo,del3,dk,ep2,ep1,d2km)

c     The change in source time is compensated eventually by a
c     change in depth

      dtokm = 0.d0

      if(im.eq.4) then

         dtokm = r(1)*6.d0
         if (zo.gt.20.d0) then
            dtokm = r(1)*7.d0
         else if (zo.gt.30.d0) then
            dtokm = r(1)*8.d0
         else if (zo.gt.260.d0) then
            dtokm = r(1)*9.d0
         else if (zo.gt.450.d0)  then
               dtokm = r(1)*10.d0
         else if (zo.gt.660.d0)  then
            dtokm = r(1)*11.d0
         endif

         if(var(4).gt.1.d-5) dtokm = dtokm - r(4)

      endif

      check = dpythag(dtokm,dk)

      call depi (elatmg,elonm,stala(istatmin),stalo(istatmin),del3,
     +                    dk,ep2,ep1,d2km)

      ilastiter = 0

      if (check.le.setcheck .or. (check.le.disper*dk .and. 
     +      (iteraz.ge.5 .or. imaxiter.ge.5 .or. 
     +       dble(maxiter)/dble(iter).lt.1.3d0   )) ) then

         ilastiter = 1

         if((dtmflag .or. .not.lastfixt .or. in.le.im) .and. 
     +      rmsold/rmso.lt.1.2d0) then

            if(ilastfixi.eq.0 .and. lastfixi .and. 
     +         in/im.gt.3                          ) then

               ilastfixi = 1

               if(thrfixi0.le.0.d0) then
                  thrfixi = 1.d0 / dble(in)
               else
                  thrfixi = thrfixi0
               endif

               infind = 0
               do 365 j=1,in

                  if(dinf(j)/dinfm .lt. thrfixi) then

                     if(typctl.gt.8) print *,j,dinf(j),dinf(j)/dinfm

                     infind = infind + 1
                     do 360 j3=1,im
                        a(j,j3) = 0.d0
360                  continue

                     do 363 j3=1,nobs
                        do 361 j4=1,4
                           if(int(dinv(j3,j4)).eq.j) then
                              if(j4.le.3) then
                               used(j3)(j4:j4) = ' '
                              else
                               idum = int( (dinv(j3,4)-
     +                                int(dinv(j3,4)))*1000.d0+0.01d0)
                               idtu(idum) = 0
                              endif
                              dinv(j3,j4)   = 0.d0
                              go to 365
                           endif
361                     continue
363                  continue
                  endif

365            continue

               if(infind.gt.0) go to 303

            endif

            go to 400

         else

            if(zo.le.0.1d0) zo=0.d0

            if(var(4).le.0.d0 .and. czo.eq.'D') then

               czo = 'B'
               if(typctl.gt.0) then
                  print *,'(3) No depth resolution, fixed at',zo,' [km]'
               endif

               if(output) then
                  write(11,'(/''No resolution for depth, fixed at''
     +                        ,f7.2)') zo
               endif

            endif

            setcheck2 = check*1.1d0

            dchang = dchang0

            go to 390

         endif
      endif

      dtmflag = .false.
      direct  = .false.

      if( ifixto.le.5  .and. (var(1).le.0.d0 .or. 
     +    (jdt.lt.in/3 .and. jdt.lt.nstat/2)) ) then

         ifixto = ifixto + 1

         rzo = sngl(zo)
         rdel = sngl(del3)
         nphas = 0
         call tauget_mod(rzo,rdel,nphas,phcd,ttc,dtdd,
     +                           dtdh,dddp,modnam)

         tom  = - dble(ttc(1))
         tome = timemin + tom

         rs(1) = rs(1) + (-tom)

         dtp   = dtp  * 1.2d0
         dts   = dts  * 1.2d0
         dtm2  = dtp + dts

         setcheck = setcheck*1.2d0
         setcheck2= setcheck*10.d0

         go to 101

      endif

      if(iter.gt.mosci) then
        moscil = mosci
      else
        moscil = iter-1
      endif

      nextiter1 = 0

c
c     check for oscillating solutions
c

      mosci2 = 0

      if(moscil.gt.2) then
         do 370 i = moscil,1,-1

           if( dabs(dzoos(i)-zo).le.rminh            .and.
     +         dabs(dtos(i)-tom).le.rmint            .and.
     +         dabs(dlaos(i)-elatm).le.rming         .and.
     +         dabs(alpha1(dloos(i)-elonm)).le.rming        ) then
               mosci2 = i
c              print *,dzoos(i)-zo,rminh,dtos(i)-tom,rmint,
c    +          dlaos(i)-elatm,dloos(i)-elonm,rming
               go to 371
           endif

370      continue
      endif

371   continue

c      Print *, 'MOSCI2: ', mosci2,'mosci: ',mosci,' iter: ', iter

      if(mosci2.gt.0) then
c
c     we have to calculate a new initial solution from all
c     oscillating solutions!
c

         nextiter = nextiter + 1
         if(nextiter.gt.8) then
            if(output) then
               write(11,'(/''Oscillating solution: after'',i3,
     +                   '' iterations stopped!'')') iter
            endif
            go to 399
         endif

         nextiter1 = 1

         if (czo.eq.'D') then

             zo  = dmean(dzoos,moscil,mosci2) 
             rs(4) = ddmax(rzos,moscil,mosci2)
             if(rs(4).le.1.d-5) rs(4) = sdzo

             if(zo.le.0.1d0 .or. nextiter.gt.moscil*2/3) then
                if(iterz.eq.1 .and. zoflag) go to 9998
                iterz = iterz + 1
                czo = 'B'
                var(4) = 0.d0
                rs(4)  = 1.d0
                if(typctl.gt.0) then
                   print *,'(2) No depth resolution, fixed at',zo,
     +                     ' [km]'
                endif

                call findrange(zmin,zmax1,dzoos,moscil,mosci2,1)

                if(output) then

                   write(11,'(/''No resolution for depth, '',
     +                         ''oscillating between:'',f6.1,
     +                         '' and'',f6.1,'' [km]''/
     +               ''depth fixed at:'',f6.1,'' km'')') zmin,zmax1,zo
                endif

             endif

         endif

         tom = dmean(dtos,moscil,mosci2)
         tome = tom  + timemin

         rs(1) = ddmax(rtos,moscil,mosci2)
         if(rs(1).lt.1.d-3) rs(1) = sdto

         elatm = dmean(dlaos,moscil,mosci2)
         elatmr= deg2rad*elatm
         elatmg  = convlat(elatm,2)
         coelatm = 90.d0 - elatm
         coelatmr= deg2rad*coelatm
         rs(2) = ddmax(rlaos,moscil,mosci2)
         if(rs(2).lt.1.d-5) rs(2) = sdlat

         elonm1= dmean(dlo1os,moscil,mosci2)
         elonm2= dmean(dlo2os,moscil,mosci2)
         elonm = rad2deg*datan2(elonm2,elonm1)

         rs(3) = ddmax(rloos,moscil,mosci2)
         if(rs(3).lt.1.d-5) rs(3) = sdlon

         if (nrms1.gt.5) then

             dtp = dtp  * 0.9d0
             dtp = dmin1(dtp,rmso*2.d0)
             if(dtp.lt.3.d0) dtp = 3.d0


             dts = dts  * 0.9d0
             dts = dmin1(dts,rmso*4.d0)
             if(dts.lt.6.d0) dts = 6.d0

             dtm2  = dtp + dts

         endif

         if(nobs.gt.1) then

            dazim1 = dazim1*0.9d0
            if(dazim1.lt.3.d0) dazim1 = 3.d0

               dpam1  = dpam1 *0.9d0
               if(dpam1.lt.1.5d0) dpam1 = 1.0d0

            endif

            rminh = rminh*1.2d0
            rming = rming*1.2d0
            rmint = rmint*1.2d0

            disper = disper*1.1d0
            if(disper.gt.0.05d0) disper=0.05d0

            setcheck = setcheck*1.5d0
            setcheck2= setcheck*10.d0

            dchang = dchang * 0.8d0

            direct = .true.

         endif
         
         if(var(4).gt.4.d0*zo .and. czo.eq.'D' .and. var(4).gt.30.d0 
     +      .and. iter.gt.5 ) then

            if(zoflag .and. iterz.ge.6) go to 9998

            czo = 'B'
            rs(4) = 1.d0
            var(4) = 0.d0
            if(mosci2.gt.0)then
               zo  = dmean(dzoos,moscil,mosci2)
            else
               zo  = dmean(dzoos,moscil,1)
            endif

            iterz = iterz + 1

            if(typctl.gt.0) then
               print *,'Bad resolution for depth; depth fixed at',zo
            endif

            if(output) then
               write(11,'(/''Bad resolution for depth; depth fixed'')')
            endif

            direct = .true.

            go to 380

         endif

         if(czo.eq.'D' .and. nextiter1.eq.0 .and. iterz.le.3 .and.
     +      moscil.ge.1) then

c
c           check for oscillation in the solutions for the focal depth
c
            mosci3 = 0

            do 377 i = moscil,1,-1

              if(dabs(dzoos(i)-zo).le.rminh) then
                 mosci3 = i
                 go to 378
              endif

377         continue

378         if((mosci3.ne.0      .and. mosci3.ne.moscil) .or. 
     +         (mosci3.eq.moscil .and. zo.le.0.1d0)         ) then
c
c           we have to calculate a new initial value for the depth and 
c           fix it 
c

              rs(4)  = 1.d0
              var(4) = 0.d0
              if (zo.gt.0.1d0) zo  = dmean(dzoos,moscil,mosci3)
              czo = 'B'

              iterz = iterz + 1

              if(typctl.gt.0) then
                print *,'Oscillating solution, ',
     +                  'depth fixed at: ',zo,' km'
              endif

              call findrange(zmin,zmax1,dzoos,moscil,mosci3,1)

              if(output) then
                 write(11,'(/''Oscillating solution between:'',f6.1,
     +                       '' and'',f6.1,'' [km] depth''/
     +             ''Depth fixed at:'',f6.1,'' km'')') zmin,zmax1,zo
            endif

            direct = .true.

         endif

      endif

380   continue
        
      if (iter.gt.mosci) then
        do 381 i = 1,mosci-1

        i2 = i + 1

        dtos(i)  = dtos(i2)
        dlaos(i) = dlaos(i2)
        dloos(i) = dloos(i2)
        dlo1os(i)= dlo1os(i2)
        dlo2os(i)= dlo2os(i2)
        dzoos(i) = dzoos(i2)
        rtos(i)  = rtos(i2)
        rlaos(i) = rlaos(i2)
        rloos(i) = rloos(i2)
        rzos(i)  = rzos(i2)

381     continue

        dtos(mosci)  = tom
        dlaos(mosci) = elatm
        dloos(mosci) = elonm
        dlo1os(mosci)= dcos(deg2rad*elonm)
        dlo2os(mosci)= dsin(deg2rad*elonm)
        dzoos(mosci) = zo
        rtos(mosci)  = dabs(r(1)) + var(1)
        rlaos(mosci) = dabs(r(2)) + var(2)
        rloos(mosci) = dabs(r(3)) + var(3)
        rzos(mosci)  = dabs(r(4)) + var(4)

      endif

      if(direct) go to 100

390   if(typctl.ge.4) then
         print*,'Iteration: ',iter,'   # of def.: ',in
         print*,'New source time  : ',tome,' +/- ',var(1)
         print*,'New epicenter lat: ',elatmg,' +/- ',sdlatg
         print*,'New epicenter lon: ',elonm,' +/- ',var(3)
         print*,'New source depth : ',zo,' +/- ',var(4)
      endif

c

      if(iter.gt.nint(maxiter*0.75) .and. imaxiter.lt.5 .and.
     +   ilastiter.eq.0) then

c
c     If we were coming close to the end of all iterations,
c     let's try it with a mean solution of the last 
c     4 solutions as new initial solution.
c

         imaxiter = imaxiter + 1
         maxiter  = maxiter + nint(maxiter*0.25/imaxiter)

         if (czo.eq.'D') then
 
            zo  = dmean(dzoos,mosci,1)
            rs(4) = ddmax(rzos,mosci,1)
            if(rs(4).le.1.d-5) rs(4) = sdzo

         endif
 
         tom = dmean(dtos,mosci,1)
         tome = tom  + timemin
         rs(1) = ddmax(rtos,mosci,1)
         if(rs(1).lt.1.d-3) rs(1) = sdto
 
         elatm = dmean(dlaos,mosci,1)
         elatmr= deg2rad*elatm
         elatmg  = convlat(elatm,2)
         coelatm = 90.d0 - elatm
         coelatmr= deg2rad*coelatm
         rs(2) = ddmax(rlaos,mosci,1)
         if(rs(2).eq.1.d-5) rs(2) = sdlat
 
         elonm1= dmean(dlo1os,mosci,1)
         elonm2= dmean(dlo2os,mosci,1)
         elonm = rad2deg*datan2(elonm2,elonm1)

         rs(3) = ddmax(rloos,mosci,1)
         if(rs(3).eq.1.d-5) rs(3) = sdlon

         dazim1 = dazim1*2.d0
         if(dazim1.gt.90.d0) dazim1 = 90.d0

         dpam1  = dpam1 *2.d0
         if(dpam1.gt.15.d0) dpam1 = 15.d0

          dtp   = dtp * 2.d0
          dts   = dts * 2.d0
         dtm2  = dtp + dts

         setcheck = setcheck*1.5d0
         setcheck2= setcheck*15.d0
         disper   = disper*2.d0
         if(disper.gt.0.05d0) disper=0.05d0

         rminh = rminh*1.5d0
         rming = rming*1.5d0
         rmint = rmint*1.5d0

      endif

      go to 100

399   continue

      last = .true.

      call findrange(tomin,tomax,dtos,mosci,mosci2,1)
      call findrange(dlamin,dlamax,dlaos,mosci,mosci2,1)
      call findrange(dlomin,dlomax,dloos,mosci,mosci2,2)
      call findrange(zmin,zmax1,dzoos,mosci,mosci2,1)

      if(output) then
         write(11,'(''Rel. source time between'',f9.2,'' and'',
     +           f9.2,'' [s]'')') tomin,tomax
      endif

      flamin = convlat(dlamin,2)
      flamax = convlat(dlamax,2)
      if(output) then
         write(11,'(''Latitude         between'',f8.2,'' and'',
     +           f8.2,'' [deg]'')') flamin,flamax

         write(11,'(''Longitude        between'',f8.2,'' and'',
     +           f8.2,'' [deg]'')') dlomin,dlomax

         if(czo.eq.'D') then
            write(11,'(''Depth            between'',f7.1,''  and'',
     +              f7.1,''  [km]''/)') zmin,zmax1
         endif

            write(11,'(//''Following the last (must not be the '',
     +                  ''best!) solution:''/)') 

      endif

400   continue
      if(output) then
         write(11,'(/''Iterations        :'',i5)') iter
      endif

      ibad = 0

      if(rmso .gt. 50.d0) ibad = ibad + 1

      if(infind.gt.0) in = in - infind
      if(output) then
         write(11,'(''Number of defining:'',i5)') in
         if(iloc) then
            if( (mod2flag .or. mod3flag .or. mod4flag) .and.
     +          (imodn(2).gt.1 .or. imodn(3).gt.1 .or. imodn(4).gt.1)
     +           ) then
                write(11,'(''First reference models  : '',a,'' and '',
     +                a)') trim(filloc),modnam
                if(mod2flag .and. imodn(2).gt.1) then
                   write(11,'(''Second reference model  : '',a)') 
     +                  modnam2
                endif
                if(mod3flag .and. imodn(3).gt.1) then
                   write(11,'(''Third reference model   : '',a)') 
     +                  modnam3
                endif
                if(mod4flag .and. imodn(4).gt.1) then
                   write(11,'(''Fourth reference model  : '',a)') 
     +                  modnam4
                endif
            else
                write(11,'(''Reference models  : '',a,'' and '',a)') 
     +                trim(filloc),modnam
            endif
         else
            if( (mod2flag .or.  mod3flag .or. mod4flag) .and.
     +          (imodn(2).gt.1 .or. imodn(3).gt.1 .or. imodn(4).gt.1)
     +           ) then
              write(11,'(   ''Main reference model    : '',a)') modnam
              if(mod2flag .and. imodn(2).gt.1) 
     +           write(11,'(''Second reference model  : '',a)') modnam2
              if(mod3flag .and. imodn(3).gt.1) 
     +           write(11,'(''Third reference model   : '',a)') modnam3
              if(mod4flag .and. imodn(4).gt.1) 
     +           write(11,'(''Fourth reference model  : '',a)') modnam4
            else
              write(11,'(''Reference model   : '',a)') modnam
            endif
         endif
      endif

      call fetoh(tome,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

c
      if(lsmu .and. ismu.gt.0) then

         f1     = dsqrt (fmus / ismu)
         var2(1) = fchi1*dpythag(var(1),f1)

         fmud = dsqrt (fmud / ismu)

         f1     = var(2) /( eps*q2(dcos(elatmr))
     +                         +q2(dsin(elatmr)) ) 
         sdlatg =  fchi1 * dpythag(f1,fmud) 

         f1 = fmud/dcos(deg2rad*elatmg)
         var2(3) = fchi1 * dpythag(var(3),f1)

         fmuh = dsqrt (fmuh / ismu)
         var2(4) = fchi1 * dpythag(var(4),fmuh)

      else

         var2(1) = fchi1 * var(1)

         f1      = var(2) /( eps*q2(dcos(elatmr))
     +                         +q2(dsin(elatmr)) ) 
         sdlatg  = fchi1 * f1
         var2(3) = fchi1 * var(3)
         var2(4) = fchi1 * var(4)

      endif

      if(output) then

         write(11,'(/''The new source parameters:'')')

         write(11,'(/''Confidence level of given uncertainties:'',
     +         f7.2,'' %'')') confl

         if(kmout) then
            if(disminst.gt.0.d0 .and. dismaxst.ge.21000.d0) then
               write(11,'(/''Location for observations at distances '',
     +               ''larger than'',f8.1,'' [km]'')') disminst
            endif
            if(disminst.le.0.d0 .and. dismaxst.lt.21000.d0) then
               write(11,'(/''Location for observations at distances '',
     +               ''smaller than'',f8.1,'' [km]'')') dismaxst
            endif
            if(disminst.gt.0.d0 .and. dismaxst.lt.21000.d0) then
               write(11,'(/''Location for observations at distances '',
     +               ''between'',f8.1,'' and'',f8.1,'' [km]'')') 
     +               disminst,dismaxst
            endif
         else
            if(disminst.gt.0.d0 .and. dismaxst.ge.180.d0) then
               write(11,'(/''Location for observations at distances '',
     +               ''larger than'',f6.1,'' [deg]'')') disminst
            endif
            if(disminst.le.0.d0 .and. dismaxst.lt.180.d0) then
               write(11,'(/''Location for observations at distances '',
     +               ''smaller than'',f6.1,'' [deg]'')') dismaxst
            endif
            if(disminst.gt.0.d0 .and. dismaxst.lt.180.d0) then
               write(11,'(/''Location for observations at distances '',
     +               ''between'',f6.1,'' and'',f6.1,'' [deg]'')') 
     +               disminst,dismaxst
            endif
         endif

         write(11,'(/''Source time  :'',i5,4i3.2,f7.3,'' +/- '',
     +           f8.3,'' [s]'')')  yy,mon,dd,hh,mi,sec,var2(1)
         write(11,'(''        or'',12x,f16.3,'' +/- '',f8.3,
     +           '' [s]'')') tome,var2(1)

         isec1 = nint(sec*1000)
         isec  = isec1/1000
         msec  = isec1-isec*1000
         write(11,'(''        or'',7x,i4,''-'',i3.3,'':'',
     +         3(i2.2,''.''),i3.3,'' +/- '',f8.3,'' [s]''//)') 
     +         yy,idoy,hh,mi,isec,msec,var2(1)
         
      endif

      if(var2(1).gt.zo/6.d0 .and. czo.eq.'D') ibad = ibad + 1

      if(sdlatg.lt.45.d0) then
         if(output) then
            write(11,'(''Epicenter lat:'',14x,f10.4,'' +/- '',f8.4,
     +              '' [deg]'')')  elatmg,sdlatg
         endif
      else
         ibad = ibad + 1
         if(sdlatg.gt.180.d0) sdlatg = 180.d0
         if(output) then
            write(11,'(''Epicenter lat:'',14x,f10.4,'' +/- '',f8.4,
     +              '' [deg] (no resolution!)'')') elatmg,sdlatg
         endif
      endif
      if(var2(3).lt.90.d0) then
         if(output) then
            write(11,'(''Epicenter lon:'',14x,f10.4,'' +/- '',f8.4,
     +              '' [deg]'')')  elonm,var2(3)
         endif
      else
         ibad = ibad + 1
         if(var2(3).gt.180.d0) var2(3) = 180.d0
         if(output) then
            write(11,'(''Epicenter lon:'',14x,f10.4,'' +/- '',f8.4,
     +              '' [deg] (no resolution!)'')') elonm,var2(3)
         endif
      endif
      if(czo.eq.'D') then

         if((var2(4).ge.zo       .and. zo.ge. 50.d0) .or. 
     +      (var2(4).ge..75d0*zo .and. zo.gt.300.d0)      ) then

            ibad = ibad + 1

            if(var2(4).gt.660.d0) var2(4) = 660.d0
            if(output) then
               write(11,'(''Source depth :'',15x,f7.2,''   +/- '',f6.2,
     +               ''   [km] (no resolution!)''/)') zo,var2(4)
            endif

         else

            if(output) then
               write(11,'(''Source depth :'',15x,f7.2,''   +/- '',f6.2,
     +               ''   [km]''/)') zo,var2(4)
            endif

         endif

      else if(czo.eq.'F' .or. czo.eq.'B') then
         if(output) then
            write(11,'(''Source depth :'',15x,f7.2,''   [km] Fixed''
     +                  /)') zo
         endif
      endif

c
c     let us now calculate the final residuals and print them out
c
405   stmean    = 0.d0
      strmean   = 0.d0
      samean    = 0.d0
      sarmean   = 0.d0
      rmsazi    = 0.d0
      spmean    = 0.d0
      sprmean   = 0.d0
      rmsp      = 0.d0
      rms       = 0.d0
      rmsisc    = 0.d0
      wisc      = 0.d0
      dtmin     = 9999.d0

c
c     misfit parameters for all input data!
c
      tmisf     = 0.d0
      tmisfl    = 0.d0
      ntmisf    = 0
      dmisf     = 0.d0
      dmisfl    = 0.d0
      ndmisf    = 0
      amisf     = 0.d0
      amisfl    = 0.d0
      namisf    = 0
      pmisf     = 0.d0
      pmisfl    = 0.d0
      npmisf    = 0
      wmisf     = 0.d0
      wmisfl    = 0.d0
      nwmisf    = 0

      nobst     = 0
      nobsa     = 0
      nobsp     = 0
      stato     = ' '

      if(magflag) then
         namp     = 0
         imsm     = 0
         dmsm     = 0.d0
         imbm     = 0
         dmbm     = 0.d0
         imlm     = 0
         dmlm     = 0.d0
      endif 

      do 450 i = 1,nobs

      epiaz(i) = -999.0
      epiaz2(i) = -999.0
      emeran(i) = -999.d0

      if(uppcas(phase(i)(1:1)).ne.'P' .and. 
     +   uppcas(phase(i)(1:1)).ne.'S' .and. 
     +   phase(i)(1:1).ne.'R' .and. phase(i)(1:1).ne.'L' ) go to 413
c
c     Mark all phases, which have been used as part of a defining
c     travel-time difference measure.
c

      do 412 j = i+1,nobs

         if(uppcas(phase(j)(1:1)).ne.'P' .and. 
     +      uppcas(phase(j)(1:1)).ne.'S' .and. 
     +      phase(j)(1:1).ne.'R' .and. phase(j)(1:1).ne.'L' ) go to 412

         if(sta(iev(i)).ne.sta(iev(j))) go to 412
         if(phaseu(i).eq.phaseu(j)) go to 412
         do 411 i3 = 1, ndt
             if(idtu(i3).eq.(j*i + j+i)) then
                used(i)(4:4) = 'D'
                used(j)(4:4) = 'D'
                go to 412
             endif
411      continue
412   continue

413   if (sta(iev(i)).ne.stato) then

         istad(iev(i)) = 0

         stato = sta(iev(i))

         call depi(stala(iev(i)),stalo(iev(i)),elatmg,elonm,
     +              del(iev(i)),delk(iev(i)),azie(iev(i)),baz(iev(i)),
     +        d2km)
         rzo   = sngl(zo)
         rdel  = sngl(del(iev(i)))
         rdelk = sngl(delk(iev(i)))

CCC      rays(iev(i),1) = stala(iev(i))
CCC      rays(iev(i),2) = stalo(iev(i))
CCC      rays(iev(i),3) = elatmg
CCC      rays(iev(i),4) = elonm
CCC      rays(iev(i),5) = delk(iev(i))
CCC      rays(iev(i),6) = del(iev(i))

         if(rdel.lt.rdmi) rdmi = rdel
         if(rdel.gt.rdma) rdma = rdel

         fla1 = sngl(deg2rad*(90.d0-elatm))
         razi = sngl(azie(iev(i)))

         nphas = 0

         if (mod2flag .and. touse(i)(7:7).eq.'2') then
           modn = modnam2
         else if (mod3flag .and. touse(i)(7:7).eq.'3') then
           modn = modnam3
         else if (mod4flag .and. touse(i)(7:7).eq.'4') then
           modn = modnam4
         else
           if(iloc) imod2 = 1
           modn = modnam
         endif

         loctt = 0

         if(imod2.eq.0 .or. del(iev(i)).gt.rmax .or. zo.gt.zmax) then

           nphas = 0
           call tauget_mod(rzo,rdel,nphas,phcd,ttc,dtdd,
     +                         dtdh,dddp,modn)

         else

           ierr = 0
           indph = istaph(iev(i))*10000 + indph0
           elatc = elatmg
           elonc = elonm

           elat2 = stala(iev(i))
           elon2 = stalo(iev(i))
           sdep = 0.d0
           if(locsta) sdep  = - stael(iev(i))

           nphas = 0
           call ttloc(rzo,rdel,czo,nphas,ttc,dtdd,dtdh,dpdh,dddp,
     +                phcd,rmax,sdep,typctl,ierr,indph,emerout)
           loctt = 1

         endif

      endif

      text(i) = ' '
      arr(i) = rdel

      if(kmout) then
         if (delk(iev(i)) .gt. dismaxst+10.d0 ) go to 450
         if (delk(iev(i)) .lt. disminst-10.d0 ) go to 450
      else
         if (del(iev(i)) .gt. dismaxst+0.1d0 ) go to 450
         if (del(iev(i)) .lt. disminst-0.1d0 ) go to 450
      endif

      epiaz2(iev(i)) = razi

      istad(iev(i)) = 1

      phid = phase(i)
      ttobs  = timemin + tt(i)

      dpa   = 0.d0
      phid1 = phase(i)
      ttres  = -9999.d0
      ttres1 =  9999.d0
      pares  =  -1000.d0
      surf = .false.

      if(sglgdis.gt.0.d0) then
         if(phid.eq.'Sg' .and. delk(iev(i)).ge.sglgdis) phid = 'Lg'
         if(phid.eq.'Lg' .and. delk(iev(i)).lt.sglgdis) phid = 'Sg'
      endif

      call fetoh(ttobs,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

      if(phid.eq.'Rg') then
         if(delk(iev(i)).le.400.d0) then
            surf = .true.
            vsurf = vrg
         else
            phid = 'LR'
         endif
      endif
                
      if(phid.eq.'Lg' ) then
         if(delk(iev(i)).le.3000.d0) then
            surf = .true.
            vsurf = vlg
         else
            phid= 'LR'
         endif
      endif

c
c     LR is also assumed for far-distant Lg and Rg
c
      if(phid.eq.'LR') then
         surf = .true.
         vsurf = vlr
      endif
                
      if(phid.eq.'LQ') then
         surf = .true.
         vsurf = vlq
      endif

      if(phid.eq.'T ') then
         surf = .true.
         vsurf = vt
      endif

      if(phid.eq.'IS') then
         surf = .true.
         vsurf = vi
      endif

      if (surf) then
         nphas       = nphas + 1
         ttc(nphas)  = delk(iev(i))/vsurf
         dtdd(nphas) = d2km/vsurf
         dtdh(nphas) = 0.d0
         phcd(nphas) = phid
         dddp(nphas) = 0.d0
         dpdh(nphas) = 0.d0
      endif
                
      icha = 0
      imin = 0

      if((phid(2:3).eq.'n ' .or. phid(2:3).eq.'g ' .or. phid(2:3).eq.
     +   'b ') .and. rdel.gt.30. ) phid(2:3)='  '

410   continue

      if(phase(i).eq.'P1') then
         nphass = 1
         if (icha.ne.0) then
            if (phid(1:1).eq.'P')   nphass = 2
            if (phid(1:2).eq.'PK' .or. phid(1:3).eq.'Pdi') nphass = 10
         endif
      else
         nphass = nphas
      endif

      dpa2= 0.d0

      do 420 j = 1,nphass

      dpa = 0.d0
c
c     Let's take the first phase from IASPEI-1991-type tables,
c     which fits the phase name of the onset
c
      phid1 = phcd(j)

      firston = .false.

      if (firstph) then

         if (j.eq.1) then
            firston = .true.
            firsts  = .false.
         else
            if (phid1(1:1).eq.'S' .and. .not. firsts) then
                firston = .true.
                firsts  = .true.
            endif
         endif

      endif

      if(phid.eq.'P1' .and. phid1(1:1).eq.'P' .and.
     +   phid1(2:2).ne.'P' .and. phid1(2:2).ne.'S' .and.
     +   (phid1(3:3).eq.' ' .or. phid1(4:4).eq.'f' .or.
     +    phid1(1:5).eq.'PKPdf')) then

        phid = phid1

        if(p(i).lt.0.9*pdif .and. phid1(1:3).eq.'Pdi') then
           if(del(iev(i)).gt.110.d0) then
             phid = 'PKPdf'
             go to 420
           endif
        endif

      endif

      if(phid.eq.'S1' .and. phid1(1:1).eq.'S' .and.
     +   phid1(2:2).ne.'S' .and. phid1(2:2).ne.'P' .and.
     +   (phid1(3:3).eq.' ' .or. phid1(1:5).eq.'SKSac')) then

        phid = phid1

      endif

      if(phid.eq.'PKPdf') then
         
         if(phid1(1:5).eq.'PKiKP' .and.
     +                       del(iev(i)).ge.90.d0) phid='PKiKP'

      endif

      if(icha.eq.0) then
         
        if(phid1(2:4).eq.'dif') then
          if(phid.eq.'P') phid='Pdif'
          if(phid.eq.'S') phid='Sdif'
        endif

        if(phid(1:5).eq.'PKPab' .and. del(iev(i)).lt.150.d0) 
     +     phid(1:5) = 'PKP  '

      endif

      if(phid1.eq.phid .or. imin.eq.1) then

c
c     checking : any ellipticity correction for this phase?
c

         phid1 = phcd(j)
         dpa   = dtdd(j)
         dtdz  = dtdh(j)
         pa    = sngl(dpa)
         dpaa  = dabs(dpa)

         if(single .and. rayokf) then

            if(dabs(p(i)-dpaa) .ge. 0.01d0) go to 420

         endif

         rzoe = rzo

         ecor = 0.

         if (.not.surf .and. .not.locgeo) then
c
c           Ellipticity corrections are yet not available for sources 
c           deeper than 700 km. Therefore we accept a small error
c           and set the source depth to 700. km
c
            if(rzoe.gt.700.) rzoe=700.

            call ellip(fla1,razi,rdel,rzoe,phid1,pa,ecor,ierr)

         endif

         phase_t = phase_type(phid1)

         statict = 0.d0
         if(firstph) then
            if (firston) then
               if(phase_t.eq.'P') statict = statp(iev(i))
               if(phase_t.eq.'S') statict = stats(iev(i))
            endif
         else
            if(phase_t.eq.'P') statict = statp(iev(i))
            if(phase_t.eq.'S') statict = stats(iev(i))
            if(surf.and. phid1(2:2).eq.'g') statict = statr(iev(i))
         endif

         th = 0.d0

         tcrust = 0.d0
         zoc    = 9999.d0

         dph = 0.d0

         if(vlflag .and. .not.surf) then

           hsta = stael(iev(i))

           vloc = 99999.d0
           if(phase_t.eq.'P') vloc = stavp(iev(i))
           if(phase_t.eq.'S') vloc = stavs(iev(i))

           if(vloc .lt. 999.d0) then

              if(imo.eq.2 .or. imo.eq.4 .and. loctt.ne.1) then
              
                 elatc  = stalae(iev(i))
                 elonc  = stalo(iev(i))
                 indr = 1
                 tcrust = crustc(phase_t,dpaa,zoc,indr,typctl)

                 if(dabs(tcrust).gt.0.d0) hsta = hsta - elev

              endif

              if(locsta .and. hsta.lt.0.d0 .and. loctt.eq.1) hsta = 0.d0

C
C   Was ist mit Deltaeffekt der Höhenkorrektur????? Check!!!!
C

              if(dabs(hsta).gt.0.d0) then

                 radkm = deg2rad*radloc(stala(iev(i)),1)
                 phin = vloc*dpaa/radkm

                 if(phin.lt.1.d0) then

                    dl = hsta / dcos(dasin(phin))

                    ddis2 = q2(dl)-q2(hsta)

                    if(ddis2.gt.0.d0) then

                       ddis = dsqrt(ddis2)/radkm

                       if(dl.gt.0.d0) then
                          th = dl/vloc - dpaa*ddis/radkm
                          dph = ddis * dddp(j)
                       else if(dl.lt.0.d0) then
                          th = dl/vloc + dpaa*ddis/radkm
                          dph = ddis * dddp(j)
                       endif

                    endif

                 endif

              endif
           endif
         endif

c
c        We have eventually to correct this phase for the local 
c        structure at the reflection point at the surface (if 
c        CRUST 1.0 or a local/regional model is available).
c
         trefl  = 0.d0
         trefls = 0.d0
         treflp = 0.d0

         if( .not.surf .and. used(i)(5:5).eq.'R' .and. loctt.eq.0) then 

            if(phid1(1:1).ne.'p' .and. phid1(1:1).ne.'s') then
               if(phid1(1:1).eq.phid1(2:2) .and. 
     +            phid1(1:1).eq.phid1(3:3))  then
                  goto  415
               endif
               if(phid1(1:2).eq.phid1(3:4) .and. 
     +            phid1(1:2).eq.phid1(5:6))  then
                  goto  415
               endif
            else
               if(phid1(2:2).eq.phid1(4:4) .and. 
     +            phid1(2:2).eq.phid1(5:5))  then
                  goto  415
               endif
               if(phid1(2:4).eq.phid1(5:6) .and. 
     +            phid1(2:4).eq.phid1(7:8))  then
                  goto  415
               endif
            endif

            chgcas = uppcas(phid1)
            phase_t = chgcas(1:1)
             
            fmult = 1.d0
            del0  = dirdel(dpaa,zo,fmult,phase_t)
            azi0 = azie(iev(i))

            if(phid1(1:1).ne.'p' .and. phid1(1:1).ne.'s') then
               if(dpa.ge.0.d0) then
                  del0 = (del(iev(i))+del0)/2.d0
               else
                  del0 = (360.d0 - del(iev(i)) - del0)/2.d0
                  azi0 = alpha2(azie(iev(i)) + 180.d0)
               endif
            endif
                
            inddel = 1
            call delazd(elatmg,elonm,azi0,del0,inddel,elatc,elonc)
    
c
c     correction for depth phases (e.g.: pP, sP...)
c
            if(phid1(1:1).eq.'p' .or. phid1(1:1).eq.'s' .and.
     +         phase_t.ne.' ' .and. del0.le.rmax)  then

               if (zo.ge.zmax) then

                   if ((phid1(1:1).eq.'p' .and. phid1(2:2).eq.'P') .or.
     +                 (phid1(1:1).eq.'s' .and. phid1(2:2).eq.'S')) then
                       indr = 2
                   else
                       indr = 3
                   endif
                   zoc    = 9999.d0

               else if(zo.lt.zmax) then

                   if ((phid1(1:1).eq.'p' .and. phid1(2:2).eq.'P') .or.
     +                 (phid1(1:1).eq.'s' .and. phid1(2:2).eq.'S')) then
                       indr = 4
                   else
                       indr = 5
                   endif
                   zoc    =  zo

               endif

               trefls = crustc(phase_t,dpaa,zoc,indr,typctl)

            endif

c
c     correction for surface multiples (e.g.: PnPn,...,PP,SS,P'P')
c
            if( phid1(1:1).eq.phid1(2:2) .or.
     +          phid1(1:2).eq.phid1(3:4)  ) then

                indr = 2
                zoc    = 9999.d0
                treflp = crustc(phase_t,dpaa,zoc,indr,typctl)
            endif

c
c      correction for converted surface multiples (e.g.: PnSn,...)
c
            conr = .false.
            if( (phid1(1:1).eq.'P' .or. phid1(2:2).eq.'P') .and.
     +          (phid1(1:1).eq.'S' .or. phid1(2:2).eq.'S') .and.
     +          (phid1(3:3).eq.' ' .or. phid1(3:3).eq.'g' .or.
     +           phid1(3:3).eq.'b' .or. phid1(3:3).eq.'n')) then
                   conr=.true.
                   phidr = phid1(2:)
            endif

            if( (phid1(1:1).eq.'P' .or. phid1(3:3).eq.'P') .and.
     +          (phid1(1:1).eq.'S' .or. phid1(3:3).eq.'S') .and.
     +           phid1(2:2).ne.'b' .and. phid1(2:2).ne.'m' .and.
     +           phid1(2:2).ne.'c' .and. phid1(2:2).ne.'k' .and.
     +          (phid1(2:2).eq.phid(4:4) .or. phid1(2:2).eq.'g' .or.
     +           phid1(2:2).eq.'n'                      )) then
                   conr=.true.
                   phidr = phid1(3:)
            endif

            if(conr) then

                zor = 0.d0
                call tauget_ray(phidr,phase_t,dpaa,modn,zor,
     +                        del0,ttray,rayok)

                azi0 = azie(iev(i))
                if(dpa.lt.0.d0) azi0 = alpha2(azie(iev(i))-180.d0)

                inddel = 1
                call delazd(stala(iev(i)),stalo(iev(i)),azi0,del0,
     +                      inddel,elatc,elonc)

                indr = 3
                zoc    = 9999.d0
                treflp = crustc(phase_t,dpaa,zoc,indr,typctl)

            endif

            trefl = trefls + treflp

            if (dabs(trefl).gt.0.d0.and.typctl.gt.6) then
               print *,'dirdel: ',phid,' azi ',azi0,' del ',del0
               print *,'  lat ',elatc,' lon ',elonc,' trefl ',trefl
            endif

         endif

415      continue

         ttt(i) = tome + ttc(j) + dble(ecor) 
     +                 + th + tcrust + trefl + statict

         if (typctl.gt.8) then
           print *,'ttt, t0, TT, ECOR, Height, Crust, Refl, Stat,Tobs'
           print
     +     *,ttt(i),tome,ttc(j),dble(ecor),th,tcrust,trefl,statict,
     +       ttobs
         endif

         ttres  = ttobs - ttt(i)

         if(ttres.gt.100.d0 .and. phase(i)(1:2).eq.'P1' .and. 
     +      phid(1:3).eq.'Pdi' .and. used(i)(1:1).eq.'T') then
            phid='PKPdf'
            go to 420
         endif

         ttres1 = dabs(dble(ttc(j) - ttc(1)))

         if(rdel.gt.110. .and. phcd(1)(1:3).eq.'Pdi' .and.
     +                         ttres1.gt.200.d0        ) then

           if(phcd(2)(1:2).eq.'PK') ttres1 = dabs(dble(ttc(j) - ttc(2)))
           if(phcd(3)(1:2).eq.'PK') ttres1 = dabs(dble(ttc(j) - ttc(3)))
           if(phcd(4)(1:2).eq.'PK') ttres1 = dabs(dble(ttc(j) - ttc(4)))

         endif

         if(imin.eq.0) then

            if(used(i)(1:1).eq.'T') then
               stmean  = stmean + ttres
               strmean = strmean + dabs(ttres)
               rms     = rms    + q2(ttres)
               rmsisc  = rmsisc + q2(ttres)*ttu(i)
               wisc    = wisc   + ttu(i)
               epiaz(i)= razi
               nobst   = nobst + 1
            endif

            pares = p(i)-dabs(dpa - dph)
            if(used(i)(3:3).eq.'S') then
               spmean  = spmean + pares
               sprmean = sprmean + dabs(pares)
               rmsp    = rmsp + q2(pares)
               epiaz(i)= razi
               nobsp   = nobsp + 1
            endif

            if(emerout .and. .not.surf) then
c
c           we use global models only
c
               if((rzv(1).lt.0.d0 .and. uppcas(phid1(1:1)).eq.'P')  .or.
     +            (rzv(2).lt.0.d0 .and. uppcas(phid1(1:1)).eq.'S')) then
                  dtdx = dabs(dpa - dph)*180.d0/pi/rearth
                  if(dtdz.lt.0.d0) then
                     emeran(i) = rad2deg*datan(dabs(dtdx/dtdz))
                  else if(dtdz.gt.0.d0) then
                     emeran(i) = 180.d0 - rad2deg*datan(dabs(dtdx/dtdz))
                  else
                     emeran(i) = 90.d0
                  endif
               else
c
c           a local/regional model had been used to locate the source 
c

                  vsource = -999.d0
                  if(uppcas(phid1(1:1)).eq.'P') vsource = rzv(1)
                  if(uppcas(phid1(1:1)).eq.'S') vsource = rzv(2)
                  if(vsource.gt.0.d0) then
                     fac = vsource*180.d0*dabs(dpa-dph)/(pi*(rearth-zo))
                     if(dabs(fac).le.1.d0) then
                        emeran(i) = rad2deg*dasin(fac)
                        if(dtdz.gt.0.d0) emeran(i) = 180.d0 - emeran(i)
                     endif
                  endif
               endif

c              print*,'1 ',phcd(j),dtdd(j),dtdx,vsource,dtdz,emeran(i)

            endif

c          if(typctl.ge.8) then
c            print *,'i,ttt,ttobs,ttres,ecor,th,tcrust,used,tts,tt2,ttu'
c            print *,i,ttt(i),ttobs,ttres,ecor,th,tcrust,used(i),tts(i),
c    +               tt2(i),ttu(i),emeran(i)
c          endif

            go to 430

         else if(imin.eq.1.and.dabs(ttres).lt.dabs(dtmin)) then

            phid2 = phid1
            dtmin = ttres
            ttresm = ttres1
            pares2 = p(i)-dabs(dpa-dph)
            dpa2  = dpa
            if(emerout .and. .not.surf) then
c
c           we use global models only
c

               if((rzv(1).lt.0.d0 .and. uppcas(phid1(1:1)).eq.'P')  .or.
     +            (rzv(2).lt.0.d0 .and. uppcas(phid1(1:1)).eq.'S')) then
                  dtdx = dabs(dpa - dph)*180.d0/pi/rearth
                  if(dtdz.lt.0.d0) then
                     emeran(i) = rad2deg*datan(dabs(dtdx/dtdz))
                  else if(dtdz.gt.0.d0) then
                     emeran(i) = 180.d0 - rad2deg*datan(dabs(dtdx/dtdz))
                  else
                     emeran(i) = 90.d0
                  endif
               else
c
c           a local/regional model had been used to locate the source 
c
                  vsource = -999.d0
                  if(uppcas(phid1(1:1)).eq.'P') vsource = rzv(1)
                  if(uppcas(phid1(1:1)).eq.'S') vsource = rzv(2)
                  if(vsource.gt.0.d0) then
                     fac = vsource*180.d0*dabs(dpa-dph)/(pi*(rearth-zo))
                     if(dabs(fac).le.1.d0) then
                        emeran(i) = rad2deg*dasin(fac)
                        if(dtdz.gt.0.d0) emeran(i) = 180.d0 - emeran(i)
                     endif
                  endif
               endif

c              print*,'2 ',phcd(j),dtdd(j),dtdx,vsource,dtdz,emeran(i)

            endif

         endif

      endif

420   continue

      if (imin.eq.1) go to 425

c
c     Try it with another phase-name from the same phase-type.
c
 
      if(single) go to 410

      call testphase (phid,icha)

      if(icha.ne.999) go to 410

      if(imin.eq.0) then
         imin  = 1
         dtmin = 9999.d0
         go to 410
      endif

425   if(dabs(dtmin).le.15.d0 .or. used(i)(1:1).eq.'T' .or.
     +   (used(i)(1:1).eq.'t' .and. used(i)(4:4).eq.'D')) then
         ttres  = dtmin
         ttres1 = ttresm
         pares  = pares2
         phid1  = phid2
         dpa    = dpa2
      else
         ttres  = -9999.d0
         ttres1 =  9999.d0
         pares  = -999.d0
         dpa    = 0.d0
         phid1  = ' '
      endif

430   continue

      phaseu(i) = phid1

      if(used(i)(1:1).eq.'t') used(i)(1:1) = ' '

      if(touse(i)(1:1).eq.'T' .and. dabs(ttres).lt.999.d0) then
         fmis    = ttres/ttu(i)
         tmisfl  = tmisfl + dabs(fmis)
         tmisf   = tmisf + q2(fmis)
         ntmisf  = ntmisf + 1
      endif

      if(touse(i)(3:3).eq.'S' .and. pares.gt.-99.d0 ) then
         fmis    = pares/ps(i)
         pmisfl  = pmisfl + dabs(fmis)
         pmisf   = pmisf + q2(fmis)
         npmisf  = npmisf + 1
      endif

      if(dpa.lt.0.d0 .or.phase(i)(1:4).eq.'P3KP') then
         azires = alpha1(azi(i) - alpha2(baz(iev(i))-180.d0))
      else
         azires = alpha1(azi(i) - baz(iev(i)))
      endif

      if(touse(i)(2:2).eq.'A' ) then
         fmis    = azires/azis(i)
         amisfl  = amisfl + dabs(fmis)
         amisf   = amisf + q2(fmis)
         namisf  = namisf + 1
      endif

      if(used(i)(2:2).eq.'A') then
         samean  = samean + azires
         sarmean = sarmean + dabs(azires)
         rmsazi  = rmsazi + q2(azires)
         epiaz(i)= razi
         nobsa   = nobsa + 1
      endif

      imod = ' '
      if(touse(i)(7:7).eq.'2' .and. mod2flag) imod = '2'
      if(touse(i)(7:7).eq.'3' .and. mod3flag) imod = '3'
      if(touse(i)(7:7).eq.'4' .and. mod4flag) imod = '4'

      do 432 iu = 1,4
      if(used(i)(iu:iu).eq.' ' .and. touse(i)(iu:iu).eq.'m') then
         used(i)(iu:iu) = 'm'
      endif
432   continue

      statw = stato
      if(touse(i)(9:9).eq.'*') then
        chgcas = lowcas(stato)
        statw = chgcas(1:5)
      endif

      write(text(i),'(a5,f8.3,f7.2,1x,a8,8x,2i3.2,f7.3,f8.3,
     +                f7.2,f8.2,f6.2,f7.2,1x,a5,a1)') 
     +                statw,rdel,razi,phase(i),hh,mi,sec,
     +                ttres,azi(i),azires,p(i),pares,used(i)(1:5),imod
      
      if(kmout)             write(text(i)(6:13),'(f8.2)') rdelk

      if(phase(i).ne.phid1) write(text(i)(30:37),'(a8)') phid1

      if(ttres.lt.-99.999d0 .or. ttres.gt.999.999d0) then
         write(text(i)(51:58),'(1x,f7.2)') ttres
      endif

      if(dabs(ttres).gt.998.d0) text(i)(51:58)='        '
      if(dabs(ttres).lt.1.d-3) text(i)(51:58)='   0.000'

      if(azires.lt.-360.d0) text(i)(66:73)='        '
      if(azires.gt.360.d0)  text(i)(66:73)='        '
      if(azi(i).le.0.d0) then
           text(i)(59:73)='               '
           if(thbaz) write (text(i)(60:65),'(F6.2)') baz(iev(i))
      endif

      if(pares.ge. 990.d0)  text(i)(80:86)='       '
      if(pares.le.-990.d0)  text(i)(80:86)='       '
      if(p(i).le.0.d0)   then
           text(i)(74:86)='             '
           if(thray) write (text(i)(75:79),'(F5.2)') dpa
      endif

      if(snr(i).gt.0.d0) then
         if(snr(i).lt.10000.d0) then
            write (text(i)(94:101),'(1x,f7.2)') snr(i)
         else
            if(snr(i).lt.100000.d0) then
               write (text(i)(94:101),'(1x,f7.1)') snr(i)
            else
               if(snr(i).lt.1000000.d0) then
                  write (text(i)(94:101),'(1x,f7.0)') snr(i)
               else
                  text(i)(94:101) = ' 999999.'
               endif
            endif
         endif
         isnr = isnr + 1
      endif

      if(amplit(i).gt.0.d0) then
         write (text(i)(102:114),'(1x,f12.2)') amplit(i)
      endif

      if(period(i).gt.0.d0) then
         write (text(i)(115:121),'(1x,f6.3)') period(i)
      endif

      onscha = '_'
      if (phase_type(phase(i)).eq.'P') then
         if(ttu(i).le.sisfe) onscha = 'e'
         if(ttu(i).le.sisfi) onscha = 'i'
      else if(phid1.ne.'LR'.and.phid1.ne.'LQ') then
         if(ttu(i).le.sisfe*2.d0) onscha = 'e'
         if(ttu(i).le.sisfi*2.d0) onscha = 'i'
      endif

      text(i)(131:131) = onscha

c     print *, statw, phid1, ttres, ttres1, touse(i)
      if(amplit(i).gt.0.d0) then

         namp = namp + 1

         if (magflag .and. touse(i)(6:6).eq.'M') then

c
c     the standard IASPEI (1967) formula:
c     Ms = log (A / T )  + 1.66 log (DELTA) + 0.3 (A in nanometer)
c
c     or the Rezapour/Pearce(BSSA 88, 43-61) formula:
c
c     Ms = log (A / T ) + 1/3 log (DELTA) + 1/2 log (sin(DELA) +
c          0.0046 DELTA + 2.370  (A in nanometer)
cc

           dmag = -9.99d0
           statmag='  '

           if(phid1.eq.'LR' .and. period(i).gt.5.d0 ) then

              d1 = del(iev(i))
              if(magtyps(1:6).eq.'IASPEI') then
                dmag = dlog10(amplit(i)/period(i)) +
     +              dlog10(d1)*1.66d0  + 0.3d0
              else if(magtyps(1:3).eq.'R-P') then
                dmag = dlog10(amplit(i)/period(i)) +
     +              dlog10(d1)/3.d0 + dlog10(d1*deg2rad)/2.d0
     +              + 0.0046d0*d1 + 2.370d0
              else
                print *,' Ms attenuation model not defined!'
                go to 449
              endif

              dmsm = dmsm + dmag
              imsm = imsm + 1
              statmag = 'MS'

           else if(dabs(ttres).lt.60.d0 .and. (phid1.eq.'Lg'    .or. 
     +            (phid1(1:1).eq.'S'.and.(phid1(2:2).eq.'g'.or.
     +             phid1(2:2).eq.'b' .or.phid1(2:2).eq.'n'.or.
     +             phid1(2:2).eq.' '))           .or.
     +           ((phid1(1:2).eq.'pS'.or.phid1(1:2).eq.'sS').and.
     +            (phid1(3:3).eq.'g' .or.phid1(3:3).eq.'b' .or.
     +             phid1(3:3).eq.'n' .or.phid1(3:3).eq.' '))
     +            ))then
c
c             we will use ML attenuation file for S-type onsets
c

              magfile = file_check(magmlfile)
              irc = 0

              if(magfile.ne.' ') then

                 if ((magtypml(1:7).eq.'Richter').and.
     +               (period(i).le.0.d0)) period(i) = 1.d0

                 if (period(i).gt.0.d0) then
                    call epmagc(real(period(i)), rdelk, rmcorr, typctl, 
     +                       irc, magfile)
                 else
                    irc = 9
                 endif

                 if(irc.eq.0) then

                    if (magtypml(1:5) .eq. 'Bath ') then
                       dmag = dlog10(amplit(i)*0.1d0) + dble(rmcorr)
                    else if (magtypml(1:7) .eq. 'Richter') then
                       dmag = dlog10(amplit(i)) + dble(rmcorr)
                    else
                       if(typctl . gt. 6 ) then 
                          print *,' Cannot find ML attenuation model '
     +                           , magtypml
                       endif
                       go to 449
                    endif

                    dmlm = dmlm + dmag
                    imlm = imlm + 1
                    statmag = 'ml'

                 else

                    if(typctl . gt. 6 ) then
                      print *,' No ML attenuation corrections found!'
                    endif
                    go to 449

                 endif
              else
                 if(typctl . gt. 6 ) then 
                    print *,' No ML attenuation model defined!'
                 endif
              endif

           else if( ttres1.lt.9.d0 .and. dabs(ttres).le.8.d0 .and. 
     +             (phid1(1:1).eq.'P' .or. phid1(1:2).eq.'pP' .or. 
     +              phid1(1:2).eq.'sP') .and.period(i).gt.0.d0) then

              magfile = ' '

              if(rdel.le.110. ) then

                if(magtypp.eq.'G-R' .and. rdel.ge.11.) 
     +                        magfile = 'MB_G-R.DAT'

                if(magtypp.eq.'V-C') magfile = 'MB_V-C.DAT'

                if(magtypp.eq.'M-R' .and. (rdel.le.100. 
     +                        .and. rdel.ge.21.)) magfile = 'MB_M-R.DAT'

              else if( rdel.gt.110 .and. rdel.le.150. .and. 
     +                (phid1(1:3).eq.'PKP'  .or. phid1(1:3).eq.'PKi' 
     +            .or. phid1(1:4).eq.'pPKP' .or. phid1(1:4).eq.'sPKP'
     +            .or. phid1(1:4).eq.'pPKi' .or. phid1(1:4).eq.'sPKi')
     +               .and. magtypp.eq.'V-C' ) then

                magfile = 'MB_V-C.DAT'

              else if( rdel.gt.150 .and. 
     +               (phid1(1:5).eq.'PKPdf' .or. phid1(1:6).eq.'pPKPdf'
     +                 .or. phid1(1:6).eq.'pPKPdf') 
     +               .and. magtypp.eq.'V-C') then

                magfile = 'MB_V-C.DAT'

              endif

              if(magfile.ne.' ') then

                 magfile = file_check(magfile)
                 irc = 0
                 call magfact(magfile, rdel, rzo, rmcorr, irc)
                 if(irc.eq.0) then

                    dmag = dlog10(amplit(i)/period(i)) + dble(rmcorr)

                    dmbm = dmbm + dmag
                    imbm = imbm + 1
                    statmag = 'mb'

                 else
                    if(typctl . gt. 6 ) then 
                       print *,' No mb attenuation corrections found!'
                    endif
                 endif

              else
                 if(typctl . gt. 6 ) then 
                    print *,' No mb attenuation corrections file!'
                 endif
              endif

           endif

449        if(dmag.gt.-9.99d0) then
              if (dmag.lt.0.d0) then
                 write (text(i)(122:129),'(1x,f4.1,1x,a2)') dmag,statmag
              else
                 write (text(i)(122:129),'(1x,f4.2,1x,a2)') dmag,statmag
              endif
           endif
        endif
      endif

      if(emerout .and. emeran(i).ge.0.d0) then
        write(text(i)(132:138),'(1x,f6.2)') emeran(i)
      endif

      if(arid(i).ne.'        ') then
        write(text(i)(139:147),'(1x,a8)') arid(i)
      endif

      if(typctl.ge.10) then
         print*,i,text(i)
      endif

450   continue

      if(nobst.gt.0) then
         stmean  = stmean/dble(nobst)
         strmean = strmean/dble(nobst)
         rms     = dsqrt(rms/dble(nobst))
         rmsisc  = dsqrt(rmsisc/wisc)
      endif 

      if(dabs(stmean).gt.0.025d0 .and. mttres) then

         if(output) then
            write(11,'(''Source time corrected for mean '', 
     +           ''travel-time residual ('',f8.3,'')'')') stmean
         endif

         if(dabs(stmean).gt.var(1) .and. itso.lt.1) then
            
            miteras  = miteras + iter
            iter     = 0
            nextiter = 0
            dchang   = dchang0
            iteraz   = 0
            itso     = itso + 1

            tom  = tom  + stmean
            tome = tom  + timemin

            rs(1) = dpythag(var(1),stmean)
            rs(2) = var(2)
            rs(3) = var(3)
            rs(4) = var(4)

            go to 100

         endif

         tome = tome + stmean
         
         call fetoh(tome,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

         if(output) then

            write(11,'(/''Source time  :'',i5,4i3.2,f7.3,'' +/- '',
     +              f8.3,'' [s]'')')  yy,mon,dd,hh,mi,sec,var2(1)
            write(11,'(''        or'',12x,f16.3,'' +/- '',f8.3,
     +              '' [s]'')') tome,var2(1)

            isec1 = nint(sec*1000)
            isec  = isec1/1000
            msec  = isec1-isec*1000
            write(11,'(''        or'',7x,i4,''-'',i3.3,'':'',
     +            3(i2.2,''.''),i3.3,'' +/- '',f8.3,'' [s]''//)') 
     +            yy,idoy,hh,mi,isec,msec,var2(1)

         endif

         go to 405

      endif

      if(nq.lt.3 .and.single) iellip =.false.

      if (iellip .and. .not.last) then

         call ellcal(elatmg,ax1,ax2,fchi2,elmax,elmin,eazi,earea,fmud)

         if(output) then

            if(earea.lt.10000.d0) then
              write(11,'( ''Epicenter error ellipse:''/
     +           ''Major half axis: '',f8.2,
     +           '' [km]  Minor half axis: '',f8.2,
     +           '' [km]'',/''Azimuth:'',f11.1,''  [deg] '',
     +           ''Area: '',f14.2,'' [km**2]''/)')
     +           elmax,elmin,eazi,earea
            else
              write(11,'( ''Epicenter error ellipse:''/
     +           ''Major half axis: '',f8.0,
     +           '' [km]  Minor half axis: '',f8.0,
     +           '' [km]'',/''Azimuth:'',f11.1,''  [deg] '',
     +           ''Area: '',2p e14.6,'' [km**2]''/)')
     +           elmax,elmin,eazi,earea
            endif

         endif

c        radsource = radloc(elatmg,1)
c        print *,'Earth radius in source region: ',radsource

      else if (iellipi.eq.1 .and. .not.iellip) then

         if(output) then

            write(11,'( ''Epicenter error ellipse calculation '',
     +         ''not possible (too less parameter resolution)''/)')

         endif

      endif

      if (ibad.ge.3 .and. .not.last) go to 455

      rlat = sngl(elatmg)
      rlon = sngl(elonm)
      call hyposat_geo( rlat,rlon, isreg, regnum, region , ierr )

      if(output) then
         write(11,'(''Flinn-Engdahl Region ('',i4,'' ): '',a/)')
     +            isreg, trim(region)
      endif

      if(ibad.ge.3) then
         ibad0 = ibad0 + 1
      else if(ibad.eq.0) then
         ibad0 = 0
      endif
c
c     We will now calculate the maximum azimuthal gaps for 
c        - all as defining used observations
c        - all observations
c
c     If possible, also the secondary azimuthal gap is calculated
c

      dazgap = 360.
      dazgap2 = 360.
      if(nstat.gt.1) then

         if (.not.gapobs) then

           call indexx(nobs,epiaz,indx)
           call azigap(epiaz,dazgap,d1azi,d2azi,dazgap2,d1azi2,
     +                 d2azi2,nobs,indx)

         else

           call indexx(nobs,epiaz2,indx)
           call azigap(epiaz2,dazgap,d1azi,d2azi,dazgap2,d1azi2,
     +                 d2azi2,nobs,indx)

         endif

      endif

c     now we count the number of stations available in the chosen 
c     distance range

      nstata = 0
      do 451 i = 1, nstat
        nstata = nstata + istad(i)
451   continue

c
c     if chosen: open the 'new' input file
c
      if(new_in) then

         inputfilen = trim(inputfile) // '_rev'
         open (unit=31,file=trim(inputfilen))
         write (31,'(a)') trim(title)

      endif

c
c     if chosen: open the ISF formatted output file
c
      if(isf_out) then

        if(trim(outputfile).eq.'hyposat-out') then
           outputisf = 'hyposat-isf'
        else
           outputisf = trim(outputfile) // '_isf'
        endif
        open (unit=12,file=outputisf)

        cdum = ' '
        cdum2 = 'BULLETIN'
        dformat = 'IMS1.0'
        dsformat = 'short'
        itest = write_data_type(12,cdum2,cdum,dformat,dsformat)

        cdum2 = ' '

        if(corid.eq.' ' .or. corid.eq.'_') then

          if(cevid(9:9).eq.' ') then
            corid = cevid(1:8)
          else if(cevid(1:1).eq.' ') then
            corid = cevid(2:9)
          else
            corid = '99999999'
          endif

        endif

        cdum = ' '
        rdum = REAL(ISF_NULL)

        itest = write_event_id(12,cevid,region)

        if(itest.eq.20) then
           itest = write_isf_error(typctl)
        endif

        itest = write_origin_head(12)

        call fetoh(tome,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

        isec1 = nint(sec*1000)
        isec  = isec1/1000
        isec2 = isec1-isec*1000

        relmax = real(elmax)
        if(relmax.lt.0.) relmax = 0.
        if(relmax.gt.9999.) relmax = 9999.

        relmin = real(elmin)
        if(relmin.lt.0.) relmin = 0.
        if(relmin.gt.9999.) relmin = 9999.

        ieazi = nint(eazi)
        if(ieazi.lt.0) ieazi = 0

        rmsisf = rms
        if(lrmsisc) rmsisf = rmsisc

        cdum = ' '

        if(czo.eq.'D') then

           itest = write_origin(12,yy,mon,dd,hh,mi,isec,isec2,cdum,
     +     real(var2(1)),real(rmsisf),real(elatmg),real(elonm),cdum,
     +     relmax,relmin,ieazi,real(zo),cdum,
     +     real(var2(4)),in,nstata,nint(dazgap),rdmi,rdma,cdum,cdum,
     +     cdum,author,corid)

        else if(czo.eq.'F' .or. czo.eq.'B') then

           itest = write_origin(12,yy,mon,dd,hh,mi,isec,isec2,cdum,
     +     real(var2(1)),real(rmsisf),real(elatmg),real(elonm),cdum,
     +     relmax,relmin,ieazi,real(zo),'f',rdum,
     +     in,nstata,nint(dazgap),rdmi,rdma,cdum,cdum,
     +     cdum,author,corid)

        endif

        if(itest.eq.20) then
           itest = write_isf_error(typctl)
        endif

        write(outputisf,'(''Confidence level of given uncertainties:'',
     +         f7.2,'' %'')') confl
        itest = write_comment(12,trim(outputisf))
        if(itest.eq.20) then
           itest = write_isf_error(typctl)
        endif

        if(diffflag .and. ndt.gt.0) then
          write(outputisf,'(i4,'' Travel-time differences used as '',
     +         ''defining'')') ndt
          itest = write_comment(12,outputisf)
          if(itest.eq.20) then
             itest = write_isf_error(typctl)
          endif

        endif

        if(iloc) then
           if( (mod2flag .or. mod3flag .or. mod4flag) .and.
     +         (imodn(2).gt.1 .or. imodn(3).gt.1 .or. imodn(4).gt.1)
     +          ) then

               write(outputisf,'(''First reference models  : '',a,
     +               '' and '',a)') trim(filloc),modnam
               itest = write_comment(12,outputisf)

               if(mod2flag .and. imodn(2).gt.1) then
                  write(outputisf,'(''Second reference model  : '',a)')
     +                 modnam2
                  itest = write_comment(12,outputisf)
               endif

               if(mod3flag .and. imodn(3).gt.1) then
                  write(outputisf,'(''Third reference model   : '',a)')
     +                 modnam3
                  itest = write_comment(12,outputisf)
               endif

               if(mod4flag .and. imodn(4).gt.1) then
                  write(outputisf,'(''Fourth reference model  : '',a)')
     +                 modnam4
                  itest = write_comment(12,outputisf)
               endif
           else
               write(outputisf,'(''Reference models  : '',a,'' and '',
     +               a)') trim(filloc),modnam
               itest = write_comment(12,outputisf)
           endif
        else
           if( (mod2flag .or.  mod3flag .or. mod4flag) .and.
     +         (imodn(2).gt.1 .or. imodn(3).gt.1 .or. imodn(4).gt.1)
     +          ) then

               write(outputisf,'(''First reference model   : '',a)') 
     +               modnam
               itest = write_comment(12,outputisf)

               if(mod2flag .and. imodn(2).gt.1) then
                  write(outputisf,'(''Second reference model  : '',a)')
     +                 modnam2
                  itest = write_comment(12,outputisf)
               endif

               if(mod3flag .and. imodn(3).gt.1) then
                  write(outputisf,'(''Third reference model   : '',a)')
     +                 modnam3
                  itest = write_comment(12,outputisf)
               endif

               if(mod4flag .and. imodn(4).gt.1) then
                  write(outputisf,'(''Fourth reference model  : '',a)')
     +                 modnam4
                  itest = write_comment(12,outputisf)
               endif
           else
               write(outputisf,'(''Reference model   : '',a)') modnam
               itest = write_comment(12,outputisf)
           endif
        endif

      endif

      if(itest.eq.20) then
         itest = write_isf_error(typctl)
      endif


      call indexx(nobs,arr,indx)

      ntext = 0
      stato = sta(iev(indx(1)))
      lenons = 131

      texth = ' '
      texth(1:45)   = ' Stat  Delta   Azi   Phase   [used]    Onset '
      texth(46:91)  = 'time    Res     Baz     Res   Rayp   Res  Used' 
      texth(92:131) = '                                       Q' 

      if(isnr.gt.0) texth(97:99) = 'SNR'

      if(emerout) then
          texth(133:138) = 'Em-Ang'
          lenons = 138
      endif

      if(larid) then
          texth(139:143) = ' ARID'
          lenons = 147
      endif

      if(namp.gt.0 ) then

        texth(104:120) = 'Amplitude  Period'

        if(imsm.gt.0) then
           dmsm = dmsm / dble(imsm)
        endif

        if(imbm.gt.0) then
           dmbm = dmbm / dble(imbm)
        endif

        if(imlm.gt.0) then
           dmlm = dmlm / dble(imlm)
        endif

        if(output) then

           if(imsm.gt.0 .and.imbm.eq.0 .and. imlm.eq.0) then
              write(11,'(''Magnitude: '',f4.1,'' (Ms, '',a,'')''/)') 
     +              dmsm,trim(magtyps)
           else if (imsm.eq.0 .and. imlm.eq.0 .and. imbm.gt.0) then
              write(11,'(''Magnitude: '',f4.1,'' (mb, '',a,'')''/)')
     +              dmbm,trim(magtypp)
           else if (imsm.eq.0 .and. imlm.gt.0 .and. imbm.eq.0) then
              write(11,'(''Magnitude: '',f4.1,'' (ML, '',a,'')''/)')
     +              dmlm,trim(magtypml)
           else if (imsm.gt.0 .and.imbm.gt.0 .and. imlm.eq.0) then
              write(11,'(''Magnitudes: '',f4.1,'' (mb, '',a,'') '',
     +              f4.1,'' (Ms, '',a,'')''/)') dmbm,trim(magtypp)
     +              ,dmsm,trim(magtyps)
           else if (imsm.gt.0 .and.imbm.eq.0 .and. imlm.gt.0) then
              write(11,'(''Magnitudes: '',f4.1,'' (ML, '',a,'') '',
     +              f4.1,'' (Ms, '',a,'')''/)') dmlm,trim(magtypml),
     +              dmsm,trim(magtyps)
           else if (imsm.eq.0 .and.imbm.gt.0 .and. imlm.gt.0) then
              write(11,'(''Magnitudes: '',f4.1,'' (mb, '',a,'') '',
     +              f4.1,'' (ML, '',a,'')''/)') dmbm,trim(magtypp),
     +              dmlm,trim(magtypml)
           else if (imsm.gt.0 .and.imbm.gt.0 .and. imlm.gt.0) then
              write(11,'(''Magnitudes: '',f4.1,'' (mb, '',a,'') '',
     +              f4.1,'' (Ms, '',a,'') '',f4.1,'' (ML, '',a,'')''/)')
     +              dmbm,trim(magtypp),dmsm,trim(magtyps),
     +              dmlm,trim(magtypml)
           endif
         endif

         if(imsm.ne.0 .or. imbm.ne.0 .or. imlm.ne.0) 
     +     texth(124:126) = 'MAG'

      endif

      write (11,'(a/)') trim(texth)

      if(isf_out .and. (imsm.gt.0 .or. imbm.gt.0 .or. imlm.gt.0)) then

        itest = write_netmag_head(12)

        if(imbm.gt.0) then
           itest = write_netmag(12,'mb', cdum, real(dmbm), rdum, imbm,
     +             trim(magtypp),corid)
        endif

        if(imsm.gt.0) then
           itest = write_netmag(12,'MS', cdum, real(dmsm), rdum, imsm,
     +             trim(magtyps),corid)
        endif

        if(imlm.gt.0) then
           itest = write_netmag(12,'Ml', cdum, real(dmlm), rdum, imlm,
     +             trim(magtyps),corid)
        endif

        if(itest.eq.20) then
           itest = write_isf_error(typctl)
        endif

      endif

      if(isf_out) then
        
        itest = write_phase_head(12)

      endif

      do 453 i=1,nobs

      if(new_in) then

c        print *,i,touse(i),' , ',used(i),text(i)
         if(touse(i)(1:1).eq.'T' .or. touse(i)(2:2).eq.'A' .or.
     +      touse(i)(3:3).eq.'S' .or. touse(i)(4:4).eq.'D'   ) then

            string = ' '

            string(1:5) = text(i)(1:5)
            string(7:14) = text(i)(22:29)
            if(text(i)(30:37).ne.'        ') 
     +         string(7:14) = text(i)(30:37)

            string(71:77) = touse(i)(1:7)

            if(touse(i)(1:1).eq.'m') then
               chgcas = lowcas(string(7:14))
               string(7:14) = chgcas(1:8)
               string(71:71) = 'x'
               string(74:74) = 'x'
            endif

            ttobs  = timemin + tt(i)
            mm = ' '
            idoy = 0
            call fetoh(ttobs,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

            write(string(16:44),'(i4,4(1x,i2),1x,f6.3,1x,f5.3)') 
     +            yy,mon,dd,hh,mi,sec,ttu(i)

            if(azi(i).ge.0.d0 .and. touse(i)(2:2).ne.'m') then
               write(string(46:57),'(f6.2,1x,f5.2)') azi(i),azis(i)
               string(72:72) = 'A'
            else if(touse(i)(2:2).eq.'m') then
               string(72:72) = 'x'
            else
               string(72:72) = '_'
            endif

            if(p(i).ge.0.d0 .and. touse(i)(3:3).ne.'m') then
               write(string(59:69),'(f5.2,1x,f5.2)') p(i),ps(i)
               string(73:73) = 'S'
            else if(touse(i)(3:3).eq.'m') then
               string(73:73) = 'x'
            else
               string(73:73) = '_'
            endif

            if(amplit(i).gt.0.d0 .and. period(i).gt.0.d0) then
               write(string(79:97),'(f6.3,1x,f12.2)')period(i),amplit(i)
            endif

            if(snr(i).gt.0.d0) then
               write(string(99:105),'(f7.2)') snr(i)
            else
               string(99:105) = '       '
            endif

            if(arid(i).ne.'        ') then
               write(string(106:114),'(1x,a8)') arid(ii)
            endif

            if((dabs(tt2(i)-tts(i)).gt.0.001d0) .and. itresw.eq.1) then
               write(string(115:120),'(1x,f5.2)') tt2(i)
            endif

            write(31,'(a)') trim(string)

         endif

      endif

      if(sta(iev(indx(i))).ne.stato.or.i.eq.nobs) then

        ntext0 = ntext

        if(i.eq.nobs) then 

          ntext        = ntext + 1
          arr(ntext)   = sngl(tt(indx(i)))
          text2(ntext) = text(indx(i))

          ntext0       = ntext
          if(sta(iev(indx(i))).ne.stato) ntext0 = ntext - 1

        endif

        call indexx(ntext0,arr,indx2)

        do 452 j=1,ntext0
        textout = trim(text2(indx2(j)))
        if(textout(1:5).ne.'     ') then
           if(output) write(11,'(a)') textout(1:lenons)
           
           if(isf_out) call isf_out_line(12,textout)
           if(itest.eq.20) then
              itest = write_isf_error(typctl) 
           endif

        endif
452     continue

        if(i.eq.nobs) then
           if(ntext0.eq.ntext-1) then
             textout = trim(text(indx(i)))
             if(textout(1:5).ne.'     ') then
               if(output) write(11,'(a)') textout(1:lenons)

               if(isf_out) call isf_out_line(12,textout)
               if(itest.eq.20) then
                  itest = write_isf_error(typctl)
               endif

             endif
           endif
           go to 453
        endif

        stato = sta(iev(indx(i)))
        ntext = 1

      else

        ntext = ntext + 1

      endif

      arr(ntext)   = sngl(tt(indx(i)))
      text2(ntext) = text(indx(i))

453   continue
      
455   continue

      if(nobsa.gt.0) then
         samean  = samean/dble(nobsa)
         sarmean = sarmean/dble(nobsa)
         rmsazi  = dsqrt(rmsazi/dble(nobsa))
      endif
      if(nobsp.gt.0) then
         spmean  = spmean/dble(nobsp)
         sprmean = sprmean/dble(nobsp)
         rmsp    = dsqrt(rmsp/dble(nobsp))
      endif

      if (ibad.ge.3 .and. .not.last) go to 470

c
      if(ndt.eq.0) go to 466
      if(output) then
        write(11,'(/''Defining travel-time differences:''/)')
        write(11,'('' Stat  Delta  Phases'',11x,''Observed   Res''/)')
      endif

      i2 = 0
      sdmean  = 0.d0
      sdrmean = 0.d0
      rmsdt   = 0.d0

      do 461 i = 1,nobs-1

      if(used(i)(4:4).ne.'D') go to 461

      do 460 j = i+1,nobs

         if(used(j)(4:4).ne.'D') go to 460

         if(sta(iev(i)).ne.sta(iev(j))) go to 460
         if(phaseu(i).eq.phaseu(j)) go to 460

         do  458 i3 = 1, ndt
            if(idtu(i3).eq.(j*i + j+i)) go to 459
458      continue
         go to 460

459      i2    = i2 + 1
         arr(i2) = sngl(del(iev(i)))

         dtth  = ttt(j) - ttt(i)
         dtobs = tt(j) - tt(i)
         dtres = dtobs - dtth

         sdmean  = sdmean  + dtres
         sdrmean = sdrmean + dabs(dtres)
         rmsdt   = rmsdt   + q2(dtres)

         fmis   = dtres/dpythag(ttu(i),ttu(j))
         dmisfl = dmisfl + dabs(fmis)
         dmisf  = dmisf + q2(fmis)
         ndmisf = ndmisf + 1
 
         art = trim(phaseu(j))//' - '//trim(phaseu(i))

c         if(typctl.gt.5) then
c            print *,sta(iev(j)),del(iev(j)),
c     +               trim(art),dtobs,dtres
c         endif

         if(output) then
            statw = sta(iev(i))
            if(touse(i)(9:9).eq.'*' .or. touse(j)(9:9).eq.'*') then
              chgcas = lowcas(statw)
              statw = chgcas(1:5)
            endif

            if (kmout) then
               write(text(i2),'(a5,f8.2,1x,a16,f9.3,f8.3)') 
     +               statw,delk(iev(i)),art,dtobs,dtres
            else
               write(text(i2),'(a5,f8.3,1x,a16,f9.3,f8.3)') 
     +               statw,del(iev(i)),art,dtobs,dtres
            endif
         endif

460   continue

461   continue

      sdmean  = sdmean / dble(ndmisf)
      sdrmean = sdrmean / dble(ndmisf)
      rmsdt   = dsqrt(rmsdt  / dble(ndmisf))

      call indexx(i2,arr,indx)

      if(output) then
         do 463 i=1,i2
         write(11,'(a)') trim(text(indx(i)))
463      continue
      endif

466   continue

      if(output) then
         write(11,'(/''Number of usable stations: '',i4)') nstata
      endif

      miteras = miteras + iter
      if(output) then
         if(miteras.gt.iter) then
            write(11,'(/''Total number of iterations: '',i5)') miteras
         endif
      endif

      if(output .and. nstat.gt.1) then
         if(.not.gapobs) then
              write(11,'(/''Maximum azimuthal gap of defining '',
     +        ''observations: '',f5.1,'' -> '',f5.1,'' [deg]'',
     +        '' = '',f5.1,'' [deg]'')') d1azi,d2azi,dazgap
         else
              write(11,'(/''Maximum azimuthal gap for all '',
     +        ''observing stations: '',f5.1,'' -> '',f5.1,'' [deg]'',
     +        '' = '',f5.1,'' [deg]'')') d1azi,d2azi,dazgap
         endif
      endif

      if(output .and. nstat.gt.2) then
         if(.not.gapobs) then
              write(11,'(/''Maximum secondary azimuthal gap of '',
     +        ''defining observations: '',f5.1,'' -> '',f5.1,
     +        '' [deg] = '',f5.1,'' [deg]'')') d1azi2,d2azi2,dazgap2
         else
              write(11,'(/''Maximum secondary azimuthal gap for '',
     +        ''all observing stations: '',f5.1,'' -> '',f5.1,
     +        '' [deg] = '',f5.1,'' [deg]'')') d1azi2,d2azi2,dazgap2
         endif
      endif

c
c     output of mean errors
c
      if(output) then
         write(11,'(/''Residuals of defining data'',10x,
     +               ''RMS    MEAN-ERROR      MEAN'')')

         if(nobst.eq.1) 
     +   write(11,'(i6,'' onset time              : '',f8.3,x,2(3x,
     +         f8.3),''  [s]'')') nobst,rms,strmean,stmean

         if(nobst.gt.1) 
     +   write(11,'(i6,'' onset times             : '',f8.3,x,2(3x,
     +         f8.3),''  [s]'')') nobst,rms,strmean,stmean

         if(nobsa.eq.1) 
     +   write(11,'(i6,'' backazimuth value       : '',f8.3,x,2(3x,
     +         f8.3),''  [deg]'')') nobsa,rmsazi,sarmean,samean

         if(nobsa.gt.1) 
     +   write(11,'(i6,'' backazimuth values      : '',f8.3,x,2(3x,
     +         f8.3),''  [deg]'')') nobsa,rmsazi,sarmean,samean

         if(nobsp.eq.1)
     +   write(11,'(i6,'' ray parameter           : '',f8.3,x,2(3x,
     +         f8.3),''  [s/deg]'')') nobsp,rmsp,sprmean,spmean

         if(nobsp.gt.1)
     +   write(11,'(i6,'' ray parameters          : '',f8.3,x,2(3x,
     +         f8.3),''  [s/deg]'')') nobsp,rmsp,sprmean,spmean

         if(ndmisf.eq.1) 
     +   write(11,'(i6,'' travel-time difference  : '',f8.3,x,2(3x,
     +         f8.3),''  [s]'')') ndmisf,rmsdt,sdrmean,sdmean

         if(ndmisf.gt.1) 
     +   write(11,'(i6,'' travel-time differences : '',f8.3,x,2(3x,
     +         f8.3),''  [s]'')') ndmisf,rmsdt,sdrmean,sdmean

         write(11,'(/''Weighted RMS of onset times (ISC type): '',
     +              f8.3,'' [s]''/)') rmsisc

         write(11,'(''Weighted misfit of input data'',9x,
     +              ''L1      L2'')')

         if(ntmisf.ge.1) then
           tmisf1 = dsqrt(tmisf/ntmisf)
           tmisfl1 = tmisfl/ntmisf
c          PRINT *,'tmisf tmisfl ntmisf tmisf1 tmisfl1'
c          PRINT *,tmisf,tmisfl,ntmisf,tmisf1,tmisfl1 
           write(11,'(i6,'' onset times             :'',2(x,f8.3))')
     +     ntmisf,tmisfl1,tmisf1
         endif

         if(namisf.ge.1) then
           amisf1 = dsqrt(amisf/namisf)
           amisfl1 = amisfl/namisf
           write(11,'(i6,'' backazimuth values      :'',2(x,f8.3))')
     +     namisf,amisfl1,amisf1
         endif

         if(npmisf.ge.1) then
           pmisf1 = dsqrt(pmisf/npmisf)
           pmisfl1 = pmisfl/npmisf
           write(11,'(i6,'' ray parameters          :'',2(x,f8.3))')
     +     npmisf,pmisfl1,pmisf1
         endif

         if(ndmisf.ge.1) then
           dmisf1 = dsqrt(dmisf/ndmisf)
           dmisfl1 = dmisfl/ndmisf
           write(11,'(i6,'' travel-time differences :'',2(x,f8.3))')
     +     ndmisf,dmisfl1,dmisf1
         endif

         nwmisf = ntmisf + namisf + npmisf + ndmisf
         wmisf = dsqrt((tmisf + amisf + pmisf + dmisf)/nwmisf)
         wmisfl = (tmisfl + amisfl + pmisfl + dmisfl)/nwmisf
         write(11,'(i6,'' misfit over all         :'',2(x,f8.3))')
     +   nwmisf,wmisfl,wmisf

      endif
c
c     output of one line with all calculated source parameters and
c     quality parameters
c

      if(output) then
         write(11,'(/''T0'',25x,''LAT'',6x,''LON'',7x,''Z'',5x,
     +      ''VPVS'',3x,''DLAT'',5x,''DLON'',6x,''DZ'',7x,''DT0'',4x,
     +      ''DVPVS DEF'',4x,''RMS'' )' )
      endif
      write(*,'(/''T0'',25x,''LAT'',6x,''LON'',7x,''Z'',5x,
     +      ''VPVS'',3x,''DLAT'',5x,''DLON'',6x,''DZ'',7x,''DT0'',4x,
     +      ''DVPVS DEF'',4x,''RMS'' )' )
      
      call fetoh(tome,idum,yy,mon,mm,dd,idoy,hh,mi,sec)
 
      isec1 = nint(sec*1000)
      isec  = isec1/1000
      isec2 = isec1-isec*1000

      if(czo.eq.'D') then

        if(output) then
           write(11,'(i4,''-'',i2.2,''-'',i2.2,3i3.2,''.'',i3.3,
     +        2f9.3,f8.2,f7.2,2f9.4,f8.2,f9.3,f7.2,i5,f9.3)') 
     +        yy,mon,dd,hh,mi,isec,isec2,elatmg,elonm
     +        ,zo,vpvsm,sdlatg,var2(3),var2(4),var2(1),sdvpvs,in,rms

        endif

        write(*,'(i4,''-'',i2.2,''-'',i2.2,3i3.2,''.'',i3.3,2f9.3,
     +           f8.2,f7.2,2f9.4,f8.2,f9.3,f7.2,i5,f9.3)') 
     +           yy,mon,dd,hh,mi,isec,isec2,elatmg,elonm
     +           ,zo,vpvsm,sdlatg,var2(3),var2(4),var2(1),sdvpvs,in,rms
        
      else if(czo.eq.'F' .or. czo.eq.'B') then

        if(output) then

              write(11,'(i4,''-'',i2.2,''-'',i2.2,3i3.2,''.'',i3.3,
     +          2f9.3,f8.2,f7.2,2f9.4,''  Fixed '',f9.3,f7.2,i5,f9.3)') 
     +          yy,mon,dd,hh,mi,isec,isec2,elatmg,elonm
     +          ,zo,vpvsm,sdlatg,var2(3),var2(1),sdvpvs,in,rms

        endif

        write(*,'(i4,''-'',i2.2,''-'',i2.2,3i3.2,''.'',i3.3,2f9.3,
     +           f8.2,f7.2,2f9.4,''  Fixed '',f9.3,f7.2,i5,f9.3)') 
     +           yy,mon,dd,hh,mi,isec,isec2,elatmg,elonm
     +           ,zo,vpvsm,sdlatg,var2(3),var2(1),sdvpvs,in,rms
      endif

      if (isf_out) then
         write(12,'(/,''STOP'')')
         close(12)
      endif

      if (new_in) close(31)

      if(isf_ref.ne.' ' .and. output .and. isf_in) then

         call depi(dlati,dloni,elatmg,elonm,del3,dk,ep2,ep1,d2km)

         write(11,'(/,''Distance to ISF Ref ( '',a10,f9.4,f10.4,'' ):'',
     +           f8.2,'' km; DZ:'',f7.1,'' km, AZI:'',f6.1,'' MA: '',
     +           f7.1,'' MI: '',f7.1,'' AZI: '',f6.1)') 
     +           author2,dlati,dloni,dk,ddepi-zo,ep2,elmax,elmin,eazi
      endif

      if(ref_eve .and. output) then

         call depi(dlati,dloni,elatmg,elonm,del3,dk,ep2,ep1,d2km)

         write(11,'(/,''Distance to Reference Event ('',f9.4,f10.4,
     +           '' ):'',f8.2,'' km; DZ:'',f7.1,'' km, AZI:'',f6.1,
     +           '' MA: '',f8.1,'' MI: '',f7.1,'' AZI: '',f6.1)') 
     +           dlati,dloni,dk,ddepi-zo,ep2,elmax,elmin,eazi
      endif
      
      if( output .and. iloc .and. modout) then

        if(imo.ge.3) then

           if(.not.kmout) then
             write(11,'(/,''CRUST 1.0 model for source-region (max. '',
     +              ''delta ='',f6.2,'' deg):'',/,''    DEPTH   '',
     +              ''    VP        VS    DISCON'')') rmax
           else
             radkm = rmax*deg2rad*radloc(elatmg,1)
             write(11,'(/,''CRUST 1.0 model for source-region (max. '',
     +              ''delta ='',f8.1,'' km):'',/,''    DEPTH   '',
     +              ''    VP        VS    DISCON'')') radkm
           endif
        
           itrue = 0
           elatc = elatmg
           elonc = elonm
           inum  = 1
           call get_mod_c10(itrue,inum,typctl,ierr)
           write(11,'(3f10.3,a6)') (z(i),v0(1,i),v0(2,i),
     +           azo(i),i=1,jmod)

        else if(imo.eq.1 .or. imo.eq.2) then

           if(.not.kmout) then
              write(11,'(/,''Local model '',a,'' used (max. delta ='',
     +              f6.2,'' deg):'',/,''    DEPTH       VP'',
     +              ''        VS'',''    DISCON'')') 
     +              trim(filloc),rmax
              write(11,'(3f10.3,a6)') (z(i),v0(1,i),v0(2,i),
     +              azo(i),i=1,jmod)
           else
              radkm = rmax*deg2rad*radloc(elatmg,1)
              write(11,'(/,''Local model '',a,'' used (max. delta ='',
     +              f8.1,'' km):'',/,''    DEPTH       VP'',
     +              ''        VS'',''    DISCON'')') 
     +              trim(filloc),radkm
              write(11,'(3f10.3,a6)') (z(i),v0(1,i),v0(2,i),
     +              azo(i),i=1,jmod)
           endif
        endif

        if(rzv(1).gt.0.d0) then
          write(11,'(/,''P velocity in source depth: '',f6.2)') rzv(1)
          if(rzv(2).gt.0.d0) then
            write(11,'(''S velocity in source depth: '',f6.2)') rzv(2)
          endif
        endif

      endif

c
c     Now let us try to locate the event better without a fixed depth.
c

c     print *,'czo ',czo,' iterz ',iterz,' zoflag ',zoflag

470   if(czo.eq.'B' .and. ((.not. zoflag) .or. iterz.lt.2)) then
         zoflag=.true.
         czo ='D'
         if(zo.le.0.1d0) zo = 33.d0
         rs(4)  = sdzo
         rs(1)  = dpythag(rs(1),25.d0)

         if(idelw.gt.1) idelw = idelw - 2

         miteras = miteras + iter
         iter = 0
         itso = 0
         nextiter = 0
         nextiter1= 0
         dchang   = dchang0
         iteraz = 0
         rmso   = 9999.d0
         rmsold   = 9999.d0
         datmax0  = 9999.d0
         dazim1 = dazim0*2.0d0
         dpam1  = dpam0*2.0d0
         dtmin  = 9999.d0
         dtp    = dtp0
         dts    = dts0
         dtm2   = dts
         check  = 9999.d0
         setcheck = setcheck1
         setcheck2 = 15.d0*setcheck

         go to 100
      endif

      go to 9999

9998  continue
      if(output) then
         write(11,'(/''Could not find a better or non-oscillating'',
     +            '' solution without fixed depth'')')
      endif
      write(*,'(/''Could not find a better or non-oscillating'',
     +            '' solution without fixed depth'')')

9999  continue

CCC      do 99999 ilst = 1,nstat
CCC
CCC      write(*,'('' ray path '',4f12.5,f8.1,f8.3)')
CCC     +      rays(ilst,1),rays(ilst,2),rays(ilst,3),rays(ilst,4),
CCC     +      rays(ilst,5),rays(ilst,6)
CCC
CCC99999 continue

      if(output) close(11)

      stop

c     end program HYPOSAT_6_0d
      end 
