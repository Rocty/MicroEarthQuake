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
      program HYPOMOD_2_0b

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      character  version*25
      parameter (version='HYPOMOD Version 2.0b     ') 

c
c     last changes:  9 October 2018
c
c----------------------------------------------------------------------
c
c     Short desciption - for more details see HYPOMOD-Manual 
c
c     This program calculates for a given seismic event all residuals
c     observed travel times, backazimuth, and slowness values.
c
c     All input and output files are identical to hyposat.
c     See HYPOSAT manual for details. However some features are just
c     ignored because we do not invert any data!
c
c     HYPOMOD 2.0 is based on HYPOSAT 6.0
c
c--------------------------------------------------------------------
c
c               Program History
c
c         see file hypomod_history!
c
c--------------------------------------------------------------------
c
c            Main Program dependencies:
c
c     calls:     delazd, depi, fetoh, fhtoe, findrange, get_station,
c                hyposat_cross, ellip, hyposat_gmi, indexx, 
c                tauget_mod, testphase, hyposat_geo, ttloc, get_mod_c10,
c                ellcal, plane, dpythag, mult_ons,
c                tauget_ray, get_mod_global
c
c     functions: alpha1, alpha2, convlat, phase_type, crustc, 
c                dirdel, q2, radloc
c                file_checkpara, read_event_id, read_origin, 
c                read_origin_head, read_phase_head, read_phase, 
c                lowcas, uppcas
c
c
c     data exchange by common blocks in include files:
c                
c                include 'ttimes.h'
c                include 'model.h'
c                include 'modelc.h'
c                include 'modelg.h'
c
c     PARAMETER settings for: mstat, mread, mvar 
c
c****6789012345678901234567890123456789012345678901234567890123456789012
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c     Functions called: variable definitions
c
      real*8           alpha1, alpha2, convlat, crustc, dirdel, q2, 
     +                 radloc, dpythag
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
     +          magfile*120,magmlfile*120,outputfile1*120
c
c     mread = maximum number of phases (hypomod-in file)
c
c     mrd2 = maximum number of observations per station
c
      parameter  (mread = 4000, mrd2 = 100)

      character phase(mread)*8,phaseu(mread)*8,phid*8,
     +          phid2*8,text2(mrd2)*155,phid1*8,phase_t*1,
     +          string*155,touse(mread)*9,touse0*9,phidr*8,
     +          o_string*155,textout*155,text(mread)*155,
     +          arid(mread)*8,statcorstr*80,texth*155

      dimension azi(mread),tt(mread),p(mread),azis(mread),tts(mread),
     +          ps(mread),period(mread),amplit(mread),
     +          tt2(mread),ttu(mread),iev(mread),indx(mread),
     +          indx2(mrd2),snr(mread),emeran(mread)

ccc  dimension  rays(mread,6)

      real*4    arr(mread),epiaz(mread),epiaz2(mread),dazgap,d1azi,
     +          d2azi,dazgap2,d1azi2,d2azi2

c
      dimension var2(4)

c
c     variables for travel-time calculations
c
      include 'ttimes.h'

      dimension dpdh(mphas) 

      character art*16,modnam*20,
     +          modnam2*20,modn*20,modnam3*20,modnam4*20, mtyp0*3

      real*4 rzo,rdel,ecor,fla1,razi,pa,rzo1,
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
     +          phisf*10, isf_ref*10, phidd*10, author2*10, 
     +          corid*8, corid2*8, cpick*1

      real*4    rdum, rpa, ramp, rper, rsnr, rlati, rloni, rdepi,
     +          rdmi, rdma, rdum0

c
c     Functions called from the ISF software libarary
c
      integer   read_event_id, read_origin, read_origin_head, 
     +          read_phase_head, read_phase, write_isf_error

c
c     other variables
c
      integer   yy,mon,dd,hh,idoy,ierr,typctl,mi,idum, isreg, regnum,
     +          y00, mon00, d00, h00

      character mm*4,name*48
      real*8    lat,lon,dlati,dloni,ddepi
      real*4    elevs, sec, rlat, rlon, smag

      character title*140, region*80, czo*1, magtypp*3, 
     +          magtyps*6, magtypml*7, statmag*2

c

      dimension ttt(mread)

      logical   vlflag  , stcorfl, iloc, surf,
     +          diffflag, single , output, modout, 
     +          magflag, mod2flag, plflag,
     +          conr, rayok, rayokf, mod3flag, mod4flag, 
     +          kmout, thbaz, thray, 
     +          isf_in, fixinp, ref_eve, 
     +          pflag, lgflag, sgflag, check_mult, aziflag,
     +          sloflag, gapobs, isf_epi, isf_dep, 
     +          lrmsisc, firston, firstph, firsts, ldepth0,
     +          old_syntax, emerout, larid

c
c     some constants and initial or default values
c

      pi      = 4.d0*datan(1.d0)
      deg2rad = pi / 180.d0
      rad2deg = 180.d0 / pi

      rearth  = 6371.d0
c
      dtmin   = 9999.d0

      ttray   = 0.d0

      cevid  = '999999999'
      corid  = '_'
      corid2 = ' '
      ievid  = 0
      AUTHOR = 'HYPOSAT'

      vlflag   = .true.
      stcorfl  = .false.
      diffflag = .true.
      iloc     = .false.
      single   = .false.
      output   = .false.
      isf_in   = .false.
      isf_epi  = .false.
      isf_dep  = .false.
      modout   = .false.
      magflag  = .false.
      mod2flag = .false.
      mod3flag = .false.
      mod4flag = .false.
      plflag   = .true.
      rayok    = .false.
      rayokf   = .false.
      kmout    = .false.
      thray    = .false.
      thbaz    = .false.
      fixinp   = .false.
      locgeo   = .false.
      locsta   = .false.
      ref_eve  = .false.
      pflag    = .false.
      lgflag   = .false.
      sgflag   = .false.
      check_mult = .false.
      aziflag  = .true.
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
      outputfile  = 'hypomod-out'
      inputfile   = 'hypomod-in'
      old_syntax = .false.
      statcorfile = ' '
      vpl = 5.8d0 
      vsl = 3.46d0
      zo1 =  0.d0
      czo = 'F'
      typctl = 4
      islow = 1
      indph0 = 3333
      epilat0 = -999.d0
      epilon0 = -999.d0
      tome0   = -2840140801.d0
      dismaxst = 21001.d0
      disminst = -1.d0

      elmax  = -999.d0
      elmin  = -999.d0
      ieazi  = -999

      dtphmax = 5.d0

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
     +                  'hyposat-parameter file'
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

      if(string(1:21).eq.'STARTING SOURCE DEPTH') then
          read (string(38:),*) zo1
          go to 1
      endif

      if(string(1:24).eq.'STARTING SOURCE LATITUDE') then
          abc = -999.0d0
          read (string(38:),*) abc
          if(abc.ge.-90.d0 .and. abc.le.90.d0) epilat0 = abc
          go to 1
      endif

      if(string(1:25).eq.'STARTING SOURCE LONGITUDE') then
          abc = -999.0d0
          read (string(38:),*) abc
          if(abc.ge.-180.d0 .and. abc.le.180.d0) epilon0 = abc
          go to 1
      endif

      if(string(1:30).eq.'AZIMUTHAL GAP FOR OBSERVATIONS') then
          intinp = 0
          gapobs = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) gapobs = .true.
          go to 1
      endif

      if(string(1:12).eq.'OUTPUT LEVEL') then
          read (string(38:),*) typctl
          if (typctl.lt.0)     typctl = 0
          if (typctl.gt.10)    typctl = 10
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

      if(string(1:16).eq.'ISC-TYPE ISF RMS') then
          intinp = 0
          lrmsisc = .false.
          read (string(38:),*) intinp
          if(intinp.eq.1) lrmsisc=.true.
          go to 1
      endif

      if(string(1:16).eq.'OUTPUT FILE NAME') then
          outputfile1 = 'hypomod-out'
          if(string(38:38).ne.' ' .and. string(38:38) .ne.'_')
     +                   read (string(38:),*) outputfile1
          if(outputfile1.ne.'hypomod-out' .or.
     +       outputfile1.ne.'hyposat-out'   ) outputfile=outputfile1
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

      if(string(1:17).eq.'LG GROUP-VELOCITY') then
         read (string(38:),*) vlg
         if(vlg.le.0.d0) vlg = vlg0
         go to 1
      endif

      if(string(1:17).eq.'RG GROUP-VELOCITY') then
          read (string(38:),*) vrg
          if(vrg.le.0.d0) vrg = vrg0
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

      if(string(1:17).eq.'T PHASE GROUP-VEL') then
          read (string(38:),*) vt
          if(vt.le.0.d0) vt = vt0
          go to 1
      endif

      if(string(1:18).eq.'IS PHASE GROUP-VEL') then
          read (string(38:),*) vi
          if(vi.le.0.d0) vi = vi0
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
         print *,'vpl = ',vpl
         print *,'vsl = ',vsl
         print *,'vrg = ',vrg
         print *,'vlg = ',vlg
         print *,'vlq = ',vlq
         print *,'vlr = ',vlr
         print *,'vt  = ',vt 
         print *,'zo1   = ',zo1
         print *,'epilat0 = ',epilat0
         print *,'epilon0 = ',epilon0
         print *,'typctl = ',typctl
         print *,'islow = ',islow
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
            indph  = 0
            elatc  = 0.d0
            elonc  = 0.d0
            sdep   = 0.d0
            rmax = 0.d0
            jmod = 0

            nphas = 0
            call ttloc(rzo,rdel,czo,nphas,ttc,dtdd,dtdh,dpdh,dddp,
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

      if (stcorfl) open(13,file=statcorfile)

      zo = zo1

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
         endif
         if(isf_dep) then
            zo1     = ddepi
            zo      = zo1
         endif

34       read (10,'(a)',end=9999) string
         itest = read_phase_head(string)
         if(itest.eq.20) go to 34

      else
         if(output) write (11,'(a,/)') trim(title)
         print *,'EVENT ',title
      endif

      isnr = 0

      sdpmean = 0.d0
      nsdp    = 0
      sdsmean = 0.d0
      nsds    = 0

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

             print *,'Reading error with onset:',string
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
      endif

      touse(ii)=touse0

      if(tt2(ii).le.0.d0) tt2(ii) = tts(ii)

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
        iev(ii) = j
        istat = j

        go to 11

      endif

10    continue

11    continue

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

      nobs  = ii
      nstat = istat

      if(nobs.eq.1) then

         rzo1  = 0.
         nphas = 0

         rdel1 = 3.
         call tauget_mod(rzo1,rdel1,nphas,phcd,ttc,dtdd,
     +                   dtdh,dddp,modnam)

         single = .true.

         if(typctl.gt.4) print *, 'Case: single array observation!'

      endif

      if (check_mult) then
         nobs0 = nobs
         call mult_ons(nobs,nobs0,iev,phase,tt,tts,tt2,azi,azis,p,ps,
     +                 touse,amplit,period,dtphmax,mread,arid)
         if(nobs.ne.nobs0) nobs = nobs0
      endif

      larid = .false.
      do 15 i = 1,nobs
      if(arid(i).ne.'        ') larid = .true.
15    continue

      tome = tome0

      zo = zo1

      elonm  = epilon0
      elatmg = epilat0
      elatm  = convlat(elatmg,1)

      ndt = 0

      if(output) then
         if(iloc) then
            if( (mod2flag .or. mod3flag .or. mod4flag) .and.
     +          (imodn(2).gt.1 .or. imodn(2).gt.1 .or. imodn(3).gt.1)
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
     +          (imodn(2).gt.1 .or. imodn(2).gt.1 .or. imodn(3).gt.1)
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
      if(output) then

         write(11,'(/''The source parameters:'')')

         if(kmout) then
            if(disminst.gt.0.d0 .and. dismaxst.ge.21000.d0) then
               write(11,'(/''For observations at distances '',
     +               ''larger than'',f8.1,'' [km]'')') disminst
            endif
            if(disminst.le.0.d0 .and. dismaxst.lt.21000.d0) then
               write(11,'(/''For observations at distances '',
     +               ''smaller than'',f8.1,'' [km]'')') dismaxst
            endif
            if(disminst.gt.0.d0 .and. dismaxst.lt.21000.d0) then
               write(11,'(/''For observations at distances '',
     +               ''between'',f8.1,'' and'',f8.1,'' [km]'')') 
     +               disminst,dismaxst
            endif
         else
            if(disminst.gt.0.d0 .and. dismaxst.ge.180.d0) then
               write(11,'(/''For observations at distances '',
     +               ''larger than'',f6.1,'' [deg]'')') disminst
            endif
            if(disminst.le.0.d0 .and. dismaxst.lt.180.d0) then
               write(11,'(/''For observations at distances '',
     +               ''smaller than'',f6.1,'' [deg]'')') dismaxst
            endif
            if(disminst.gt.0.d0 .and. dismaxst.lt.180.d0) then
               write(11,'(/''For observations at distances '',
     +               ''between'',f6.1,'' and'',f6.1,'' [deg]'')') 
     +               disminst,dismaxst
            endif
         endif

         write(11,'(/''Source time  :'',i5,4i3.2,f7.3,
     +           '' [s]'')')  yy,mon,dd,hh,mi,sec
         write(11,'(''        or'',12x,f16.3,
     +           '' [s]'')') tome

         isec1 = nint(sec*1000)
         isec  = isec1/1000
         msec  = isec1-isec*1000
         write(11,'(''        or'',7x,i4,''-'',i3.3,'':'',
     +         3(i2.2,''.''),i3.3,'' [s]''//)') 
     +         yy,idoy,hh,mi,isec,msec
         
      endif

      if(output) then
        write(11,'(''Epicenter lat:'',14x,f10.4,
     +           '' [deg]'')')  elatmg
        write(11,'(''Epicenter lon:'',14x,f10.4,
     +           '' [deg]'')')  elonm
        write(11,'(''Source depth :'',15x,f7.2,
     +               ''   [km]''/)') zo
      endif

c
c     let us now calculate the final residuals and print them out
c
      stmean    = 0.d0
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

      rzv(1)    = -999.d0
      rzv(2)    = -999.d0

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

      ttu(i) = tts(i)

      ttobs  = tt(i)

      phid1 = phase(i)
      phid2 = phase(i)
      phid  = phase(i)

      call fetoh(ttobs,idum,yy,mon,mm,dd,idoy,hh,mi,sec)

      epiaz(i) = -999.0
      epiaz2(i) = -999.0
      emeran(i) = -999.d0

      cdum = ' '

      if(diffflag) then

         chgcas = uppcas(phase(i)(1:1))
         cdum = chgcas(1:1)
         if( cdum.ne.'P' .and. cdum.ne.'S' .and. cdum.ne.'R' .and.
     +       cdum.ne.'L' ) go to 413
c
c        Mark all possible phases, which can be used as part of a defining
c        travel-time difference measure.
c
         do 412 j = i+1,nobs

            chgcas = uppcas(phase(j)(1:1))
            cdum = chgcas(1:1)

            if( cdum.ne.'P' .and. cdum.ne.'S' .and. cdum.ne.'R' .and.
     +          cdum.ne.'L' ) go to 412
            if(sta(iev(i)).ne.sta(iev(j))) go to 412
            if(phase(i).eq.phase(j)) go to 412

            if(touse(i)(4:4).eq.'D' .and. touse(j)(4:4).eq.'D') then
               ndt = ndt + 1
            endif

412      continue

      endif

413   if (sta(iev(i)).ne.stato) then

         istad(iev(i)) = 0

         stato = sta(iev(i))

         call depi(stala(iev(i)),stalo(iev(i)),elatmg,elonm,
     +              del(iev(i)),delk(iev(i)),azie(iev(i)),baz(iev(i)),
     +        d2km)
         rzo   = sngl(zo)
         rdel  = sngl(del(iev(i)))
         rdelk = sngl(delk(iev(i)))

         if(rdel.lt.rdmi) rdmi = rdel
         if(rdel.gt.rdma) rdma = rdel

         fla1 = sngl(deg2rad*(90.d0-elatm))
         razi = sngl(azie(iev(i)))

         nphas = 0

         if (mod2flag .and. touse(i)(7:7).eq.'2') then
           imodn(2) = 2
           modn = modnam2
         else if (mod3flag .and. touse(i)(7:7).eq.'3') then
           imodn(3) = 2
           modn = modnam3
         else if (mod4flag .and. touse(i)(7:7).eq.'4') then
           imodn(4) = 2
           modn = modnam4
         else
           if(iloc) imod2 = 1
           imodn(1) = 2
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

      dpa   = 0.d0
      phid1 = phase(i)
      ttres  = -9999.d0
      ttres1 =  9999.d0
      pares  =  -1000.d0
      surf = .false.

      if(sglgdis.gt.0.d0 .and. .not.fixinp) then
         if(phid.eq.'Sg' .and. delk(iev(i)).ge.sglgdis) phid = 'Lg'
         if(phid.eq.'Lg' .and. delk(iev(i)).lt.sglgdis) phid = 'Sg'
      endif

      if(phid.eq.'Rg') then
         if(delk(iev(i)).le.400.d0 .or.fixinp) then
            surf = .true.
            vsurf = vrg
         else
            phid = 'LR'
         endif
      endif
                
      if(phid.eq.'Lg' ) then
         if(delk(iev(i)).le.3000.d0 .or.fixinp) then
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
     +   'b ') .and. rdel.gt.30. .and. .not.fixinp) phid(2:3)='  '

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
     +   (phid1(3:3).eq.' ' .or. phid1(2:4).eq.'dif' .or.
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

      if(phid.eq.'PKPdf'.and. .not.fixinp) then
         
         if(phid1(1:5).eq.'PKiKP' .and.
     +                       del(iev(i)).ge.90.d0) phid='PKiKP'

      endif

      if(icha.eq.0) then
         
        if(phid1(2:4).eq.'dif') then
          if(phid.eq.'P') phid='Pdif'
          if(phid.eq.'S') phid='Sdif'
        endif

        if(phid(1:5).eq.'PKPab' .and. del(iev(i)).lt.150.d0
     +     .and. .not.fixinp ) phid(1:5) = 'PKP  '

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
                 tcrust = crustc(phase_t,dpaa,zo,indr,typctl)

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
         trefl = 0.d0
         trefls = 0.d0
         treflp = 0.d0

         if(.not.surf .and. touse(i)(5:5).eq.'R' .and. loctt.eq.0) then 

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
                   zoc = zo

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
                treflp   = crustc(phase_t,dpaa,zoc,indr,typctl)

            endif

            trefl = trefls + treflp

            if (dabs(trefl).gt.0.d0.and.typctl.gt.6) then
               print *,'dirdel: ',phid,' azi ',azi0,' del ',del0
               print *,'  lat ',elatc,' lon ',elonc,' trefl ',trefl
            endif

         endif

415      continue

         ttt(i) = tome+ttc(j)+dble(ecor)+th+tcrust+trefl+statict

         if (typctl.gt.8) then
           print *,'ttt, t0, TT, ECOR, Height, Crust, Refl, Stat,Tobs'
           print
     +     *,ttt(i),tome,ttc(j),dble(ecor),th,tcrust,trefl,statict,
     +       ttobs
         endif

         ttres  = ttobs - ttt(i)

         if(ttres.gt.100.d0 .and. phase(i)(1:2).eq.'P1' .and. 
     +      phid(1:3).eq.'Pdi' .and. touse(i)(1:1).eq.'T') then
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

c          if(typctl.ge.8) then
c            print *,'i,ttt,ttobs,ttres,ecor,th,tcrust,touse,tts,tt2,ttu'
c            print *,i,ttt(i),ttobs,ttres,ecor,th,tcrust,touse(i),tts(i),
c    +               tt2(i),ttu(i)
c          endif

            if(touse(i)(1:1).eq.'T') then
               stmean  = stmean + ttres
               strmean = strmean + dabs(ttres)
               rms     = rms    + q2(ttres)
               rmsisc  = rmsisc + q2(ttres)*ttu(i)
               wisc    = wisc   + ttu(i)
               epiaz(i)= razi
               nobst   = nobst + 1
            endif

            pares = p(i)-dabs(dpa - dph)
            if(touse(i)(3:3).eq.'S') then
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
            pares2 = p(i)-dabs(dpa)
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

      if (imin.eq.1 .or. fixinp) go to 425

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

425   if(dabs(dtmin).le.15.d0 .or. touse(i)(1:1).eq.'T') then
         ttres  = dtmin
         ttres1 = ttresm
         pares  = pares2
         phid1  = phid2
         dpa    = dpa2
      else
         touse(i)(1:1)='_'
         touse(i)(4:4)='_'
         ttres  = -9999.d0
         ttres1 =  9999.d0
         pares  = -999.d0
         dpa    = 0.d0
         phid1  = ' '
      endif

430   continue

      phaseu(i) = phid1

      if(touse(i)(1:1).eq.'T' .and. ttres.gt.-999.d0 .and. 
     +   ttres.lt.9999.d0) then
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

      statw = stato
      if(touse(i)(9:9).eq.'*') then
        chgcas = lowcas(stato)
        statw = chgcas(1:5)
      endif

      write(text(i),'(a5,f8.3,f7.2,1x,a8,8x,2i3.2,f7.3,f8.3,
     +                f7.2,f8.2,f6.2,f7.2,1x,a5,a1)') 
     +                statw,rdel,razi,phase(i),hh,mi,sec,
     +                ttres,azi(i),azires,p(i),pares,touse(i)(1:5),imod
      
      if(kmout)             write(text(i)(6:13),'(f8.2)') rdelk

      if(phase(i).ne.phid1) write(text(i)(30:37),'(a8)') phid1

      if(ttres.lt.-99.999d0 .or. ttres.gt.999.999d0) then
         write(text(i)(51:58),'(1x,f7.2)') ttres
      endif

      if(dabs(ttres).gt.998.d0)  text(i)(51:58)='        '
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

           dmag = 0.d0
           statmag='  '

c           print *, phid1, ttres, ttres1

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

449        if(dmag.ne.0.d0) then
              if (dmag.gt. -0.99d0) then
                 write (text(i)(122:129),'(1x,f4.2,1x,a2)') dmag,statmag
              else
                 write (text(i)(122:129),'(1x,f4.1,1x,a2)') dmag,statmag
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

      rlat = sngl(elatmg)
      rlon = sngl(elonm)
      call hyposat_geo( rlat,rlon, isreg, regnum, region , ierr )

      if(output) then
         write(11,'(''Flinn-Engdahl Region ('',i4,'' ): '',a/)')
     +            isreg, trim(region)
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

      do 453 i=1,nobs

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

        endif
452     continue

        if(i.eq.nobs) then
           if(ntext0.eq.ntext-1) then
             textout = trim(text(indx(i)))
             if(textout(1:5).ne.'     ') then
               if(output) write(11,'(a)') textout(1:lenons)

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

      if(touse(i)(4:4).ne.'D') go to 461

      do 460 j = i+1,nobs

         if(touse(j)(4:4).ne.'D') go to 460

         if(sta(iev(i)).ne.sta(iev(j))) go to 460
         if(phaseu(i).eq.phaseu(j)) go to 460

         i2    = i2 + 1
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

      ndt = ndmisf

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
         write(11,'(/''Number of usable stations: '',i4)') nstat
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
      in = nobst + nobsa + nobsp + ndt
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

      if(output) then

         write(11,'(i4,''-'',i2.2,''-'',i2.2,3i3.2,''.'',i3.3,
     +        2f9.3,f8.2,f7.2,2f9.4,''  Fixed '',f9.3,f7.2,i5,f9.3)') 
     +        yy,mon,dd,hh,mi,isec,isec2,elatmg,elonm
     +        ,zo,0.,0.,var2(3),var2(1),0.,in,rms

      endif

      write(*,'(i4,''-'',i2.2,''-'',i2.2,3i3.2,''.'',i3.3,2f9.3,
     +         f8.2,f7.2,2f9.4,''  Fixed '',f9.3,f7.2,i5,f9.3)') 
     +         yy,mon,dd,hh,mi,isec,isec2,elatmg,elonm
     +         ,zo,0.,0.,var2(3),var2(1),0.,in,rms

      if(isf_ref.ne.' ' .and. output .and. isf_in) then

         call depi(dlati,dloni,elatmg,elonm,del3,dk,ep2,ep1,d2km)

         write(11,'(/,''Distance to ISF Ref ( '',a10,f9.4,f10.4,'' ):'',
     +           f8.2,'' km; DZ:'',f7.1,'' km, AZI:'',f6.1)')
     +           author2,dlati,dloni,dk,ddepi-zo,ep2
      endif

      if(ref_eve .and. output) then

         call depi(dlati,dloni,elatmg,elonm,del3,dk,ep2,ep1,d2km)

         write(11,'(/,''Distance to Reference Event ('',f9.4,f10.4,
     +           '' ):'',f8.2,'' km; DZ:'',f7.1,'' km, AZI:'',f6.1)') 
     +           dlati,dloni,dk,ddepi-zo,ep2
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

9999  continue

      if(output) close(11)

      stop

c     end program HYPOSAT_6_0
      end 
