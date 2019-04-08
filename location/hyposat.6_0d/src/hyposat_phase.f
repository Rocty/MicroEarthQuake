      function phase_type(phid0)
c
c           last corrections :  16 February 1997
c                               correcting range of characters to check
c
c                               November 3, 1997
c                               corrected for AB and AC-branches
c
c                               October 8, 2002
c                               corrected for P1 and S1 phase names
c                               and new IASPEI phase names like P4 ...
c
c                               27 August 2005
c
      character phid0*8, phase_type*1, phid*8
      integer icph

      phid = phid0
      phase_type=' '

      icph = len_trim(phid)

      if(icph.le.1) go to 10

      if(phid(icph:icph).eq.'1') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'2') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'3') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'4') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'5') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'6') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'7') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'8') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'9') then
          icph=icph-1
      endif

5     if(icph.le.1) go to 10

      if(icph.ge.4) then
         if(phid(icph-2:icph).eq.'dif') then
            icph=icph-3
            go to 5
         endif
         if(phid(icph-2:icph).eq.'pre'.or.phid(icph-2:icph).eq.
     +      'PRE') then
            icph=icph-3
            go to 5
         endif
      endif

      if(icph.ge.3) then

         if(phid(icph-1:icph).eq.'ab') then
            icph=icph-2
            go to 5
         endif
         if(phid(icph-1:icph).eq.'ac') then
            icph=icph-2
            go to 5
         endif
         if(phid(icph-1:icph).eq.'bc') then
            icph=icph-2
            go to 5
         endif
         if(phid(icph-1:icph).eq.'df') then
            icph=icph-2
            go to 5
         endif

      endif

      if(phid(icph:icph).eq.'n') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'g') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'b') then
         icph=icph-1
         if(icph.le.1) go to 10
      endif

      if(phid(icph:icph).eq."'") then
         icph=icph-1
         if(icph.le.1) go to 10
      endif
      if(phid(icph:icph).eq.'*') then
         icph=icph-1
      endif

10    if(icph.le.1) icph=1

      if(phid(icph:icph).eq.'P') phase_type='P'
      if(phid(icph:icph).eq.'S') phase_type='S'

      return
      end
c
c     subroutine testphase (phid0,icha)
c
c     some changes added for surface-reflections
c
c     3 August 1997, Johannes Schweitzer, NORSAR
c
c     13 February 2002: PKPab + PKPdif added
c
c     corrections 04 September 2005
c
c     2 Febrayry 2018: P'P'P' + S'S'S' added
c  
      subroutine testphase (phid0,icha)
      integer icha
      character phid*8,phid0*8,s*1
      logical flag

      flag = .false.
      phid = phid0
      s    = ' '

      if(phid0(1:1).eq.'p' .or. phid0(1:1).eq.'s') then
         flag = .true.
         s    = phid0(1:1)
         phid = phid0(2:)
      endif
      if(phid(1:1).eq.'P' .and. icha.ge.19) go to 50
      if(phid(1:1).eq.'S' .and. icha.ge.17) go to 50

      if(phid(1:4).eq.'PKP ') then
         phid='PKPdf'
         icha = icha+1
         goto 100
      endif
      if(phid(1:5).eq.'PKPdf') then
         phid='PKPdif'
         icha = icha+1
         goto 100
      endif
      if(phid(1:6).eq.'PKPdif') then
         phid='PKPbc'
         icha = icha+1
         goto 100
      endif
      if(phid(1:5).eq.'PKPbc') then
         phid='PKPab'
         icha = icha+1
         goto 100
      endif
      if(phid(1:5).eq.'PKPab') then
         phid='Pdif'
         icha = icha+1
         goto 100
      endif
      if(phid(1:4).eq.'Pdif') then
         phid='P'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'P') then
         phid='Pn'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'PmP') then
         phid='Pn'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'Pn') then
         phid='Pb'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'PbP') then
         phid='Pb'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'Pb') then
         phid='Pg'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'Pg') then
         phid='PKPdf'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'P1') then
         phid='Pg'
         icha = icha+1
         goto 100
      endif
      if(phid(1:5).eq."P'P' ") then
         phid="P'P'ab"
         icha = icha+1
         goto 100
      endif
      if(phid(1:6).eq."P'P'ab") then
         phid="P'P'bc"
         icha = icha+1
         goto 100
      endif
      if(phid(1:6).eq."P'P'bc") then
         phid="P'P'df"
         icha = icha+1
         goto 100
      endif
      if(phid(1:7).eq."P'P'P' ") then
         phid="P'P'P'ab"
         icha = icha+1
         goto 100
      endif
      if(phid(1:8).eq."P'P'P'ab") then
         phid="P'P'P'bc"
         icha = icha+1
         goto 100
      endif
      if(phid(1:8).eq."P'P'P'bc") then
         phid="P'P'P'df"
         icha = icha+1
         goto 100
      endif
      if(phid(1:6).eq.'PKKP  ') then
         phid='PKKPab'
         icha = icha+1
         goto 100
      endif
      if(phid(1:6).eq.'PKKPab') then
         phid='PKKPbc'
         icha = icha+1
         goto 100
      endif
      if(phid(1:6).eq.'PKKPbc') then
         phid='PKKPdf'
         icha = icha+1
         goto 100
      endif

      if(phid(1:4).eq.'SKS ') then
         phid='SKSdf'
         icha = icha+1
         goto 100
      endif
      if(phid(1:5).eq.'SKSdf') then
         phid='SKSac'
         icha = icha+1
         goto 100
      endif
      if(phid(1:5).eq.'SKSac') then
         phid='Sdif'
         icha = icha+1
         goto 100
      endif
      if(phid(1:4).eq.'Sdif') then
         phid='S'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'S') then
         phid='Sn'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'SmS') then
         phid='Sn'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'Sn') then
         phid='Sb'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'SbS') then
         phid='Sb'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'Sb') then
         phid='Sg'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'S1') then
         phid='Sg'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'Lg') then
         phid='Sg'
         icha = icha+1
         goto 100
      endif
      if(phid.eq.'Rg') then
         phid='SKSdf'
         icha = icha+1
         goto 100
      endif
      if(phid(1:5).eq."S'S' ") then
         phid="S'S'ac"
         icha = icha+1
         goto 100
      endif
      if(phid(1:6).eq."S'S'ac") then
         phid="S'S'df"
         icha = icha+1
         goto 100
      endif
      if(phid(1:7).eq."S'S'S' ") then
         phid="S'S'S'ac"
         icha = icha+1
         goto 100
      endif
      if(phid(1:8).eq."S'S'S'ac") then
         phid="S'S'S'df"
         icha = icha+1
         goto 100
      endif
      if(phid(1:5).eq.'SKKS ') then
         phid='SKKSac'
         icha = icha+1
         goto 100
      endif
      if(phid(1:6).eq.'SKKSac') then
         phid='SKKSdf'
         icha = icha+1
         goto 100
      endif
      if(phid(1:2).eq.'Sg') then
         phid='SKSdf'
         icha = icha+1
         goto 100
      endif

50    icha = 999
      return

100   if(flag) then
         phid0= s // phid(1:7)
      else
         phid0 = phid
      endif

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mult_ons(nobs,nobs0,iev,phase,tt,tts,tt2,azi,azis,
     +           p,ps,touse,amp,per,dtphmax,mread,arid)

      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      character*(*) touse(*), phase(*), arid(*)

      character*8 ph0*8, tu0*9

      real*8 dmerg1, dmerg2

      dimension iev(*),tt(*),tts(*),tt2(*),azi(*),azis(*),p(*),ps(*),
     +          amp(*),per(*)

      logical first

c
c     last changes 26 October 2016
c

      nobs0 = nobs

      do 9000 i = 1, nobs-1

         if (touse(i)(1:1) .eq. 'm') go to 8205

         ind = 1

         ph0 = phase(i)
         t0 = tt(i)
         ts0 = tts(i)
         t20 = tt2(i)
         tu0 = touse(i)
         p0 = p(i)
         ps0 = ps(i)
         indp = 0
         if(p0.gt.0.d0) indp = 1

         do 8200 j = i+1,nobs 

         if(iev(i).ne.iev(j)) go to 8200

         first = .false.


         if( (ph0(1:1).eq.'P'.and.phase(j).eq.'P1') .or.
     +       (ph0(1:1).eq.'S'.and.phase(j).eq.'S1') .or.
     +       (phase(j)(1:1).eq.'P'.and.ph0.eq.'P1') .or.
     +       (phase(j)(1:1).eq.'S'.and.ph0.eq.'S1')  ) first = .true.

         if( ph0.ne.phase(j) .and.  .not.first ) go to 8200

         dt0 = dabs(t0-tt(j))

         if(dt0.le.dtphmax) then 

               touse(j)(1:1) = 'm'
               if(touse(j)(4:4).eq.'D') touse(j)(4:4) = 'm'
               touse(j)(5:5) = ' '

               if(dt0.ge.0.01d0) then
                    
                  t0 = dmerg1(ind,t0,tt(j))

                  ts0 = dmerg2(tts(j),ts0,dt0)
                  t20 = dmerg2(tt2(j),t20,dt0)

                  touse(i)(1:1) = 'm'
                  if(touse(i)(4:4).eq.'D') touse(i)(4:4) = 'm'
                  touse(i)(5:5) = ' '

                  ind = ind + 1

               endif

               if(p(j).lt.0.d0 .or. touse(j)(3:3).eq.'m') go to 8200

               dp = dabs(p0-p(j))

               if(dp.le.0.0001d0) then 
                  if(ind.gt.1) touse(i)(3:3) = 'm'
                  go to 8150
               endif

               if(p0.lt.0.0001d0) then
                  p0 = p(j)
                  ps0 = ps(j)
                  indp = 1
                  if(ind.gt.1) then
                     touse(i)(3:3) = 'm'
                  else
                     p(i) = p0
                     ps(i) = ps(j)
                     touse(i)(3:3) = touse(j)(3:3)
                  endif
               else
                  p0 = dmerg1(indp,p0,p(j))
                  ps0 = dmerg2(ps0,ps(j),dp)
                  indp = indp + 1
                  touse(i)(3:3) = 'm'
               endif

8150           if(touse(j)(3:3).eq.'S') tu0(3:3) = 'S'
               touse(j)(3:3) = 'm'

               if (first .and. ph0(1:1).eq.'P') ph0= 'P1'
               if (first .and. ph0(1:1).eq.'S') ph0= 'S1'

         endif

8200     continue


         if (ind.gt.1) then
            nobs0 = nobs0 + 1

            if (nobs0 .le. mread) then
               iev(nobs0) = iev(i)
               phase(nobs0) = ph0
               tt(nobs0) = t0
               tts(nobs0) = ts0
               tt2(nobs0) = t20
               azi(nobs0) = -999.d0
               azis(nobs0) = -999.d0
               p(nobs0) = p0
               ps(nobs0) = ps0
               touse(nobs0) = tu0
               touse(nobs0)(2:2) = ' '
               touse(nobs0)(9:9) = '*'
               amp(nobs0) = -999.d0
               per(nobs0) =  -999.d0
               touse(nobs0)(6:6) = ' '
               arid(nobs0) = '  merged'
c              if(i.eq.1) then
c                print*,nobs, nobs0, iold
c                print*, iev(nobs0), phase(nobs0), tt(nobs0),
c    +           tts(nobs0), ind,indp
c              endif
            else
               print *, ' Merging of onsets extend maximum number '
     +              ,'of allowed onsets!'
            endif
         endif

8205     if (touse(i)(2:2).eq.'m' .or. azi(i).lt.0.d0) go to 9000

         do 8300 j = i+1,nobs

            if(iev(i).eq.iev(j)) then

               if(azi(j).lt.0.d0 .or. touse(j)(2:2).eq.'m') go to 8300

               if(dabs(azi(i)-azi(j)).le.0.1d0) then
                  if(touse(j)(2:2).eq.'A') touse(j)(2:2) = 'm'
               endif
            endif

8300     continue

9000  continue

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function dmerg1(n,v1,v2)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      real*8 dmerg1

      dmerg1 = ( v1 * dble(n) + v2 ) / dble(n + 1)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function dmerg2(v1,v2,v3)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)

      real*8 dmerg2, dpythag

      dmerg0 = dpythag(v1,v2)

      dmerg2 = dpythag(dmerg0,v3)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
