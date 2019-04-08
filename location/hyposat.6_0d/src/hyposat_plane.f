cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine plane(stala,stalo,iev,data,ndat,azi,dazi,
     +                 ray,dray,phas,touse,phase,jref,typctl)
C
C     author: Johannes Schweitzer, NORSAR
C             January 2002
C
      IMPLICIT real*8 (A-H,O-Z)

      real*8 dpythag
      dimension stala(*),stalo(*),data(*),iev(*)
      character*(*) touse(*),phase(*)
      character phas*8

      integer ndat,typctl

      include 'lsq.h'

c
      pi=4.0d0*datan(1.0d0)
      rad2phi=180.d0/pi
      phi2rad=pi/180.d0

      azi  = 0.d0
      dazi = 180.d0
      ray  = 0.d0
      dray = 100.d0

      j= 0
      d2k0 = 0.d0
      is   = 1
      jref = 0

51    DO 100 i=is,ndat

      if(touse(i)(8:8).ne.'P') go to 100

      if(j.eq.0) then
       phas = phase(i)
       stla = stala(iev(i))
       stlo = stalo(iev(i))
       dat0 = data(i)
       jref = i
       go to 90
      endif
      if(phase(i).ne.phas) go to 100
      call depi(stala(iev(i)),stalo(iev(i)),stla,stlo,dg,dk,
     +          azi0,baz,d2k)
      if(dk.le.0.d0) go to 100
      ggl(j,1) = dk*dsin(phi2rad*azi0)
      ggl(j,2) = dk*dcos(phi2rad*azi0)
      ddl(j) = data(i) - dat0
      d2k0 = d2k0 + d2k
90    j = j + 1
100   continue

      n=j-1

      if(n.lt.3) then

      j  = 0
      is = jref + 1
      jref = 0

      if (is.ge.ndat-2 .or. n.lt.0) then
         print *,'Tried, but cannot find enough data for a ',
     +          'plane-wave approximation'
         go to 900
      endif
      go to 51

      endif

      m=2
      d2k0 = d2k0 / dble(n)
c
      if(typctl .ge. 8) then
      write(*,
     +      '(//''The kernel for the plane wave inversion'')')
       do 150 i=1,n
           write(*,'(i3,2f10.1)') i,ggl(i,1),ggl(i,2)
150      continue
      endif
c
c     the least-squares fit
c
      call dlsq(n,m)
      if(lsqerr.gt.0) then
         print *,'ERROR for plane wave approximation'
         jref = 0
         go to 900
      endif
c

      phi = datan2(aal(1),aal(2))
      sph = dsin(phi)
      cph = dcos(phi)

      v  = sph / aal(1)

      a1 = sph*v
      a2 = cph*cph/aal(2)

      dazi  = rad2phi*dpythag(a2*vvl(1),a1*vvl(2))

      dray  = d2k0*dpythag(sph*vvl(1),cph*vvl(2))

      phi = phi + pi
      azi = alpha2(phi*rad2phi)

      ray = d2k0 / v

c
      if(typctl.gt.5 ) then

         write (*,'('' PLANE: Plane wave fit of '',i3,
     +       '' measured onsets:'')') n+1

         write (*,'(''ray parameter'',f7.2,
     +      ''+/-'',f7.2,''sec/deg azimuth:'',f7.2,''+/-'',f7.2)')
     +      ray,dray,azi,dazi

      end if

      if(typctl.gt.8 ) then
      write(*,
     +         '(//''Model parameter and standard deviations'')')
         do 200 i=1,m
              write(*,'(i3,3f8.4)') i,aal(i),vvl(i)
200        continue

        write(*,'(//''Observed data, calculated data, and residuals'')')
        do 250 i=1,n
           dt=ddl(i)-rrl(i)
           write(*,'(i3,3f10.4)') i,ddl(i),dt,rrl(i)
250     continue
      end if
c
900   continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
