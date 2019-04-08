CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function alpha1(a)
c  all angels are between -180 and + 180 degrees
      real*8 a,alpha1
      if(a.gt. 180.d0)  then
         alpha1 = a  - 360.d0
      else if(a.lt.-180.d0) then
         alpha1 = a  + 360.d0
      else
         alpha1 = a
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      function alpha2(a)
c  all angels are between 0 and 360 degrees
      real*8 a,alpha2
      if(a.gt. 360.d0)  then
         alpha2 = a  - 360.d0
      else if(a.lt.0.d0) then
         alpha2 = a  + 360.d0
      else
         alpha2 = a
      endif
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine findrange(amin,amax,a,n1,n2,ind)
c
c       fixed for longitudes: Nov 5, 1997 JS
c
      implicit real*8 (a-h,o-z)
      dimension a(*)
      integer ind,n1,n2

      amin = 9.d99
      amax = -9.d99

      if(ind.eq.1) then

       do 10 i=n1,n2,-1

         if(a(i).lt.amin) amin = a(i)
         if(a(i).gt.amax) amax = a(i)

10       continue

      else if (ind.eq.2) then
c
c       We have to handle the longitude values especially!
c

         deg2rad = datan(1.d0)/45.d0
       bmin = 9.d99
       bmax = -9.d99
       cmin = 9.d99
       cmax = -9.d99

       do 20 i=n1,n2,-1

         p1 = deg2rad*a(i)
         p2 = dcos(p1)
         p3 = dsin(p1)

         if(p2.lt.bmin) bmin = p2
         if(p2.gt.bmax) bmax = p2

         if(p3.lt.cmin) cmin = p3
         if(p3.gt.cmax) cmax = p3

20       continue

         amin = alpha1(datan2(cmin,bmax)/deg2rad)
         amax = alpha1(datan2(cmax,bmin)/deg2rad)

      endif
            
        return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       function bilinear to interpolate the value of a
c       point inside a given area of 4 cornerpoints and
c       the values at these cornerpoints respectively.
c 
c       Johannes Schweitzer, March 1999
c
c       assuming a cartesian grid is spanned by:
c
c              x(4),y(4)  .       x(3),y(3)
c                         .
c               ....... x1,y1 ...........
c                         .
c              x(1),y(1)  .       x(2),y(2)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function bilinear(x,y,z,x1,y1,indbi)
      implicit real*8 (a-h,o-z)

      dimension x(4), y(4), z(4)

      save

      if(indbi.ne.0) go to 10

      dx = x(2)-x(1)
      if(dx.ne.0.d0) then
        t = (x1-x(1))/dx
      else
        t = 1.d0
      endif

      dy = y(4)-y(1)
      if(dy.ne.0.d0) then
        u = (y1-y(1))/dy
      else
        u  = 1.d0
      endif


10      bilinear = (1.d0 - t)*(1.d0 - u) *z(1) +
     +                     t *(1.d0 - u) *z(2) +
     +                     t *        u  *z(3) +
     +             (1.d0 - t)*        u  *z(4)
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      function q2(a)
      real*8 a,q2
      q2 = a*a
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function dmean(a,n1,n2)
      real*8 a(*),dmean, b
      integer n1,n2
      if(n1.ne.n2) then
         b = 0.d0
         do 10 i = n1,n2,-1
10       b = b + a(i)
         dmean = b /dble(n1-n2+1)
      else
c dmean = 0.d0
       dmean =  a(n1)
      endif
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      function ddmax(a,n1,n2)
      real*8 a(*),ddmax
      integer n1,n2
      ddmax = -9.d99
      do 10 i = n1,n2,-1
      if(a(i).gt.ddmax) ddmax = a(i)
10    continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION UPPCAS(WORD)
      character*(*) word
      character*30 uppcas,tword
      integer lent

      lent = len(word) 
      if(lent.gt.30) then
         print *,' UPPCAS input too long ',word
         STOP
      endif
      uppcas = ' '
      tword = word(1:lent)
      do 100 i=1,lent
      ii= ichar(tword(i:i))
      if(ii.ge.97 .and. ii.le.122) then
         ii=ii-32
         tword(i:i)=char(ii)
      endif
100   continue
      uppcas(1:lent) = tword(1:lent)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION LOWCAS(WORD)
      character*(*) word
      character*30 lowcas,tword
      integer lent

      lent = len(word)
      if(lent.gt.30) then
         print *,' LOWCAS input too long ',word
         STOP
      endif
      lowcas = ' '
      tword = word(1:lent)
      do 100 i=1,lent
      ii= ichar(tword(i:i))
      if(ii.ge.65 .and. ii.le.90) then
         ii=ii+32
         tword(i:i)=char(ii)
      endif
100   continue
      lowcas(1:lent) = tword(1:lent)
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
