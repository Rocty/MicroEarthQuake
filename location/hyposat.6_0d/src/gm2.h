c
c  include file gm2.h
c
c     to define the common block between the GMI routines
c
      parameter (nn = 5000, mm = 4)
      common  /gm2/ b(nn,nn),w(mm),v(mm,mm)
