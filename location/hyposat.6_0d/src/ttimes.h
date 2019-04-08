c
c    include file ttimes.h
c

c
c     mphas  = max number of phases for tt-calculations (see in
c              tauget_mod and hyposat_loc )
c
c              The PARAMETER maxp in ttlim.h must have the same value!
c
      integer   mphas

      parameter (mphas = 120)

      character phcd(mphas)*8

      real*8 ttc(mphas),dtdd(mphas),dtdh(mphas),dddp(mphas)

