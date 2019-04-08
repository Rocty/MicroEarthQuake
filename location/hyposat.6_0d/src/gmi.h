c
c    include file gmi.h
c
c     common block and variable definition for GMI and LSQ-fit

c
c     mvar = maximum number of varibales for lsq-fit or gmi
c
c
c     mread should be 2*mread, as defined in main routine
c
c     Warning!!!!
c
c     Whenever changing the following parameter settings: 
c
c     change also the parameters in gm2.h !!!!
c
      parameter  (mvar = 4 , mread2 = 5000)

      dimension a(mread2,mvar),var(mvar),
     +          dat(mread2),r(mvar),res(mread2),dats(mread2),rs(mvar),
     +          ax1(2),ax2(2), dinf(mread2)

      logical   iellip

      common  /gmi/ a,var,dat,r,res,dats,rs,ax1,ax2,dinf,
     +              dinfm,iellip


