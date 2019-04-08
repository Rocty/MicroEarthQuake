c
c    include file model.h
c

c     common block with information about the local/regional model
c     (used also for CRUST 1.0)
c
c     maxla = maximum number of layers for a local/regional velocity
c             model
c

      parameter (maxla = 101)
      DIMENSION V0(2,maxla),z(maxla),rzv(2)
      character mtyp*3
      integer   jmod,iread,imo
      logical   locgeo,locsta

      COMMON /MODEL/  v0,z,rzv,elev,elatc,elonc,zmax,elat2,elon2,elev2,
     +                jmod,iread,imo,locgeo,locsta,mtyp

      real*8 v0,z,rzv,elev,elatc,elonc,zmax,elat2,elon2,elev2
