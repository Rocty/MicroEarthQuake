c
c    include file modelg.h
c

c     common block with information about the local/regional model
c     (used also for CRUST 1.0)
c
c     maxla = defined in model.h
c

      dimension v0g(2,maxla),zg(maxla)

      integer jmodg

c     character azog(maxla)*4

c     COMMON /MODELG/  v0g,zg,elevg,zmaxg,jmodg,azog
      COMMON /MODELG/  v0g,zg,elevg,zmaxg,jmodg

