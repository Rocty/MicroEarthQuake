c
c    include file modelc.h
c

c     common block with information about the local/regional model
c     (used also for CRUST 1.0)
c
c     maxla = defined in model.h
c

      character filloc*120,azo(maxla)*4

      COMMON /MODELC/  filloc,azo

