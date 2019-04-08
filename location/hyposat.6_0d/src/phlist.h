c
c     include file phlist.h
c
c     for common block phlist
c
c     if changing np, change also np in subroutine reflex(1)
c
      parameter (np=45)

      character phlist(np)*8

      common  /phasel/ phlist

