c
c     include file ref.h
c
c     for common block ref
c

      DIMENSION RR(2),TT(2),V(2,maxla),G(2,maxla),V2(2,maxla),
     *          ttp(np,3,3),ppp(np,3,3),ion(np,3,3),
     *          PAD(2),VHQ(2),FA(2,maxla),pa(2),del(3),
     *          ndisc(maxla)

      character phase*8

      logical   conv

      COMMON    /REF/ rr,tt,pa,fa,ttp,ppp,PIM,aa,del,
     *                VHQ,PAD,PI,V,G,V2,ndisc,
     *                ion,IB,IQQ,IQL,jh1,phase,conv

