CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE indexx(n,arr,indx)
C  (C) Copr. 1986-92 Numerical Recipes Software
      INTEGER n,indx(n),M,NSTACK
      REAL*4 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL*4 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
c j.s.  if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(jstack.gt.NSTACK) stop  'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SVDCMP(M,N,INDERR)
C  (C) Copr. 1986-92 Numerical Recipes Software
      IMPLICIT real*8 (A-H,O-Z)

      include 'gm2.h'

      DIMENSION RV1(NN)
      G=0.0d0
      SCALE=0.0d0
      ANORM=0.0d0
      INDERR = 0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0d0
        S=0.0d0
        SCALE=0.0d0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+DABS(B(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0d0) THEN
            DO 12 K=I,M
              B(K,I)=B(K,I)/SCALE
              S=S+B(K,I)*B(K,I)
12          CONTINUE
            F=B(I,I)
            G=-DSIGN(DSQRT(S),F)
            H=F*G-S
            B(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0d0
                DO 13 K=I,M
                  S=S+B(K,I)*B(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  B(K,J)=B(K,J)+F*B(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              B(K,I)=SCALE*B(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0d0
        S=0.0d0
        SCALE=0.0d0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+DABS(B(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0d0) THEN
            DO 18 K=L,N
              B(I,K)=B(I,K)/SCALE
              S=S+B(I,K)*B(I,K)
18          CONTINUE
            F=B(I,L)
            G=-DSIGN(DSQRT(S),F)
            H=F*G-S
            B(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=B(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0d0
                DO 21 K=L,N
                  S=S+B(J,K)*B(I,K)
21              CONTINUE
                DO 22 K=L,N
                  B(J,K)=B(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              B(I,K)=SCALE*B(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=DMAX1(ANORM,(DABS(W(I))+DABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0d0) THEN
            DO 26 J=L,N
              V(J,I)=(B(I,J)/B(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0d0
              DO 27 K=L,N
                S=S+B(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0d0
            V(J,I)=0.0d0
31        CONTINUE
        ENDIF
        V(I,I)=1.0d0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            B(I,J)=0.0d0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0d0) THEN
          G=1.0d0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0d0
              DO 34 K=L,M
                S=S+B(K,I)*B(K,J)
34            CONTINUE
              F=(S/B(I,I))*G
              DO 35 K=I,M
                B(K,J)=B(K,J)+F*B(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            B(J,I)=B(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            B(J,I)=0.0d0
38        CONTINUE
        ENDIF
        B(I,I)=B(I,I)+1.0d0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,60
          DO 41 L=K,1,-1
            NM=L-1
            IF ((DABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((DABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0d0
          S=1.0d0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((DABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=dpythag(f,g)
              W(I)=H
              H=1.0d0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=B(J,NM)
                Z=B(J,I)
                B(J,NM)=(Y*C)+(Z*S)
                B(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0d0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.100) THEN
             INDERR=100
             RETURN
          ENDIF
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0d0*H*Y)
          G=dpythag(f,1.0d0)
          F=((X-Z)*(X+Z)+H*((Y/(F+DSIGN(G,F)))-H))/X
          C=1.0d0
          S=1.0d0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=dpythag(F,H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=dpythag(F,H)
            W(J)=Z
            IF (Z.NE.0.0d0) THEN
              Z=1.0d0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=B(NM,J)
              Z=B(NM,I)
              B(NM,J)= (Y*C)+(Z*S)
              B(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0d0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION dpythag(a,b)
C  (C) Copr. 1986-92 Numerical Recipes Software
      real*8 a,b,dpythag
      real*8 absa,absb,q2
      absa=dabs(a)
      absb=dabs(b)
      if(absa.gt.absb)then
        dpythag=absa*dsqrt(1.0d0+q2(absb/absa))
      else
        if(absb.eq.0.0d0)then
          dpythag=0.0d0
        else
          dpythag=absb*dsqrt(1.0d0+q2(absa/absb))
        endif
      endif
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function getchi(dfree0,confi0)
c
c     uses subroutines and functions of the Numerical Recipes 
c     Software to calculate the chi**2. values for a specific
c     given confidence level and degree of freedom.
c

      real*8 dfree,dfree0,confi,getchi,confi0

      real*8 rtbis,xacc,xa,xe

      confi = 1.0d0 - (confi0 / 100.d0)

      dfree = dfree0 / 2.0d0

      call zbrak(xa,xe,dfree,confi)

      xacc=0.5d-6*(xa+xe)

      getchi = 2.d0 * rtbis(xa,xe,xacc,dfree,confi)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE zbrak(xa,xe,dfree,confi)
C  (C) Copr. 1986-92 Numerical Recipes Software
c
c     here adopted for calculating chi**2 for a specific
c     confidence level.
c
c     Johannes Schweitzer, NORSAR, October 2000
c
      real*8 confi,xa,xe,dfree,gammq
      INTEGER i,n
      real*8 dx,fc,fp,x,x1,x2
      n = 50
      x1= .5d0
      x2= 25.0d0
      x = x1
      dx=(x2-x1)/n
      fp= confi - gammq(dfree,x)

      do 11 i=1,n
        x=x+dx
        fc= confi - gammq(dfree,x)
        if(fc*fp.lt.0.d0) then
          xa=x-dx
          xe=x
        go to 1
        endif
        fp=fc
11    continue
1     continue
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION rtbis(x1,x2,xacc,dfree,confi)
C  (C) Copr. 1986-92 Numerical Recipes Software
c
c     Here adopted to get chi**2 for a specific confidence level
c     Johannes Schweitzer, October 2000, NORSAR
c
      real*8 rtbis,x1,x2,xacc,gammq,dfree,confi
      real*8 dx,f,fmid,xmid
      fmid= confi - gammq(dfree,x2)
      f=    confi - gammq(dfree,x1)
      if(f.lt.0.d0)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,100
        dx=dx*.5d0
        xmid=rtbis+dx
        fmid= confi - gammq(dfree,xmid)
        if(fmid.le.0.d0)rtbis=xmid
        if(dabs(dx).lt.xacc .or. fmid.eq.0.d0) return
11    continue
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION gammq(a0,x0)
      real*8 a0,a,gammq,x0,x
C  (C) Copr. 1986-92 Numerical Recipes Software
CU    USES gcf,gser
      real*8 gammcf,gamser,gln
      a = a0
      x = x0
      if(x.lt.a+1.d0)then
        call gser(gamser,a,x,gln)
        gammq=1.d0-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION gammln(xx)
C  (C) Copr. 1986-92 Numerical Recipes Software
      real*8 gammln,xx
      real*8 ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*dlog(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+dlog(stp*ser/x)
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gser(gamser,a,x,gln)
C  (C) Copr. 1986-92 Numerical Recipes Software
      INTEGER n, itmax
      real*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
CU    USES gammln
      real*8 ap,del,sum,gammln

      gln=gammln(a)
      if(x.le.0.d0)then
        if(x.lt.0.d0) print *,'GAMMQ: x < 0 in gser'
        gamser=0.d0
        return
      endif
      ap=a
      sum=1.d0/a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.d0
        del=del*x/ap
        sum=sum+del
        if(dabs(del).lt.dabs(sum)*EPS)goto 1
11    continue
      print *,'GAMMQ: a too large, ITMAX too small in gser'
1     gamser=sum*dexp(-x+a*dlog(x)-gln)
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gcf(gammcf,a,x,gln)
C  (C) Copr. 1986-92 Numerical Recipes Software
      INTEGER ITMAX, I
      real*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
CU    USES gammln
      real*8 an,b,c,d,del,h,gammln

      gln=gammln(a)
      b=x+1.d0-a
      c=1.d0/FPMIN
      d=1.d0/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.d0
        d=an*d+b
        if(dabs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(dabs(c).lt.FPMIN)c=FPMIN
        d=1.d0/d
        del=d*c
        h=h*del
        if(dabs(del-1.d0).lt.EPS)goto 1
11    continue
      print *, 'GAMMQ: a too large, ITMAX too small in gcf'
1     gammcf=dexp(-x+a*dlog(x)-gln)*h
      return
      END
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
