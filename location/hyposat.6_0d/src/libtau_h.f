      subroutine tabin(in,modnam)
c
c     6. februar 2002 changes added from libtau.src (of 21 FEB 1996)
c     j.s.
c
c     October 2002
c     Phase names with diff changed to dif as recommanded by IASPEI
c     working group on phase names.
c
      save

      include 'ttlim.h'
      character*(*) modnam
cjs     character requ*80
c     logical log

      character*8 phcd,phdif(9)
      real*8 pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp,dn
c
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
c
      common/tabc/dn,us(2),pt(jout),tau(4,jout),xlim(2,jout),
     1 xbrn(jbrn,3),dbrn(jbrn,2),xn,pn,tn,hn,jndx(jbrn,2),idel(jbrn,3),
     2 mbr1,mbr2
c
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
c
      common/pcdc/phcd(jbrn)
c
      character  ic*120, file_check*120
 
cjs   data tauc,xc/jtsm*0d0,jxsm*0d0/

      do i = 1,jtsm
      tauc(i) = 0.d0
      enddo
      do i = 1,jxsm
      xc(i) = 0.d0
      enddo
c
      nin=in
      phdif(1)='P'
      phdif(2)='S'
      phdif(3)='pP'
      phdif(4)='sP'
      phdif(5)='pS'
      phdif(6)='sS'
      phdif(7)='PKPab'
      phdif(8)='pPKPab'
      phdif(9)='sPKPab'
c++

      nb=index(modnam,' ')-1
      if(nb.le.0) nb=len(modnam)
c
      ic = modnam(1:nb) // '.hed'
      ic = file_check(ic)

      call assign(nin,1,ic)
c++ 

      read(nin,'(a)')
      read(nin,'(a)')
      read(nin,*) jsrctab ,jsegtab ,jbrntab ,jouttab ,jtsmtab

      if(jsrctab.gt.jsrc .or. jsegtab.gt.jseg .or. jbrntab.gt.jbrn .or.
     1   jouttab.gt.jout .or. jtsmtab.gt.jtsm) then

         print *,'LIBTAU variable dimensions too small for this table,'
         print *,'!! check dimensions in INCLUDE files for tables i/o '
         stop
      endif

      read(nin,'(a)')
      read(nin,'(a)')
      read(nin,*) nasgr,nl,len2,xn,pn,tn,mt,nseg,nbrn,ku,km

      read(nin,'(a)')
      read(nin,'(a)')
      read(nin,'(7i8,3g15.5)') ((nafl(i,j),j=1,3),
     1  (indx(i,j),j=1,2),(kndx(i,j),j=1,2),
     2  (fcs(i,j),j=1,3),i=1,jsegtab)

      if(jsegtab.lt.jseg) then
        do i = jsegtab+1,jseg
           do j = 1,2
              nafl(i,j) = 0
              indx(i,j) = 0
              kndx(i,j) = 0
              fcs(i,j) = 0.
            enddo
            nafl(i,3) = 0
            fcs(i,3) = 0.
         enddo
      endif

      read(nin,'(a)')
      read(nin,'(a)')
      read(nin,'(4g25.17,2i8)') ((pm(i,j),j=1,2),(zm(i,j),j=1,2),
     1      (ndex(i,j),j=1,2), i=1,jsrctab)

      if(jsrctab.lt.jsrc) then
        do i = jsrctab+1,jsrc
           do j = 1,2
              pm(i,j) = 0.d0
              zm(i,j) = 0.d0
              ndex(i,j) = 0
            enddo
         enddo
      endif

      read(nin,'(a)')
      read(nin,'(a)')
      jtsm0tab = jtsmtab + 1
      read(nin,'(2g25.17)')((pu(i,j),j=1,2),i=1,jtsm0tab)

      if(jtsm0tab.lt.jtsm0) then
        do i = jtsm0tab+1,jtsm0
           do j = 1,2
              pu(i,j) = 0.d0
            enddo
         enddo
      endif

      read(nin,'(a)')
      read(nin,'(a)')
      jxsmtab = jbrntab
      read(nin,'(2g25.17)')((pux(i,j),j=1,2),i=1,jxsmtab)

      if(jxsmtab.lt.jxsm) then
        do i = jxsmtab+1,jxsm
           do j = 1,2
              pux(i,j) = 0.d0
            enddo
         enddo
      endif

      read(nin,'(a)')
      read(nin,'(a)') text
      do 777 i = 1,jbrntab
      read(nin,'(a8,x,2i8,4g25.17)') phcd(i),(jndx(i,j),j=1,2),
     1      (px(i,j),j=1,2),(xt(i,j),j=1,2)
777   continue

      if(jbrntab.lt.jbrn) then
        do i = jbrntab+1,jbrn
           phcd(i) = '        '
           do j = 1,2
              jndx(i,j) = 0
              px(i,j) = 0.d0
              xt(i,j) = 0.d0
            enddo
         enddo
      endif

      read(nin,'(a)')
      read(nin,'(a)')
      read(nin,'(2g25.17)') (pt(i),taut(i),i=1,jouttab)

      read(nin,'(a)')
      read(nin,'(a)')
      read(nin,'(5g25.17)') ((coef(j,i),j=1,5),i=1,jouttab)

      if(jouttab.lt.jout) then
        do i = jouttab+1,jout
           pt(i) = 0.d0
           taut(i) = 0.d0
           do j = 1,5
              coef(i,j) = 0.d0
            enddo
         enddo
      endif

c     read(nin) nasgr,nl,len2,xn,pn,tn,mt,nseg,nbrn,ku,km,fcs,nafl,
c    1 indx,kndx
c     read(nin) pm,zm,ndex
c     read(nin) pu,pux
c     read(nin) phcd,px,xt,jndx
c     read(nin) pt,taut
c     read(nin) coef

      close(nin)
c     call retrns(nin)
c
c++
      nb=index(modnam,' ')-1
      if(nb.le.0) nb=len(modnam)

      ic =  modnam(1:nb) // '.tbl'
      ic = file_check (ic) 

      call assign(nin,1,ic)
c+++

      do 11 nph=1,2
 11   pu(ku(nph)+1,nph)=pm(1,nph)

c
c     write(10,*)'nasgr nl len2',nasgr,nl,len2
c     write(10,*)'nseg nbrn mt ku km',nseg,nbrn,mt,ku,km
c     write(10,*)'xn pn tn',xn,pn,tn
c     write(10,200)(i,(ndex(i,j),pm(i,j),zm(i,j),j=1,2),i=1,mt(2))
c200  format(/(1x,i3,i7,2f12.6,i7,2f12.6))
c     write(10,201)(i,(pu(i,j),j=1,2),i=1,ku(2)+1)
c201  format(/(1x,i3,2f12.6))
c     write(10,201)(i,(pux(i,j),j=1,2),i=1,km(2))
c     write(10,202)(i,(nafl(i,j),j=1,3),(indx(i,j),j=1,2),(kndx(i,j),
c    1 j=1,2),(fcs(i,j),j=1,3),i=1,nseg)
c202  format(/(1x,i3,7i5,3f5.0))
c     cn=180./3.1415927
c     write(10,203)(i,(jndx(i,j),j=1,2),(px(i,j),j=1,2),(cn*xt(i,j),
c    1 j=1,2),phcd(i),i=1,nbrn)
c203  format(/(1x,i3,2i5,2f12.6,2f12.2,2x,a))
c     write(10,204)(i,pt(i),taut(i),(coef(j,i),j=1,5),i=1,jout)
c204  format(/(1x,i5,0p2f12.6,1p5d10.2))
c
      tn=1./tn
c     dn=3.1415927/(180.*pn*xn)
      dn=datan(1.d0)/(45.d0*dble(pn)*dble(xn))
      odep=-1.
      ki=0
      msrc(1)=0
      msrc(2)=0
      k=1
      do 3 i=1,nbrn
      jidx(i)=jndx(i,2)
      do 4 j=1,2
 4    dbrn(i,j)=-1.d0
 8    if(jndx(i,2).le.indx(k,2)) go to 7
      k=k+1
      go to 8
 7    if(nafl(k,2).gt.0) go to 9
      ind=nafl(k,1)
      l=0
      do 10 j=jndx(i,1),jndx(i,2)
      l=l+1
 10   tp(l,ind)=pt(j)
 9    if(nafl(k,1).gt.0.and.(phcd(i)(1:1).eq.'P'.or.
     1 phcd(i)(1:1).eq.'S')) go to 3
      do 5 j=1,9
      if(phcd(i).eq.phdif(j)) go to 6
 5    continue
      go to 3
 6    dbrn(i,1)=1.d0
      phdif(j)=' '
 3    continue
c     write(10,205)(i,phcd(i),(dbrn(i,j),j=1,2),jidx(i),i=1,nbrn)
c205  format(/(1x,i5,2x,a,2f8.2,i5))
c     write(10,206)(i,(tp(i,j),j=1,2),i=1,jbrnu)
c206  format(/(1x,i5,2f12.6))
      return
      end

cjs      subroutine depset(dep,usrc)
      subroutine depset(dep)
      save 
      include 'ttlim.h'
c     logical dop,dos,segmsk,prnt
      logical dop,dos
cjs      real*4 usrc(2)
c     real*8 pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
      real*8 us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp,dn
      common/tabc/dn,us(2),pt(jout),tau(4,jout),xlim(2,jout),
     1 xbrn(jbrn,3),dbrn(jbrn,2),xn,pn,tn,hn,jndx(jbrn,2),idel(jbrn,3),
     2 mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
c     common/prtflc/segmsk(jseg),prnt(2)
cjs   data segmsk,prnt/jseg*.true.,2*.false./

      if(amax1(dep,.011).ne.odep) go to 1
      dop=.false.
      dos=.false.
      do 2 i=1,nseg
c     if(.not.segmsk(i).or.iidx(i).gt.0) go to 2
      if(iidx(i).gt.0) go to 2
      if(iabs(nafl(i,1)).le.1) dop=.true.
      if(iabs(nafl(i,1)).ge.2) dos=.true.
 2    continue
      if(.not.dop.and..not.dos) return
      go to 3
c
 1    nph0=0
      int0(1)=0
      int0(2)=0
      mbr1=nbrn+1
      mbr2=0
      dop=.false.
      dos=.false.
      do 4 i=1,nseg
c     if(.not.segmsk(i)) go to 4
      if(iabs(nafl(i,1)).le.1) dop=.true.
      if(iabs(nafl(i,1)).ge.2) dos=.true.
 4    continue
      do 5 i=1,nseg
      if(nafl(i,2).gt.0.or.odep.lt.0.) go to 5
      ind=nafl(i,1)
      k=0
      do 15 j=indx(i,1),indx(i,2)
      k=k+1
 15   pt(j)=tp(k,ind)
 5    iidx(i)=-1
      do 6 i=1,nbrn
 6    jndx(i,2)=-1
      if(ki.le.0) go to 7
      do 8 i=1,ki
      j=kk(i)
 8    pt(j)=pk(i)
      ki=0
c   Sample the model at the source depth.
 7    odep=amax1(dep,.011)
      rdep=dep
      if(rdep.lt..011) rdep=0.
      zs=dble(amin1(alog(amax1(1.-rdep*xn,1e-30)),0.))
      hn=1./(pn*(1.-rdep*xn))
c     if(prnt(1).or.prnt(2)) write(10,100)dep
c100  format(/1x,'Depth =',f7.2/)
c
 3    if(nph0.gt.1) go to 12
      if(dop) call depcor(1)
      if(dos) call depcor(2)
      go to 14
 12   if(dos) call depcor(2)
      if(dop) call depcor(1)
c
c   Interpolate all tau branches.
c
 14   j=1
      do 9 i=1,nseg
c     if(.not.segmsk(i)) go to 9
      nph=iabs(nafl(i,1))
c     print *,'i iidx nph msrc nafl =',i,iidx(i),nph,msrc(nph),nafl(i,1)
      if(iidx(i).gt.0.or.(msrc(nph).le.0.and.nafl(i,1).gt.0)) go to 9
      iidx(i)=1
      if(nafl(i,2).le.0) int=nafl(i,1)
      if(nafl(i,2).gt.0.and.nafl(i,2).eq.iabs(nafl(i,1)))
     1  int=nafl(i,2)+2
      if(nafl(i,2).gt.0.and.nafl(i,2).ne.iabs(nafl(i,1)))
     1  int=iabs(nafl(i,1))+4
      if(nafl(i,2).gt.0.and.nafl(i,2).ne.nafl(i,3)) int=nafl(i,2)+6
 11   if(jndx(j,1).ge.indx(i,1)) go to 10
      j=j+1
      go to 11
 10   idel(j,3)=nafl(i,1)
c     print *,'spfit:  j int =',j,int 
      call spfit(j,int)
      mbr1=min0(mbr1,j)
      mbr2=max0(mbr2,j)
      if(j.ge.nbrn) go to 9
      j=j+1
c     print *,'j jidx indx jndx =',j,jidx(j),indx(i,2),jndx(j,2)
      if(jidx(j).le.indx(i,2).and.jndx(j,2).gt.0) go to 10
 9    continue
c     write(10,*)'mbr1 mbr2',mbr1,mbr2
c     write(10,*)'msrc isrc odep zs us',msrc,isrc,odep,sngl(zs),
c    1 sngl(us(1)),sngl(us(2))
c     write(10,200)ki,(i,iidx(i),kk(i),pk(i),i=1,nseg)
c200  format(/10x,i5/(1x,3i5,f12.6))
cjs      usrc(1)=sngl(us(1))/pn
cjs      usrc(2)=sngl(us(2))/pn
      return
      end
      subroutine depcor(nph)
      save
      include 'ttlim.h'
      character*8 phcd
c     logical noend,noext,segmsk,prnt
      logical noend,noext
      real*8 pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp,ua,taua
      real*8 tup(jrec),umod,zmod,tauus1(2),tauus2(2),xus1(2),
     1 xus2(2),ttau,tx,sgn,umin,dtol,u0,u1,z0,z1,fac,du,dn
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
      common/tabc/dn,us(2),pt(jout),tau(4,jout),xlim(2,jout),
     1 xbrn(jbrn,3),dbrn(jbrn,2),xn,pn,tn,hn,jndx(jbrn,2),idel(jbrn,3),
     2 mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
      common/pdec/ua(5,2),taua(5,2),deplim,ka
c     common/prtflc/segmsk(jseg),prnt(2)
c     equivalence (tauc,tup)
cjs   data tol,dtol,deplim,ka,lpower/.01,1.d-6,1.1,4,7/
      data tol,dtol,lpower/.01,1.d-6,7/

c     deplim = 2.0
      deplim = 1.1
      ka     = 4
c
c/      do 111 i =1,jtsm
c/      tauu(i,1) = 0.d0
c/111   tauu(i,2) = 0.d0
c/      do 112 i =1,jrec
c/112   tup(i) = 0.d0
c/      do 113 i=1,2
c/      tauus1 (i) = 0.d0
c/      tauus2 (i) = 0.d0
c/      xus1 (i) = 0.d0
c/      xus2 (i) = 0.d0
c/113   continue
c/      do 114 i = 1,jxsm
c/      xu(i,1) = 0.d0
c/      xu(i,2) = 0.d0
c/114   xc(i)   = 0.d0
c/      du = 0.d0
c/      fac = 0.d0
c
c     umin=umod(zs,isrc,nph) 
      umin=umod(zs,isrc,nph,xn) 
c     write(10,*)'depcor:  nph nph0',nph,nph0
      if(nph.eq.nph0) go to 1
      nph0=nph
c     us(nph)=umod(zs,isrc,nph)
      us(nph)=umod(zs,isrc,nph,xn)
c   If we are in a high slowness zone, find the slowness of the lid.
      umin=us(nph)
      ks=isrc(nph)
c     write(10,*)'ks us',ks,sngl(umin)
      do 2 i=1,ks
      if(pm(i,nph).gt.umin) go to 2
      umin=pm(i,nph)
 2    continue
c   Find where the source slowness falls in the ray parameter array.
      n1=ku(nph)+1
      do 3 i=2,n1
      if(pu(i,nph).gt.umin) go to 4
 3    continue
      k2=n1
      if(pu(n1,nph).eq.umin) go to 50
cj.s. call abort('Source slowness too large.')
      print *,'Source slowness too large.'
      Stop
 4    k2=i
c50   write(10,*)'k2 umin',k2,sngl(umin)
c
c   Read in the appropriate depth correction values.
c
 50   noext=.false.
      sgn=1.d0
      if(msrc(nph).eq.0) msrc(nph)=1
c   See if the source depth coincides with a model sample
      ztol=xn*tol/(1.-xn*odep)
      if(dabs(zs-zm(ks+1,nph)).gt.dble(ztol)) go to 5
      ks=ks+1
      go to 6
 5    if(dabs(zs-zm(ks,nph)).gt.dble(ztol)) go to 7
c   If so flag the fact and make sure that the right integrals are
c   available.
 6    noext=.true.
      if(msrc(nph).eq.ks) go to 8
      call bkin(nin,ndex(ks,nph),ku(nph)+km(nph),tup)
      go to 11
c   If it is necessary to interpolate, see if appropriate integrals
c   have already been read in.
 7    if(msrc(nph).ne.ks+1) go to 9
      ks=ks+1
      sgn=-1.d0
      go to 8
 9    if(msrc(nph).eq.ks) go to 8
c   If not, read in integrals for the model depth nearest the source
c   depth.
      if(dabs(zm(ks,nph)-zs).le.dabs(zm(ks+1,nph)-zs)) go to 10
      ks=ks+1
      sgn=-1.d0
 10   call bkin(nin,ndex(ks,nph),ku(nph)+km(nph),tup)
c   Move the depth correction values to a less temporary area.
 11   do 31 i=1,ku(nph)
 31   tauu(i,nph)=tup(i)
      k=ku(nph)
      do 12 i=1,km(nph)
      k=k+1
      xc(i)=tup(k)
 12   xu(i,nph)=tup(k)
c     write(10,*)'bkin',ks,sngl(sgn),sngl(tauu(1,nph)),sngl(xu(1,nph))
c
c   Fiddle pointers.
c
 8    msrc(nph)=ks
c     write(10,*)'msrc sgn',msrc(nph),sngl(sgn)
      noend=.false.
      if(dabs(umin-pu(k2-1,nph)).le.dtol*umin) k2=k2-1
      if(dabs(umin-pu(k2,nph)).le.dtol*umin) noend=.true.
      if(msrc(nph).le.1.and.noext) msrc(nph)=0
      k1=k2-1
      if(noend) k1=k2
c     write(10,*)'noend noext k2 k1',noend,noext,k2,k1
      if(noext) go to 14
c
c   Correct the integrals for the depth interval [zm(msrc),zs].
c
      ms=msrc(nph)
Cjs, 26.02.2010
C      if(sgn)15,16,16
C 16   u0=pm(ms,nph)
C      z0=zm(ms,nph)
C      u1=us(nph)
C      z1=zs
C      go to 17
C 15   u0=us(nph)
C      z0=zs
C      u1=pm(ms,nph)
C      z1=zm(ms,nph)
C 17   mu=1

      if(sgn.lt.0.D0) then
         u0=us(nph)
         z0=zs
         u1=pm(ms,nph)
         z1=zm(ms,nph)
      else
         u0=pm(ms,nph)
         z0=zm(ms,nph)
         u1=us(nph)
         z1=zs
      endif

      mu = 1
         
c     write(10,*)'u0 z0',sngl(u0),sngl(z0)
c     write(10,*)'u1 z1',sngl(u1),sngl(z1)
      do 18 k=1,k1
      call tauint(pu(k,nph),u0,u1,z0,z1,ttau,tx)
      tauc(k)=tauu(k,nph)+sgn*ttau
      if(dabs(pu(k,nph)-pux(mu,nph)).gt.dtol) go to 18
      xc(mu)=xu(mu,nph)+sgn*tx
c     write(10,*)'up x:  k mu',k,mu,sngl(xu(mu,nph)),sngl(xc(mu))
      mu=mu+1
 18   continue
      go to 39
c   If there is no correction, copy the depth corrections to working
c   storage.
 14   mu=1
      do 40 k=1,k1
      tauc(k)=tauu(k,nph)
      if(dabs(pu(k,nph)-pux(mu,nph)).gt.dtol) go to 40
      xc(mu)=xu(mu,nph)
c     write(10,*)'up x:  k mu',k,mu,sngl(xu(mu,nph)),sngl(xc(mu))
      mu=mu+1
 40   continue
c
c   Calculate integrals for the ray bottoming at the source depth.
c
 39   xus1(nph)=0.d0
      xus2(nph)=0.d0
      mu=mu-1
      if(dabs(umin-us(nph)).gt.dtol.and.dabs(umin-pux(mu,nph)).le.dtol)
     1  mu=mu-1
c   This loop may be skipped only for surface focus as range is not
c   available for all ray parameters.
      if(msrc(nph).le.0) go to 1
      is=isrc(nph)
      tauus2(nph)=0.d0
      if(dabs(pux(mu,nph)-umin).gt.dtol.or.dabs(us(nph)-umin).gt.dtol)
     1  go to 48
c   If we happen to be right at a discontinuity, range is available.
      tauus1(nph)=tauc(k1)
      xus1(nph)=xc(mu)
c     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)),
c    1 sngl(xus1(nph)),'  *'
      go to 33
c   Integrate from the surface to the source.
 48   tauus1(nph)=0.d0
      j=1
      if(is.lt.2) go to 42
      do 19 i=2,is
      call tauint(umin,pm(j,nph),pm(i,nph),zm(j,nph),zm(i,nph),ttau,tx)
      tauus1(nph)=tauus1(nph)+ttau
      xus1(nph)=xus1(nph)+tx
 19   j=i
c     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)),
c    1 sngl(xus1(nph))
 42   if(dabs(zm(is,nph)-zs).le.dtol) go to 33
c   Unless the source is right on a sample slowness, one more partial
c   integral is needed.
      call tauint(umin,pm(is,nph),us(nph),zm(is,nph),zs,ttau,tx)
      tauus1(nph)=tauus1(nph)+ttau
      xus1(nph)=xus1(nph)+tx
c     write(10,*)'is ks tauus1 xus1',is,ks,sngl(tauus1(nph)),
c    1 sngl(xus1(nph))
 33   if(pm(is+1,nph).lt.umin) go to 41
c   If we are in a high slowness zone, we will also need to integrate
c   down to the turning point of the shallowest down-going ray.
      u1=us(nph)
      z1=zs
      do 35 i=is+1,mt(nph)
      u0=u1
      z0=z1
      u1=pm(i,nph)
      z1=zm(i,nph)
      if(u1.lt.umin) go to 36
      call tauint(umin,u0,u1,z0,z1,ttau,tx)
      tauus2(nph)=tauus2(nph)+ttau
 35   xus2(nph)=xus2(nph)+tx
c36   write(10,*)'is ks tauus2 xus2',is,ks,sngl(tauus2(nph)),
c    1 sngl(xus2(nph))
 36   z1=zmod(umin,i-1,nph)
      if(dabs(z0-z1).le.dtol) go to 41
c   Unless the turning point is right on a sample slowness, one more
c   partial integral is needed.
      call tauint(umin,u0,umin,z0,z1,ttau,tx)
      tauus2(nph)=tauus2(nph)+ttau
      xus2(nph)=xus2(nph)+tx
c     write(10,*)'is ks tauus2 xus2',is,ks,sngl(tauus2(nph)),
c    1 sngl(xus2(nph))
c
c   Take care of converted phases.
c
 41   iph=mod(nph,2)+1
      xus1(iph)=0.d0
      xus2(iph)=0.d0
      tauus1(iph)=0.d0
      tauus2(iph)=0.d0
      go to (59,61),nph
 61   if(umin.gt.pu(ku(1)+1,1)) go to 53
c
c   If we are doing an S-wave depth correction, we may need range and
c   tau for the P-wave which turns at the S-wave source slowness.  This
c   would bd needed for sPg and SPg when the source is in the deep mantle.
c
      do 44 j=1,nbrn
      if((phcd(j)(1:2).ne.'sP'.and.phcd(j)(1:2).ne.'SP').or.
     1 px(j,2).le.0.d0) go to 44
c     write(10,*)'Depcor:  j phcd px umin =',j,' ',phcd(j),px(j,1),
c    1 px(j,2),umin
      if(umin.ge.px(j,1).and.umin.lt.px(j,2)) go to 45
 44   continue
      go to 53
c
c   If we are doing an P-wave depth correction, we may need range and
c   tau for the S-wave which turns at the P-wave source slowness.  This
c   would be needed for pS and PS.
c
 59   do 60 j=1,nbrn
      if((phcd(j)(1:2).ne.'pS'.and.phcd(j)(1:2).ne.'PS').or.
     1 px(j,2).le.0.d0) go to 60
c     write(10,*)'Depcor:  j phcd px umin =',j,' ',phcd(j),px(j,1),
c    1 px(j,2),umin
      if(umin.ge.px(j,1).and.umin.lt.px(j,2)) go to 45
 60   continue
      go to 53
c
c   Do the integral.
 45   j=1
c     write(10,*)'Depcor:  do pS or sP integral - iph =',iph
      do 46 i=2,mt(iph)
      if(umin.ge.pm(i,iph)) go to 47
      call tauint(umin,pm(j,iph),pm(i,iph),zm(j,iph),zm(i,iph),ttau,tx)
      tauus1(iph)=tauus1(iph)+ttau
      xus1(iph)=xus1(iph)+tx
 46   j=i
 47   z1=zmod(umin,j,iph)
      if(dabs(zm(j,iph)-z1).le.dtol) go to 53
c   Unless the turning point is right on a sample slowness, one more
c   partial integral is needed.
      call tauint(umin,pm(j,iph),umin,zm(j,iph),z1,ttau,tx)
      tauus1(iph)=tauus1(iph)+ttau
      xus1(iph)=xus1(iph)+tx
c     write(10,*)'is ks tauusp xusp',j,ks,sngl(tauus1(iph)),
c    1 sngl(xus1(iph))
c
 53   ua(1,nph)=-1.d0
c     if(odep.ge.deplim.or.odep.le..1) go to 43
      if(odep.ge.deplim) go to 43
      do 57 i=1,nseg
c     if(.not.segmsk(i)) go to 57
      if(nafl(i,1).eq.nph.and.nafl(i,2).eq.0.and.iidx(i).le.0) go to 58
 57   continue
      go to 43
c
c   If the source is very shallow, we will need to insert some extra
c   ray parameter samples into the up-going branches.
c
 58   du=dble(amin1(1.e-5+(odep-.4)*2.e-5,1.e-5))
c     write(10,*)'Add:  nph is ka odep du us =',nph,is,ka,odep,
c    1 sngl(du),sngl(us(nph))
      lp=lpower
      k=0
      do 56 l=ka,1,-1
      k=k+1
      ua(k,nph)=us(nph)-(l**lp)*du
      lp=lp-1
      taua(k,nph)=0.d0
      j=1
      if(is.lt.2) go to 54
      do 55 i=2,is
      call tauint(ua(k,nph),pm(j,nph),pm(i,nph),zm(j,nph),zm(i,nph),
     1 ttau,tx)
      taua(k,nph)=taua(k,nph)+ttau
 55   j=i
c     write(10,*)'l k ua taua',l,k,sngl(ua(k,nph)),sngl(taua(k,nph))
 54   if(dabs(zm(is,nph)-zs).le.dtol) go to 56
c   Unless the source is right on a sample slowness, one more partial
c   integral is needed.
      call tauint(ua(k,nph),pm(is,nph),us(nph),zm(is,nph),zs,ttau,tx)
      taua(k,nph)=taua(k,nph)+ttau
c     write(10,*)'l k ua taua',l,k,sngl(ua(k,nph)),sngl(taua(k,nph))
 56   continue
      go to 43
c
c   Construct tau for all branches.
c
 1    mu=mu+1
 43   j=1
c     write(10,*)'mu',mu
c     write(10,*)'killer loop:'
      do 20 i=1,nseg
c     if(.not.segmsk(i)) go to 20
c     write(10,*)'i iidx nafl nph',i,iidx(i),nafl(i,1),nph
      if(iidx(i).gt.0.or.iabs(nafl(i,1)).ne.nph.or.(msrc(nph).le.0.and.
     1 nafl(i,1).gt.0)) go to 20
c
      iph=nafl(i,2)
      kph=nafl(i,3)
c   Handle up-going P and S.
      if(iph.le.0) iph=nph
      if(kph.le.0) kph=nph
      sgn=isign(1,nafl(i,1))
      i1=indx(i,1)
      i2=indx(i,2)
c     write(10,*)'i1 i2 sgn iph',i1,i2,sngl(sgn),iph
      m=1
      do 21 k=i1,i2
      if(pt(k).gt.umin) go to 22
23    if(dabs(pt(k)-pu(m,nph)).le.dtol) go to 21
      m=m+1
      go to 23
 21   tau(1,k)=taut(k)+sgn*tauc(m)
      k=i2
c     write(10,*)'k m',k,m
      go to 24
c22   write(10,*)'k m',k,m
 22   if(dabs(pt(k-1)-umin).le.dtol) k=k-1
      ki=ki+1
      kk(ki)=k
      pk(ki)=pt(k)
      pt(k)=umin
      fac=dble(fcs(i,1))
c     write(10,*)'ki fac',ki,sngl(fac)
      tau(1,k)=fac*(tauus1(iph)+tauus2(iph)+tauus1(kph)+tauus2(kph))+
     1 sgn*tauus1(nph)
c     write(10,*)'&&&&& nph iph kph tauus1 tauus2 tau =',
c    1 nph,iph,kph,sngl(tauus1(1)),sngl(tauus1(2)),sngl(tauus2(1)),
c    2 sngl(tauus2(2)),sngl(tau(1,k))
 24   m=1
 26   if(jndx(j,1).ge.indx(i,1)) go to 25
      j=j+1
      go to 26
 25   jndx(j,2)=min0(jidx(j),k)
      if(jndx(j,1).lt.jndx(j,2)) go to 37
      jndx(j,2)=-1
      go to 20
c37   write(10,*)'j jndx jidx',j,jndx(j,1),jndx(j,2),jidx(j),' ',
c    1 phcd(j)
 37   do 30 l=1,2
 28   if(dabs(pux(m,nph)-px(j,l)).le.dtol) go to 27
      if(m.ge.mu) go to 29
      m=m+1
      go to 28
 27   xbrn(j,l)=xt(j,l)+sgn*xc(m)
c     write(10,*)'x up:  j l m  ',j,l,m
      go to 30
 29   xbrn(j,l)=fac*(xus1(iph)+xus2(iph)+xus1(kph)+xus2(kph))+
     1 sgn*xus1(nph)
c     write(10,*)'x up:  j l end',j,l
c     write(10,*)'&&&&& nph iph kph xus1 xus2 xbrn =',
c    1 nph,iph,kph,sngl(xus1(1)),sngl(xus1(2)),sngl(xus2(1)),
c    2 sngl(xus2(2)),sngl(xbrn(j,l))
 30   continue
      if(j.ge.nbrn) go to 20
      j=j+1
      if(jndx(j,1).le.k) go to 25
 20   continue
      return
      end

cj.s. real*8 function umod(zs,isrc,nph)
      function umod(zs,isrc,nph,xn)
      real*8 umod
      save 
      include 'ttlim.h'
      character*31 msg
cjs   real*8 pm,zm,us,pt,tau,xlim,xbrn,dbrn
      real*8 pm,zm
      real*8 zs,uend,dtol,zmod
      dimension isrc(2)
      common/umdc/pm(jsrc,2),zm(jsrc,2),ndex(jsrc,2),mt(2)
cjs      common/tabc/us(2),pt(jout),tau(4,jout),xlim(2,jout),xbrn(jbrn,3),
cjs     1 dbrn(jbrn,2),xn,pn,tn,dn,hn,jndx(jbrn,2),idel(jbrn,3),mbr1,mbr2
      data dtol/1.d-6/
c
      m1=mt(nph)
      do 1 i=2,m1
      if(zm(i,nph).le.zs) go to 2
 1    continue
      dep=sngl((1.d0-dexp(zs))/dble(xn))
      write(msg,100)dep
 100  format('Source depth (',f6.1,') too deep.')
c     write(*,100)dep
cj.s. call abort(msg)
      print *,msg
      stop
 2    if(dabs(zs-zm(i,nph)).le.dtol.and.dabs(zm(i,nph)-zm(i+1,nph)).le.
     1 dtol) go to 3
      j=i-1
      isrc(nph)=j
      umod=pm(j,nph)+(pm(i,nph)-pm(j,nph))*(dexp(zs-zm(j,nph))-1.d0)/
     1 (dexp(zm(i,nph)-zm(j,nph))-1.d0)
      return
 3    isrc(nph)=i
      umod=pm(i+1,nph)
      return
c
      entry zmod(uend,js,nph)
      i=js+1
      zmod=zm(js,nph)+dlog(dmax1((uend-pm(js,nph))*(dexp(zm(i,nph)-
     1 zm(js,nph))-1.d0)/(pm(i,nph)-pm(js,nph))+1.d0,1.d-30))
      return
      end
c
      subroutine bkin(nin,nrec,len,buf)

c
c $$$$$ calls no other routines $$$$$
c
c   Bkin reads a block of len real*8 words into array buf(len)
c   from record nrec of the direct access unformatted file connected to
c   logical unit lu.
c
c   Here version changed to read input file with ASCII formatted data 
c   Johannes Schweizter, NORSAR, March 2017
c
      save

      include 'ttlim.h'

      integer nin,nrec,len,j2

      parameter (j2 = jtsm+jxsm)
      real*8 buf(j2),tmp(jsrc,j2)

      real*4 zso
      character modnamo*20
      logical first
      common /bkin0/first,modnamo,zso

c
      if(first) then
         first = .false.

         read(nin,*) jsrctab,j2tab

         if(jsrctab.gt.jsrc .or. j2tab.gt.j2) then

            print *,'LIBTAU variable dimensions too small,'
            print *,'!! check dimensions in INCLUDE files '
            stop
         endif

         do i = 1,jsrctab

            read(nin,'(2i8)',end=2) nrecr,lenr
            read(nin,'(5g25.17)') (tmp(nrecr,i0),i0=1,lenr)

            if(lenr.lt.j2) then
               do i2 = lenr+1,j2
                  tmp(nrecr,i2) = 0.d0
               enddo
            endif
         enddo

2        nrecm = nrecr
         close(nin)

      endif

      if(nrec.le.0 .or. nrec.gt.nrecm) go to 5

      do i = 1,len
         buf(i) = tmp(nrec,i)
      enddo

      return

c   If the record doesn't exist, zero fill the buffer.
5     do i=1,len
         buf(i)=0.d0
      enddo

      return
      end

      subroutine tauint(ptk,ptj,pti,zj,zi,tau,x)
      save
c
c $$$$$ calls warn $$$$$
c
c   Tauint evaluates the intercept (tau) and distance (x) integrals  for
c   the spherical earth assuming that slowness is linear between radii
c   for which the model is known.  The partial integrals are performed
c   for ray slowness ptk between model radii with slownesses ptj and pti
c   with equivalent flat earth depths zj and zi respectively.  The partial
c   integrals are returned in tau and x.  Note that ptk, ptj, pti, zj, zi,
c   tau, and x are all real*8.
c
      character*71 msg
      real*8 ptk,ptj,pti,zj,zi,tau,x,x90
      real*8 xx,b,sqk,sqi,sqj,sqb
c
      x90 = 2.d0 * datan(1.d0)
      if(dabs(zj-zi).le.1.d-9) go to 13
      if(dabs(ptj-pti).gt.1.d-9) go to 10
      if(dabs(ptk-pti).le.1.d-9) go to 13
      b=dabs(zj-zi)
      sqj=dsqrt(dabs(ptj*ptj-ptk*ptk))
      tau=b*sqj
      x=b*ptk/sqj
      go to 4
 10   if(ptk.gt.1.d-9.or.pti.gt.1.d-9) go to 1
c   Handle the straight through ray.
      tau=ptj
c     x=1.5707963267948966d0
      x = x90
      go to 4
 1    b=ptj-(pti-ptj)/(dexp(zi-zj)-1.d0)
      if(ptk.gt.1.d-9) go to 2
      tau=-(pti-ptj+b*dlog(pti/ptj)-b*dlog(dmax1((ptj-b)*pti/
     1 ((pti-b)*ptj),1.d-30)))
      x=0.d0
      go to 4
 2    if(ptk.eq.pti) go to 3
      if(ptk.eq.ptj) go to 11
      sqk=ptk*ptk
      sqi=dsqrt(dabs(pti*pti-sqk))
      sqj=dsqrt(dabs(ptj*ptj-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(sqb.gt.1.d-30) go to 5
      xx=0.d0
      x=ptk*(dsqrt(dabs((pti+b)/(pti-b)))-dsqrt(dabs((ptj+b)/
     1 (ptj-b))))/b
      go to 6
 5    if(b*b.lt.sqk) go to 7
      xx=dlog(dmax1((ptj-b)*(sqb*sqi+b*pti-sqk)/((pti-b)*
     1 (sqb*sqj+b*ptj-sqk)),1.d-30))
      x=ptk*xx/sqb
      go to 6
 7    xx=dasin(dmax1(dmin1((b*pti-sqk)/(ptk*dabs(pti-b)),1.d0),-1.d0))-
     1 dasin(dmax1(dmin1((b*ptj-sqk)/(ptk*dabs(ptj-b)),1.d0),-1.d0))
      x=-ptk*xx/sqb
 6    tau=-(sqi-sqj+b*dlog((pti+sqi)/(ptj+sqj))-sqb*xx)
      go to 4
 3    sqk=pti*pti
      sqj=dsqrt(dabs(ptj*ptj-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(b*b.lt.sqk) go to 8
      xx=dlog(dmax1((ptj-b)*(b*pti-sqk)/((pti-b)*(sqb*sqj+b*ptj-sqk)),
     1 1.d-30))
      x=pti*xx/sqb
      go to 9
c8    xx=dsign(1.5707963267948966d0,b-pti)-dasin(dmax1(dmin1((b*ptj-
 8    xx=dsign(x90,b-pti)-dasin(dmax1(dmin1((b*ptj-
     1 sqk)/(pti*dabs(ptj-b)),1.d0),-1.d0))
      x=-pti*xx/sqb
 9    tau=-(b*dlog(pti/(ptj+sqj))-sqj-sqb*xx)
      go to 4
 11   sqk=ptj*ptj
      sqi=dsqrt(dabs(pti*pti-sqk))
      sqb=dsqrt(dabs(b*b-sqk))
      if(b*b.lt.sqk) go to 12
      xx=dlog(dmax1((ptj-b)*(sqb*sqi+b*pti-sqk)/((pti-b)*(b*ptj-sqk)),
     1 1.d-30))
      x=ptj*xx/sqb
      go to 14
 12   xx=dasin(dmax1(dmin1((b*pti-sqk)/(ptj*dabs(pti-b)),1.d0),-1.d0))-
     1 dsign(x90,b-ptj)
c    1 dsign(1.5707963267948966d0,b-ptj)
      x=-ptj*xx/sqb
 14   tau=-(b*dlog((pti+sqi)/ptj)+sqi-sqb*xx)
c
c   Handle various error conditions.
c
 4    if(x.ge.-1.d-10) go to 15
      write(msg,100)ptk,ptj,pti,tau,x
 100  format('Bad range: ',1p5d12.4)
      call warn(msg)
 15   if(tau.ge.-1.d-10) go to 16
      write(msg,101)ptk,ptj,pti,tau,x
 101  format('Bad tau: ',1p5d12.4)
      call warn(msg(1:69))
 16   return
c   Trap null integrals and handle them properly.
 13   tau=0.d0
      x=0.d0
      return
      end

      subroutine spfit(jb,int)
      save 
      include 'ttlim.h'
      character*3 disc
      character*8 phcd
c     logical newgrd,makgrd,segmsk,prnt
      logical newgrd,makgrd
c     logical log
c     real*8 pm,zm,us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
      real*8 us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,
     1 xu,px,xt,taut,coef,tauc,xc,tcoef,tp
c     real*8 pmn,dmn,dmx,hm,shm,thm,p0,p1,tau0,tau1,x0,x1,pe,
c    1 pe0,spe0,scpe0,pe1,spe1,scpe1,dpe,dtau,dbrnch,cn,x180,x360,dtol,
c    2 ptol,xmin,difpkp,xbot,dn
      real*8 pmn,dmn,dmx,hm,shm,thm,p0,p1,tau0,tau1,x0,x1,pe,
     1 pe0,spe0,scpe0,pe1,spe1,scpe1,dpe,dtau,dbrnch,x180,x360,dtol,
     2 ptol,xmin,difpkp,xbot,dn
      common/tabc/dn,us(2),pt(jout),tau(4,jout),xlim(2,jout),
     1 xbrn(jbrn,3),dbrn(jbrn,2),xn,pn,tn,hn,jndx(jbrn,2),idel(jbrn,3),
     2 mbr1,mbr2
      common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
     1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
     2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
     3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
     4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
     5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)
cjs   data dbrnch,cn,x180,x360,dtol,ptol/2.5307274d0,57.295779d0,
cjs  1 3.1415927d0,6.2831853d0,1.d-6,2.d-6/
cjs   data difpkp/3.1415926d0/
cjs   data dbrnch,dtol,ptol/2.5307274d0,1d-6,2d-6/
      data dtol,ptol/1d-6,2d-6/
c
      x180 = 4.d0*datan(1.0D0)
      x360 = 2.d0*x180
c     cn   = 180.d0 / x180
      difpkp = x180 - 1.d-7
      dbrnch = x180*145.d0/180.d0
c
      xmin=3.92403d-3
c     if(prnt(1)) write(10,102)
      i1=jndx(jb,1)
      i2=jndx(jb,2)
c     write(10,*)'Spfit:  jb i1 i2 pt =',jb,i1,i2,sngl(pt(i1)),
c    1 sngl(pt(i2))
      if(i2-i1.gt.1.or.dabs(pt(i2)-pt(i1)).gt.ptol) go to 14
      jndx(jb,2)=-1
      return
 14   newgrd=.false.
      makgrd=.false.
      if(dabs(px(jb,2)-pt(i2)).gt.dtol) newgrd=.true.
c     write(10,*)'Spfit:  px newgrd =',sngl(px(jb,2)),newgrd
      if(.not.newgrd) go to 10
      k=mod(int-1,2)+1
      if(int.ne.int0(k)) makgrd=.true.
c     write(10,*)'Spfit:  int k int0 makgrd =',int,k,int0(k),makgrd
      if(int.gt.2) go to 12
c     call query('Enter xmin:',log)
c     read *,xmin
c     xmin=xmin*xn
      xmin=dble(xn*amin1(amax1(2.*odep,2.),25.))
c     write(10,*)'Spfit:  xmin =',xmin,xmin/xn
c     call pdecu(i1,i2,xbrn(jb,1),xbrn(jb,2),xmin,int,i3)
      call pdecu(i1,i2,xbrn(jb,1),xbrn(jb,2),xmin,int,ijo)
      i2 = ijo
      jndx(jb,2)=i2
 12   nn=i2-i1+1
      if(makgrd) call tauspl(1,nn,pt(i1),tcoef(1,1,k))
c     write(10,301,iostat=ios)jb,k,nn,int,newgrd,makgrd,
c    1 xbrn(jb,1),xbrn(jb,2),(i,pt(i-1+i1),tau(1,i-1+i1),
c    2 (tcoef(j,i,k),j=1,5),i=1,nn)
c301  format(/1x,4i3,2l3,2f12.8/(1x,i5,0p2f12.8,1p5d10.2))
      call fitspl(1,nn,tau(1,i1),xbrn(jb,1),xbrn(jb,2),tcoef(1,1,k))
      int0(k)=int
      go to 11
 10   call fitspl(i1,i2,tau,xbrn(jb,1),xbrn(jb,2),coef)
 11   pmn=pt(i1)
      dmn=xbrn(jb,1)
      dmx=dmn
      mxcnt=0
      mncnt=0
c     call appx(i1,i2,xbrn(jb,1),xbrn(jb,2))
c     write(10,300)(i,pt(i),(tau(j,i),j=1,3),i=i1,i2)
c300  format(/(1x,i5,4f12.6))
      pe=pt(i2)
      p1=pt(i1)
      tau1=tau(1,i1)
      xbot=tau(2,i1)
      x1=xbot
      pe1=pe-p1
      spe1=dsqrt(dabs(pe1))
      scpe1=pe1*spe1
      j=i1
      is=i1+1
      do 2 i=is,i2
      p0=p1
      p1=pt(i)
      tau0=tau1
      tau1=tau(1,i)
      x0=x1
      x1=tau(2,i)
      dpe=p0-p1
      dtau=tau1-tau0
      pe0=pe1
      pe1=pe-p1
      spe0=spe1
      spe1=dsqrt(dabs(pe1))
      scpe0=scpe1
      scpe1=pe1*spe1
      tau(4,j)=(2.d0*dtau-dpe*(x1+x0))/(.5d0*(scpe1-scpe0)-1.5d0*spe1*
     1 spe0*(spe1-spe0))
      tau(3,j)=(dtau-dpe*x0-(scpe1+.5d0*scpe0-1.5d0*pe1*spe0)*tau(4,j))/
     1 (dpe*dpe)
      tau(2,j)=(dtau-(pe1*pe1-pe0*pe0)*tau(3,j)-(scpe1-scpe0)*tau(4,j))/
     1 dpe
      tau(1,j)=tau0-scpe0*tau(4,j)-pe0*(pe0*tau(3,j)+tau(2,j))
      xlim(1,j)=dmin1(x0,x1)
      xlim(2,j)=dmax1(x0,x1)
      if(xlim(1,j).ge.dmn) go to 5
      dmn=xlim(1,j)
      pmn=pt(j)
      if(x1.lt.x0) pmn=pt(i)
 5    disc=' '
      if(dabs(tau(3,j)).le.1.d-30) go to 4
      shm=-.375d0*tau(4,j)/tau(3,j)
      hm=shm*shm
      if(shm.le.0.d0.or.(hm.le.pe1.or.hm.ge.pe0)) go to 4
cjs.s 7    thm=tau(2,j)+shm*(2d0*shm*tau(3,j)+1.5d0*tau(4,j))
      thm=tau(2,j)+shm*(2.d0*shm*tau(3,j)+1.5d0*tau(4,j))
      xlim(1,j)=dmin1(xlim(1,j),thm)
      xlim(2,j)=dmax1(xlim(2,j),thm)
      if(thm.ge.dmn) go to 6
      dmn=thm
      pmn=pe-hm
 6    disc='max'
      if(tau(4,j).lt.0.d0) disc='min'
      if(disc.eq.'max') mxcnt=mxcnt+1
      if(disc.eq.'min') mncnt=mncnt+1
c4    if(prnt(1)) write(10,100,iostat=ios)disc,j,pt(j),
c    1 (tau(k,j),k=1,4),(cn*xlim(k,j),k=1,2)
 4    continue
c100  format(1x,a,i5,f10.6,1p4e10.2,0p2f7.2)
      dmx=dmax1(dmx,xlim(2,j))
 2    j=i
c     if(prnt(1)) write(10,100,iostat=ios)'   ',j,pt(j)
      xbrn(jb,1)=dmn
      xbrn(jb,2)=dmx
      xbrn(jb,3)=pmn
      idel(jb,1)=1
      idel(jb,2)=1
      if(xbrn(jb,1).gt.x180) idel(jb,1)=2
      if(xbrn(jb,2).gt.x180) idel(jb,2)=2
      if(xbrn(jb,1).gt.x360) idel(jb,1)=3
      if(xbrn(jb,2).gt.x360) idel(jb,2)=3
      if(int.gt.2) go to 1
      phcd(jb)=phcd(jb)(1:1)
      i=jb
      do 8 j=1,nbrn
      i=mod(i,nbrn)+1
      if(phcd(i)(1:1).eq.phcd(jb).and.phcd(i)(2:2).ne.'P'.and.
     1 (pe.ge.px(i,1).and.pe.le.px(i,2))) go to 9
 8    continue
      go to 1
 9    phcd(jb)=phcd(i)
      if(dabs(pt(i2)-pt(jndx(i,1))).le.dtol) phcd(jb)=phcd(i-1)
c1    if(prnt(1).and.prnt(2)) write(10,102)
 1    continue
c102  format()
      if(dbrn(jb,1).le.0.d0) go to 3
      dbrn(jb,1)=xbot
      dbrn(jb,2)=dbrnch
      if(index(phcd(jb),'ab').gt.0) dbrn(jb,2)=difpkp

c     if(prnt(2)) write(10,101,iostat=ios)phcd(jb),
c    1 (jndx(jb,k),k=1,2),(cn*xbrn(jb,k),k=1,2),xbrn(jb,3),
c    2 (cn*dbrn(jb,k),k=1,2),(idel(jb,k),k=1,3),int,newgrd,makgrd
c101  format(1x,a,2i5,2f8.2,f8.4,2f8.2,4i3,2l2)

      go to 15

c3    if(prnt(2)) write(10,103,iostat=ios)phcd(jb),
c    1 (jndx(jb,k),k=1,2),(cn*xbrn(jb,k),k=1,2),xbrn(jb,3),
c    2 (idel(jb,k),k=1,3),int,newgrd,makgrd
c103  format(1x,a,2i5,2f8.2,f8.4,16x,4i3,2l2)

 3    continue
 15   if(mxcnt.gt.mncnt.or.mncnt.gt.mxcnt+1)
     1 call warn('Bad interpolation on '//phcd(jb))
      return
      end
      subroutine pdecu(i1,i2,x0,x1,xmin,int,len)
      save 
      include 'ttlim.h'
      real*8 us,pt,tau,xlim,xbrn,dbrn,ua,taua,dn
      real*8 x0,x1,xmin,dx,dx2,sgn,rnd,xm,axm,x,h1,h2,hh,xs
      common/tabc/dn,us(2),pt(jout),tau(4,jout),xlim(2,jout),
     1 xbrn(jbrn,3),dbrn(jbrn,2),xn,pn,tn,hn,jndx(jbrn,2),idel(jbrn,3),
     2 mbr1,mbr2
      common/pdec/ua(5,2),taua(5,2),deplim,ka
c
c     write(10,*)'Pdecu:  ua =',sngl(ua(1,int))
      if(ua(1,int).le.0.d0) go to 17
c     write(10,*)'Pdecu:  fill in new grid'
      k=i1+1
      do 18 i=1,ka
      pt(k)=ua(i,int)
      tau(1,k)=taua(i,int)
 18   k=k+1
      pt(k)=pt(i2)
      tau(1,k)=tau(1,i2)
      go to 19
c
 17   is=i1+1
      ie=i2-1
      xs=x1
      do 11 i=ie,i1,-1
      x=xs
      if(i.ne.i1) go to 12
      xs=x0
      go to 14
 12   h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
      xs=-(h2*tau(1,i-1)-(h2+h1)*tau(1,i)+h1*tau(1,i+1))/hh
 14   if(dabs(x-xs).le.xmin) go to 15
 11   continue
      len=i2
      return
 15   ie=i
      if(dabs(x-xs).gt..75d0*xmin.or.ie.eq.i2) go to 16
      xs=x
      ie=ie+1
 16   n=max0(idint(dabs(xs-x0)/xmin+.8d0),1)
      dx=(xs-x0)/n
      dx2=dabs(.5d0*dx)
      sgn=dsign(1.d0,dx)
      rnd=0.d0
      if(sgn.gt.0.d0) rnd=1.d0
      xm=x0+dx
      k=i1
      m=is
      axm=1d10
      do 1 i=is,ie
      if(i.lt.ie) go to 8
      x=xs
      go to 5
 8    h1=pt(i-1)-pt(i)
      h2=pt(i+1)-pt(i)
      hh=h1*h2*(h1-h2)
      h1=h1*h1
      h2=-h2*h2
      x=-(h2*tau(1,i-1)-(h2+h1)*tau(1,i)+h1*tau(1,i+1))/hh
 5    if(sgn*(x-xm).le.dx2) go to 2
      if(k.lt.m) go to 3
      do 4 j=m,k
 4    pt(j)=-1.d0
 3    m=k+2
      k=i-1
      axm=1d10
cj.s.  7    xm=xm+dx*idint((x-xm-dx2)/dx+rnd)
      xm=xm+dx*idint((x-xm-dx2)/dx+rnd)
 2    if(dabs(x-xm).ge.axm) go to 1
      axm=dabs(x-xm)
      k=i-1
 1    continue
      if(k.lt.m) go to 9
      do 6 j=m,k
 6    pt(j)=-1.d0
 9    k=i1
      do 10 i=is,i2
      if(pt(i).lt.0.d0) go to 10
      k=k+1
      pt(k)=pt(i)
      tau(1,k)=tau(1,i)
 10   continue
 19   len=k
c     write(10,300)(i,pt(i),tau(1,i),i=i1,len)
c300  format(/(1x,i5,0pf12.6,1pd15.4))
      return
      end
      subroutine tauspl(i1,i2,pt,coef)
c
c $$$$$ calls only library routines $$$$$
c
c   Given ray parameter grid pt;i (pt sub i), i=i1,i1+1,...,i2, tauspl
c   determines the i2-i1+3 basis functions for interpolation I such
c   that:
c
c      tau(p) = a;1,i + Dp * a;2,i + Dp**2 * a;3,i + Dp**(3/2) * a;4,i
c
c   where Dp = pt;n - p, pt;i <= p < pt;i+1, and the a;j,i's are
c   interpolation coefficients.  Rather than returning the coefficients,
c   a;j,i, which necessarily depend on tau(pt;i), i=i1,i1+1,...,i2 and
c   x(pt;i) (= -d tau(p)/d p | pt;i), i=i1,i2, tauspl returns the
c   contribution of each basis function and its derivitive at each
c   sample.  Each basis function is non-zero at three grid points,
c   therefore, each grid point will have contributions (function values
c   and derivitives) from three basis functions.  Due to the basis
c   function normalization, one of the function values will always be
c   one and is not returned in array coef with the other values.
c   Rewritten on 23 December 1983 by R. Buland.
c
      save
      real*8 pt(i2),coef(5,i2)
      real*8 del(5),sdel(5),deli(5),d3h(4),d1h(4),dih(4),
     1 d(4),ali,alr,b3h,b1h,bih,th0p,th2p,th3p,th2m
c
      n2=i2-i1-1
      if(n2.le.-1) return
      is=i1+1
c
c   To achieve the requisite stability, proceed by constructing basis
c   functions G;i, i=0,1,...,n+1.  G;i will be non-zero only on the
c   interval [p;i-2,p;i+2] and will be continuous with continuous first
c   and second derivitives.  G;i(p;i-2) and G;i(p;i+2) are constrained
c   to be zero with zero first and second derivitives.  G;i(p;i) is
c   normalized to unity.
c
c   Set up temporary variables appropriate for G;-1.  Note that to get
c   started, the ray parameter grid is extrapolated to yeild p;i, i=-2,
c   -1,0,1,...,n.
      del(2)=pt(i2)-pt(i1)+3.d0*(pt(is)-pt(i1))
      sdel(2)=dsqrt(dabs(del(2)))
      deli(2)=1.d0/sdel(2)
      m=2
      do 1 k=3,5
      del(k)=pt(i2)-pt(i1)+(5-k)*(pt(is)-pt(i1))
      sdel(k)=dsqrt(dabs(del(k)))
      deli(k)=1.d0/sdel(k)
      d3h(m)=del(k)*sdel(k)-del(m)*sdel(m)
      d1h(m)=sdel(k)-sdel(m)
      dih(m)=deli(k)-deli(m)
 1    m=k
      l=i1-1
      if(n2.le.0) go to 10
c   Loop over G;i, i=0,1,...,n-3.
      do 2 i=1,n2
      m=1
c   Update temporary variables for G;i-1.
      do 3 k=2,5
      del(m)=del(k)
      sdel(m)=sdel(k)
      deli(m)=deli(k)
      if(k.ge.5) go to 3
      d3h(m)=d3h(k)
      d1h(m)=d1h(k)
      dih(m)=dih(k)
 3    m=k
      l=l+1
      del(5)=pt(i2)-pt(l+1)
      sdel(5)=dsqrt(dabs(del(5)))
      deli(5)=1.d0/sdel(5)
      d3h(4)=del(5)*sdel(5)-del(4)*sdel(4)
      d1h(4)=sdel(5)-sdel(4)
      dih(4)=deli(5)-deli(4)
c   Construct G;i-1.
      ali=1.d0/(.125d0*d3h(1)-(.75d0*d1h(1)+.375d0*dih(1)*del(3))*
     1 del(3))
      alr=ali*(.125d0*del(2)*sdel(2)-(.75d0*sdel(2)+.375d0*del(3)*
     1 deli(2)-sdel(3))*del(3))
      b3h=d3h(2)+alr*d3h(1)
      b1h=d1h(2)+alr*d1h(1)
      bih=dih(2)+alr*dih(1)
      th0p=d1h(1)*b3h-d3h(1)*b1h
      th2p=d1h(3)*b3h-d3h(3)*b1h
      th3p=d1h(4)*b3h-d3h(4)*b1h
      th2m=dih(3)*b3h-d3h(3)*bih
c   The d;i's completely define G;i-1.
      d(4)=ali*((dih(1)*b3h-d3h(1)*bih)*th2p-th2m*th0p)/((dih(4)*b3h-
     1 d3h(4)*bih)*th2p-th2m*th3p)
      d(3)=(th0p*ali-th3p*d(4))/th2p
      d(2)=(d3h(1)*ali-d3h(3)*d(3)-d3h(4)*d(4))/b3h
      d(1)=alr*d(2)-ali
c   Construct the contributions G;i-1(p;i-2) and G;i-1(p;i).
c   G;i-1(p;i-1) need not be constructed as it is normalized to unity.
      coef(1,l)=(.125d0*del(5)*sdel(5)-(.75d0*sdel(5)+.375d0*deli(5)*
     1 del(4)-sdel(4))*del(4))*d(4)
      if(i.ge.3) coef(2,l-2)=(.125d0*del(1)*sdel(1)-(.75d0*sdel(1)+
     1 .375d0*deli(1)*del(2)-sdel(2))*del(2))*d(1)
c   Construct the contributions -dG;i-1(p)/dp | p;i-2, p;i-1, and p;i.
      coef(3,l)=-.75d0*(sdel(5)+deli(5)*del(4)-2.d0*sdel(4))*d(4)
      if(i.ge.2) coef(4,l-1)=-.75d0*((sdel(2)+deli(2)*del(3)-
     1 2.d0*sdel(3))*d(2)-(d1h(1)+dih(1)*del(3))*d(1))
      if(i.ge.3) coef(5,l-2)=-.75d0*(sdel(1)+deli(1)*del(2)-
     1 2.d0*sdel(2))*d(1)
 2    continue
c   Loop over G;i, i=n-2,n-1,n,n+1.  These cases must be handled
c   seperately because of the singularities in the second derivitive
c   at p;n.
 10   do 4 j=1,4
      m=1
c   Update temporary variables for G;i-1.
      do 5 k=2,5
      del(m)=del(k)
      sdel(m)=sdel(k)
      deli(m)=deli(k)
      if(k.ge.5) go to 5
      d3h(m)=d3h(k)
      d1h(m)=d1h(k)
      dih(m)=dih(k)
 5    m=k
      l=l+1
      del(5)=0.d0
      sdel(5)=0.d0
      deli(5)=0.d0
c   Construction of the d;i's is different for each case.  In cases
c   G;i, i=n-1,n,n+1, G;i is truncated at p;n to avoid patching across
c   the singularity in the second derivitive.
      if(j.lt.4) go to 6

c   For G;n+1 constrain G;n+1(p;n) to be .25.
      d(1)=2.d0/(del(1)*sdel(1))
      go to 9

c   For G;i, i=n-2,n-1,n, the condition dG;i(p)/dp|p;i = 0 has been
c   substituted for the second derivitive continuity condition that
c   can no longer be satisfied.
 6    alr=(sdel(2)+deli(2)*del(3)-2.d0*sdel(3))/(d1h(1)+dih(1)*del(3))

      d(2)=1.d0/(.125d0*del(2)*sdel(2)-(.75d0*sdel(2)+.375d0*deli(2)*
     1 del(3)-sdel(3))*del(3)-(.125d0*d3h(1)-(.75d0*d1h(1)+.375d0*
     2 dih(1)*del(3))*del(3))*alr)

      d(1)=alr*d(2)

cjs 26.02.2010
c
c     if(j-2)8,7,9
c

      if((j-2).lt.0) then
c   No additional constraints are required for G;n-2.
         d(3)=-((d3h(2)-d1h(2)*del(4))*d(2)+(d3h(1)-d1h(1)*del(4))*
     1        d(1))/(d3h(3)-d1h(3)*del(4))
         d(4)=(d3h(3)*d(3)+d3h(2)*d(2)+d3h(1)*d(1))/(del(4)*sdel(4))
         goto 9
      else if((j-2).gt.0) then
         goto 9
      endif

c   For G;n-1 constrain G;n-1(p;n) to be .25.
cjs 7     d(3)=(2.d0+d3h(2)*d(2)+d3h(1)*d(1))/(del(3)*sdel(3))
      d(3)=(2.d0+d3h(2)*d(2)+d3h(1)*d(1))/(del(3)*sdel(3))

c   Construct the contributions G;i-1(p;i-2) and G;i-1(p;i).
 9    if(j.le.2) then
         coef(1,l)=(.125d0*del(3)*sdel(3)-(.75d0*sdel(3)+.375d0*
     1   deli(3)*del(4)-sdel(4))*del(4))*d(3)-(.125d0*d3h(2)-(.75d0*
     2   d1h(2)+.375d0*dih(2)*del(4))*del(4))*d(2)-(.125d0*d3h(1)-(
     3   .75d0*d1h(1)+.375d0*dih(1)*del(4))*del(4))*d(1)

c   Construct the contributions -dG;i-1(p)/dp | p;i-2, p;i-1, and p;i.
         coef(3,l)=-.75d0*((sdel(3)+deli(3)*del(4)-
     1   2.d0*sdel(4))*d(3)-(d1h(2)+dih(2)*del(4))*d(2)-(d1h(1)+
     2   dih(1)*del(4))*d(1))

      endif

c   Construct the contributions G;i-1(p;i-2) and G;i-1(p;i).
      if(l-i1.gt.1) then
         coef(2,l-2)=(.125d0*del(1)*sdel(1)-(.75d0*sdel(1)+
     1   .375d0*deli(1)*del(2)-sdel(2))*del(2))*d(1)

c   Construct the contributions -dG;i-1(p)/dp | p;i-2, p;i-1, and p;i.
         coef(5,l-2)=-.75d0*(sdel(1)+deli(1)*del(2)-
     1   2.d0*sdel(2))*d(1)
      endif

      if(j.le.3.and.l-i1.gt.0) coef(4,l-1)=0.d0

 4    continue

      return
      end
      subroutine fitspl(i1,i2,tau,x1,xn,coef)
c
c $$$$$ calls only library routines $$$$$
c
c   Given ray parameter grid p;i (p sub i), i=1,2,...,n, corresponding
c   tau;i values, and x;1 and x;n (x;i = -dtau/dp|p;i); tauspl finds
c   interpolation I such that:  tau(p) = a;1,i + Dp * a;2,i + Dp**2 *
c   a;3,i + Dp**(3/2) * a;4,i where Dp = p;n - p and p;i <= p < p;i+1.
c   Interpolation I has the following properties:  1) x;1, x;n, and
c   tau;i, i=1,2,...,n are fit exactly, 2) the first and second
c   derivitives with respect to p are continuous everywhere, and
c   3) because of the paramaterization d**2 tau/dp**2|p;n is infinite.
c   Thus, interpolation I models the asymptotic behavior of tau(p)
c   when tau(p;n) is a branch end due to a discontinuity in the
c   velocity model.  Note that array a must be dimensioned at least
c   a(4,n) though the interpolation coefficients will be returned in
c   the first n-1 columns.  The remaining column is used as scratch
c   space and returned as all zeros.  Programmed on 16 August 1982 by
c   R. Buland.
c
      save 
      real*8 tau(4,i2),x1,xn,coef(5,i2),a(2,i2),ap(3),
     1 b(i2),alr,gn
cj.s.     1 b(100),alr,g1,gn
c
cjs 26.02.2010      if(i2-i1)13,1,2

      if((i2-i1).lt.0) then
         return
      else if((i2-i1).gt.0) then
         n=0
      else
         tau(2,i1)=x1
      endif

      do 3 i=i1,i2
      n=n+1
      b(n)=tau(1,i)
      do 3 j=1,2
 3    a(j,n)=coef(j,i)
      do 4 j=1,3
 4    ap(j)=coef(j+2,i2)
      n1=n-1
c
c   Arrays ap(*,1), a, and ap(*,2) comprise n+2 x n+2 penta-diagonal
c   matrix A.  Let x1, tau, and xn comprise corresponding n+2 vector b.
c   Then, A * g = b, may be solved for n+2 vector g such that
c   interpolation I is given by I(p) = sum(i=0,n+1) g;i * G;i(p).
c
c   Eliminate the lower triangular portion of A to form A'.  A
c   corresponding transformation applied to vector b is stored in
c   a(4,*).
      alr=a(1,1)/coef(3,i1)
      a(1,1)=1.d0-coef(4,i1)*alr
      a(2,1)=a(2,1)-coef(5,i1)*alr
      b(1)=b(1)-x1*alr
      j=1
      do 5 i=2,n
      alr=a(1,i)/a(1,j)
      a(1,i)=1.d0-a(2,j)*alr
      b(i)=b(i)-b(j)*alr
 5    j=i
      alr=ap(1)/a(1,n1)
      ap(2)=ap(2)-a(2,n1)*alr
      gn=xn-b(n1)*alr
      alr=ap(2)/a(1,n)
c   Back solve the upper triangular portion of A' for coefficients g;i.
c   When finished, storage g(2), a(4,*), g(5) will comprise vector g.
      gn=(gn-b(n)*alr)/(ap(3)-a(2,n)*alr)
      b(n)=(b(n)-gn*a(2,n))/a(1,n)
      j=n
      do 6 i=n1,1,-1
      b(i)=(b(i)-b(j)*a(2,i))/a(1,i)
 6    j=i
cj.s.      g1=(x1-coef(4,i1)*b(1)-coef(5,i1)*b(2))/coef(3,i1)
c
      tau(2,i1)=x1
      is=i1+1
      ie=i2-1
      j=1
      do 7 i=is,ie
      j=j+1
 7    tau(2,i)=coef(3,i)*b(j-1)+coef(4,i)*b(j)+coef(5,i)*b(j+1)
      tau(2,i2)=xn
      return
      end

      subroutine trtm(delta,n,tt,dtdd,dtdh,dddp,phnm)
      save 
      include 'ttlim.h'
      character*(*) phnm(maxp)
      character*8 ctmp(maxp)
      real*8 tt(maxp),dtdd(maxp),dtdh(maxp),dddp(maxp)
      dimension tmp1(maxp),tmp2(maxp),tmp3(maxp),tmp4(maxp),iptr(maxp)
      real*8 us,pt,tau,xlim,xbrn,dbrn
      real*8 x(3),cn,dtol,pi,pi2,dn
      common/tabc/dn,us(2),pt(jout),tau(4,jout),xlim(2,jout),
     1 xbrn(jbrn,3),dbrn(jbrn,2),xn,pn,tn,hn,jndx(jbrn,2),idel(jbrn,3),
     2 mbr1,mbr2
c     data cn,dtol,atol,pi,pi2/.017453292519943296d0,1.d-6,.005,
c    1 3.1415926535897932d0,6.2831853071795865d0/

      data dtol,atol/1.d-6,.005/
c
      pi = 4.d0*datan(1.d0)
      pi2 = 2.d0*pi
      cn = pi / 180.d0

      n=0
      if(mbr2.le.0) return
      x(1)=dmod(dabs(cn*dble(delta)),pi2)
      if(x(1).gt.pi) x(1)=pi2-x(1)
      x(2)=pi2-x(1)
      x(3)=x(1)+pi2
      if(dabs(x(1)).gt.dtol) go to 9
      x(1)=dtol
      x(3)=-10.d0
 9    if(dabs(x(1)-pi).gt.dtol) go to 7
      x(1)=pi-dtol
      x(2)=-10.d0
 7    do 1 j=mbr1,mbr2
 1    if(jndx(j,2).gt.0) call findtt(j,x,n,tmp1,tmp2,tmp3,tmp4,ctmp)
c
c js 26.02.2010
c
c     if(n-1)3,4,5
c
      if(n-1.lt.0) then
         return
      else if(n-1.gt.0) then
        call r4sort(n,tmp1,iptr)
      else
        iptr(1)=1
      endif
      k=0

      do 2 i=1,n
      j=iptr(i)
      if(k.le.0) go to 8
      if(phnm(k).eq.ctmp(j).and.abs(sngl(tt(k))-tmp1(j)).le.atol) 
     +      go to 2
 8    k=k+1
      tt(k)  =dble(tmp1(j))
      dtdd(k)=dble(tmp2(j))
      dtdh(k)=dble(tmp3(j))
      dddp(k)=dble(tmp4(j))
      phnm(k)=ctmp(j)
 2    continue
      n=k
      return
      end

      subroutine findtt(jb,x0,n,tt,dtdd,dtdh,dddp,phnm)
      save 
      include 'ttlim.h'
      character*(*) phnm(maxp)
      character*8 phcd
      character*67 msg
      dimension tt(maxp),dtdd(maxp),dtdh(maxp),dddp(maxp)
      real*8 us,pt,tau,xlim,xbrn,dbrn,dn,dsgn
      real*8 x,x0(3),p0,p1,arg,dp,dps,delp,tol,ps,deps
      common/tabc/dn,us(2),pt(jout),tau(4,jout),xlim(2,jout),
     1 xbrn(jbrn,3),dbrn(jbrn,2),xn,pn,tn,hn,jndx(jbrn,2),idel(jbrn,3),
     2 mbr1,mbr2
      common/pcdc/phcd(jbrn)
      data tol/3.d-6/,deps/1.d-10/
c
      nph=iabs(idel(jb,3))
      hsgn=isign(1,idel(jb,3))*hn
      dsgn=(-1.d0)**dble(idel(jb,1))*dn
      dpn=-1./tn
      do 10 ij=idel(jb,1),idel(jb,2)
      x=x0(ij)
      dsgn=-dsgn
      if(x.lt.xbrn(jb,1).or.x.gt.xbrn(jb,2)) go to 12
      j=jndx(jb,1)
      is=j+1
      ie=jndx(jb,2)
      do 1 i=is,ie
      if(x.le.xlim(1,j).or.x.gt.xlim(2,j)) go to 8
      le=n
      p0=pt(ie)-pt(j)
      p1=pt(ie)-pt(i)
      delp=dmax1(tol*(pt(i)-pt(j)),1.d-3)
      if(dabs(tau(3,j)).gt.1.d-30) go to 2
      dps=(x-tau(2,j))/(1.5d0*tau(4,j))
      dp=dsign(dps*dps,dps)
      dp0=sngl(dp)
      if(dp.lt.p1-delp.or.dp.gt.p0+delp) go to 9
      if(n.ge.maxp) go to 13
      n=n+1
      ps=pt(ie)-dp
      tt(n)=tn*sngl(tau(1,j)+dp*(tau(2,j)+dps*tau(4,j))+ps*x)
      dtdd(n)=sngl(dsgn*ps)
cjs   dtdh(n)=hsgn*sqrt(abs(sngl(us(nph)*us(nph)-ps*ps)))
      dtdh(n)=hsgn*sngl(dsqrt(dabs(us(nph)*us(nph)-ps*ps)))
      dddp(n)=dpn*sngl(.75d0*tau(4,j)/dmax1(dabs(dps),deps))
      phnm(n)=phcd(jb)
      in=index(phnm(n),'ab')
      if(in.le.0) go to 8
      if(ps.le.xbrn(jb,3)) phnm(n)(in:)='bc'
      go to 8
 2    do 4 jj=1,2
      go to (5,6),jj
 5    arg=9.d0*tau(4,j)*tau(4,j)+32.d0*tau(3,j)*(x-tau(2,j))
      if(arg.ge.0.d0) go to 3
      write(msg,100)arg
 100  format('Bad sqrt argument:',1pd11.2,'.')
      call warn(msg(1:30))
 3    dps=-(3.d0*tau(4,j)+dsign(dsqrt(dabs(arg)),tau(4,j)))/(8.d0*
     1 tau(3,j))
      dp=dsign(dps*dps,dps)
      dp0=sngl(dp)
      go to 7
 6    dps=(tau(2,j)-x)/(2.d0*tau(3,j)*dps)
      dp=dsign(dps*dps,dps)
 7    if(dp.lt.p1-delp.or.dp.gt.p0+delp) go to 4
      if(n.ge.maxp) go to 13
      n=n+1
      ps=pt(ie)-dp
      tt(n)=tn*sngl(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+
     1      dps*tau(4,j))+ps*x)
      dtdd(n)=sngl(dsgn*ps)
      dtdh(n)=hsgn*sngl(dsqrt(dabs(us(nph)*us(nph)-ps*ps)))
      dddp(n)=dpn*sngl(2.d0*tau(3,j)+.75d0*tau(4,j)/
     1        dmax1(dabs(dps),deps))
      phnm(n)=phcd(jb)
      in=index(phnm(n),'ab')
      if(in.le.0) go to 4
      if(ps.le.xbrn(jb,3)) phnm(n)(in:)='bc'
 4    continue
 9    if(n.gt.le) go to 8
      write(msg,101)phcd(jb),x,dp0,dp,p1,p0
 101  format('Failed to find phase:  ',a,f8.1,4f7.4)
      call warn(msg)
 8    j=i
 1    continue
c
 12   if(x.lt.dbrn(jb,1).or.x.gt.dbrn(jb,2)) go to 10
      if(n.ge.maxp) go to 13
      j=jndx(jb,1)
      i=jndx(jb,2)
      dp=pt(i)-pt(j)
      dps=dsqrt(dabs(dp))
      n=n+1
      tt(n)=tn*sngl(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+
     1      dps*tau(4,j))+pt(j)*x)
      dtdd(n)=sngl(dsgn*pt(j))
      dtdh(n)=hsgn*sngl(dsqrt(dabs
     1        (us(nph)*us(nph)-pt(j)*pt(j))))
      dddp(n)=dpn*sngl(2.d0*tau(3,j)+.75d0*tau(4,j)/
     1        dmax1(dabs(dps),deps))
      ln=index(phcd(jb),'ab')-1
      if(ln.le.0) ln=index(phcd(jb),' ')-1
      if(ln.le.0) ln=len(phcd(jb))
      phnm(n)=phcd(jb)(1:ln)//'dif'
 10   continue
      return
 13   write(msg,102)maxp
 102  format('More than ',i3,' arrivals found.')
      call warn(msg(1:28))
      return
      end
      subroutine assign(lu,mode,ia)
c
c $$$$$ calls no other routine $$$$$
c
c   Subroutine assign opens (connects) logical unit lu to the disk file
c   named by the character string ia with mode mode.  If iabs(mode) = 1,
c   then open the file for reading.  If iabs(mode) = 2, then open the
c   file for writing.  If iabs(mode) = 3, then open a scratch file for
c   writing.  If mode > 0, then the file is formatted.  If mode < 0,
c   then the file is unformatted.  All files opened by assign are
c   assumed to be sequential.  Programmed on 3 December 1979 by
c   R. Buland.
c
c     save
      character*(*) ia
      logical exst
      integer lu,mode
c
c js  na = index(ia,' ')-1
c js  if(na.le.0) na=len(ia)
c
c     Because of WINDOWS implementation, we have to accept blanks in 
c     file and path names. FUNCTION LEN is searching the first blank
c     from the left.
c     Solution: We check the file name length from the right and
c               remove trailing blanks with FUNCTION TRIM .
c

      if(mode.ge.0) nf=1
      if(mode.lt.0) nf=2
      ns=iabs(mode)
      if(ns.le.0.or.ns.gt.3) ns=3
      go to (1,2),nf
 1    go to (11,12,13),ns
 11   open(lu,file=trim(ia),status='old',form='formatted')
      rewind lu
      return
 12   inquire(file=trim(ia),exist=exst)
      if(exst) go to 11
 13   open(lu,file=trim(ia),status='new',form='formatted')
      return
 2    go to (21,22,23),ns
 21   open(lu,file=trim(ia),status='old',form='unformatted')
      rewind lu
      return
 22   inquire(file=trim(ia),exist=exst)
      if(exst) go to 21
 23   open(lu,file=trim(ia),status='new',form='unformatted')
      return
      end
      subroutine r4sort(n,rkey,iptr)
c
c $$$$$ calls no other routine $$$$$
c
c   R4sort sorts the n elements of array rkey so that rkey(i), 
c   i = 1, 2, 3, ..., n are in asending order.  R4sort is a trivial
c   modification of ACM algorithm 347:  "An efficient algorithm for
c   sorting with minimal storage" by R. C. Singleton.  Array rkey is
c   sorted in place in order n*alog2(n) operations.  Coded on
c   8 March 1979 by R. Buland.  Modified to handle real*4 data on
c   27 September 1983 by R. Buland.
c
c     save
      dimension rkey(*),iptr(*),il(10),iu(10)
c   Note:  il and iu implement a stack containing the upper and
c   lower limits of subsequences to be sorted independently.  A
c   depth of k allows for n<=2**(k+1)-1.
      if(n.le.0) return
      do 1 i=1,n
 1    iptr(i)=i
      if(n.le.1) return
      r=.375
      m=1
      i=1
      j=n
c
c   The first section interchanges low element i, middle element ij,
c   and high element j so they are in order.
c
 5    if(i.ge.j) go to 70
 10   k=i
c   Use a floating point modification, r, of Singleton's bisection
c   strategy (suggested by R. Peto in his verification of the
c   algorithm for the ACM).
      if(r.gt..58984375) go to 11
      r=r+.0390625
      go to 12
 11   r=r-.21875
 12   ij=i+nint((j-i)*r)
      if(rkey(iptr(i)).le.rkey(iptr(ij))) go to 20
      it=iptr(ij)
      iptr(ij)=iptr(i)
      iptr(i)=it
 20   l=j
      if(rkey(iptr(j)).ge.rkey(iptr(ij))) go to 39
      it=iptr(ij)
      iptr(ij)=iptr(j)
      iptr(j)=it
      if(rkey(iptr(i)).le.rkey(iptr(ij))) go to 39
      it=iptr(ij)
      iptr(ij)=iptr(i)
      iptr(i)=it
 39   tmpkey=rkey(iptr(ij))
      go to 40
c
c   The second section continues this process.  K counts up from i and
c   l down from j.  Each time the k element is bigger than the ij
c   and the l element is less than the ij, then interchange the
c   k and l elements.  This continues until k and l meet.
c
 30   it=iptr(l)
      iptr(l)=iptr(k)
      iptr(k)=it
 40   l=l-1
      if(rkey(iptr(l)).gt.tmpkey) go to 40
 50   k=k+1
      if(rkey(iptr(k)).lt.tmpkey) go to 50
      if(k.le.l) go to 30
c
c   The third section considers the intervals i to l and k to j.  The
c   larger interval is saved on the stack (il and iu) and the smaller
c   is remapped into i and j for another shot at section one.
c
      if(l-i.le.j-k) go to 60
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 80
 60   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 80
c
c   The fourth section pops elements off the stack (into i and j).  If
c   necessary control is transfered back to section one for more
c   interchange sorting.  If not we fall through to section five.  Note
c   that the algorighm exits when the stack is empty.
c
 70   m=m-1
      if(m.eq.0) return
      i=il(m)
      j=iu(m)
 80   if(j-i.ge.11) go to 10
      if(i.eq.1) go to 5
      i=i-1
c
c   The fifth section is the end game.  Final sorting is accomplished
c   (within each subsequence popped off the stack) by rippling out
c   of order elements down to their proper positions.
c
 90   i=i+1
      if(i.eq.j) go to 70
      if(rkey(iptr(i)).le.rkey(iptr(i+1))) go to 90
      k=i
      kk=k+1
      ib=iptr(kk)
 100  iptr(kk)=iptr(k)
      kk=k
      k=k-1
      if(rkey(ib).lt.rkey(iptr(k))) go to 100
      iptr(kk)=ib
      go to 90
      end
c
      subroutine warn(msg)
      character*(*) msg
      write(*,100) msg
 100  format(1x,a)
      return
      end
c

      subroutine oneray(phnm,dtdd,xcor,tcor,onerayf)
c
c   Given a phase code, phnm, oneray returns the distance, xcor, in
c   degrees and the travel time, tcor, in seconds corresponding to ray
c   parameter dtdd in km/(km*s).
c
      save
      real*8 cn
c     parameter(cn=57.295779d0)
      include 'ttlim.h'
      character*(*) phnm
      character*8 phcd,phsyn
c     real*8 us,pt,tau,xlim,xbrn,dbrn,zs,pk,pu,pux,tauu,xu,
c    1 px,xt,taut,tauc,xc,coef,tcoef,tp
      real*8 us,pt,tau,xlim,xbrn,dbrn
      real*8 x,dp,dps,ps,dn
      common/tabc/dn,us(2),pt(jout),tau(4,jout),xlim(2,jout),
     1 xbrn(jbrn,3),dbrn(jbrn,2),xn,pn,tn,hn,jndx(jbrn,2),idel(jbrn,3),
     2 mbr1,mbr2
c     common/brkc/zs,pk(jseg),pu(jtsm0,2),pux(jxsm,2),tauu(jtsm,2),
c    1 xu(jxsm,2),px(jbrn,2),xt(jbrn,2),taut(jout),coef(5,jout),
c    2 tauc(jtsm),xc(jxsm),tcoef(5,jbrna,2),tp(jbrnu,2),odep,
c    3 fcs(jseg,3),nin,nph0,int0(2),ki,msrc(2),isrc(2),nseg,nbrn,ku(2),
c    4 km(2),nafl(jseg,3),indx(jseg,2),kndx(jseg,2),iidx(jseg),
c    5 jidx(jbrn),kk(jseg)
      common/pcdc/phcd(jbrn)

      logical onerayf
c

      cn = 45.d0 / datan(1.d0)

      onerayf=.true.
      phsyn=phnm
c
c changed and added to get the correct PKP branch J.S.
c
c     j=index(phsyn,'bc')
c     if(j.gt.0) phsyn(j:j+1)='ab'
      jp=index(phsyn,'bc')
      ja=index(phsyn,'ab')
      if(jp.gt.0) phsyn(jp:jp+1)='ab'

c
c   Find the branch.
c     print *,'mbr1 mbr2',mbr1,mbr2
      do jb=mbr1,mbr2
c       print *,'jb phcd xbrn',jb,'  ',phcd(jb),xbrn(jb,1)
        if(phcd(jb).eq.phsyn) then
c
c   Got the branch.  See if the ray parameter is OK.
          is=jndx(jb,1)+1
          ie=jndx(jb,2)
          ps=dabs(dble(dtdd))/dn
c         print *,'jb is ie dtdd dn ps',jb,is,ie,dtdd,dn,ps
          if(ps.ge.pt(is-1).and.ps.le.pt(ie)) then
c
c   The ray parameter is OK.  Find the right ray parameter interval.
            do i=is,ie
c             print *,'i pt',i,pt(i)
c
c   Got the ray parameter interval.  Interpolate.
              if(ps.le.pt(i)) then
c
c  added to get the right PKP branch J.S.
c
              if(jp.gt.0 .and. ps.gt.xbrn(jb,3)) go to 10
              if(ja.gt.0 .and. ps.le.xbrn(jb,3)) go to 10

                j=i-1
                dp=pt(ie)-ps
                dps=dsqrt(dabs(dp))
                x=tau(2,j)+2.d0*dp*tau(3,j)+1.5d0*dps*tau(4,j)
c               print *,'j pt dp dps x',j,pt(ie),dp,dps,x
                tcor=sngl(dble(tn)*(tau(1,j)+dp*(tau(2,j)+dp*tau(3,j)+
     1           dps*tau(4,j))+ps*x))
                xcor=sngl(cn*x)
c               print *,'oneray xcor tcor',onerayf,xcor,tcor
                return
              endif
            enddo
          endif
        endif
      enddo
c
c   Didn't get it.
cj.s.  onerayf=.false.
10    onerayf=.false.
      xcor=0.
      tcor=0.
c     print *,'oneray xcor tcor',onerayf,xcor,tcor
      return
      end
c
