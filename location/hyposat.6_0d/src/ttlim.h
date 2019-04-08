c   The "j" parameters (1 st line) are intended to be user settable:
c        jsrc   Maximum number of discrete model slowness samples above
c               the maximum source depth of interest.
c        jseg   Maximum number of different types of travel-times
c               considered.
c        jbrn   Maximum number of different travel-time branches to be
c               searched.
c        jout   Maximum length of all travel-time branches strung
c               together.
c        jtsm   Maximum length of the tau depth increments.
c        jxsm   Maximum number of x-values needed for the depth
c               increments.
c        jbrnu  Maximum length of the up-going branches.
c        jbrna  Maximum length of branches which may need
c               re-interpolation.
c
c        maxp   maximum number of phases at one distance allowed
c               The PARAMETER mphas in ttimes.h must have the same value!!
c
      parameter(jsrc=250,jseg=75,jbrn=150,jout=3000,jtsm=450,jxsm=jbrn,
     1 jbrnu=jbrn,jbrna=jbrn,maxp=120)
c   A few derived parameters are also needed.
      parameter(jrec=jtsm+jxsm,jtsm0=jtsm+1)
