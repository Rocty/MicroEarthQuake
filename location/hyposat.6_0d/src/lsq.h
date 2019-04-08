c
C       linear equation sytem: GGL*AAL=DDL 
C               (n equations for m unknowns, n>=M)
C       kernel matrix GGL (n times m), model AAL (m times 1),
C       data DDL (n times 1)
C       VVL = STANDARD DEVIATIONS (sqrt of variances) of AAL
C       RRL = RESIDUALS = GGL*AAL-DDL
c       lsqerr = set to .ne. 0 in case of an error

      PARAMETER (NNL=2000,MML=2)
      DIMENSION GGL(NNL,MML),AAL(MML),DDL(NNL),VVL(MML),RRL(NNL)
      integer   lsqerr
      common  /lsq/ ggl,ddl,aal,vvl,rrl,lsqerr
