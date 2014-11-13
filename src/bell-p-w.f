      SUBROUTINE DPSS(NMAX, KMAX, N, W, V, SIG, TOTIT, SINES, VOLD,   
     * U, SCR1, IFAULT)
C
C  CALCULATES DISCRETE PROLATE SPHEROIDAL SEQUENCES FOR USE AS DATA
C  TAPERS.
C
C  FORTRAN 77 
C
C  SUBMITTED BY BRAD BELL, DON PERCIVAL AND ANDREW WALDEN.
C  Comments/queries to
C  dbp@apl.washington.edu  OR  a.walden@ic.ac.uk
C
C  This software may be freely used for non-commercial purposes and can
C  be freely distributed.
C
C  Equation numbers and comments refer to the article:
C
C  Bell, B., Percival, D.B. and Walden, A.T. "Calculating Thomson's
C    Spectral Multitapers by Inverse Iteration", J. Comput. and Graph.
C    Stat., 1993.
C
C  Calls auxiliary routines SYTOEP and SPOL also given below.
C
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     NMAX    integer    input:  maximum possible length of taper
C
C     KMAX    integer    input:  dpss orders 0 to KMAX required
C
C     N       integer    input:  length of sequence to generate
C
C     W       real       input:  half-bandwidth, W < 1/2
C
C     V(NMAX, KMAX+1)       
C            real array output:  columns contain tapers
C
C     SIG(KMAX+1)
C            real array output:  eigenvalues are 1+SIG(j)
C
C     TOTIT   integer   output:  total number of iterations
C
C     SINES(0:N-1)
C            real array          work array
C
C     VOLD(N) real array         work array
C
C     U(N)    real array         work array
C
C     SCR1(N) real array         work array
C
C     IFAULT  integer   output:   0 indicates success
C                                 1 if W > 1/2
C                                 2 if N < 2
C                                 3 if NMAX < N; matrix too small
C                                 4 if KMAX < 0
C                                 5 failure in SYTOEP
C                                 6 > MAXIT its required for some order
C                                 (Output values are undefined for
C                                 IFAULT in the range 1 to 5.)
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER J, K, K1, KMAX, M, N, TOTIT, ISIG, ILOW,
     *  IHIG, IT, MAXIT, IFAIL, IFAULT
      DOUBLE PRECISION V(NMAX,0:KMAX), SIG(0:KMAX), U(N),  
     * VOLD(N), SINES(0:N-1), SCR1(N), ZERO, ONE, TWO, FOUR, 
     * EPS, DELTA, PI, W, ROOTN, PROJ, SNORM, SSNORM, DIFF,
     * SUM, HALF, RONE 
      DATA ZERO, HALF, ONE,TWO,FOUR,EPS
     * /0.0D0, 0.5D0, 1.0D0,2.0D0,4.0D0,0.5D-6/
C
C  initialize max number of iterations flag
C
      IFAIL=0 
C
C check input parameters
C      
      IFAULT=1  
      IF(W .GT. HALF) RETURN
      IFAULT=2
      IF(N .LT. 2) RETURN
      IFAULT=3
      IF(NMAX .LT. N) RETURN
      IFAULT=4
      IF(KMAX .LT. 0 .OR. KMAX .GT. N-1) RETURN
C
C  set up SINES so that S in eqn. (1) is given by
C  S(n,m)=SINES(n-m) for n not equal to m. 
C
      PI=FOUR*ATAN(ONE)
      DO 5 M=1,N-1
        SINES(M)=SIN(TWO*PI*W*M)/(PI*M)
 5    CONTINUE
C
C  set total iteration counter and constant
C
      TOTIT=0
      ROOTN=SQRT(FLOAT(N))
      RONE=ONE/ROOTN
C
C  major loop over dpss orders 0 to KMAX
C
C  modify SINES(0) so that B_k in Section 2.2 is given by
C  B_k(n,m)=SINES(n-m)
C
      DO 200 K=0,KMAX
         IF(K .EQ. 0) THEN
              SINES(0)=TWO*W-ONE
         ELSE
              SINES(0)=TWO*W-(ONE+SIG(K-1))
         END IF
C
C  define suitable starting vector for inverse iteration;
C  see Section 2.2.
C
         ISIG=1
         K1=K+1
         DO 15 J=1, K1
                   ILOW=((J-1)*N/K1)+1
                   IHIG=(J*N/K1)
                   DO 20 JJ=ILOW,IHIG
                        U(JJ)=ISIG*RONE
 20                CONTINUE
                   ISIG=ISIG*(-1)
 15      CONTINUE
         IF(MOD(K,2).GT.0 .AND. MOD(N,2).GT.0) U((N/2)+1)=ZERO
C
C  maximum number of iterations
C
         MAXIT=(K+3)*ROOTN 
         IT=0
C
C  carry out inverse iteration
C
         DO 180 IT=1,MAXIT
C
C  copy U into old V; VOLD = previous iterate
C
         DO 50 J=1,N
 50        VOLD(J)=U(J)
C
C  solve symmetric Toeplitz matrix equation B_k*U=VOLD for U
C
         CALL SYTOEP(N, SINES, VOLD, U, SCR1, IFAIL)
C
C  check no problems 
C 
         IFAULT=5
         IF(IFAIL .NE. 0) RETURN
C
C  new vector must be orthogonal to previous eigenvectors
C
         IF(K .GT. 0) THEN
              DO 80 K1=0,K-1
C
C  projection of U onto V(*,K1):
C
              PROJ=ZERO
              DO 85 J=1,N
 85             PROJ=PROJ + U(J)*V(J,K1)
C
C subtract projection
C
                DO 90 J=1,N
 90               U(J)=U(J) - PROJ*V(J,K1)
C
 80           CONTINUE
         END IF
C
C  normalize
C
         SNORM=ZERO
         DO 100 J=1,N
100        SNORM=SNORM+U(J)*U(J)
         SSNORM=SQRT(SNORM)
         DO 105 J=1,N
105               U(J)=U(J)/SSNORM
C
C  check for convergence
C
         SUM=ZERO
         DIFF=ZERO
         DO 120 J=1,N
C         
C  first previous-current:
C
           DIFF=DIFF+(VOLD(J)-U(J))**2
C
C  next, previous+current
C            
 120     SUM=SUM+(VOLD(J)+U(J))**2
         DELTA=SQRT(MIN(DIFF,SUM))
         IF(DELTA .LE. EPS) GOTO 190
180      CONTINUE
C
C  if here, max number of iterations exceeded for this order dpss
C  
         IT=MAXIT
         IFAIL=1
190      CONTINUE
         TOTIT=TOTIT+IT
         IF(SUM .LT. DIFF) THEN
               IF(K .EQ. 0) THEN
                    SIG(0)= - ONE/SSNORM
               ELSE
                    SIG(K)=SIG(K-1) - ONE/SSNORM
               END IF
         ELSE
               IF(K .EQ. 0) THEN
                    SIG(0)= ONE/SSNORM
               ELSE
                    SIG(K)=SIG(K-1) + ONE/SSNORM
               END IF
         END IF
C
C  ensure tapers satisfy Slepian convention
C         
         CALL SPOL(N, U, K, IFAULT)
         DO 220 J=1,N
220        V(J,K)=U(J)
C
200   CONTINUE
C
C   if one order of dpss did not converge set IFAULT to 6
C
      IFAULT=0
      IF( IFAIL .EQ. 1) IFAULT=6 
      RETURN
      END


      SUBROUTINE SYTOEP(N, R, G, F, W, IFAULT)
C
C  FINDS FILTER CORRESPONDING TO A SYMMETRIC TOEPLITZ MATRIX
C  WITH FIRST ROW R(.) AND CROSSCORRELATION VECTOR G(.)
C
C  I.E., "R" F = G
C
C  To be used with DPSS and SPOL. See
C  Bell, B., Percival, D.B. and Walden, A.T. "Calculating Thomson's
C    Spectral Multitapers by Inverse Iteration", J. Comput. and Graph.
C    Stat., 1993.
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     N       integer    input:  dimension of Toeplitz matrix and
C                                 cross-correlation vector
C
C     R(N)    real       input:  autocovariances from lag 0 to N-1
C
C     G(N)    real       input:  cross-correlation vector
C
C     F(N)    real      output:  required filter
C
C     W(N)    real       input:  work array
C
C     IFAULT  integer   output:  0 indicates successful
C                                1 if N < 1
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C   THIS PROGRAM IS A SUBSTANTIALLY CORRECTED AND MODIFIED VERSION OF 
C   "EUREKA" IN ROBINSON, E.A. (1967) MULTICHANNEL TIME SERIES ANALYSIS 
C   WITH DIGITAL COMPUTER PROGRAMS, HOLDEN-DAY.
C
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DOUBLE PRECISION R(0:N-1), G(N), F(N), W(N), V, D, Q,
     *  HOLD, ZERO, ONE
      DATA ZERO,ONE /0.0D0,1.0D0/
C
C  check for special "matrix" sizes
C
      IFAULT=1
      IF(N .LT. 1) RETURN
      V=R(0)
      F(1)=G(1)/V
      IFAULT=0
      IF(N .EQ. 1) RETURN
C
      D=R(1)
      W(1)=ONE
      Q=F(1)*R(1)
      DO 5 L=2,N
        W(L)=-D/V
        IF(L .GT. 2) THEN
             L1=(L-2)/2
             L2=L1+1
             IF(L .NE. 3) THEN
                  DO 10 J=2,L2
                  HOLD=W(J)
                  K=L-J+1
                  W(J)=W(J)+W(L)*W(K)
 10               W(K)=W(K)+W(L)*HOLD
             END IF
             IF((2*L1 .NE. L-2) .OR. L .EQ. 3) W(L2+1)=
     *       W(L2+1)+W(L)*W(L2+1)
        END IF
      V=V+W(L)*D
      F(L)=(G(L)-Q)/V
      L3=L-1
      DO 15 J=1,L3
           K=L-J+1
 15   F(J)=F(J)+F(L)*W(K)
      IF(L .EQ. N) RETURN
      D=ZERO
      Q=ZERO
      DO 5 I=1,L
          K=L-I+2
          D=D+W(I)*R(K-1)
  5   Q=Q+F(I)*R(K-1)
      RETURN
      END



      SUBROUTINE SPOL(N, V, K, IFAULT)
C
C  SCALES THE DISCRETE PROLATE SPHEROIDAL SEQUENCE AND SETS THE 
C  POLARITY TO AGREE WITH SLEPIAN'S CONVENTION.
C
C  To be used with DPSS and SYTOEP. See
C  Bell, B., Percival, D.B. and Walden, A.T. "Calculating Thomson's
C    Spectral Multitapers by Inverse Iteration", J. Comput. and Graph.
C    Stat., 1993. (Section 1.2.)
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C  N      Integer             input        length of dpss sequence
C
C  V      Real(N)             input        eigenvector (dpss) with unit energy
C                            output        unit energy dpss conforming to
C                                          Slepian's polarity convention
C
C  K      Integer             input        the order of the dpss 0=<K=<N-1
C
C  IFAULT Integer            output        0 indicates successful
C                                          1 indicates N < 1
C
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      INTEGER N, I, L, K, IFAULT
      DOUBLE PRECISION V(N), RN, STD, DSUM, DWSUM, ZERO, ONE, 
     * TWO, ZFLOAT
C
      DATA ZERO, ONE, TWO/0.0D0, 1.0D0, 2.0D0/
C
C
      IFAULT=1
      IF(N .LT. 1) RETURN
      IFAULT=0
      RN=DBLE(N)
      DSUM=ZERO
      DWSUM=ZERO
      DO 5 L=1,N
        DSUM=DSUM+V(L)
  5   DWSUM=DWSUM+V(L)*(RN-ONE-TWO*DBLE(L-1))
      IF(((MOD(K,2).EQ.0).AND.(DSUM.LT.ZERO)).OR.
     *     ((MOD(K,2).EQ.1).AND.(DWSUM.LT.ZERO))) THEN
      DO 10 L=1,N
 10          V(L)=-V(L)
      END IF
      RETURN
      END
