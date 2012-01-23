C***********************************************************************
C
C     File of the REAL              Level 2 BLAS routines:  
C
C      SGEMV, SGBMV, SSYMV, SSBMV, SSPMV, STRMV, STBMV, STPMV,
C      SGER , SSYR , SSPR ,
C      SSYR2, SSPR2,
C      STRSV, STBSV, STPSV.
C
C     See: 
C
C        Dongarra J. J., Du Croz J. J., Hammarling S. and Hanson R. J.. 
C        A proposal for an extended set of Fortran Basic Linear Algebra
C        Subprograms. Technical Memorandum No.41 (revision 1),
C        Mathematics and Computer Science Division, Argone National
C        Laboratory, 9700 South Cass Avenue, Argonne, Illinois 60439,
C        USA, or NAG Technical Report TR4/85, Numerical Algorithms Group
C        Inc., 1101 31st Street, Suite 100, Downers Grove, Illinois
C        60606-1263, USA.
C
C***********************************************************************
C
      SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      CHARACTER *1 TRANS
      INTEGER M,N,LDA,INCX,INCY
      REAL ALPHA,A(LDA,*),X(*),BETA,Y(*)
*
*  Purpose
*  =======
*
*  SGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, 
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows: 
*
*              TRANS = 'N'  y := alpha*A*x + beta*y.
*
*              TRANS = 'T'  y := alpha*A'*x + beta*y.
*
*              TRANS = 'C'  y := alpha*A'*x + beta*y
*.
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           max(m,1).
*           Unchanged on exit.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the 
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL
*           On entry, BETA specifies the scalar beta. When BETA is 
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*  Note that TRANS, M, N and LDA must be such that the value of the
*  LOGICAL variable OK in the following statement is true.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-August-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER KX,KY,LENX,LENY 
      REAL ONE,ZERO 
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP
      LOGICAL OK,LSAME
      OK = (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. ((M.GT.0) .AND. (N.GT.0) .AND.
     .     (LDA.GE.M))
*
*     Quick return if possible.
*
      IF (((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE)) .OR. .NOT. OK) RETURN
*
*     Set LENX and LENY, the lengths of the vectors x and y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
*
      ELSE
          LENX = M
          LENY = N
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          IF (BETA.NE.ONE) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10,I = 1,LENY
                     Y(I) = ZERO
   10             CONTINUE
*
              ELSE
                  DO 20,I = 1,LENY
                     Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
*
          END IF
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (LENX-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (LENY-1)*INCY
          END IF
*
          IF (BETA.NE.ONE) THEN
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30,I = 1,LENY
                     Y(IY) = ZERO
                     IY = IY + INCY
   30             CONTINUE
*
              ELSE
                  DO 40,I = 1,LENY
                     Y(IY) = BETA*Y(IY) 
                     IY = IY + INCY
   40             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     DO 50,I = 1,M
                        Y(I) = Y(I) + TEMP*A(I,J) 
   50                CONTINUE 
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              DO 80,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IY = KY
                     DO 70,I = 1,M
                        Y(IY) = Y(IY) + TEMP*A(I,J)
                        IY = IY + INCY
   70                CONTINUE 
                 END IF
*
                 JX = JX + INCX
   80         CONTINUE
          END IF
*
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 100,J = 1,N
                 TEMP = ZERO
                 DO 90,I = 1,M
                    TEMP = TEMP + A(I,J)*X(I)
   90            CONTINUE
                 Y(J) = Y(J) + ALPHA*TEMP
  100         CONTINUE
*
          ELSE
              JY = KY
              DO 120,J = 1,N
                 TEMP = ZERO
                 IX = KX
                 DO 110,I = 1,M
                    TEMP = TEMP + A(I,J)*X(IX)
                    IX = IX + INCX
  110            CONTINUE
                 Y(JY) = Y(JY) + ALPHA*TEMP
                 JY = JY + INCY
  120         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SGEMV .
*
      END 
      SUBROUTINE SGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      CHARACTER *1 TRANS
      INTEGER M,N,KL,KU,LDA,INCX,INCY
      REAL ALPHA,A(LDA,*),X(*),BETA,Y(*)
*
*  Purpose
*  =======
*
*  SGBMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y, 
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n band matrix, with kl sub-diagonals and ku super-diagonals. 
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows: 
*
*              TRANS = 'N'  y := alpha*A*x + beta*y.
*
*              TRANS = 'T'  y := alpha*A'*x + beta*y.
*
*              TRANS = 'C'  y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  KL     - INTEGER.
*           On entry, KL specifies the number of sub-diagonals of the
*           matrix A. KL must satisfy  0 .le. KL.
*           Unchanged on exit.
*
*  KU     - INTEGER.
*           On entry, KU specifies the number of super-diagonals of the
*           matrix A. KU must satisfy   0 .le. KU. 
*           Unchanged on exit.
*
*  Users may find that efficiency of their application is enhanced by
*  adjusting the values of m and n so that KL .ge. max(0,m-n) and
*  KU .ge. max(0,n-m) or KL and KU so that KL .lt. m and KU .lt. n.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading ( kl + ku + 1 ) by n part of the
*           array A must contain the matrix of coefficients, supplied
*           column by column, with the leading diagonal of the matrix in
*           row ( ku + 1 ) of the array, the first super-diagonal
*           starting at position 2 in row ku, the first sub-diagonal
*           starting at position 1 in row ( ku + 2 ), and so on.
*           This placement of the data can be realized with the
*           following loops: 
*               DO 20 J =1,N
*                    K=KU+1-J 
*                    DO 10 I =MAX(1,J-KU),MIN(M,J+KL)
*                         A(K+I,J)=matrix entry of row I, column J. 
*     10             CONTINUE 
*     20        CONTINUE
*           Elements in the array A that do not correspond to elements
*           in the band matrix (such as the top left ku by ku triangle)
*           are not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           ( kl + ku + 1 ).
*           Unchanged on exit.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the 
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            . 
*           On entry, BETA specifies the scalar beta. When BETA is 
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry, the incremented array Y must contain the 
*           vector y. On exit, Y is overwritten by the updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 27-Sept-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTRINSIC MAX,MIN
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER K,KUP1,KX,KY,LENX,LENY
      REAL ONE,ZERO 
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP
      LOGICAL OK,LSAME
      OK = (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (M.GT.0) .AND. (N.GT.0) .AND.
     .     (KL.GE.0) .AND. (KU.GE.0) .AND.
     .     (LDA.GE. (KL+KU+1))
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN 
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
*
      ELSE
          LENX = M
          LENY = N
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the band part of A.
*
*     First form  y := beta*y  and set up the start points in  X  and  Y
*     if the increments are not both unity.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          IF (BETA.NE.ONE) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10,I = 1,LENY
                     Y(I) = ZERO
   10             CONTINUE
*
              ELSE
                  DO 20,I = 1,LENY
                     Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
*
          END IF
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (LENX-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (LENY-1)*INCY
          END IF
*
          IF (BETA.NE.ONE) THEN
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30,I = 1,LENY
                     Y(IY) = ZERO
                     IY = IY + INCY
   30             CONTINUE
*
              ELSE
                  DO 40,I = 1,LENY
                     Y(IY) = BETA*Y(IY) 
                     IY = IY + INCY
   40             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      IF (ALPHA.EQ.ZERO) RETURN
      KUP1 = KU + 1 
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     K = KUP1 - J
                     DO 50,I = MAX(1,J-KU),MIN(M,J+KL)
                        Y(I) = Y(I) + TEMP*A(K+I,J)
   50                CONTINUE 
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              DO 80,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IY = KY
                     K = KUP1 - J
                     DO 70,I = MAX(1,J-KU),MIN(M,J+KL)
                        Y(IY) = Y(IY) + TEMP*A(K+I,J)
                        IY = IY + INCY
   70                CONTINUE 
                 END IF
*
                 JX = JX + INCX
                 IF (J.GT.KU) KY = KY + INCY
   80         CONTINUE
          END IF
*
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 100,J = 1,N
                 TEMP = ZERO
                 K = KUP1 - J 
                 DO 90,I = MAX(1,J-KU),MIN(M,J+KL)
                    TEMP = TEMP + A(K+I,J)*X(I)
   90            CONTINUE
                 Y(J) = Y(J) + ALPHA*TEMP
  100         CONTINUE
*
          ELSE
              JY = KY
              DO 120,J = 1,N
                 TEMP = ZERO
                 IX = KX
                 K = KUP1 - J 
                 DO 110,I = MAX(1,J-KU),MIN(M,J+KL)
                    TEMP = TEMP + A(K+I,J)*X(IX)
                    IX = IX + INCX
  110            CONTINUE
                 Y(JY) = Y(JY) + ALPHA*TEMP
                 JY = JY + INCY
                 IF (J.GT.KU) KX = KX + INCX
  120         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SGBMV .
*
      END 
      SUBROUTINE SSYMV(UPLO,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      CHARACTER *1 UPLO
      INTEGER N,LDA,INCX,INCY 
      REAL ALPHA,A(LDA,*),X(*),BETA,Y(*)
*
*  Purpose
*  =======
*
*  SSYMV  performs the matrix-vector  operation 
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows: 
*
*              UPLO = 'U'          Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L'          Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U', the leading n by n 
*           upper triangular part of the array A must contain the upper 
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced.
*           Before entry with UPLO = 'L', the leading n by n
*           lower triangular part of the array A must contain the lower 
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(n,1).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            . 
*           On entry, BETA specifies the scalar beta. When BETA is 
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 27-Sept-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER KX,KY 
      REAL ONE,ZERO 
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP1,TEMP2
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0) .AND.
     .     (LDA.GE.N)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN 
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          IF (BETA.NE.ONE) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10,I = 1,N
                     Y(I) = ZERO
   10             CONTINUE
*
              ELSE
                  DO 20,I = 1,N
                     Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
*
          END IF
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
*
          IF (BETA.NE.ONE) THEN
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30,I = 1,N
                     Y(IY) = ZERO
                     IY = IY + INCY
   30             CONTINUE
*
              ELSE
                  DO 40,I = 1,N
                     Y(IY) = BETA*Y(IY) 
                     IY = IY + INCY
   40             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  y  when A is stored in upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 DO 50,I = 1,J - 1
                    Y(I) = Y(I) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(I)
   50            CONTINUE
                 Y(J) = Y(J) + TEMP1*A(J,J) + ALPHA*TEMP2
   60         CONTINUE
*
          ELSE
              IX = KX - INCX
              DO 80,J = 1,N
                 TEMP1 = ALPHA*X(IX+INCX)
                 TEMP2 = ZERO 
                 IX = KX
                 IY = KY
                 DO 70,I = 1,J - 1
                    Y(IY) = Y(IY) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(IX)
                    IX = IX + INCX
                    IY = IY + INCY
   70            CONTINUE
                 Y(IY) = Y(IY) + TEMP1*A(J,J) + ALPHA*TEMP2 
   80         CONTINUE
          END IF
*
      ELSE
*
*        Form  y  when A is stored in lower triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 100,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 Y(J) = Y(J) + TEMP1*A(J,J)
                 DO 90,I = J + 1,N
                    Y(I) = Y(I) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(I)
   90            CONTINUE
                 Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 120,J = 1,N
                 TEMP1 = ALPHA*X(JX)
                 TEMP2 = ZERO 
                 Y(JY) = Y(JY) + TEMP1*A(J,J)
                 IX = JX
                 IY = JY
                 DO 110,I = J + 1,N
                    IX = IX + INCX
                    IY = IY + INCY
                    Y(IY) = Y(IY) + TEMP1*A(I,J)
                    TEMP2 = TEMP2 + A(I,J)*X(IX)
  110            CONTINUE
                 Y(JY) = Y(JY) + ALPHA*TEMP2
                 JX = JX + INCX
                 JY = JY + INCY
  120         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSYMV .
*
      END 
      SUBROUTINE SSBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      CHARACTER *1 UPLO
      INTEGER N,K,LDA,INCX,INCY
      REAL ALPHA,A(LDA,*),X(*),BETA,Y(*)
*
*  Purpose
*  =======
*
*  SSBMV  performs the matrix-vector  operation 
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric band matrix, with k super-diagonals.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the band matrix A is being supplied as
*           follows: 
*
*              UPLO = 'U'          The upper triangular part of A is
*                                  being supplied.
*
*              UPLO = 'L'          The lower triangular part of A is
*                                  being supplied.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry, K specifies the number of super-diagonals of the
*           matrix A. K must satisfy  0 .le. K .lt. n.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the symmetric matrix, supplied column by 
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           Before entry with UPLO = 'L', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the symmetric matrix, supplied column by 
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the 
*           array A is not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the 
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            . 
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  Y      - REAL             array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the 
*           vector y. On exit, Y is overwritten by the updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTRINSIC MAX,MIN
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER KPLUS1,KX,KY,L
      REAL ONE,ZERO 
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP1,TEMP2
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0) .AND.
     .     (K.GE.0) .AND. (K.LT.N) .AND. (LDA.GE. (K+1))
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN 
*
*     Start the operations. In this version the elements of the array A
*     are accessed sequentially with one pass through A.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          IF (BETA.NE.ONE) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10,I = 1,N
                     Y(I) = ZERO
   10             CONTINUE
*
              ELSE
                  DO 20,I = 1,N
                     Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
*
          END IF
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
*
          IF (BETA.NE.ONE) THEN
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30,I = 1,N
                     Y(IY) = ZERO
                     IY = IY + INCY
   30             CONTINUE
*
              ELSE
                  DO 40,I = 1,N
                     Y(IY) = BETA*Y(IY) 
                     IY = IY + INCY
   40             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  y  when upper triangle of A is stored.
*
          KPLUS1 = K + 1
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 I = MAX(1,J-K)
                 DO 50,L = KPLUS1 + I - J,K
                    Y(I) = Y(I) + TEMP1*A(L,J)
                    TEMP2 = TEMP2 + A(L,J)*X(I)
                    I = I + 1 
   50            CONTINUE
                 Y(J) = Y(J) + TEMP1*A(KPLUS1,J) + ALPHA*TEMP2
   60         CONTINUE
*
          ELSE
              IX = KX - INCX
              DO 80,J = 1,N
                 TEMP1 = ALPHA*X(IX+INCX)
                 TEMP2 = ZERO 
                 IX = KX
                 IY = KY
                 DO 70,L = 1 + MAX(KPLUS1-J,0),K
                    Y(IY) = Y(IY) + TEMP1*A(L,J)
                    TEMP2 = TEMP2 + A(L,J)*X(IX)
                    IX = IX + INCX
                    IY = IY + INCY
   70            CONTINUE
                 Y(IY) = Y(IY) + TEMP1*A(KPLUS1,J) + ALPHA*TEMP2
                 IF (J.GT.K) THEN
                     KX = KX + INCX
                     KY = KY + INCY
                 END IF
*
   80         CONTINUE
          END IF
*
      ELSE
*
*        Form  y  when lower triangle of A is stored.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 100,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 Y(J) = Y(J) + TEMP1*A(1,J)
                 I = J + 1
                 DO 90,L = 2,1 + MIN(K,N-J)
                    Y(I) = Y(I) + TEMP1*A(L,J)
                    TEMP2 = TEMP2 + A(L,J)*X(I)
                    I = I + 1 
   90            CONTINUE
                 Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 120,J = 1,N
                 TEMP1 = ALPHA*X(JX)
                 TEMP2 = ZERO 
                 Y(JY) = Y(JY) + TEMP1*A(1,J)
                 IX = JX
                 IY = JY
                 DO 110,L = 2,1 + MIN(K,N-J)
                    IX = IX + INCX
                    IY = IY + INCY
                    Y(IY) = Y(IY) + TEMP1*A(L,J)
                    TEMP2 = TEMP2 + A(L,J)*X(IX)
  110            CONTINUE
                 Y(JY) = Y(JY) + ALPHA*TEMP2
                 JX = JX + INCX
                 JY = JY + INCY
  120         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSBMV .
*
      END 
      SUBROUTINE SSPMV(UPLO,N,ALPHA,AP,X,INCX,BETA,Y,INCY)
      CHARACTER *1 UPLO
      INTEGER N,INCX,INCY
      REAL ALPHA,AP(*),X(*),BETA,Y(*)
*
*  Purpose
*  =======
*
*  SSPMV  performs the matrix-vector operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows: 
*
*              UPLO = 'U'          The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L'          The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  AP     - REAL             array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with UPLO = 'U', the array AP must
*           contain the upper triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) 
*           and a( 2, 2 ) respectively, and so on.
*           Before entry with UPLO = 'L', the array AP must
*           contain the lower triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) 
*           and a( 3, 1 ) respectively, and so on.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  BETA   - REAL            . 
*           On entry, BETA specifies the scalar beta. When BETA is 
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 27-Sept-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER K,KK,KX,KY
      REAL ONE,ZERO 
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      REAL TEMP1,TEMP2
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN 
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
*     First form  y := beta*y  and set up the start points in X and Y if
*     the increments are not both unity.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          IF (BETA.NE.ONE) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10,I = 1,N
                     Y(I) = ZERO
   10             CONTINUE
*
              ELSE
                  DO 20,I = 1,N
                     Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
*
          END IF
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
*
          IF (BETA.NE.ONE) THEN
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30,I = 1,N
                     Y(IY) = ZERO
                     IY = IY + INCY
   30             CONTINUE
*
              ELSE
                  DO 40,I = 1,N
                     Y(IY) = BETA*Y(IY) 
                     IY = IY + INCY
   40             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      IF (ALPHA.EQ.ZERO) RETURN
      K = 1
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  y  when AP contains the upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 DO 50,I = 1,J - 1
                    Y(I) = Y(I) + TEMP1*AP(K)
                    TEMP2 = TEMP2 + AP(K)*X(I)
                    K = K + 1 
   50            CONTINUE
                 Y(J) = Y(J) + TEMP1*AP(K) + ALPHA*TEMP2
                 K = K + 1
   60         CONTINUE
*
          ELSE
              IX = KX - INCX
              DO 80,J = 1,N
                 TEMP1 = ALPHA*X(IX+INCX)
                 TEMP2 = ZERO 
                 IX = KX
                 IY = KY
                 KK = K
                 DO 70,K = KK,KK + J - 2
                    Y(IY) = Y(IY) + TEMP1*AP(K)
                    TEMP2 = TEMP2 + AP(K)*X(IX)
                    IX = IX + INCX
                    IY = IY + INCY
   70            CONTINUE
                 Y(IY) = Y(IY) + TEMP1*AP(K) + ALPHA*TEMP2
                 K = K + 1
   80         CONTINUE
          END IF
*
      ELSE
*
*        Form  y  when AP contains the upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 100,J = 1,N
                 TEMP1 = ALPHA*X(J)
                 TEMP2 = ZERO 
                 Y(J) = Y(J) + TEMP1*AP(K)
                 K = K + 1
                 DO 90,I = J + 1,N
                    Y(I) = Y(I) + TEMP1*AP(K)
                    TEMP2 = TEMP2 + AP(K)*X(I)
                    K = K + 1 
   90            CONTINUE
                 Y(J) = Y(J) + ALPHA*TEMP2
  100         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 120,J = 1,N
                 TEMP1 = ALPHA*X(JX)
                 TEMP2 = ZERO 
                 Y(JY) = Y(JY) + TEMP1*AP(K)
                 IX = JX
                 IY = JY
                 KK = K + 1
                 DO 110,K = KK,KK + N - (J+1)
                    IX = IX + INCX
                    IY = IY + INCY
                    Y(IY) = Y(IY) + TEMP1*AP(K)
                    TEMP2 = TEMP2 + AP(K)*X(IX)
  110            CONTINUE
                 Y(JY) = Y(JY) + ALPHA*TEMP2
                 JX = JX + INCX
                 JY = JY + INCY
  120         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSPMV .
*
      END 
      SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,LDA,INCX
      REAL A(LDA,*),X(*)
*
*  Purpose
*  =======
*
*  STRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x, 
*
*  where x is n element vector and A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U'          A is an upper triangular matrix.
*
*              UPLO = 'L'          A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows: 
*
*              TRANS = 'N' x := A*x.
*
*              TRANS = 'T' x := A'*x. 
*
*              TRANS = 'C' x := A'*x. 
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U'          A is assumed to be unit triangular.
*
*              DIAG = 'N'          A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U', the leading n by n 
*           upper triangular part of the array A must contain the upper 
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L', the leading n by n
*           lower triangular part of the array A must contain the lower 
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(n,1).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x. 
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0) .AND. (LDA.GE.N) 
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := A*x. 
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         DO 10,I = 1,J - 1
                            X(I) = X(I) + X(J)*A(I,J)
   10                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*A(J,J)
                     END IF
*
   20             CONTINUE
*
              ELSE
                  JX = KX
                  DO 40,J = 1,N
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         DO 30,I = 1,J - 1
                            X(IX) = X(IX) + X(JX)*A(I,J)
                            IX = IX + INCX
   30                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                     END IF
*
                     JX = JX + INCX
   40             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         DO 50,I = N,J + 1,-1
                            X(I) = X(I) + X(J)*A(I,J)
   50                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*A(J,J)
                     END IF
*
   60             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80,J = N,1,-1
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         DO 70,I = N,J + 1,-1
                            X(IX) = X(IX) + X(JX)*A(I,J)
                            IX = IX - INCX
   70                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                     END IF
*
                     JX = JX - INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := A'*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 100,J = N,1,-1
                     IF (NOUNIT) X(J) = X(J)*A(J,J)
                     DO 90,I = J - 1,1,-1
                        X(J) = X(J) + A(I,J)*X(I) 
   90                CONTINUE 
  100             CONTINUE
*
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 120,J = N,1,-1
                     IX = JX
                     IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                     DO 110,I = J - 1,1,-1
                        IX = IX - INCX
                        X(JX) = X(JX) + A(I,J)*X(IX)
  110                CONTINUE 
                     JX = JX - INCX
  120             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140,J = 1,N
                     IF (NOUNIT) X(J) = X(J)*A(J,J)
                     DO 130,I = J + 1,N 
                        X(J) = X(J) + A(I,J)*X(I) 
  130                CONTINUE 
  140             CONTINUE
*
              ELSE
                  JX = KX
                  DO 160,J = 1,N
                     IX = JX
                     IF (NOUNIT) X(JX) = X(JX)*A(J,J)
                     DO 150,I = J + 1,N 
                        IX = IX + INCX
                        X(JX) = X(JX) + A(I,J)*X(IX)
  150                CONTINUE 
                     JX = JX + INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STRMV .
*
      END 
      SUBROUTINE STBMV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,K,LDA,INCX
      REAL A(LDA,*),X(*)
*
*  Purpose
*  =======
*
*  STBMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x, 
*
*  where x is n element vector and A is an n by n unit, or non-unit,
*  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U' A is an upper triangular matrix.
*
*              UPLO = 'L' A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows: 
*
*              TRANS = 'N' x := A*x.
*
*              TRANS = 'T' x := A'*x. 
*
*              TRANS = 'C' x := A'*x. 
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U' A is assumed to be unit triangular.
*
*              DIAG = 'N' A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with UPLO = 'U', K specifies the number of
*           super-diagonals of the matrix A. 
*           On entry with UPLO = 'L', K specifies the number of
*           sub-diagonals of the matrix A.
*           K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           Before entry with UPLO = 'L', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the 
*           array A is not referenced.
*           Note that when DIAG = 'U' the elements of the array A
*           corresponding to the diagonal elements of the matrix are not
*           referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x. 
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 5-November-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTRINSIC MAX,MIN
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,KPLUS1,KX
      INTEGER L
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0) .AND. (K.GE.0) .AND.
     .     (LDA.GE. (K+1))
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX   too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*         Form  x := A*x.
*
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 20,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         I = MAX(1,J-K) 
                         DO 10,L = KPLUS1 + I - J,K
                            X(I) = X(I) + X(J)*A(L,J)
                            I = I + 1
   10                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*A(KPLUS1,J)
                     END IF
*
   20             CONTINUE
*
              ELSE
                  JX = KX
                  DO 40,J = 1,N
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         DO 30,L = 1 + MAX(KPLUS1-J,0),K
                            X(IX) = X(IX) + X(JX)*A(L,J)
                            IX = IX + INCX
   30                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*A(KPLUS1,J)
                     END IF
*
                     JX = JX + INCX
                     IF (J.GT.K) KX = KX + INCX
   40             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         I = MIN(N,J+K) 
                         DO 50,L = 1 + I - J,2,-1 
                            X(I) = X(I) + X(J)*A(L,J)
                            I = I - 1
   50                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*A(1,J)
                     END IF
*
   60             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80,J = N,1,-1
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         DO 70,L = 1 + MIN(K,N-J),2,-1
                            X(IX) = X(IX) + X(JX)*A(L,J)
                            IX = IX - INCX
   70                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*A(1,J)
                     END IF
*
                     JX = JX - INCX
                     IF ((N-J).GE.K) KX = KX - INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := A'*x.
*
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 100,J = N,1,-1
                     I = J
                     IF (NOUNIT) X(J) = X(J)*A(KPLUS1,J)
                     DO 90,L = K,1 + MAX(KPLUS1-J,0),-1
                        I = I - 1
                        X(J) = X(J) + A(L,J)*X(I) 
   90                CONTINUE 
  100             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 120,J = N,1,-1
                     KX = KX - INCX
                     IX = KX
                     IF (NOUNIT) X(JX) = X(JX)*A(KPLUS1,J)
                     DO 110,L = K,1 + MAX(KPLUS1-J,0),-1
                        X(JX) = X(JX) + A(L,J)*X(IX)
                        IX = IX - INCX
  110                CONTINUE 
                     JX = JX - INCX
  120             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140,J = 1,N
                     I = J
                     IF (NOUNIT) X(J) = X(J)*A(1,J)
                     DO 130,L = 2,1 + MIN(K,N-J)
                        I = I + 1
                        X(J) = X(J) + A(L,J)*X(I) 
  130                CONTINUE 
  140             CONTINUE
*
              ELSE
                  JX = KX
                  DO 160,J = 1,N
                     KX = KX + INCX
                     IX = KX
                     IF (NOUNIT) X(JX) = X(JX)*A(1,J)
                     DO 150,L = 2,1 + MIN(K,N-J)
                        X(JX) = X(JX) + A(L,J)*X(IX)
                        IX = IX + INCX
  150                CONTINUE 
                     JX = JX + INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STBMV .
*
      END 
      SUBROUTINE STPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,INCX
      REAL AP(*),X(*)
*
*  Purpose
*  =======
*
*  STPMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x, 
*
*  where x is n element vector and A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U'          A is an upper triangular matrix.
*
*              UPLO = 'L'          A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows: 
*
*              TRANS = 'N' x := A*x.
*
*              TRANS = 'T' x := A'*x. 
*
*              TRANS = 'C' x := A'*x. 
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U'          A is assumed to be unit triangular.
*
*              DIAG = 'N'          A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  AP     - REAL             array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U', the array AP must
*           contain the upper triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
*           respectively, and so on.
*           Before entry with UPLO = 'L', the array AP must
*           contain the lower triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
*           respectively, and so on.
*           Note that when  DIAG = 'U', the diagonal elements of
*           A are not referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x. 
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*  Note that UPLO, TRANS, DIAG and N must be such that the value of the
*  LOGICAL variable OK in the following statement is true.
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 2-October-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,K,KK
      INTEGER KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0)
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of AP are
*     accessed sequentially with one pass through AP.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x:= A*x.
*
          IF (LSAME(UPLO,'U')) THEN
              K = 1 
              IF (INCX.EQ.1) THEN
                  DO 20,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         DO 10,I = 1,J - 1
                            X(I) = X(I) + X(J)*AP(K)
                            K = K + 1
   10                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*AP(K)
                         K = K + 1
*
                     ELSE
                         K = K + J
                     END IF
*
   20             CONTINUE
*
              ELSE
                  JX = KX
                  DO 40,J = 1,N
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         KK = K
                         DO 30,K = KK,KK + J - 2
                            X(IX) = X(IX) + X(JX)*AP(K)
                            IX = IX + INCX
   30                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*AP(K)
                         K = K + 1
*
                     ELSE
                         K = K + J
                     END IF
*
                     JX = JX + INCX
   40             CONTINUE
              END IF
*
          ELSE
              K = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 60,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         DO 50,I = N,J + 1,-1
                            X(I) = X(I) + X(J)*AP(K)
                            K = K - 1
   50                    CONTINUE
                         IF (NOUNIT) X(J) = X(J)*AP(K)
                         K = K - 1
*
                     ELSE
                         K = K - (N-J+1)
                     END IF
*
   60             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 80,J = N,1,-1
                     IF (X(JX).NE.ZERO) THEN
                         IX = KX
                         KK = K
                         DO 70,K = KK,KK - (N- (J+1)),-1
                            X(IX) = X(IX) + X(JX)*AP(K)
                            IX = IX - INCX
   70                    CONTINUE
                         IF (NOUNIT) X(JX) = X(JX)*AP(K)
                         K = K - 1
*
                     ELSE
                         K = K - (N-J+1)
                     END IF
*
                     JX = JX - INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := A'*x.
*
          IF (LSAME(UPLO,'U')) THEN
              K = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 100,J = N,1,-1
                     IF (NOUNIT) X(J) = X(J)*AP(K)
                     K = K - 1
                     DO 90,I = J - 1,1,-1
                        X(J) = X(J) + AP(K)*X(I)
                        K = K - 1
   90                CONTINUE 
  100             CONTINUE
*
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 120,J = N,1,-1
                     IX = JX
                     IF (NOUNIT) X(JX) = X(JX)*AP(K)
                     KK = K - 1
                     DO 110,K = KK,KK - J + 2,-1
                        IX = IX - INCX
                        X(JX) = X(JX) + AP(K)*X(IX)
  110                CONTINUE 
                     JX = JX - INCX
  120             CONTINUE
              END IF
*
          ELSE
              K = 1 
              IF (INCX.EQ.1) THEN
                  DO 140,J = 1,N
                     IF (NOUNIT) X(J) = X(J)*AP(K)
                     K = K + 1
                     DO 130,I = J + 1,N 
                        X(J) = X(J) + AP(K)*X(I)
                        K = K + 1
  130                CONTINUE 
  140             CONTINUE
*
              ELSE
                  JX = KX
                  DO 160,J = 1,N
                     IX = JX
                     IF (NOUNIT) X(JX) = X(JX)*AP(K)
                     KK = K + 1
                     DO 150,K = KK,KK + N - (J+1) 
                        IX = IX + INCX
                        X(JX) = X(JX) + AP(K)*X(IX)
  150                CONTINUE 
                     JX = JX + INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STPMV .
*
      END 
      SUBROUTINE STRSV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,LDA,INCX
      REAL A(LDA,*),X(*)
*
*  Purpose
*  =======
*
*  STRSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular matrix.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U'          A is an upper triangular matrix.
*
*              UPLO = 'L'          A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows: 
*
*              TRANS = 'N' A*x = b.
*
*              TRANS = 'T' A'*x = b.
*
*              TRANS = 'C' A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U'          A is assumed to be unit triangular.
*
*              DIAG = 'N'          A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U', the leading n by n 
*           upper triangular part of the array A must contain the upper 
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L', the leading n by n
*           lower triangular part of the array A must contain the lower 
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(n,1).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0) .AND. (LDA.GE.N) 
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := inv( A )*x. 
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 20,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/A(J,J)
                         DO 10,I = J - 1,1,-1
                            X(I) = X(I) - X(J)*A(I,J)
   10                    CONTINUE
                     END IF
*
   20             CONTINUE
*
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40,J = N,1,-1
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                         IX = JX
                         DO 30,I = J - 1,1,-1
                            IX = IX - INCX
                            X(IX) = X(IX) - X(JX)*A(I,J)
   30                    CONTINUE
                     END IF
*
                     JX = JX - INCX
   40             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/A(J,J)
                         DO 50,I = J + 1,N
                            X(I) = X(I) - X(J)*A(I,J)
   50                    CONTINUE
                     END IF
*
   60             CONTINUE
*
              ELSE
                  JX = KX
                  DO 80,J = 1,N
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                         IX = JX
                         DO 70,I = J + 1,N
                            IX = IX + INCX
                            X(IX) = X(IX) - X(JX)*A(I,J)
   70                    CONTINUE
                     END IF
*
                     JX = JX + INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := inv( A' )*x.
*
          IF (LSAME(UPLO,'U')) THEN
              IF (INCX.EQ.1) THEN
                  DO 100,J = 1,N
                     DO 90,I = 1,J - 1
                        X(J) = X(J) - A(I,J)*X(I) 
   90                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/A(J,J)
  100             CONTINUE
*
              ELSE
                  JX = KX
                  DO 120,J = 1,N
                     IX = KX
                     DO 110,I = 1,J - 1 
                        X(JX) = X(JX) - A(I,J)*X(IX)
                        IX = IX + INCX
  110                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                     JX = JX + INCX
  120             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140,J = N,1,-1
                     DO 130,I = N,J + 1,-1
                        X(J) = X(J) - A(I,J)*X(I) 
  130                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/A(J,J)
  140             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160,J = N,1,-1
                     IX = KX
                     DO 150,I = N,J + 1,-1
                        X(JX) = X(JX) - A(I,J)*X(IX)
                        IX = IX - INCX
  150                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/A(J,J)
                     JX = JX - INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STRSV .
*
      END 
      SUBROUTINE STBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,K,LDA,INCX
      REAL A(LDA,*),X(*)
*
*  Purpose
*  =======
*
*  STBSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular band matrix, with ( k + 1 )
*  diagonals.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U'          A is an upper triangular matrix.
*
*              UPLO = 'L'          A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows: 
*
*              TRANS = 'N' A*x = b.
*
*              TRANS = 'T' A'*x = b.
*
*              TRANS = 'C' A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U' A is assumed to be unit triangular.
*
*              DIAG = 'N' A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with UPLO = 'U', K specifies the number of
*           super-diagonals of the matrix A. 
*           On entry with UPLO = 'L', K specifies the number of
*           sub-diagonals of the matrix A.
*           K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           Before entry with UPLO = 'L', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the 
*           array A is not referenced.
*           Note that when DIAG = 'U' the elements of the array A
*           corresponding to the diagonal elements of the matrix are not
*           referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the leading dimension of A as
*           declared in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 7-November-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTRINSIC MAX,MIN
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,KPLUS1,KX
      INTEGER L
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0) .AND. (K.GE.0) .AND.
     .     (LDA.GE. (K+1))
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed by sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := inv( A )*x. 
*
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 20,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                         I = J
                         DO 10,L = K,1 + MAX(KPLUS1-J,0),-1 
                            I = I - 1
                            X(I) = X(I) - X(J)*A(L,J)
   10                    CONTINUE
                     END IF
*
   20             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40,J = N,1,-1
                     KX = KX - INCX
                     IX = KX
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                         DO 30 L = K,1 + MAX(KPLUS1-J,0),-1 
                            X(IX) = X(IX) - X(JX)*A(L,J)
                            IX = IX - INCX
   30                    CONTINUE
                     END IF
*
                     JX = JX - INCX
   40             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/A(1,J)
                         I = J
                         DO 50,L = 2,1 + MIN(K,N-J)
                            I = I + 1
                            X(I) = X(I) - X(J)*A(L,J)
   50                    CONTINUE
                     END IF
*
   60             CONTINUE
*
              ELSE
                  JX = KX
                  DO 80,J = 1,N
                     KX = KX + INCX
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                         IX = KX
                         DO 70,L = 2,1 + MIN(K,N-J)
                            X(IX) = X(IX) - X(JX)*A(L,J)
                            IX = IX + INCX
   70                    CONTINUE
                     END IF
*
                     JX = JX + INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := inv( A')*x. 
*
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
                  DO 100,J = 1,N
                     I = MAX(1,J-K)
                     DO 90,L = KPLUS1 + I - J,K
                        X(J) = X(J) - A(L,J)*X(I) 
                        I = I + 1
   90                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
  100             CONTINUE
*
              ELSE
                  JX = KX
                  DO 120,J = 1,N
                     IX = KX
                     DO 110,L = 1 + MAX(KPLUS1-J,0),K
                        X(JX) = X(JX) - A(L,J)*X(IX)
                        IX = IX + INCX
  110                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                     JX = JX + INCX
                     IF (J.GT.K) KX = KX + INCX
  120             CONTINUE
              END IF
*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 140,J = N,1,-1
                     I = MIN(N,J+K)
                     DO 130,L = 1 + I - J,2,-1
                        X(J) = X(J) - A(L,J)*X(I) 
                        I = I - 1
  130                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/A(1,J)
  140             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160,J = N,1,-1
                     IX = KX
                     DO 150,L = 1 + MIN(K,N-J),2,-1
                        X(JX) = X(JX) - A(L,J)*X(IX)
                        IX = IX - INCX
  150                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                     JX = JX - INCX
                     IF ((N-J).GE.K) KX = KX - INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STBSV .
*
      END 
      SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
      CHARACTER *1 UPLO,TRANS,DIAG
      INTEGER N,INCX
      REAL AP(*),X(*)
*
*  Purpose
*  =======
*
*  STPSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular matrix.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows: 
*
*              UPLO = 'U' A is an upper triangular matrix.
*
*              UPLO = 'L' A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows: 
*
*              TRANS = 'N' A*x = b.
*
*              TRANS = 'T' A'*x = b.
*
*              TRANS = 'C' A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows: 
*
*              DIAG = 'U' A is assumed to be unit triangular.
*
*              DIAG = 'N' A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  AP     - REAL             array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U', the array AP must
*           contain the upper triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
*           respectively, and so on.
*           Before entry with UPLO = 'L', the array AP must
*           contain the lower triangular matrix packed sequentially,
*           column by column, so that AP( 1 ) contains a( 1, 1 ),
*           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
*           respectively, and so on.
*           Note that when  DIAG = 'U', the diagonal elements of
*           A are not referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 11-November-1985. 
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      LOGICAL NOUNIT
      INTEGER I,IX,J,JX,K,KK
      INTEGER KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND.
     .     (LSAME(TRANS,'N') .OR. LSAME(TRANS,'T') .OR.
     .     LSAME(TRANS,'C')) .AND. (LSAME(DIAG,'U') .OR.
     .     LSAME(DIAG,'N')) .AND. (N.GT.0)
*
*     Quick return if possible.
*
      IF ( .NOT. OK) RETURN
      NOUNIT = LSAME(DIAG,'N')
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of AP are
*     accessed sequentially with one pass through AP.
*
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  x := inv( A )*x. 
*
          IF (LSAME(UPLO,'U')) THEN
              K = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 20,J = N,1,-1
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/AP(K)
                         K = K - 1
                         DO 10,I = J - 1,1,-1
                            X(I) = X(I) - X(J)*AP(K)
                            K = K - 1
   10                    CONTINUE
*
                     ELSE
                         K = K - J
                     END IF
*
   20             CONTINUE
*
              ELSE
                  JX = KX + (N-1)*INCX
                  DO 40,J = N,1,-1
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/AP(K)
                         IX = JX
                         KK = K - 1
                         DO 30,K = KK,KK - J + 2,-1
                            IX = IX - INCX
                            X(IX) = X(IX) - X(JX)*AP(K)
   30                    CONTINUE
*
                     ELSE
                         K = K - J
                     END IF
*
                     JX = JX - INCX
   40             CONTINUE
              END IF
*
          ELSE
              K = 1 
              IF (INCX.EQ.1) THEN
                  DO 60,J = 1,N
                     IF (X(J).NE.ZERO) THEN
                         IF (NOUNIT) X(J) = X(J)/AP(K)
                         K = K + 1
                         DO 50,I = J + 1,N
                            X(I) = X(I) - X(J)*AP(K)
                            K = K + 1
   50                    CONTINUE
*
                     ELSE
                         K = K + N - J + 1
                     END IF
*
   60             CONTINUE
*
              ELSE
                  JX = KX
                  DO 80,J = 1,N
                     IF (X(JX).NE.ZERO) THEN
                         IF (NOUNIT) X(JX) = X(JX)/AP(K)
                         IX = JX
                         KK = K + 1
                         DO 70,K = KK,KK + N - (J+1)
                            IX = IX + INCX
                            X(IX) = X(IX) - X(JX)*AP(K)
   70                    CONTINUE
*
                     ELSE
                         K = K + N - J + 1
                     END IF
*
                     JX = JX + INCX
   80             CONTINUE
              END IF
*
          END IF
*
      ELSE
*
*        Form  x := inv( A' )*x.
*
          IF (LSAME(UPLO,'U')) THEN
              K = 1 
              IF (INCX.EQ.1) THEN
                  DO 100,J = 1,N
                     DO 90,I = 1,J - 1
                        X(J) = X(J) - AP(K)*X(I)
                        K = K + 1
   90                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/AP(K)
                     K = K + 1
  100             CONTINUE
*
              ELSE
                  JX = KX
                  DO 120,J = 1,N
                     IX = KX
                     KK = K
                     DO 110,K = KK,KK + J - 2
                        X(JX) = X(JX) - AP(K)*X(IX)
                        IX = IX + INCX
  110                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/AP(K)
                     K = K + 1
                     JX = JX + INCX
  120             CONTINUE
              END IF
*
          ELSE
              K = (N* (N+1))/2
              IF (INCX.EQ.1) THEN
                  DO 140,J = N,1,-1
                     DO 130,I = N,J + 1,-1
                        X(J) = X(J) - AP(K)*X(I)
                        K = K - 1
  130                CONTINUE 
                     IF (NOUNIT) X(J) = X(J)/AP(K)
                     K = K - 1
  140             CONTINUE
*
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160,J = N,1,-1
                     IX = KX
                     KK = K
                     DO 150,K = KK,KK - (N- (J+1)),-1
                        X(JX) = X(JX) - AP(K)*X(IX)
                        IX = IX - INCX
  150                CONTINUE 
                     IF (NOUNIT) X(JX) = X(JX)/AP(K)
                     K = K - 1
                     JX = JX - INCX
  160             CONTINUE
              END IF
*
          END IF
*
      END IF
*
      RETURN
*
*     End of STPSV .
*
      END 
      SUBROUTINE SGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      INTEGER M,N,INCX,INCY,LDA
      REAL ALPHA,X(*),Y(*),A(LDA,*)
*
*  Purpose
*  =======
*
*  SGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix. 
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(1,m).
*           Unchanged on exit.
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-August-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,J,JY,KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      REAL TEMP
      LOGICAL OK
      OK = (M.GT.0) .AND. (N.GT.0) .AND. (LDA.GE.M)
*
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN
          DO 20,J = 1,N
             IF (Y(J).NE.ZERO) THEN
                 TEMP = ALPHA*Y(J)
                 DO 10,I = 1,M
                    A(I,J) = A(I,J) + X(I)*TEMP
   10            CONTINUE
             END IF 
*
   20     CONTINUE
*
      ELSE
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (M-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              JY = 1
*
          ELSE
              JY = 1 - (N-1)*INCY
          END IF
*
          DO 40,J = 1,N
             IF (Y(JY).NE.ZERO) THEN
                 TEMP = ALPHA*Y(JY)
                 IX = KX
                 DO 30,I = 1,M
                    A(I,J) = A(I,J) + X(IX)*TEMP
                    IX = IX + INCX
   30            CONTINUE
             END IF 
*
             JY = JY + INCY
   40     CONTINUE
      END IF
*
      RETURN
*
*     End of SGER  .
*
      END 
      SUBROUTINE SSYR(UPLO,N,ALPHA,X,INCX,A,LDA)
      CHARACTER *1 UPLO
      INTEGER N,INCX,LDA
      REAL ALPHA,X(*),A(LDA,*)
*
*  Purpose
*  =======
*
*  SSYR   performs the symmetric rank 1 operation
*
*     A := alpha*x*x' + A,
*
*  where alpha is a real scalar, x is an n element vector and A is an
*  n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows: 
*
*              UPLO = 'U' Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U', the leading n by n 
*           upper triangular part of the array A must contain the upper 
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L', the leading n by n
*           lower triangular part of the array A must contain the lower 
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(1,n).
*           Unchanged on exit.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 27-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,J,JX,KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      REAL TEMP
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0) .AND.
     .     (LDA.GE.N)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when A is stored in upper triangle.
*
          IF (INCX.EQ.1) THEN 
              DO 20,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     DO 10,I = 1,J
                        A(I,J) = A(I,J) + X(I)*TEMP
   10                CONTINUE 
                 END IF
*
   20         CONTINUE
*
          ELSE
              JX = KX
              DO 40,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IX = KX
                     DO 30,I = 1,J
                        A(I,J) = A(I,J) + X(IX)*TEMP
                        IX = IX + INCX
   30                CONTINUE 
                 END IF
*
                 JX = JX + INCX
   40         CONTINUE
          END IF
*
      ELSE
*
*        Form  A  when A is stored in lower triangle.
*
          IF (INCX.EQ.1) THEN 
              DO 60,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     DO 50,I = J,N
                        A(I,J) = A(I,J) + X(I)*TEMP
   50                CONTINUE 
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              DO 80,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IX = JX
                     DO 70,I = J,N
                        A(I,J) = A(I,J) + X(IX)*TEMP
                        IX = IX + INCX
   70                CONTINUE 
                 END IF
*
                 JX = JX + INCX
   80         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSYR  .
*
      END 
      SUBROUTINE SSPR(UPLO,N,ALPHA,X,INCX,AP)
      CHARACTER *1 UPLO
      INTEGER N,INCX
      REAL ALPHA,X(*),AP(*)
*
*  Purpose
*  =======
*
*  SSPR    performs the symmetric rank 1 operation
*
*     A := alpha*x*x' + A,
*
*  where alpha is a real scalar, x is an n element vector and A is an
*  n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows: 
*
*              UPLO = 'U' The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L' The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  AP     - REAL             array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U', the array AP must
*           contain the upper triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) 
*           and a( 2, 2 ) respectively, and so on. On exit, the array
*           AP is overwritten by the upper triangular part of the
*           updated matrix.
*           Before entry with UPLO = 'L', the array AP must
*           contain the lower triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) 
*           and a( 3, 1 ) respectively, and so on. On exit, the array
*           AP is overwritten by the lower triangular part of the
*           updated matrix.
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,J,JX,K,KK
      INTEGER KX
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      REAL TEMP
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set the start point in X if the increment is not unity.
*
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX 
*
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
      K = 1
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when upper triangle is stored in AP.
*
          IF (INCX.EQ.1) THEN 
              DO 20,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     DO 10,I = 1,J
                        AP(K) = AP(K) + X(I)*TEMP 
                        K = K + 1
   10                CONTINUE 
*
                 ELSE
                     K = K + J
                 END IF
*
   20         CONTINUE
*
          ELSE
              JX = KX
              DO 40,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IX = KX
                     KK = K
                     DO 30,K = KK,KK + J - 1
                        AP(K) = AP(K) + X(IX)*TEMP
                        IX = IX + INCX
   30                CONTINUE 
*
                 ELSE
                     K = K + J
                 END IF
*
                 JX = JX + INCX
   40         CONTINUE
          END IF
*
      ELSE
*
*        Form  A  when lower triangle is stored in AP.
*
          IF (INCX.EQ.1) THEN 
              DO 60,J = 1,N
                 IF (X(J).NE.ZERO) THEN 
                     TEMP = ALPHA*X(J)
                     DO 50,I = J,N
                        AP(K) = AP(K) + X(I)*TEMP 
                        K = K + 1
   50                CONTINUE 
*
                 ELSE
                     K = K + N - J + 1
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              DO 80,J = 1,N
                 IF (X(JX).NE.ZERO) THEN
                     TEMP = ALPHA*X(JX) 
                     IX = JX
                     KK = K
                     DO 70,K = KK,KK + N - J
                        AP(K) = AP(K) + X(IX)*TEMP
                        IX = IX + INCX
   70                CONTINUE 
*
                 ELSE
                     K = K + N - J + 1
                 END IF
*
                 JX = JX + INCX
   80         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSPR  .
*
      END 
      SUBROUTINE SSYR2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
      CHARACTER *1 UPLO
      INTEGER N,INCX,INCY,LDA 
      REAL ALPHA,X(*),Y(*),A(LDA,*)
*
*  Purpose
*  =======
*
*  SSYR2  performs the symmetric rank 2 operation
*
*     A := alpha*x*y' + alpha*y*x' + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an n
*  by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows: 
*
*              UPLO = 'U' Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*  A      - REAL             array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U', the leading n by n 
*           upper triangular part of the array A must contain the upper 
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L', the leading n by n
*           lower triangular part of the array A must contain the lower 
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least max(1,n).
*           Unchanged on exit.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 27-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER KX,KY 
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      REAL TEMP1,TEMP2
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0) .AND.
     .     (LDA.GE.N)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set up the start points in X and Y if the increments are not both 
*     unity.
*
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
*
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when A is stored in the upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 20,J = 1,N
                 IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(J) 
                     TEMP2 = ALPHA*X(J) 
                     DO 10,I = 1,J
                        A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   10                CONTINUE 
                 END IF
*
   20         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 40,J = 1,N
                 IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(JY)
                     TEMP2 = ALPHA*X(JX)
                     IX = KX
                     IY = KY
                     DO 30,I = 1,J
                        A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                        IX = IX + INCX
                        IY = IY + INCY
   30                CONTINUE 
                 END IF
*
                 JX = JX + INCX
                 JY = JY + INCY
   40         CONTINUE
          END IF
*
      ELSE
*
*        Form  A  when A is stored in the upper triangle.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(J) 
                     TEMP2 = ALPHA*X(J) 
                     DO 50,I = J,N
                        A(I,J) = A(I,J) + X(I)*TEMP1 + Y(I)*TEMP2
   50                CONTINUE 
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 80,J = 1,N
                 IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(JY)
                     TEMP2 = ALPHA*X(JX)
                     IX = JX
                     IY = JY
                     DO 70,I = J,N
                        A(I,J) = A(I,J) + X(IX)*TEMP1 + Y(IY)*TEMP2
                        IX = IX + INCX
                        IY = IY + INCY
   70                CONTINUE 
                 END IF
*
                 JX = JX + INCX
                 JY = JY + INCY
   80         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSYR2 .
*
      END 
      SUBROUTINE SSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
      CHARACTER *1 UPLO
      INTEGER N,INCX,INCY
      REAL ALPHA,X(*),Y(*),AP(*)
*
*  Purpose
*  =======
*
*  SSPR2  performs the symmetric rank 2 operation
*
*     A := alpha*x*y' + alpha*y*x' + A,
*
*  where alpha is a scalar, x and y are n element vectors and A is an
*  n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the matrix A is supplied in the packed
*           array AP as follows: 
*
*              UPLO = 'U' The upper triangular part of A is
*                                  supplied in AP.
*
*              UPLO = 'L' The lower triangular part of A is
*                                  supplied in AP.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - REAL            . 
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X.
*           Unchanged on exit.
*
*  Y      - REAL             array of dimension at least 
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y.
*           Unchanged on exit.
*
*  AP     - REAL             array of DIMENSION at least
*           ( ( n*( n + 1 ) )/2 ).
*           Before entry with  UPLO = 'U', the array AP must
*           contain the upper triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 ) 
*           and a( 2, 2 ) respectively, and so on. On exit, the array
*           AP is overwritten by the upper triangular part of the
*           updated matrix.
*           Before entry with UPLO = 'L', the array AP must
*           contain the lower triangular part of the symmetric matrix
*           packed sequentially, column by column, so that AP( 1 )
*           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 ) 
*           and a( 3, 1 ) respectively, and so on. On exit, the array
*           AP is overwritten by the lower triangular part of the
*           updated matrix.
*
*
*
*
*
*  Level 2 Blas routine.
*
*  -- Written on 30-September-1985.
*     Sven Hammarling, Nag Central Office.
C     REVISED 860623
C     REVISED YYMMDD
C     BY R. J. HANSON, SANDIA NATIONAL LABS.
*
      INTEGER I,IX,IY,J,JX,JY 
      INTEGER K,KK,KX,KY
      REAL ZERO
      PARAMETER (ZERO=0.0E+0) 
      REAL TEMP1,TEMP2
      LOGICAL OK,LSAME
      OK = (LSAME(UPLO,'U') .OR. LSAME(UPLO,'L')) .AND. (N.GT.0)
*
*     Quick return if possible.
*
      IF ( .NOT. OK .OR. (ALPHA.EQ.ZERO)) RETURN
*
*     Set up the start points in X and Y if the increments are not both 
*     unity.
*
      IF ((INCX.NE.1) .OR. (INCY.NE.1)) THEN
          IF (INCX.GT.0) THEN 
              KX = 1
*
          ELSE
              KX = 1 - (N-1)*INCX
          END IF
*
          IF (INCY.GT.0) THEN 
              KY = 1
*
          ELSE
              KY = 1 - (N-1)*INCY
          END IF
*
      END IF
*
*     Start the operations. In this version the elements of the array AP
*     are accessed sequentially with one pass through AP.
*
      K = 1
      IF (LSAME(UPLO,'U')) THEN
*
*        Form  A  when upper triangle is stored in AP.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 20,J = 1,N
                 IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(J) 
                     TEMP2 = ALPHA*X(J) 
                     DO 10,I = 1,J
                        AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                        K = K + 1
   10                CONTINUE 
*
                 ELSE
                     K = K + J
                 END IF
*
   20         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 40,J = 1,N
                 IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(JY)
                     TEMP2 = ALPHA*X(JX)
                     IX = KX
                     IY = KY
                     KK = K
                     DO 30,K = KK,KK + J - 1
                        AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                        IX = IX + INCX
                        IY = IY + INCY
   30                CONTINUE 
*
                 ELSE
                     K = K + J
                 END IF
*
                 JX = JX + INCX
                 JY = JY + INCY
   40         CONTINUE
          END IF
*
      ELSE
*
*        Form  A  when lower triangle is stored in AP.
*
          IF ((INCX.EQ.1) .AND. (INCY.EQ.1)) THEN 
              DO 60,J = 1,N
                 IF ((X(J).NE.ZERO) .OR. (Y(J).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(J) 
                     TEMP2 = ALPHA*X(J) 
                     DO 50,I = J,N
                        AP(K) = AP(K) + X(I)*TEMP1 + Y(I)*TEMP2
                        K = K + 1
   50                CONTINUE 
*
                 ELSE
                     K = K + N - J + 1
                 END IF
*
   60         CONTINUE
*
          ELSE
              JX = KX
              JY = KY
              DO 80,J = 1,N
                 IF ((X(JX).NE.ZERO) .OR. (Y(JY).NE.ZERO)) THEN
                     TEMP1 = ALPHA*Y(JY)
                     TEMP2 = ALPHA*X(JX)
                     IX = JX
                     IY = JY
                     KK = K
                     DO 70,K = KK,KK + N - J
                        AP(K) = AP(K) + X(IX)*TEMP1 + Y(IY)*TEMP2
                        IX = IX + INCX
                        IY = IY + INCY
   70                CONTINUE 
*
                 ELSE
                     K = K + N - J + 1
                 END IF
*
                 JX = JX + INCX
                 JY = JY + INCY
   80         CONTINUE
          END IF
*
      END IF
*
      RETURN
*
*     End of SSPR2 .
*
      END 
C      LOGICAL FUNCTION LSAME(CA,CB)
C     TEST IF TWO CHARACTERS ARE ESSENTIALLY THE SAME.
C     THE CHARACTER CB IS ONE OF THE FORTRAN SET. 
C     (LOWER AND UPPER CASE LETTERS ARE EQUIVALENT.)
C     THIS IS A SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
C
C     THIS SUBPROGRAM IS MACHINE-DEPENDENT.
C     VERSION FOR CDC SYSTEMS USING 6-12 BIT REPRESENTATIONS.
C      CHARACTER CA(*)
C      CHARACTER *1 CB
C      INTEGER ICIRFX
C      DATA ICIRFX/62/
C     SEE IF THE FIRST CHAR. IN STRING CA EQUALS STRING CB. 
C      LSAME = CA(1) .EQ. CB .AND. CA(1) .NE. CHAR(ICIRFX)
C      IF (LSAME) RETURN
C     THE CHARS. ARE NOT IDENTICAL.  NOW CHECK THEM FOR EQUIVALENCE.
C     LOOK FOR THE 'ESCAPE' CHARACTER, CIRCUMFLEX, FOLLOWED BY
C     THE LETTER.
C      IVAL = ICHAR(CA(2))
C      IF (IVAL.GE.ICHAR('A') .AND. IVAL.LE.ICHAR('Z')) THEN 
C          LSAME = CA(1) .EQ. CHAR(ICIRFX) .AND. CA(2) .EQ. CB
C      END IF
C*
C      RETURN
C     END 
      LOGICAL FUNCTION LSAME(CA,CB)
C     TEST IF TWO CHARACTERS ARE ESSENTIALLY THE SAME.
C     THE CHARACTER CB IS ONE OF THE FORTRAN SET. 
C     (LOWER AND UPPER CASE LETTERS ARE EQUIVALENT.)
C     THIS IS A SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
C
C     THIS SUBPROGRAM IS MACHINE-DEPENDENT.
C     VERSION FOR ANY ASCII MACHINE.
      CHARACTER *1 CA
      CHARACTER *1 CB
      INTEGER IOFF
      DATA IOFF/32/ 
C     SEE IF THE  CHAR. IN STRING CA EQUALS STRING CB.
      LSAME = CA .EQ. CB
      IF (LSAME) RETURN
C      THE CHARS. ARE NOT IDENTICAL.  NOW CHECK THEM FOR EQUIVALENCE.
      ISHIFT = ICHAR(CA) - IOFF
      IF (ISHIFT.GE.ICHAR('A') .AND. ISHIFT.LE.ICHAR('Z')) THEN
          LSAME = ISHIFT .EQ. ICHAR(CB) 
      END IF

      RETURN
      END 
C
C     LOGICAL FUNCTION LSAME(CA,CB)
C     TEST IF TWO CHARACTERS ARE ESSENTIALLY THE SAME.
C     THE CHARACTER CB IS ONE OF THE FORTRAN SET. 
C     (LOWER AND UPPER CASE LETTERS ARE EQUIVALENT.)
C     THIS IS A SUBPROGRAM FOR THE LEVEL TWO BLAS.
C     REVISED 860623
C     REVISED YYMMDD
C     AUTH=R. J. HANSON, SANDIA NATIONAL LABS.
C
C     THIS SUBPROGRAM IS MACHINE-DEPENDENT.
C     VERSION FOR ANY EBCDIC MACHINE.
C     CHARACTER *1 CA
C     CHARACTER *1 CB
C     INTEGER IOFF
C     DATA IOFF/64/ 
C     SEE IF THE  CHAR. IN STRING CA EQUALS STRING CB.
C     LSAME = CA .EQ. CB
C     IF (LSAME) RETURN
C     THE CHARS. ARE NOT IDENTICAL.  NOW CHECK THEM FOR EQUIVALENCE.
C     ISHIFT = ICHAR(CA) + IOFF
C     IF (ISHIFT.GE.ICHAR('A') .AND. ISHIFT.LE.ICHAR('I')) THEN
C         LSAME = ISHIFT .EQ. ICHAR(CB) 
C     END IF
C
C     IF (ISHIFT.GE.ICHAR('J') .AND. ISHIFT.LE.ICHAR('R')) THEN
C         LSAME = ISHIFT .EQ. ICHAR(CB) 
C     END IF
C
C     IF (ISHIFT.GE.ICHAR('S') .AND. ISHIFT.LE.ICHAR('Z')) THEN
C         LSAME = ISHIFT .EQ. ICHAR(CB) 
C     END IF
C
C     RETURN
C     END 
 
