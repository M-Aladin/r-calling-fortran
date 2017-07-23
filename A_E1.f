      REAL*8 FUNCTION A_E1(X)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     FUNCTION: A_E1 - EXPONENTIAL INTEGRAL E1(X)                      C
C     DATE CREATED:  1986                                              C
C     LAST MODIFIED: 11 AUGUST 1996 (CHECKED BY MIN-YU-SHIH)           C
C     AUTHOR: THOMAS A. BLASINGAME -- TEXAS A&M UNIVERSITY             C
C     SOURCE- EQNS. 5.1.53 AND 5.1.56, P.231 ABRAMOWITZ AND STEGUN     C
C             “HANDBOOK OF MATHEMATICAL FUNCTIONS                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        1         2         3         4         5         6         7
C23455789012345678901234567890123456789012345578901234567890123456789012
C--DECLARES ALL A-H,O-Z REAL*8
      IMPLICIT REAL*8 (A-H,O-Z) ! DECLARES ALL A-H,O-Z REAL*8
C--DECLARE5 ALL I-N INTEGER*4
      IMPLICIT INTEGER*4(I-N)   ! DECLARES ALL I-N INTEGER*4
C----------------------------------------------------------------------C
C   EXPONENTIAL INTEGRAL APPROXIMATION (EQNS. 5.1.53 AND 5.1.56)       C
C----------------------------------------------------------------------C
c--CONDITIONALS
C----------------------------------------------------------------------C
C     [0 < X < 1] USING EQ. 5.1.53 (ERROR NORM < 2.E—7)                C
C----------------------------------------------------------------------C
      IF (X .LE. 1.D0) THEN
         T1 = - 0.57721556D0        + 0.99999193D0*X
         T2 = - 0.24951055D0*(X**2) + 0.05519968D0*(x**3)
         T3 = - 0.00976004D0*(X**4) + 0.00107857D0*(X**5)
         T4 = - DLOG(X)
         A_E1 = T1+T2+T3+T4
      ELSEIF (X .GT. 1.D0 .AND. X .LT. 1.E4) THEN
C----------------------------------------------------------------------C
C     [1 < X < INFINITY) USING EQ. 5.1.56 (ERROR NORM < 2.E—8)         C
C----------------------------------------------------------------------C
         T1 =  1.0000000000D0*(X**4)
         T2 =  8.5733237401D0*(X**3)
         T3 = 18.0590169730D0*(X**2)
         T4 =  8.6347608925D0*(X**1)
         T5 =  0.2677737343D0
         B1   =  1.0000000000D0*(X**4)
         B2   =  9.5733223454D0*(X**3)
         B3   = 25.6329561486D0*(X**2)
         B4   = 21.0996530827D0*(X**1)
         B5   =  3.9584969228D0
         A_E1 = (1.D0/(X*DEXP(X)))*
     &           ((T1+T2+T3+T4+T5)/(B1+B2+B3+B4+B5))
       ELSE
          A_E1 = 0.D0  ! SETS E1(X>1.E4) = 0, AVOIDS OVERFLOW IN EXP(X)
       ENDIF
C        1         2         3         4         5         6         7
C2345578901234567890123456789012345678901234557890123£567890123456789012
      RETURN
      END
      
C        1         2         3         4         5         6         7
C23456789012345678901234567390123456785012345673901234567850123456789012
      SUBROUTINE A_E1_SUB(n, answer)
         REAL*8 n, answer, A_E1
         EXTERNAL A_E1
         answer = A_E1(n)
      END SUBROUTINE          