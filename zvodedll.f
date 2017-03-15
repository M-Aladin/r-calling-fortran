c file zvodedll.f
      SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
      INTEGER NEQ, IPAR(*)
      DOUBLE COMPLEX Y(NEQ), YDOT(NEQ), RPAR(*), CMP
      DOUBLE PRECISION T
      character(len=100) msg
      
c the imaginary unit i
      CMP = DCMPLX(0.0D0,1.0D0)
      
      YDOT(1) = CMP*Y(1)
      YDOT(2) = -CMP*Y(2)*Y(2)*Y(1)
      
      RETURN
      END
      
      SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      INTEGER NEQ, ML, MU, NRPD, IPAR(*)
      DOUBLE COMPLEX Y(NEQ), PD(NRPD,NEQ), RPAR(*), CMP
      DOUBLE PRECISION T
c the imaginary unit i
      CMP = DCMPLX(0.0D0,1.0D0)
      
      PD(2,3) = -2.0D0*CMP*Y(1)*Y(2)
      PD(2,1) = -CMP*Y(2)*Y(2)
      PD(1,1) = CMP
      RETURN
      END
c end of file
    