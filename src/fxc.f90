
!**************************************************************

FUNCTION FXC(N)
      IMPLICIT REAL*8 (A-Z)
      LOGICAL*4 FRSTME
      DATA FRSTME /.TRUE./
      SAVE

      FRSTME =.TRUE.
      If (FRSTME) then
         FRSTME=.FALSE.
         ONTHRD=1.D0/3.D0
         PI=4.D0*DATAN(1.D0)
         ALPHA=(4.D0/(9.D0*PI))**ONTHRD
         RSMUX=-2.D0/(PI*ALPHA)
         VOLFAC=3.D0/(4.D0*PI)
         ! PERDEW-ZUNGER PARAMETERIZATION OF CEPERLEY-ALDER PARAMAGNETIC
         ! CORRELATION POTENTIAL (PHYS. REV. B 23, 5048 (1981))
         ! NOTE FACTOR OF 2 FOR RY UNITS
         GAMMA=2.D0*(-.1423D0)
         BETA1=1.0529D0
         BETA2=.3334D0
         A=2.D0*.0311D0
         B=2.D0*(-.048D0)
         C=2.D0*.0020D0
         D=2.D0*(-.0116D0)
         XBET1=(7.D0/6.D0)*BETA1
         XBET2=(4.D0/3.D0)*BETA2
         XBA=B-A/3.D0
         XC=(2.D0/3.D0)*C
         XDC=(2.D0*D-C)/3.D0
      End If

      ! IN CASE FUNCTION IS GIVEN A NEGATIVE DENSITY:
      NN=DMAX1(N,1.D-30)
      RS=(VOLFAC/NN)**ONTHRD
      IF(RS.GE.1.D0) THEN
         SQRTRS=DSQRT(RS)
         Y=(1.D0+BETA1*SQRTRS+BETA2*RS)**2
         MUC=GAMMA*(1.D0+XBET1*SQRTRS+XBET2*RS)/Y
      ELSE
         LNRS=DLOG(RS)
         MUC=A*LNRS+XBA+XC*RS*LNRS+XDC*RS
      END IF
      MUX=RSMUX/RS
      FXC=MUX+MUC

      RETURN
END FUNCTION
