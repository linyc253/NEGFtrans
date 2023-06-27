!******************************************************************************
!
! File : anderson_mixing.f90
!   by : Alan Tackett
!   on : 07/19/99
!  for : Potential and Wij Mixing
!
!  This module contains routines to implement the extended Anderson
!  mixing method as outlined in V. Eyert, J. Comp. Phys. 124,  271(1996).
!
!  The method is defined by eq's: 2.1, 8.2, 7.7
!
!******************************************************************************

Module anderson_mixing


Implicit NONE!!!!!!

Type Anderson_Context  !** Anderson Mixing context
    Real*8    :: w0        !** w0^2 used in 8.2 to ensure linear independence 
    Real*8    :: NewMix    !** Amount of new vectors to mix, ie beta in paper.
    Integer :: Nmax      !** Max number of vectors to keep 
    Integer :: N         !** Current number of vectors in list
    Integer :: Slot      !** New fill Slot
    Integer :: VecSize   !** Size of each vector
    Integer :: Err_Unit  !** Error unit
    Integer :: Lwork
    Integer, Pointer :: IPIV(:)
    Complex*16, Pointer :: Work(:)
    Complex*16, Pointer :: Matrix(:,:)
    Complex*16, Pointer :: DupMatrix(:,:)
    Complex*16, Pointer :: Gamma(:)  !** Gamma as defined in 7.6
    Complex*16, Pointer :: DF(:,:)   !** Delta F
    Complex*16, Pointer :: Fprev(:)
    Complex*16, Pointer :: DX(:,:)
    Complex*16, Pointer :: Xprev(:)
    Complex*16, Pointer :: Xcurr(:)
End Type
    

!******************************************************************************
Contains
!******************************************************************************


!******************************************************************************
!
! Anderson_Mix - Performs the actual mixing of the input vector with the
!                history and retuns the result.
!
!   AC - Anderson context
!   X  - Current vector on input and new quess on output
!   F  - F(X) - X. Nonlinear mixing of input vector
!
!******************************************************************************

Subroutine Anderson_Mix(AC, X, F)
  Type  (Anderson_Context), Intent(INOUT) :: AC
  Complex*16,                  Intent(INOUT) :: X(:)
  Complex*16,                  Intent(IN)    :: F(:)

  Integer :: i, slot
  Complex*16 :: term
  Real*8    :: tmp

  !** First determine where to store the new correction vectors ***
  AC%slot = AC%slot + 1
  If (AC%Slot>AC%Nmax) AC%Slot = 1

  If ((AC%N < 0) .OR. (AC%Nmax == 0)) then  !** Simple mixing for 1st time ***
     AC%Xprev = X
     X = X + AC%NewMix*F
  else
     slot = AC%Slot
!Write(AC%Err_Unit,*) 'And_Mix: Slot=',Slot, ' * N=',AC%N

     AC%DF(:,slot) = F - AC%Fprev   !** Make new DF vector
     AC%DX(:,slot) = X - AC%Xprev   !** Make new DX vector

     Do i=1, Min(AC%N+1,AC%Nmax)     !*** Add row/col to matrix
        term = DOT_PRODUCT(AC%DF(:,i), AC%DF(:,slot))
        AC%Matrix(i,slot) = term
        if (i /= slot) AC%Matrix(slot,i) = CONJG(term)

        AC%Gamma(i) = DOT_PRODUCT(AC%DF(:,i), F)
     End Do

     AC%Matrix(slot,slot) = (1+AC%w0) * AC%Matrix(slot,slot)

     AC%DupMatrix = AC%Matrix

     slot = Min(AC%N+1,AC%Nmax) 
     Call ZGESV(Slot, 1, AC%DupMatrix(1,1), AC%Nmax, AC%IPIV(1), &
          AC%Gamma(1), AC%Nmax, i)

     If (i /= 0) then
        Write(AC%Err_Unit,*) 'Anderson_Mix: Error in ZHESV. Error=',i
        tmp = 0
        tmp = 1.0/tmp
        STOP
     end If


     AC%Xprev = X

     !*** Now calculate the new vector ***
     X = X + AC%NewMix*F

     Do i=1, Min(AC%N+1,AC%Nmax) 
        X = X - AC%Gamma(i)*(AC%DX(:,i) + AC%NewMix*AC%DF(:,i))
     End Do
  End If
  

  AC%Fprev = F
  AC%Xcurr = X

  AC%N = AC%N + 1
  If (AC%N > AC%Nmax) AC%N = AC%Nmax

  Return
End Subroutine

!******************************************************************************
!
!  Anderson_ResetMix - Resets the mixing history to None
!
!     AC - Anderson context to reset
!
!******************************************************************************

Subroutine Anderson_ResetMix(AC)
  Type  (Anderson_Context), Intent(INOUT) :: AC

  AC%N = -1
  AC%Slot = -1

  Return
End Subroutine


!******************************************************************************
!
!  FreeAnderson - Frees all the data associated with the AC data structure
!    
!      AC -Pointer to the Anderson context to free
!
!******************************************************************************

Subroutine FreeAnderson(AC)
   Type (Anderson_Context), Pointer :: AC  

   DeAllocate(AC%Gamma, AC%Work, AC%Fprev, AC%DF, AC%Matrix, &
        AC%IPIV, AC%DX, AC%Xprev, AC%DupMatrix)

   DeAllocate(AC)

   Return
end Subroutine


!******************************************************************************
!
!  InitAnderson - Initializes and Anderson_Context data structure for use
!
!   AC       - Anderson context created and returned
!   Err_Unit - Output error unit
!   Nmax     - Max number of vectors to keep
!   VecSize  - Size of each vector
!   w0       - Fudge for diagonal to keep linear independece ~1E-4.
!   NewMix   - Mixing factor
!
!******************************************************************************

Subroutine InitAnderson(AC, Err_Unit, Nmax, VecSize, w0, NewMix)
   Type (Anderson_Context), Pointer     :: AC
   Integer,                 Intent(IN)  :: Err_Unit
   Integer,                 Intent(IN)  :: Nmax
   Integer,                 Intent(IN)  :: VecSize
   Real*8,                    Intent(IN)  :: w0
   Real*8,                    Intent(IN)  :: NewMix

   Integer :: i
   Real*8    :: tmp

   Allocate(AC)   !** Allocate the pointer

   AC%Nmax = Nmax          !*** Store the contants
   AC%VecSize = VecSize
   AC%w0 = w0
   AC%NewMix = NewMix
   AC%Err_Unit = Err_Unit

   AC%N = -1                !** Init the rest of the structure
   AC%Slot = -1
   AC%Lwork = 2*Nmax
   Allocate(AC%Gamma(Nmax), AC%Work(AC%Lwork), AC%Fprev(VecSize), &
        AC%DF(VecSize, Nmax), AC%Matrix(Nmax, Nmax), AC%IPIV(Nmax), &
        AC%DX(VecSize, Nmax), AC%Xprev(VecSize), AC%Xcurr(VecSize), &
        AC%DupMatrix(Nmax, Nmax), STAT=i)

   If (i /= 0) then
      Write(Err_Unit,*) 'InitAnderson: Allocate Error! Error=',i
      Write(Err_Unit, *) 'InitAnderson: Nmax=',Nmax, ' * VecSize=',VecSize
      tmp = 0
      tmp = 1.0/tmp
      STOP
   End If


   AC%Matrix = 0

   Return
End Subroutine


End Module

