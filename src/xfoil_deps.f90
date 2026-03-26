! This file is part of Construct2D.
!
! Construct2D is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Construct2D is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Construct2D.  If not, see <http://www.gnu.org/licenses/>.
!
! Contains xfoil paneling subroutine (pangen) and dependencies.
! Converted from fixed-form Fortran 77 to free-form Fortran 90 module.
! Only minor changes made to xfoil source code.
!
! Copyright (C) 2013 -- 2018 Daniel Prosser (this modified version of
! XFoil code)
! Original copyright (C) 2000 Mark Drela (original XFoil code)

module xfoil_deps

  implicit none
  
  private
  public :: pangen, scalc, segspl, lefind, seval

contains

  !-----------------------------------------------------------------------------
  subroutine SPLIND(X, XS, S, N, XS1, XS2)
  !-----------------------------------------------------------------------------
  ! Calculates spline coefficients for X(S).
  ! Specified 1st derivative and/or usual zero 2nd derivative end conditions
  ! are used.  To evaluate the spline at some value of S, use SEVAL and/or
  ! DEVAL.
  !
  ! S        independent variable array (input)
  ! X        dependent variable array   (input)
  ! XS       dX/dS array                (calculated)
  ! N        number of points           (input)
  ! XS1,XS2  optional endpoint derivatives (input)
  !          If present, specified first derivatives are used.
  !          If absent, usual zero second derivative end conditions are used.
  !-----------------------------------------------------------------------------
    integer,          intent(in)    :: N
    double precision, intent(in)    :: X(N), S(N)
    double precision, intent(in), optional :: XS1, XS2
    double precision, intent(inout) :: XS(N)

    integer, parameter :: NMAX = 1000
    double precision   :: A(NMAX), B(NMAX), C(NMAX)
    integer            :: I
    double precision   :: DSM, DSP

    if (N > NMAX) stop 'SPLIND: array overflow, increase NMAX'

    do I = 2, N-1
      DSM  = S(I) - S(I-1)
      DSP  = S(I+1) - S(I)
      B(I) = DSP
      A(I) = 2.0d0*(DSM+DSP)
      C(I) = DSM
      XS(I) = 3.0d0*((X(I+1)-X(I))*DSM/DSP + (X(I)-X(I-1))*DSP/DSM)
    end do

    if (.not. present(XS1)) then
      !----- set zero second derivative end condition
      A(1)  = 2.0d0
      C(1)  = 1.0d0
      XS(1) = 3.0d0*(X(2)-X(1)) / (S(2)-S(1))
    else
      !----- set specified first derivative end condition
      A(1)  = 1.0d0
      C(1)  = 0.0d0
      XS(1) = XS1
    end if

    if (.not. present(XS2)) then
      B(N)  = 1.0d0
      A(N)  = 2.0d0
      XS(N) = 3.0d0*(X(N)-X(N-1)) / (S(N)-S(N-1))
    else
      A(N)  = 1.0d0
      B(N)  = 0.0d0
      XS(N) = XS2
    end if

    !---- solve for derivative array XS
    call TRISOL(A, B, C, XS, N)

  end subroutine SPLIND


  !-----------------------------------------------------------------------------
  double precision function D2VAL(SS, X, XS, S, N)
  !-----------------------------------------------------------------------------
  ! Calculates d2X/dS2(SS).
  ! XS array must have been calculated by SPLIND.
  !-----------------------------------------------------------------------------
    integer,          intent(in) :: N
    double precision, intent(in) :: SS, X(N), XS(N), S(N)

    integer          :: ILOW, I, IMID
    double precision :: DS, T, CX1, CX2

    ILOW = 1
    I    = N
    do while (I - ILOW > 1)
      IMID = (I + ILOW) / 2
      if (SS < S(IMID)) then
        I    = IMID
      else
        ILOW = IMID
      end if
    end do

    DS    = S(I) - S(I-1)
    T     = (SS - S(I-1)) / DS
    CX1   = DS*XS(I-1) - X(I) + X(I-1)
    CX2   = DS*XS(I)   - X(I) + X(I-1)
    D2VAL = (6.0d0*T - 4.0d0)*CX1 + (6.0d0*T - 2.0d0)*CX2
    D2VAL = D2VAL / DS**2

  end function D2VAL


  !-----------------------------------------------------------------------------
  subroutine SCALC(X, Y, S, N)
  !-----------------------------------------------------------------------------
  ! Calculates the arc length array S for a 2-D array of points (X,Y).
  !-----------------------------------------------------------------------------
    integer,          intent(in)  :: N
    double precision, intent(in)  :: X(N), Y(N)
    double precision, intent(out) :: S(N)

    integer :: I

    S(1) = 0.0d0
    do I = 2, N
      S(I) = S(I-1) + sqrt((X(I)-X(I-1))**2 + (Y(I)-Y(I-1))**2)
    end do

  end subroutine SCALC


  !-----------------------------------------------------------------------------
  subroutine SEGSPL(X, XS, S, N)
  !-----------------------------------------------------------------------------
  ! Splines X(S) array just like SPLIND, but allows derivative discontinuities
  ! at segment joints.  Segment joints are defined by identical successive S
  ! values.
  !-----------------------------------------------------------------------------
    integer,          intent(in)    :: N
    double precision, intent(in)    :: X(N), S(N)
    double precision, intent(inout) :: XS(N)

    integer :: ISEG, ISEG0, NSEG

    if (S(1) == S(2)  ) stop 'SEGSPL:  First input point duplicated'
    if (S(N) == S(N-1)) stop 'SEGSPL:  Last  input point duplicated'

    ISEG0 = 1
    do ISEG = 2, N-2
      if (S(ISEG) == S(ISEG+1)) then
        NSEG  = ISEG - ISEG0 + 1
        call SPLIND(X(ISEG0), XS(ISEG0), S(ISEG0), NSEG)
        ISEG0 = ISEG + 1
      end if
    end do

    NSEG = N - ISEG0 + 1
    call SPLIND(X(ISEG0), XS(ISEG0), S(ISEG0), NSEG)

  end subroutine SEGSPL


  !-----------------------------------------------------------------------------
  double precision function CURV(SS, X, XS, Y, YS, S, N)
  !-----------------------------------------------------------------------------
  ! Calculates curvature of splined 2-D curve at S = SS.
  !
  ! S        arc length array of curve
  ! X, Y     coordinate arrays of curve
  ! XS, YS   derivative arrays (calculated earlier by SPLIND)
  !-----------------------------------------------------------------------------
    integer,          intent(in) :: N
    double precision, intent(in) :: SS, X(N), XS(N), Y(N), YS(N), S(N)

    integer          :: ILOW, I, IMID
    double precision :: DS, T
    double precision :: CX1, CX2, CY1, CY2
    double precision :: XD, XDD, YD, YDD, SD

    ILOW = 1
    I    = N
    do while (I - ILOW > 1)
      IMID = (I + ILOW) / 2
      if (SS < S(IMID)) then
        I    = IMID
      else
        ILOW = IMID
      end if
    end do

    DS  = S(I) - S(I-1)
    T   = (SS - S(I-1)) / DS

    CX1 = DS*XS(I-1) - X(I) + X(I-1)
    CX2 = DS*XS(I)   - X(I) + X(I-1)
    XD  = X(I) - X(I-1) + (1.0d0 - 4.0d0*T + 3.0d0*T*T)*CX1 + T*(3.0d0*T - 2.0d0)*CX2
    XDD = (6.0d0*T - 4.0d0)*CX1 + (6.0d0*T - 2.0d0)*CX2

    CY1 = DS*YS(I-1) - Y(I) + Y(I-1)
    CY2 = DS*YS(I)   - Y(I) + Y(I-1)
    YD  = Y(I) - Y(I-1) + (1.0d0 - 4.0d0*T + 3.0d0*T*T)*CY1 + T*(3.0d0*T - 2.0d0)*CY2
    YDD = (6.0d0*T - 4.0d0)*CY1 + (6.0d0*T - 2.0d0)*CY2

    SD   = sqrt(XD*XD + YD*YD)
    SD   = max(SD, 0.001d0*DS)
    CURV = (XD*YDD - YD*XDD) / SD**3

  end function CURV


  !-----------------------------------------------------------------------------
  subroutine LEFIND(SLE, X, XP, Y, YP, S, N)
  !-----------------------------------------------------------------------------
  ! Locates leading edge spline-parameter value SLE.
  !
  ! The defining condition is
  !    (X-XTE, Y-YTE) . (X', Y') = 0   at  S = SLE
  !
  ! i.e. the surface tangent is normal to the chord line connecting
  ! X(SLE), Y(SLE) and the TE point.
  !-----------------------------------------------------------------------------
    integer,          intent(in)  :: N
    double precision, intent(in)  :: X(*), XP(*), Y(*), YP(*), S(*)
    double precision, intent(out) :: SLE

    integer          :: I, ITER
    double precision :: DSEPS, XTE, YTE
    double precision :: DXTE, DYTE, DX, DY, DOTP
    double precision :: XLE, YLE, DXDS, DYDS, DXDD, DYDD
    double precision :: XCHORD, YCHORD, RES, RESS, DSLE

    !---- convergence tolerance
    DSEPS = (S(N) - S(1)) * 1.0d-5

    !---- set trailing edge point coordinates
    XTE = 0.5d0*(X(1) + X(N))
    YTE = 0.5d0*(Y(1) + Y(N))

    !---- get first guess for SLE
    do I = 3, N-2
      DXTE = X(I) - XTE
      DYTE = Y(I) - YTE
      DX   = X(I+1) - X(I)
      DY   = Y(I+1) - Y(I)
      DOTP = DXTE*DX + DYTE*DY
      if (DOTP < 0.0d0) exit
    end do

    SLE = S(I)

    !---- check for sharp LE case
    if (S(I) == S(I-1)) return

    !---- Newton iteration to get exact SLE value
    do ITER = 1, 50
      XLE  = SEVAL(SLE, X, XP, S, N)
      YLE  = SEVAL(SLE, Y, YP, S, N)
      DXDS = DEVAL(SLE, X, XP, S, N)
      DYDS = DEVAL(SLE, Y, YP, S, N)
      DXDD = D2VAL(SLE, X, XP, S, N)
      DYDD = D2VAL(SLE, Y, YP, S, N)

      XCHORD = XLE - XTE
      YCHORD = YLE - YTE

      !------ drive dot product between chord line and LE tangent to zero
      RES  = XCHORD*DXDS + YCHORD*DYDS
      RESS = DXDS*DXDS + DYDS*DYDS &
           + XCHORD*DXDD + YCHORD*DYDD

      !------ Newton delta for SLE
      DSLE = -RES/RESS

      DSLE = max( DSLE, -0.02d0*abs(XCHORD+YCHORD) )
      DSLE = min( DSLE,  0.02d0*abs(XCHORD+YCHORD) )
      SLE  = SLE + DSLE
      if (abs(DSLE) < DSEPS) return
    end do

    write(*,*) 'LEFIND:  LE point not found.  Continuing...'
    SLE = S(I)

  end subroutine LEFIND


  !-----------------------------------------------------------------------------
  double precision function SEVAL(SS, X, XS, S, N)
  !-----------------------------------------------------------------------------
  ! Calculates X(SS).
  ! XS array must have been calculated by SPLIND.
  !-----------------------------------------------------------------------------
    integer,          intent(in) :: N
    double precision, intent(in) :: SS, X(N), XS(N), S(N)

    integer          :: ILOW, I, IMID
    double precision :: DS, T, CX1, CX2

    ILOW = 1
    I    = N
    do while (I - ILOW > 1)
      IMID = (I + ILOW) / 2
      if (SS < S(IMID)) then
        I    = IMID
      else
        ILOW = IMID
      end if
    end do

    DS    = S(I) - S(I-1)
    T     = (SS - S(I-1)) / DS
    CX1   = DS*XS(I-1) - X(I) + X(I-1)
    CX2   = DS*XS(I)   - X(I) + X(I-1)
    SEVAL = T*X(I) + (1.0d0-T)*X(I-1) + (T-T*T)*((1.0d0-T)*CX1 - T*CX2)

  end function SEVAL


  !-----------------------------------------------------------------------------
  double precision function DEVAL(SS, X, XS, S, N)
  !-----------------------------------------------------------------------------
  ! Calculates dX/dS(SS).
  ! XS array must have been calculated by SPLIND.
  !-----------------------------------------------------------------------------
    integer,          intent(in) :: N
    double precision, intent(in) :: SS, X(N), XS(N), S(N)

    integer          :: ILOW, I, IMID
    double precision :: DS, T, CX1, CX2

    ILOW = 1
    I    = N
    do while (I - ILOW > 1)
      IMID = (I + ILOW) / 2
      if (SS < S(IMID)) then
        I    = IMID
      else
        ILOW = IMID
      end if
    end do

    DS    = S(I) - S(I-1)
    T     = (SS - S(I-1)) / DS
    CX1   = DS*XS(I-1) - X(I) + X(I-1)
    CX2   = DS*XS(I)   - X(I) + X(I-1)
    DEVAL = X(I) - X(I-1) + (1.0d0 - 4.0d0*T + 3.0d0*T*T)*CX1 + T*(3.0d0*T - 2.0d0)*CX2
    DEVAL = DEVAL / DS

  end function DEVAL


  !-----------------------------------------------------------------------------
  subroutine TRISOL(A, B, C, D, KK)
  !-----------------------------------------------------------------------------
  ! Solves KK long tri-diagonal system:
  !
  !     A C          D
  !     B A C        D
  !       B A .      .
  !         . . C    .
  !           B A    D
  !
  ! The righthand side D is replaced by the solution.  A and C are destroyed.
  !-----------------------------------------------------------------------------
    integer,          intent(in)    :: KK
    double precision, intent(in)    :: B(KK)
    double precision, intent(inout) :: A(KK), C(KK), D(KK)

    integer :: K, KM

    do K = 2, KK
      KM    = K - 1
      C(KM) = C(KM) / A(KM)
      D(KM) = D(KM) / A(KM)
      A(K)  = A(K) - B(K)*C(KM)
      D(K)  = D(K) - B(K)*D(KM)
    end do

    D(KK) = D(KK) / A(KK)

    do K = KK-1, 1, -1
      D(K) = D(K) - C(K)*D(K+1)
    end do

  end subroutine TRISOL


  !-----------------------------------------------------------------------------
  subroutine PANGEN(XNEW, YNEW, NPAN, XBIN, YBIN, NB, CVPAR, CTERAT, &
                    CTRRAT, XSREF1, XSREF2, XPREF1, XPREF2)
  !-----------------------------------------------------------------------------
  ! Set paneling distribution from buffer airfoil geometry, thus creating the
  ! current airfoil.
  !
  ! If REFINE=True, bunch points at x=XSREF on top side and at x=XPREF on
  ! bottom side by setting a fictitious local curvature of CTRRAT*(LE curv).
  !-----------------------------------------------------------------------------
    integer,          intent(in)                :: NPAN, NB
    double precision, dimension(NB),   intent(in)  :: XBIN, YBIN
    double precision, dimension(NPAN), intent(out) :: XNEW, YNEW
    double precision, intent(in) :: CVPAR, CTERAT, CTRRAT
    double precision, intent(in) :: XSREF1, XSREF2, XPREF1, XPREF2

    double precision, allocatable :: W1(:), W2(:), W3(:), W4(:), W5(:), W6(:)
    double precision, allocatable :: SB(:), SNEW(:), S(:)
    double precision, allocatable :: XB(:), YB(:), XBP(:), YBP(:), X(:), Y(:)

    integer :: IPFAC, N, NN, NCORN, NK, K, IBLE, NFRAC1, NN1, NN2
    integer :: work_size, max_surface_pts, stat_alloc
    integer :: I, IB, J, IND
    logical :: converged
    double precision :: SBREF, SBLE, CVLE, XBLE, YBLE, XBTE, YBTE, CHBSQ
    double precision :: FRAC, SBK, CVK, CVSUM, CVAVG, CVMAX, CC, CVTE
    double precision :: SMOOL, SMOOSQ, DSM, DSP, DSO
    double precision :: DSAVG, DSAVG1, DSAVG2, RDSTE, RTF
    double precision :: XBCORN, YBCORN, SBCORN, XOC
    double precision :: CV1, CV2, CV3, CVS1, CVS2, CVS3
    double precision :: CAVM, CAVP, CAVM_S1, CAVM_S2, CAVP_S2, CAVP_S3
    double precision :: FM, FP, REZ, RLX, DMAX, DS, DDS, DSRAT

    if (NB < 2) then
      write(*,*) 'PANGEN: Buffer airfoil not available.'
      return
    end if

    !---- Number of temporary nodes for panel distribution calculation
    !       exceeds the specified panel number by factor of IPFAC.
    IPFAC = 3
    IPFAC = 5

    !---- Allocate working arrays from current problem size.
    NN = IPFAC*(NPAN-1) + 1
    work_size = max(NB, NN)
    max_surface_pts = NPAN + NB

    allocate(W1(work_size), W2(work_size), W3(work_size), W4(work_size), &
             W5(work_size), W6(work_size),                                &
             SB(NB), SNEW(NN),                                            &
             XB(NB), YB(NB), XBP(NB), YBP(NB),                            &
             X(max_surface_pts), Y(max_surface_pts), S(max_surface_pts),  &
             stat=stat_alloc)
    if (stat_alloc /= 0) then
      write(*,*) 'PANGEN: Failed to allocate working arrays.'
      return
    end if

    !---- Populate XB and YB arrays
    XB(1:NB) = XBIN
    YB(1:NB) = YBIN

    !---- Number of temporary nodes for panel distribution calculation
    !       exceeds the specified panel number by factor of IPFAC.
    !---- number of airfoil panel points
    N = NPAN

    !---- set arc length spline parameter
    call SCALC(XB, YB, SB, NB)

    !---- spline raw airfoil coordinates
    call SEGSPL(XB, XBP, SB, NB)
    call SEGSPL(YB, YBP, SB, NB)

    !---- normalizing length (~ chord)
    SBREF = 0.5d0*(SB(NB) - SB(1))

    !---- set up curvature array
    do I = 1, NB
      W5(I) = abs( CURV(SB(I), XB, XBP, YB, YBP, SB, NB) ) * SBREF
    end do

    !---- locate LE point arc length value and the normalized curvature there
    call LEFIND(SBLE, XB, XBP, YB, YBP, SB, NB)
    CVLE = abs( CURV(SBLE, XB, XBP, YB, YBP, SB, NB) ) * SBREF

    !---- check for doubled point (sharp corner) at LE
    IBLE = 0
    do I = 1, NB-1
      if (SBLE == SB(I) .and. SBLE == SB(I+1)) then
        IBLE = I
        write(*,*)
        write(*,*) 'Sharp leading edge'
        exit
      end if
    end do

    !---- set LE, TE points
    XBLE  = SEVAL(SBLE, XB, XBP, SB, NB)
    YBLE  = SEVAL(SBLE, YB, YBP, SB, NB)
    XBTE  = 0.5d0*(XB(1) + XB(NB))
    YBTE  = 0.5d0*(YB(1) + YB(NB))
    CHBSQ = (XBTE-XBLE)**2 + (YBTE-YBLE)**2

    !---- set average curvature over 2*NK+1 points within Rcurv of LE point
    NK    = 3
    CVSUM = 0.0d0
    do K = -NK, NK
      FRAC  = dble(K) / dble(NK)
      SBK   = SBLE + FRAC*SBREF / max(CVLE, 20.0d0)
      CVK   = abs( CURV(SBK, XB, XBP, YB, YBP, SB, NB) ) * SBREF
      CVSUM = CVSUM + CVK
    end do
    CVAVG = CVSUM / dble(2*NK+1)

    !---- dummy curvature for sharp LE
    if (IBLE /= 0) CVAVG = 10.0d0

    !---- set curvature attraction coefficient actually used
    CC = 6.0d0 * CVPAR

    !---- set artificial curvature at TE to bunch panels there
    CVTE   = CVAVG * CTERAT
    W5(1)  = CVTE
    W5(NB) = CVTE

    !**** smooth curvature array for smoother panel size distribution ****
    !---- set smoothing length = 1 / averaged LE curvature, but
    !-    no more than 5% of chord and no less than 1/4 average panel spacing
    SMOOL  = max( 1.0d0/max(CVAVG, 20.0d0), 0.25d0/dble(NPAN/2) )
    SMOOSQ = (SMOOL*SBREF)**2

    !---- set up tri-diagonal system for smoothed curvatures
    W2(1) = 1.0d0
    W3(1) = 0.0d0
    do I = 2, NB-1
      DSM = SB(I) - SB(I-1)
      DSP = SB(I+1) - SB(I)
      DSO = 0.5d0*(SB(I+1) - SB(I-1))
      if (DSM == 0.0d0 .or. DSP == 0.0d0) then
        !------- leave curvature at corner point unchanged
        W1(I) = 0.0d0
        W2(I) = 1.0d0
        W3(I) = 0.0d0
      else
        W1(I) =  SMOOSQ * (          - 1.0d0/DSM) / DSO
        W2(I) =  SMOOSQ * ( 1.0d0/DSP + 1.0d0/DSM) / DSO  +  1.0d0
        W3(I) =  SMOOSQ * (-1.0d0/DSP             ) / DSO
      end if
    end do
    W1(NB) = 0.0d0
    W2(NB) = 1.0d0

    !---- fix curvature at LE point by modifying equations adjacent to LE
    do I = 2, NB-1
      if (SB(I) == SBLE .or. I == IBLE .or. I == IBLE+1) then
        !------- if node falls right on LE point, fix curvature there
        W1(I) = 0.0d0
        W2(I) = 1.0d0
        W3(I) = 0.0d0
        W5(I) = CVLE
      else if (SB(I-1) < SBLE .and. SB(I) > SBLE) then
        !------- modify equation at node just before LE point
        DSM = SB(I-1) - SB(I-2)
        DSP = SBLE    - SB(I-1)
        DSO = 0.5d0*(SBLE - SB(I-2))
        W1(I-1) =  SMOOSQ * (          - 1.0d0/DSM) / DSO
        W2(I-1) =  SMOOSQ * ( 1.0d0/DSP + 1.0d0/DSM) / DSO  +  1.0d0
        W3(I-1) =  0.0d0
        W5(I-1) =  W5(I-1) + SMOOSQ*CVLE/(DSP*DSO)
        !------- modify equation at node just after LE point
        DSM = SB(I)   - SBLE
        DSP = SB(I+1) - SB(I)
        DSO = 0.5d0*(SB(I+1) - SBLE)
        W1(I) =  0.0d0
        W2(I) =  SMOOSQ * ( 1.0d0/DSP + 1.0d0/DSM) / DSO  +  1.0d0
        W3(I) =  SMOOSQ * (-1.0d0/DSP             ) / DSO
        W5(I) =  W5(I) + SMOOSQ*CVLE/(DSM*DSO)
        exit  ! done modifying adjacent nodes
      end if
    end do

    !---- set artificial curvature at bunching points and fix it there
    do I = 2, NB-1
      !------ chord-based x/c coordinate
      XOC = ((XB(I)-XBLE)*(XBTE-XBLE) &
           + (YB(I)-YBLE)*(YBTE-YBLE)) / CHBSQ
      if (SB(I) < SBLE) then
        !------- check if top side point is in refinement area
        if (XOC > XSREF1 .and. XOC < XSREF2) then
          W1(I) = 0.0d0
          W2(I) = 1.0d0
          W3(I) = 0.0d0
          W5(I) = CVLE*CTRRAT
        end if
      else
        !------- check if bottom side point is in refinement area
        if (XOC > XPREF1 .and. XOC < XPREF2) then
          W1(I) = 0.0d0
          W2(I) = 1.0d0
          W3(I) = 0.0d0
          W5(I) = CVLE*CTRRAT
        end if
      end if
    end do

    !---- solve for smoothed curvature array W5
    if (IBLE == 0) then
      call TRISOL(W2, W1, W3, W5, NB)
    else
      I = 1
      call TRISOL(W2(I), W1(I), W3(I), W5(I), IBLE)
      I = IBLE + 1
      call TRISOL(W2(I), W1(I), W3(I), W5(I), NB-IBLE)
    end if

    !---- find max curvature
    CVMAX = 0.0d0
    do I = 1, NB
      CVMAX = max( CVMAX, abs(W5(I)) )
    end do

    !---- normalize curvature array
    do I = 1, NB
      W5(I) = W5(I) / CVMAX
    end do

    !---- spline curvature array
    call SEGSPL(W5, W6, SB, NB)

    !---- Set initial guess for node positions uniform in s.
    !     More nodes than specified (by factor of IPFAC) are temporarily used
    !     for more reliable convergence.
    NN = IPFAC*(N-1) + 1

    !---- ratio of lengths of panel at TE to one away from the TE
    RDSTE = 0.667d0
    RTF   = (RDSTE - 1.0d0)*2.0d0 + 1.0d0

    if (IBLE == 0) then
      DSAVG   = (SB(NB)-SB(1)) / (dble(NN-3) + 2.0d0*RTF)
      SNEW(1) = SB(1)
      do I = 2, NN-1
        SNEW(I) = SB(1) + DSAVG * (dble(I-2) + RTF)
      end do
      SNEW(NN) = SB(NB)
      NN1 = 0   ! DP mod: initialised to suppress compiler warning
    else
      NFRAC1  = (N * IBLE) / NB
      NN1     = IPFAC*(NFRAC1-1) + 1
      DSAVG1  = (SBLE-SB(1)) / (dble(NN1-2) + RTF)
      SNEW(1) = SB(1)
      do I = 2, NN1
        SNEW(I) = SB(1) + DSAVG1 * (dble(I-2) + RTF)
      end do
      NN2    = NN - NN1 + 1
      DSAVG2 = (SB(NB)-SBLE) / (dble(NN2-2) + RTF)
      do I = 2, NN2-1
        SNEW(I-1+NN1) = SBLE + DSAVG2 * (dble(I-2) + RTF)
      end do
      SNEW(NN) = SB(NB)
    end if

    !---- Newton iteration loop for new node positions
    converged = .false.
    do I = 1, 20
      !------ set up tri-diagonal system for node position deltas
      CV1  = SEVAL(SNEW(1), W5, W6, SB, NB)
      CV2  = SEVAL(SNEW(2), W5, W6, SB, NB)
      CVS1 = DEVAL(SNEW(1), W5, W6, SB, NB)
      CVS2 = DEVAL(SNEW(2), W5, W6, SB, NB)

      CAVM = sqrt(CV1**2 + CV2**2)
      if (CAVM == 0.0d0) then
        CAVM_S1 = 0.0d0
        CAVM_S2 = 0.0d0
      else
        CAVM_S1 = CVS1 * CV1/CAVM
        CAVM_S2 = CVS2 * CV2/CAVM
      end if

      do J = 2, NN-1
        DSM  = SNEW(J) - SNEW(J-1)
        DSP  = SNEW(J) - SNEW(J+1)
        CV3  = SEVAL(SNEW(J+1), W5, W6, SB, NB)
        CVS3 = DEVAL(SNEW(J+1), W5, W6, SB, NB)

        CAVP = sqrt(CV3**2 + CV2**2)
        if (CAVP == 0.0d0) then
          CAVP_S2 = 0.0d0
          CAVP_S3 = 0.0d0
        else
          CAVP_S2 = CVS2 * CV2/CAVP
          CAVP_S3 = CVS3 * CV3/CAVP
        end if

        FM = CC*CAVM + 1.0d0
        FP = CC*CAVP + 1.0d0

        REZ = DSP*FP + DSM*FM

        !-------- lower, main, and upper diagonals
        W1(J) =      -FM  +  CC*               DSM*CAVM_S1
        W2(J) =  FP + FM  +  CC*(DSP*CAVP_S2 + DSM*CAVM_S2)
        W3(J) = -FP       +  CC* DSP*CAVP_S3

        !-------- residual, requiring (1 + C*curv)*deltaS to be equal on both sides
        W4(J) = -REZ

        CV1     = CV2
        CV2     = CV3
        CVS1    = CVS2
        CVS2    = CVS3
        CAVM    = CAVP
        CAVM_S1 = CAVP_S2
        CAVM_S2 = CAVP_S3
      end do

      !------ fix endpoints (at TE)
      W2(1)  = 1.0d0
      W3(1)  = 0.0d0
      W4(1)  = 0.0d0
      W1(NN) = 0.0d0
      W2(NN) = 1.0d0
      W4(NN) = 0.0d0

      if (RTF /= 1.0d0) then
        !------- fudge equations adjacent to TE to get TE panel length ratio RTF
        J = 2
        W4(J) = -((SNEW(J) - SNEW(J-1)) + RTF*(SNEW(J) - SNEW(J+1)))
        W1(J) = -1.0d0
        W2(J) =  1.0d0 + RTF
        W3(J) =        - RTF

        J = NN-1
        W4(J) = -((SNEW(J) - SNEW(J+1)) + RTF*(SNEW(J) - SNEW(J-1)))
        W3(J) = -1.0d0
        W2(J) =  1.0d0 + RTF
        W1(J) =        - RTF
      end if

      !------ fix sharp LE point
      if (IBLE /= 0) then
        J     = NN1
        W1(J) = 0.0d0
        W2(J) = 1.0d0
        W3(J) = 0.0d0
        W4(J) = SBLE - SNEW(J)
      end if

      !------ solve for changes W4 in node position arc length values
      call TRISOL(W2, W1, W3, W4, NN)

      !------ find under-relaxation factor to keep nodes from changing order
      RLX  = 1.0d0
      DMAX = 0.0d0
      do J = 1, NN-1
        DS    = SNEW(J+1) - SNEW(J)
        DDS   = W4(J+1) - W4(J)
        DSRAT = 1.0d0 + RLX*DDS/DS
        if (DSRAT > 4.0d0) RLX = (4.0d0-1.0d0)*DS/DDS
        if (DSRAT < 0.2d0) RLX = (0.2d0-1.0d0)*DS/DDS
        DMAX = max(abs(W4(J)), DMAX)
      end do

      !------ update node position
      do J = 2, NN-1
        SNEW(J) = SNEW(J) + RLX*W4(J)
      end do

      if (abs(DMAX) < 1.0d-3) then
        converged = .true.
        exit
      end if
    end do

    if (.not. converged) write(*,*) 'Paneling convergence failed.  Continuing anyway...'

    !---- set new panel node coordinates
    do I = 1, N
      IND  = IPFAC*(I-1) + 1
      S(I) = SNEW(IND)
      X(I) = SEVAL(SNEW(IND), XB, XBP, SB, NB)
      Y(I) = SEVAL(SNEW(IND), YB, YBP, SB, NB)
    end do

    !---- go over buffer airfoil again, checking for corners (double points)
    NCORN = 0
    outer: do IB = 1, NB-1
      if (SB(IB) == SB(IB+1)) then
        !------- found one !
        NCORN  = NCORN + 1
        XBCORN = XB(IB)
        YBCORN = YB(IB)
        SBCORN = SB(IB)

        !------- find current-airfoil panel which contains corner
        do I = 1, N
          !--------- keep stepping until first node past corner
          if (S(I) <= SBCORN) cycle

          !---------- move remainder of panel nodes to make room for additional node
          do J = N, I, -1
            X(J+1) = X(J)
            Y(J+1) = Y(J)
            S(J+1) = S(J)
          end do
          N = N + 1

          if (N > max_surface_pts) &
            stop 'PANEL: Too many panels. Increase working array size.'

          X(I) = XBCORN
          Y(I) = YBCORN
          S(I) = SBCORN

          !---------- shift nodes adjacent to corner to keep panel sizes comparable
          if (I-2 >= 1) then
            S(I-1) = 0.5d0*(S(I) + S(I-2))
            X(I-1) = SEVAL(S(I-1), XB, XBP, SB, NB)
            Y(I-1) = SEVAL(S(I-1), YB, YBP, SB, NB)
          end if

          if (I+2 <= N) then
            S(I+1) = 0.5d0*(S(I) + S(I+2))
            X(I+1) = SEVAL(S(I+1), XB, XBP, SB, NB)
            Y(I+1) = SEVAL(S(I+1), YB, YBP, SB, NB)
          end if

          !---------- go on to next input geometry point to check for corner
          cycle outer
        end do
      end if
    end do outer

    !---- Output new geometry
    XNEW(1:NPAN) = X(1:NPAN)
    YNEW(1:NPAN) = Y(1:NPAN)

    deallocate(W1, W2, W3, W4, W5, W6, SB, SNEW, XB, YB, XBP, YBP, X, Y, S)

  end subroutine PANGEN

end module xfoil_deps
