!+
! PROGRAM AbbottVonDoenhoff                            ! pdas/naca456/avd.f90
! ---------------------------------------------------------------------------
! PURPOSE - Make tables of the data in Abbott and von Doenhoff book
!   entitled "Theory of Airfoil Sections".
!   The output will be HTML 5 files for viewing in a browser.
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 20Nov01  0.1   RLC   Original coding (used program in \atmos as guide )
! 21Nov01  0.2   RLC   Added the actual airfoil calculations
! 23Nov01  0.3   RLC   Cambered sections both interp and non-interp
! 27Nov01  0.4   RLC   Rearranged tables at head of each section
! 30Nov01  0.5   RLC   Changed output to four figures; added leRadius
! 02Dec01  0.6   RLC   Renamed some procedures. Added remainder of profiles
! 07Dec01  0.7   RLC   Finally got leRadius, leSlope right
! 20Jan02  0.8   RLC   Put #top at top of page
! 26Jan02  0.9   RLC   Separated sections into 3 groups
! 07Apr02  1.0   RLC   Added the missing profiles and sections. Noted charset
! 08Feb03  1.1   RLC   Made the resulting files XHTML-1.1 STRICT 
! 16Feb03  1.2   RLC   Fixed items that caused non-compliance with XHTML 1.1
! 15Jan09  1.3   RLC   Additional fixes for XHTML 1.1 and file names are *.html
! 16Jan09  1.4   RLC   Added WriteBanner,WriteCrumb,WriteFooter
! 17Jan09  1.5   RLC   Added MainPage; removed Welcome and debug file 
! 07Nov10  1.6   RLC   Revised to create HTML 5
! 18Jul17  1.7   RLC   Correcting some HTML 5 errors
! 24Apr21  1.8   RLC   Minor revisions, moved to httpsdocs from httpdocs
! 31May22  1.9   RLC   Removed td.links from WriteHead (obsolete)
!
! NOTES- If you should recompile this source code, you should use the
!   appropriate compiler flag making the default real kind as 64-bit.
!   With gfortran this is -fdefault-real-8
!-------------------------------------------------------------------------------

INCLUDE 'splprocs.f90'
INCLUDE 'epspsi.f90'
INCLUDE 'nacax.f90'
! NOTE - If you want to modify and recompile this program, you will need the
! three files above as well as avd.f90. These are also used by naca456.f90.

!+
MODULE AbbottVonDoenhoffSubs
! ------------------------------------------------------------------------------
! PURPOSE - Subroutines called by Program AbbottVonDoenhoff.
!   Make tables of the data in Abbott and von Doenhoff.
!   This module is not used by the naca456 program.  

IMPLICIT NONE

  PUBLIC:: WriteFourDigitProfile
  PUBLIC:: WriteModifiedFourDigitProfile
  PUBLIC:: WriteSixSeriesProfile

  PUBLIC:: WriteTwoDigitMeanLine
  PUBLIC:: WriteThreeDigitMeanLine
  PUBLIC:: WriteThreeDigitReflexMeanLine
  PUBLIC:: WriteSixSeriesMeanLine
  PUBLIC:: WriteSixSeriesModifiedMeanLine

  PUBLIC:: WriteFourDigitSection
  PUBLIC:: WriteModifiedFourDigitSection
  PUBLIC:: WriteFiveDigitSection
  PUBLIC:: WriteFiveDigitReflexSection
  PUBLIC:: WriteSixSeriesSection

  PRIVATE:: WriteOneProfile
  PRIVATE:: WriteOneMeanLine
  PRIVATE:: WriteOneSection

  PRIVATE:: WriteColumns, WriteColumns2
  PUBLIC:: WriteTextTable, WriteTextTableRows

  PUBLIC:: WriteBanner,WriteCrumb,WriteFooter,WriteHead
!-------------------------------------------------------------------------------

CONTAINS


!+
SUBROUTINE WriteFourDigitProfile(efu,toc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a four-digit profile
USE NACAauxilary,ONLY: Thickness4,LeadingEdgeRadius4
! Uses module procedure WriteOneProfile

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: toc      ! thickness / chord   (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   ! profile name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: rle ! radius of leading edge, fraction of chord
  REAL,DIMENSION(SIZE(x)):: yt,ytp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
  CALL Thickness4(toc, x,yt,ytp)   ! compute yt and ytp at each x-value
  rle=LeadingEdgeRadius4(toc)
  CALL WriteOneProfile(efu,name,id, x,yt,ytp, rle)
  RETURN
END Subroutine WriteFourDigitProfile   ! ---------------------------------------

!+
SUBROUTINE WriteModifiedFourDigitProfile(efu,leIndex,xmaxt,toc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a modified four-digit profile
USE NACAauxilary,ONLY: Thickness4M,LeadingEdgeRadius4M

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: leIndex  ! leading edge index
                             ! coded to allow leIndex to be a REAL
  REAL,INTENT(IN):: xmaxt    ! x-coor of maximum thickness
  REAL,INTENT(IN):: toc      ! thickness / chord   (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   ! profile name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: rle
  REAL,DIMENSION(SIZE(x)):: yt,ytp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
! compute yt and ytp at each x-value:
  CALL Thickness4M(toc, leIndex, xmaxt, x,yt,ytp)
  rle=LeadingEdgeRadius4M(toc,leIndex)
  CALL WriteOneProfile(efu,name,id, x,yt,ytp, rle)
  RETURN
END Subroutine WriteModifiedFourDigitProfile   ! -------------------------------

!+
SUBROUTINE WriteSixSeriesProfile(efu,family,toc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a six-series profile.
USE NACAauxilary,ONLY: Thickness6,LeadingEdgeRadius6

  INTEGER,INTENT(IN):: efu   ! external file unit
  INTEGER,INTENT(IN):: family  ! =1 63; =2 64; =3 65; =4 66; =5 67
                               ! =6 63A; =7 64A; =8 65A
  REAL,INTENT(IN):: toc      ! thickness / chord (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   ! profile name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: rle   ! leading edge radius, fraction of chord
  REAL,DIMENSION(SIZE(x)):: yt,ytp   ! tables of y and dy/dx
!-------------------------------------------------------------------------------
! compute yt and ytp at each x-value:
  CALL Thickness6(family,toc, x,yt,ytp)
  rle=LeadingEdgeRadius6(family,toc)
  CALL WriteOneProfile(efu, name,id, x,yt,ytp, rle)
  RETURN
END Subroutine WriteSixSeriesProfile   ! ---------------------------------------

!+
SUBROUTINE WriteTwoDigitMeanLine(efu,cmax,xmaxc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a four-digit meanline.
USE NACAauxilary,ONLY: MeanLine2

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cmax   ! max. camber (fraction of chord)
  REAL,INTENT(IN):: xmaxc  ! x-location of max. camber (fraction of chord) 
  CHARACTER(LEN=*),INTENT(IN):: name   ! mean line name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
! compute ymean and ymeanp at each x-value:
  CALL MeanLine2(cmax,xmaxc, x,ymean,ymeanp)
  CALL WriteOneMeanLine(efu, name,id, x,ymean,ymeanp)
  RETURN
END Subroutine WriteTwoDigitMeanLine   ! --------------------------------------

!+
SUBROUTINE WriteThreeDigitMeanLine(efu,cl,xmaxc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a three-digit meanline.
! EXAMPLE - A 210 mean line has design CL=0.3 and max camber at 5 percent
!   (first digit is 2/3 CL in tenths; 
!    second digit is twice x of max camber in percent chord;
!    third digit, zero, indicates non-reflexed )
! REF - Eq. 6.6, p.115 of Abbott and von Doenhoff

USE NACAauxilary,ONLY: MeanLine3

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cl    ! design lift coefficient
  REAL,INTENT(IN):: xmaxc !  x-location of max. camber (fraction of chord) 
  CHARACTER(LEN=*),INTENT(IN):: name ! mean line name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
! compute ymean and ymeanp at each x-value:
  CALL MeanLine3(cl,xmaxc, x,ymean,ymeanp)
  CALL WriteOneMeanLine(efu, name,id, x,ymean,ymeanp)
  RETURN
END Subroutine WriteThreeDigitMeanLine   ! -------------------------------------

!+
SUBROUTINE WriteThreeDigitReflexMeanLine(efu,cl,xmaxc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a 3-digit reflex meanline
! REF - p.8 of NASA Technical Memorandum 4741 and NACA Report 537.
USE NACAauxilary,ONLY: MeanLine3Reflex

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cl    ! design lift coefficient (same as non-reflex)
  REAL,INTENT(IN):: xmaxc !  x-location of max. camber (fraction of chord) 
  CHARACTER(LEN=*),INTENT(IN):: name   ! mean line name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
! compute ymean and ymeanp at each x-value:
  CALL MeanLine3Reflex(cl,xmaxc, x,ymean,ymeanp)
  CALL WriteOneMeanLine(efu, name,id, x,ymean,ymeanp)
  RETURN
END Subroutine WriteThreeDigitReflexMeanLine   ! -------------------------------

!+
SUBROUTINE WriteSixSeriesMeanLine(efu,a,cl,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a six-series meanline
USE NACAauxilary,ONLY: MeanLine6

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: a    ! x-coor of aft-position of pressure rooftop
  REAL,INTENT(IN):: cl   ! design lift coefficient
  CHARACTER(LEN=*),INTENT(IN):: name   ! mean line name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
! compute ymean and ymeanp at each x-value:
  CALL MeanLine6(a,cl, x,ymean,ymeanp)
  CALL WriteOneMeanLine(efu, name,id, x,ymean,ymeanp)
  RETURN
END Subroutine WriteSixSeriesMeanLine   ! --------------------------------------

!+
SUBROUTINE WriteSixSeriesModifiedMeanLine(efu,cl,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a six-series modified meanline
! NOTE - a is not input; a is always 0.8 
USE NACAauxilary,ONLY: MeanLine6M

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cl       ! design lift coefficient
  CHARACTER(LEN=*),INTENT(IN):: name   ! mean line name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp
!-------------------------------------------------------------------------------
! compute ymean and ymeanp at each x-value:
  CALL MeanLine6M(cl, x,ymean,ymeanp)
  CALL WriteOneMeanLine(efu, name,id, x,ymean,ymeanp)
  RETURN
END Subroutine WriteSixSeriesModifiedMeanLine   ! ------------------------------

!+
SUBROUTINE WriteFourDigitSection(efu,cmax,xmaxc,toc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a four-digit airfoil section
!   Three quantities (cmax,xmaxc and toc define the airfoil)
USE NACAauxilary,ONLY: Thickness4,LeadingEdgeRadius4,MeanLine2
USE NACAauxilary,ONLY: CombineThicknessAndCamber,InterpolateCombinedAirfoil

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cmax ! max. camber (fraction of chord)
  REAL,INTENT(IN):: xmaxc  ! x-location of max. camber (fraction of chord) 
  REAL,INTENT(IN):: toc  ! thickness / chord (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   ! section name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: leSlope
  REAL:: rle   ! leading edge radius, fraction of chord
  REAL,DIMENSION(SIZE(x)):: yt   ! table of half-thickness
  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
  REAL,DIMENSION(SIZE(x)):: xupper,yupper, xlower,ylower
  REAL,DIMENSION(SIZE(x)):: yup,ylo

!-------------------------------------------------------------------------------
  CALL Thickness4(toc, x,yt)
  rle=LeadingEdgeRadius4(toc)
  leSlope=2.0*cmax/xmaxc
  CALL MeanLine2(cmax,xmaxc, x,ymean,ymeanp)
  CALL CombineThicknessAndCamber(x,yt,ymean,ymeanp, &
    xupper,yupper, xlower,ylower)
  CALL InterpolateCombinedAirfoil(x,yt,ymean,ymeanp, yup,ylo)
  CALL WriteOneSection(efu,name,id, x,yup, x,ylo, &
    xupper,yupper, xlower,ylower, rle,leSlope)
  RETURN
END Subroutine WriteFourDigitSection   ! ---------------------------------------

!+
SUBROUTINE WriteModifiedFourDigitSection(efu,cmax,xmaxc, &
  leIndex,xmaxt,toc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a modified four-digit airfoil.
! NOTE - leIndex is coded as REAL, although it is usually integer.

USE NACAauxilary,ONLY: Thickness4M,LeadingEdgeRadius4M,MeanLine2
USE NACAauxilary,ONLY: CombineThicknessAndCamber,InterpolateCombinedAirfoil

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cmax   ! max. camber (fraction of chord)
  REAL,INTENT(IN):: xmaxc  ! x-location of max. camber (fraction of chord) 
  REAL,INTENT(IN):: leIndex   ! leading edge index
  REAL,INTENT(IN):: xmaxt  ! x-location of max. thickness (fraction of chord) 
  REAL,INTENT(IN):: toc   ! thickness / chord (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   ! section name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: leSlope
  REAL:: rle   ! leading edge radius, fraction of chord
  REAL,DIMENSION(SIZE(x)):: yt
  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
  REAL,DIMENSION(SIZE(x)):: xupper,yupper, xlower,ylower
  REAL,DIMENSION(SIZE(x)):: yup,ylo
!-------------------------------------------------------------------------------
  CALL Thickness4M(leIndex,xmaxt,toc, x,yt)
  rle=LeadingEdgeRadius4M(leIndex,toc)
  leSlope=2.0*cmax/xmaxc
  CALL MeanLine2(cmax,xmaxc, x,ymean,ymeanp)
  CALL CombineThicknessAndCamber(x,yt,ymean,ymeanp, &
    xupper,yupper, xlower,ylower)
  CALL InterpolateCombinedAirfoil(x,yt,ymean,ymeanp, yup,ylo)
  CALL WriteOneSection(efu,name,id, x,yup, x,ylo, &
    xupper,yupper, xlower,ylower, rle,leSlope)
  RETURN
END Subroutine WriteModifiedFourDigitSection   ! -------------------------------

!+
SUBROUTINE WriteFiveDigitSection(efu,cl,xmaxc,toc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a five-digit airfoil.
USE NACAauxilary,ONLY: Thickness4,LeadingEdgeRadius4,MeanLine3,GetRk1
USE NACAauxilary,ONLY: CombineThicknessAndCamber,InterpolateCombinedAirfoil

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cl  ! design lift coefficient
  REAL,INTENT(IN):: xmaxc ! x-coor of max camber (fraction of chord)
  REAL,INTENT(IN):: toc  ! thickness / chord (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   ! section name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: k1    ! from GetRk1
  REAL:: leSlope
  REAL:: r    ! from GetRk1
  REAL:: rle   ! leading edge radius, fraction of chord (from LeadingEdgeRadius)
  REAL,DIMENSION(SIZE(x)):: yt             ! table of HALF-thickness
  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx for mean line
  REAL,DIMENSION(SIZE(x)):: xupper,yupper  ! tables of (x,y) for upper surface
  REAL,DIMENSION(SIZE(x)):: xlower,ylower  ! tables of (x,y) for lower surface
  REAL,DIMENSION(SIZE(x)):: yup,ylo        ! interpolated values, upper and lower 
!-------------------------------------------------------------------------------
  CALL Thickness4(toc, x,yt)
  rle=LeadingEdgeRadius4(toc)
  CALL GetRk1(xmaxc,r,k1)
  leSlope=r*r*(3.0-r)*k1*cl/1.8
  CALL MeanLine3(cl,xmaxc, x,ymean,ymeanp)
  CALL CombineThicknessAndCamber(x,yt,ymean,ymeanp, &
    xupper,yupper, xlower,ylower)
  CALL InterpolateCombinedAirfoil(x,yt,ymean,ymeanp, yup,ylo)
  CALL WriteOneSection(efu, name,id, x,yup, x,ylo, &
    xupper,yupper, xlower,ylower, rle,leSlope)
  RETURN
END Subroutine WriteFiveDigitSection   ! ---------------------------------------

!+
SUBROUTINE WriteFiveDigitReflexSection(efu,cl,xmaxc,toc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a 5-digit reflex airfoil.
USE NACAauxilary,ONLY: Thickness4,LeadingEdgeRadius4,MeanLine3Reflex,GetRk1k2
USE NACAauxilary,ONLY: CombineThicknessAndCamber,InterpolateCombinedAirfoil

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cl   ! design lift coefficient
  REAL,INTENT(IN):: xmaxc   ! x-coor of max camber (fraction of chord)
  REAL,INTENT(IN):: toc   ! thickness / chord (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   !  section name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: k1,k21
  REAL:: leSlope
  REAL:: mr3,r,r3
  REAL:: rle   ! leading edge radius, fraction of chord
  REAL,DIMENSION(SIZE(x)):: yt
  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
  REAL,DIMENSION(SIZE(x)):: xupper,yupper, xlower,ylower
  REAL,DIMENSION(SIZE(x)):: yup,ylo
!-------------------------------------------------------------------------------
  CALL Thickness4(toc, x,yt)
  rle=LeadingEdgeRadius4(toc)
  CALL GetRk1k2(xmaxc,r,k1,k21)
  r3=r**3
  mr3=(1.0-r)**3
  leSlope=(3.0*r*r-k21*mr3-r3)*k1*cl/1.8

  CALL MeanLine3reflex(cl,xmaxc, x,ymean,ymeanp)
  CALL CombineThicknessAndCamber(x,yt,ymean,ymeanp, &
    xupper,yupper, xlower,ylower)
  CALL InterpolateCombinedAirfoil(x,yt,ymean,ymeanp, yup,ylo)
  CALL WriteOneSection(efu, name,id, x,yup, x,ylo, &
    xupper,yupper, xlower,ylower, rle,leSlope)
  RETURN
END Subroutine WriteFiveDigitReflexSection   ! ---------------------------------

!+
SUBROUTINE WriteSixSeriesSection(efu,family,a,cl,toc,name,id,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a six-series airfoil.
USE NACAauxilary,ONLY: Thickness6,MeanLine6,MeanLine6M,LeadingEdgeRadius6
USE NACAauxilary,ONLY: CombineThicknessAndCamber,InterpolateCombinedAirfoil

  INTEGER,INTENT(IN):: efu   ! external file unit
  INTEGER,INTENT(IN):: family  ! =1 63; =2 64; =3 65; =4 66; =5 67
                               ! =6 63A; =7 64A; =8 65A
  REAL,INTENT(IN):: a ! x-coor of aft-position of pressure rooftop
  REAL,INTENT(IN):: cl   ! design lift coefficient
  REAL,INTENT(IN):: toc   ! thickness / chord (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   !  section name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: leSlope
  REAL:: rle   ! leading edge radius, fraction of chord
  REAL,DIMENSION(SIZE(x)):: yt
  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
  REAL,DIMENSION(SIZE(x)):: xupper,yupper, xlower,ylower
  REAL,DIMENSION(SIZE(x)):: yup,ylo
!-------------------------------------------------------------------------------
  CALL Thickness6(family,toc, x,yt)
  IF (family < 6) THEN
    CALL MeanLine6(a,cl, x,ymean,ymeanp)
  ELSE
    CALL MeanLine6M(cl,x,ymean,ymeanp)
  END IF
  rle=LeadingEdgeRadius6(family,toc)
  leSlope=ymeanp(2)   ! we know the 2nd element from LoadX is 0.05   (MAGIC)
  CALL CombineThicknessAndCamber(x,yt,ymean,ymeanp, &
    xupper,yupper, xlower,ylower)
  CALL InterpolateCombinedAirfoil(x,yt,ymean,ymeanp, yup,ylo)
  CALL WriteOneSection(efu, name,id, x,yup, x,ylo, &
    xupper,yupper, xlower,ylower, rle,leSlope)
  RETURN
END Subroutine WriteSixSeriesSection   ! ---------------------------------------

!+
SUBROUTINE WriteOneMeanLine(efu,name,id,x,y,yp)
! ------------------------------------------------------------------------------
! PURPOSE - Write the description of one mean line to the HTML file.
  INTEGER,INTENT(IN):: efu   ! external file unit
  CHARACTER(LEN=*),INTENT(IN):: name   ! mean line name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x,y,yp   ! tables of x,y,dy/dz

!-------------------------------------------------------------------------------
  WRITE(efu,*) '<table>' 
  WRITE(efu,*) &
   '<tr><th colspan="3"><a id="'//id//'">NACA Mean Line </a>'//name//'</th></tr>'

  WRITE(efu,*) '<tr><td class="center" colspan="3">'
  WRITE(efu,*) '(Stations and ordinates given<br />'
  WRITE(efu,*) 'in per cent of airfoil chord)</td></tr>'

  WRITE(efu,*) '<tr><td class="center">x</td><td class="center">y</td>'// &
                   '<td class="center">dy/dx</td></tr>'

  CALL WriteColumns2(efu, 100.0*x,100.0*y,yp)
  WRITE(efu,*) "</table>"
  WRITE(efu,*) '<p>Go to the <a href="#TopMeanLines">top</a> of the page</p>'

  RETURN
END Subroutine WriteOneMeanLine   ! --------------------------------------------

!+
SUBROUTINE WriteOneProfile(efu,name,id,x,y,yp, rle)
! ------------------------------------------------------------------------------
! PURPOSE - Write the description of one thickness to the HTML file.
  INTEGER,INTENT(IN):: efu   ! external file unit
  CHARACTER(LEN=*),INTENT(IN):: name   ! profile name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: x,y,yp   ! tables of x,y,dy/dz
  REAL,INTENT(IN):: rle   ! leading edge radius, fraction of chord

  CHARACTER(LEN=*),PARAMETER:: FMT1 = &
    "('<tr><td colspan=""3"">L.E. radius = ', F5.3, ' percent chord</td></tr>')"
!-------------------------------------------------------------------------------
  WRITE(efu,*) '<table class="avd">'
  WRITE(efu,*) &
   '<tr><th colspan="3"><a id="'//id//'">NACA Profile </a>'//name//'</th></tr>'

  WRITE(efu,*) '<tr><td class="center" colspan="3">'
  WRITE(efu,*) '(Stations and ordinates given<br />'
  WRITE(efu,*) 'in per cent of airfoil chord)</td></tr>'
  WRITE(efu,*) '<tr><td class="center">x</td><td class="center">y</td>'// &
    '<td class="center">dy/dx</td></tr>'

  CALL WriteColumns2(efu, 100.0*x,100.0*y,yp)
  WRITE(efu,FMT1) 100.0*rle
  WRITE(efu,*) "</table>"
  WRITE(efu,*) '<p>Go to the <a href="#TopProfiles">top</a> of the page</p>'

  RETURN
END Subroutine WriteOneProfile   ! ---------------------------------------------

!+
SUBROUTINE WriteOneSection(efu,name,id,xup,yup,xlo,ylo, &
  xupper,yupper, xlower,ylower, rle,leSlope)
! ------------------------------------------------------------------------------
! PURPOSE - Write the description of one airfoil section to the HTML file.
! NOTE - Writes a table with one row and two columns. Each of the two cells
!   of this table contains a table with four columns and enough rows to
!   display xup, etc. along with column headers.

  INTEGER,INTENT(IN):: efu   ! external file unit
  CHARACTER(LEN=*),INTENT(IN):: name   ! section name
  CHARACTER(LEN=*),INTENT(IN):: id   ! HTML identification code
  REAL,INTENT(IN),DIMENSION(:):: xup,yup,xlo,ylo
  REAL,INTENT(IN),DIMENSION(:):: xupper,yupper,xlower,ylower
  REAL,INTENT(IN):: rle   ! leading edge radius, fraction of chord
  REAL,INTENT(IN):: leSlope   ! slope of mean line at leading edge
                              ! for 6-series mean lines, this is the slope
                              ! at 0.5 per cent chord


  CHARACTER(LEN=*),PARAMETER:: FMT1 = &
    "('<tr><td class=""center"" colspan=""4"">L.E. radius = ', F5.3, ' percent c</td></tr>')"
  CHARACTER(LEN=*),PARAMETER:: FMT2 = &
    "('<tr><td class=""center"" colspan=""4"">slope of mean line at LE = ', F7.4,'</td></tr>')"

!-------------------------------------------------------------------------------
  
  WRITE(efu,*) '<table>'     
  WRITE(efu,*) '<tr>'

  WRITE(efu,*) '<td style="border-style: none;">'
  WRITE(efu,*) '<table>'
  WRITE(efu,*) &
   '<tr><th colspan="4"><a id="'//id//'">NACA Section </a>'//name// '</th></tr>'

  WRITE(efu,*) '<tr><td class="center" colspan="4">'
  WRITE(efu,*) '(Stations and ordinates given<br />'
  WRITE(efu,*) 'in per cent of airfoil chord)</td></tr>'

  WRITE(efu,*) '<tr><td colspan="2" class="center">Upper Surface</td>'
  WRITE(efu,*)     '<td colspan="2" class="center">Lower Surface</td></tr>'

  WRITE(efu,*) &
   '<tr><td class="center">Station</td><td class="center">Ordinate</td>'
  WRITE(efu,*) &
   '<td class="center">Station</td><td class="center">Ordinate</td></tr>'

  CALL WriteColumns2(efu, 100.0*xup,100.0*yup,100.0*xlo,100.0*ylo)
  WRITE(efu,FMT1) 100.0*rle
  WRITE(efu,FMT2) leSlope
  WRITE(efu,*) "</table>"
  WRITE(efu,*) "</td>"

  WRITE(efu,*) '<td style="border-style: none;">'
  WRITE(efu,*) '<table>'
  WRITE(efu,*) '<tr><th colspan="4">NACA '//name//'</th></tr>'

  WRITE(efu,*) '<tr><td class="center" colspan="4">'
  WRITE(efu,*) '(Stations and ordinates given<br />'
  WRITE(efu,*) 'in per cent of airfoil chord)</td></tr>'


  WRITE(efu,*) '<tr><td class="center" colspan="2">Upper Surface</td>'
  WRITE(efu,*)     '<td class="center" colspan="2">Lower Surface</td></tr>'

  WRITE(efu,*) &
   '<tr><td class="center">Station</td><td class="center">Ordinate</td>'
  WRITE(efu,*) &
       '<td class="center">Station</td><td class="center">Ordinate</td></tr>'

  CALL WriteColumns2(efu, 100.0*xupper,100.0*yupper,100.0*xlower,100.0*ylower)
  WRITE(efu,FMT1) 100.0*rle
  WRITE(efu,FMT2) leSlope


  WRITE(efu,*) "</table>"
  WRITE(efu,*) "</td>"

  WRITE(efu,*) '</tr>'
  WRITE(efu,*) '</table>'
  WRITE(efu,*) '<p>Go to the <a href="#TopSections">top</a> of the page</p>'


  RETURN
END Subroutine WriteOneSection   ! ---------------------------------------------


!+
SUBROUTINE WriteColumns(efu,a,b,c,d,e,f)
! ------------------------------------------------------------------------------
! PURPOSE - Write the data in the arrays to the HTML file.
!   Print 4 figures after decimal point

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN),OPTIONAL,DIMENSION(:):: a,b,c,d,e,f

  CHARACTER(LEN=*),PARAMETER:: FMT='("<td>",F8.4,"</td>")'
  INTEGER:: k
!  REAL:: xx
!-------------------------------------------------------------------------------
  IF (.NOT.Present(a)) RETURN
  DO k=1,SIZE(a)
    WRITE(efu,*) '<tr>'
    WRITE(efu,FMT) a(k)
    IF (Present(b)) WRITE(efu,FMT) b(k)
    IF (Present(c)) WRITE(efu,FMT) c(k)
    IF (Present(d)) WRITE(efu,FMT) d(k)
    IF (Present(e)) WRITE(efu,FMT) e(k)
    IF (Present(f)) WRITE(efu,FMT) f(k)
    WRITE(efu,*) "</tr>"
  END DO

  RETURN
END Subroutine WriteColumns   ! ------------------------------------------------


!+
SUBROUTINE WriteColumns2(efu,a,b,c,d,e,f)
! ------------------------------------------------------------------------------
! PURPOSE - Write the data in the arrays to the HTML file.
!   Print 4 figures after decimal point

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN),OPTIONAL,DIMENSION(:):: a,b,c,d,e,f

  CHARACTER(LEN=8):: as,bs,cs,ds,es,fs
  CHARACTER(LEN=132):: buffer
!  CHARACTER(LEN=*),PARAMETER:: FMT='("<td>",F8.4,"</td>")'
  CHARACTER(LEN=*),PARAMETER:: FMT='(F8.4)'
  INTEGER:: k
!  REAL:: xx
!-------------------------------------------------------------------------------
  IF (.NOT.Present(a)) RETURN
  DO k=1,SIZE(a)
    WRITE(as,FMT) a(k)
    IF (Present(b)) WRITE(bs,FMT) b(k)
    IF (Present(c)) WRITE(cs,FMT) c(k)
    IF (Present(d)) WRITE(ds,FMT) d(k)
    IF (Present(e)) WRITE(es,FMT) e(k)
    IF (Present(f)) WRITE(fs,FMT) f(k)

    buffer='<tr><td>'//Trim(AdjustL(as))//'</td>'
    IF (Present(b)) buffer=Trim(buffer)//'<td>'//Trim(AdjustL(bs))//'</td>'
    IF (Present(c)) buffer=Trim(buffer)//'<td>'//Trim(AdjustL(cs))//'</td>'
    IF (Present(d)) buffer=Trim(buffer)//'<td>'//Trim(AdjustL(ds))//'</td>'
    IF (Present(e)) buffer=Trim(buffer)//'<td>'//Trim(AdjustL(es))//'</td>'
    IF (Present(f)) buffer=Trim(buffer)//'<td>'//Trim(AdjustL(fs))//'</td>'
    buffer=Trim(buffer)//'</tr>'
    WRITE(efu,'(A)') Trim(buffer)
  END DO

  RETURN
END Subroutine WriteColumns2   ! ------------------------------------------------

!+
SUBROUTINE WriteTextTable(efu,entries,ncol)
! ------------------------------------------------------------------------------
! PURPOSE - Write a table with ncol columns in HTML format.
!  Each cell has an <a> link and the name in the link is the
!  value of the entry prefaced with a #. This is the HTML convention.

  INTEGER,INTENT(IN):: efu   ! external file unit
  CHARACTER(LEN=*),INTENT(IN),DIMENSION(:):: entries ! array of strings
  INTEGER,INTENT(IN):: ncol  ! number of columns in the table


  CHARACTER(LEN=*),PARAMETER:: FMT = '("<td class=""links""><a href=",A,">",A,"</a></td>")'
  INTEGER:: k,k1,k2
  INTEGER:: n
  CHARACTER(LEN=132):: s
!-------------------------------------------------------------------------------
  WRITE(efu,*) '<table>'//'   <!-- links to positions on this page-->'

  n=SIZE(entries)
  k2=0
  DO
    WRITE(efu,*) '<tr>'
    k1=k2+1
    k2=MIN(k2+ncol,n)
    DO k=k1,k2
      s='"#' // Trim(AdjustL(entries(k))) // '"'
      WRITE(efu,FMT) Trim(s),Trim(entries(k))
    END DO
    WRITE(efu,*) '</tr>'
    IF (k2==n) EXIT
  END DO
  WRITE(efu,*) '</table>'
  RETURN
END Subroutine WriteTextTable   ! ----------------------------------------------

!+
SUBROUTINE WriteTextTableRows(efu,header,entries)
! ------------------------------------------------------------------------------
! PURPOSE - Write two rows of a table in HTML format.
!  The number of columns in the table is the SIZE of entries.
!  The first line is a header that spans all the columns.
!  If the header has trimmed length of zero, don't write it.
!  Each cell has an <a> link and the name in the link is the
!  value of the entry prefaced with a #. This is the HTML convention.

  INTEGER,INTENT(IN):: efu   ! external file unit
  CHARACTER(LEN=*),INTENT(IN):: header
  CHARACTER(LEN=*),INTENT(IN),DIMENSION(:):: entries


  CHARACTER(LEN=*),PARAMETER:: FMT1 = &
   '("<tr><td class=""links"" colspan=",I1,">",A,"</td></tr>")'
  CHARACTER(LEN=*),PARAMETER:: FMT2 = &
   '("<td class=""links""><a href=",A,"</a></td>")'
  INTEGER:: k
  CHARACTER(LEN=132):: s
!-------------------------------------------------------------------------------
!  WRITE(efu,*) '<table cellpadding=2 '// &
!    'summary="links to positions on this page">'
!!!  WRITE(HTML,*) '<tr><td colspan=8>Four digit profiles</td></tr>'

  k=MAX(2, SIZE(entries)) 
  IF (LEN_TRIM(header) > 0) WRITE(efu,FMT1) k, header
  WRITE(efu,*) '<tr>'
  DO k=1,SIZE(entries)
    s=AdjustL(entries(k))
    s='"#' // Trim(s) // '">' // Trim(s) 
    WRITE(efu,FMT2) Trim(s)
  END DO
  WRITE(efu,*) '</tr>'
  RETURN
END Subroutine WriteTextTableRows   ! ------------------------------------------

!+
SUBROUTINE WriteBanner(efu)
! ------------------------------------------------------------------------------
! PURPOSE - Write the PDAS banner into the HTML file
  INTEGER,INTENT(IN):: efu
!-------------------------------------------------------------------------------
  WRITE(efu,*) '<div class="newbanner">'
  WRITE(efu,*) 'Public Domain Aeronautical Software (PDAS) &nbsp;</div>'  
  RETURN
END Subroutine WriteBanner   ! -------------------------------------------------




!+
SUBROUTINE WriteCrumb(efu, crumbTitle)
! ------------------------------------------------------------------------------
! PURPOSE- Write the crumb line showing the path to the current page
  INTEGER,INTENT(IN):: efu
  CHARACTER(LEN=*),INTENT(IN):: crumbTitle
!-------------------------------------------------------------------------------
  WRITE(efu,*) '<nav class="crumb">'
  WRITE(efu,*) '<a href="index.html">PDAS home</a> &gt; '
  WRITE(efu,*) '<a href="contents16.html">Contents</a> &gt; '
  WRITE(efu,*) '<a href="naca456.html">NACA Airfoil</a> &gt; '
  WRITE(efu,*) '<a href="avd.html">AVD Appendices</a> &gt; '
  WRITE(efu,*) Trim(crumbTitle)
  WRITE(efu,*) '</nav>'
  RETURN
END Subroutine WriteCrumb   ! --------------------------------------------------

!+
SUBROUTINE WriteFooter(efu)
! ------------------------------------------------------------------------------
! PURPOSE - Write the footer, picking up address1 that contains the cuurent date
  INTEGER,INTENT(IN):: efu   ! external file unit for writing output
  CHARACTER(LEN=*),PARAMETER:: DASH = "-"
  CHARACTER(LEN=*),PARAMETER:: MONTH='JanFebMarAprMayJunJulAugSepOctNovDec'
  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(I4,A1,I2.2,A1,I2.2)'
  CHARACTER(LEN=*),PARAMETER:: FMT2 = '(I4,1X,A3,1X,I2.0)'
  CHARACTER(LEN=15):: s1,s2
  CHARACTER(LEN=*),PARAMETER:: T1 = '<a href="order.html">Get a copy</a>' // &
    ' of all the programs from'
  CHARACTER(LEN=*),PARAMETER:: T2 =   &
   '<cite>Public Domain Computer Programs for the Aeronautical Engineer.</cite><br />'
  INTEGER,DIMENSION(8):: v

!-------------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)
  WRITE(s1,FMT1) v(1),DASH,v(2),DASH,v(3)    ! current date in ISO 8601 format
  WRITE(s2,FMT2) v(1), MONTH(3*v(2)-2:3*v(2)), v(3)   ! human readable 
  WRITE(efu,*) '<footer>'
  WRITE(efu,*) T1
  WRITE(efu,*) T2

   
  WRITE(efu,'(5A)') &
    'Last updated: <time datetime="', TRIM(s1), '">',TRIM(s2), '</time>'
  WRITE(efu,*) 'by Ralph Carmichael, pdaerowebmaster AT gmail DOT com.'
  WRITE(efu,*) '</footer>'
  RETURN
END Subroutine WriteFooter   ! -------------------------------------------------

!+
SUBROUTINE WriteHead(efu,filename,title)
! ------------------------------------------------------------------------------
! PURPOSE - 
  INTEGER,INTENT(IN):: efu
  CHARACTER(LEN=*),INTENT(IN):: filename
  CHARACTER(LEN=*),INTENT(IN):: title
!-------------------------------------------------------------------------------
  WRITE(efu,*) '<!DOCTYPE html>'
  WRITE(efu,*) '<html lang="en-us">'
  WRITE(efu,*) '<!-- www.pdas.com/'//filename//'   called from avd page -->'
  WRITE(efu,*) '<head>'
  WRITE(efu,*) '<title>'//title//'</title>'
  WRITE(efu,*) '<meta charset="utf-8" />'
  WRITE(efu,*) '<meta name="generator" content="avd program (Fortran)" />'
  WRITE(efu,*) '<meta name="author" content="Ralph Carmichael" />'
  WRITE(efu,*) '<meta name="keywords" content= "NACA airfoil, AIAA-2001-5235,'// &
    '  profiles, meanlines, public domain, PDAS" />'
  WRITE(efu,*) '<link rel="stylesheet" href="pdas.css" type="text/css" />'
  WRITE(efu,*) '<link rel="stylesheet" href="pdasp.css" type="text/css" media="print" />'

  WRITE(efu,*) '<style>'
!  WRITE(efu,*) 'body {margin: 5px; color: black; background: white;'
!  WRITE(efu,*) '  font-family: Arial,helvetica,sans-serif;}'

  WRITE(efu,*) 'th {border: 1px solid black; text-align: center; ' // &
     'border-collapse: collapse; }'
  WRITE(efu,*) 'td {border: 1px solid black; text-align: right; ' // &
     'border-collapse: collapse; }'
!  WRITE(efu,*) 'td.center {text-align: center; }'
!  WRITE(efu,*) 'td.links { padding: 2px; border-style: none}'
!  WRITE(efu,*) 'h1 { font-family: arial,helvetica,univers,sans-serif;'
!  WRITE(efu,*) '     font-size: 2em; font-weight: bold; color: blue;'
!  WRITE(efu,*) '     text-align: center; margin-top: 0.5em; }'

  WRITE(efu,*) '</style></head>'
  
  RETURN
END Subroutine WriteHead   ! ---------------------------------------------------



END Module AbbottVonDoenhoffSubs   ! ===========================================

!+
MODULE MainProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Collect 

USE AbbottVonDoenhoffSubs
IMPLICIT NONE

  INTEGER,PARAMETER:: HTML=1

!  CHARACTER(LEN=*),PARAMETER:: ADDRESS2 = &
!    ' by Ralph Carmichael <a href="mailto:pdaerowebmaster@gmail.com">pdaerowebmaster AT gmail DOT com</a>'



  CHARACTER(LEN=*),PARAMETER:: PAUSENOTE = &
    'Let the page finish loading before clicking on a section.<br />'  

  PUBLIC:: MainPage
  PUBLIC:: Profiles
  PUBLIC:: Meanlines
  PUBLIC:: Sections45
  PUBLIC:: Sections6
  PUBLIC:: Sections6A
!  PUBLIC:: Welcome
!-------------------------------------------------------------------------------


CONTAINS

!+
SUBROUTINE MainPage()
! ------------------------------------------------------------------------------
! PURPOSE - Write the main page with links to the other data pages 
!  There are links to profiles.html, meanline.html, sections45.html,
!  sections6.html, sections6a.html

!-------------------------------------------------------------------------------
  OPEN(UNIT=HTML, FILE='avd.html', STATUS='REPLACE', ACTION='WRITE')
  CALL WriteHead(HTML, 'avd.html', 'Airfoil Tables')
  WRITE(HTML,*) '<body>'
  WRITE(HTML,*) '<div class="crumb">'
  WRITE(HTML,*) '<a href="index.html">PDAS home</a> &gt; '
  WRITE(HTML,*) '<a href="contents16.html">Contents</a> &gt; '
  WRITE(HTML,*) '<a href="naca456.html">NACA Airfoil</a> &gt; '
  WRITE(HTML,*) 'AVD Appendices </div>'

!  WRITE(HTML,*) &
!    '<div class="newbanner">Public Domain Aeronautical Software (PDAS) &nbsp;</div>'  
  CALL WriteBanner(HTML)
  WRITE(HTML,*) '<div id="header"><h1>Tables from Appendices I,II,III of '
  WRITE(HTML,*) '<cite>Theory of Airfoil Sections</cite></h1></div>'
  WRITE(HTML,*) '<p>The referenced book has three large appendices'
  WRITE(HTML,*) 'containing tables of coordinates of a selection of NACA'
  WRITE(HTML,*) 'profiles, mean lines and sections.</p>'
  WRITE(HTML,*) '<p>This page and the linked pages reproduces these data'
  WRITE(HTML,*) 'recomputed using the procedures of the PDAS program naca456'
  WRITE(HTML,*) 'module.</p>'
  WRITE(HTML,*) '<ul>'
  WRITE(HTML,*) '<li><a href="profiles.html">Appendix I</a>- Profiles</li>'
  WRITE(HTML,*) '<li><a href="meanlines.html">Appendix II</a>- Mean Lines</li>'
  WRITE(HTML,*) '<li><a href="sections45.html">' // &
    'Appendix III</a>- 4 and 5 digit Sections</li>'
  WRITE(HTML,*) '<li><a href="sections6.html">' // &
    'Appendix III</a>- 6-Series Sections</li>'
  WRITE(HTML,*) '<li><a href="sections6a.html">' // &
    'Appendix III</a>- 6A-Series Sections</li>'
  WRITE(HTML,*) '</ul>'

  CALL WriteFooter(HTML)
  WRITE(HTML,*) '<nav>'
  WRITE(HTML,*) '<a href="index.html">PDAS home</a> &gt; '
  WRITE(HTML,*) '<a href="contents16.html">Contents</a> &gt; '
  WRITE(HTML,*) '<a href="naca456.html">NACA Airfoil</a> &gt; '
  WRITE(HTML,*) 'AVD Appendices </nav>'

  CALL WriteBanner(HTML)
  WRITE(HTML,*) " </body></html>"  
  CLOSE(UNIT=HTML)

  RETURN
END Subroutine MainPage   ! ----------------------------------------------------

!+
SUBROUTINE Profiles(x)
! ------------------------------------------------------------------------------
! PURPOSE - Create a HTML file with the tables of profiles - Appendix I
USE AbbottVonDoenhoffSubs

  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coor for calculation

!  CHARACTER(LEN=*),PARAMETER:: DESCRIPTION = &
!    '<meta name="description" content= "NACA Airfoil Profiles" />'

!-------------------------------------------------------------------------------
  OPEN(UNIT=HTML,FILE='profiles.html',STATUS='REPLACE',ACTION='WRITE')
  CALL WriteHead(HTML, 'profiles.html', 'Appendix I - Tables of Profiles')

  WRITE(HTML,*) '<body class="profiles">'
  WRITE(HTML,*) '<p><a id="TopProfiles"></a></p>'
  CALL WriteCrumb(HTML, 'Profiles')

  CALL WriteBanner(HTML)
  WRITE(HTML,*) '<div id="header"><h1>Appendix I - Profiles</h1></div>'

  WRITE(HTML,*) '<p><a href="avd.html">Go back to the AVD main page</a><br />'
!  WRITE(HTML,*) PAUSENOTE
  WRITE(HTML,*) '</p>'

  WRITE(HTML,*) '<p>Four digit profiles<br />'
!  CALL WriteTextTable(HTML, (/  &
!   "0006", "0008", "0010", "0012", "0015", "0018", "0021" /), 7)
  WRITE(HTML,*) '<a href="#p0006">0006</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p0008">0008</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p0010">0010</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p0012">0012</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p0015">0015</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p0018">0018</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p0021">0021</a></p>'
    
  
  WRITE(HTML,*) '<p>Modified four digit profiles<br />'
!  CALL WriteTextTable(HTML, (/ &
!   "0008-34", "0010-34", "0010-35", "0010-64",           &
!   "0010-66", "0010-66", "0012-34", "0012-64" /), 4)
  WRITE(HTML,*) '<a href="#p000834">0008-34</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p001034">0010-34</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p001035">0010-35</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p001064">0010-64</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p001065">0010-65</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p001066">0010-66</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p001234">0012-34</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p001264">0012-64</a></p>'

  
  WRITE(HTML,*) '<p>16-series profiles<br />'
!  CALL WriteTextTable(HTML, (/ &
!   "16-006", "16-009", "16-012", "16-015", "16-018", "16-021" /), 6)
  WRITE(HTML,*) '<a href="#p16006">16-006</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p16009">16-009</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p16012">16-012</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p16015">16-015</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p16018">16-018</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p16021">16-021</a></p>'
   
   
  WRITE(HTML,*) '<p>63 series profiles<br />'
!  CALL WriteTextTable(HTML,(/ &
!   "63-006", "63-008", "63-009", "63-010",               &
!   "63-012", "63-015", "63-018", "63-021"/), 4)
  WRITE(HTML,*) '<a href="#p63006">63-006</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p63008">63-008</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p63009">63-009</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p63010">63-010</a> &nbsp;<br />'
  WRITE(HTML,*) '<a href="#p63012">63-012</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p63015">63-015</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p63018">63-018</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p63021">63-021</a></p>'
    
  WRITE(HTML,*) '<p>63A series profiles<br />'
!  CALL WriteTextTable(HTML,  (/ &
!   "63A006", "63A008", "63A010", "63A012", "63A015" /), 5) 
  WRITE(HTML,*) '<a href="#p63A006">63A-006</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p63A008">63A-000</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p63A010">63A-010</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p63A012">63A-012</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p63A015">63A-015</a></p>'

  
  WRITE(HTML,*) '<p>64 series profiles<br />'
!  CALL WriteTextTable(HTML, (/ &
!   "64-006", "64-008", "64-009", "64-010",               &
!   "64-012", "64-015", "64-018", "64-021"/), 4)
  WRITE(HTML,*) '<a href="#p64006">64-006</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p64008">64-008</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p64009">64-009</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p64010">64-010</a> &nbsp;<br />'
  WRITE(HTML,*) '<a href="#p64012">64-012</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p64015">64-015</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p64018">64-018</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p64021">64-021</a></p>'

  WRITE(HTML,*) '<p>64A series profiles<br />'
!  CALL WriteTextTable(HTML, (/ &
!   "64A006", "64A008", "64A010", "64A012", "64A015" /), 5) 
  WRITE(HTML,*) '<a href="#p64A006">64A-006</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p64A008">64A-008</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p64A010">64A-010</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p64A012">64A-012</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p64A015">64A-015</a></p>'
   
  
  WRITE(HTML,*) '<p>65 series profiles<br />'
!  CALL WriteTextTable(HTML, (/ &
!   "65-006", "65-008", "65-009", "65-010",               &
!   "65-012", "65-015", "65-018", "65-021"/), 4)
  WRITE(HTML,*) '<a href="#p65006">65-006</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p65008">65-008</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p65009">65-009</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p65010">65-010</a> &nbsp;<br />'
  WRITE(HTML,*) '<a href="#p65012">65-012</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p65015">65-015</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p65018">65-018</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p65021">65-021</a></p>'

  WRITE(HTML,*) '<p>65A series profiles<br />'
!  CALL WriteTextTable(HTML,(/ &
!   "65A006", "65A008", "65A010", "65A012", "65A015" /), 5) 
  WRITE(HTML,*) '<a href="#p65A006">65A-006</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p65A008">65A-008</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p65A010">65A-010</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p65A012">65A-012</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p65A015">65A-015</a></p>'
  
  WRITE(HTML,*) '<p>66 series profiles<br />'
!  CALL WriteTextTable(HTML, (/ &
!   "66-006", "66-008", "66-009", "66-010",               &
!   "66-012", "66-015", "66-018", "66-021"/), 4)
  WRITE(HTML,*) '<a href="#p66006">66-006</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p66008">66-008</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p66009">66-009</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p66010">66-010</a> &nbsp;<br />'
  WRITE(HTML,*) '<a href="#p66012">66-012</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p66015">66-015</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p66018">66-018</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p66021">66-021</a></p>'

  WRITE(HTML,*) '<p>67 series profiles<br />'
!  CALL WriteTextTable(HTML, (/"67-012", "67-015" /), 2)
  WRITE(HTML,*) '<a href="#p67012">67-012</a> &nbsp;'
  WRITE(HTML,*) '<a href="#p67015">67-015</a></p>'

  
  CALL WriteFourDigitProfile(HTML, 0.06, '0006','p0006', x)
  CALL WriteFourDigitProfile(HTML, 0.08, '0008','p0008', x)
  CALL WriteFourDigitProfile(HTML, 0.10, '0010','p0010', x)
  CALL WriteFourDigitProfile(HTML, 0.12, '0012','p0012', x)
  CALL WriteFourDigitProfile(HTML, 0.15, '0015','p0015', x)
  CALL WriteFourDigitProfile(HTML, 0.18, '0018','p0018', x)
  CALL WriteFourDigitProfile(HTML, 0.21, '0021','p0021', x)

  CALL WriteModifiedFourDigitProfile(HTML, 3.0, 0.4, 0.08, '0008-34','p000834', x)
  CALL WriteModifiedFourDigitProfile(HTML, 3.0, 0.4, 0.10, '0010-34','p001034', x)
  CALL WriteModifiedFourDigitProfile(HTML, 3.0, 0.5, 0.10, '0010-35','p001035', x)
  CALL WriteModifiedFourDigitProfile(HTML, 6.0, 0.4, 0.10, '0010-64','p001064', x)
  CALL WriteModifiedFourDigitProfile(HTML, 6.0, 0.5, 0.10, '0010-65','p001065', x)
  CALL WriteModifiedFourDigitProfile(HTML, 6.0, 0.6, 0.10, '0010-66','p001066', x)
  CALL WriteModifiedFourDigitProfile(HTML, 3.0, 0.4, 0.12, '0012-34','p001234', x)
  CALL WriteModifiedFourDigitProfile(HTML, 6.0, 0.4, 0.12, '0012-64','p001264', x)

  CALL WriteModifiedFourDigitProfile(HTML, 4.0, 0.5, 0.06, '16-006','p16006', x)
  CALL WriteModifiedFourDigitProfile(HTML, 4.0, 0.5, 0.09, '16-009','p16009', x)
  CALL WriteModifiedFourDigitProfile(HTML, 4.0, 0.5, 0.12, '16-012','p16012', x)
  CALL WriteModifiedFourDigitProfile(HTML, 4.0, 0.5, 0.15, '16-015','p16015', x)
  CALL WriteModifiedFourDigitProfile(HTML, 4.0, 0.5, 0.18, '16-018','p16018', x)
  CALL WriteModifiedFourDigitProfile(HTML, 4.0, 0.5, 0.21, '16-021','p16021', x)

  CALL WriteSixSeriesProfile(HTML, 1, 0.06, '63-006','p63006', x)
  CALL WriteSixSeriesProfile(HTML, 1, 0.08, '63-008','p63008', x)
  CALL WriteSixSeriesProfile(HTML, 1, 0.09, '63-009','p63009', x)
  CALL WriteSixSeriesProfile(HTML, 1, 0.10, '63-010','p63010', x)
  CALL WriteSixSeriesProfile(HTML, 1, 0.12, '63-012','p63012', x)
  CALL WriteSixSeriesProfile(HTML, 1, 0.15, '63-015','p63015', x)
  CALL WriteSixSeriesProfile(HTML, 1, 0.18, '63-018','p63018', x)
  CALL WriteSixSeriesProfile(HTML, 1, 0.21, '63-021','p63021', x)

  CALL WriteSixSeriesProfile(HTML, 6, 0.06, '63A006','p63A006', x)
  CALL WriteSixSeriesProfile(HTML, 6, 0.08, '63A008','p63A008', x)
  CALL WriteSixSeriesProfile(HTML, 6, 0.10, '63A010','p63A010', x)
  CALL WriteSixSeriesProfile(HTML, 6, 0.12, '63A012','p63A012', x)
  CALL WriteSixSeriesProfile(HTML, 6, 0.15, '63A015','p63A015', x)

  CALL WriteSixSeriesProfile(HTML, 2, 0.06, '64-006','p64006', x)
  CALL WriteSixSeriesProfile(HTML, 2, 0.08, '64-008','p64008', x)
  CALL WriteSixSeriesProfile(HTML, 2, 0.09, '64-009','p64009', x)
  CALL WriteSixSeriesProfile(HTML, 2, 0.10, '64-010','p64010', x)
  CALL WriteSixSeriesProfile(HTML, 2, 0.12, '64-012','p64012', x)
  CALL WriteSixSeriesProfile(HTML, 2, 0.15, '64-015','p64015', x)
  CALL WriteSixSeriesProfile(HTML, 2, 0.18, '64-018','p64018', x)
  CALL WriteSixSeriesProfile(HTML, 2, 0.21, '64-021','p64021', x)

  CALL WriteSixSeriesProfile(HTML, 7, 0.06, '64A006','p64A006', x)
  CALL WriteSixSeriesProfile(HTML, 7, 0.08, '64A008','p64A008', x)
  CALL WriteSixSeriesProfile(HTML, 7, 0.10, '64A010','p64A010', x)
  CALL WriteSixSeriesProfile(HTML, 7, 0.12, '64A012','p64A012', x)
  CALL WriteSixSeriesProfile(HTML, 7, 0.15, '64A015','p64A015', x)

  CALL WriteSixSeriesProfile(HTML, 3, 0.06, '65-006','p65006', x)
  CALL WriteSixSeriesProfile(HTML, 3, 0.08, '65-008','p65008', x)
  CALL WriteSixSeriesProfile(HTML, 3, 0.09, '65-009','p65009', x)
  CALL WriteSixSeriesProfile(HTML, 3, 0.10, '65-010','p65010', x)
  CALL WriteSixSeriesProfile(HTML, 3, 0.12, '65-012','p65012', x)
  CALL WriteSixSeriesProfile(HTML, 3, 0.15, '65-015','p65015', x)
  CALL WriteSixSeriesProfile(HTML, 3, 0.18, '65-018','p65018', x)
  CALL WriteSixSeriesProfile(HTML, 3, 0.21, '65-021','p65021', x)

  CALL WriteSixSeriesProfile(HTML, 8, 0.06, '65A006','p65A006', x)
  CALL WriteSixSeriesProfile(HTML, 8, 0.08, '65A008','p65A008', x)
  CALL WriteSixSeriesProfile(HTML, 8, 0.10, '65A010','p65A010', x)
  CALL WriteSixSeriesProfile(HTML, 8, 0.12, '65A012','p65A012', x)
  CALL WriteSixSeriesProfile(HTML, 8, 0.15, '65A015','p65A015', x)

  CALL WriteSixSeriesProfile(HTML, 4, 0.06, '66-006','p66006', x)
  CALL WriteSixSeriesProfile(HTML, 4, 0.08, '66-008','p66008', x)
  CALL WriteSixSeriesProfile(HTML, 4, 0.09, '66-009','p66009', x)
  CALL WriteSixSeriesProfile(HTML, 4, 0.10, '66-010','p66010', x)
  CALL WriteSixSeriesProfile(HTML, 4, 0.12, '66-012','p66012', x)
  CALL WriteSixSeriesProfile(HTML, 4, 0.15, '66-015','p66015', x)
  CALL WriteSixSeriesProfile(HTML, 4, 0.18, '66-018','p66018', x)
  CALL WriteSixSeriesProfile(HTML, 4, 0.21, '66-021','p66021', x)

  CALL WriteSixSeriesProfile(HTML, 5, 0.12, '67-012','p67012', x)
  CALL WriteSixSeriesProfile(HTML, 5, 0.15, '67-015','p67015', x)

  WRITE(HTML,*) '<div>Go back to the <a href="avd.html">AVD</a> page.</div>'

  CALL WriteFooter(HTML)
  CALL WriteCrumb(HTML, 'Profiles')
  CALL WriteBanner(HTML)
  WRITE(HTML,*) '</body></html>'
  CLOSE(UNIT=HTML)
  RETURN
END Subroutine Profiles   ! ----------------------------------------------------

!+
SUBROUTINE MeanLines(x)
! ------------------------------------------------------------------------------
! PURPOSE - Create a HTML file with the tables of mean lines - Appendix II
USE AbbottVonDoenhoffSubs
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coor for calculation
!-------------------------------------------------------------------------------
  OPEN(UNIT=HTML,FILE='meanlines.html',STATUS='REPLACE',ACTION='WRITE')
  CALL WriteHead(HTML, 'meanlines.html', 'Appendix II - Mean Lines')
 
  WRITE(HTML,*) '<body class="meanlines">'
  WRITE(HTML,*) '<p id="TopMeanLines"></p>'
  CALL WriteCrumb(HTML, 'Meanlines')
  CALL WriteBanner(HTML)
  WRITE(HTML,*) ''
  WRITE(HTML,*) '<div id="header">'
  WRITE(HTML,*) '<h1>Appendix II - Mean Lines</h1></div>'

  WRITE(HTML,*) '<p><a href="avd.html">Go back to the AVD main page</a><br />'
  WRITE(HTML,*) PAUSENOTE
  WRITE(HTML,*) '</p>'
 
  WRITE(HTML,*) "<p>Two digit mean lines<br />"
  WRITE(HTML,*) '<a href="#ml62">62</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml63">63</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml64">64</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml65">65</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml66">66</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml67">67</a></p>'

  WRITE(HTML,*) "<p>Three digit mean lines<br />"
  WRITE(HTML,*) '<a href="#ml210">210</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml220">220</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml230">230</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml240">240</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml250">250</a></p>'

  WRITE(HTML,*) "<p>Three digit reflex mean lines<br />"
  WRITE(HTML,*) '<a href="#ml211">211</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml221">221</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml231">231</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml241">241</a> &nbsp;'
  WRITE(HTML,*) '<a href="#ml251">251</a></p>'

  WRITE(HTML,*) "<p>Six-series mean lines<br />"
  WRITE(HTML,*) '<a href="#mla0">a=0.0</a> &nbsp;'
  WRITE(HTML,*) '<a href="#mla1">a=0.1</a> &nbsp;'
  WRITE(HTML,*) '<a href="#mla2">a=0.2</a> &nbsp;'
  WRITE(HTML,*) '<a href="#mla3">a=0.3</a> &nbsp;'
  WRITE(HTML,*) '<a href="#mla4">a=0.4</a> &nbsp;'
  WRITE(HTML,*) '<a href="#mla5">a=0.5</a><br />'
  WRITE(HTML,*) '<a href="#mla6">a=0.6</a> &nbsp;'
  WRITE(HTML,*) '<a href="#mla7">a=0.7</a> &nbsp;'
  WRITE(HTML,*) '<a href="#mla8">a=0.8</a> &nbsp;'
  WRITE(HTML,*) '<a href="#mla9">a=0.9</a> &nbsp;'
  WRITE(HTML,*) '<a href="#mla10">a=1.0</a></p>'
    
  WRITE(HTML,*) "<p>Six-series modified mean line<br />"
  WRITE(HTML,*) '<a href="#mla8m">a=0.8 (modified)</a></p>'

  CALL WriteTwoDigitMeanLine(HTML, 0.06, 0.2, '62','ml62', x)
  CALL WriteTwoDigitMeanLine(HTML, 0.06, 0.3, '63','ml63', x)
  CALL WriteTwoDigitMeanLine(HTML, 0.06, 0.4, '64','ml64', x)
  CALL WriteTwoDigitMeanLine(HTML, 0.06, 0.5, '65','ml65', x)
  CALL WriteTwoDigitMeanLine(HTML, 0.06, 0.6, '66','ml66', x)
  CALL WriteTwoDigitMeanLine(HTML, 0.06, 0.7, '67','ml67', x)

  CALL WriteThreeDigitMeanLine(HTML, 0.3, 0.05, '210','ml210', x)
  CALL WriteThreeDigitMeanLine(HTML, 0.3, 0.10, '220','ml220', x)
  CALL WriteThreeDigitMeanLine(HTML, 0.3, 0.15, '230','ml230', x)
  CALL WriteThreeDigitMeanLine(HTML, 0.3, 0.20, '240','ml240', x)
  CALL WriteThreeDigitMeanLine(HTML, 0.3, 0.25, '250','ml250', x)
  CALL WriteThreeDigitMeanLine(HTML, 0.6, 0.15, '430','ml430', x)

  CALL WriteThreeDigitReflexMeanLine(HTML, 0.3, 0.05, '211','ml211', x)
  CALL WriteThreeDigitReflexMeanLine(HTML, 0.3, 0.10, '221','ml221', x)
  CALL WriteThreeDigitReflexMeanLine(HTML, 0.3, 0.15, '231','ml231', x)
  CALL WriteThreeDigitReflexMeanLine(HTML, 0.3, 0.20, '241','ml241', x)
  CALL WriteThreeDigitReflexMeanLine(HTML, 0.3, 0.25, '251','ml251', x)

  CALL WriteSixSeriesMeanLine(HTML,0.0,1.0,'a=0.0','mla0',x)
  CALL WriteSixSeriesMeanLine(HTML,0.1,1.0,'a=0.1','mla1',x)
  CALL WriteSixSeriesMeanLine(HTML,0.2,1.0,'a=0.2','mla2',x)
  CALL WriteSixSeriesMeanLine(HTML,0.3,1.0,'a=0.3','mla3',x)
  CALL WriteSixSeriesMeanLine(HTML,0.4,1.0,'a=0.4','mla4',x)
  CALL WriteSixSeriesMeanLine(HTML,0.5,1.0,'a=0.5','mla5',x)
  CALL WriteSixSeriesMeanLine(HTML,0.6,1.0,'a=0.6','mla6',x)
  CALL WriteSixSeriesMeanLine(HTML,0.7,1.0,'a=0.7','mla7',x)
  CALL WriteSixSeriesMeanLine(HTML,0.8,1.0,'a=0.8','mla8',x)
  CALL WriteSixSeriesMeanLine(HTML,0.9,1.0,'a=0.9','mla9',x)
  CALL WriteSixSeriesMeanLine(HTML,1.0,1.0,'a=1.0','mla10',x)

  CALL WriteSixSeriesModifiedMeanLine(HTML,1.0,'a=0.8 (modified)','mla8m',x)

  WRITE(HTML,*) '<div>Go back to the <a href="avd.html">AVD</a> page.</div>'

  CALL WriteFooter(HTML)
  CALL WriteCrumb(HTML, 'Meanlines')
  CALL WriteBanner(HTML)
  WRITE(HTML,*) '</body></html>'
  CLOSE(UNIT=HTML)
  RETURN
END Subroutine MeanLines   ! ---------------------------------------------------

!+
SUBROUTINE Sections45(x)
! ------------------------------------------------------------------------------
! PURPOSE - Create a HTML file with the tables of profiles - Appendix I
!   This subroutine writes the 4 and 5 digit sections.
USE AbbottVonDoenhoffSubs
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coor for calculation
!-------------------------------------------------------------------------------

  OPEN(UNIT=HTML,FILE='sections45.html',STATUS='REPLACE',ACTION='WRITE')
  CALL WriteHead(HTML, 'sections45.html', 'Appendix III - 4 and 5 Digit Sections')

  WRITE(HTML,*) '<body class="sections45">'
  WRITE(HTML,*) '<p id="TopSections"></p>'
  CALL WriteCrumb(HTML, '4- and 5-Digit Sections')
  CALL WriteBanner(HTML)
  WRITE(HTML,*) '<div id="header"><h1>Appendix III - 4 and 5 Digit Sections</h1></div>'

  WRITE(HTML,*) '<p><a href="avd.html">Go back to the AVD main page</a><br />'
  WRITE(HTML,*) PAUSENOTE
  WRITE(HTML,*) '</p>'

  WRITE(HTML,*) '<p>Four digit cambered sections<br />'
  WRITE(HTML,*) '<a href="#s1408">1408</a>'
  WRITE(HTML,*) '<a href="#s1410">1410</a>'
  WRITE(HTML,*) '<a href="#s1412">1412</a>'
  WRITE(HTML,*) '<a href="#s2408">2408</a>'
  WRITE(HTML,*) '<a href="#s2410">2410</a><br />'
  WRITE(HTML,*) '<a href="#s2412">2412</a>'
  WRITE(HTML,*) '<a href="#s2415">2415</a>'
  WRITE(HTML,*) '<a href="#s2418">2418</a>'
  WRITE(HTML,*) '<a href="#s2421">2421</a>'
  WRITE(HTML,*) '<a href="#s2424">2424</a><br />'
  WRITE(HTML,*) '<a href="#s4412">4412</a>'
  WRITE(HTML,*) '<a href="#s4415">4415</a>'
  WRITE(HTML,*) '<a href="#s4418">4418</a>'
  WRITE(HTML,*) '<a href="#s4421">4421</a>'
  WRITE(HTML,*) '<a href="#s4424">4424</a></p>'

!  CALL WriteTextTable(HTML, (/                                          &
!    "1408", "1410", "1412", "2408", "2410",                             &
!    "2412", "2415", "2418", "2421", "2424",                             &
!    "4412", "4415", "4418", "4421", "4424" /), 5)

  WRITE(HTML,*) '<p>Five-digit cambered sections<br />'
!  CALL WriteTextTable(HTML, (/ &
!    "23012", "23015", "23018", "23021", "23024" /), 5)
  WRITE(HTML,*) '<a href="#s23012">23012</a>'
  WRITE(HTML,*) '<a href="#s23015">23015</a>'
  WRITE(HTML,*) '<a href="#s23018">23018</a>'
  WRITE(HTML,*) '<a href="#s23021">23021</a>'
  WRITE(HTML,*) '<a href="#s23024">23024</a></p>'
  

  CALL WriteFourDigitSection(HTML, 0.01,0.4,0.08, '1408','s1408', x)
  CALL WriteFourDigitSection(HTML, 0.01,0.4,0.10, '1410','s1410', x)
  CALL WriteFourDigitSection(HTML, 0.01,0.4,0.12, '1412','s1412', x)
  CALL WriteFourDigitSection(HTML, 0.02,0.4,0.08, '2408','s2408', x)
  CALL WriteFourDigitSection(HTML, 0.02,0.4,0.10, '2410','s2410', x)
  CALL WriteFourDigitSection(HTML, 0.02,0.4,0.12, '2412','s2412', x)
  CALL WriteFourDigitSection(HTML, 0.02,0.4,0.15, '2415','s2415', x)
  CALL WriteFourDigitSection(HTML, 0.02,0.4,0.18, '2418','s2418', x)
  CALL WriteFourDigitSection(HTML, 0.02,0.4,0.21, '2421','s2421', x)
  CALL WriteFourDigitSection(HTML, 0.02,0.4,0.24, '2424','s2424', x)
  CALL WriteFourDigitSection(HTML, 0.04,0.4,0.12, '4412','s4412', x)
  CALL WriteFourDigitSection(HTML, 0.04,0.4,0.15, '4415','s4415', x)
  CALL WriteFourDigitSection(HTML, 0.04,0.4,0.18, '4418','s4418', x)
  CALL WriteFourDigitSection(HTML, 0.04,0.4,0.21, '4421','s4421', x)
  CALL WriteFourDigitSection(HTML, 0.04,0.4,0.24, '4424','s4424', x)

  CALL WriteFiveDigitSection(HTML, 0.3,0.15,0.12, '23012','s23012', x)
  CALL WriteFiveDigitSection(HTML, 0.3,0.15,0.15, '23015','s23015', x)
  CALL WriteFiveDigitSection(HTML, 0.3,0.15,0.18, '23018','s23018', x)
  CALL WriteFiveDigitSection(HTML, 0.3,0.15,0.21, '23021','s23021', x)
  CALL WriteFiveDigitSection(HTML, 0.3,0.15,0.24, '23024','s23024', x)

  WRITE(HTML,*) '<p>Go back to the <a href="avd.html">AVD</a> page.</p>'

  CALL WriteFooter(HTML)
  CALL WriteCrumb(HTML, '4- and 5-Digit Sections')
  CALL WriteBanner(HTML)
  WRITE(HTML,*) '</body></html>'
  CLOSE(UNIT=HTML)
  RETURN
END Subroutine Sections45   ! ----------------------------------------------------

!+
SUBROUTINE Sections6(x)
! ------------------------------------------------------------------------------
! PURPOSE - Create a HTML file with the tables of sections - Appendix III
USE AbbottVonDoenhoffSubs
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coor for calculation
!-------------------------------------------------------------------------------
  OPEN(UNIT=HTML,FILE='sections6.html',STATUS='REPLACE',ACTION='WRITE')
  CALL WriteHead(HTML, 'sections6.html', 'Appendix III - 6-Series Sections')

  WRITE(HTML,*) '<body class="sections6">'
  WRITE(HTML,*) '<p id="TopSections"></p>'
  CALL WriteCrumb(HTML, '6-Series Sections')
  CALL WriteBanner(HTML)
  WRITE(HTML,*) '<div id="header"><h1>Appendix III - 6-Series Sections</h1></div>'

  WRITE(HTML,*) '<p><a href="avd.html">Go back to the AVD main page</a><br />'
  WRITE(HTML,*) PAUSENOTE
  WRITE(HTML,*) '</p>'


  WRITE(HTML,*) '<p>63-series cambered sections<br />'
!  CALL WriteTextTable(HTML, (/ &
!    "63-206", "63-209", "63-210", "63-212", "63-412", &
!    "63-215", "63-415", "63-615", "63-218", "63-418", &
!    "63-618", "63-221", "63-421" /), 4)
  WRITE(HTML,*) '<a href="#s63206">63-206</a>'
  WRITE(HTML,*) '<a href="#s63209">63-209</a>'
  WRITE(HTML,*) '<a href="#s63210">63-210</a>'
  WRITE(HTML,*) '<a href="#s63212">63-212</a>'
  WRITE(HTML,*) '<a href="#s63412">63-412</a><br />'
  WRITE(HTML,*) '<a href="#s63215">63-215</a>'
  WRITE(HTML,*) '<a href="#s63415">63-415</a>'
  WRITE(HTML,*) '<a href="#s63615">63-615</a>'
  WRITE(HTML,*) '<a href="#s63218">63-218</a>'
  WRITE(HTML,*) '<a href="#s63418">63-418</a><br />'
  WRITE(HTML,*) '<a href="#s63618">63-618</a>'
  WRITE(HTML,*) '<a href="#s63221">63-221</a>'
  WRITE(HTML,*) '<a href="#s63421">63-421</a></p>'

  WRITE(HTML,*) '<p>64-series cambered sections<br />'
!  CALL WriteTextTable(HTML, (/ &
!    "64-108", "64-110", "64-206", "64-208", "64-209", &
!    "64-210", "64-112", "64-212", "64-412", "64-215", &
!    "64-415", "64-218", "64-418", "64-618", "64-421" /), 4)
  WRITE(HTML,*) '<a href="#s64108">64-108</a>'
  WRITE(HTML,*) '<a href="#s64110">64-110</a>'
  WRITE(HTML,*) '<a href="#s64206">64-206</a>'
  WRITE(HTML,*) '<a href="#s64208">64-208</a>'
  WRITE(HTML,*) '<a href="#s64209">64-209</a><br />'
  WRITE(HTML,*) '<a href="#s64210">64-210</a>'
  WRITE(HTML,*) '<a href="#s64112">64-112</a>'
  WRITE(HTML,*) '<a href="#s64212">64-212</a>'
  WRITE(HTML,*) '<a href="#s64412">64-412</a>'
  WRITE(HTML,*) '<a href="#s64215">64-215</a><br />'
  WRITE(HTML,*) '<a href="#s64415">64-415</a>'
  WRITE(HTML,*) '<a href="#s64218">64-218</a>'
  WRITE(HTML,*) '<a href="#s64418">64-418</a>'
  WRITE(HTML,*) '<a href="#s64618">64-618</a>'
  WRITE(HTML,*) '<a href="#s64421">64-421</a></p>'
    

  WRITE(HTML,*) '<p>65-series cambered sections<br />'
!  CALL WriteTextTable(HTML,(/ &
!    "65-206     ", "65-209     ", "65-210     ", "65-410     ", "65-212     ", &
!    "65-212a=0.5", "65-412     ", "65-215     ", "65-415     ", "65-415a=0.5", &
!    "65-218     ", "65-418     ", "65-418a=0.5", "65-618     ", "65-618a=0.5", &
!    "65-221     ", "65-421     ", "65-421a=0.5" /), 4)
  WRITE(HTML,*) '<a href="#s65206">65-206</a>'
  WRITE(HTML,*) '<a href="#s65209">65-209</a>'
  WRITE(HTML,*) '<a href="#s65210">65-210</a>'
  WRITE(HTML,*) '<a href="#s65410">65-410</a>'
  WRITE(HTML,*) '<a href="#s65212">65-212</a><br />'
  WRITE(HTML,*) '<a href="#s65212a5">65-212 a=0.5</a>'
  WRITE(HTML,*) '<a href="#s65412">65-412</a>'
  WRITE(HTML,*) '<a href="#s65215">65-215</a>'
  WRITE(HTML,*) '<a href="#s65415">65-415</a>'
  WRITE(HTML,*) '<a href="#s65415a5">65-415 a=0.5</a><br />'
  WRITE(HTML,*) '<a href="#s65218">65-218</a>'
  WRITE(HTML,*) '<a href="#s65418">65-418</a>'
  WRITE(HTML,*) '<a href="#s65418a5">65-418 a=0.5</a>'
  WRITE(HTML,*) '<a href="#s65618">65-618</a>'
  WRITE(HTML,*) '<a href="#s65618a5">65-618 a=0.5</a><br />'
  WRITE(HTML,*) '<a href="#s65221">65-221</a>'
  WRITE(HTML,*) '<a href="#s65421">65-421</a>'
  WRITE(HTML,*) '<a href="#s65421a5">65-421 a=0.5</a></p>'

  WRITE(HTML,*) "<p>66-series cambered sections<br />"
!  CALL WriteTextTable(HTML,(/ &
!    "66-206", "66-209", "66-210", "66-212", "66-215", &
!    "66-415", "66-218", "66-418", "66-221" /), 4)
  WRITE(HTML,*) '<a href="#s66206">66-206</a>'
  WRITE(HTML,*) '<a href="#s66209">66-209</a>'
  WRITE(HTML,*) '<a href="#s66210">66-210</a>'
  WRITE(HTML,*) '<a href="#s66212">66-212</a><br />'
  WRITE(HTML,*) '<a href="#s66215">66-215</a>'
  WRITE(HTML,*) '<a href="#s66415">66-415</a>'
  WRITE(HTML,*) '<a href="#s66218">66-218</a>'
  WRITE(HTML,*) '<a href="#s66418">66-418</a>'
  WRITE(HTML,*) '<a href="#s66221">66-221</a></p>'

  WRITE(HTML,*) '<p>67-series cambered sections<br />'
!  CALL WriteTextTable(HTML,(/ "67-212", "67-215" /), 2 )
  WRITE(HTML,*) '<a href="#s67212">67-212</a>'
  WRITE(HTML,*) '<a href="#s67215">67-215</a></p>'

  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.2, 0.06, '63-206','s63206', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.2, 0.09, '63-209','s63209', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.2, 0.10, '63-210','s63210', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.2, 0.12, '63-212','s63212', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.4, 0.12, '63-412','s63412', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.2, 0.15, '63-215','s63215', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.4, 0.15, '63-415','s63415', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.6, 0.15, '63-615','s63615', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.2, 0.18, '63-218','s63218', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.4, 0.18, '63-418','s63418', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.6, 0.18, '63-618','s63618', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.2, 0.21, '63-221','s63221', x)
  CALL WriteSixSeriesSection(HTML, 1, 1.0, 0.4, 0.21, '63-421','s63421', x)

  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.1, 0.08, '64-108','s64108', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.1, 0.10, '64-110','s64110', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.2, 0.06, '64-206','s64206', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.2, 0.08, '64-208','s64208', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.2, 0.09, '64-209','s64209', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.2, 0.10, '64-210','s64210', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.1, 0.12, '64-112','s64112', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.2, 0.12, '64-212','s64212', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.4, 0.12, '64-412','s64412', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.2, 0.15, '64-215','s64215', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.4, 0.15, '64-415','s64415', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.2, 0.18, '64-218','s64218', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.4, 0.18, '64-418','s64418', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.6, 0.18, '64-618','s64618', x)
  CALL WriteSixSeriesSection(HTML, 2, 1.0, 0.4, 0.21, '64-421','s64421', x)
  
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.1, 0.08, '65-108','s65108', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.1, 0.10, '65-110','s65110', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.2, 0.06, '65-206','s65206', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.2, 0.08, '65-208','s65208', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.2, 0.09, '65-209','s65209', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.2, 0.10, '65-210','s65210', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.4, 0.10, '65-410','s65410', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.1, 0.12, '65-112','s65112', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.2, 0.12, '65-212','s65212', x)
  CALL WriteSixSeriesSection(HTML, 3, 0.5, 0.2, 0.12, '65-212a=0.5','s65212a5', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.4, 0.12, '65-412','s65412', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.2, 0.15, '65-215','s65215', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.4, 0.15, '65-415','s65415', x)
  CALL WriteSixSeriesSection(HTML, 3, 0.5, 0.4, 0.15, '65-415a=0.5','s65415a5', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.2, 0.18, '65-218','s65218', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.4, 0.18, '65-418','s65418', x)
  CALL WriteSixSeriesSection(HTML, 3, 0.5, 0.4, 0.18, '65-418a=0.5','s65418a5', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.6, 0.18, '65-618','s65618', x)
  CALL WriteSixSeriesSection(HTML, 3, 0.5, 0.6, 0.18, '65-618a=0.5','s65618a5', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.2, 0.21, '65-221','s65221', x)
  CALL WriteSixSeriesSection(HTML, 3, 1.0, 0.4, 0.21, '65-421','s65421', x)
  CALL WriteSixSeriesSection(HTML, 3, 0.5, 0.4, 0.21, '65-421a=0.5','s65421a5', x)
   
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.1, 0.08, '66-108','s66108', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.1, 0.10, '66-110','s66110', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.2, 0.06, '66-206','s66206', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.2, 0.08, '66-208','s66208', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.2, 0.09, '66-209','s66209', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.2, 0.10, '66-210','s66210', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.1, 0.12, '66-112','s66112', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.2, 0.12, '66-212','s66212', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.4, 0.12, '66-412','s66412', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.2, 0.15, '66-215','s66215', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.4, 0.15, '66-415','s66415', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.2, 0.18, '66-218','s66218', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.4, 0.18, '66-418','s66418', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.6, 0.18, '66-618','s66618', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.2, 0.21, '66-221','s66221', x)
  CALL WriteSixSeriesSection(HTML, 4, 1.0, 0.4, 0.21, '66-421','s66421', x)
  
  CALL WriteSixSeriesSection(HTML, 5, 1.0, 0.2, 0.12, '67-212','s67212', x)
  CALL WriteSixSeriesSection(HTML, 5, 1.0, 0.2, 0.15, '67-215','s67215', x)
 
  WRITE(HTML,*) '<p>Go back to the <a href="avd.html">AVD</a> page.</p>'

  CALL WriteFooter(HTML)
  CALL WriteCrumb(HTML, '6-Series Sections')
  CALL WriteBanner(HTML)
  WRITE(HTML,*) '</body></html>'
  CLOSE(UNIT=HTML)
  RETURN
END Subroutine Sections6   ! ----------------------------------------------------

!+
SUBROUTINE Sections6A(x)
! ------------------------------------------------------------------------------
! PURPOSE - Create a HTML file with the tables of sections - Appendix III
USE AbbottVonDoenhoffSubs
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coor for calculation
!-------------------------------------------------------------------------------
  OPEN(UNIT=HTML,FILE='sections6a.html',STATUS='REPLACE',ACTION='WRITE')
  CALL WriteHead(HTML, 'sections6a.html', 'Appendix III - 6A-Series Sections')

  WRITE(HTML,*) '<body class="sections6A">'
  WRITE(HTML,*) '<p id="TopSections"></p>'
  CALL WriteCrumb(HTML, '6A-Series Sections')
  CALL WriteBanner(HTML)
  WRITE(HTML,*) '<div id="header"><h1>Appendix III - 6A-Series Sections</h1></div>'

  WRITE(HTML,*) '<p><a href="avd.html">Go back to the AVD main page</a><br />'
  WRITE(HTML,*) 'Only a few of these 6A sections were published in the book.'
  WRITE(HTML,*) 'I have added those that I know were used on a'
  WRITE(HTML,*) 'production airplane.<br />'
  WRITE(HTML,*) PAUSENOTE
  WRITE(HTML,*) '</p>'
  
  WRITE(HTML,*) '<p>63A-series cambered airfoils<br />'
!  CALL WriteTextTable(HTML, (/"63A112", "63A210", "63A412", &
!    "63A415", "63A418", "63A421", "63A615" /), 4)
  WRITE(HTML,*) '<a href="#s63A112">63A112</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s63A210">63A210</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s63A412">63A412</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s63A415">63A415</a> &nbsp;<br />'
  WRITE(HTML,*) '<a href="#s63A418">63A418</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s63A421">63A421</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s63A615">63A615</a> &nbsp;</p>'
 
  WRITE(HTML,*) '<p>64A-series cambered airfoils<br />'
!  CALL WriteTextTable(HTML, (/ "64A204", "64A106",              &
!    "64A210", "64A410", "64A212", "64A412", "64A612", "64A114",  &
!    "64A215", "64A415", "64A218",  "64A221" /), 4)
  WRITE(HTML,*) '<a href="#s64A204">64A204</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s64A106">64A106</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s64A210">64A210</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s64A410">64A410</a> &nbsp;<br />'
  WRITE(HTML,*) '<a href="#s64A212">64A212</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s64A412">64A412</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s64A612">64A612</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s64A114">64A114</a> &nbsp;<br />'
  WRITE(HTML,*) '<a href="#s64A215">64A215</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s64A415">64A415</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s64A218">64A218</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s64A221">64A221</a> &nbsp;</p>'

  WRITE(HTML,*) '<p>65A-series cambered airfoils<br />'
!  CALL WriteTextTable(HTML, (/ "65A112", "65A210", "65A212", "65A412", &
!    "65A215", "65A415", "65A615", "65A418", "65A421" /), 4 )
  WRITE(HTML,*) '<a href="#s65A112">65A112</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s65A210">65A210</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s65A212">65A212</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s65A412">65A412</a> &nbsp;<br />'
  WRITE(HTML,*) '<a href="#s65A215">65A215</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s65A415">65A415</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s65A615">65A615</a> &nbsp;<br />'
  WRITE(HTML,*) '<a href="#s65A418">65A418</a> &nbsp;'
  WRITE(HTML,*) '<a href="#s65A421">65A421</a> &nbsp;</p>'
  
  CALL WriteSixSeriesSection(HTML, 6, 0.8, 0.2, 0.10, '63A210','s63A210', x)
  CALL WriteSixSeriesSection(HTML, 6, 0.8, 0.1, 0.12, '63A112','s63A112', x)
  CALL WriteSixSeriesSection(HTML, 6, 0.8, 0.4, 0.12, '63A412','s63A412', x)
  CALL WriteSixSeriesSection(HTML, 6, 0.8, 0.4, 0.15, '63A415','s63A415', x)
  CALL WriteSixSeriesSection(HTML, 6, 0.8, 0.6, 0.15, '63A615','s63A615', x)
  CALL WriteSixSeriesSection(HTML, 6, 0.8, 0.4, 0.18, '63A418','s63A418', x)
  CALL WriteSixSeriesSection(HTML, 6, 0.8, 0.4, 0.21, '63A421','s63A421', x)

  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.2, 0.04, '64A204','s64A204', x)
  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.1, 0.06, '64A106','s64A106', x)  
!  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.2, 0.04, '64A204','s63A204', x)  
!  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.2, 0.04, '64A206','s63A206', x)

  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.2, 0.10, '64A210','s64A210', x)
  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.4, 0.10, '64A410','s64A410', x)
  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.2, 0.12, '64A212','s64A212', x)
  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.4, 0.12, '64A412','s64A412', x)
  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.6, 0.12, '64A612','s64A612', x)
  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.1, 0.14, '64A114','s64A114', x)
  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.2, 0.15, '64A215','s64A215', x)
  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.4, 0.15, '64A415','s64A415', x)
  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.2, 0.18, '64A218','s64A218', x)
  CALL WriteSixSeriesSection(HTML, 7, 0.8, 0.2, 0.21, '64A221','s64A221', x)
  
  CALL WriteSixSeriesSection(HTML, 8, 0.8, 0.2, 0.10, '65A210','s65A210', x)
  CALL WriteSixSeriesSection(HTML, 8, 0.8, 0.2, 0.12, '65A112','s65A112', x)
  CALL WriteSixSeriesSection(HTML, 8, 0.8, 0.2, 0.12, '65A212','s65A212', x)
  CALL WriteSixSeriesSection(HTML, 8, 0.8, 0.4, 0.12, '65A412','s65A412', x)
  CALL WriteSixSeriesSection(HTML, 8, 0.8, 0.2, 0.15, '65A215','s65A215', x)
  CALL WriteSixSeriesSection(HTML, 8, 0.8, 0.4, 0.15, '65A415','s65A415', x)
  CALL WriteSixSeriesSection(HTML, 8, 0.8, 0.6, 0.15, '65A615','s65A615', x)
  CALL WriteSixSeriesSection(HTML, 8, 0.8, 0.4, 0.18, '65A418','s65A418', x)
  CALL WriteSixSeriesSection(HTML, 8, 0.8, 0.4, 0.21, '65A421','s65A421', x)
  
  WRITE(HTML,*) '<p>Go back to the <a href="avd.html">AVD</a> page.</p>'

  CALL WriteFooter(HTML)
  CALL WriteCrumb(HTML, '6A-Series Sections')
  CALL WriteBanner(HTML)
  WRITE(HTML,*) '</body></html>'
  CLOSE(UNIT=HTML)
  RETURN
END Subroutine Sections6A   ! --------------------------------------------------


END Module MainProcedures   ! ==================================================


!+
PROGRAM AbbottVonDoenhoff                                 ! pdas/naca456/avd.f90
! ------------------------------------------------------------------------------
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

USE ISO_FORTRAN_ENV
USE NACAauxilary,ONLY: LoadX
USE MainProcedures,ONLY: MainPage, Profiles, MeanLines, Sections45, Sections6, Sections6A
!USE AbbottVonDoenhoffSubs, ONLY: address1  !,ADDRESS2,ADDRESS3,WriteBanner
IMPLICIT NONE
  
  CHARACTER(LEN=*),PARAMETER:: AUTHOR = &
    "Ralph L. Carmichael, Public Domain Aeronautical Software"
  CHARACTER(LEN=15):: dateTimeStr
  CHARACTER(LEN=*),PARAMETER:: GREETING="avd - format airfoil data"
  CHARACTER(LEN=*),PARAMETER:: MODIFIER=" "  ! your name if you change it

  INTEGER,PARAMETER:: denCode=1   ! tells LoadX to use coarse spacing
  INTEGER:: nx                    ! length of x  (from LoadX)
  REAL,DIMENSION(100):: x         ! x-coor for calculation (from LoadX)
                                  ! 100 is more than enough
                                  
  CHARACTER(LEN=*),PARAMETER:: VERSION=" 1.8 (2021 Apr 26)"
!-------------------------------------------------------------------------------
  WRITE(*,*) GREETING
  WRITE(*,*) AUTHOR
  IF (MODIFIER /= ' ') WRITE(*,*) 'Modified by '//MODIFIER
  WRITE(*,*) VERSION
  WRITE(*,'(3A/A)') 'This file was compiled by ', compiler_version(),   &
                    ' using the options ',        compiler_options()


  CALL MainPage()
  CALL LoadX(denCode, nx,x)

  CALL Profiles(x(1:nx))
  CALL MeanLines(x(1:nx))
  CALL Sections45(x(1:nx))
  CALL Sections6(x(1:nx))
  CALL Sections6A(x(1:nx))

  WRITE(*,*) 'Normal termination of program avd.'
  WRITE(*,*) 'Files avd.html, profiles.html, meanlines.html, sections45.html,'
  WRITE(*,*) '   sections6.html, and sections6a.html added to your directory'
  STOP
END Program AbbottVonDoenhoff   ! ==============================================
