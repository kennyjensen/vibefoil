!+
! PROGRAM AbbottVonDoenhoff                         ! pdas/naca456/avdpdf.f90
! ---------------------------------------------------------------------------
! PURPOSE - Make tables of the data in Abbott and von Doenhoff book
!   entitled "Theory of Airfoil Sections".
!   The output will be a LaTeX file
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 2021Apr29  0.1   RLC   Copied avd.f90 as a guide
! 2021May01  0.2   RLC   Added LaTeX coding
! 2021May05  0.3   RLC   Added index
! 2021May21  0.4   RLC   Removed HTML id variable
!
! NOTES- If you should recompile this source code, you should use the
!   appropriate compiler flag making the default real kind as 64-bit.
!   With gfortran this is -fdefault-real-8
! NOTE - If you want to modify and recompile this program, you will need the
! three files below as well as avdpdf.f90. These are also used by naca456.f90.
!-------------------------------------------------------------------------------
INCLUDE 'splprocs.f90'
INCLUDE 'epspsi.f90'
INCLUDE 'nacax.f90'

!+
MODULE AbbottVonDoenhoffSubsPDF
! ------------------------------------------------------------------------------
! PURPOSE - Subroutines called by Program PDFVonDoenhoff.
!   Make tables of the data in Abbott and von Doenhoff.
!   This module is not used by the naca456 program.  
!-------------------------------------------------------------------------------
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

  PRIVATE:: WriteColumns

CONTAINS

!+
SUBROUTINE WriteFourDigitProfile(efu,toc,name,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a four-digit profile
USE NACAauxilary,ONLY: Thickness4,LeadingEdgeRadius4
! CALLed from Profiles
! CALLs module procedures Thickness4 and WriteOneProfile

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: toc      ! thickness / chord   (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   ! profile name
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: rle ! radius of leading edge, fraction of chord
  REAL,DIMENSION(SIZE(x)):: yt,ytp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
  CALL Thickness4(toc, x,yt,ytp)   ! compute yt and ytp at each x-value
  rle=LeadingEdgeRadius4(toc)
  CALL WriteOneProfile(efu,name,x,yt,ytp,rle)
  RETURN
END Subroutine WriteFourDigitProfile   ! ---------------------------------------

!+
SUBROUTINE WriteModifiedFourDigitProfile(efu,leIndex,xmaxt,toc,name,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a modified four-digit profile
USE NACAauxilary,ONLY: Thickness4M,LeadingEdgeRadius4M

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: leIndex  ! leading edge index
                             ! coded to allow leIndex to be a REAL
  REAL,INTENT(IN):: xmaxt    ! x-coor of maximum thickness
  REAL,INTENT(IN):: toc      ! thickness / chord   (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   ! profile name
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: rle
  REAL,DIMENSION(SIZE(x)):: yt,ytp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
  CALL Thickness4M(toc,leIndex,xmaxt,x,yt,ytp)  ! compute yt and ytp at each x-value:
  rle=LeadingEdgeRadius4M(toc,leIndex)
  CALL WriteOneProfile(efu,name,x,yt,ytp,rle)
  RETURN
END Subroutine WriteModifiedFourDigitProfile   ! -------------------------------

!+
SUBROUTINE WriteSixSeriesProfile(efu,family,toc,name, x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a six-series profile.
USE NACAauxilary,ONLY: Thickness6,LeadingEdgeRadius6

  INTEGER,INTENT(IN):: efu   ! external file unit
  INTEGER,INTENT(IN):: family  ! =1 63; =2 64; =3 65; =4 66; =5 67
                               ! =6 63A; =7 64A; =8 65A
  REAL,INTENT(IN):: toc      ! thickness / chord (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   ! profile name
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: rle   ! leading edge radius, fraction of chord
  REAL,DIMENSION(SIZE(x)):: yt,ytp   ! tables of y and dy/dx
!-------------------------------------------------------------------------------
  CALL Thickness6(family,toc, x,yt,ytp)   ! compute yt and ytp at each x-value:
  rle=LeadingEdgeRadius6(family,toc)
  CALL WriteOneProfile(efu,name,x,yt,ytp,rle)
  RETURN
END Subroutine WriteSixSeriesProfile   ! ---------------------------------------

!+
SUBROUTINE WriteTwoDigitMeanLine(efu,cmax,xmaxc,name,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a four-digit meanline.
USE NACAauxilary,ONLY: MeanLine2

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cmax   ! max. camber (fraction of chord)
  REAL,INTENT(IN):: xmaxc  ! x-location of max. camber (fraction of chord) 
  CHARACTER(LEN=*),INTENT(IN):: name   ! mean line name
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
  CALL MeanLine2(cmax,xmaxc, x,ymean,ymeanp)   ! compute ymean and ymeanp at each x-value:
  CALL WriteOneMeanLine(efu,name,x,ymean,ymeanp)
  RETURN
END Subroutine WriteTwoDigitMeanLine   ! --------------------------------------

!+
SUBROUTINE WriteThreeDigitMeanLine(efu,cl,xmaxc,name,x)
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
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
  CALL MeanLine3(cl,xmaxc,x,ymean,ymeanp)   ! compute ymean and ymeanp at each x-value:
  CALL WriteOneMeanLine(efu,name,x, ymean,ymeanp)
  RETURN
END Subroutine WriteThreeDigitMeanLine   ! -------------------------------------

!+
SUBROUTINE WriteThreeDigitReflexMeanLine(efu,cl,xmaxc,name,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a 3-digit reflex meanline
! REF - p.8 of NASA Technical Memorandum 4741 and NACA Report 537.
USE NACAauxilary,ONLY: MeanLine3Reflex

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cl    ! design lift coefficient (same as non-reflex)
  REAL,INTENT(IN):: xmaxc !  x-location of max. camber (fraction of chord) 
  CHARACTER(LEN=*),INTENT(IN):: name   ! mean line name
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
  CALL MeanLine3Reflex(cl,xmaxc,x,ymean,ymeanp)   ! compute ymean and ymeanp at each x-value:
  CALL WriteOneMeanLine(efu,name,x,ymean,ymeanp)
  RETURN
END Subroutine WriteThreeDigitReflexMeanLine   ! -------------------------------

!+
SUBROUTINE WriteSixSeriesMeanLine(efu,a,cl,name,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a six-series meanline
USE NACAauxilary,ONLY: MeanLine6

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: a    ! x-coor of aft-position of pressure rooftop
  REAL,INTENT(IN):: cl   ! design lift coefficient
  CHARACTER(LEN=*),INTENT(IN):: name   ! mean line name
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
!-------------------------------------------------------------------------------
  CALL MeanLine6(a,cl,x,ymean,ymeanp)   ! compute ymean and ymeanp at each x-value:
  CALL WriteOneMeanLine(efu, name, x,ymean,ymeanp)
  RETURN
END Subroutine WriteSixSeriesMeanLine   ! --------------------------------------

!+
SUBROUTINE WriteSixSeriesModifiedMeanLine(efu,cl,name,x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a six-series modified meanline
! NOTE - a is not input; a is always 0.8 
USE NACAauxilary,ONLY: MeanLine6M

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cl
  CHARACTER(LEN=*),INTENT(IN):: name   ! mean line name
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp
!-------------------------------------------------------------------------------
  CALL MeanLine6M(cl,x,ymean,ymeanp)   ! compute ymean and ymeanp at each x-value:
  CALL WriteOneMeanLine(efu,name,x,ymean,ymeanp)
  RETURN
END Subroutine WriteSixSeriesModifiedMeanLine   ! ------------------------------

!+
SUBROUTINE WriteFourDigitSection(efu,cmax,xmaxc,toc,name, x)
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
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: leSlope
  REAL:: rle   ! leading edge radius, fraction of chord
  REAL,DIMENSION(SIZE(x)):: yt   ! table of half-thickness
  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
  REAL,DIMENSION(SIZE(x)):: xupper,yupper, xlower,ylower
  REAL,DIMENSION(SIZE(x)):: yup,ylo

!-------------------------------------------------------------------------------
  CALL Thickness4(toc, x,yt)   ! compute yt at each x-value:
  rle=LeadingEdgeRadius4(toc)
  leSlope=2.0*cmax/xmaxc
  CALL MeanLine2(cmax,xmaxc, x,ymean,ymeanp)   ! compute ymean and ymeanp at each x-value:
  CALL CombineThicknessAndCamber(x,yt,ymean,ymeanp, &
    xupper,yupper, xlower,ylower)
  CALL InterpolateCombinedAirfoil(x,yt,ymean,ymeanp, yup,ylo)
  CALL WriteOneSection(efu,name,  x,yup, x,ylo, &
    xupper,yupper, xlower,ylower, rle,leSlope)
  RETURN
END Subroutine WriteFourDigitSection   ! ---------------------------------------

!+
SUBROUTINE WriteModifiedFourDigitSection(efu,cmax,xmaxc, &
  leIndex,xmaxt,toc,name, x)
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
  CALL WriteOneSection(efu,name,  x,yup, x,ylo, &
    xupper,yupper, xlower,ylower, rle,leSlope)
  RETURN
END Subroutine WriteModifiedFourDigitSection   ! -------------------------------

!+
SUBROUTINE WriteFiveDigitSection(efu,cl,xmaxc,toc,name, x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a five-digit airfoil.
USE NACAauxilary,ONLY: Thickness4,LeadingEdgeRadius4,MeanLine3,GetRk1
USE NACAauxilary,ONLY: CombineThicknessAndCamber,InterpolateCombinedAirfoil

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cl  ! design lift coefficient
  REAL,INTENT(IN):: xmaxc ! x-coor of max camber (fraction of chord)
  REAL,INTENT(IN):: toc  ! thickness / chord (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   ! section name
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coordinates for calculation

  REAL:: k1
  REAL:: leSlope
  REAL:: r
  REAL:: rle   ! leading edge radius, fraction of chord
  REAL,DIMENSION(SIZE(x)):: yt
  REAL,DIMENSION(SIZE(x)):: ymean,ymeanp   ! yables of y and dy/dx
  REAL,DIMENSION(SIZE(x)):: xupper,yupper, xlower,ylower
  REAL,DIMENSION(SIZE(x)):: yup,ylo
!-------------------------------------------------------------------------------
  CALL Thickness4(toc, x,yt)
  rle=LeadingEdgeRadius4(toc)
  CALL GetRk1(xmaxc,r,k1)
  leSlope=r*r*(3.0-r)*k1*cl/1.8
  CALL MeanLine3(cl,xmaxc, x,ymean,ymeanp)
  CALL CombineThicknessAndCamber(x,yt,ymean,ymeanp, &
    xupper,yupper, xlower,ylower)
  CALL InterpolateCombinedAirfoil(x,yt,ymean,ymeanp, yup,ylo)
  CALL WriteOneSection(efu, name,  x,yup, x,ylo, &
    xupper,yupper, xlower,ylower, rle,leSlope)
  RETURN
END Subroutine WriteFiveDigitSection   ! ---------------------------------------

!+
SUBROUTINE WriteFiveDigitReflexSection(efu,cl,xmaxc,toc,name, x)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate and print the coordinates of a 5-digit reflex airfoil.
USE NACAauxilary,ONLY: Thickness4,LeadingEdgeRadius4,MeanLine3Reflex,GetRk1k2
USE NACAauxilary,ONLY: CombineThicknessAndCamber,InterpolateCombinedAirfoil

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN):: cl   ! design lift coefficient
  REAL,INTENT(IN):: xmaxc   ! x-coor of max camber (fraction of chord)
  REAL,INTENT(IN):: toc   ! thickness / chord (fraction, not percent)
  CHARACTER(LEN=*),INTENT(IN):: name   !  section name
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
  CALL WriteOneSection(efu, name,  x,yup, x,ylo, &
    xupper,yupper, xlower,ylower, rle,leSlope)
  RETURN
END Subroutine WriteFiveDigitReflexSection   ! ---------------------------------

!+
SUBROUTINE WriteSixSeriesSection(efu,family,a,cl,toc,name, x)
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
  CALL WriteOneSection(efu, name,  x,yup, x,ylo, &
    xupper,yupper, xlower,ylower, rle,leSlope)
  RETURN
END Subroutine WriteSixSeriesSection   ! ---------------------------------------

!+
SUBROUTINE WriteOneMeanLine(efu,name,x,y,yp)
! ------------------------------------------------------------------------------
! PURPOSE - Write the description of one mean line to the TeX file.
  INTEGER,INTENT(IN):: efu   ! external file unit
  CHARACTER(LEN=*),INTENT(IN):: name   ! mean line name     
  REAL,INTENT(IN),DIMENSION(:):: x,y,yp   ! tables of x,y,dy/dz
  CHARACTER(LEN=80):: temp

!-------------------------------------------------------------------------------
  WRITE(efu,*) "\phantomsection\label{ml" // TRIM(name) // "}"
  WRITE(efu,*) "\begin{Large}"
  WRITE(efu,*) "NACA Mean Line " // TRIM(name)
  WRITE(efu,*) "\end{Large}"
  WRITE(efu,*) " "  
  WRITE(efu,*) '\vspace{8mm}'
  
  WRITE(efu,*) '\begin{tabular}{|c|c|c|}  \hline' 
  WRITE(efu,*)  'x & y & dy/dx \\'
  WRITE(efu,*) '\hline'
  CALL WriteColumns(efu, 100.0*x,100.0*y,yp)
 
  WRITE(efu,*) '\hline'
  WRITE(efu,*) '\end{tabular}'
  WRITE(efu,*) '\vspace{8mm}'

  WRITE(efu,'(/A/)') '(Stations and ordinates given in per cent of airfoil chord)'

  WRITE(efu,*) '\newpage'
  RETURN
END Subroutine WriteOneMeanLine   ! --------------------------------------------

!+
SUBROUTINE WriteOneProfile(efu,name, x,y,yp, rle)
! ------------------------------------------------------------------------------
! PURPOSE - Write the description of one thickness to the HTML file.
  INTEGER,INTENT(IN):: efu   ! external file unit
  CHARACTER(LEN=*),INTENT(IN):: name   ! profile name
  REAL,INTENT(IN),DIMENSION(:):: x,y,yp   ! tables of x,y,dy/dz
  REAL,INTENT(IN):: rle   ! leading edge radius, fraction of chord

  CHARACTER(LEN=*),PARAMETER:: FMT = '(//A,F7.4,A)'
!-------------------------------------------------------------------------------
  WRITE(efu,*) "\phantomsection \label{p" // TRIM(name) // "}"
  WRITE(efu,*) "\begin{Large}"
  WRITE(efu,*) "NACA Profile " // TRIM(name) 
  WRITE(efu,*) "\end{Large}"
  WRITE(efu,*) " "
  WRITE(efu,*) "\vspace{8mm}"
  
  WRITE(efu,*) "\begin{tabular}{|c|c|c|} \hline "
  WRITE(efu,*) " x  &  y  &  dy/dx \\"
  WRITE(efu,*) "\hline"

  CALL WriteColumns(efu, 100.0*x,100.0*y,yp)
  WRITE(efu,*) "\hline"
  WRITE(efu,*) "\end{tabular}"
  WRITE(efu,*) "\vspace{8mm}"
  WRITE(efu,'(//A)') "Stations and ordinates given in per cent of airfoil chord "
  WRITE(efu,FMT) "L.E. Radius= ", 100.0*rle, " percent chord"

  WRITE(efu,*) "\newpage"
  RETURN
END Subroutine WriteOneProfile   ! ---------------------------------------------

!+
SUBROUTINE WriteOneSection(efu,name, xup,yup,xlo,ylo, &
  xupper,yupper, xlower,ylower, rle,leSlope)
! ------------------------------------------------------------------------------
! PURPOSE - Write the description of one airfoil section to the HTML file.
! NOTE - Writes a table with one row and two columns. Each of the two cells
!   of this table contains a table with four columns and enough rows to
!   display xup, etc. along with column headers.

  INTEGER,INTENT(IN):: efu   ! external file unit
  CHARACTER(LEN=*),INTENT(IN):: name   ! section name
  REAL,INTENT(IN),DIMENSION(:):: xup,yup,xlo,ylo
  REAL,INTENT(IN),DIMENSION(:):: xupper,yupper,xlower,ylower
  REAL,INTENT(IN):: rle   ! leading edge radius, fraction of chord
  REAL,INTENT(IN):: leSlope   ! slope of mean line at leading edge
                              ! for 6-series mean lines, this is the slope
                              ! at 0.5 per cent chord


  CHARACTER(LEN=*),PARAMETER:: FMT = "(//A,F7.4,A)"

!-------------------------------------------------------------------------------
  WRITE(efu,*) "\phantomsection \label{s" // TRIM(name) // "}"  
  WRITE(efu,*) "\begin{Large}"
  WRITE(efu,*) "NACA Section " // TRIM(name)
  WRITE(efu,*) "\end{Large}"
  WRITE(efu,*) " "
  WRITE(efu,*) '\vspace{8mm}'

  WRITE(efu,*) "\begin{tabular}{|r|r|r|r|} \hline "
  WRITE(efu,*) &
    '\multicolumn{2}{|c|}{Upper surface} & \multicolumn{2}{|c|}{Lower surface} \\'
  WRITE(efu,*) '\hline'
  WRITE(efu,*) 'Station & Ordinate & Station & Ordinate \\'
  WRITE(efu,*) "\hline"
  CALL WriteColumns(efu, 100.0*xup,100.0*yup,100.0*xlo,100.0*ylo)
  WRITE(efu,*) "\hline "
  WRITE(efu,*) "\end{tabular}"
  WRITE(efu,*) "\hspace{4mm}"
  
! the second table follows immediately with no paragraph break
  WRITE(efu,*) "\begin{tabular}{|r|r|r|r|} \hline "
  WRITE(efu,*) &
    '\multicolumn{2}{|c|}{Upper surface} & \multicolumn{2}{|c|}{Lower surface} \\'
  WRITE(efu,*) '\hline'
  WRITE(efu,*) 'Station & Ordinate & Station & Ordinate \\'
  WRITE(efu,*) "\hline"
  CALL WriteColumns(efu, 100.0*xupper,100.0*yupper,100.0*xlower,100.0*ylower)
  WRITE(efu,*) "\hline "
  WRITE(efu,*) "\end{tabular}"

  WRITE(efu,*) "\vspace{8mm}"
  WRITE(efu,'(/A/)') 'Stations and ordinates given in per cent of airfoil chord'
  
  WRITE(efu,FMT) "L.E. Radius= ", 100.0*rle, " percent chord"

  WRITE(efu,FMT) " slope of mean line at LE = ",leSlope
  WRITE(efu,*) "\newpage"
  
  RETURN
END Subroutine WriteOneSection   ! ---------------------------------------------

!+
SUBROUTINE WriteColumns(efu,a,b,c,d,e,f)
! ------------------------------------------------------------------------------
! PURPOSE - Write the data in the arrays to the LaTeX file.
!   Print 4 figures after decimal point separated by ampersands

  INTEGER,INTENT(IN):: efu   ! external file unit
  REAL,INTENT(IN),OPTIONAL,DIMENSION(:):: a,b,c,d,e,f

  CHARACTER(LEN=12):: as,bs,cs,ds,es,fs
  CHARACTER(LEN=132):: buffer
  CHARACTER(LEN=*),PARAMETER:: FMT='(F8.4)'
  INTEGER:: k
!-------------------------------------------------------------------------------
  IF (.NOT.Present(a)) RETURN
  DO k=1,SIZE(a)
    WRITE(as,FMT) a(k)
    IF (Present(b)) WRITE(bs,FMT) b(k)
    IF (Present(c)) WRITE(cs,FMT) c(k)
    IF (Present(d)) WRITE(ds,FMT) d(k)
    IF (Present(e)) WRITE(es,FMT) e(k)
    IF (Present(f)) WRITE(fs,FMT) f(k)

    buffer=Trim(AdjustL(as))
    IF (Present(b)) buffer=Trim(buffer) // " & " // Trim(AdjustL(bs))
    IF (Present(c)) buffer=Trim(buffer) // " & " // Trim(AdjustL(cs))
    IF (Present(d)) buffer=Trim(buffer) // " & " // Trim(AdjustL(ds))
    IF (Present(e)) buffer=Trim(buffer) // " & " // Trim(AdjustL(es))
    IF (Present(f)) buffer=Trim(buffer) // " & " // Trim(AdjustL(fs))
    WRITE(efu,'(A)') Trim(buffer) // ' \\'
  END DO

  RETURN
END Subroutine WriteColumns   ! ------------------------------------------------

END Module AbbottVonDoenhoffSubsPDF   ! ========================================

!+
MODULE LaTeXProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Collect 

USE AbbottVonDoenhoffSubsPDF
IMPLICIT NONE

  PUBLIC:: FrontMatter
  PUBLIC:: MakeIndex
  PUBLIC:: Profiles
  PUBLIC:: Meanlines
  PUBLIC:: Sections45
  PUBLIC:: Sections6
  PUBLIC:: Sections6A
  PUBLIC:: BackMatter
!-------------------------------------------------------------------------------


CONTAINS

!+
SUBROUTINE FrontMatter(efu)
! ------------------------------------------------------------------------------
! PURPOSE - Write the title page and introduction to the LaTeX file

  INTEGER,INTENT(IN):: efu   ! external file unit of the TeX file  
!-------------------------------------------------------------------------------
  WRITE(efu,*) '\documentclass[11pt]{book}'
  WRITE(efu,*) '\usepackage{layout}'
  WRITE(efu,*) '\usepackage[hidelinks]{hyperref}'
  WRITE(efu,*) '\setlength{\textwidth}{6.5in}'
  WRITE(efu,*) '\setlength{\marginparwidth}{0pt}'
  WRITE(efu,*) '\setlength{\oddsidemargin}{5pt}'
  WRITE(efu,*) '\setlength{\evensidemargin}{5pt}'
  WRITE(efu,*) '\pagestyle{plain}'

  WRITE(efu,*) '\begin{document}'
  
!!! COVER PAGE *****************************************************************
  WRITE(efu,*) '\begin{center}'
  WRITE(efu,*) '\vspace{8mm}'
  WRITE(efu,*) '\Huge'
  WRITE(efu,*) '\textbf{Recomputed Tables of Coordinates} \\'
  WRITE(efu,*) '\textbf{of NACA Airfoil} \\'
  WRITE(efu,*) '\textbf{Profiles, Meanlines, and Sections}\\'
  WRITE(efu,*) '\vspace{8mm}'
  WRITE(efu,*) '\large'
  WRITE(efu,*) 'Based on NACA Report 824 \\and the book \textit{Theory of Airfoil Sections}\\'
  WRITE(efu,*) 'by Ira H. Abbott, Alfred E. von Doenhoff, and Louis S. Stivers\\'
  WRITE(efu,*) '\vspace{70mm}'
  WRITE(efu,*) '\LARGE'
  WRITE(efu,*) 'by Ralph L. Carmichael\\'    
  WRITE(efu,*) '\vspace{8mm}'  

  WRITE(efu,*) '\Large'
  WRITE(efu,*) 'Public Domain Aeronautical Software'
  WRITE(efu,*) '\\Santa Cruz, California \\'
  WRITE(efu,*) '\end{center}'
  WRITE(efu,*) '\newpage'
  WRITE(efu,*) 'This edition compiled on \today \\'
  WRITE(efu,*) 'The latest version of this document may be downloaded from '
  WRITE(efu,*) ' https://www.pdas.com/avd.pdf \\'
  WRITE(efu,*) 'The latest version of the LaTeX source file may be downloaded from '
  WRITE(efu,*) ' https://www.pdas.com/tex/avdpdf.tex \\'
  WRITE(efu,*) 'The latest version of the Fortran source code that writes the '
  WRITE(efu,*) ' LaTeX source file is included in the archive file '
  WRITE(efu,*) ' https://www.pdas.com/packages/naca456.zip \\'
  WRITE(efu,*) '\newpage'
  
!!!  INTRODUCTION **************************************************************
  WRITE(efu,*) '\begin{center} \Large Introduction \normalsize \end{center}'
  WRITE(efu,*) 'The referenced book has three large appendices'
  WRITE(efu,*) 'containing tables of coordinates of a selection of NACA'
  WRITE(efu,*) 'profiles, mean lines and sections.'
  WRITE(efu,*) 'This page and the linked pages reproduces these data'
  WRITE(efu,*) 'recomputed using the procedures of the PDAS program naca456'
  WRITE(efu,*) 'module.  \\'
  
  WRITE(efu,*) 'There are several differences between this document and the'
  WRITE(efu,*) 'referenced book.'
  WRITE(efu,*)  '\begin{itemize}'
  WRITE(efu,*)  '\item Results are shown to 4 decimal places.'
  WRITE(efu,*)  '\item A few extra sections of popular airfoils have been added.'
  WRITE(efu,*)  '\item Two tables are shown for each section:'
  WRITE(efu,*)  '\begin{itemize}'
  WRITE(efu,*)  '\item The left table presents data interpolated at the same x-location.'
  WRITE(efu,*)  '\item The right table presents data with the thickness added ' 
  WRITE(efu,*)  'perpendicular to the mean line (as in the reference book).'
  WRITE(efu,*)  '\end{itemize}'
  WRITE(efu,*)  '\end{itemize}'
  WRITE(efu,*)  

  WRITE(efu,*) '\newpage'
!!!  WRITE(efu,*) '\layout'
  WRITE(efu,*) '\newpage'
  
  CALL MakeIndex(efu)
  
  RETURN
END Subroutine FrontMatter   ! -------------------------------------------------

!+
SUBROUTINE MakeIndex(efu)
! ------------------------------------------------------------------------------
! PURPOSE - Write the LaTeX statements to make a two-column index.
!   Be sure to return to one-column mode at the end.

  INTEGER,INTENT(IN):: efu   ! external file unit of the LaTeX file
  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(/A)'
!-------------------------------------------------------------------------------
  WRITE(efu,*) ' '
  WRITE(efu,*) '\twocolumn'
  WRITE(efu,FMT1) '\begin{theindex}'
  WRITE(efu,*) 'Profiles'
  WRITE(efu,*) '\item 0006, \hyperref[p0006]{\pageref{p0006}}'
  WRITE(efu,*) '\item 0008, \hyperref[p0008]{\pageref{p0008}}'
  WRITE(efu,*) '\item 0010, \hyperref[p0010]{\pageref{p0010}}'
  WRITE(efu,*) '\item 0012, \hyperref[p0012]{\pageref{p0012}}'
  WRITE(efu,*) '\item 0015, \hyperref[p0015]{\pageref{p0015}}'
  WRITE(efu,*) '\item 0018, \hyperref[p0018]{\pageref{p0018}}'
  WRITE(efu,*) '\item 0021, \hyperref[p0021]{\pageref{p0021}}'
  WRITE(efu,*) '\item 0008-34, \hyperref[p0008-34]{\pageref{p0008-34}}'
  WRITE(efu,*) '\item 0010-34, \hyperref[p0010-34]{\pageref{p0010-34}}'
  WRITE(efu,*) '\item 0010-35, \hyperref[p0010-35]{\pageref{p0010-35}}'
  WRITE(efu,*) '\item 0010-64, \hyperref[p0010-64]{\pageref{p0010-64}}'
  WRITE(efu,*) '\item 0010-65, \hyperref[p0010-65]{\pageref{p0010-65}}'
  WRITE(efu,*) '\item 0010-66, \hyperref[p0010-66]{\pageref{p0010-66}}'
  WRITE(efu,*) '\item 0012-34, \hyperref[p0012-34]{\pageref{p0012-34}}'
  WRITE(efu,*) '\item 0012-64, \hyperref[p0012-64]{\pageref{p0012-64}}'
  WRITE(efu,*) '\item 16-006, \hyperref[p16-006]{\pageref{p16-006}}'  
  WRITE(efu,*) '\item 16-009, \hyperref[p16-009]{\pageref{p16-009}}'  
  WRITE(efu,*) '\item 16-012, \hyperref[p16-012]{\pageref{p16-012}}'  
  WRITE(efu,*) '\item 16-015, \hyperref[p16-015]{\pageref{p16-015}}'  
  WRITE(efu,*) '\item 16-018, \hyperref[p16-018]{\pageref{p16-018}}'  
  WRITE(efu,*) '\item 16-021, \hyperref[p16-021]{\pageref{p16-021}}'  

  WRITE(efu,*) ' '
  WRITE(efu,*) '6- and 6A-series profiles' 
  
  WRITE(efu,*) '\item 63-006, \hyperref[p63-006]{\pageref{p63-006}}'  
  WRITE(efu,*) '\item 63-008, \hyperref[p63-008]{\pageref{p63-008}}'  
  WRITE(efu,*) '\item 63-009, \hyperref[p63-009]{\pageref{p63-009}}'  
  WRITE(efu,*) '\item 63-010, \hyperref[p63-010]{\pageref{p63-010}}'  
  WRITE(efu,*) '\item 63-012, \hyperref[p63-012]{\pageref{p63-012}}'  
  WRITE(efu,*) '\item 63-015, \hyperref[p63-015]{\pageref{p63-015}}'  
  WRITE(efu,*) '\item 63-018, \hyperref[p63-018]{\pageref{p63-018}}'  
  WRITE(efu,*) '\item 63-021, \hyperref[p63-021]{\pageref{p63-021}}'  

  WRITE(efu,*) '\item 63A006, \hyperref[p63A006]{\pageref{p63A006}}'  
  WRITE(efu,*) '\item 63A008, \hyperref[p63A008]{\pageref{p63A008}}'  
  WRITE(efu,*) '\item 63A010, \hyperref[p63A010]{\pageref{p63A010}}'  
  WRITE(efu,*) '\item 63A012, \hyperref[p63A012]{\pageref{p63A012}}'  
  WRITE(efu,*) '\item 63A015, \hyperref[p63A015]{\pageref{p63A015}}'  

  WRITE(efu,*) '\item 64-006, \hyperref[p64-006]{\pageref{p64-006}}'  
  WRITE(efu,*) '\item 64-008, \hyperref[p64-008]{\pageref{p64-008}}'  
  WRITE(efu,*) '\item 64-009, \hyperref[p64-009]{\pageref{p64-009}}'  
  WRITE(efu,*) '\item 64-010, \hyperref[p64-010]{\pageref{p64-010}}'  
  WRITE(efu,*) '\item 64-012, \hyperref[p64-012]{\pageref{p64-012}}'  
  WRITE(efu,*) '\item 64-015, \hyperref[p64-015]{\pageref{p64-015}}'  
  WRITE(efu,*) '\item 64-018, \hyperref[p64-018]{\pageref{p64-018}}'  
  WRITE(efu,*) '\item 64-021, \hyperref[p64-021]{\pageref{p64-021}}'  

  WRITE(efu,*) '\item 64A006, \hyperref[p64A006]{\pageref{p64A006}}'  
  WRITE(efu,*) '\item 64A008, \hyperref[p64A008]{\pageref{p64A008}}'  
  WRITE(efu,*) '\item 64A010, \hyperref[p64A010]{\pageref{p64A010}}'  
  WRITE(efu,*) '\item 64A012, \hyperref[p64A012]{\pageref{p64A012}}'  
  WRITE(efu,*) '\item 64A015, \hyperref[p64A015]{\pageref{p64A015}}'  

  WRITE(efu,*) '\item 65-006, \hyperref[p65-006]{\pageref{p65-006}}'  
  WRITE(efu,*) '\item 65-008, \hyperref[p65-008]{\pageref{p65-008}}'  
  WRITE(efu,*) '\item 65-009, \hyperref[p65-009]{\pageref{p65-009}}'  
  WRITE(efu,*) '\item 65-010, \hyperref[p65-010]{\pageref{p65-010}}'  
  WRITE(efu,*) '\item 65-012, \hyperref[p65-012]{\pageref{p65-012}}'  
  WRITE(efu,*) '\item 65-015, \hyperref[p65-015]{\pageref{p65-015}}'  
  WRITE(efu,*) '\item 65-018, \hyperref[p65-018]{\pageref{p65-018}}'  
  WRITE(efu,*) '\item 65-021, \hyperref[p65-021]{\pageref{p65-021}}'  

  WRITE(efu,*) '\item 65A006, \hyperref[p65A006]{\pageref{p65A006}}'  
  WRITE(efu,*) '\item 65A008, \hyperref[p65A008]{\pageref{p65A008}}'  
  WRITE(efu,*) '\item 65A010, \hyperref[p65A010]{\pageref{p65A010}}'  
  WRITE(efu,*) '\item 65A012, \hyperref[p65A012]{\pageref{p65A012}}'  
  WRITE(efu,*) '\item 65A015, \hyperref[p65A015]{\pageref{p65A015}}'  

  WRITE(efu,*) '\item 66-006, \hyperref[p66-006]{\pageref{p66-006}}'  
  WRITE(efu,*) '\item 66-008, \hyperref[p66-008]{\pageref{p66-008}}'  
  WRITE(efu,*) '\item 66-009, \hyperref[p66-009]{\pageref{p66-009}}'  
  WRITE(efu,*) '\item 66-010, \hyperref[p66-010]{\pageref{p66-010}}'  
  WRITE(efu,*) '\item 66-012, \hyperref[p66-012]{\pageref{p66-012}}'  
  WRITE(efu,*) '\item 66-015, \hyperref[p66-015]{\pageref{p66-015}}'  
  WRITE(efu,*) '\item 66-018, \hyperref[p66-018]{\pageref{p66-018}}'  
  WRITE(efu,*) '\item 66-021, \hyperref[p66-021]{\pageref{p66-021}}'  


  WRITE(efu,*) '\item 67-012, \hyperref[p67-012]{\pageref{p67-012}}'  
  WRITE(efu,*) '\item 67-015, \hyperref[p67-015]{\pageref{p67-015}}'  
  
   
  WRITE(efu,FMT1) 'Meanlines'
  WRITE(efu,*) '\item 62, \hyperref[ml62]{\pageref{ml62}}'  
  WRITE(efu,*) '\item 63, \hyperref[ml63]{\pageref{ml63}}'  
  WRITE(efu,*) '\item 64, \hyperref[ml64]{\pageref{ml64}}'  
  WRITE(efu,*) '\item 65, \hyperref[ml65]{\pageref{ml65}}'  
  WRITE(efu,*) '\item 66, \hyperref[ml66]{\pageref{ml66}}'  
  WRITE(efu,*) '\item 67, \hyperref[ml67]{\pageref{ml67}}'  
  WRITE(efu,*) '\item 210, \hyperref[ml210]{\pageref{ml210}}'  
  WRITE(efu,*) '\item 220, \hyperref[ml220]{\pageref{ml220}}'  
  WRITE(efu,*) '\item 230, \hyperref[ml230]{\pageref{ml230}}'  
  WRITE(efu,*) '\item 240, \hyperref[ml240]{\pageref{ml240}}'  
  WRITE(efu,*) '\item 250, \hyperref[ml250]{\pageref{ml250}}'  
  WRITE(efu,*) '\item 430, \hyperref[ml430]{\pageref{ml430}}'  
  WRITE(efu,*) '\item 211, \hyperref[ml211]{\pageref{ml211}}'  
  WRITE(efu,*) '\item 221, \hyperref[ml221]{\pageref{ml221}}'  
  WRITE(efu,*) '\item 231, \hyperref[ml231]{\pageref{ml231}}'  
  WRITE(efu,*) '\item 241, \hyperref[ml241]{\pageref{ml241}}'  
  WRITE(efu,*) '\item 251, \hyperref[ml251]{\pageref{ml251}}'  
  
  WRITE(efu,*) '\item a=0.0, \hyperref[mla=0.0]{\pageref{mla=0.0}}'
  WRITE(efu,*) '\item a=0.1, \hyperref[mla=0.1]{\pageref{mla=0.1}}'
  WRITE(efu,*) '\item a=0.2, \hyperref[mla=0.2]{\pageref{mla=0.2}}'
  WRITE(efu,*) '\item a=0.3, \hyperref[mla=0.3]{\pageref{mla=0.3}}'
  WRITE(efu,*) '\item a=0.4, \hyperref[mla=0.4]{\pageref{mla=0.4}}'
  WRITE(efu,*) '\item a=0.5, \hyperref[mla=0.5]{\pageref{mla=0.5}}'
  WRITE(efu,*) '\item a=0.6, \hyperref[mla=0.6]{\pageref{mla=0.6}}'
  WRITE(efu,*) '\item a=0.7, \hyperref[mla=0.7]{\pageref{mla=0.7}}'
  WRITE(efu,*) '\item a=0.8, \hyperref[mla=0.8]{\pageref{mla=0.8}}'
  WRITE(efu,*) '\item a=0.9, \hyperref[mla=0.9]{\pageref{mla=0.9}}'
  WRITE(efu,*) '\item a=1.0, \hyperref[mla=1.0]{\pageref{mla=1.0}}'
  WRITE(efu,*) '\item a=0.8(modified), \hyperref[mla=0.8(modified)]{\pageref{mla=0.8(modified)}}'


  WRITE(efu,FMT1) 'Sections'  
  WRITE(efu,*) '\item 1408, \hyperref[s1408]{\pageref{s1408}}'        
  WRITE(efu,*) '\item 1410, \hyperref[s1410]{\pageref{s1410}}'    
  WRITE(efu,*) '\item 1412, \hyperref[s1412]{\pageref{s1412}}'    
  WRITE(efu,*) '\item 2408, \hyperref[s2408]{\pageref{s2408}}'    
  WRITE(efu,*) '\item 2410, \hyperref[s2410]{\pageref{s2410}}'    
  WRITE(efu,*) '\item 2412, \hyperref[s2412]{\pageref{s2412}}'    
  WRITE(efu,*) '\item 2415, \hyperref[s2415]{\pageref{s2415}}'    
  WRITE(efu,*) '\item 2418, \hyperref[s2418]{\pageref{s2418}}'    
  WRITE(efu,*) '\item 2421, \hyperref[s2421]{\pageref{s2421}}'    
  WRITE(efu,*) '\item 2424, \hyperref[s2424]{\pageref{s2424}}'    
  WRITE(efu,*) '\item 4412, \hyperref[s4412]{\pageref{s4412}}'    
  WRITE(efu,*) '\item 4415, \hyperref[s4415]{\pageref{s4415}}'    
  WRITE(efu,*) '\item 4418, \hyperref[s4418]{\pageref{s4418}}'    
  WRITE(efu,*) '\item 4421, \hyperref[s4421]{\pageref{s4421}}'    
  WRITE(efu,*) '\item 4424, \hyperref[s4424]{\pageref{s4424}}'    
  
  WRITE(efu,*) '\item 23012, \hyperref[s23012]{\pageref{s23012}}'    
  WRITE(efu,*) '\item 23015, \hyperref[s23015]{\pageref{s23015}}'     
  WRITE(efu,*) '\item 23018, \hyperref[s23018]{\pageref{s23018}}'    
  WRITE(efu,*) '\item 23021, \hyperref[s23021]{\pageref{s23021}}'     
  WRITE(efu,*) '\item 23024, \hyperref[s23024]{\pageref{s23024}}'    
  
  
  
  WRITE(efu,*) '\item 63-206, \hyperref[s63-206]{\pageref{s63-206}}'        
  WRITE(efu,*) '\item 63-209, \hyperref[s63-209]{\pageref{s63-209}}'   
  WRITE(efu,*) '\item 63-210, \hyperref[s63-209]{\pageref{s63-209}}'  
  WRITE(efu,*) '\item 63-212, \hyperref[s63-209]{\pageref{s63-209}}'    
  WRITE(efu,*) '\item 63-412, \hyperref[s63-209]{\pageref{s63-209}}'    
  WRITE(efu,*) '\item 63-215, \hyperref[s63-209]{\pageref{s63-209}}'    
  WRITE(efu,*) '\item 63-415, \hyperref[s63-209]{\pageref{s63-209}}'    
  WRITE(efu,*) '\item 63-218, \hyperref[s63-209]{\pageref{s63-209}}'    
  WRITE(efu,*) '\item 63-418, \hyperref[s63-209]{\pageref{s63-209}}'    
  WRITE(efu,*) '\item 63-618, \hyperref[s63-209]{\pageref{s63-209}}'    
  WRITE(efu,*) '\item 63-221, \hyperref[s63-209]{\pageref{s63-209}}'    
  WRITE(efu,*) '\item 63-421, \hyperref[s63-209]{\pageref{s63-209}}'    

      
  WRITE(efu,*) '\item 63A210, \hyperref[s63A210]{\pageref{s63A210}}'       
  WRITE(efu,*) '\item 63A112, \hyperref[s63A112]{\pageref{s63A112}}'       
  WRITE(efu,*) '\item 63A412, \hyperref[s63A412]{\pageref{s63A412}}'       
  WRITE(efu,*) '\item 63A415, \hyperref[s63A415]{\pageref{s63A415}}'       
  WRITE(efu,*) '\item 63A615, \hyperref[s63A615]{\pageref{s63A615}}'       
  WRITE(efu,*) '\item 63A418, \hyperref[s63A418]{\pageref{s63A418}}'       
  WRITE(efu,*) '\item 63A421, \hyperref[s63A421]{\pageref{s63A421}}'       

  WRITE(efu,*) '\item 64-108, \hyperref[s64-108]{\pageref{s64-108}}'        
  WRITE(efu,*) '\item 64-110, \hyperref[s64-110]{\pageref{s64-110}}'          
  WRITE(efu,*) '\item 64-206, \hyperref[s64-206]{\pageref{s64-206}}'        
  WRITE(efu,*) '\item 64-208, \hyperref[s64-208]{\pageref{s64-208}}'   
  WRITE(efu,*) '\item 64-209, \hyperref[s64-209]{\pageref{s64-209}}'  
  WRITE(efu,*) '\item 64-210, \hyperref[s64-210]{\pageref{s64-210}}'    
  WRITE(efu,*) '\item 64-112, \hyperref[s64-112]{\pageref{s64-112}}'    
  WRITE(efu,*) '\item 64-212, \hyperref[s64-212]{\pageref{s64-212}}'    
  WRITE(efu,*) '\item 64-412, \hyperref[s64-412]{\pageref{s64-412}}'    
  WRITE(efu,*) '\item 64-215, \hyperref[s64-215]{\pageref{s64-215}}'    
  WRITE(efu,*) '\item 64-218, \hyperref[s64-218]{\pageref{s64-218}}'    
  WRITE(efu,*) '\item 64-418, \hyperref[s64-418]{\pageref{s64-418}}'    
  WRITE(efu,*) '\item 64-618, \hyperref[s64-618]{\pageref{s63-618}}'    
  WRITE(efu,*) '\item 64-421, \hyperref[s64-421]{\pageref{s64-421}}'        




  WRITE(efu,*) '\item 64A204, \hyperref[s64A204]{\pageref{s64A204}}'       
  WRITE(efu,*) '\item 64A106, \hyperref[s64A106]{\pageref{s64A106}}'       
  WRITE(efu,*) '\item 64A210, \hyperref[s64A210]{\pageref{s64A218}}'       
  WRITE(efu,*) '\item 64A410, \hyperref[s64A410]{\pageref{s64A410}}'       
  WRITE(efu,*) '\item 64A212, \hyperref[s64A212]{\pageref{s64A212}}'       
  WRITE(efu,*) '\item 64A612, \hyperref[s64A612]{\pageref{s64A612}}'       
  WRITE(efu,*) '\item 64A114, \hyperref[s64A114]{\pageref{s64A114}}'       
  WRITE(efu,*) '\item 64A215, \hyperref[s64A215]{\pageref{s64A215}}'       
  WRITE(efu,*) '\item 64A415, \hyperref[s64A415]{\pageref{s64A415}}'       
  WRITE(efu,*) '\item 64A218, \hyperref[s64A218]{\pageref{s64A218}}'       
  WRITE(efu,*) '\item 64A221, \hyperref[s64A221]{\pageref{s64A221}}'       
  
  WRITE(efu,*) '\item 65-108, \hyperref[s65-108]{\pageref{s65-108}}'        
  WRITE(efu,*) '\item 65-110, \hyperref[s65-110]{\pageref{s65-110}}'          
  WRITE(efu,*) '\item 65-206, \hyperref[s65-206]{\pageref{s65-206}}'        
  WRITE(efu,*) '\item 65-208, \hyperref[s65-208]{\pageref{s65-208}}'   
  WRITE(efu,*) '\item 65-209, \hyperref[s65-209]{\pageref{s65-209}}'  
  WRITE(efu,*) '\item 65-210, \hyperref[s65-210]{\pageref{s65-210}}'   
  WRITE(efu,*) '\item 65-410, \hyperref[s65-410]{\pageref{s65-210}}'      
  WRITE(efu,*) '\item 65-112, \hyperref[s65-112]{\pageref{s65-112}}'    
  WRITE(efu,*) '\item 65-212, \hyperref[s65-212]{\pageref{s65-212}}'  
  WRITE(efu,*) '\item 65-212a=0.5, \hyperref[s65-212a=0.5]{\pageref{s65-212a=0.5}}'    
  WRITE(efu,*) '\item 65-412, \hyperref[s65-412]{\pageref{s65-412}}'    
  WRITE(efu,*) '\item 65-215, \hyperref[s65-215]{\pageref{s65-215}}'  
  WRITE(efu,*) '\item 65-415, \hyperref[s65-415]{\pageref{s65-415}}'  
  WRITE(efu,*) '\item 65-415a=0.5, \hyperref[s65-415a=0.5]{\pageref{s65-415a=0.5}}'      
  WRITE(efu,*) '\item 65-218, \hyperref[s65-218]{\pageref{s65-218}}'    
  WRITE(efu,*) '\item 65-418, \hyperref[s65-418]{\pageref{s65-418}}' 
  WRITE(efu,*) '\item 65-418a=0.5, \hyperref[s65-418a=0.5]{\pageref{s65-418a=0.5}}'      
  WRITE(efu,*) '\item 65-618, \hyperref[s65-618]{\pageref{s65-618}}'  
  WRITE(efu,*) '\item 65-618a=0.5, \hyperref[s65-618a=0.5]{\pageref{s65-618a=0.5}}'     
  WRITE(efu,*) '\item 65-421, \hyperref[s65-421]{\pageref{s65-421}}'        
  WRITE(efu,*) '\item 65-421a=0.5, \hyperref[s65-421a=0.5]{\pageref{s65-421a=0.5}}' 
  
  
  WRITE(efu,*) '\item 65A210, \hyperref[s65A210]{\pageref{s65A210}}'       
  WRITE(efu,*) '\item 65A112, \hyperref[s65A112]{\pageref{s65A112}}'       
  WRITE(efu,*) '\item 65A212, \hyperref[s65A212]{\pageref{s65A212}}'       
  WRITE(efu,*) '\item 65A412, \hyperref[s65A412]{\pageref{s65A412}}'       
  WRITE(efu,*) '\item 65A215, \hyperref[s65A215]{\pageref{s65A215}}'       
  WRITE(efu,*) '\item 65A415, \hyperref[s65A415]{\pageref{s65A415}}'       
  WRITE(efu,*) '\item 65A615, \hyperref[s65A615]{\pageref{s65A615}}'       
  WRITE(efu,*) '\item 65A418, \hyperref[s65A418]{\pageref{s65A418}}'       
  WRITE(efu,*) '\item 65A421, \hyperref[s65A421]{\pageref{s65A421}}'       
  
  WRITE(efu,*) '\item 66-108, \hyperref[s66-108]{\pageref{s66-108}}'        
  WRITE(efu,*) '\item 66-110, \hyperref[s66-110]{\pageref{s66-110}}'          
  WRITE(efu,*) '\item 66-206, \hyperref[s66-206]{\pageref{s66-206}}'        
  WRITE(efu,*) '\item 66-208, \hyperref[s66-208]{\pageref{s66-208}}'   
  WRITE(efu,*) '\item 66-209, \hyperref[s66-209]{\pageref{s66-209}}'  
  WRITE(efu,*) '\item 66-210, \hyperref[s66-210]{\pageref{s66-210}}'    
  WRITE(efu,*) '\item 66-112, \hyperref[s66-112]{\pageref{s66-112}}'    
  WRITE(efu,*) '\item 66-212, \hyperref[s66-212]{\pageref{s66-212}}'    
  WRITE(efu,*) '\item 66-412, \hyperref[s66-412]{\pageref{s66-412}}'    
  WRITE(efu,*) '\item 66-215, \hyperref[s66-215]{\pageref{s66-215}}'   
  WRITE(efu,*) '\item 66-415, \hyperref[s66-415]{\pageref{s66-415}}'    
  WRITE(efu,*) '\item 66-218, \hyperref[s66-218]{\pageref{s66-218}}'    
  WRITE(efu,*) '\item 66-418, \hyperref[s66-418]{\pageref{s66-418}}'    
  WRITE(efu,*) '\item 66-618, \hyperref[s66-618]{\pageref{s63-618}}'    
  WRITE(efu,*) '\item 66-221, \hyperref[s66-221]{\pageref{s66-221}}'      
  WRITE(efu,*) '\item 66-421, \hyperref[s66-421]{\pageref{s66-421}}'        
 
  WRITE(efu,*) '\item 67-215, \hyperref[s67-212]{\pageref{s67-212}}'     
  WRITE(efu,*) '\item 67-215, \hyperref[s67-215]{\pageref{s67-215}}'       
   
  WRITE(efu,*) '\end{theindex}'
  WRITE(efu,*) '\onecolumn'
  WRITE(efu,*) ' '
  WRITE(efu,*) '\newpage'
  RETURN
END Subroutine MakeIndex   ! ---------------------------------------------------    

!+
SUBROUTINE Profiles(efu, x)
! ------------------------------------------------------------------------------
! PURPOSE - Add the tables of profiles - Appendix I - to the TeX file

  INTEGER,INTENT(IN):: efu
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coor for calculation

!-------------------------------------------------------------------------------

  
  CALL WriteFourDigitProfile(efu, 0.06, '0006', x)
  CALL WriteFourDigitProfile(efu, 0.08, '0008', x)
  CALL WriteFourDigitProfile(efu, 0.10, '0010', x)
  CALL WriteFourDigitProfile(efu, 0.12, '0012', x)
  CALL WriteFourDigitProfile(efu, 0.15, '0015', x)
  CALL WriteFourDigitProfile(efu, 0.18, '0018', x)
  CALL WriteFourDigitProfile(efu, 0.21, '0021', x)

  CALL WriteModifiedFourDigitProfile(efu, 3.0, 0.4, 0.08, '0008-34', x)
  CALL WriteModifiedFourDigitProfile(efu, 3.0, 0.4, 0.10, '0010-34', x)
  CALL WriteModifiedFourDigitProfile(efu, 3.0, 0.5, 0.10, '0010-35', x)
  CALL WriteModifiedFourDigitProfile(efu, 6.0, 0.4, 0.10, '0010-64', x)
  CALL WriteModifiedFourDigitProfile(efu, 6.0, 0.5, 0.10, '0010-65', x)
  CALL WriteModifiedFourDigitProfile(efu, 6.0, 0.6, 0.10, '0010-66', x)
  CALL WriteModifiedFourDigitProfile(efu, 3.0, 0.4, 0.12, '0012-34', x)
  CALL WriteModifiedFourDigitProfile(efu, 6.0, 0.4, 0.12, '0012-64', x)

  CALL WriteModifiedFourDigitProfile(efu, 4.0, 0.5, 0.06, '16-006', x)
  CALL WriteModifiedFourDigitProfile(efu, 4.0, 0.5, 0.09, '16-009', x)
  CALL WriteModifiedFourDigitProfile(efu, 4.0, 0.5, 0.12, '16-012', x)
  CALL WriteModifiedFourDigitProfile(efu, 4.0, 0.5, 0.15, '16-015', x)
  CALL WriteModifiedFourDigitProfile(efu, 4.0, 0.5, 0.18, '16-018', x)
  CALL WriteModifiedFourDigitProfile(efu, 4.0, 0.5, 0.21, '16-021', x)

  CALL WriteSixSeriesProfile(efu, 1, 0.06, '63-006', x)
  CALL WriteSixSeriesProfile(efu, 1, 0.08, '63-008', x)
  CALL WriteSixSeriesProfile(efu, 1, 0.09, '63-009', x)
  CALL WriteSixSeriesProfile(efu, 1, 0.10, '63-010', x)
  CALL WriteSixSeriesProfile(efu, 1, 0.12, '63-012', x)
  CALL WriteSixSeriesProfile(efu, 1, 0.15, '63-015', x)
  CALL WriteSixSeriesProfile(efu, 1, 0.18, '63-018', x)
  CALL WriteSixSeriesProfile(efu, 1, 0.21, '63-021', x)

  CALL WriteSixSeriesProfile(efu, 6, 0.06, '63A006', x)
  CALL WriteSixSeriesProfile(efu, 6, 0.08, '63A008', x)
  CALL WriteSixSeriesProfile(efu, 6, 0.10, '63A010', x)
  CALL WriteSixSeriesProfile(efu, 6, 0.12, '63A012', x)
  CALL WriteSixSeriesProfile(efu, 6, 0.15, '63A015', x)

  CALL WriteSixSeriesProfile(efu, 2, 0.06, '64-006', x)
  CALL WriteSixSeriesProfile(efu, 2, 0.08, '64-008', x)
  CALL WriteSixSeriesProfile(efu, 2, 0.09, '64-009', x)
  CALL WriteSixSeriesProfile(efu, 2, 0.10, '64-010', x)
  CALL WriteSixSeriesProfile(efu, 2, 0.12, '64-012', x)
  CALL WriteSixSeriesProfile(efu, 2, 0.15, '64-015', x)
  CALL WriteSixSeriesProfile(efu, 2, 0.18, '64-018', x)
  CALL WriteSixSeriesProfile(efu, 2, 0.21, '64-021', x)

  CALL WriteSixSeriesProfile(efu, 7, 0.06, '64A006', x)
  CALL WriteSixSeriesProfile(efu, 7, 0.08, '64A008', x)
  CALL WriteSixSeriesProfile(efu, 7, 0.10, '64A010', x)
  CALL WriteSixSeriesProfile(efu, 7, 0.12, '64A012', x)
  CALL WriteSixSeriesProfile(efu, 7, 0.15, '64A015', x)

  CALL WriteSixSeriesProfile(efu, 3, 0.06, '65-006', x)
  CALL WriteSixSeriesProfile(efu, 3, 0.08, '65-008', x)
  CALL WriteSixSeriesProfile(efu, 3, 0.09, '65-009', x)
  CALL WriteSixSeriesProfile(efu, 3, 0.10, '65-010', x)
  CALL WriteSixSeriesProfile(efu, 3, 0.12, '65-012', x)
  CALL WriteSixSeriesProfile(efu, 3, 0.15, '65-015', x)
  CALL WriteSixSeriesProfile(efu, 3, 0.18, '65-018', x)
  CALL WriteSixSeriesProfile(efu, 3, 0.21, '65-021', x)

  CALL WriteSixSeriesProfile(efu, 8, 0.06, '65A006', x)
  CALL WriteSixSeriesProfile(efu, 8, 0.08, '65A008', x)
  CALL WriteSixSeriesProfile(efu, 8, 0.10, '65A010', x)
  CALL WriteSixSeriesProfile(efu, 8, 0.12, '65A012', x)
  CALL WriteSixSeriesProfile(efu, 8, 0.15, '65A015', x)

  CALL WriteSixSeriesProfile(efu, 4, 0.06, '66-006', x)
  CALL WriteSixSeriesProfile(efu, 4, 0.08, '66-008', x)
  CALL WriteSixSeriesProfile(efu, 4, 0.09, '66-009', x)
  CALL WriteSixSeriesProfile(efu, 4, 0.10, '66-010', x)
  CALL WriteSixSeriesProfile(efu, 4, 0.12, '66-012', x)
  CALL WriteSixSeriesProfile(efu, 4, 0.15, '66-015', x)
  CALL WriteSixSeriesProfile(efu, 4, 0.18, '66-018', x)
  CALL WriteSixSeriesProfile(efu, 4, 0.21, '66-021', x)

  CALL WriteSixSeriesProfile(efu, 5, 0.12, '67-012', x)
  CALL WriteSixSeriesProfile(efu, 5, 0.15, '67-015', x)

  RETURN
END Subroutine Profiles   ! ----------------------------------------------------

!+
SUBROUTINE MeanLines(efu,x)
! ------------------------------------------------------------------------------
! PURPOSE - Add the tables of mean lines - Appendix II - to the TeX file
  INTEGER,INTENT(IN):: efu  ! external file unit of the TEX file
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coor for calculation
!-------------------------------------------------------------------------------

  CALL WriteTwoDigitMeanLine(efu, 0.06, 0.2, '62', x)
  CALL WriteTwoDigitMeanLine(efu, 0.06, 0.3, '63', x)
  CALL WriteTwoDigitMeanLine(efu, 0.06, 0.4, '64', x)
  CALL WriteTwoDigitMeanLine(efu, 0.06, 0.5, '65', x)
  CALL WriteTwoDigitMeanLine(efu, 0.06, 0.6, '66', x)
  CALL WriteTwoDigitMeanLine(efu, 0.06, 0.7, '67', x)

  CALL WriteThreeDigitMeanLine(efu, 0.3, 0.05, '210', x)
  CALL WriteThreeDigitMeanLine(efu, 0.3, 0.10, '220', x)
  CALL WriteThreeDigitMeanLine(efu, 0.3, 0.15, '230', x)
  CALL WriteThreeDigitMeanLine(efu, 0.3, 0.20, '240', x)
  CALL WriteThreeDigitMeanLine(efu, 0.3, 0.25, '250', x)
  CALL WriteThreeDigitMeanLine(efu, 0.6, 0.15, '430', x)

  CALL WriteThreeDigitReflexMeanLine(efu, 0.3, 0.05, '211', x)
  CALL WriteThreeDigitReflexMeanLine(efu, 0.3, 0.10, '221', x)
  CALL WriteThreeDigitReflexMeanLine(efu, 0.3, 0.15, '231', x)
  CALL WriteThreeDigitReflexMeanLine(efu, 0.3, 0.20, '241', x)
  CALL WriteThreeDigitReflexMeanLine(efu, 0.3, 0.25, '251', x)

  CALL WriteSixSeriesMeanLine(efu,0.0,1.0,'a=0.0',x)
  CALL WriteSixSeriesMeanLine(efu,0.1,1.0,'a=0.1',x)
  CALL WriteSixSeriesMeanLine(efu,0.2,1.0,'a=0.2',x)
  CALL WriteSixSeriesMeanLine(efu,0.3,1.0,'a=0.3',x)
  CALL WriteSixSeriesMeanLine(efu,0.4,1.0,'a=0.4',x)
  CALL WriteSixSeriesMeanLine(efu,0.5,1.0,'a=0.5',x)
  CALL WriteSixSeriesMeanLine(efu,0.6,1.0,'a=0.6',x)
  CALL WriteSixSeriesMeanLine(efu,0.7,1.0,'a=0.7',x)
  CALL WriteSixSeriesMeanLine(efu,0.8,1.0,'a=0.8',x)
  CALL WriteSixSeriesMeanLine(efu,0.9,1.0,'a=0.9',x)
  CALL WriteSixSeriesMeanLine(efu,1.0,1.0,'a=1.0',x)

  CALL WriteSixSeriesModifiedMeanLine(efu,1.0,'a=0.8(modified)',x)

  RETURN
END Subroutine MeanLines   ! ---------------------------------------------------

!+
SUBROUTINE Sections45(efu,x)
! ------------------------------------------------------------------------------
! PURPOSE - Add the tables of profiles - Appendix I
!   This subroutine writes the 4 and 5 digit sections.
  INTEGER,INTENT(IN):: efu  ! external file unit of the TEX file
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coor for calculation
!-------------------------------------------------------------------------------

  CALL WriteFourDigitSection(efu, 0.01,0.4,0.08, '1408', x)
  CALL WriteFourDigitSection(efu, 0.01,0.4,0.10, '1410', x)
  CALL WriteFourDigitSection(efu, 0.01,0.4,0.12, '1412', x)
  CALL WriteFourDigitSection(efu, 0.02,0.4,0.08, '2408', x)
  CALL WriteFourDigitSection(efu, 0.02,0.4,0.10, '2410', x)
  CALL WriteFourDigitSection(efu, 0.02,0.4,0.12, '2412', x)
  CALL WriteFourDigitSection(efu, 0.02,0.4,0.15, '2415', x)
  CALL WriteFourDigitSection(efu, 0.02,0.4,0.18, '2418', x)
  CALL WriteFourDigitSection(efu, 0.02,0.4,0.21, '2421', x)
  CALL WriteFourDigitSection(efu, 0.02,0.4,0.24, '2424', x)
  CALL WriteFourDigitSection(efu, 0.04,0.4,0.12, '4412', x)
  CALL WriteFourDigitSection(efu, 0.04,0.4,0.15, '4415', x)
  CALL WriteFourDigitSection(efu, 0.04,0.4,0.18, '4418', x)
  CALL WriteFourDigitSection(efu, 0.04,0.4,0.21, '4421', x)
  CALL WriteFourDigitSection(efu, 0.04,0.4,0.24, '4424', x)

  CALL WriteFiveDigitSection(efu, 0.3,0.15,0.12, '23012', x)
  CALL WriteFiveDigitSection(efu, 0.3,0.15,0.15, '23015', x)
  CALL WriteFiveDigitSection(efu, 0.3,0.15,0.18, '23018', x)
  CALL WriteFiveDigitSection(efu, 0.3,0.15,0.21, '23021', x)
  CALL WriteFiveDigitSection(efu, 0.3,0.15,0.24, '23024', x)

  RETURN
END Subroutine Sections45   ! ----------------------------------------------------

!+
SUBROUTINE Sections6(efu,x)
! ------------------------------------------------------------------------------
! PURPOSE - Add the tables of sections - Appendix III
  INTEGER,INTENT(IN):: efu  ! external file unit of the TEX file
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coor for calculation
!-------------------------------------------------------------------------------

  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.2, 0.06, '63-206', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.2, 0.09, '63-209', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.2, 0.10, '63-210', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.2, 0.12, '63-212', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.4, 0.12, '63-412', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.2, 0.15, '63-215', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.4, 0.15, '63-415', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.6, 0.15, '63-615', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.2, 0.18, '63-218', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.4, 0.18, '63-418', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.6, 0.18, '63-618', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.2, 0.21, '63-221', x)
  CALL WriteSixSeriesSection(efu, 1, 1.0, 0.4, 0.21, '63-421', x)

  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.1, 0.08, '64-108', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.1, 0.10, '64-110', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.2, 0.06, '64-206', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.2, 0.08, '64-208', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.2, 0.09, '64-209', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.2, 0.10, '64-210', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.1, 0.12, '64-112', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.2, 0.12, '64-212', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.4, 0.12, '64-412', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.2, 0.15, '64-215', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.4, 0.15, '64-415', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.2, 0.18, '64-218', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.4, 0.18, '64-418', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.6, 0.18, '64-618', x)
  CALL WriteSixSeriesSection(efu, 2, 1.0, 0.4, 0.21, '64-421', x)
  
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.1, 0.08, '65-108', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.1, 0.10, '65-110', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.2, 0.06, '65-206', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.2, 0.08, '65-208', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.2, 0.09, '65-209', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.2, 0.10, '65-210', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.4, 0.10, '65-410', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.1, 0.12, '65-112', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.2, 0.12, '65-212', x)
  CALL WriteSixSeriesSection(efu, 3, 0.5, 0.2, 0.12, '65-212a=0.5', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.4, 0.12, '65-412', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.2, 0.15, '65-215', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.4, 0.15, '65-415', x)
  CALL WriteSixSeriesSection(efu, 3, 0.5, 0.4, 0.15, '65-415a=0.5', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.2, 0.18, '65-218', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.4, 0.18, '65-418', x)
  CALL WriteSixSeriesSection(efu, 3, 0.5, 0.4, 0.18, '65-418a=0.5', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.6, 0.18, '65-618', x)
  CALL WriteSixSeriesSection(efu, 3, 0.5, 0.6, 0.18, '65-618a=0.5', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.2, 0.21, '65-221', x)
  CALL WriteSixSeriesSection(efu, 3, 1.0, 0.4, 0.21, '65-421', x)
  CALL WriteSixSeriesSection(efu, 3, 0.5, 0.4, 0.21, '65-421a=0.5', x)
   
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.1, 0.08, '66-108', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.1, 0.10, '66-110', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.2, 0.06, '66-206', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.2, 0.08, '66-208', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.2, 0.09, '66-209', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.2, 0.10, '66-210', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.1, 0.12, '66-112', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.2, 0.12, '66-212', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.4, 0.12, '66-412', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.2, 0.15, '66-215', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.4, 0.15, '66-415', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.2, 0.18, '66-218', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.4, 0.18, '66-418', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.6, 0.18, '66-618', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.2, 0.21, '66-221', x)
  CALL WriteSixSeriesSection(efu, 4, 1.0, 0.4, 0.21, '66-421', x)
  
  CALL WriteSixSeriesSection(efu, 5, 1.0, 0.2, 0.12, '67-212', x)
  CALL WriteSixSeriesSection(efu, 5, 1.0, 0.2, 0.15, '67-215', x)
 
  RETURN
END Subroutine Sections6   ! ----------------------------------------------------

!+
SUBROUTINE Sections6A(efu,x)
! ------------------------------------------------------------------------------
! PURPOSE - Create a LaTeX file with the tables of sections - Appendix III
  INTEGER,INTENT(IN):: efu  ! external file unit of the TEX file
  REAL,INTENT(IN),DIMENSION(:):: x   ! table of x-coor for calculation
!-------------------------------------------------------------------------------
  
  CALL WriteSixSeriesSection(efu, 6, 0.8, 0.2, 0.10, '63A210', x)
  CALL WriteSixSeriesSection(efu, 6, 0.8, 0.1, 0.12, '63A112', x)
  CALL WriteSixSeriesSection(efu, 6, 0.8, 0.4, 0.12, '63A412', x)
  CALL WriteSixSeriesSection(efu, 6, 0.8, 0.4, 0.15, '63A415', x)
  CALL WriteSixSeriesSection(efu, 6, 0.8, 0.6, 0.15, '63A615', x)
  CALL WriteSixSeriesSection(efu, 6, 0.8, 0.4, 0.18, '63A418', x)
  CALL WriteSixSeriesSection(efu, 6, 0.8, 0.4, 0.21, '63A421', x)

  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.2, 0.04, '64A204', x)
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.1, 0.06, '64A106', x)  
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.2, 0.10, '64A210', x)
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.4, 0.10, '64A410', x)
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.2, 0.12, '64A212', x)
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.4, 0.12, '64A412', x)
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.6, 0.12, '64A612', x)
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.1, 0.14, '64A114', x)
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.2, 0.15, '64A215', x)
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.4, 0.15, '64A415', x)
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.2, 0.18, '64A218', x)
  CALL WriteSixSeriesSection(efu, 7, 0.8, 0.2, 0.21, '64A221', x)
  CALL WriteSixSeriesSection(efu, 8, 0.8, 0.2, 0.10, '65A210', x)
  CALL WriteSixSeriesSection(efu, 8, 0.8, 0.2, 0.12, '65A112', x)
  CALL WriteSixSeriesSection(efu, 8, 0.8, 0.2, 0.12, '65A212', x)
  CALL WriteSixSeriesSection(efu, 8, 0.8, 0.4, 0.12, '65A412', x)
  CALL WriteSixSeriesSection(efu, 8, 0.8, 0.2, 0.15, '65A215', x)
  CALL WriteSixSeriesSection(efu, 8, 0.8, 0.4, 0.15, '65A415', x)
  CALL WriteSixSeriesSection(efu, 8, 0.8, 0.6, 0.15, '65A615', x)
  CALL WriteSixSeriesSection(efu, 8, 0.8, 0.4, 0.18, '65A418', x)
  CALL WriteSixSeriesSection(efu, 8, 0.8, 0.4, 0.21, '65A421', x)
  
  RETURN
END Subroutine Sections6A   ! --------------------------------------------------

!+
SUBROUTINE BackMatter(efu)
! ------------------------------------------------------------------------------
! PURPOSE - Write all of the material following the tables

  INTEGER,INTENT(IN):: efu   ! external file unit of the TeX file  
!-------------------------------------------------------------------------------

  WRITE(efu,*) "\end{document}"

  RETURN
END Subroutine BackMatter   ! --------------------------------------------------

END Module LaTeXProcedures   ! =================================================


!+
PROGRAM AbbottVonDoenhoffPDF                           ! pdas/naca456/avdpdf.f90
! ------------------------------------------------------------------------------
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

USE ISO_FORTRAN_ENV
USE NACAauxilary,ONLY: LoadX
USE LaTeXProcedures,ONLY: FrontMatter, BackMatter, Profiles, MeanLines, &
    Sections45, Sections6, Sections6A
IMPLICIT NONE
  
  CHARACTER(LEN=*),PARAMETER:: AUTHOR = &
    "Ralph L. Carmichael, Public Domain Aeronautical Software"
  CHARACTER(LEN=*),PARAMETER:: GREETING="avdpdf - create avd.pdf"
  CHARACTER(LEN=*),PARAMETER:: MODIFIER=" "  ! your name if you change it

  INTEGER,PARAMETER:: denCode=1   ! tells LoadX to use coarse spacing
  INTEGER:: nx                    ! length of x  (from LoadX)
  REAL,DIMENSION(100):: x         ! x-coor for calculation (from LoadX)
                                  ! 100 is more than enough
  INTEGER,PARAMETER:: TEX = 1
  CHARACTER(LEN=*),PARAMETER:: VERSION=" 0.42 (2022 May 17)"
  REAL:: time1,time2
!-------------------------------------------------------------------------------
  CALL CPU_TIME(time1)
  WRITE(*,*) GREETING
  WRITE(*,*) AUTHOR
  IF (MODIFIER /= ' ') WRITE(*,*) 'Modified by '//MODIFIER
  WRITE(*,*) VERSION
  WRITE(*,'(3A/A)') 'This file was compiled by ', compiler_version(),   &
                    ' using the options ',        compiler_options()

  OPEN(UNIT=TEX, FILE='avd.tex', STATUS='REPLACE', ACTION='WRITE')
  
  CALL FrontMatter(TEX)
  CALL LoadX(denCode, nx,x)

  CALL Profiles(TEX, x(1:nx))
  CALL MeanLines(TEX, x(1:nx))
  CALL Sections45(TEX, x(1:nx))
  CALL Sections6(TEX, x(1:nx))
  CALL Sections6A(TEX, x(1:nx))

  CALL BackMatter(TEX) 
  CLOSE(UNIT=TEX)

  WRITE(*,*) 'Normal termination of program avdpdf.'
  WRITE(*,*) 'File avd.tex added to your directory'
  CALL SYSTEM('pdflatex avd.tex')
  CALL CPU_TIME(time2)
  WRITE(*,'(A,I0,A)') ' CPU time= ', NINT(1000.0*(time2-time1)), ' milliseconds'
  STOP
END Program AbbottVonDoenhoffPDF   ! ===========================================
