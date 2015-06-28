C***********************************************************************
C
C	TITLE: 	GEN_RLNINT  SINGLE LINEAR INTERPOLATION
C
C	PREPARED BY:	SJ Fournier
C			IPSS for DOFASCO
C
C***********************************************************************
C
C	DEVELOPED ON: DEC PDP 11/70
C
C		  IN: DEC FORTRAN 77
C
C***********************************************************************
C
C	DATE :		12-August-87
C
C	REVISIONS
C	
C***********************************************************************
C
C 		GEN_GEN_RLNINT FUNCTION DESCRITPION
C
C***********************************************************************
C 
C	GEN_RLNINT is a function which will, given X,
C	compute F(X).  It requires a table of X values and a table of
C	F(X) values.  The X values can be in ascending or descending order.
C	Input X must be in the same units as the values in the
C	X table.
C
C	All values except LEN are to be in single precision real format.
C	If input X falls outside the X table range, one of two things can 
C	occur depending on the shape of the X array.  If there are no equal 
C	values of X then F(X) will be extrapolated along a straight line 
C	extending thru the last two points.  If the last two X points are the 
C	same, F(X) will be returned as the last value in the F(X) array.  In 
C	this way, a clamp function may be used to prevent extrapolation
C
C	NOTE: GEN_RLNINT IS A LOCAL FUNCTION.
C***********************************************************************


C	Function Format (FORTRAN):
C		FX = GEN_RLNINT(X,XARY,FXARY,ARRY_LEN)
C		where
C			X	= Number to be interpolated
C			XARY	= Name of X array
C			FXARY	= Name of F(X) array
C			ARRY_LEN	= Dimension of two arrays - Integer format

	REAL*8	FUNCTION GEN_RLNINT(X,XARY,FXARY,ARRY_LEN)

	INTEGER*2	ARRY_LEN
	INTEGER*2	IBEG,IEND,ISTEP
	INTEGER*2	J,JX

	REAL*8	X,XARY(ARRY_LEN),FXARY(ARRY_LEN)

	IF(XARY(1) .LT. XARY(ARRY_LEN)) THEN
		IBEG = 2
		IEND = ARRY_LEN
		ISTEP = 1
	ELSE
		IBEG = ARRY_LEN-1
		IEND = 1
		ISTEP = -1
	END IF

D	WRITE(5,*)XARY(1),XARY(ARRY_LEN)
D	WRITE(5,*)IBEG,IEND,ISTEP

C	Save JX as the end element in the array (either 1 or ARRY_LEN)
	JX = IBEG - ISTEP

	DO 10 J = IBEG,IEND,ISTEP

C		If the input value less than or equal to the Jth value
C		of the XARY go to interpolate
		IF(X .LE. XARY(J)) GO TO 50

C		Is this not the last loop
		IF(J .NE. IEND) THEN

C			If the Jth element is not equal to the J+ISTEPth element
C			Save JX as the last element before identical X elements
C			It is the normal case for JX to be saved on each loop
			IF(XARY(J) .NE. XARY(J+ISTEP)) THEN
				JX = J
			END IF
		END IF
10	CONTINUE

C	X input is beyond the bounds of the array
C	Reset J to the end of the array
	J = IEND

C	Are the last two X elements the same value
C	JX=IEND-ISTEP if the elements are different
	IF(JX .NE. IEND-ISTEP) THEN

C		Clamp the result to the last f(x)
		GEN_RLNINT = FXARY(IEND)
		RETURN
	END IF

50	CONTINUE

D	WRITE(5,*)JX,XARY(JX),FXARY(JX),J,XARY(J),FXARY(J)

C	Calculate the interpolation
C	... first check if first two X array elements are the same and if so,
C	    move to the third element of the X and Y arrays to calculate the
C	    slope ...
	IF ( XARY(J) .EQ. XARY(JX) ) THEN
		GEN_RLNINT = FXARY(JX) + (FXARY(J+ISTEP) - FXARY(JX))
     &			*(X-XARY(JX))/(XARY(J+ISTEP) - XARY(JX))
		GO TO 60
	END IF


	GEN_RLNINT = FXARY(JX) + (FXARY(J) - FXARY(JX))
     &		*(X-XARY(JX))/(XARY(J) - XARY(JX))

60	CONTINUE

	RETURN
	END
