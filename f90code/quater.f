cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Objetive: Program to generate quaternion transform power spectrum, from a starting frequency (fmin) until a maximum frequency (fmax). (Output in ps.dgt) 
c Input Parameters: 
c        file: Input file with data column
c        fmin: starting frequency
c        fmax: maximum frequency (fmax)
c        days: number of days recalculted in based of a power 2
c        sampling: in second
c Input data file: Radial velocity measurement (excerpt). 1 column. Units: m/s
c Output data file: ps.dtg
c        7 columns
c        column 1: Frequencies in microHz
c        column 2: Classical Power
c	 column 3: Quaternion Power

c  Author: Rafael Garrido Haba
c  Supervised: Jose Ramon Rodon
c  2022 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c----------------
c  Variables
c----------------

        double precision dpi,wo,wp,th,gr,st,rci1,rci2,rci3,rmi
        real AM(8388608),TRADAT(8388608),FDAT(8388608), tradatg(8388608)
        CHARACTER*200 OBJ,OBJF,FICHE1,FICHE2
c -------------------
c Inicialization
c ------------------

        dpi=datan(1.d0)*8
        amedi=0.d0
c        obj='golfsel2.dat'
c        tpas=20.
	tradatg=0.0

        write(*,*) 'file, initial frequency, final frequency, days,
     1  sampling(sec)'
c	write(*,*) 'fmin,fmax,rti' 
c       rti max 2616 days of continuos 20 sec golf@soho data
        read(*,*) obj, anuino, anufio,rti,tpas
c	read(*,*) anuino,anufio,rti!,fmax!,res
        n=rti*86400d0/(1*tpas)
        nw=log(n*1.)/log(2.)



	
c-------------------------------------------------
c
c   Input data normalization (mean 0 variance 1)
c
c-------------------------------------------------

	OPEN(UNIT=1,FILE=OBJ)
	
        DO I=1,n
	   READ(1,*,END=1811,err=1811) am(i)
           AMEDI=AMEDI+AM(I)
	END DO
	
	goto 1812
1811    write(*,*) i,po,am(i)
	stop
1812	N=I-1
	CLOSE(UNIT=1)

	AMEDI=AMEDI/(N)
        sig = 0.d0

	DO I=1,N
          sig = sig + (AM(I)-AMEDI)**2
        END DO

	sig = sqrt(sig/N)

c	DO I=1,N
c	  am((I))=(AM(I)-AMEDI)/sig
c	END DO
	
c       Nyquist
	fny=1000000./(2*tpas)
c       Fourier fundamental
	delnu=1.d6/((n)*tpas)
	ji=floor(anuino/delnu)
	anuin=ji*delnu
	jf=floor(anufio/delnu)
c       NP: Number of frequencies which must to calculate the program
	NP=jf-ji
	write(*,*) 'Frequencies, data, days, Ny'
	WRITE(*,*) NP,n,n*tpas/86400.,fny
c       Solution: Number of days for analysis.
	F=ANUIN
	BMAX=0.d0

	DO J=1,NP
	  rci=0.d0
	  cci=0.0
c         Sums for a determinated frequency.
   	  rmi=0.0
	  rci1=0.0
	  rci2=0.0
          rci3=0.0
	  DO I=1,N
  	    th1=dmod((i)*1.d0*(j+ji)/(n),1.d0)!+gora
c           Sin/Cos for standart fourier
  	    gr=th1*dpi
	    c1=dcos(gr)
	    s1=dsin(gr)
c           Change of variable
	    ff=exp(asinh(th1/2.0))
c           two frequencies
	    wo=dpi*(ff-1.d0/ff)/2
	    wp=dpi*(ff+1.d0/ff)/2
	
            z1=dcos(wp)*dcos(wo)
 	    z2=dsin(wp)*dcos(wo)
	    z3=dcos(wp)*dsin(wo)
	    z4=dsin(wp)*dsin(wo)

 	    rci=rci+am((i))*(c1)
	    cci=cci+am((i))*(s1)
c           rmi: Real part of quaternion transform 
            rmi=rmi+am((I))*z1
c           rci1, rci2, rci3: Three imaginary parts
            rci2=rci2+am((i))*z2
            rci3=rci3+am((i))*z3
            rci1=rci1+am((i))*z4
	  END DO

c------------------------------
c
c Normalization
c
c------------------------------	

	  rci=rci*2.d0/n !real clasica
	  cci=cci*2.d0/n!imag clasica
          rmi=rmi*2.d0/n !real nueva	  
	  rci2=rci2*2.d0/n !imag nueva
	  rci3=rci3*2.d0/n!imag nueva
	  rci1=rci1*2.d0/n!imag nueva

c         B: the classical power
	  b=((rci)**2+(cci)**2)
c         Bg: The quaternionic power
          Bg=rmi**2+(rci2**2+rci1**2+rci3**2)

          F=(ANUIN+(J-1)*DELNU)
	    IF (bg.gt.BMAX) then
	      FMAX=f
	      BMAX=Bg
	      jmax=j
	    end IF

	  FDAT(J)=SNGL(F)
	  TRADAT(J)=SNGL(B)
	  TRADATg(J)=SNGL(Bg)
	END DO

c In file ps.dgt frequency, squared modulus
        fiche1='ps.dgt'

        OPEN(UNIT=41,FILE=FICHE1)
        DO I=1,NP
	  WRITE(41,*) FDAT((I)),tradat(i),TRADATg(I)
	END DO
     	WRITE(*,*) "data, fourier freq, freq max, max valor"
	WRITE(*,*) il,delnu,fdat(jmax),BMAX
	close(41)
22	STOP
	END
c-------------------------------------------------------------c
c                                                             c
c  Subroutine sffteu( x, y, n, m, itype )                     c
c                                                             c
c  This routine is a slight modification of a complex split   c
c  radix FFT routine presented by C.S. Burrus.  The original  c
c  program header is shown below.                             c
c                                                             c
c  Arguments:                                                 c
c     x - real array containing real parts of transform       c
c              sequence (in/out)                              c
c     y - real array containing imag parts of transform       c
c              sequence (in/out)                              c
c     n - integer length of transform (in)                    c
c     m - integer such that n = 2**m  (in)                    c
c     itype - integer job specifier (in)                      c
c              itype .ne. -1 --> foward transform             c
c              itype .eq. -1 --> backward transform           c
c                                                             c
c  The forward transform computes                             c
c     X(k) = sum_{j=0}^{N-1} x(j)*exp(-2ijk*pi/N)             c
c                                                             c
c  The backward transform computes                            c
c     x(j) = (1/N) * sum_{k=0}^{N-1} X(k)*exp(2ijk*pi/N)      c
c                                                             c
c                                                             c
c  Requires standard FORTRAN functions - sin, cos             c
c                                                             c
c  Steve Kifowit, 9 July 1997                                 c
c                                                             c
C-------------------------------------------------------------C
C  A Duhamel-Hollman Split-Radix DIF FFT                      C
C  Reference:  Electronics Letters, January 5, 1984           C
C  Complex input and output in data arrays X and Y            C
C  Length is N = 2**M                                         C
C                                                             C
C  C.S. Burrus          Rice University         Dec 1984      C
C-------------------------------------------------------------C
c
      SUBROUTINE SFFTEU( X, Y, N, M, ITYPE )

c       INTEGER*8  N, M, ITYPE
c      REAL*16  X(*), Y(*)
c      INTEGER*8  I, J, K, N1, N2, N4, IS, ID, I0, I1, I2, I3
c      REAL*16  TWOPI, E, A, A3, CC1, SS1, CC3, SS3
c      REAL*16  R1, R2, S1, S2, S3, XT

      INTEGER  N, M, ITYPE
      REAL X(*), Y(*)
      INTEGER  I, J, K, N1, N2, N4, IS, ID, I0, I1, I2, I3
      REAL  TWOPI, E, A, A3, CC1, SS1, CC3, SS3
      REAL  R1, R2, S1, S2, S3, XT
      INTRINSIC  SIN, COS
      PARAMETER  ( TWOPI = 6.28318530717958647690)
      IF ( N .EQ. 1 ) RETURN
c
      IF ( ITYPE .EQ. -1 ) THEN
	 DO 1, I = 1, N
	    Y(I) = - Y(I)
 1       CONTINUE
      ENDIF
c
      N2 = 2 * N
      DO 10, K = 1, M-1
	 N2 = N2 / 2
	 N4 = N2 / 4
	 E = TWOPI / N2
	 A = 0.0
	 DO 20, J = 1, N4
	    A3 = 3 * A
	    CC1 = COS( A )
	    SS1 = SIN( A )
	    CC3 = COS( A3 )
	    SS3 = SIN( A3 )
	    A = J * E
	    IS = J
	    ID = 2 * N2
 40         DO 30, I0 = IS, N-1, ID
	       I1 = I0 + N4
	       I2 = I1 + N4
	       I3 = I2 + N4
	       R1 = X(I0) - X(I2)
	       X(I0) = X(I0) + X(I2)
	       R2 = X(I1) - X(I3)
	       X(I1) = X(I1) + X(I3)
	       S1 = Y(I0) - Y(I2)
	       Y(I0) = Y(I0) + Y(I2)
	       S2 = Y(I1) - Y(I3)
	       Y(I1) = Y(I1) + Y(I3)
	       S3 = R1 - S2
	       R1 = R1 + S2
	       S2 = R2 - S1
	       R2 = R2 + S1
	       X(I2) = R1 * CC1 - S2 * SS1
	       Y(I2) = - S2 * CC1 - R1 * SS1
	       X(I3) = S3 * CC3 + R2 * SS3
	       Y(I3) = R2 * CC3 - S3 * SS3
 30         CONTINUE
	    IS = 2 * ID - N2 + J
	    ID = 4 * ID
	    IF ( IS .LT. N ) GOTO 40
 20      CONTINUE
 10   CONTINUE
c
C--------LAST STAGE, LENGTH-2 BUTTERFLY ----------------------C
c
      IS = 1
      ID = 4
 50   DO 60, I0 = IS, N, ID
	 I1 = I0 + 1
	 R1 = X(I0)
	 X(I0) = R1 + X(I1)
	 X(I1) = R1 - X(I1)
	 R1 = Y(I0)
	 Y(I0) = R1 + Y(I1)
	 Y(I1) = R1 - Y(I1)
 60   CONTINUE
      IS = 2 * ID - 1
      ID = 4 * ID
      IF ( IS .LT. N ) GOTO 50
c
C-------BIT REVERSE COUNTER-----------------------------------C
c
 100  J = 1
      N1 = N - 1
      DO 104, I = 1, N1
	 IF ( I .GE. J ) GOTO 101
	 XT = X(J)
	 X(J) = X(I)
	 X(I) = XT
	 XT = Y(J)
	 Y(J) = Y(I)
	 Y(I) = XT
 101     K = N / 2
 102     IF ( K .GE. J ) GOTO 103
	 J = J - K
	 K = K / 2
	 GOTO 102
 103     J = J + K
 104  CONTINUE
c
      IF ( ITYPE .EQ. -1 ) THEN
	 DO 2, I = 1, N
	    X(I) = X(I) / N
	    Y(I) = - Y(I) / N
 2       CONTINUE
      ENDIF
c
      RETURN
c
c ... End of subroutine SFFTEU ...
c
      END

