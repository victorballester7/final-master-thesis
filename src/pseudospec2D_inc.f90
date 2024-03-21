!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines for computing spatial derivatives and nonlinear
! terms in incompressible HD, MHD and Hall-MHD equations in 2D
! using a pseudo-spectral method. You should use the FFTPLANS
! and MPIVARS modules (see the file 'fftp2D_mod.f90') in each
! program that call any of the subroutines in this file.
!
! NOTATION: index 'i' is 'x'
!           index 'j' is 'y'
!
! 2003 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

!*****************************************************************
SUBROUTINE derivk2(a,b,dir)
!-----------------------------------------------------------------
!
! Two-dimensional derivative of the matrix 'a'
!
! Parameters
!     a  : input matrix
!     b  : at the output contains the derivative da/dk_dir
!     dir: =1 derivative in the x-direction
!          =2 derivative in the y-direction
!
   USE ali
   USE kes
   USE var
   USE grid
!      USE mpivars
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: a,b
   INTEGER :: dir
   INTEGER :: i,j

!
! Derivative in the x-direction
!
   IF (dir.eq.1) THEN
      DO j = 1,n
         DO i = 1,n/2+1
            IF ((ka2(i,j).le.kmax).and.(ka2(i,j).ge.tiny)) THEN
               b(i,j) = im*ka(i)*a(i,j)
            ELSE
               b(i,j) = 0.0d0
            ENDIF
         END DO
      END DO
!
! Derivative in the y-direction
!
   ELSE
      DO j = 1,n
         DO i = 1,n/2+1
            IF ((ka2(i,j).le.kmax).and.(ka2(i,j).ge.tiny)) THEN
               b(i,j) = im*ka(j)*a(i,j)
            ELSE
               b(i,j) = 0.0d0
            ENDIF
         END DO
      END DO
   ENDIF

   RETURN
END SUBROUTINE derivk2

!*****************************************************************
SUBROUTINE laplak2(a,b)
!-----------------------------------------------------------------
!
! Two-dimensional Laplacian of the matrix 'a'
!
! Parameters
!     a: input matrix
!     b: at the output contains the Laplacian d2a/dka2
!
   USE kes
   USE grid
   USE ali
!      USE mpivars
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: a,b
   INTEGER :: i,j

   DO j = 1,n
      DO i = 1,n/2+1
         IF ((ka2(i,j).le.kmax).and.(ka2(i,j).ge.tiny)) THEN
            b(i,j) = -ka2(i,j)*a(i,j)
         ELSE
            b(i,j) = 0.0d0
         ENDIF
      END DO
   END DO

   RETURN
END SUBROUTINE laplak2

!*****************************************************************
SUBROUTINE poisson(a,b,c)
!-----------------------------------------------------------------
!
! Poisson bracket of the scalar fields A and B
! in real space.
!
! Parameters
!     a: input matrix
!     b: input matrix
!     c: Poisson bracket {A,B} [output]
!
!      USE mpivars
   USE kes
   USE ali
   USE grid
   USE fft
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: a,b,c
   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: c1,c2
   DOUBLE PRECISION, DIMENSION(n,n)    :: r1,r2,r3
   INTEGER :: i,j

!
! Computes dA/dx.dB/dy
!
   CALL derivk2(a,c1,1)
   CALL derivk2(b,c2,2)
!     CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!     CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
   CALL rfftwnd_f77_one_complex_to_real(plancr, c1, r1)
   CALL rfftwnd_f77_one_complex_to_real(plancr, c2, r2)
   DO j = 1,n
      DO i = 1,n
         r3(i,j) = r1(i,j)*r2(i,j)
      END DO
   END DO

!
! Computes dA/dy.dB/dx
!
   CALL derivk2(a,c1,2)
   CALL derivk2(b,c2,1)
!     CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
!     CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
   CALL rfftwnd_f77_one_complex_to_real(plancr, c1, r1)
   CALL rfftwnd_f77_one_complex_to_real(plancr, c2, r2)
   DO j = 1,n
      DO i = 1,n
         r3(i,j) = (r3(i,j)-r1(i,j)*r2(i,j))/dble(n)**4
      END DO
   END DO
!     CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)
   CALL rfftwnd_f77_one_real_to_complex(planrc, r3, c)
   DO j = 1,n
      DO i = 1,n/2+1
         IF ((ka2(i,j).ge.kmax).and.(ka2(i,j).le.tiny)) THEN
            c(i,j) = 0.0d0
         ENDIF
      END DO
   END DO
   RETURN
END SUBROUTINE poisson

!*****************************************************************
SUBROUTINE energy(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the mean kinetic or magnetic energy in 2D,
! and the mean square current density or vorticity.
! The output is valid only in the first node.
!
! Parameters
!     a  : input matrix with the scalar field
!     b  : at the output contains the energy
!     kin: =2 computes the square of the scalar field
!          =1 computes the energy
!          =0 computes the current or vorticity
!
   USE kes
   USE grid
   USE ali
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: a
   DOUBLE PRECISION  :: b
   DOUBLE PRECISION  :: bloc
   DOUBLE PRECISION  :: tmp,dbl
   INTEGER :: kin
   INTEGER :: i,j

   bloc = 0.0d0
   dbl  = 2.0d0

   tmp = 1./dble(n)**4

!
! Computes the square of the scalar field
!
! Computes the energy
!
! Computes the current or vorticity
!
   DO j = 1,n
      DO i = 1,n/2+1
         dbl  = 2.0d0
         IF (i.eq.1) dbl=1.0d0
         IF ((ka2(i,j).le.kmax).and.(ka2(i,j).ge.tiny)) THEN
            bloc = bloc+dbl*ka2(i,j)**kin*abs(a(i,j))**2*tmp
         ENDIF
      END DO
   END DO
   b=bloc
!
! Computes the reduction between nodes
!
!      CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                      MPI_COMM_WORLD,ierr)

   RETURN
END SUBROUTINE energy

!###############################
!*****************************************************************
SUBROUTINE inerprod(a,b,kin,rslt)
!-----------------------------------------------------------------
! Parameters
!     a  : first  input matrix
!     b  : second input matrix
!     kin: = multiplies by the laplacian to this power

   USE kes
   USE grid
!      USE mpivars
   USE ali

   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: a,b
   DOUBLE PRECISION  :: tmp,tmq,dbl
   DOUBLE PRECISION  :: rslt
   INTEGER :: kin
   INTEGER :: i,j

   tmp = 0.0d0
   tmq = 1./dble(n)**4

   DO j = 1,n
      DO i = 1,n/2+1
         dbl  = 2.0d0
         IF (i.eq.1) dbl=1.0d0
         IF ((ka2(i,j).le.kmax).and.(ka2(i,j).ge.tiny)) THEN
            tmp = tmp + dbl*ka2(i,j)**kin*dble(b(i,j)*conjg(a(i,j))) *tmq
         ENDIF
      END DO
   END DO
   rslt=tmp
!      CALL MPI_REDUCE(tmp,rslt,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                      MPI_COMM_WORLD,ierr)

   RETURN
END SUBROUTINE inerprod
!################################


!*****************************************************************
SUBROUTINE hdcheck(a,b,t,eng,ens,dir)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in HD 2D
!
! Parameters
!     a  : streamfunction
!     b  : external force
!     t  : time
!
   USE kes
   USE grid
!      USE mpivars

   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: a,b
   DOUBLE PRECISION    :: eng,ens,feng,fens,dens
   DOUBLE PRECISION    :: t
   CHARACTER*100 :: dir

!
! Computes the mean energy and enstrophy
!
   CALL energy(a, eng,1)
   CALL energy(a, ens,2)
   CALL energy(a,dens,3)

!
! Computes the energy injection rate
!
   CALL inerprod(b,a,1,feng)
   CALL inerprod(b,a,2,fens)
!
! Creates external files to store the results
!
!      IF (myrank.eq.0) THEN
   OPEN(1,file=trim(dir)//'energy_bal.txt',position='append')
   WRITE(1,20) t,eng,ens,feng
20 FORMAT( D22.14,D22.14,D22.14,D22.14 )
   CLOSE(1)
   OPEN(1,file=trim(dir)//'enstrophy_bal.txt',position='append')
   WRITE(1,21) t,ens,dens,fens
21 FORMAT( D22.14,D22.14,D22.14,D22.14 )
   CLOSE(1)
!      ENDIF
   RETURN
END SUBROUTINE hdcheck

!*****************************************************************
SUBROUTINE mhdcheck(a,b,c,d,t,eng,dir)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in MHD 2D
!
! Parameters
!     a  : streamfunction
!     b  : vector potential
!     c  : external kinetic force
!     d  : external magnetic force
!     t  : time
!
   USE kes
   USE grid
!      USE mpivars

   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: a,b,c,d
   DOUBLE PRECISION  :: eng
   DOUBLE PRECISION  :: enk, denk, fenk, henk
   DOUBLE PRECISION  :: enm, denm, fenm, henm
   DOUBLE PRECISION  :: ens, dens, fens, hens
   DOUBLE PRECISION  :: crs, dcrs, fcrs ,ecrs, hcrs
   DOUBLE PRECISION  :: asq, dasq, fasq, hasq
   DOUBLE PRECISION  :: t
   CHARACTER*100 :: dir

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ENERGY BALANCE
   CALL energy(a, enk, 1)
   CALL energy(a,denk, 2)
   CALL energy(a,henk,-1)
   CALL energy(b, enm, 1)
   CALL energy(b,denm, 2)
   CALL energy(b,henm,-1)
   CALL inerprod(c,a,1,fenk)
   CALL inerprod(d,b,1,fenm)
   eng = enk+enm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     ENSTROPHY BALANCE
   ens=denk
   CALL energy(a,dens,3)
   CALL inerprod(c,a,2,fens)
   CALL energy(a,hens,0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     CROSS HELICITY BALANCE
   CALL inerprod(b,a, 1, crs)
   CALL inerprod(b,a, 2,dcrs)
   CALL inerprod(b,a,-1,hcrs)
   CALL inerprod(a,d, 1,fcrs)
   CALL inerprod(b,c, 1,ecrs)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SQUARE VECTOR POTENTIAL
   CALL energy(b, asq,0)
   dasq=enm
   CALL inerprod(b,d,0,fasq)
   CALL energy(b,hasq,-2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Creates external files to store the results
!
!      IF (myrank.eq.0) THEN
   OPEN(1,file=trim(dir)//'energy_bal.txt',position='append')
   WRITE(1,10) t,eng,enk,enm,denk,denm,fenk,fenm,henk,henm
10 FORMAT(E22.14,E22.14,E22.14,E22.14,E22.14,E22.14,E22.14,E22.14,&
      E22.14,E22.14)
   CLOSE(1)
   OPEN(1,file=trim(dir)//'cross_bal.txt',position='append')
   WRITE(1,11) t,crs,dcrs,fcrs,ecrs,hcrs
11 FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14,E22.14)
   CLOSE(1)
   OPEN(1,file=trim(dir)//'enstrophy_bal.txt',position='append')
   WRITE(1,12) t,ens,dens,fens,hens
12 FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14)
   CLOSE(1)
   OPEN(1,file=trim(dir)//'sqr_vecpot_bal.txt',position='append')
   WRITE(1,13) t,asq,dasq,fasq,hasq
13 FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14)
   CLOSE(1)
!      ENDIF

   RETURN
END SUBROUTINE mhdcheck

!*****************************************************************
SUBROUTINE spectrum(a,ext,kin,dir)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D.
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction or vector potential
!     ext: the extension used when writting the file
!     kin: =1 computes the kinetic spectrum
!          =0 computes the magnetic spectrum
!
   USE kes
   USE grid
!      USE mpivars
   IMPLICIT NONE

   DOUBLE PRECISION, DIMENSION(n/2+1)        :: Ek
   DOUBLE COMPLEX, DIMENSION(n/2+1,n)          :: a
   DOUBLE PRECISION        :: tmp,dbl
   INTEGER     :: kin
   INTEGER     :: kmn
   INTEGER     :: i,j
   CHARACTER*3 :: ext
   CHARACTER*100 :: dir

!
! Sets Ek to zero
!
   DO i = 1,n/2+1
      Ek(i) = 0.0d0
   END DO
!
! Computes the energy spectrum
!
   tmp = 1.0d0/dble(n)**4
   DO j = 1,n
      DO i = 1,n/2+1
         dbl = 2.0d0
         IF (i.eq.1) dbl=1.0d0
         kmn = int(sqrt(ka2(i,j))+.5d0)
         IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
            Ek(kmn) = Ek(kmn)+dbl*ka2(i,j)*abs(a(i,j))**2*tmp
         ENDIF
      END DO
   END DO
!
! Computes the reduction between nodes
! and exports the result to a file
!
!      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                      MPI_COMM_WORLD,ierr)
!      IF (myrank.eq.0) THEN
   IF (kin.eq.1) THEN
      OPEN(1,file=trim(dir)//'kspectrum/kspectrum.' // ext // '.txt')
   ELSE
      OPEN(1,file=trim(dir)//'mspectrum/mspectrum.' // ext // '.txt')
   ENDIF
   WRITE(1,20) Ek
20 FORMAT( E23.15 )
   CLOSE(1)
!      ENDIF

   RETURN
END SUBROUTINE spectrum

!*****************************************************************
SUBROUTINE vectrans(a,b,c,ext1,ext2,dir)
!-----------------------------------------------------------------
!
! Computes the square vector potential transfer in
! Fourier space in 2D MHD. The output is written
! to a file by the first node.
!
! Parameters
!     a  : streamfunction
!     b  : vector potential
!     ext: the extension used when writting the file
!
   USE kes
   USE grid
!      USE mpivars
   IMPLICIT NONE

   DOUBLE PRECISION, DIMENSION(n/2+1)        :: Ek
   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: a,b,c,d
   DOUBLE PRECISION        :: tmp
   INTEGER     :: kmn,dbl
   INTEGER     :: i,j
   CHARACTER*3 :: ext1,ext2
   CHARACTER*100 :: dir

!
! Sets Ek to zero
!
   DO i = 1,n/2+1
      Ek(i) = 0.
   END DO
!
! Computes the square vector potential flux
!
   CALL poisson(b,c,d)
   tmp = 1.0d0/dble(n)**4
   DO j = 1,n
      DO i = 1,n/2+1
         dbl=2
         IF (i.eq.1) dbl=1
         kmn = int(sqrt(ka2(i,j))+.5)
         IF (kmn.le.n/2+1) THEN
            Ek(kmn) = Ek(kmn)+dbl*dble(a(i,j)*conjg(d(i,j)))*tmp
         ENDIF
      END DO
   END DO
!
! Computes the reduction between nodes
! and exports the result to a file
!
!      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                      MPI_COMM_WORLD,ierr)
!      IF (myrank.eq.0) THEN
   OPEN(1,file=trim(dir)//'vectrans/vectrans_' // ext1 // '.' // ext2 // '.txt')
   WRITE(1,20) Ek
20 FORMAT( E23.15 )
   CLOSE(1)
!      ENDIF

   RETURN
END SUBROUTINE vectrans

!*****************************************************************
SUBROUTINE Eprof(a,ext,dir)
!-----------------------------------------------------------------
!
! Computes the energy profile in circles of radius 1, 2, 3, ...
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction or vector potential
!     ext: the extension used when writting the file
!
   USE kes
   USE grid
   USE fft
   USE ali
!      USE mpivars
   IMPLICIT NONE

   DOUBLE PRECISION, DIMENSION(n/2+1)        :: E_R
   DOUBLE COMPLEX, DIMENSION(n/2+1,n)          :: a
   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: c1
   DOUBLE PRECISION, DIMENSION(n,n)    :: r1
   INTEGER     :: r
   INTEGER, DIMENSION(n/2+1) :: E_Num ! Number of points in each circle
   INTEGER     :: i,j
   CHARACTER*3 :: ext
   CHARACTER*100 :: dir

!
! Sets Ek to zero
!
   DO i = 1,n/2+1
      E_R(i) = 0.0d0
      E_Num(i) = 0.0d0
   END DO
!
! Computes the contribution for u_x = partial_y psi
!
   CALL derivk2(a,c1,2)
   CALL rfftwnd_f77_one_complex_to_real(plancr, c1, r1)
   DO i = 1,n
      DO j = 1,n
         ! if (i,j) is the center of the circle, skip it
         IF (i.eq.(n+1)/2 .AND. j.eq.(n+1)/2) CYCLE
         r = int(sqrt(real((i-(n+1)/2)**2+(j-(n+1)/2)**2)))
         E_R(r) = E_R(r) + r1(i,j)**2
         E_Num(r) = E_Num(r) + 1
      END DO
   END DO

! Computes the contribution for u_y = -partial_x psi
   CALL derivk2(a,c1,1)
   CALL rfftwnd_f77_one_complex_to_real(plancr, c1, r1)
   DO i = 1,n
      DO j = 1,n
         ! if (i,j) is the center of the circle, skip it
         IF (i.eq.(n+1)/2 .AND. j.eq.(n+1)/2) CYCLE
         r = int(sqrt(real((i-(n+1)/2)**2+(j-(n+1)/2)**2)))
         E_R(r) = E_R(r) + r1(i,j)**2
      END DO
   END DO

   DO i = 1,n/2+1
      IF (E_Num(i).gt.0) THEN
         E_R(i) = E_R(i)/dble(E_Num(i))
      ENDIF
   END DO

   OPEN(1,file=trim(dir)//'Eprofile/Eprofile.' // ext // '.txt')
   WRITE(1,20) E_R
20 FORMAT( E23.15 )
   CLOSE(1)

   RETURN
END SUBROUTINE Eprof
!*****************************************************************
