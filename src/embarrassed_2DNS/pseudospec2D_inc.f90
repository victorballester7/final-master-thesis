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
            IF ((ka2(i,j).le.kmax2).and.(ka2(i,j).ge.tiny)) THEN
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
            IF ((ka2(i,j).le.kmax2).and.(ka2(i,j).ge.tiny)) THEN
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
         IF ((ka2(i,j).le.kmax2).and.(ka2(i,j).ge.tiny)) THEN
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
         IF ((ka2(i,j).ge.kmax2).and.(ka2(i,j).le.tiny)) THEN
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

   tmp = 1./dble(n)**4 ! Normalization factor n^2 for FFT and we have to square it because the field is raised to the power of 2

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
         IF ((ka2(i,j).le.kmax2).and.(ka2(i,j).ge.tiny)) THEN
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
         IF ((ka2(i,j).le.kmax2).and.(ka2(i,j).ge.tiny)) THEN
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
SUBROUTINE hdcheck(a,b,t,inu,nu,imu,mu,eng,ens,node,dir)
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
   DOUBLE PRECISION    :: eng,deng,heng,feng
   DOUBLE PRECISION    :: ens,dens,hens,fens
   DOUBLE PRECISION    :: t,nu,mu
   INTEGER :: inu,imu
   CHARACTER*3         :: node
   CHARACTER*100       :: dir


!
! Computes the mean energies and enstrophies for different terms (kinetic, dissipation, hypo-Dissipation)
!
   CALL energy(a, eng,1    ) ! Kinetic energy
   CALL energy(a,deng,1+inu) ! Dissipation Energy
   CALL energy(a,heng,1-imu) ! Hypo-dissipation Energy

   CALL energy(a, ens,2    ) ! Enstrophy
   CALL energy(a,dens,2+inu) ! Dissipation Enstrophy
   CALL energy(a,hens,2-imu) ! Hypo-dissipation Enstrophy

!
! Computes the energy injection rate
!
   CALL inerprod(b,a,1,feng)
   CALL inerprod(b,a,2,fens)
!
! Creates external files to store the results
!
!      IF (myrank.eq.0) THEN
   OPEN(1,file=trim(dir)//'/energy_bal.'//node//'.txt',position='append')
   WRITE(1,20) t,eng,nu*deng,mu*heng,feng
20 FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14 )
   CLOSE(1)
   OPEN(1,file=trim(dir)//'/enstrophy_bal.'//node//'.txt',position='append')
   WRITE(1,21) t,ens,nu*dens,mu*hens,fens
21 FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14 )
   CLOSE(1)
!      ENDIF
   RETURN
END SUBROUTINE hdcheck


!*****************************************************************
SUBROUTINE spectrum(a,ext,kin,node,dir)
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
   CHARACTER*3 :: node
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
      OPEN(1,file=trim(dir)//'/kspectrum/kspectrum.'//node//'.' // ext // '.txt')
   ELSE
      OPEN(1,file=trim(dir)//'/mspectrum/mspectrum.'//node//'.' // ext // '.txt')
   ENDIF
   WRITE(1,20) Ek
20 FORMAT( E23.15 )
   CLOSE(1)
!      ENDIF

   RETURN
END SUBROUTINE spectrum

!*****************************************************************
SUBROUTINE vectrans(a,b,c,ext1,ext2,node,dir)
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
   CHARACTER*3 :: ext1,ext2,node
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
   OPEN(1,file=trim(dir)//'/vectrans/vectrans_'// ext1 //'.'//node//'.' // ext2 // '.txt')
   WRITE(1,20) Ek
20 FORMAT( E23.15 )
   CLOSE(1)
!      ENDIF

   RETURN
END SUBROUTINE vectrans


!*****************************************************************
SUBROUTINE shift(a,b,dx,dy)
!-----------------------------------------------------------------
!
! Two-dimensional displacement
!
! Parameters
!     a  : input matrix
!     b  : output
!
   USE ali
   USE kes
   USE var
   USE grid
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: a,b
   DOUBLE PRECISION        :: dx,dy
   INTEGER :: i,j

   DO i = 1,n/2+1
      DO j = 1,n
         b(i,j) = a(i,j) *exp(-im*(ka(i)*dx+ka(j)*dy))
      END DO
   END DO

   RETURN
END SUBROUTINE shift


!*****************************************************************
SUBROUTINE initialcond(a,seed)
!-----------------------------------------------------------------
! a is the streamfunction
! seed is the seed for the random number generator
   USE var
   USE kes
   USE ali
   USE grid
   USE random
   USE fft
   IMPLICIT NONE
   DOUBLE COMPLEX, DIMENSION(n/2+1,n)     :: a
   DOUBLE PRECISION, DIMENSION(n,n)   :: R1
   INTEGER :: i,j,jj,ki
   INTEGER :: seed
   DOUBLE PRECISION :: phase1,phase2,phase3,phase4

   DO j = 1,n
      DO i = 1,n
         R1(i,j) = 0.0d0
      END DO
   END DO
   DO ki=1,5
      DO j = 1,n
         DO i = 1,n
            phase1=randu(seed)
            phase2=randu(seed)
            phase3=randu(seed)
            phase4=randu(seed)
            R1(i,j) = R1(i,j)                                &
               + 1.0d0*sin(2.0d0*ki *pi*(dble(i)-1)/dble(n)+phase1) &
               + 1.0d0*sin(2.0d0*ki *pi*(dble(j)-1)/dble(n)+phase2) &
               + 1.0d0*sin(2.0d0*ki *pi*(dble(i)-1)/dble(n)+phase3) &
               * 1.0d0*sin(2.0d0*ki *pi*(dble(j)-1)/dble(n)+phase4)
         END DO
      END DO
   ENDDO
   CALL rfftwnd_f77_one_real_to_complex(planrc, R1, a)
   RETURN

END SUBROUTINE initialcond
!-----------------------------------------------------------------

!*****************************************************************
SUBROUTINE forcing(iflow,f0,kup,kdn,seed,fk)
!-----------------------------------------------------------------

   USE var
   USE kes
   USE ali
   USE grid
   USE random
   USE fft
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n)     :: fk,c1,c2
   DOUBLE PRECISION, DIMENSION(n,n)   :: r1
   INTEGER :: i,j,jj,iflow,ikup
   INTEGER :: seed
   DOUBLE PRECISION        :: f0,kup,kdn,zz,zx,zy,kf
   DOUBLE PRECISION        :: tmp,tmp1,tmp2
   DOUBLE PRECISION        :: phase1,phase2,dkup

   kf=0.5d0*dble(kup+kdn)
   !if (myrank.eq.0) print*,"DBG pseudo: * iflow=",iflow,f0,kf
!!!!!!!  STEADY FORCING  !!!!!!!!!!!
   IF (iflow.eq.0) THEN
      DO j = 1,n
         DO i = 1,n
            r1(i,j) = sin(2*kup*pi*(dble(i)-1)/dble(n)) &
               * sin(2*kdn*pi*(dble(j)-1)/dble(n))
         END DO
      END DO
      CALL rfftwnd_f77_one_real_to_complex(planrc,r1,fk)
      CALL energy(fk,tmp1,1)
      tmp=f0/sqrt(tmp1)! ==> sum k^2 fk^2=f0^2
      DO i = 1,n/2+1
         DO j = 1,n
            IF ((ka2(i,j).le.kmax2).and.(ka2(i,j).ge.tiny)) THEN
               fk(i,j) = tmp*fk(i,j)
            ELSE
               fk(i,j) = 0.0d0
            ENDIF
         END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! (4) DISK of FIRE !!!!
   ELSE IF (iflow.eq.3) THEN
      !if (myrank.eq.0) print*,"DBG pseudo: iflow=",iflow,f0,kf
      ikup = 2*int(kf)+1
      dkup = dble(kf)/2
      DO j = 1,n
         DO i = 1,n
            zz=(2*pi*(dble(i-n/2)-1)/dble(n))**2
            zz=(2*pi*(dble(j-n/2)-1)/dble(n))**2+zz
            zz=0.5d0*kup*kup*zz
            if (zz.gt.80.d0) zz=80.0d0
            r1(i,j) = f0*exp(-zz)
         END DO
      END DO
      CALL rfftwnd_f77_one_real_to_complex(planrc,r1,c1)
      DO i = 1,n/2+1
         DO j = 1,n
            fk(i,j) = 0.0
         END DO
      END DO
      DO jj=1,10
         phase1 = 2*pi*randu(seed)
         tmp1= pi/kdn/2*sin(phase1)
         tmp2= pi/kdn/2*cos(phase1)
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = 1,n/2+1
            DO j = 1,n
               fk(i,j) = fk(i,j) +c2(i,j)
            END DO
         END DO
         tmp1= -tmp1
         tmp2= -tmp2
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = 1,n/2+1
            DO j = 1,n
               fk(i,j) = fk(i,j) -c2(i,j)
            END DO
         END DO
      END DO  !! jj
      DO i = 1,n/2+1
         DO j = 1,n
            IF (ka2(i,j).ge.0.1) THEN
               fk(i,j) = fk(i,j)/ka2(i,j)
            ELSE
               fk(i,j) = 0.0d0
            ENDIF
         END DO
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! (4) RING of FIRE !!!!
   ELSE IF (iflow.eq.4) THEN
      !if (myrank.eq.0) print*,"DBG pseudo: iflow=",iflow,f0,kf
      ikup = 2*int(kf)+1
      dkup = dble(kf)/2
      DO j = 1,n
         DO i = 1,n
            zz=(2*pi*(dble(i-n/2)-1)/dble(n))**2
            zz=(2*pi*(dble(j-n/2)-1)/dble(n))**2+zz
            zz=0.5d0*kup*kup*zz
            if (zz.gt.80.d0) zz=80.0d0
            r1(i,j) = f0*exp(-zz)
         END DO
      END DO
      CALL rfftwnd_f77_one_real_to_complex(planrc,r1,c1)
      DO i = 1,n/2+1
         DO j = 1,n
            fk(i,j) = 0.0
         END DO
      END DO
      DO jj=1,10
         phase1 = 2*pi*randu(seed)
         tmp1= pi/kdn/2*sin(phase1)
         tmp2= pi/kdn/2*cos(phase1)
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = 1,n/2+1
            DO j = 1,n
               fk(i,j) = fk(i,j) +c2(i,j)
            END DO
         END DO
         tmp1= -tmp1
         tmp2= -tmp2
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = 1,n/2+1
            DO j = 1,n
               fk(i,j) = fk(i,j) -c2(i,j)
            END DO
         END DO
      END DO  !! jj
      DO i = 1,n/2+1
         DO j = 1,n
            IF (ka2(i,j).ge.0.1) THEN
               fk(i,j) = fk(i,j)/ka2(i,j)
            ELSE
               fk(i,j) = 0.0d0
            ENDIF
         END DO
      END DO
   ENDIF
   RETURN
END SUBROUTINE forcing


!*****************************************************************
SUBROUTINE EnergyEnstropy_profiles(a,p,ext,node,dir)
!-----------------------------------------------------------------
!
! Computes the energy profile in circles of radius 1, 2, 3, ...
! The output is written to a file containing only one column, corresponding to the energy at each radius.
!
! Parameters
!     a  : streamfunction or vector potential
!     w  : vorticity
!     p  : p-th moment of the field (p=2 for the usual energy and enstrophy)
!     ext: the extension used when writting the file
!
   USE kes
   USE grid
   USE fft
   USE ali
!      USE mpivars
   IMPLICIT NONE

   INTEGER    :: r_max


   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)        :: E_R, W_R
   DOUBLE COMPLEX, DIMENSION(n/2+1,n)          :: a
   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: c1, c2
   DOUBLE PRECISION, DIMENSION(n,n)    :: r1,r2, r3
   INTEGER     :: r, p
   INTEGER, ALLOCATABLE, DIMENSION(:) :: Num  ! Number of points in each circle
   INTEGER     :: i,j
   CHARACTER*3 :: ext,node
   CHARACTER*100 :: dir
   CHARACTER*3 :: p_str
   DOUBLE PRECISION :: tmp

   WRITE(p_str,'(i0)') p

   r_max = int(n/sqrt(2.0)) ! The maximum radius is the diagonal of the domain (I computed it, there's no need to add 1 to be conservative)

   ALLOCATE(E_R(r_max),W_R(r_max),Num(r_max))

!
! Sets Ek to zero
!
   DO i = 1,r_max
      E_R(i) = 0.0d0
      W_R(i) = 0.0d0
      Num(i) = 0.0d0
   END DO



! !
! ! Computes the contribution for u_x = partial_y psi
! !
   ! normalization of c1 for FFT
   CALL derivk2(a,c1,2)
   CALL rfftwnd_f77_one_complex_to_real(plancr, c1, r1)

   ! Computes the contribution for u_y = -partial_x psi
   CALL derivk2(a,c2,1)
   CALL rfftwnd_f77_one_complex_to_real(plancr, c2, r2)

   CALL laplak2(a,c1)     ! make W = vorticity
   CALL rfftwnd_f77_one_complex_to_real(plancr, c1, r3)

   DO i = 1,n
      DO j = 1,n
         r = int(sqrt(real((i-(n/2+0.5))**2+(j-(n/2+0.5))**2)))
         r = r+1 ! to avoid the zero radius
         E_R(r) = E_R(r) + (r1(i,j)**2 + r2(i,j)**2)**(p/2)
         W_R(r) = W_R(r) + r3(i,j)**p
         Num(r) = Num(r) + 1
      END DO
   END DO

   tmp=dble(n)**(2*p)
   DO i = 1,r_max
      IF (Num(i).gt.0) THEN
         E_R(i) = E_R(i)/dble(Num(i))/tmp ! n^2p = (n^2)^p, where the n^2 is the normalization factor for the FFT and the other ^p is because the field is raised to the power of p
         W_R(i) = W_R(i)/dble(Num(i))/tmp
      ENDIF
   END DO

   OPEN(1,file=trim(dir)//'/EnergyProf/Energy.' // trim(p_str)// '.'//node//'.'// ext // '.txt')
   WRITE(1,20) E_R
20 FORMAT( E23.15 )
   CLOSE(1)
   OPEN(1,file=trim(dir)//'/EnstrophyProf/Enstrophy.'// trim(p_str)//'.' //node//'.' // ext // '.txt')
   WRITE(1,24) W_R
24 FORMAT( E23.15 )
   CLOSE(1)
   RETURN
END SUBROUTINE EnergyEnstropy_profiles
!*****************************************************************

!*****************************************************************
SUBROUTINE outputfields(a,f,ext,node,dir)
!-----------------------------------------------------------------
!
! Computes the energy profile in circles of radius 1, 2, 3, ...
! The output is written to a file containing only one column, corresponding to the energy at each radius.
!
! Parameters
!     a  : streamfunction
!     p  : p-th moment of the field (p=2 for the usual energy and enstrophy)
!     ext: the extension used when writting the file
!
   USE kes
   USE grid
   USE fft
   USE ali
   !  USE mpivars
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: a,f, C1
   DOUBLE PRECISION, DIMENSION(n,n)    :: R1
   INTEGER     :: i,j
   CHARACTER*3 :: ext,node
   CHARACTER*100 :: dir

   DO i = 1,n/2+1
      DO j = 1,n
         C1(i,j) = a(i,j)/dble(n)**2 ! normalizes for the FFT
      END DO
   END DO
   CALL rfftwnd_f77_one_complex_to_real(plancr,C1,R1)
   OPEN(1,file=trim(dir) // '/output/hd2Dps.' // node // '.'// ext // '.out',form='unformatted')
   WRITE(1) R1
   CLOSE(1)
   DO i = 1,n/2+1
      DO j = 1,n
         C1(i,j) = a(i,j)*ka2(i,j)/dble(n)**2
      END DO
   END DO
   CALL rfftwnd_f77_one_complex_to_real(plancr,C1,R1)
   OPEN(1,file=trim(dir) // '/output/hd2Dww.' // node // '.' // ext // '.out',form='unformatted')
   WRITE(1) R1
   CLOSE(1)
   DO i = 1,n/2+1
      DO j = 1,n
         C1(i,j) = f(i,j)*ka2(i,j)/dble(n)**2
      END DO
   END DO
   CALL rfftwnd_f77_one_complex_to_real(plancr,C1,R1)
   OPEN(1,file=trim(dir) // '/output/hd2Dfw.' // node // '.' // ext // '.out',form='unformatted')
   WRITE(1) R1
   CLOSE(1)
   DO i = 1,n/2+1
      DO j = 1,n
         C1(i,j) = f(i,j)/dble(n)**2
      END DO
   END DO
   CALL rfftwnd_f77_one_complex_to_real(plancr,C1,R1)
   OPEN(1,file=trim(dir) // '/output/hd2Dfp.' // node // '.' // ext // '.out',form='unformatted')
   WRITE(1) R1
   CLOSE(1)
   RETURN
END SUBROUTINE outputfields
!*****************************************************************

!*****************************************************************
SUBROUTINE CFL_condition(cfl,c1,inu,nu,dt)
!-----------------------------------------------------------------

!        Parameters
!     cfl :cfl factor
!      c1 : stream fun
   !  USE mpivars
   USE kes
   USE ali
   USE grid
   USE fft
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n/2+1,n) :: c1,c3,c4
   DOUBLE PRECISION, DIMENSION(n,n)    :: r1,r2,r3
   INTEGER :: i,j,inu
   DOUBLE PRECISION        :: tmp,dt,nu,cfl
   DOUBLE PRECISION        :: tmp1,kcut,nrm

   kcut=(dble(n)/3.0d0)
   nrm=(dble(n))**2
   CALL derivk2(c1,c3,1)
   CALL derivk2(c1,c4,2)
   CALL rfftwnd_f77_one_complex_to_real(plancr,c3,r1)
   CALL rfftwnd_f77_one_complex_to_real(plancr,c4,r2)
   DO j = 1,n
      DO i = 1,n
         r3(i,j) = r1(i,j)*r1(i,j)+r2(i,j)*r2(i,j) ! u^2+v^2
      END DO
   END DO
   tmp1=maxval(r3)
   ! think of dx = 1/kcut
   ! 1st CFL condition (advection): CFL * dx/dt >= max_speed
   ! 2nd CFL condition (diffusion): nu * dt/dx^(2*inu) <= CFL  ==> CFL * dx/dt >= nu * kcut^(2*inu-1)

   ! so we take CFL * dx/dt >= max_speed + nu * kcut^(2*inu-1)
   tmp=sqrt(tmp1)/nrm+nu*kcut**(2*inu-1)
   dt = cfl/(kcut*tmp)
   !! if (myrank.eq.0) print*,"*",myrank,dt,cfl,tmp,sqrt(tmp1)/nrm,nu*kcut**(2*inu-1)

   RETURN
END SUBROUTINE CFL_condition
