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
   USE mpivars
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
   INTEGER :: dir
   INTEGER :: i,j

!
! Derivative in the x-direction
!
   IF (dir.eq.1) THEN
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.tiny)) THEN
               b(j,i) = im*ka(i)*a(j,i)
            ELSE
               b(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
!
! Derivative in the y-direction
!
   ELSE
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.tiny)) THEN
               b(j,i) = im*ka(j)*a(j,i)
            ELSE
               b(j,i) = 0.0d0
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
   USE mpivars
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
   INTEGER :: i,j

   DO i = ista,iend
      DO j = 1,n
         b(j,i) = -ka2(j,i)*a(j,i)
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
   USE mpivars
   USE kes
   USE ali
   USE grid
   USE fft
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b,c
   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: c1,c2
   DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: r1,r2,r3
   INTEGER :: i,j

!
! Computes dA/dx.dB/dy
!
   CALL derivk2(a,c1,1)
   CALL derivk2(b,c2,2)
   CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
   CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
   DO j = jsta,jend
      DO i = 1,n
         r3(i,j) = r1(i,j)*r2(i,j)
      END DO
   END DO

!
! Computes dA/dy.dB/dx
!
   CALL derivk2(a,c1,2)
   CALL derivk2(b,c2,1)
   CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
   CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
   DO j = jsta,jend
      DO i = 1,n
         r3(i,j) = (r3(i,j)-r1(i,j)*r2(i,j))/dble(n)**4
      END DO
   END DO
   CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)
   DO i = ista,iend
      DO j = 1,n
         IF ((ka2(j,i).gt.kmax2).or.(ka2(j,i).lt.tiny)) THEN
            c(j,i) = 0.0d0
         ENDIF
      END DO
   END DO
   RETURN
END SUBROUTINE poisson


!*****************************************************************
SUBROUTINE energy(a,b,kin)
!-----------------------------------------------------------------
!
! Computes the energy, vorticity or square of the scalar field
!
! Parameters
!     a  : input matrix with the scalar field
!     b  : at the output contains the energy
!     kin: =2 computes the square of the scalar field (~streamfunction)
!          =1 computes the energy
!          =0 computes the vorticity
!
   USE kes
   USE grid
   USE mpivars
   USE ali
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a
   DOUBLE PRECISION  :: b
   DOUBLE PRECISION  :: bloc
   DOUBLE PRECISION  :: tmp
   INTEGER :: kin,two
   INTEGER :: i,j

   bloc = 0.0d0
   tmp = 1.0d0/dble(n)**4 ! Normalization factor n^2 for FFT and we have to square it because the field is raised to the power of 2

   DO i = ista,iend
      two = 2
      if (i.eq.1) two=1
      DO j = 1,n
         IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.tiny)) THEN
            bloc = bloc+two* ka2(j,i)**kin *abs(a(j,i))**2*tmp
         ENDIF
      END DO
   END DO
!
! Computes the reduction between nodes
!
   CALL MPI_REDUCE(bloc,b,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(b ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   RETURN
END SUBROUTINE energy

!*****************************************************************
SUBROUTINE inerprod(a,b,kin,rslt)
!-----------------------------------------------------------------
! Parameters
!     a  : first  input matrix
!     b  : second input matrix
!     kin: = multiplies by the laplacian to this power

   USE kes
   USE grid
   USE mpivars
   USE ali

   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
   DOUBLE PRECISION  :: tmp,tmq
   DOUBLE PRECISION  :: rslt
   INTEGER :: kin,two
   INTEGER :: i,j

   tmp = 0.0d0
   tmq = 1./dble(n)**4 ! Normalization factor n^2 for FFT and we have to square it because we have two fields


   DO i = ista,iend
      two = 2
      if (i.eq.1) two=1
      DO j = 1,n
         IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.tiny)) THEN
            tmp = tmp+two*(ka2(j,i)**kin)*dble(b(j,i)*conjg(a(j,i)))*tmq
         ENDIF
      END DO
   END DO
   CALL MPI_REDUCE(tmp,rslt,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      MPI_COMM_WORLD,ierr)
   CALL MPI_BCAST(rslt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

   RETURN
END SUBROUTINE inerprod
!################################


!*****************************************************************
SUBROUTINE hdcheck(a,b,t,inu,nu,imu,mu,eng,ens,dir)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in HD 2D
!
! Parameters
!     a  : streamfunction
!     b  : external force
!     t  : time
!     inu: power of the laplacian for the viscosity
!     nu : viscosity
!     imu: power of the laplacian for the hypo-viscosity
!     mu : hypo-viscosity
!     eng: at the output contains the kinetic energy
!     ens: at the output contains the enstrophy
!    dir: directory where the results are stored
!
   USE kes
   USE grid
   USE mpivars

   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
   DOUBLE PRECISION    :: eng,deng,heng,feng
   DOUBLE PRECISION    :: ens,dens,hens,fens
   DOUBLE PRECISION    :: t,nu,mu
   DOUBLE PRECISION    :: tmq,tmp
   INTEGER :: i,j,inu,imu
   CHARACTER*100 :: dir

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
   IF (myrank.eq.0) THEN
      OPEN(1,file=trim(dir)//'/energy_bal.txt',position='append')
      WRITE(1,20) t,eng,nu*deng,mu*heng,feng
20    FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14 )
      CLOSE(1)
      OPEN(1,file=trim(dir)//'/enstrophy_bal.txt',position='append')
      WRITE(1,21) t,ens,nu*dens,mu*hens,fens
21    FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14 )
      CLOSE(1)
   ENDIF
   RETURN
END SUBROUTINE hdcheck

!*****************************************************************
SUBROUTINE spectrum(a,ext,dir)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D.
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction
!     ext: the extension used when writting the file
!     dir: directory where the results are stored

   USE kes
   USE grid
   USE mpivars
   IMPLICIT NONE

   DOUBLE PRECISION, DIMENSION(n/2+1)        :: Ek,Ektot
   DOUBLE COMPLEX, DIMENSION(n,ista:iend)    :: a,b
   DOUBLE PRECISION        :: tmp,two
   INTEGER     :: kin
   INTEGER     :: kmn
   INTEGER     :: i,j
   CHARACTER*4 :: ext
   CHARACTER*100 :: dir


   tmp = 1.0d0/dble(n)**4

   ! Computes the velocity energy spectrum
   DO i = 1,n/2+1
      Ek(i) = 0.0d0
   END DO
   DO i = ista,iend
      two=2.0d0
      IF (i.eq.1) two=1.0d0
      DO j = 1,n
         kmn = int(sqrt(ka2(j,i))+.5d0) ! we add 0.5 to prevent kmn from being 0
         IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
            Ek(kmn) = Ek(kmn)+two*ka2(j,i)*abs(a(j,i))**2*tmp
         ENDIF
      END DO
   END DO
   CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      MPI_COMM_WORLD,ierr)
   IF (myrank.eq.0) THEN
      OPEN(1,file=trim(dir)//'/kspectrum/kspectrum.' // ext // '.txt')
      WRITE(1,20) Ektot
      CLOSE(1)
   ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
20 FORMAT( E23.15 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   RETURN
END SUBROUTINE spectrum


!*****************************************************************
SUBROUTINE vectrans(a,b,c,ext1,ext2,dir)
!-----------------------------------------------------------------
!
! Computes the square vector potential transfer in
! Fourier space in 2DHD. The output is written
! to a file by the first node.
!
! Parameters
!     a  : streamfunction
!     b  : advected field
!     c  : advecting SF
!     ext: the extension used when writting the file
!
   USE kes
   USE grid
   USE mpivars
   IMPLICIT NONE

   DOUBLE PRECISION, DIMENSION(n/2+1)        :: Ek,Ektot
   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b,c,d
   DOUBLE PRECISION        :: tmp,flx
   INTEGER     :: kmn,kin,two
   INTEGER     :: i,j
   CHARACTER*3 :: ext1
   CHARACTER*4 :: ext2
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
   tmp = 1./dble(n)**4
   CALL poisson(b,c,d)
   DO i = ista,iend
      two = 2
      if (i.eq.1) two=1
      DO j = 1,n
         kmn = int(sqrt(ka2(j,i))+.5)
         IF (kmn.le.n/2+1) THEN
            Ek(kmn) = Ek(kmn)+two*dble(a(j,i)*conjg(d(j,i)))*tmp
         ENDIF
      END DO
   END DO
!
! Computes the reduction between nodes
! and exports the result to a file
!
   CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      MPI_COMM_WORLD,ierr)
   IF (myrank.eq.0) THEN
      OPEN(1,file=trim(dir)//'/vectrans/vectrans_' // ext1 // '.' // ext2 // '.txt')
      flx=0.0d0
      DO i=1,n/2+1
         WRITE(1,20) Ektot(i),flx
         flx=flx+Ektot(i)
      ENDDO
20    FORMAT( E23.15,E23.15 )
      CLOSE(1)
   ENDIF

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
!     dx : displacement in the x-direction
!     dy : displacement in the y-direction

   USE ali
   USE kes
   USE var
   USE grid
   USE mpivars
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
   DOUBLE PRECISION        :: dx,dy
   INTEGER :: i,j

   DO i = ista,iend
      DO j = 1,n
         b(j,i) = a(j,i) *exp(-im*(ka(i)*dx+ka(j)*dy))
      END DO
   END DO

   RETURN
END SUBROUTINE shift


!*****************************************************************
SUBROUTINE CFL_condition(cfl,c1,inu,nu,dt)
!-----------------------------------------------------------------

! Parameters
!     cfl : cfl factor
!      c1 : stream fun
!     inu : power of the laplacian
!      nu : viscosity
!      dt : time step

   USE mpivars
   USE kes
   USE ali
   USE grid
   USE fft
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: c1,c2,c3,c4
   DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: r1,r2,r3
   INTEGER :: i,j,inu
   DOUBLE PRECISION        :: tmp,dt,b0,mu,nu,cfl
   DOUBLE PRECISION        :: tmp1,tmp2,kcut,nrm

   kcut=(dble(n)/3.0d0)
   nrm=(dble(n))**2
   CALL derivk2(c1,c3,1)
   CALL derivk2(c1,c4,2)
   CALL fftp2d_complex_to_real(plancr,c3,r1,MPI_COMM_WORLD)
   CALL fftp2d_complex_to_real(plancr,c4,r2,MPI_COMM_WORLD)
   DO j = jsta,jend
      DO i = 1,n
         r3(i,j) = r1(i,j)*r1(i,j)+r2(i,j)*r2(i,j)
      END DO
   END DO
   tmp=maxval(r3)
   call MPI_REDUCE(tmp,tmp1,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
   call MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   ! think of dx = 1/kcut
   ! 1st CFL condition (advection): CFL * dx/dt >= max_speed
   ! 2nd CFL condition (diffusion): nu * dt/dx^(2*inu) <= CFL  ==> CFL * dx/dt >= nu * kcut^(2*inu-1)

   ! so we take CFL * dx/dt >= max_speed + nu * kcut^(2*inu-1)
   tmp=sqrt(tmp1)/nrm+nu*kcut**(2*inu-1)
   dt = cfl/(kcut*tmp)

   RETURN
END SUBROUTINE CFL_condition

!*****************************************************************
SUBROUTINE forcing(iflow,f0,kup,kdn,seed,myseed,fk)
!-----------------------------------------------------------------
! Sets the forcing term for the streamfunction
! Parameters
! iflow : type of forcing
! f0    : amplitude of the forcing
! kup   : 1st wavenumber of the forcing
! kdn   : 2nd wavenumber of the forcing
! seed: the seed for the random number generator (global)
! myseed: the seed for the random number generator (local for each core)
! fk    : where the forcing is stored
   USE mpivars
   USE var
   USE kes
   USE ali
   USE grid
   USE random
   USE fft
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n,ista:iend)     :: fk,c1,c2
   DOUBLE PRECISION, DIMENSION(n,jsta:jend)   :: r1
   INTEGER :: i,j,jj,iflow
   INTEGER :: seed,myseed
   DOUBLE PRECISION        :: f0,kup,kdn,zz
   DOUBLE PRECISION        :: tmp1,tmp2,amp,radius
   DOUBLE PRECISION        :: phase1,phase2

   !!! (4) RING of FIRE !!!!
   IF (iflow.eq.4) THEN
      DO j = jsta,jend
         DO i = 1,n
            zz=(2*pi*(dble(i-n/2)-1)/dble(n))**2
            zz=(2*pi*(dble(j-n/2)-1)/dble(n))**2+zz
            zz=0.5d0*kup*kup*zz
            if (zz.gt.80.d0) zz=80.0d0
            r1(i,j) = f0*exp(-zz)
         END DO
      END DO
      CALL fftp2d_real_to_complex(planrc,r1,c1,MPI_COMM_WORLD)
      DO i = ista,iend
         DO j = 1,n
            fk(j,i) = 0.0
         END DO
      END DO
      DO jj=1,10
         phase1 = 2*pi*randu(seed)
         radius = pi/kdn
         CALL MPI_BCAST(phase1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         tmp1= radius*sin(phase1)
         tmp2= radius*cos(phase1)
         amp = randu(seed)
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = ista,iend
            DO j = 1,n
               fk(j,i) = fk(j,i) +amp*c2(j,i)
            END DO
         END DO
         phase1 = 2*pi*randu(seed)
         tmp1= radius*sin(phase1)
         tmp2= radius*cos(phase1)
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = ista,iend
            DO j = 1,n
               fk(j,i) = fk(j,i) -amp*c2(j,i)
            END DO
         END DO
      END DO  !! jj
      DO i = ista,iend
         DO j = 1,n
            IF (ka2(j,i).ge.0.1) THEN
               fk(j,i) = fk(j,i)/ka2(j,i)
            ELSE
               fk(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
      !!! (4) DISK of FIRE !!!!
   ELSE IF (iflow.eq.5) THEN
      DO j = jsta,jend
         DO i = 1,n
            zz=(2*pi*(dble(i-n/2)-1)/dble(n))**2
            zz=(2*pi*(dble(j-n/2)-1)/dble(n))**2+zz
            zz=0.5d0*kup*kup*zz
            if (zz.gt.80.d0) zz=80.0d0
            r1(i,j) = f0*exp(-zz)
         END DO
      END DO
      CALL fftp2d_real_to_complex(planrc,r1,c1,MPI_COMM_WORLD)
      DO i = ista,iend
         DO j = 1,n
            fk(j,i) = 0.0
         END DO
      END DO
      DO jj=1,10
         phase1 = 2*pi*randu(seed)
         radius = sqrt(abs(randu(seed)))*pi/kdn ! we don't want the forcing to be uniform in r, but more concentrated at the outer of the disk.
         tmp1= radius*sin(phase1)
         tmp2= radius*cos(phase1)
         amp = randu(seed)
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = ista,iend
            DO j = 1,n
               fk(j,i) = fk(j,i) +amp*c2(j,i)
            END DO
         END DO
         phase1 = 2*pi*randu(seed)
         radius = sqrt(abs(randu(seed)))*pi/kdn ! we don't want the forcing to be uniform in r, but more concentrated at the outer of the disk.
         tmp1= radius*sin(phase1)
         tmp2= radius*cos(phase1)
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = ista,iend
            DO j = 1,n
               fk(j,i) = fk(j,i) -amp*c2(j,i)
            END DO
         END DO
      END DO  !! jj
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.tiny)) THEN
               fk(j,i) = fk(j,i)/ka2(j,i) ! we want streamfunction forcing
            ELSE
               fk(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
   ELSE IF (iflow.gt.5) THEN
      DO i = ista,iend
         DO j = 1,n
            fk(j,i) = 0.0d0
         END DO
      END DO
   ENDIF
   RETURN
END SUBROUTINE forcing

!*****************************************************************
SUBROUTINE initialcond(iflow,u0,kup,kdn,seed,myseed,a)
!-----------------------------------------------------------------
! Sets the initial condition for the streamfunction
! Parameters
! iflow : type of forcing
! u0    : amplitude of the initial condition
! kup   : 1st wavenumber of the forcing
! kdn   : 2nd wavenumber of the forcing
! seed: the seed for the random number generator (global)
! myseed: the seed for the random number generator (local for each core)
! a     : streamfunction
   USE var
   USE kes
   USE ali
   USE grid
   USE random
   USE fft
   USE mpivars
   IMPLICIT NONE
   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a
   INTEGER :: seed, myseed, iflow
   DOUBLE PRECISION :: u0,kup,kdn,tmp1,enerk

   CALL forcing(iflow,u0,kup,kdn,seed,myseed,a)
   CALL energy(a,enerk,1)
   tmp1=u0/sqrt(enerk) ! we normalize the energy injection rate, not the forcing amplitude
   CALL MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   a = tmp1*a
   RETURN

END SUBROUTINE initialcond
!-----------------------------------------------------------------

!*****************************************************************
SUBROUTINE EnergyEnstropy_profiles(a,ext,dir)
!-----------------------------------------------------------------
!
! Computes the energy profile in circles of radius 1, 2, 3, ...
! The output is written to a file containing only one column, corresponding to the energy at each radius.
!
! Parameters
!     a  : streamfunction
!     ext: the extension used when writting the file
!    dir: directory where the results are stored
   USE kes
   USE grid
   USE fft
   USE ali
   USE mpivars
   IMPLICIT NONE

   ! INTEGER, PARAMETER :: r_max = int(n/sqrt(2.0)) ! The maximum radius is the diagonal of the domain (I computed it, there's no need to add 1 to be conservative)
   INTEGER, PARAMETER :: r_max = int(n/2) ! The maximum radius is the diagonal of the domain (I computed it, there's no need to add 1 to be conservative)
   INTEGER, PARAMETER :: p_max = 8 ! if you want to change this, remeber tochange the number in the format statement in the write statement

   DOUBLE PRECISION, DIMENSION(r_max*p_max)        :: E_R, W_R, E_R_total, W_R_total
   DOUBLE COMPLEX, DIMENSION(n,ista:iend)          :: a
   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: c1, c2
   DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: r1,r2, r3
   INTEGER     :: r, p
   INTEGER, DIMENSION(r_max) :: Num, Num_total  ! Number of points in each circle
   INTEGER     :: i,j
   CHARACTER*4 :: ext
   CHARACTER*100 :: dir
   DOUBLE PRECISION :: tmp

!
! Sets Ek to zero
!
   DO r = 1,r_max
      DO p = 1,p_max
         E_R((r-1)*p_max+p) = 0.0d0
         W_R((r-1)*p_max+p) = 0.0d0
      END DO
      Num(r) = 0.0d0
   END DO

! !
! ! Computes the contribution for u_x = partial_y psi
! !
   ! normalization of c1 for FFT
   CALL derivk2(a,c1,2)
   CALL fftp2d_complex_to_real(plancr, c1, r1,MPI_COMM_WORLD)

   ! Computes the contribution for u_y = -partial_x psi
   CALL derivk2(a,c2,1)
   CALL fftp2d_complex_to_real(plancr, c2, r2,MPI_COMM_WORLD)

   CALL laplak2(a,c1)     ! make W = vorticity
   CALL fftp2d_complex_to_real(plancr, c1, r3,MPI_COMM_WORLD)

   DO j = jsta,jend
      DO i = 1,n
         r = int(sqrt(real((i-(n/2+0.5))**2+(j-(n/2+0.5))**2)))
         r = r+1 ! to avoid the zero radius
         if (r.gt.r_max) cycle
         DO p = 1,p_max
            E_R((r-1)*p_max+p) = E_R((r-1)*p_max+p) + (r1(i,j)**2 + r2(i,j)**2)**(dble(p)/2.0)
            W_R((r-1)*p_max+p) = W_R((r-1)*p_max+p) + abs(r3(i,j))**p
         END DO
         Num(r) = Num(r) + 1
      END DO
   END DO

   ! Computes the reduction between nodes
   !
   CALL MPI_REDUCE(E_R,E_R_total,p_max*r_max,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

   CALL MPI_REDUCE(W_R,W_R_total,p_max*r_max,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)

   CALL MPI_REDUCE(Num,Num_total,r_max,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)

   IF (myrank.eq.0) THEN
      DO p = 1,p_max
         tmp=dble(n)**(2*p)
         DO r = 1,r_max
            IF (Num_total(r).gt.0) THEN
               E_R_total((r-1)*p_max + p) = (E_R_total((r-1)*p_max + p)/dble(Num_total(r))/tmp)**(1.0/dble(p)) ! n^2p = (n^2)^p, where the n^2 is the normalization factor for the FFT and the other ^p is because the field is raised to the power of p
               W_R_total((r-1)*p_max + p) = (W_R_total((r-1)*p_max + p)/dble(Num_total(r))/tmp)**(1.00/dble(p))
            ENDIF
         END DO
      END DO

      OPEN(1,file=trim(dir)//'/EnergyProf/Energy.' // ext // '.txt')
      DO r = 1,r_max
         WRITE(1,20) (E_R_total((r-1)*p_max+p),p=1,p_max)
20       FORMAT(8E22.14)
      END DO
      CLOSE(1)

      OPEN(1,file=trim(dir)//'/EnstrophyProf/Enstrophy.' // ext // '.txt')
      DO r = 1,r_max
         WRITE(1,20) (W_R_total((r-1)*p_max+p),p=1,p_max)
      END DO
      CLOSE(1)
   ENDIF

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
!     a : streamfunction
!    f  : forcing term in the streamfunction
!     p : p-th moment of the field (p=2 for the usual energy and enstrophy)
!    ext: the extension used when writting the file
!   node: node name of the core
!    dir: directory where the results are stored

   USE kes
   USE grid
   USE fft
   USE ali
   USE mpivars
   IMPLICIT NONE

   DOUBLE COMPLEX, DIMENSION(n,ista:iend)          :: a,f
   DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: C1
   DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: R1
   INTEGER     :: i,j
   CHARACTER*3 :: ext,node
   CHARACTER*100 :: dir

   DO i = ista,iend
      DO j = 1,n
         C1(j,i) = a(j,i)/dble(n)**2 ! normalizes for the FFT
      END DO
   END DO
   CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
   OPEN(1,file=trim(dir) // '/output/hd2Dps.' // node // '.'// ext // '.out',form='unformatted')
   WRITE(1) R1
   CLOSE(1)
   DO i = ista,iend
      DO j = 1,n
         C1(j,i) = a(j,i)*ka2(j,i)/dble(n)**2
      END DO
   END DO
   CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
   OPEN(1,file=trim(dir) // '/output/hd2Dww.' // node // '.'// ext // '.out',form='unformatted')
   WRITE(1) R1
   CLOSE(1)
   if (ext.eq.'001') then
      DO i = ista,iend
         DO j = 1,n
            C1(j,i) = f(j,i)*ka2(j,i)/dble(n)**2
         END DO
      END DO
      CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      OPEN(1,file=trim(dir) // '/output/hd2Dfw.' // node // '.' // ext // '.out',form='unformatted')
      WRITE(1) R1
      CLOSE(1)
      DO i = ista,iend
         DO j = 1,n
            C1(j,i) = f(j,i)/dble(n)**2
         END DO
      END DO
      CALL fftp2d_complex_to_real(plancr,C1,R1,MPI_COMM_WORLD)
      OPEN(1,file=trim(dir) // '/output/hd2Dfp.' // node // '.' // ext // '.out',form='unformatted')
      WRITE(1) R1
      CLOSE(1)
   endif
   RETURN
END SUBROUTINE outputfields
!*****************************************************************




