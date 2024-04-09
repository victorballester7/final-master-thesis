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
      SUBROUTINE filterk(a,b,q)
!-----------------------------------------------------------------
!
!     a  : input matrix
!     b  : at the output
!
      USE ali
      USE kes
      USE var
      USE grid
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      INTEGER :: i,j
      DOUBLE PRECISION  :: q,tmp

         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
                  b(j,i) = a(j,i)*exp(-0.5d0*ka2(j,i)/q**2)
               ELSE
                  b(j,i) = 0.0d0
               ENDIF
            END DO
         END DO

      RETURN
      END SUBROUTINE filterk

!*****************************************************************
      SUBROUTINE volfrac(a,ext)
!-----------------------------------------------------------------
!
!     a  : input matrix
!     b  : at the output
      USE ali
      USE kes
      USE var
      USE grid
      USE mpivars
      USE fft
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(n,jsta:jend) :: r
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      INTEGER :: i,j,ip
      DOUBLE PRECISION  :: q,tmp,tmp1,tmp2
      DOUBLE PRECISION  :: n0,n1,n2,n3,n4
      DOUBLE PRECISION  :: m0,m1,m2,m3,m4
      CHARACTER*4 :: ext


      !! if (myrank.eq.0) print*,"DBG volfrac 1 ****"
      IF (myrank.eq.0) THEN
         OPEN(1,file='volume_frac.'//ext//'.txt',position='append')
      ENDIF

      DO ip=1,8
        q=2.0d0**ip
        CALL filterk(a,b,q)
        CALL fftp2d_complex_to_real(plancr,b,r,MPI_COMM_WORLD)

        tmp=0.0d0
        DO j = jsta,jend
           DO i = 1,n
           tmp = tmp +r(i,j)**2    
           END DO
        END DO
        tmp = tmp/dble(n)**2
        !! if (myrank.eq.0) print*,"DBG volfrac 2, ",ip
        CALL MPI_REDUCE(tmp,tmp1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr) 
        CALL MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        !! if (myrank.eq.0) print*,"DBG volfrac 3, ",ip
        n0=0.0D0
        n1=0.0D0
        n2=0.0D0
        n3=0.0D0
        n4=0.0D0

        DO j = jsta,jend
           DO i = 1,n
           tmp = r(i,j)**2   
           if (tmp.gt.0.5*tmp1) n0=n0+1
           if (tmp.gt.1.0*tmp1) n1=n1+1
           if (tmp.gt.2.0*tmp1) n2=n2+1
           if (tmp.gt.4.0*tmp1) n3=n3+1
           if (tmp.gt.8.0*tmp1) n4=n4+1 
           END DO
        END DO
        !! if (myrank.eq.0) print*,"DBG volfrac 4, ",ip
        CALL MPI_REDUCE(n0,m0,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(n1,m1,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(n2,m2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(n3,m3,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(n4,m4,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
        n0=m0/dble(n)**2
        n1=m1/dble(n)**2
        n2=m2/dble(n)**2
        n3=m3/dble(n)**2
        n4=m4/dble(n)**2
        !! if (myrank.eq.0) print*,"DBG volfrac 5, ",ip
        IF (myrank.eq.0) THEN
         WRITE(1,27) q,n0,n1,n2,n3,n4
        ENDIF
      END DO !ip
   27 FORMAT( 6E22.14 )
      CLOSE(1)    
      RETURN
      END SUBROUTINE volfrac


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
               IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
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
               IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
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
            IF ((ka2(j,i).gt.kmax2).or.(ka2(j,i).lt.kmin2)) THEN
               c(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
      RETURN
      END SUBROUTINE poisson

!*****************************************************************
      SUBROUTINE parseval(a,b)
!-----------------------------------------------------------------
      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: c1,c2
      DOUBLE PRECISION, DIMENSION(n,jsta:jend)    :: r1,r2
      DOUBLE PRECISION  :: tmp,tmp1,tmp2,tmp3,tmp4
      INTEGER :: i,j
      
      CALL inerprod(a,b,0,tmp)
      CALL energy(a,tmp3,0) 
      CALL energy(b,tmp4,0)
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
               c1(j,i) = a(j,i)
               c2(j,i) = b(j,i) 
            ELSE
               c1(j,i) = 0.
               c2(j,i) = 0. 
            ENDIF
         END DO
      END DO
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      tmp1 = 0.
      DO j = jsta,jend
         DO i = 1,n
            tmp1 = tmp1+r2(i,j)*r1(i,j)/dble(n)**6
         END DO
      END DO
      CALL MPI_REDUCE(tmp1,tmp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      if (myrank.eq.0) print*,"PARCEVAL1: ",tmp
      if (myrank.eq.0) print*,"PARCEVAL2: ",tmp2 
      if (myrank.eq.0) print*,"PARCEVAL3: ",tmp3
      if (myrank.eq.0) print*,"PARCEVAL4: ",tmp4

      END SUBROUTINE parseval


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
      tmp = 1.0d0/dble(n)**4

!
! Computes the square of the scalar field
!
      DO i = ista,iend
        two = 2
        if (i.eq.1) two=1
        DO j = 1,n
           IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
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
      USE mpivars
      USE ali

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION  :: tmp,tmq
      DOUBLE PRECISION  :: rslt
      INTEGER :: kin,two
      INTEGER :: i,j

      tmp = 0.0d0
      tmq = 1./dble(n)**4

      DO i = ista,iend
        two = 2
        if (i.eq.1) two=1
        DO j = 1,n
            IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
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
      SUBROUTINE hdcheck(a,b,t,inu,nu,imu,mu,eng,ens) 
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
      USE mpivars

      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION    :: eng,deng,heng,feng
      DOUBLE PRECISION    :: ens,dens,hens,fens
      DOUBLE PRECISION    :: t,nu,mu
      DOUBLE PRECISION    :: tmq,tmp
      INTEGER :: i,j,inu,imu

!
! Computes the mean energy and enstrophy
!
      CALL energy(a, eng,1    )
      CALL energy(a,deng,1+inu)
      CALL energy(a,heng,1-imu)

      CALL energy(a, ens,2    )
      CALL energy(a,dens,2+inu)
      CALL energy(a,hens,2-imu)


!
! Computes the energy injection rate
!
      CALL inerprod(b,a,1,feng)
      CALL inerprod(b,a,2,fens)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='energy_bal.txt',position='append')
         WRITE(1,20) t,eng,nu*deng,mu*heng,feng
   20    FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14 )
         CLOSE(1)
         OPEN(1,file='enstrophy_bal.txt',position='append')
         WRITE(1,21) t,ens,nu*dens,mu*hens,fens
   21    FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14 )
         CLOSE(1)
      ENDIF      
      RETURN
      END SUBROUTINE hdcheck

!*****************************************************************
      SUBROUTINE spectrum(a,ext)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction 
!     ext: the extension used when writting the file
!
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


      tmp = 1.0d0/dble(n)**4
!eudospec2D_inc.f90!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Computes the VELOCITY energy spectrum!!!!!!!!!!!!!!
      DO i = 1,n/2+1
         Ek(i) = 0.0d0
      END DO
      DO i = ista,iend
         two=2.0d0
         IF (i.eq.1) two=1.0d0
         DO j = 1,n
            kmn = int(sqrt(ka2(j,i))+.5d0)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               Ek(kmn) = Ek(kmn)+two*ka2(j,i)*abs(a(j,i))**2*tmp
            ENDIF
         END DO
      END DO
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
            OPEN(1,file='kspectrum.' // ext // '.txt')
         WRITE(1,20) Ektot
         CLOSE(1)
      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   20    FORMAT( E23.15 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      RETURN
      END SUBROUTINE spectrum


!*****************************************************************
      SUBROUTINE vectrans(a,b,c,ext1,ext2)
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
         OPEN(1,file='vectrans_' // ext1 // '.' // ext2 // '.txt')
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
      SUBROUTINE structure(a,ext)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D. 
! The output is written to a file by the first node.
!
! Parameters
!     a  : streamfunction 
!     ext: the extension used when writting the file
!

      USE mpivars
      USE kes
      USE ali
      USE grid
      USE fft
      IMPLICIT NONE
      CHARACTER*4   :: ext
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a
      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: cx,cy
      DOUBLE PRECISION, DIMENSION(n,jsta:jend) :: rx,ry
      DOUBLE PRECISION, DIMENSION(8,n/2)       :: Stp,Stptot
      DOUBLE PRECISION, DIMENSION(8,n/2)       :: Stk,Stktot
      DOUBLE PRECISION, DIMENSION(8,n/2)       :: Stpa,Stpatot
      DOUBLE PRECISION, DIMENSION(8,n/2)       :: Stka,Stkatot
      DOUBLE PRECISION, DIMENSION(10,1001)        :: Pdfk,Pdfktot
      DOUBLE PRECISION, DIMENSION(10,1001)        :: Pdfp,Pdfptot 
      DOUBLE PRECISION  :: tmp,tmp1,tmp2,tmp3,tmp4,Umax,DU,dux,duy
      INTEGER :: i,j,ir,imom,irx,idux,iduy,iu

      CALL derivk2(a,cy,1)
      CALL derivk2(a,cx,2)
      CALL fftp2d_complex_to_real(plancr,cx,rx,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,cy,ry,MPI_COMM_WORLD)
      DO ir = 1,n/2
      DO imom=1,8
         Stp(imom,ir) =0.
         Stk(imom,ir) =0.
         Stpa(imom,ir) =0.
         Stka(imom,ir) =0.
      END DO
      END DO
      Umax=50.0d0
      DU=2*Umax/1000
      DO iu = 1,1001
      DO ir=1,10
         Pdfk(ir,iu) =0.
         Pdfp(ir,iu) =0.
      END DO
      END DO

      tmp=1.0d0/dble(n)**2
      DO j = jsta,jend
        DO i = 1,n
          DO ir = 1,n/2
             irx = i+ir
             if (irx.gt.n) irx = irx - n
             DO imom=1,8
             Stp(imom,ir)=Stp(imom,ir)+((rx(irx,j)-rx(i,j))*tmp)**(imom)
             Stk(imom,ir)=Stk(imom,ir)+((ry(irx,j)-ry(i,j))*tmp)**(imom)
             Stpa(imom,ir)=Stpa(imom,ir)+abs((rx(irx,j)-rx(i,j))*tmp)**(imom)
             Stka(imom,ir)=Stka(imom,ir)+abs((ry(irx,j)-ry(i,j))*tmp)**(imom)
             END DO
          END DO
        END DO
      END DO
      CALL MPI_REDUCE(Stp,Stptot,8*(n/2),MPI_DOUBLE_PRECISION,MPI_SUM, &
                                         0, MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(Stk,Stktot,8*(n/2),MPI_DOUBLE_PRECISION,MPI_SUM, &
                                         0, MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(Stpa,Stpatot,8*(n/2),MPI_DOUBLE_PRECISION,MPI_SUM, &
                                         0, MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(Stka,Stkatot,8*(n/2),MPI_DOUBLE_PRECISION,MPI_SUM, &
                                         0, MPI_COMM_WORLD,ierr)
      !!!!!!!!!
      DO j = jsta,jend
        DO i = 1,n
          DO imom = 2,10
             ir   = n/2**(imom-1) 
             irx = i+ir
             if (irx.gt.n) irx = irx - n
             dux=((rx(irx,j)-rx(i,j))*tmp)
             duy=((ry(irx,j)-ry(i,j))*tmp)
             idux=int((Umax+dux)/DU+1.5d0)
             iduy=int((Umax+duy)/DU+1.5d0) 
             if (idux.le.1) idux=1
             if (iduy.le.1) iduy=1
             if (idux.ge.1001) idux=1001
             if (iduy.ge.1001) iduy=1001
             Pdfp(imom,idux) = Pdfp(imom,idux) +1.0D0
             Pdfk(imom,iduy) = Pdfk(imom,iduy) +1.0D0
          END DO
        END DO
      END DO
      CALL MPI_REDUCE(Pdfp,Pdfptot,10010,MPI_DOUBLE_PRECISION,MPI_SUM, &
                                         0, MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(Pdfk,Pdfktot,10010,MPI_DOUBLE_PRECISION,MPI_SUM, &
                                         0, MPI_COMM_WORLD,ierr)
      DO iu = 1,1001
      Pdfktot(1,iu) =-Umax+DU*(iu-1)
      Pdfptot(1,iu) =-Umax+DU*(iu-1)
      ENDDO
      !!!!!!!!!
      IF (myrank.eq.0) THEN
         OPEN(1,file='Structure_p.' // ext // '.txt')
         DO, ir=1,n/2
         WRITE(1,33) ( Stptot(imom,ir)/dble(n)**2, imom=1,8 )
         END DO
         CLOSE(1)
         OPEN(1,file='Structure_k.' // ext // '.txt')
         DO, ir=1,n/2
         WRITE(1,33) ( Stktot(imom,ir)/dble(n)**2, imom=1,8 )
         END DO
         CLOSE(1)
         OPEN(1,file='Structure_ap.' // ext // '.txt')
         DO, ir=1,n/2
         WRITE(1,33) ( Stpatot(imom,ir)/dble(n)**2, imom=1,8 )
         END DO
         CLOSE(1)
         OPEN(1,file='Structure_ak.' // ext // '.txt')
         DO, ir=1,n/2
         WRITE(1,33) ( Stkatot(imom,ir)/dble(n)**2, imom=1,8 )
         END DO
         CLOSE(1)
   33    FORMAT( 8E23.15 )
         OPEN(1,file='DUpdf_p.' // ext // '.txt')
         DO, ir=1,1001
         WRITE(1,44) ( Pdfptot(imom,ir), imom=1,10 )
         END DO
         CLOSE(1)
         OPEN(1,file='DUpdf_k.' // ext // '.txt')
         DO, ir=1,1001
         WRITE(1,44) ( Pdfktot(imom,ir), imom=1,10 )
         END DO
         CLOSE(1)
   44    FORMAT( 10E23.15 ) 
      ENDIF

      RETURN
      END SUBROUTINE structure

!*****************************************************************
      SUBROUTINE CFL_condition(cfl,c1,inu,nu,dt)
!-----------------------------------------------------------------

!        Parameters
!     cfl :cfl factor
!      c1 : stream fun
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
!%%%
      kcut=(dble(n)/3.0d0)
      tmp=sqrt(tmp1)/nrm+nu*kcut**(2*inu-1)
      dt = cfl/(kcut*tmp)
      !! if (myrank.eq.0) print*,"*",myrank,dt,cfl,tmp,sqrt(tmp1)/nrm,nu*kcut**(2*inu-1)
 
      RETURN
      END SUBROUTINE CFL_condition

!*****************************************************************
      SUBROUTINE forcing(iflow,f0,kup,kdn,dt,myseed,fk)
!-----------------------------------------------------------------

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
      INTEGER :: i,j,jj,iflow,ikup
      INTEGER :: seed,myseed
      DOUBLE PRECISION        :: f0,kup,kdn,dt,zz,zx,zy,kf
      DOUBLE PRECISION        :: tmp,tmp1,tmp2
      DOUBLE PRECISION        :: phase1,phase2,dkup

      kf=0.5d0*dble(kup+kdn)
      !if (myrank.eq.0) print*,"DBG pseudo: * iflow=",iflow,f0,kf
!!!!!!!  STEADY FORCING  !!!!!!!!!!!
      IF (iflow.eq.0) THEN
         DO j = jsta,jend
            DO i = 1,n
               r1(i,j) = sin(2*kup*pi*(dble(i)-1)/dble(n)) &
                       * sin(2*kdn*pi*(dble(j)-1)/dble(n))
            END DO
         END DO
         CALL fftp2d_real_to_complex(planrc,r1,fk,MPI_COMM_WORLD)
         CALL energy(fk,tmp1,1)
         CALL MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         tmp=f0/sqrt(tmp1)
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
                  fk(j,i) = tmp*fk(j,i)
               ELSE
                  fk(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
!!!!!!! (1)  RANDOM HOMOGENEOUS FORCING  !!!!!!!!!!!
      ELSE IF (iflow.eq.1) THEN
          !! print*,t,myrank,iflow,myseed
          DO i = ista,iend
           DO j = 1,n
             IF ((ka2(j,i).le.kup*kup).and.(ka2(j,i).ge.kdn*kdn)) THEN
                phase1 = 2*pi*randu(myseed)
                fk(j,i) = exp(im*phase1)
             ELSE
                fk(j,i) = 0.
             ENDIF
            END DO
          END DO
          CALL MPI_BCAST(seed ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          IF (myrank.eq.0) THEN
          i=1
          DO j = 1,n
            jj=n-j+2
            if (jj.gt.n) jj=1
            fk(j,i) = conjg(fk(jj,i))
          END DO
          ENDIF
      CALL energy(fk,tmp1,1)
      CALL MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      tmp=f0/sqrt(0.5d0*tmp1*dt)
      DO i = ista,iend
         DO j = 1,n
            IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
              fk(j,i) = tmp*fk(j,i)
            ELSE
              fk(j,i) = 0.0d0
            ENDIF
         END DO
      END DO
!!!!!!!  (2)  POINT VORTEX FORCING   !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  PSI = [cos(x)+cos(y)]^(2n+1) with n>>1  (n=kup)
      ELSE IF (iflow.eq.2) THEN
         ikup = 2*int(kf)+1
         dkup = dble(kf)
         DO j = jsta,jend
            DO i = 1,n
               r1(i,j) = 0.5d0*cos(2*pi*(dble(i)-1)/dble(n)) &
                       + 0.5d0*cos(2*pi*(dble(j)-1)/dble(n))
               !! r1(i,j) = r1(i,j)**ikup
               !!!!!!!!!!!!!!!!!!!
               tmp     = r1(i,j)
               tmp1    = dkup*dlog(r1(i,j)**2)  
               tmp2    = tmp1
               if (tmp2.lt.-80.d0) tmp2=-80.0d0
               r1(i,j) = r1(i,j)*exp(tmp2)
            END DO
         END DO
         CALL fftp2d_real_to_complex(planrc,r1,fk,MPI_COMM_WORLD)
         CALL energy(fk,tmp1,1)
         CALL MPI_BCAST(tmp1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         tmp=f0/sqrt(tmp1)
         DO i = ista,iend
            DO j = 1,n
               IF ((ka2(j,i).le.5*kmax2).and.(ka2(j,i).ge.kmin2)) THEN
                  fk(j,i) = -tmp*fk(j,i) !!/ka2(j,i)
               ELSE
                  fk(j,i) = 0.0d0
               ENDIF
            END DO
         END DO
!!!!!!!  (3) POINT VORTEX FORCING  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  PSI = exp(-kup**2(x**2+y**2) )
      ELSE IF (iflow.eq.3) THEN
         !if (myrank.eq.0) print*,"DBG pseudo: iflow=",iflow,f0,kf
         ikup = 2*int(kf)+1
         dkup = dble(kf)/2
         DO j = jsta,jend
            DO i = 1,n
               zz=(2*pi*(dble(i-n/2)-1)/dble(n))**2
               zz=(2*pi*(dble(j-n/2)-1)/dble(n))**2+zz
               zz=0.5d0*dkup*dkup*zz
               if (zz.gt.80.d0) zz=80.0d0
               r1(i,j) = f0*exp(-zz)
            END DO
         END DO
         CALL fftp2d_real_to_complex(planrc,r1,c1,MPI_COMM_WORLD)
         tmp1= pi/dkup/2
         tmp2= pi/dkup/2
         CALL shift(c1,fk,tmp1,tmp2)
         tmp1= -pi/dkup/2
         tmp2= -pi/dkup/2
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = ista,iend
            DO j = 1,n
            fk(j,i) = fk(j,i) +c2(j,i)          
            END DO
         END DO
         tmp1= +pi/dkup/2
         tmp2= -pi/dkup/2
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = ista,iend
            DO j = 1,n
            fk(j,i) = fk(j,i) -c2(j,i)
            END DO
         END DO
         tmp1= -pi/dkup/2
         tmp2= +pi/dkup/2
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = ista,iend
            DO j = 1,n
            fk(j,i) = fk(j,i) -c2(j,i)
            END DO
         END DO
         DO i = ista,iend
            DO j = 1,n
            IF (ka2(j,i).ge.0.1) THEN
                  fk(j,i) = fk(j,i)/ka2(j,i)
               ELSE
                  fk(j,i) = 0.0d0
            ENDIF
            END DO
         END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!! (4) RING of FIRE !!!!
         ELSE IF (iflow.eq.4) THEN
         !if (myrank.eq.0) print*,"DBG pseudo: iflow=",iflow,f0,kf
         ikup = 2*int(kf)+1
         dkup = dble(kf)/2
         DO j = jsta,jend
            DO i = 1,n
               zz=(2*pi*(dble(i-n/2)-1)/dble(n))**2
               zz=(2*pi*(dble(j-n/2)-1)/dble(n))**2+zz
               zz=0.5d0*dkup*dkup*zz
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
         phase1 = 2*pi*randu(myseed)
         CALL MPI_BCAST(phase1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         tmp1= pi/dkup/2*sin(phase1)
         tmp2= pi/dkup/2*cos(phase1)
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = ista,iend
            DO j = 1,n
            fk(j,i) = fk(j,i) +c2(j,i)
            END DO
         END DO
         tmp1= -tmp1
         tmp2= -tmp2
         CALL shift(c1,c2,tmp1,tmp2)
         DO i = ista,iend
            DO j = 1,n
            fk(j,i) = fk(j,i) -c2(j,i)
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
!!!!!!!  (5) VORTEX LINE FORCING  !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  PSI = exp(-(kup*x)**2)+exp(-(kup*y)**2) 
      ELSE IF (iflow.gt.4) THEN
         ikup = 2*int(kf)+1
         dkup = dble(kf)/1.5
         DO j = jsta,jend
            DO i = 1,n
               zx=(2*pi*(dble(i-n/2)-1)/dble(n))**2
               zy=(2*pi*(dble(j-n/2)-1)/dble(n))**2
               zx=0.5d0*dkup*dkup*zx
               zy=0.5d0*dkup*dkup*zy
               if (zx.gt.80.d0) zx=80.0d0
               if (zy.gt.80.d0) zy=80.0d0
               r1(i,j) = f0*(exp(-zx)-exp(-zy))
            END DO
         END DO
         CALL fftp2d_real_to_complex(planrc,r1,fk,MPI_COMM_WORLD)   
      ENDIF      
      RETURN
      END SUBROUTINE forcing

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
      USE mpivars
      IMPLICIT NONE

      DOUBLE COMPLEX, DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION        :: dx,dy
      INTEGER :: i,j

         DO i = ista,iend
            DO j = 1,n
!!               IF ((ka2(j,i).le.kmax2).and.(ka2(j,i).ge.kmin2)) THEN
                  b(j,i) = a(j,i) *exp(-im*(ka(i)*dx+ka(j)*dy)) 
!!               ELSE
!!                  b(j,i) = 0.0d0
!!               ENDIF
            END DO
         END DO

      RETURN
      END SUBROUTINE shift


!*****************************************************************
      SUBROUTINE fractal05(f)
!-----------------------------------------------------------------
!
! Two-dimensional derivative of the matrix 'a'
!
! Parameters
!     a  : output matrix
      USE mpivars
      IMPLICIT NONE

      INTEGER, DIMENSION(1024) :: f
      INTEGER :: i,j,npieces,ipiece,piecesize,icut,ncut,zero

         DO i = 1,1024
         f(i) = 1
         END DO

         DO icut=0,4
            npieces = 4**(icut+1) 
            piecesize = int(1024/npieces) 
            DO ipiece = 0,npieces-1
                zero=1
                if (2*int(ipiece/2).lt.ipiece) zero=0
                DO i = 1,piecesize
                   j=ipiece*piecesize+i
                   f(j)=f(j)*zero
                   !! if (myrank.eq.0) print*,"sbrt",j,i,zero,f(j)
                ENDDO ! i
            ENDDO ! ipiece 
         ENDDO ! icut   

      RETURN
      END SUBROUTINE fractal05




