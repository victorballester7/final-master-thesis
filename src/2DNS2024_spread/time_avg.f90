      PROGRAM TIME_AVERAGE
!=================================================================
!=================================================================

      IMPLICIT NONE

! Integration parameters
!     ord  : order of the Runge-Kutta method used

      INTEGER, PARAMETER :: ngrd = 512
      INTEGER, PARAMETER :: nsta = 1902
      INTEGER, PARAMETER :: nend = 1920
      INTEGER :: i,j,k,ifile,NL,ru,rd,rc,rt
      DOUBLE PRECISION :: epsilon,pi
      DOUBLE PRECISION :: x,y,z,r,r1,r2
      DOUBLE PRECISION :: eng,ens 
      DOUBLE PRECISION :: tmp,tmp1,tmp2,tmp3,tmp4
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: SFP   
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: SFK   
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: SFAP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: SFAK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: ESP   
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: ETR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: VTR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:)   :: ETMP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: STMP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: VTMP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: PDFP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: PDFK
      CHARACTER     :: th,c,d,u
      CHARACTER*4   :: ext


!filename1 = 'kspectrum.' 
!filename2 = 'Structure_k.'  
!filename3 = 'Structure_p.'  
!filename4 = 'vectrans_euu.'  
!filename5 = 'vectrans_vrt.'
      pi = 3.14159265359d0
      NL=int(ngrd/2)
      ALLOCATE( ESP(2,NL)   )
      ALLOCATE( ETR(3,NL)   )
      ALLOCATE( VTR(3,NL)   )
      ALLOCATE( SFP(9,NL)   )
      ALLOCATE( SFK(9,NL)   )
      ALLOCATE( SFAP(9,NL)   )
      ALLOCATE( SFAK(9,NL)   )
      ALLOCATE( ETMP(NL+1)  )
      ALLOCATE( VTMP(2,NL+1))
      ALLOCATE( STMP(8,NL)  )
      ALLOCATE( PDFP(10,1001))
      ALLOCATE( PDFK(10,1001))
      

      DO i=1,NL
      ESP(1,i) = i
      ETR(1,i) = i
      VTR(1,i) = i
      SFP(1,i) = (i)*pi/(NL-1)
      SFK(1,i) = (i)*pi/(NL-1)
      SFAP(1,i) = (i)*pi/(NL-1)
      SFAK(1,i) = (i)*pi/(NL-1)
      ESP(2,i) = 0.0D0
      ETR(2,i) = 0.0D0
      VTR(2,i) = 0.0D0
       DO j=2,9
          SFP(j,i) =0
          SFK(j,i) =0
       ENDDO   
      ENDDO

      DO ifile=nsta,nend
      ru = mod(ifile,10   )
      rd = mod(ifile,100  )-ru
      rc = mod(ifile,1000 )-rd-ru
      rt = ifile-rc-rd-ru
      ru = 48 + ru
      rd = 48 + rd/10
      rc = 48 + rc/100
      rt = 48 + rt/1000
      th = char(rt)
       c = char(rc)
       d = char(rd)
       u = char(ru)
       ext = th // c // d // u
       !!!!!!!!!
       print*,'kspectrum.'//ext//'.txt'
       OPEN(1,file='kspectrum.' // ext //  '.txt',status='unknown')
       DO i=1,NL
       READ(1,*) ETMP(i)
       ESP(2,i) = ETMP(i) + ESP(2,i)
       ENDDO
       CLOSE(1)
       !!!!!!!!!
       !print*,'vectrans_euu.'// ext//'.txt'
       OPEN(1,file='vectrans_euu.' // ext //  '.txt',status='unknown')
       DO i=1,NL
       READ(1,*) VTMP(1,i),VTMP(2,i)
       ETR(2,i) = VTMP(1,i) + ETR(2,i)
       ETR(3,i) = VTMP(2,i) + ETR(3,i)
       ENDDO
       CLOSE(1)
       !!!!!!!!!
       !print*,'vectrans_vrt.' // ext //  '.txt'
       OPEN(1,file='vectrans_vrt.' // ext //  '.txt',status='unknown')
       DO i=1,NL
       READ(1,*) VTMP(1,i),VTMP(2,i)
       VTR(2,i) = VTMP(1,i) + VTR(2,i)
       VTR(3,i) = VTMP(2,i) + VTR(3,i)
       ENDDO
       CLOSE(1)
       !!!!!!!!!
       !print*,'Structure_p.' // ext //  '.txt'
       OPEN(1,file='Structure_p.' // ext //  '.txt',status='unknown')
       DO i=1,NL
       READ(1,*) (STMP(j,i), j=1,8)
       DO j=1,8
       SFP(j+1,i) = STMP(j,i) + SFP(j+1,i)
       ENDDO
       ENDDO
       CLOSE(1)
       !!!!!!
       !print*,'Structure_k.' // ext //  '.txt'
       OPEN(1,file='Structure_k.' // ext //  '.txt',status='unknown')
       DO i=1,NL
       READ(1,*) (STMP(j,i), j=1,8)
       DO j=1,8
       SFK(j+1,i) = STMP(j,i) + SFK(j+1,i)
       ENDDO
       ENDDO
       CLOSE(1)
       !!!!!!
       OPEN(1,file='Structure_ap.' // ext //  '.txt',status='unknown')
       DO i=1,NL
       READ(1,*) (STMP(j,i), j=1,8)
       DO j=1,8
       SFAP(j+1,i) = STMP(j,i) + SFAP(j+1,i)
       ENDDO
       ENDDO
       CLOSE(1)
       !!!!!!
       !print*,'Structure_ak.' // ext //  '.txt'
       OPEN(1,file='Structure_ak.' // ext //  '.txt',status='unknown')
       DO i=1,NL
       READ(1,*) (STMP(j,i), j=1,8)
       DO j=1,8
       SFAK(j+1,i) = STMP(j,i) + SFAK(j+1,i)
       ENDDO
       ENDDO
       CLOSE(1)
       !!!!!!
       OPEN(1,file='DUpdf_p.' // ext //  '.txt',status='unknown')
       DO i=1,1001
       READ(1,*) (PDFP(j,i), j=1,10)
       PDFP(1,i) = PDFP(1,i)
       DO j=2,10
       PDFP(j,i) = PDFP(j,i) + PDFP(j,i)
       ENDDO
       ENDDO
       CLOSE(1)
       !!!!!!
       OPEN(1,file='DUpdf_k.' // ext //  '.txt',status='unknown')
       DO i=1,1001
       READ(1,*) (PDFK(j,i), j=1,10)
       PDFK(1,i) = PDFK(1,i)
       DO j=2,10
       PDFK(j,i) = PDFK(j,i) + PDFK(j,i)
       ENDDO
       ENDDO
       CLOSE(1)
       !!!!!!
       !print*,"ifile=",ifile
      ENDDO

      tmp=1.0d0/(nend-nsta+1) 
      DO i=1,NL
      ESP(2,i) = ESP(2,i)*tmp
      ETR(2,i) = ETR(2,i)*tmp
      ETR(3,i) = ETR(3,i)*tmp
      VTR(2,i) = VTR(2,i)*tmp
      VTR(3,i) = VTR(3,i)*tmp
      DO j=2,9
        SFK(j,i) = SFK(j,i)*tmp
        SFP(j,i) = SFP(j,i)*tmp
        SFAK(j,i) = SFAK(j,i)*tmp
        SFAP(j,i) = SFAP(j,i)*tmp
      ENDDO
      ENDDO

      DO j=2,10
       tmp=0.0
       DO i=1,1001
       tmp=tmp+PDFP(j,i)
       ENDDO
       tmp=10.0/tmp
       DO i=1,1001
       PDFP(j,i)=PDFP(j,i)*tmp
       ENDDO
       !!!!!
       tmp=0.0
       DO i=1,1001
       tmp=tmp+PDFK(j,i)
       ENDDO
       tmp=10.0/tmp
       DO i=1,1001
       PDFK(j,i)=PDFK(j,i)*tmp
       ENDDO
      ENDDO 

      OPEN(10,file='kspectrum.avrg.txt',position='append')
      DO i=1,NL
      WRITE(10,11)  (ESP(j,i), j=1,2  )
      ENDDO
      CLOSE(10)
 11   FORMAT( E22.14,E22.14 )
      OPEN(10,file='vectrans_euu.avrg.txt',position='append')
      DO i=1,NL
      WRITE(10,12)  (ETR(j,i), j=1,3  )
      ENDDO
      CLOSE(10)
      OPEN(10,file='vectrans_vrt.avrg.txt',position='append')
      DO i=1,NL
      WRITE(10,12)  (VTR(j,i), j=1,3  )
      ENDDO
      CLOSE(10)
 12   FORMAT( E22.14,E22.14,E22.14 ) 
      OPEN(10,file='Structure_p.avrg.txt',position='append')
      DO i=1,NL
      WRITE(10,13)  (SFP(j,i), j=1,9  )
      ENDDO
      CLOSE(10)
      OPEN(10,file='Structure_k.avrg.txt',position='append')
      DO i=1,NL
      WRITE(10,13)  (SFK(j,i), j=1,9  )
      ENDDO
      CLOSE(10)
      OPEN(10,file='Structure_ap.avrg.txt',position='append')
      DO i=1,NL
      WRITE(10,13)  (SFAP(j,i), j=1,9  )
      ENDDO
      CLOSE(10)
      OPEN(10,file='Structure_ak.avrg.txt',position='append')
      DO i=1,NL
      WRITE(10,13)  (SFAK(j,i), j=1,9  )
      ENDDO
      CLOSE(10)
 13   FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14,E22.14, E22.14,E22.14,E22.14 )
      OPEN(10,file='DUpdf_p.avrg.txt',position='append')
      DO i=1,1001
      WRITE(10,14)  (PDFP(j,i), j=1,10  )
      ENDDO
      CLOSE(10)
      OPEN(10,file='DUpdf_k.avrg.txt',position='append')
      DO i=1,1001
      WRITE(10,14)  (PDFK(j,i), j=1,10  )
      ENDDO
      CLOSE(10)
  14  FORMAT( E22.14,E22.14,E22.14,E22.14,E22.14,E22.14,E22.14,E22.14,E22.14, E22.14 )
      END PROGRAM TIME_AVERAGE


