! Minimal finite-dimensional host for the production metallic orthogonalizer.
MODULE kinds
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.D0)
END MODULE
MODULE io_global
  INTEGER, PARAMETER :: stdout=6
  LOGICAL, PARAMETER :: ionode=.TRUE.
END MODULE
MODULE constants
  USE kinds
  REAL(DP), PARAMETER :: pi=3.1415926535897932384626433832795_DP
  REAL(DP), PARAMETER :: sqrtpm1=1.D0/SQRT(pi), sqrt2=SQRT(2.D0)
END MODULE
MODULE klist
  USE kinds
  LOGICAL :: lgauss=.TRUE., ltetra=.FALSE.
  REAL(DP) :: degauss=0.04_DP, wk(2)=1.D0
  INTEGER :: ngauss=-1
END MODULE
MODULE noncollin_module
  LOGICAL :: noncolin=.TRUE.
  INTEGER :: npol=2
END MODULE
MODULE wvfct
  USE kinds
  INTEGER :: npwx=2, nbnd=4
  REAL(DP) :: wg(4,2)=0.D0, et(4,2)
END MODULE
MODULE ener
  USE kinds
  REAL(DP) :: ef=0.D0
END MODULE
MODULE becmod
  USE kinds
  TYPE bec_type
     INTEGER :: unused
  END TYPE
  TYPE(bec_type) :: becp
CONTAINS
  SUBROUTINE calbec(offload,n,vkb,psi,bec,m)
    INTEGER :: offload,n,m
    COMPLEX(DP) :: vkb(:,:),psi(:,:)
    TYPE(bec_type) :: bec
    ERROR STOP 'ultrasoft branch is outside this norm-conserving probe'
  END SUBROUTINE
END MODULE
MODULE uspp
  USE kinds
  COMPLEX(DP) :: vkb(4,1)
  LOGICAL :: okvan=.FALSE.
END MODULE
MODULE mp_bands
  LOGICAL :: use_bgrp_in_hpsi=.FALSE.
  INTEGER :: inter_bgrp_comm=0, intra_bgrp_comm=0
END MODULE
MODULE mp
  USE kinds
  INTERFACE mp_sum
    MODULE PROCEDURE sum_real, sum_complex
  END INTERFACE
CONTAINS
  SUBROUTINE sum_real(a,comm)
    REAL(DP) :: a(:,:)
    INTEGER :: comm
  END SUBROUTINE
  SUBROUTINE sum_complex(a,comm)
    COMPLEX(DP) :: a(:,:)
    INTEGER :: comm
  END SUBROUTINE
END MODULE
MODULE xc_lib
CONTAINS
  LOGICAL FUNCTION exx_is_active()
    exx_is_active=.FALSE.
  END FUNCTION
END MODULE
MODULE control_flags
  LOGICAL :: gamma_only=.FALSE.
  INTEGER :: offload_type=0
END MODULE
MODULE gvect
  INTEGER :: gstart=1
END MODULE
MODULE control_lr
  USE kinds
  REAL(DP) :: alpha_pv=1.D0
  INTEGER :: nbnd_occ(2)=4
END MODULE
MODULE dfpt_tetra_mod
  USE kinds
CONTAINS
  REAL(DP) FUNCTION dfpt_tetra_beta(j,i,k)
    INTEGER :: j,i,k
    ERROR STOP 'tetrahedra are outside this smearing probe'
  END FUNCTION
END MODULE
SUBROUTINE start_clock(name)
  CHARACTER(*) :: name
END SUBROUTINE
SUBROUTINE stop_clock(name)
  CHARACTER(*) :: name
END SUBROUTINE
SUBROUTINE errore(routine,message,code)
  CHARACTER(*) :: routine,message
  INTEGER :: code
  WRITE(*,*) routine,message,code
  STOP 1
END SUBROUTINE
SUBROUTINE divide(comm,n,first,last)
  INTEGER :: comm,n,first,last
  first=1
  last=n
END SUBROUTINE
SUBROUTINE s_psi_acc(npwx,npw,m,psi,spsi)
  USE kinds
  INTEGER :: npwx,npw,m
  COMPLEX(DP) :: psi(2*npwx,m),spsi(2*npwx,m)
  spsi=psi
END SUBROUTINE
SUBROUTINE dgemm()
  ERROR STOP 'real gamma-only BLAS is outside this spinor probe'
END SUBROUTINE
SUBROUTINE mydger()
  ERROR STOP 'real gamma-only BLAS is outside this spinor probe'
END SUBROUTINE
SUBROUTINE zgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  USE kinds
  CHARACTER :: ta,tb
  INTEGER :: m,n,k,lda,ldb,ldc,i,j,l
  COMPLEX(DP) :: alpha,beta,a(lda,*),b(ldb,*),c(ldc,*),aa,bb,value
  DO j=1,n
     DO i=1,m
        value=0.D0
        DO l=1,k
           IF (ta=='N') THEN
              aa=a(i,l)
           ELSE IF (ta=='C') THEN
              aa=CONJG(a(l,i))
           ELSE
              ERROR STOP 'unsupported first BLAS transpose'
           ENDIF
           IF (tb=='N') THEN
              bb=b(l,j)
           ELSE IF (tb=='C') THEN
              bb=CONJG(b(j,l))
           ELSE
              ERROR STOP 'unsupported second BLAS transpose'
           ENDIF
           value=value+aa*bb
        ENDDO
        c(i,j)=alpha*value+beta*c(i,j)
     ENDDO
  ENDDO
END SUBROUTINE
PROGRAM probe
  USE kinds
  USE klist
  USE control_lr
  USE wvfct
  USE hubbard_nc_response, ONLY: hubbard_kramers_partner_nc,hub_spin_index,hubbard_dv_from_dns_nc
  IMPLICIT NONE
  CHARACTER(16) :: action
  COMPLEX(DP) :: evq(4,4),dvpsi(4,4),dpsi(4,4),mat(4,4)
  COMPLEX(DP) :: blocks(2,2,4,1),partner(2,2,4,1)
  INTEGER :: i,j,s,t,m,n
  REAL(DP) :: re,im,strength(1)
  CALL GET_COMMAND_ARGUMENT(1,action)
  IF (TRIM(action)=='weight') THEN
     READ(*,*) ngauss,degauss,alpha_pv,nbnd_occ
     READ(*,*) et(:,1)
     READ(*,*) et(:,2)
     DO j=1,4
        DO i=1,4
           READ(*,*) re,im
           evq(i,j)=CMPLX(re,im,DP)
        ENDDO
     ENDDO
     DO j=1,4
        DO i=1,4
           READ(*,*) re,im
           dvpsi(i,j)=CMPLX(re,im,DP)
        ENDDO
     ENDDO
     dpsi=0.D0
     CALL orthogonalize(dvpsi,evq,1,2,dpsi,2,.FALSE.)
     ! Only source occupied columns are defined by the production routine.
     dvpsi(:,nbnd_occ(1)+1:4)=0.D0
     mat=dvpsi
  ELSE IF (TRIM(action)=='partner' .OR. TRIM(action)=='dudarev') THEN
     IF (TRIM(action)=='dudarev') READ(*,*) strength
     DO j=1,4
        DO i=1,4
           READ(*,*) re,im
           mat(i,j)=CMPLX(re,im,DP)
        ENDDO
     ENDDO
     DO s=1,2
        DO t=1,2
           DO m=1,2
              DO n=1,2
                 blocks(m,n,hub_spin_index(s,t),1)=mat(m+2*(s-1),n+2*(t-1))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     IF (TRIM(action)=='partner') THEN
        CALL hubbard_kramers_partner_nc(2,1,blocks,partner)
     ELSE
        CALL hubbard_dv_from_dns_nc(2,1,strength,blocks,partner)
     ENDIF
     DO s=1,2
        DO t=1,2
           DO m=1,2
              DO n=1,2
                 mat(m+2*(s-1),n+2*(t-1))=partner(m,n,hub_spin_index(s,t),1)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ELSE
     ERROR STOP 'unknown probe action'
  ENDIF
  DO j=1,4
     DO i=1,4
        WRITE(*,'(2ES26.17)') REAL(mat(i,j),DP),AIMAG(mat(i,j))
     ENDDO
  ENDDO
END PROGRAM
