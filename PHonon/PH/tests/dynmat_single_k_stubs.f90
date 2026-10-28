! Minimal environment for the production dynmat_us single-k regression.
! A translated one-plane-wave separable projector has exactly zero curvature.
MODULE kinds
  INTEGER, PARAMETER :: dp=KIND(1.0D0)
END MODULE
MODULE constants
  USE kinds
  REAL(dp), PARAMETER :: tpi=6.2831853071795864769_dp
END MODULE
MODULE ions_base
  USE kinds
  INTEGER :: nat=1,nsp=1,ityp(1)=1
  REAL(dp) :: tau(3,1)=0
END MODULE
MODULE uspp
  USE kinds
  INTEGER :: nkb=1
  COMPLEX(dp) :: vkb(1,1)=1
END MODULE
MODULE scf
  USE kinds
  TYPE density
    REAL(dp) :: of_r(1,1)=0
  END TYPE
  TYPE(density) :: rho
END MODULE
MODULE fft_base
  TYPE descriptor
    INTEGER :: nnr=1,nl(1)=1
  END TYPE
  TYPE(descriptor) :: dfftp
END MODULE
MODULE fft_interfaces
  USE kinds
  USE fft_base
CONTAINS
  SUBROUTINE fwfft(label,data,desc)
    CHARACTER(*) :: label
    COMPLEX(dp) :: data(:)
    TYPE(descriptor) :: desc
  END SUBROUTINE
END MODULE
MODULE buffers
  USE kinds
  INTEGER :: reads=0
  COMPLEX(dp) :: direct(2,2)
CONTAINS
  SUBROUTINE get_buffer(data,length,unit,record)
    COMPLEX(dp) :: data(:,:)
    INTEGER :: length,unit,record
    reads=reads+1
    data(:,1)=direct(:,record)
  END SUBROUTINE
END MODULE
MODULE gvect
  USE kinds
  INTEGER :: ngm=1,igtongl(1)=1
  REAL(dp) :: g(3,1)=0
END MODULE
MODULE wvfct
  USE kinds
  INTEGER :: npwx=1,nbnd=1
  REAL(dp) :: wg(1,2)=0,et(1,2)=0
END MODULE
MODULE lsda_mod
  LOGICAL :: lsda=.FALSE.
  INTEGER :: current_spin=1,isk(2)=1,nspin=4
END MODULE
MODULE vlocal
  USE kinds
  REAL(dp) :: vloc(1,1)=0
END MODULE
MODULE klist
  USE kinds
  INTEGER :: ngk(2)=1,igk_k(1,2)=1
  REAL(dp) :: xk(3,2)=0
END MODULE
MODULE wavefunctions
  USE kinds
  COMPLEX(dp) :: evc(2,1)
END MODULE
MODULE cell_base
  USE kinds
  REAL(dp) :: omega=1,tpiba2=1
END MODULE
MODULE uspp_param
  INTEGER :: nh(1)=1,nhm=1
END MODULE
MODULE noncollin_module
  LOGICAL :: noncolin=.TRUE.,domag=.TRUE.
  INTEGER :: npol=2
END MODULE
MODULE becmod
  USE kinds
  TYPE bec_type
    COMPLEX(dp), ALLOCATABLE :: nc(:,:,:),k(:,:)
  END TYPE
CONTAINS
  SUBROUTINE allocate_bec_type(nkb,nbnd,bec)
    INTEGER :: nkb,nbnd
    TYPE(bec_type) :: bec
    ALLOCATE(bec%nc(nkb,2,nbnd),bec%k(nkb,nbnd))
    bec%nc=0; bec%k=0
  END SUBROUTINE
  SUBROUTINE deallocate_bec_type(bec)
    TYPE(bec_type) :: bec
    DEALLOCATE(bec%nc,bec%k)
  END SUBROUTINE
  SUBROUTINE beccopy(source,dest,nkb,nbnd)
    TYPE(bec_type) :: source,dest
    INTEGER :: nkb,nbnd
    dest%nc=source%nc; dest%k=source%k
  END SUBROUTINE
  SUBROUTINE calbec(npw,beta,psi,bec)
    INTEGER :: npw
    COMPLEX(dp) :: beta(:,:),psi(:,:)
    TYPE(bec_type) :: bec
    INTEGER :: spin
    DO spin=1,2
      bec%nc(1,spin,1)=DOT_PRODUCT(beta(:,1),psi(spin:spin,1))
    END DO
  END SUBROUTINE
END MODULE
MODULE modes
  USE kinds
  COMPLEX(dp) :: u(3,3)=0
END MODULE
MODULE dynmat
  USE kinds
  COMPLEX(dp) :: dyn(3,3)=0
END MODULE
MODULE phus
  USE becmod
  TYPE(bec_type) :: alphap(3,2)
END MODULE
MODULE units_lr
  INTEGER :: iuwfc=1,lrwfc=4
END MODULE
MODULE io_global
  INTEGER :: stdout=6
END MODULE
MODULE mp_pools
  INTEGER :: my_pool_id=0,inter_pool_comm=0
END MODULE
MODULE mp_bands
  INTEGER :: intra_bgrp_comm=0
END MODULE
MODULE mp
  USE kinds
CONTAINS
  SUBROUTINE mp_sum(data,comm)
    COMPLEX(dp) :: data(:,:)
    INTEGER :: comm
  END SUBROUTINE
END MODULE
MODULE Coul_cut_2D
  LOGICAL :: do_cutoff_2D=.FALSE.
END MODULE
MODULE Coul_cut_2D_ph
  USE kinds
CONTAINS
  SUBROUTINE cutoff_dynmat0(dyn,rho)
    COMPLEX(dp) :: dyn(:,:),rho(:)
  END SUBROUTINE
END MODULE
MODULE lrus
  USE becmod
  TYPE(bec_type) :: becp1(2)
END MODULE
MODULE qpoint
  INTEGER :: nksq=1,ikks(2)=[1,2]
END MODULE
MODULE control_lr
  INTEGER :: nbnd_occ(2)=1
  LOGICAL :: lgamma=.TRUE.
END MODULE
MODULE uspp_init
  USE kinds
CONTAINS
  SUBROUTINE init_us_2(npw,igk,xk,beta)
    INTEGER :: npw,igk(*)
    REAL(dp) :: xk(*)
    COMPLEX(dp) :: beta(:,:)
  END SUBROUTINE
END MODULE
SUBROUTINE start_clock(label)
  CHARACTER(*) :: label
END SUBROUTINE
SUBROUTINE stop_clock(label)
  CHARACTER(*) :: label
END SUBROUTINE
SUBROUTINE compute_deff_nc(deff,energy)
  USE kinds
  COMPLEX(dp) :: deff(1,1,1,4)
  REAL(dp) :: energy
  ! A nontrivial Hermitian spin-dependent separable potential.
  deff(1,1,1,:)=[CMPLX(2.0_dp,0.0_dp,dp),CMPLX(0.2_dp,0.3_dp,dp), &
                CMPLX(0.2_dp,-0.3_dp,dp),CMPLX(0.7_dp,0.0_dp,dp)]
END SUBROUTINE
SUBROUTINE compute_deff(deff,energy)
  USE kinds
  REAL(dp) :: deff(1,1,1),energy
  deff=1
END SUBROUTINE
SUBROUTINE addusdynmat(dyn)
  USE kinds
  COMPLEX(dp) :: dyn(3,3)
END SUBROUTINE
SUBROUTINE rotate_pattern_add(nat,u,dyn,work)
  USE kinds
  INTEGER :: nat
  COMPLEX(dp) :: u(3*nat,3*nat),dyn(3*nat,3*nat),work(3*nat,3*nat)
  dyn=dyn+MATMUL(CONJG(TRANSPOSE(u)),MATMUL(work,u))
END SUBROUTINE

PROGRAM probe
  USE kinds
  USE buffers
  USE klist
  USE modes
  USE dynmat
  USE phus
  USE lrus
  USE qpoint
  USE wavefunctions
  USE wvfct
  USE noncollin_module
  IMPLICIT NONE
  INTEGER :: scenario,ik,alpha
  CHARACTER(20) :: argument
  CALL get_command_argument(1,argument)
  READ(argument,*) scenario
  direct(:,1)=[CMPLX(0.3_dp,0.4_dp,dp),CMPLX(0.2_dp,-0.1_dp,dp)]
  direct(:,2)=[CMPLX(-0.1_dp,0.3_dp,dp),CMPLX(0.5_dp,0.2_dp,dp)]
  xk(:,1)=[0.21_dp,-0.17_dp,0.32_dp]
  xk(:,2)=[-0.13_dp,0.27_dp,0.11_dp]
  nksq=1
  IF (scenario==2) nksq=2
  IF (scenario==3) domag=.FALSE.
  wg=1.0_dp/REAL(nksq,dp)
  DO ik=1,2
    CALL allocate_bec_type(1,1,becp1(ik))
    becp1(ik)%nc(1,:,1)=direct(:,ik)
    DO alpha=1,3
      CALL allocate_bec_type(1,1,alphap(alpha,ik))
      alphap(alpha,ik)%nc(1,:,1)=CMPLX(0.0_dp,xk(alpha,ik),dp)*direct(:,ik)
    END DO
  END DO
  DO alpha=1,3
    u(alpha,alpha)=1
  END DO
  ! Reproduce dnsq_bare's exit state: the cached spinor is the auxiliary
  ! J*conjg(psi), but all stored becp1/alphap belong to the direct state.
  evc(:,1)=[CONJG(direct(2,1)),-CONJG(direct(1,1))]
  IF (scenario==3) evc(:,1)=direct(:,1)
  CALL dynmat_us()
  WRITE(*,'(I0,1X,ES24.16)') reads,MAXVAL(ABS(dyn))
END PROGRAM
