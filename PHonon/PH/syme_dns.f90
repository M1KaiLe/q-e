!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE syme_dns (ldim, npe, dns)
  !-----------------------------------------------------------------------
  !! DFPT+U: This routine symmetrizes the first order variation of 
  !! the occupation matrices dns due to the perturbation with the 
  !! electric field.
  !
  !! Written  by S. de Gironcoli and A. Floris.  
  !! Modified by I. Timrov (01.10.2018).
  !
  USE kinds,             ONLY : DP
  USE constants,         ONLY : tpi
  USE ions_base,         ONLY : nat, ityp
  USE ldaU,              ONLY : Hubbard_l, is_hubbard, nwfcU, Hubbard_lmax
  USE lsda_mod,          ONLY : lsda, nspin
  USE lr_symm_base,      ONLY : nsymq, irgq, minus_q, irotmq, rtau
  USE uspp_param,        ONLY : upf
  USE symm_base,         ONLY : d1, d2, d3, nsym, irt, s, invs

  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldim, npe
  COMPLEX(DP), INTENT(INOUT) :: dns(ldim,ldim,nspin,nat,npe)
  !
  ! ... local variables
  !
  INTEGER :: nt, n, l, ip, jp, na, nb, is, m1, m2, &
             m0, m00, isym, irot
  COMPLEX(DP), ALLOCATABLE :: dnr(:,:,:,:,:), dnraux(:,:,:,:,:)
  !
  IF ((nsym==1) .AND. (.NOT.minus_q)) RETURN
  !
  ! Initialization
  !
  ! D_Sl for l=1, l=2 and l=3 are already initialized, for l=0 D_S0 is 1
  !
  ALLOCATE (dnraux(ldim,ldim,nspin,nat,npe))  
  ALLOCATE (dnr(ldim,ldim,nspin,nat,npe))  
  !
  dnraux = (0.d0, 0.d0)
  dnr    = (0.d0, 0.d0)
  !
  ! Impose hermiticity of dns_{m1,m2, is,na,ip} for m1<->m2
  ! and put it in dnr zeroing dns
  !
  DO ip = 1, npe
     DO na = 1, nat  
        nt = ityp(na)
        DO is = 1, nspin  
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                 dnr(m1, m2, is, na, ip) = 0.5d0 * ( dns(m1, m2, is, na, ip) &
                                                   + dns(m2, m1, is, na, ip) )
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  !
  dns = (0.d0, 0.d0)
  !
  ! Symmetrize with -q if present (output overwritten on dnr)
  ! We use the s matrix in cryst. coord.
  !
  IF (minus_q) THEN
     DO ip = 1, npe
        DO na = 1, nat  
           nt = ityp(na)  
           IF (is_hubbard(nt)) THEN  
              DO is = 1, nspin  
                 DO m1 = 1, 2 * Hubbard_l(nt) + 1  
                    DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                       nb = irt (irotmq, na)  
                       DO m0 = 1, 2 * Hubbard_l(nt) + 1  
                          DO m00 = 1, 2 * Hubbard_l(nt) + 1  
                             DO jp = 1, npe
                                IF (Hubbard_l(nt).EQ.0) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + & 
                                   dnr(m0,m00,is,nb,jp) * s(ip,jp,irotmq)
                                ELSE IF (Hubbard_l(nt).EQ.1) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                   d1(m0 ,m1,irotmq) * d1(m00,m2,irotmq) *     &
                                   dnr(m0,m00,is,nb,jp) * s(ip,jp,irotmq) 
                                ELSE IF (Hubbard_l(nt).EQ.2) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                   d2(m0 ,m1,irotmq) * d2(m00,m2,irotmq) *     &
                                   dnr(m0,m00,is,nb,jp) * s(ip,jp,irotmq) 
                                ELSE IF (Hubbard_l(nt).EQ.3) THEN
                                   dns(m1,m2,is,na,ip) = dns(m1,m2,is,na,ip) + &
                                   d3(m0 ,m1,irotmq) * d3(m00,m2,irotmq) *     &
                                   dnr(m0,m00,is,nb,jp) * s(ip,jp,irotmq) 
                                ELSE
                                   CALL errore ('syme_dns', &
                                        'angular momentum not implemented', &
                                        ABS(Hubbard_l(nt)) )
                                ENDIF
                             ENDDO
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     dnr(:,:,:,:,:) = 0.5d0 * ( dnr(:,:,:,:,:) + CONJG(dns(:,:,:,:,:)) )
     dns = (0.d0, 0.d0)
  ENDIF
  !
  dnraux = (0.d0, 0.d0)
  !
  ! Symmetryze dnr -> dns
  !
  DO isym = 1, nsym
     !
     dns = (0.d0, 0.d0)
     !
     irot = isym
     !
     DO ip = 1, npe
        DO na = 1, nat  
           nt = ityp(na)  
           IF (is_hubbard(nt)) THEN  
              DO m1 = 1, 2 * Hubbard_l(nt) + 1  
                 DO m2 = 1, 2 * Hubbard_l(nt) + 1  
                    nb = irt (irot, na)  
                    DO m0 = 1, 2 * Hubbard_l(nt) + 1  
                       DO m00 = 1, 2 * Hubbard_l(nt) + 1  
                          DO jp = 1, npe
                             IF (Hubbard_l(nt).EQ.0) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) +  &
                                dnr(m0,m00,is,nb,jp) * s(ip,jp,irot) 
                             ELSE IF (Hubbard_l(nt).EQ.1) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                d1(m0 ,m1,irot) * d1(m00,m2,irot) *       &
                                dnr(m0,m00,:,nb,jp) * s(ip,jp,irot) 
                             ELSE IF (Hubbard_l(nt).EQ.2) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                d2(m0 ,m1,irot) * d2(m00,m2,irot) *       &
                                dnr(m0,m00,:,nb,jp) * s(ip,jp,irot) 
                             ELSE IF (Hubbard_l(nt).EQ.3) THEN
                                dns(m1,m2,:,na,ip) = dns(m1,m2,:,na,ip) + &
                                d3(m0 ,m1,irot) * d3(m00,m2,irot) *       &
                                dnr(m0,m00,:,nb,jp) * s(ip,jp,irot)
                             ELSE
                                CALL errore ('syme_dns', &
                                     'angular momentum not implemented', &
                                     ABS(Hubbard_l(nt)) )
                             ENDIF
                          ENDDO
                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     dnraux = dnraux + dns / nsym
  ENDDO
  dns = dnraux
  !
  DEALLOCATE (dnr)
  DEALLOCATE (dnraux)
  !
  RETURN
  !
END SUBROUTINE syme_dns
!---------------------------------------------------------------------

!---------------------------------------------------------------------
SUBROUTINE syme_dns_nc(ldim, npe, dns)
  !-------------------------------------------------------------------
  !! Symmetrize the noncollinear Hubbard occupation response to a
  !! homogeneous electric field.  The electric perturbation is a polar
  !! vector, while the Hubbard matrix is rotated jointly in orbital and
  !! spin space.
  !
  USE kinds,        ONLY : DP
  USE ions_base,    ONLY : nat, ityp
  USE ldaU,         ONLY : Hubbard_l, is_hubbard, d_spin_ldau
  USE symm_base,    ONLY : d1, d2, d3, nsym, irt, s, t_rev
  USE hubbard_nc_response, ONLY : hubbard_rotate_nc, hubbard_adjoint_nc, &
                                  hubbard_spin_inverse_nc
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ldim, npe
  COMPLEX(DP), INTENT(INOUT) :: dns(ldim,ldim,4,nat,npe)
  !
  INTEGER :: isym, ip, jp, na, nb, nt, ldim_nt
  REAL(DP) :: dorb(ldim,ldim)
  COMPLEX(DP) :: dspin_inv(2,2)
  COMPLEX(DP) :: block_in(ldim,ldim,4), block_out(ldim,ldim,4)
  COMPLEX(DP), ALLOCATABLE :: dnr(:,:,:,:,:), dnraux(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: adjoint(:,:,:,:)
  !
  ALLOCATE(dnr(ldim,ldim,4,nat,npe), dnraux(ldim,ldim,4,nat,npe))
  ALLOCATE(adjoint(ldim,ldim,4,nat))
  dnr = dns
  !
  ! A static electric perturbation produces a Hermitian occupation response.
  ! Enforce this before the unitary space-group average.
  !
  DO ip = 1, npe
     CALL hubbard_adjoint_nc(ldim, nat, dnr(:,:,:,:,ip), adjoint)
     dnr(:,:,:,:,ip) = 0.5_DP * (dnr(:,:,:,:,ip) + adjoint)
  END DO
  !
  dnraux = (0.0_DP, 0.0_DP)
  DO isym = 1, nsym
     IF (t_rev(isym) /= 0) CALL errore('syme_dns_nc', &
          'antiunitary electric-field symmetry is outside the supported scope', 1)
     DO ip = 1, npe
        DO na = 1, nat
           nt = ityp(na)
           IF (.NOT. is_hubbard(nt)) CYCLE
           nb = irt(isym,na)
           ldim_nt = 2 * Hubbard_l(nt) + 1
           dorb = 0.0_DP
           SELECT CASE(Hubbard_l(nt))
           CASE(0)
              dorb(1,1) = 1.0_DP
           CASE(1)
              dorb(1:ldim_nt,1:ldim_nt) = &
                   TRANSPOSE(d1(1:ldim_nt,1:ldim_nt,isym))
           CASE(2)
              dorb(1:ldim_nt,1:ldim_nt) = &
                   TRANSPOSE(d2(1:ldim_nt,1:ldim_nt,isym))
           CASE(3)
              dorb(1:ldim_nt,1:ldim_nt) = &
                   TRANSPOSE(d3(1:ldim_nt,1:ldim_nt,isym))
           CASE DEFAULT
              CALL errore('syme_dns_nc','angular momentum not implemented', &
                   ABS(Hubbard_l(nt)))
           END SELECT
           CALL hubbard_spin_inverse_nc(d_spin_ldau(:,:,isym), dspin_inv)
           DO jp = 1, npe
              block_in = (0.0_DP, 0.0_DP)
              block_in(1:ldim_nt,1:ldim_nt,:) = &
                   dnr(1:ldim_nt,1:ldim_nt,:,nb,jp)
              CALL hubbard_rotate_nc(ldim_nt, &
                   dorb(1:ldim_nt,1:ldim_nt), dspin_inv, &
                   block_in(1:ldim_nt,1:ldim_nt,:), &
                   block_out(1:ldim_nt,1:ldim_nt,:))
              dnraux(1:ldim_nt,1:ldim_nt,:,na,ip) = &
                   dnraux(1:ldim_nt,1:ldim_nt,:,na,ip) + &
                   s(ip,jp,isym) * block_out(1:ldim_nt,1:ldim_nt,:) / &
                   REAL(nsym,DP)
           END DO
        END DO
     END DO
  END DO
  dns = dnraux
  !
  DEALLOCATE(dnr, dnraux, adjoint)
  !
END SUBROUTINE syme_dns_nc
!---------------------------------------------------------------------
