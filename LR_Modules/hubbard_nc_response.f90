!
! Copyright (C) 2026 Quantum ESPRESSO/Perturbo developers
! This file is distributed under the terms of the GNU General Public License.
!
MODULE hubbard_nc_response
  !
  !! Small, state-free kernels shared by PH and qe2pert for noncollinear
  !! Dudarev DFT+U response matrices.  The flattened spin order is the one
  !! used by PW: (up,up), (up,down), (down,up), (down,down).
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  PRIVATE
  !
  INTEGER, PARAMETER, PUBLIC :: hub_spin_uu = 1
  INTEGER, PARAMETER, PUBLIC :: hub_spin_ud = 2
  INTEGER, PARAMETER, PUBLIC :: hub_spin_du = 3
  INTEGER, PARAMETER, PUBLIC :: hub_spin_dd = 4
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: hubbard_nc_format = 'PERTURBO_DFPTU_NC_V1'
  !
  PUBLIC :: hub_spin_index, hub_spin_pair, hub_spin_transpose
  PUBLIC :: hubbard_branch_indices_nc
  PUBLIC :: hubbard_kramers_indices_nc
  PUBLIC :: hubbard_kramers_partner_nc
  PUBLIC :: hubbard_response_branch_inplace_nc
  PUBLIC :: hubbard_dv_from_dns_nc
  PUBLIC :: hubbard_dns_from_dv_nc
  PUBLIC :: hubbard_energy_pair_nc
  PUBLIC :: hubbard_qpair_inner_nc
  PUBLIC :: hubbard_reverse_magnetization_nc
  PUBLIC :: hubbard_time_reverse_nc
  PUBLIC :: hubbard_time_reverse_inplace_nc
  PUBLIC :: hubbard_adjoint_nc
  PUBLIC :: hubbard_spin_inverse_nc
  PUBLIC :: hubbard_rotate_nc
  PUBLIC :: hubbard_nc_diag_matrix, hubbard_nc_diag_dyn
  PUBLIC :: hubbard_nc_diag_response5, hubbard_nc_diag_response6
  PUBLIC :: hubbard_nc_diag_qmap
  !
CONTAINS
  !
  PURE INTEGER FUNCTION hub_spin_index(is1, is2) RESULT(is)
    INTEGER, INTENT(IN) :: is1, is2
    is = 2 * (is1 - 1) + is2
  END FUNCTION hub_spin_index
  !
  PURE SUBROUTINE hub_spin_pair(is, is1, is2)
    INTEGER, INTENT(IN) :: is
    INTEGER, INTENT(OUT) :: is1, is2
    is1 = (is - 1) / 2 + 1
    is2 = MOD(is - 1, 2) + 1
  END SUBROUTINE hub_spin_pair
  !
  PURE INTEGER FUNCTION hub_spin_transpose(is) RESULT(ist)
    INTEGER, INTENT(IN) :: is
    INTEGER :: is1, is2
    CALL hub_spin_pair(is, is1, is2)
    ist = hub_spin_index(is2, is1)
  END FUNCTION hub_spin_transpose
  !
  PURE COMPLEX(DP) FUNCTION hubbard_energy_pair_nc(left, right) RESULT(value)
    !! Contract two matrices in PW's occupation-array convention.  The
    !! derivative of E_U is stored so that dE_U=sum_ab v_ab*dN_ab; no
    !! orbital or spin transpose belongs in this contraction.
    COMPLEX(DP), INTENT(IN) :: left(:,:,:), right(:,:,:)
    value = SUM(left * right)
  END FUNCTION hubbard_energy_pair_nc
  !
  PURE COMPLEX(DP) FUNCTION hubbard_qpair_inner_nc(left_q, right_q) RESULT(value)
    !! Use dN_ba(-q)=conjg(dN_ab(q)) to contract a -q/q response pair in
    !! the same-index Frobenius convention.
    COMPLEX(DP), INTENT(IN) :: left_q(:,:,:), right_q(:,:,:)
    value = SUM(CONJG(left_q) * right_q)
  END FUNCTION hubbard_qpair_inner_nc
  !
  PURE SUBROUTINE hubbard_branch_indices_nc(branch, m1, m2, is1, is2, &
                                             mbra, mket, sbra, sket, sign)
    !! Projector indices for the magnetic two-Sternheimer construction.
    !! apply_trev constructs the auxiliary time-reversed wavefunctions, while
    !! the physical density response still combines their spin components as
    !! J R^T J^dagger, consistently with incdrhoscf_nc.
    INTEGER, INTENT(IN) :: branch, m1, m2, is1, is2
    INTEGER, INTENT(OUT) :: mbra, mket, sbra, sket, sign
    IF (branch == 1) THEN
       mbra = m1
       mket = m2
       sbra = is1
       sket = is2
       sign = 1
    ELSE
       CALL hubbard_kramers_indices_nc(m1, m2, is1, is2, &
                                       mbra, mket, sbra, sket, sign)
    ENDIF
  END SUBROUTINE hubbard_branch_indices_nc
  !
  PURE SUBROUTINE hubbard_kramers_indices_nc(m1, m2, is1, is2, &
                                              mbra, mket, sbra, sket, sign)
    !! Indices for a spinor time-reversed transpose.  Both the nonmagnetic
    !! fixed-q bra partner and the magnetic -B branch require J R^T J^dagger,
    !! J=i*sigma_y, rather than an ordinary combined transpose.
    INTEGER, INTENT(IN) :: m1, m2, is1, is2
    INTEGER, INTENT(OUT) :: mbra, mket, sbra, sket, sign
    !
    mbra = m2
    mket = m1
    sbra = 3 - is2
    sket = 3 - is1
    sign = MERGE(1, -1, is1 == is2)
  END SUBROUTINE hubbard_kramers_indices_nc
  !
  SUBROUTINE hubbard_kramers_partner_nc(ldim, nat, matrix_in, matrix_out)
    !! Apply J A^T J^dagger to a combined orbital-spin response matrix.
    INTEGER, INTENT(IN) :: ldim, nat
    COMPLEX(DP), INTENT(IN) :: matrix_in(ldim, ldim, 4, nat)
    COMPLEX(DP), INTENT(OUT) :: matrix_out(ldim, ldim, 4, nat)
    INTEGER :: na
    !
    DO na = 1, nat
       matrix_out(:,:,hub_spin_uu,na) = &
            TRANSPOSE(matrix_in(:,:,hub_spin_dd,na))
       matrix_out(:,:,hub_spin_ud,na) = &
            -TRANSPOSE(matrix_in(:,:,hub_spin_ud,na))
       matrix_out(:,:,hub_spin_du,na) = &
            -TRANSPOSE(matrix_in(:,:,hub_spin_du,na))
       matrix_out(:,:,hub_spin_dd,na) = &
            TRANSPOSE(matrix_in(:,:,hub_spin_uu,na))
    END DO
  END SUBROUTINE hubbard_kramers_partner_nc
  !
  SUBROUTINE hubbard_response_branch_inplace_nc(ldim, nat, matrix)
    !! Transform a finite-q response operator for the magnetic -B branch.
    !! Since dV(-q)=dV(q)^dagger, Theta dV(-q) Theta^-1 is
    !! J dV(q)^T J^dagger, without complex conjugation.
    INTEGER, INTENT(IN) :: ldim, nat
    COMPLEX(DP), INTENT(INOUT) :: matrix(ldim, ldim, 4, nat)
    COMPLEX(DP) :: tmp(ldim, ldim, 4, nat)
    !
    CALL hubbard_kramers_partner_nc(ldim, nat, matrix, tmp)
    matrix = tmp
  END SUBROUTINE hubbard_response_branch_inplace_nc
  !
  SUBROUTINE hubbard_dv_from_dns_nc(ldim, nat, u_atom, dns, dvhub)
    !! Linearized Dudarev potential: dV_ab = -U dN_ba.
    INTEGER, INTENT(IN) :: ldim, nat
    REAL(DP), INTENT(IN) :: u_atom(nat)
    COMPLEX(DP), INTENT(IN) :: dns(ldim, ldim, 4, nat)
    COMPLEX(DP), INTENT(OUT) :: dvhub(ldim, ldim, 4, nat)
    INTEGER :: na, m1, m2, is
    !
    dvhub = (0.0_DP, 0.0_DP)
    DO na = 1, nat
       DO is = 1, 4
          DO m2 = 1, ldim
             DO m1 = 1, ldim
                dvhub(m1,m2,is,na) = -u_atom(na) * &
                     dns(m2,m1,hub_spin_transpose(is),na)
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE hubbard_dv_from_dns_nc
  !
  SUBROUTINE hubbard_dns_from_dv_nc(ldim, nat, u_atom, dvhub, dns)
    !! Inverse Dudarev map used after mixing dV_U with the local response
    !! potential. Atoms with U=0 remain zero and are not active mixer entries.
    INTEGER, INTENT(IN) :: ldim, nat
    REAL(DP), INTENT(IN) :: u_atom(nat)
    COMPLEX(DP), INTENT(IN) :: dvhub(ldim, ldim, 4, nat)
    COMPLEX(DP), INTENT(OUT) :: dns(ldim, ldim, 4, nat)
    INTEGER :: na, m1, m2, is
    !
    dns = (0.0_DP, 0.0_DP)
    DO na = 1, nat
       IF (ABS(u_atom(na)) <= TINY(1.0_DP)) CYCLE
       DO is = 1, 4
          DO m2 = 1, ldim
             DO m1 = 1, ldim
                dns(m1,m2,is,na) = -dvhub(m2,m1, &
                     hub_spin_transpose(is),na) / u_atom(na)
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE hubbard_dns_from_dv_nc
  !
  SUBROUTINE hubbard_reverse_magnetization_nc(ldim, nat, matrix)
    !! Reverse the three magnetic components while preserving charge.
    !! This is the transformation used for the -B Sternheimer branch.
    INTEGER, INTENT(IN) :: ldim, nat
    COMPLEX(DP), INTENT(INOUT) :: matrix(ldim, ldim, 4, nat)
    COMPLEX(DP) :: tmp(ldim, ldim, 4)
    INTEGER :: na
    !
    DO na = 1, nat
       tmp = matrix(:,:,:,na)
       matrix(:,:,hub_spin_uu,na) =  tmp(:,:,hub_spin_dd)
       matrix(:,:,hub_spin_ud,na) = -tmp(:,:,hub_spin_ud)
       matrix(:,:,hub_spin_du,na) = -tmp(:,:,hub_spin_du)
       matrix(:,:,hub_spin_dd,na) =  tmp(:,:,hub_spin_uu)
    END DO
  END SUBROUTINE hubbard_reverse_magnetization_nc
  !
  SUBROUTINE hubbard_time_reverse_nc(ldim, nat, matrix_in, matrix_out)
    !! Antiunitary Theta=i*sigma_y*K action on an operator.  This routine is
    !! kept separate from magnetization reversal because the magnetic PH
    !! double solve already supplies time-reversed wavefunctions.
    INTEGER, INTENT(IN) :: ldim, nat
    COMPLEX(DP), INTENT(IN) :: matrix_in(ldim, ldim, 4, nat)
    COMPLEX(DP), INTENT(OUT) :: matrix_out(ldim, ldim, 4, nat)
    !
    matrix_out(:,:,hub_spin_uu,:) =  CONJG(matrix_in(:,:,hub_spin_dd,:))
    matrix_out(:,:,hub_spin_ud,:) = -CONJG(matrix_in(:,:,hub_spin_du,:))
    matrix_out(:,:,hub_spin_du,:) = -CONJG(matrix_in(:,:,hub_spin_ud,:))
    matrix_out(:,:,hub_spin_dd,:) =  CONJG(matrix_in(:,:,hub_spin_uu,:))
  END SUBROUTINE hubbard_time_reverse_nc
  !
  SUBROUTINE hubbard_time_reverse_inplace_nc(ldim, nat, matrix)
    !! Antiunitary in-place transform for a static operator when switching
    !! the magnetic Sternheimer branch.
    INTEGER, INTENT(IN) :: ldim, nat
    COMPLEX(DP), INTENT(INOUT) :: matrix(ldim, ldim, 4, nat)
    COMPLEX(DP) :: tmp(ldim, ldim, 4, nat)
    !
    CALL hubbard_time_reverse_nc(ldim, nat, matrix, tmp)
    matrix = tmp
  END SUBROUTINE hubbard_time_reverse_inplace_nc
  !
  SUBROUTINE hubbard_adjoint_nc(ldim, nat, matrix_in, matrix_out)
    INTEGER, INTENT(IN) :: ldim, nat
    COMPLEX(DP), INTENT(IN) :: matrix_in(ldim, ldim, 4, nat)
    COMPLEX(DP), INTENT(OUT) :: matrix_out(ldim, ldim, 4, nat)
    INTEGER :: na, m1, m2, is
    !
    DO na = 1, nat
       DO is = 1, 4
          DO m2 = 1, ldim
             DO m1 = 1, ldim
                matrix_out(m1,m2,is,na) = &
                     CONJG(matrix_in(m2,m1,hub_spin_transpose(is),na))
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE hubbard_adjoint_nc
  !
  PURE SUBROUTINE hubbard_spin_inverse_nc(dspin, dspin_inv)
    !! Inverse of a unitary SU(2) rotation.  Keeping this operation next to
    !! the orbital-spin kernel makes the group-operation direction explicit.
    COMPLEX(DP), INTENT(IN) :: dspin(2,2)
    COMPLEX(DP), INTENT(OUT) :: dspin_inv(2,2)
    !
    dspin_inv = TRANSPOSE(CONJG(dspin))
  END SUBROUTINE hubbard_spin_inverse_nc
  !
  SUBROUTINE hubbard_rotate_nc(ldim, dorb, dspin, matrix_in, matrix_out)
    !! Apply the same combined orbital/SU(2) convention used by new_ns_nc.
    INTEGER, INTENT(IN) :: ldim
    REAL(DP), INTENT(IN) :: dorb(ldim, ldim)
    COMPLEX(DP), INTENT(IN) :: dspin(2,2)
    COMPLEX(DP), INTENT(IN) :: matrix_in(ldim, ldim, 4)
    COMPLEX(DP), INTENT(OUT) :: matrix_out(ldim, ldim, 4)
    INTEGER :: m1, m2, m3, m4, is1, is2, is3, is4, is, js
    !
    matrix_out = (0.0_DP, 0.0_DP)
    DO is1 = 1, 2
       DO is2 = 1, 2
          is = hub_spin_index(is1, is2)
          DO is3 = 1, 2
             DO is4 = 1, 2
                js = hub_spin_index(is3, is4)
                DO m1 = 1, ldim
                   DO m2 = 1, ldim
                      DO m3 = 1, ldim
                         DO m4 = 1, ldim
                            matrix_out(m1,m2,is) = matrix_out(m1,m2,is) + &
                                 CONJG(dspin(is1,is3)) * dorb(m1,m3) * &
                                 matrix_in(m3,m4,js) * dspin(is2,is4) * &
                                 dorb(m2,m4)
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE hubbard_rotate_nc
  !
  SUBROUTINE hubbard_nc_diag_matrix(tag, iq, branch, matrix)
    !! Print compact, ionode-only diagnostics for a 4-block Hubbard matrix.
    !! The output is intentionally a single machine-readable line per call;
    !! full matrices remain in the versioned response files.
    USE io_global, ONLY : ionode, stdout
    CHARACTER(*), INTENT(IN) :: tag
    INTEGER, INTENT(IN) :: iq, branch
    COMPLEX(DP), INTENT(IN) :: matrix(:,:,:,:)
    REAL(DP) :: norm_total, norm_max, norm_block(4)
    REAL(DP) :: norm_charge, norm_mz, norm_mx, norm_my
    INTEGER :: is
    !
    IF (.NOT. ionode) RETURN
    DO is = 1, 4
       norm_block(is) = SQRT(SUM(ABS(matrix(:,:,is,:))**2))
    END DO
    norm_total = SQRT(SUM(ABS(matrix)**2))
    norm_max = MAXVAL(ABS(matrix))
    norm_charge = SQRT(SUM(ABS(matrix(:,:,hub_spin_uu,:) + &
                               matrix(:,:,hub_spin_dd,:))**2))
    norm_mz = SQRT(SUM(ABS(matrix(:,:,hub_spin_uu,:) - &
                           matrix(:,:,hub_spin_dd,:))**2))
    norm_mx = SQRT(SUM(ABS(matrix(:,:,hub_spin_ud,:) + &
                           matrix(:,:,hub_spin_du,:))**2))
    norm_my = SQRT(SUM(ABS(matrix(:,:,hub_spin_du,:) - &
                           matrix(:,:,hub_spin_ud,:))**2))
    WRITE(stdout,'(5x,"DFPTU_NC_MATRIX tag=",a," iq=",i4," branch=",i2, &
         &" norm=",es16.8," max=",es16.8)') TRIM(tag), iq, branch, &
         norm_total, norm_max
    WRITE(stdout,'(5x,"DFPTU_NC_BLOCKS iq=",i4," branch=",i2, &
         &" uu=",es16.8," ud=",es16.8," du=",es16.8," dd=",es16.8, &
         &" charge=",es16.8," mz=",es16.8," mx=",es16.8," my=",es16.8)') &
         iq, branch, norm_block(1), norm_block(2), norm_block(3), norm_block(4), &
         norm_charge, norm_mz, norm_mx, norm_my
  END SUBROUTINE hubbard_nc_diag_matrix
  !
  SUBROUTINE hubbard_nc_diag_response5(tag, iq, branch, matrix)
    !! Diagnostic for dns_branch(ldim,ldim,4,nat,npert).
    USE io_global, ONLY : ionode, stdout
    CHARACTER(*), INTENT(IN) :: tag
    INTEGER, INTENT(IN) :: iq, branch
    COMPLEX(DP), INTENT(IN) :: matrix(:,:,:,:,:)
    REAL(DP) :: norm_total, norm_max, norm_block(4)
    INTEGER :: is
    !
    IF (.NOT. ionode) RETURN
    DO is = 1, 4
       norm_block(is) = SQRT(SUM(ABS(matrix(:,:,is,:,:))**2))
    END DO
    norm_total = SQRT(SUM(ABS(matrix)**2))
    norm_max = MAXVAL(ABS(matrix))
    WRITE(stdout,'(5x,"DFPTU_NC_RESPONSE tag=",a," iq=",i4," branch=",i2, &
         &" norm=",es16.8," max=",es16.8," uu=",es16.8," ud=",es16.8, &
         &" du=",es16.8," dd=",es16.8)') TRIM(tag), iq, branch, norm_total, &
         norm_max, norm_block(1), norm_block(2), norm_block(3), norm_block(4)
  END SUBROUTINE hubbard_nc_diag_response5
  !
  SUBROUTINE hubbard_nc_diag_response6(tag, iq, branch, matrix)
    !! Diagnostic for dns_branch(ldim,ldim,4,nat,3,nat).
    USE io_global, ONLY : ionode, stdout
    CHARACTER(*), INTENT(IN) :: tag
    INTEGER, INTENT(IN) :: iq, branch
    COMPLEX(DP), INTENT(IN) :: matrix(:,:,:,:,:,:)
    REAL(DP) :: norm_total, norm_max, norm_block(4)
    INTEGER :: is
    !
    IF (.NOT. ionode) RETURN
    DO is = 1, 4
       norm_block(is) = SQRT(SUM(ABS(matrix(:,:,is,:,:,:))**2))
    END DO
    norm_total = SQRT(SUM(ABS(matrix)**2))
    norm_max = MAXVAL(ABS(matrix))
    WRITE(stdout,'(5x,"DFPTU_NC_RESPONSE tag=",a," iq=",i4," branch=",i2, &
         &" norm=",es16.8," max=",es16.8," uu=",es16.8," ud=",es16.8, &
         &" du=",es16.8," dd=",es16.8)') TRIM(tag), iq, branch, norm_total, &
         norm_max, norm_block(1), norm_block(2), norm_block(3), norm_block(4)
  END SUBROUTINE hubbard_nc_diag_response6
  !
  SUBROUTINE hubbard_nc_diag_dyn(tag, iq, matrix)
    !! Print norm and Hermitian defect of a phonon dynamical matrix.
    USE io_global, ONLY : ionode, stdout
    CHARACTER(*), INTENT(IN) :: tag
    INTEGER, INTENT(IN) :: iq
    COMPLEX(DP), INTENT(IN) :: matrix(:,:)
    REAL(DP) :: norm_total, norm_max, norm_herm, rel_herm
    !
    IF (.NOT. ionode) RETURN
    norm_total = SQRT(SUM(ABS(matrix)**2))
    norm_max = MAXVAL(ABS(matrix))
    norm_herm = SQRT(SUM(ABS(matrix - CONJG(TRANSPOSE(matrix)))**2))
    rel_herm = norm_herm / MAX(1.0_DP, norm_total)
    WRITE(stdout,'(5x,"DFPTU_NC_DYN tag=",a," iq=",i4," norm=",es16.8, &
         &" max=",es16.8," herm_defect=",es16.8," rel_herm=",es16.8)') &
         TRIM(tag), iq, norm_total, norm_max, norm_herm, rel_herm
  END SUBROUTINE hubbard_nc_diag_dyn
  !
  SUBROUTINE hubbard_nc_diag_qmap(iq, ikks_in, ikqs_in, ikmks_in)
    !! Print a deterministic checksum of the record maps used by the active
    !! Sternheimer branches.  The time-reversed map is absent for nonmagnetic
    !! noncollinear calculations.
    USE io_global, ONLY : ionode, stdout
    INTEGER, INTENT(IN) :: iq
    INTEGER, INTENT(IN) :: ikks_in(:), ikqs_in(:)
    INTEGER, INTENT(IN), OPTIONAL :: ikmks_in(:)
    INTEGER :: ik, bad
    INTEGER(KIND=8) :: hash
    !
    IF (.NOT. ionode) RETURN
    bad = 0
    hash = 0
    IF (PRESENT(ikmks_in)) THEN
       DO ik = 1, SIZE(ikks_in)
          IF (ikks_in(ik) < 1 .OR. ikqs_in(ik) < 1 .OR. ikmks_in(ik) < 1) bad = bad + 1
          hash = MOD(hash + 17_8*ikks_in(ik) + 31_8*ikqs_in(ik) + &
               47_8*ikmks_in(ik), 2147483629_8)
       END DO
       WRITE(stdout,'(5x,"DFPTU_NC_QMAP iq=",i4," nksq=",i6," bad=",i4, &
            &" hash=",i12," first=",3i6," last=",3i6)') iq, SIZE(ikks_in), bad, hash, &
            ikks_in(1), ikqs_in(1), ikmks_in(1), &
            ikks_in(SIZE(ikks_in)), ikqs_in(SIZE(ikqs_in)), ikmks_in(SIZE(ikmks_in))
    ELSE
       DO ik = 1, SIZE(ikks_in)
          IF (ikks_in(ik) < 1 .OR. ikqs_in(ik) < 1) bad = bad + 1
          hash = MOD(hash + 17_8*ikks_in(ik) + 31_8*ikqs_in(ik), 2147483629_8)
       END DO
       WRITE(stdout,'(5x,"DFPTU_NC_QMAP_DIRECT iq=",i4," nksq=",i6," bad=",i4, &
            &" hash=",i12," first=",2i6," last=",2i6)') iq, SIZE(ikks_in), bad, hash, &
            ikks_in(1), ikqs_in(1), ikks_in(SIZE(ikks_in)), ikqs_in(SIZE(ikqs_in))
    ENDIF
  END SUBROUTINE hubbard_nc_diag_qmap
  !
END MODULE hubbard_nc_response
