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
  PUBLIC :: hubbard_dv_from_dns_nc
  PUBLIC :: hubbard_reverse_magnetization_nc
  PUBLIC :: hubbard_time_reverse_nc
  PUBLIC :: hubbard_time_reverse_inplace_nc
  PUBLIC :: hubbard_adjoint_nc
  PUBLIC :: hubbard_spin_inverse_nc
  PUBLIC :: hubbard_rotate_nc
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
  PURE SUBROUTINE hubbard_branch_indices_nc(branch, m1, m2, is1, is2, &
                                             mbra, mket, sbra, sket)
    !! Projector indices for the magnetic two-Sternheimer construction.
    !! The -B solution supplies the bra/ket-transposed half of delta N.
    INTEGER, INTENT(IN) :: branch, m1, m2, is1, is2
    INTEGER, INTENT(OUT) :: mbra, mket, sbra, sket
    IF (branch == 1) THEN
       mbra = m1
       mket = m2
       sbra = is1
       sket = is2
    ELSE
       mbra = m2
       mket = m1
       sbra = is2
       sket = is1
    ENDIF
  END SUBROUTINE hubbard_branch_indices_nc
  !
  PURE SUBROUTINE hubbard_kramers_indices_nc(m1, m2, is1, is2, &
                                              mbra, mket, sbra, sket, sign)
    !! Indices for the bra-response partner in a time-reversal-symmetric
    !! spinor calculation.  At fixed q this partner is J R^T J^dagger,
    !! J=i*sigma_y, rather than the ordinary transpose used without SOC.
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
    !! In-place wrapper used while switching the magnetic Sternheimer branch.
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
END MODULE hubbard_nc_response
