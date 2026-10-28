! Minimal state for compiling the actual Hubbard-response symmetrizers.
MODULE kinds
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
END MODULE kinds
MODULE constants
  USE kinds
  IMPLICIT NONE
  REAL(DP), PARAMETER :: tpi = 6.283185307179586476925286766559_DP
END MODULE constants
MODULE io_global
  IMPLICIT NONE
  INTEGER, PARAMETER :: stdout = 6
  LOGICAL, PARAMETER :: ionode = .TRUE.
END MODULE io_global
MODULE ions_base
  IMPLICIT NONE
  INTEGER :: nat = 2, ityp(2) = 1
END MODULE ions_base
MODULE ldaU
  USE kinds
  IMPLICIT NONE
  INTEGER :: Hubbard_l(1) = 1, Hubbard_lmax = 1, nwfcU = 12
  LOGICAL :: is_hubbard(1) = .TRUE.
  COMPLEX(DP) :: d_spin_ldau(2,2,48)
END MODULE ldaU
MODULE lsda_mod
  IMPLICIT NONE
  LOGICAL :: lsda = .FALSE.
  INTEGER :: nspin = 4
END MODULE lsda_mod
MODULE lr_symm_base
  USE kinds
  IMPLICIT NONE
  INTEGER :: nsymq, irgq(48), irotmq
  LOGICAL :: minus_q
  REAL(DP) :: rtau(3,48,2)
END MODULE lr_symm_base
MODULE modes
  USE kinds
  IMPLICIT NONE
  COMPLEX(DP) :: t(3,3,48,6), tmq(3,3,6)
END MODULE modes
MODULE qpoint
  USE kinds
  IMPLICIT NONE
  REAL(DP) :: xq(3)
END MODULE qpoint
MODULE uspp_param
  IMPLICIT NONE
  INTEGER :: upf(1)
END MODULE uspp_param
MODULE symm_base
  USE kinds
  IMPLICIT NONE
  INTEGER :: nsym, irt(48,2), s(3,3,48), invs(48), t_rev(48)
  REAL(DP) :: d1(3,3,48), d2(5,5,48), d3(7,7,48)
END MODULE symm_base
SUBROUTINE errore(routine, message, code)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: routine, message
  INTEGER, INTENT(IN) :: code
  WRITE(*,*) TRIM(routine), ': ', TRIM(message), code
  STOP 1
END SUBROUTINE errore
