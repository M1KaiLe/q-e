PROGRAM qexsd_copy_magnetic_driver
  USE kinds, ONLY : dp
  USE qes_types_module, ONLY : atomic_species_type
  USE actual_copy, ONLY : qexsd_copy_atomic_species
  IMPLICIT NONE
  TYPE(atomic_species_type) :: obj
  INTEGER :: nsp, isp, failures = 0
  CHARACTER(LEN=8) :: atm(2)
  CHARACTER(LEN=256) :: psfile(2), pseudo_dir
  REAL(dp) :: mass(2), mag(2), theta(2), phi(2), pi

  pi = ACOS(-1.0_DP)
  obj%ntyp = 2
  ALLOCATE(obj%species(2))
  obj%species(1)%name = 'Fe'
  obj%species(2)%name = 'Co'
  obj%species(1)%mass = 55.845_DP
  obj%species(2)%mass = 58.933_DP
  obj%species(1)%pseudo_file = 'Fe.upf'
  obj%species(2)%pseudo_file = 'Co.upf'
  obj%pseudo_dir_ispresent = .TRUE.
  obj%pseudo_dir = '/probe/pseudo'
  DO isp = 1, 2
     obj%species(isp)%mass_ispresent = .TRUE.
     obj%species(isp)%starting_magnetization_ispresent = .TRUE.
     obj%species(isp)%spin_teta_ispresent = .TRUE.
     obj%species(isp)%spin_phi_ispresent = .TRUE.
  END DO
  obj%species(1)%starting_magnetization = .35_DP
  obj%species(2)%starting_magnetization = -.25_DP
  obj%species(1)%spin_teta = pi/3
  obj%species(2)%spin_teta = 2*pi/3
  obj%species(1)%spin_phi = -.47_DP
  obj%species(2)%spin_phi = 2.31_DP
  CALL qexsd_copy_atomic_species(obj,nsp,atm,mass,mag,theta,phi,psfile,pseudo_dir)
  CALL check('explicit_nonzero_radians', &
       ALL(mag == [.35_DP,-.25_DP]) .AND. ALL(theta == [pi/3,2*pi/3]) .AND. &
       ALL(phi == [-.47_DP,2.31_DP]))

  ! Model a second XML read with all magnetic elements absent. The inactive
  ! object members deliberately retain unrelated values, as the reader may.
  CALL absent()
  CALL qexsd_copy_atomic_species(obj,nsp,atm,mass,mag,theta,phi)
  CALL check('second_read_missing_fields', &
       ALL(mag == 0._DP) .AND. ALL(theta == 0._DP) .AND. ALL(phi == 0._DP))

  obj%species(1)%starting_magnetization_ispresent = .TRUE.
  obj%species(1)%starting_magnetization = .17_DP
  obj%species(1)%spin_phi_ispresent = .TRUE.
  obj%species(1)%spin_phi = pi/4
  obj%species(2)%spin_teta_ispresent = .TRUE.
  obj%species(2)%spin_teta = pi/2
  mag = 91._DP; theta = 92._DP; phi = 93._DP
  CALL qexsd_copy_atomic_species(obj,nsp,atm,mass,mag,theta,phi)
  CALL check('mixed_species_presence', ALL(mag == [.17_DP,0._DP]) .AND. &
       ALL(theta == [0._DP,pi/2]) .AND. ALL(phi == [pi/4,0._DP]))

  CALL absent()
  mag = 91._DP
  CALL qexsd_copy_atomic_species(obj,nsp,atm,mass,starting_magnetization=mag)
  CALL check('only_magnetization_optional', ALL(mag == 0._DP))
  theta = 92._DP
  CALL qexsd_copy_atomic_species(obj,nsp,atm,mass,angle1=theta)
  CALL check('only_theta_optional', ALL(theta == 0._DP))
  phi = 93._DP
  CALL qexsd_copy_atomic_species(obj,nsp,atm,mass,angle2=phi)
  CALL check('only_phi_optional', ALL(phi == 0._DP))

  CALL qexsd_copy_atomic_species(obj,nsp,atm,mass)
  CALL check('all_optionals_omitted_CP_call', nsp == 2 .AND. &
       ALL(atm == ['Fe','Co']) .AND. ALL(mass == [55.845_DP,58.933_DP]))
  CALL check('nonmagnetic_outputs_preserved', &
       TRIM(psfile(1)) == 'Fe.upf' .AND. TRIM(psfile(2)) == 'Co.upf' .AND. &
       TRIM(pseudo_dir) == '/probe/pseudo')

  IF (failures /= 0) STOP 1
CONTAINS
  SUBROUTINE absent()
    INTEGER :: i
    DO i = 1, 2
       obj%species(i)%starting_magnetization_ispresent = .FALSE.
       obj%species(i)%spin_teta_ispresent = .FALSE.
       obj%species(i)%spin_phi_ispresent = .FALSE.
       obj%species(i)%starting_magnetization = 711._DP
       obj%species(i)%spin_teta = 722._DP
       obj%species(i)%spin_phi = 733._DP
    END DO
  END SUBROUTINE absent

  SUBROUTINE check(name, ok)
    CHARACTER(LEN=*), INTENT(IN) :: name
    LOGICAL, INTENT(IN) :: ok
    IF (ok) THEN
       WRITE(*,'(A)') 'PASS '//name
    ELSE
       WRITE(*,'(A)') 'FAIL '//name
       failures = failures + 1
    END IF
  END SUBROUTINE check
END PROGRAM qexsd_copy_magnetic_driver
