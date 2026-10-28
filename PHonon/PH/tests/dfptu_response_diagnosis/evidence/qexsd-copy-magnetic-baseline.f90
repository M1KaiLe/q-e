  SUBROUTINE qexsd_copy_atomic_species (atomic_species, nsp, atm, amass, &
       starting_magnetization, angle1, angle2, psfile, pseudo_dir)
    !---------------------------------------------------------------------------    !
    USE qes_types_module, ONLY : atomic_species_type
    !
    IMPLICIT NONE
    !
    TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_species
    INTEGER, INTENT(out) :: nsp
    CHARACTER(LEN=*), INTENT(out) :: atm(:)
    REAL(dp), INTENT(out) :: amass(:)
    REAL(dp), OPTIONAL, INTENT(out) :: starting_magnetization(:), &
         angle1(:), angle2(:)
    CHARACTER(LEN=*), OPTIONAL, INTENT(out) :: psfile(:), pseudo_dir
    !
    INTEGER :: isp
    !
    nsp = atomic_species%ntyp
    DO isp = 1, nsp
       amass(isp) = 0.d0
       IF (atomic_species%species(isp)%mass_ispresent) &
            amass(isp) = atomic_species%species(isp)%mass
       atm(isp) = TRIM ( atomic_species%species(isp)%name )
       IF ( PRESENT (psfile) ) THEN
          psfile(isp) = TRIM ( atomic_species%species(isp)%pseudo_file)
       END IF
       IF ( PRESENT (starting_magnetization) ) THEN
          IF ( atomic_species%species(isp)%starting_magnetization_ispresent) THEN
             starting_magnetization(isp) = atomic_species%species(isp)%starting_magnetization
          END IF
       END IF
       IF ( PRESENT (angle1) ) THEN
          IF ( atomic_species%species(isp)%spin_teta_ispresent ) THEN
             angle1(isp) =  atomic_species%species(isp)%spin_teta
          END IF
       END IF
       IF ( PRESENT (angle2) ) THEN
          IF ( atomic_species%species(isp)%spin_phi_ispresent ) THEN
             angle2(isp) = atomic_species%species(isp)%spin_phi
          END IF
       END IF
    END DO
    !
    ! ... this is where PP files were originally found (if available)
    !
    IF ( PRESENT (pseudo_dir) ) THEN
       IF ( atomic_species%pseudo_dir_ispresent) THEN
          pseudo_dir = TRIM(atomic_species%pseudo_dir)
       ELSE
          pseudo_dir = ' '
       END IF
    END IF
    !
  END SUBROUTINE qexsd_copy_atomic_species
