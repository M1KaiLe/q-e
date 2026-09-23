MODULE ph_nc_response_conv
  ! Physical-grid residuals for NC magnetic response, independent of storage.
  USE kinds, ONLY : DP
  USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nc_response_sums, nc_response_metric, nc_response_stop
CONTAINS
  PURE SUBROUTINE nc_response_sums(nnr,npe,physical,vin,vout,sums)
    INTEGER, INTENT(IN) :: nnr,npe
    LOGICAL, INTENT(IN) :: physical(nnr)
    COMPLEX(DP), INTENT(IN) :: vin(nnr,4,npe),vout(nnr,4,npe)
    REAL(DP), INTENT(OUT) :: sums(3,2,npe)
    INTEGER :: ir,j,c,block
    sums=0.0_DP
    DO j=1,npe
      DO c=1,4
        block=MIN(c,2)
        DO ir=1,nnr
          ! Test the mask before accessing padding, which may be uninitialized.
          IF (.NOT.physical(ir)) CYCLE
          sums(1,block,j)=sums(1,block,j)+ABS(vout(ir,c,j)-vin(ir,c,j))**2
          sums(2,block,j)=sums(2,block,j)+ABS(vin(ir,c,j))**2
          sums(3,block,j)=sums(3,block,j)+ABS(vout(ir,c,j))**2
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE

  PURE SUBROUTINE nc_response_metric(nphysical,npe,sums,error2,rms,reference,scaled2)
    INTEGER, INTENT(IN) :: npe
    REAL(DP), INTENT(IN) :: nphysical,sums(3,2,npe)
    REAL(DP), INTENT(OUT) :: error2,rms(2,npe),reference(2,npe),scaled2(2,npe)
    ! Unit reference amplitude in the code's response-potential units.
    ! Above it the tolerance is relative; below it the tolerance is absolute.
    REAL(DP), PARAMETER :: reference_floor=1.0_DP
    IF (nphysical<=0.0_DP .OR. .NOT.ALL(ieee_is_finite(sums))) THEN
      error2=HUGE(1.0_DP)
      rms=error2; reference=error2; scaled2=error2
      RETURN
    ENDIF
    rms=SQRT(sums(1,:,:)/nphysical)
    reference=SQRT(MAX(sums(2,:,:),sums(3,:,:))/nphysical)
    scaled2=(rms/MAX(reference_floor,reference))**2
    error2=MAXVAL(scaled2)
  END SUBROUTINE

  PURE LOGICAL FUNCTION nc_response_stop(error2,threshold,inner_converged)
    REAL(DP), INTENT(IN) :: error2,threshold
    LOGICAL, INTENT(IN) :: inner_converged
    nc_response_stop=ieee_is_finite(error2) .AND. ieee_is_finite(threshold) .AND. threshold>0.0_DP .AND. &
        error2>=0.0_DP .AND. error2<threshold .AND. inner_converged
  END FUNCTION
END MODULE
