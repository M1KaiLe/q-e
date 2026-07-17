PROGRAM test_hubbard_nc_response
  USE kinds, ONLY : DP
  USE hubbard_nc_response
  IMPLICIT NONE

  INTEGER, PARAMETER :: ldim = 2, nat = 1
  REAL(DP), PARAMETER :: tol = 1.0E-12_DP
  INTEGER :: failures, is, is1, is2, m1, m2
  INTEGER :: mbra, mket, sbra, sket, ksign
  REAL(DP) :: u_atom(nat), d1(ldim,ldim), d2(ldim,ldim), deff(ldim,ldim)
  COMPLEX(DP) :: dns(ldim,ldim,4,nat), dv(ldim,ldim,4,nat)
  COMPLEX(DP) :: a(ldim,ldim,4,nat), b(ldim,ldim,4,nat)
  COMPLEX(DP) :: r1(ldim,ldim,4), r2(ldim,ldim,4), reff(ldim,ldim,4)
  COMPLEX(DP) :: s1(2,2), s2(2,2), seff(2,2), sinv(2,2)
  REAL(DP) :: angle

  failures = 0
  DO is1 = 1, 2
     DO is2 = 1, 2
        is = hub_spin_index(is1,is2)
        CALL check_int(is,2*(is1-1)+is2,'spin flatten',failures)
        CALL hub_spin_pair(is,m1,m2)
        CALL check_int(m1,is1,'spin row inverse',failures)
        CALL check_int(m2,is2,'spin column inverse',failures)
        CALL check_int(hub_spin_transpose(is), &
             hub_spin_index(is2,is1),'spin transpose',failures)
     END DO
  END DO

  CALL hubbard_branch_indices_nc(1,1,2,1,2,mbra,mket,sbra,sket)
  CALL check_int(mbra,1,'direct branch bra orbital',failures)
  CALL check_int(mket,2,'direct branch ket orbital',failures)
  CALL check_int(sbra,1,'direct branch bra spin',failures)
  CALL check_int(sket,2,'direct branch ket spin',failures)
  CALL hubbard_branch_indices_nc(2,1,2,1,2,mbra,mket,sbra,sket)
  CALL check_int(mbra,2,'time-reversed branch bra orbital',failures)
  CALL check_int(mket,1,'time-reversed branch ket orbital',failures)
  CALL check_int(sbra,2,'time-reversed branch bra spin',failures)
  CALL check_int(sket,1,'time-reversed branch ket spin',failures)

  CALL hubbard_kramers_indices_nc(1,2,1,2,mbra,mket,sbra,sket,ksign)
  CALL check_int(mbra,2,'Kramers bra orbital',failures)
  CALL check_int(mket,1,'Kramers ket orbital',failures)
  CALL check_int(sbra,1,'Kramers bra spin',failures)
  CALL check_int(sket,2,'Kramers ket spin',failures)
  CALL check_int(ksign,-1,'Kramers off-diagonal sign',failures)

  DO is = 1, 4
     DO m2 = 1, ldim
        DO m1 = 1, ldim
           dns(m1,m2,is,1) = CMPLX(100*is+10*m2+m1, &
                7*is-3*m1+m2,kind=DP)
        END DO
     END DO
  END DO
  u_atom = 2.75_DP

  CALL hubbard_kramers_partner_nc(ldim,nat,dns,a)
  b(:,:,hub_spin_uu,1) =  TRANSPOSE(dns(:,:,hub_spin_dd,1))
  b(:,:,hub_spin_ud,1) = -TRANSPOSE(dns(:,:,hub_spin_ud,1))
  b(:,:,hub_spin_du,1) = -TRANSPOSE(dns(:,:,hub_spin_du,1))
  b(:,:,hub_spin_dd,1) =  TRANSPOSE(dns(:,:,hub_spin_uu,1))
  CALL check_array(a,b,'Kramers fixed-q partner',failures)
  CALL hubbard_kramers_partner_nc(ldim,nat,a,b)
  CALL check_array(b,dns,'Kramers partner squared',failures)

  CALL hubbard_dv_from_dns_nc(ldim,nat,u_atom,dns,dv)
  DO is = 1, 4
     DO m2 = 1, ldim
        DO m1 = 1, ldim
           CALL check_complex(dv(m1,m2,is,1), &
                -u_atom(1)*dns(m2,m1,hub_spin_transpose(is),1), &
                'Dudarev response',failures)
        END DO
     END DO
  END DO

  a = dns
  CALL hubbard_reverse_magnetization_nc(ldim,nat,a)
  b(:,:,hub_spin_uu,:) =  dns(:,:,hub_spin_dd,:)
  b(:,:,hub_spin_ud,:) = -dns(:,:,hub_spin_ud,:)
  b(:,:,hub_spin_du,:) = -dns(:,:,hub_spin_du,:)
  b(:,:,hub_spin_dd,:) =  dns(:,:,hub_spin_uu,:)
  CALL check_array(a,b,'magnetization reversal blocks',failures)
  CALL hubbard_reverse_magnetization_nc(ldim,nat,a)
  CALL check_array(a,dns,'magnetization reversal squared',failures)

  CALL hubbard_time_reverse_nc(ldim,nat,dns,a)
  b(:,:,hub_spin_uu,:) =  CONJG(dns(:,:,hub_spin_dd,:))
  b(:,:,hub_spin_ud,:) = -CONJG(dns(:,:,hub_spin_du,:))
  b(:,:,hub_spin_du,:) = -CONJG(dns(:,:,hub_spin_ud,:))
  b(:,:,hub_spin_dd,:) =  CONJG(dns(:,:,hub_spin_uu,:))
  CALL check_array(a,b,'time reversal blocks',failures)
  b = a
  CALL hubbard_time_reverse_inplace_nc(ldim,nat,b)
  CALL check_array(b,dns,'in-place time reversal squared',failures)
  CALL hubbard_time_reverse_nc(ldim,nat,a,b)
  CALL check_array(b,dns,'time reversal squared',failures)

  CALL hubbard_time_reverse_nc(ldim,nat,dns,a)
  CALL hubbard_dv_from_dns_nc(ldim,nat,u_atom,a,b)
  CALL hubbard_time_reverse_nc(ldim,nat,dv,a)
  CALL check_array(b,a,'Dudarev time-reversal covariance',failures)

  CALL hubbard_adjoint_nc(ldim,nat,dns,a)
  DO is = 1, 4
     DO m2 = 1, ldim
        DO m1 = 1, ldim
           CALL check_complex(a(m1,m2,is,1), &
                CONJG(dns(m2,m1,hub_spin_transpose(is),1)), &
                'adjoint elements',failures)
        END DO
     END DO
  END DO
  CALL hubbard_adjoint_nc(ldim,nat,a,b)
  CALL check_array(b,dns,'adjoint squared',failures)

  r1 = dns(:,:,:,1)
  d1 = 0.0_DP
  d1(1,1) = 1.0_DP
  d1(2,2) = 1.0_DP
  s1 = (0.0_DP,0.0_DP)
  s1(1,1) = (1.0_DP,0.0_DP)
  s1(2,2) = (1.0_DP,0.0_DP)
  CALL hubbard_rotate_nc(ldim,d1,s1,r1,r2)
  CALL check_block(r2,r1,'identity rotation',failures)

  angle = 0.37_DP
  d1 = RESHAPE([COS(angle),SIN(angle),-SIN(angle),COS(angle)], &
       SHAPE(d1))
  d2 = RESHAPE([COS(0.21_DP),SIN(0.21_DP), &
       -SIN(0.21_DP),COS(0.21_DP)],SHAPE(d2))
  s1 = (0.0_DP,0.0_DP)
  s1(1,1) = CMPLX(COS(0.13_DP),0.0_DP,kind=DP)
  s1(1,2) = CMPLX(0.0_DP,SIN(0.13_DP),kind=DP)
  s1(2,1) = CMPLX(0.0_DP,SIN(0.13_DP),kind=DP)
  s1(2,2) = CMPLX(COS(0.13_DP),0.0_DP,kind=DP)
  CALL hubbard_spin_inverse_nc(s1,sinv)
  seff = MATMUL(sinv,s1)
  CALL check_complex(seff(1,1),(1.0_DP,0.0_DP),'spin inverse 11',failures)
  CALL check_complex(seff(1,2),(0.0_DP,0.0_DP),'spin inverse 12',failures)
  CALL check_complex(seff(2,1),(0.0_DP,0.0_DP),'spin inverse 21',failures)
  CALL check_complex(seff(2,2),(1.0_DP,0.0_DP),'spin inverse 22',failures)
  CALL hubbard_rotate_nc(ldim,d1,s1,r1,r2)
  CALL hubbard_rotate_nc(ldim,TRANSPOSE(d1),sinv,r2,reff)
  CALL check_block(reff,r1,'inverse orbital-spin rotation',failures)
  s2 = (0.0_DP,0.0_DP)
  s2(1,1) = CMPLX(COS(0.29_DP),0.0_DP,kind=DP)
  s2(1,2) = CMPLX(SIN(0.29_DP),0.0_DP,kind=DP)
  s2(2,1) = CMPLX(-SIN(0.29_DP),0.0_DP,kind=DP)
  s2(2,2) = CMPLX(COS(0.29_DP),0.0_DP,kind=DP)
  CALL hubbard_rotate_nc(ldim,d1,s1,r1,r2)
  CALL hubbard_rotate_nc(ldim,d2,s2,r2,reff)
  deff = MATMUL(d2,d1)
  seff = MATMUL(s2,s1)
  CALL hubbard_rotate_nc(ldim,deff,seff,r1,r2)
  CALL check_block(reff,r2,'combined rotation group product',failures)

  IF (TRIM(hubbard_nc_format) /= 'PERTURBO_DFPTU_NC_V1') THEN
     WRITE(*,'(A)') 'FAIL response format marker'
     failures = failures + 1
  ENDIF
  IF (failures /= 0) ERROR STOP 1
  WRITE(*,'(A)') 'PASS hubbard_nc_response'

CONTAINS

  SUBROUTINE check_int(actual,expected,label,nfail)
    INTEGER, INTENT(IN) :: actual, expected
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nfail
    IF (actual /= expected) THEN
       WRITE(*,'(A,2(1X,I0))') 'FAIL '//TRIM(label),actual,expected
       nfail = nfail + 1
    ENDIF
  END SUBROUTINE check_int

  SUBROUTINE check_complex(actual,expected,label,nfail)
    COMPLEX(DP), INTENT(IN) :: actual, expected
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nfail
    IF (ABS(actual-expected) > tol) THEN
       WRITE(*,'(A,1X,ES12.4)') 'FAIL '//TRIM(label),ABS(actual-expected)
       nfail = nfail + 1
    ENDIF
  END SUBROUTINE check_complex

  SUBROUTINE check_array(actual,expected,label,nfail)
    COMPLEX(DP), INTENT(IN) :: actual(ldim,ldim,4,nat)
    COMPLEX(DP), INTENT(IN) :: expected(ldim,ldim,4,nat)
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nfail
    IF (MAXVAL(ABS(actual-expected)) > tol) THEN
       WRITE(*,'(A,1X,ES12.4)') 'FAIL '//TRIM(label), &
            MAXVAL(ABS(actual-expected))
       nfail = nfail + 1
    ENDIF
  END SUBROUTINE check_array

  SUBROUTINE check_block(actual,expected,label,nfail)
    COMPLEX(DP), INTENT(IN) :: actual(ldim,ldim,4)
    COMPLEX(DP), INTENT(IN) :: expected(ldim,ldim,4)
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nfail
    IF (MAXVAL(ABS(actual-expected)) > tol) THEN
       WRITE(*,'(A,1X,ES12.4)') 'FAIL '//TRIM(label), &
            MAXVAL(ABS(actual-expected))
       nfail = nfail + 1
    ENDIF
  END SUBROUTINE check_block

END PROGRAM test_hubbard_nc_response
