PROGRAM hubbard_magnetic_probe
  USE kinds
  USE ldaU
  USE ions_base
  USE symm_base
  USE lr_symm_base
  USE modes
  USE qpoint
  IMPLICIT NONE
  INTEGER :: which, npe, nops, iminus, op, ia, ip, is, m1, m2
  COMPLEX(DP), ALLOCATABLE :: dns(:,:,:,:,:)
  READ(*,*) which, npe, nsymq, nops, iminus
  nsym = nsymq
  minus_q = iminus /= 0
  irotmq = iminus
  READ(*,*) xq
  d1 = 0.0_DP
  d2 = 0.0_DP
  d3 = 0.0_DP
  d_spin_ldau = (0.0_DP,0.0_DP)
  t = (0.0_DP,0.0_DP)
  tmq = (0.0_DP,0.0_DP)
  rtau = 0.0_DP
  s = 0
  t_rev = 0
  irt = 0
  DO op = 1, nops
     READ(*,*) t_rev(op), irt(op,:)
     READ(*,*) d1(:,:,op)
     READ(*,*) d_spin_ldau(:,:,op)
     READ(*,*) rtau(:,op,:)
     READ(*,*) s(:,:,op)
     READ(*,*) t(1:npe,1:npe,op,1)
  END DO
  READ(*,*) tmq(1:npe,1:npe,1)
  ALLOCATE(dns(3,3,4,nat,npe))
  READ(*,*) dns
  IF (which == 1) THEN
     CALL sym_dns_nc(3,npe,1,dns)
  ELSE
     CALL syme_dns_nc(3,npe,dns)
  ENDIF
  DO ip = 1, npe
     DO ia = 1, nat
        DO is = 1, 4
           DO m2 = 1, 3
              DO m1 = 1, 3
                 WRITE(*,'(2ES27.17E3)') REAL(dns(m1,m2,is,ia,ip),DP), AIMAG(dns(m1,m2,is,ia,ip))
              END DO
           END DO
        END DO
     END DO
  END DO
END PROGRAM hubbard_magnetic_probe
