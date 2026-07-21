PROGRAM test_hubbard_nc_response
  USE kinds, ONLY : DP
  USE hubbard_nc_response
  IMPLICIT NONE

  INTEGER, PARAMETER :: ldim = 2, nat = 1, hdim = 2*ldim
  REAL(DP), PARAMETER :: tol = 1.0E-12_DP
  INTEGER :: failures, is, is1, is2, is3, is4, js, m1, m2, m3, m4
  INTEGER :: i, j, icell
  INTEGER :: mbra, mket, sbra, sket, bsign, ksign
  REAL(DP) :: u_atom(nat), d1(ldim,ldim), d2(ldim,ldim), deff(ldim,ldim)
  COMPLEX(DP) :: dns(ldim,ldim,4,nat), dv(ldim,ldim,4,nat)
  COMPLEX(DP) :: a(ldim,ldim,4,nat), b(ldim,ldim,4,nat)
  COMPLEX(DP) :: r1(ldim,ldim,4), r2(ldim,ldim,4), reff(ldim,ldim,4)
  COMPLEX(DP) :: rref(ldim,ldim,4)
  COMPLEX(DP) :: pair_left(ldim,ldim,4), pair_right(ldim,ldim,4)
  COMPLEX(DP) :: s1(2,2), s2(2,2), seff(2,2), sinv(2,2)
  COMPLEX(DP) :: x(ldim,2), y(ldim,2), tx(ldim,2), ty(ldim,2)
  COMPLEX(DP) :: n0(ldim,ldim,4), dn(ldim,ldim,4), v0(ldim,ldim,4)
  COMPLEX(DP) :: aq(ldim,ldim,4), bq(ldim,ldim,4)
  COMPLEX(DP) :: aq_adj(ldim,ldim,4), bq_adj(ldim,ldim,4)
  COMPLEX(DP) :: ibq(ldim,ldim,4), ibq_adj(ldim,ldim,4)
  COMPLEX(DP) :: dir_a(ldim,ldim,4), dir_b(ldim,ldim,4)
  COMPLEX(DP) :: full_n0(hdim,hdim), full_dn(hdim,hdim), full_v0(hdim,hdim)
  COMPLEX(DP) :: fd_first, fd_pair, mixed_re, mixed_im, phase
  COMPLEX(DP) :: wrong_pair, wrong_qpair
  COMPLEX(DP) :: jspin(hdim,hdim), phi_k(hdim,hdim), phi_q(hdim,hdim)
  COMPLEX(DP) :: dphi_k0(hdim,hdim), dphi_q0(hdim,hdim)
  COMPLEX(DP) :: phi_k_aux(hdim,hdim), phi_q_aux(hdim,hdim)
  COMPLEX(DP) :: dphi_k_minus(hdim,hdim), dphi_q_minus(hdim,hdim)
  COMPLEX(DP) :: dphi_k_aux(hdim,hdim), dphi_q_aux(hdim,hdim)
  COMPLEX(DP) :: dphi_k_wrong(hdim,hdim), dphi_q_wrong(hdim,hdim)
  COMPLEX(DP) :: v_static(hdim,hdim), dvq0(hdim,hdim), dv_q(hdim,hdim)
  COMPLEX(DP) :: dv_minus(hdim,hdim), v_aux(hdim,hdim), dv_aux(hdim,hdim)
  COMPLEX(DP) :: v_double(hdim,hdim), dv_conj_wrong(hdim,hdim)
  COMPLEX(DP) :: dv_u_wrong(hdim,hdim)
  COMPLEX(DP) :: v_block(ldim,ldim,4,nat), work_block(ldim,ldim,4,nat)
  COMPLEX(DP) :: psi_minus(hdim), psi_aux(hdim)
  COMPLEX(DP) :: proj_minus(hdim), dproj_minus(hdim)
  COMPLEX(DP) :: proj_aux(hdim), dproj_aux(hdim), dproj_wrong(hdim)
  COMPLEX(DP) :: rhs_minus(hdim), rhs_aux(hdim), rhs_ref(hdim)
  COMPLEX(DP) :: rhs_wrong(hdim), rhs_term_minus(hdim), rhs_term_aux(hdim)
  COMPLEX(DP) :: u_mode, u_hub
  REAL(DP) :: angle, explicit_err, inverse_difference, fd_step

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

  CALL hubbard_branch_indices_nc(1,1,2,1,2,mbra,mket,sbra,sket,bsign)
  CALL check_int(mbra,1,'direct branch bra orbital',failures)
  CALL check_int(mket,2,'direct branch ket orbital',failures)
  CALL check_int(sbra,1,'direct branch bra spin',failures)
  CALL check_int(sket,2,'direct branch ket spin',failures)
  CALL check_int(bsign,1,'direct branch sign',failures)
  DO is1 = 1, 2
     DO is2 = 1, 2
        CALL hubbard_branch_indices_nc(2,1,2,is1,is2, &
             mbra,mket,sbra,sket,bsign)
        CALL check_int(mbra,2,'time-reversed branch bra orbital',failures)
        CALL check_int(mket,1,'time-reversed branch ket orbital',failures)
        CALL check_int(sbra,3-is2,'time-reversed branch bra spin',failures)
        CALL check_int(sket,3-is1,'time-reversed branch ket spin',failures)
        CALL check_int(bsign,MERGE(1,-1,is1 == is2), &
             'time-reversed branch sign',failures)
     END DO
  END DO

  CALL hubbard_kramers_indices_nc(1,2,1,2,mbra,mket,sbra,sket,ksign)
  CALL check_int(mbra,2,'Kramers bra orbital',failures)
  CALL check_int(mket,1,'Kramers ket orbital',failures)
  CALL check_int(sbra,1,'Kramers bra spin',failures)
  CALL check_int(sket,2,'Kramers ket spin',failures)
  CALL check_int(ksign,-1,'Kramers off-diagonal sign',failures)

  ! A Hermitian purely imaginary off-diagonal pair distinguishes PW's
  ! same-index energy convention from an accidental joint transpose.
  pair_left = (0.0_DP,0.0_DP)
  pair_left(1,2,hub_spin_ud) = (0.0_DP,1.0_DP)
  pair_left(2,1,hub_spin_du) = (0.0_DP,-1.0_DP)
  pair_right = pair_left
  CALL check_complex(hubbard_energy_pair_nc(pair_left,pair_right), &
       (-2.0_DP,0.0_DP),'same-index Hubbard energy pair',failures)
  CALL check_complex(hubbard_qpair_inner_nc(pair_left,pair_right), &
       (2.0_DP,0.0_DP),'same-index finite-q Frobenius pair',failures)

  ! Derive the contraction from the scalar Dudarev functional, independently
  ! of the production kernels.  PW stores dE/dN_ab at the same array index as
  ! dN_ab, even though that derivative contains the transposed N_ba.
  u_hub = CMPLX(2.25_DP,0.0_DP,kind=DP)
  fd_step = 0.125_DP
  DO j = 1, hdim
     DO i = 1, hdim
        full_n0(i,j) = CMPLX(0.031_DP*(2*i-j), &
             0.023_DP*(i+3*j),kind=DP)
        full_dn(i,j) = CMPLX(0.017_DP*(i-2*j), &
             0.029_DP*(2*i+j),kind=DP)
     END DO
  END DO
  full_n0 = 0.5_DP * (full_n0 + TRANSPOSE(CONJG(full_n0)))
  full_dn = 0.5_DP * (full_dn + TRANSPOSE(CONJG(full_dn)))
  DO i = 1, hdim
     full_n0(i,i) = full_n0(i,i) + CMPLX(0.11_DP*i,0.0_DP,kind=DP)
     full_v0(i,i) = 0.5_DP * u_hub - u_hub * full_n0(i,i)
     DO j = 1, hdim
        IF (i /= j) full_v0(i,j) = -u_hub * full_n0(j,i)
     END DO
  END DO
  CALL full_to_blocks(full_n0,n0)
  CALL full_to_blocks(full_dn,dn)
  CALL full_to_blocks(full_v0,v0)
  fd_first = (dudarev_energy(u_hub,n0+fd_step*dn) - &
       dudarev_energy(u_hub,n0-fd_step*dn)) / (2.0_DP*fd_step)
  CALL check_complex(hubbard_energy_pair_nc(v0,dn),fd_first, &
       'Dudarev finite-difference same-index derivative',failures)
  wrong_pair = joint_transpose_pair(v0,dn,.FALSE.)
  CALL check_separated(wrong_pair,fd_first, &
       'old joint-transpose energy counterexample',failures)

  ! A four-cell Fourier sum turns arbitrary complex q matrices into Hermitian
  ! real-space perturbations.  The mixed Dudarev curvature recovers both the
  ! real and imaginary parts of Tr(A^dagger B) without assuming the answer.
  DO is = 1, 4
     DO m2 = 1, ldim
        DO m1 = 1, ldim
           aq(m1,m2,is) = CMPLX(0.019_DP*(m1+2*m2+3*is), &
                0.013_DP*(2*m1-m2+is),kind=DP)
           bq(m1,m2,is) = CMPLX(0.017_DP*(3*m1-m2+2*is), &
                0.011_DP*(-m1+3*m2+is),kind=DP)
        END DO
     END DO
  END DO
  CALL independent_adjoint(aq,aq_adj)
  CALL independent_adjoint(bq,bq_adj)
  ibq = CMPLX(0.0_DP,1.0_DP,kind=DP) * bq
  CALL independent_adjoint(ibq,ibq_adj)
  mixed_re = (0.0_DP,0.0_DP)
  mixed_im = (0.0_DP,0.0_DP)
  DO icell = 0, 3
     SELECT CASE (icell)
     CASE (0)
        phase = CMPLX(1.0_DP,0.0_DP,kind=DP)
     CASE (1)
        phase = CMPLX(0.0_DP,1.0_DP,kind=DP)
     CASE (2)
        phase = CMPLX(-1.0_DP,0.0_DP,kind=DP)
     CASE DEFAULT
        phase = CMPLX(0.0_DP,-1.0_DP,kind=DP)
     END SELECT
     dir_a = phase*aq + CONJG(phase)*aq_adj
     dir_b = phase*bq + CONJG(phase)*bq_adj
     mixed_re = mixed_re + dudarev_mixed_difference( &
          u_hub,n0,dir_a,dir_b,fd_step)
     dir_b = phase*ibq + CONJG(phase)*ibq_adj
     mixed_im = mixed_im + dudarev_mixed_difference( &
          u_hub,n0,dir_a,dir_b,fd_step)
  END DO
  fd_pair = CMPLX(-REAL(mixed_re,DP)/(8.0_DP*REAL(u_hub,DP)), &
       REAL(mixed_im,DP)/(8.0_DP*REAL(u_hub,DP)),kind=DP)
  CALL check_complex(hubbard_qpair_inner_nc(aq,bq),fd_pair, &
       'Dudarev Fourier mixed-difference q pair',failures)
  wrong_qpair = joint_transpose_pair(aq,bq,.TRUE.)
  CALL check_separated(wrong_qpair,fd_pair, &
       'old joint-transpose finite-q counterexample',failures)

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

  b = 0.5_DP * (dns + a)
  CALL hubbard_kramers_partner_nc(ldim,nat,b,a)
  CALL check_array(a,b,'pure-TR ground-state average',failures)

  a = dns
  CALL hubbard_response_branch_inplace_nc(ldim,nat,a)
  CALL hubbard_kramers_partner_nc(ldim,nat,dns,b)
  CALL check_array(a,b,'finite-q response branch transpose',failures)
  CALL hubbard_time_reverse_nc(ldim,nat,dns,b)
  IF (MAXVAL(ABS(a-b)) <= tol) THEN
     WRITE(*,'(A)') 'FAIL finite-q transpose confused with complex conjugation'
     failures = failures + 1
  ENDIF
  CALL hubbard_response_branch_inplace_nc(ldim,nat,a)
  CALL check_array(a,dns,'finite-q response branch squared',failures)

  a = (0.0_DP,0.0_DP)
  DO is1 = 1, 2
     DO is2 = 1, 2
        is = hub_spin_index(is1,is2)
        DO m2 = 1, ldim
           DO m1 = 1, ldim
              CALL hubbard_branch_indices_nc(2,m1,m2,is1,is2, &
                   mbra,mket,sbra,sket,bsign)
              a(m1,m2,is,1) = REAL(bsign,DP) * &
                   dns(mbra,mket,hub_spin_index(sbra,sket),1)
           END DO
        END DO
     END DO
  END DO
  CALL hubbard_kramers_partner_nc(ldim,nat,dns,b)
  CALL check_array(a,b,'magnetic branch-2 Kramers transpose',failures)

  ! Independent wavefunction-level check.  apply_trev produces tx=J*x^*
  ! and ty=J*y^*.  Mapping their auxiliary bilinear back with J R^T J^dagger
  ! must give the missing bra response y^* x^T in the physical spin basis.
  x(:,1) = [CMPLX(0.7_DP,-0.2_DP,kind=DP), &
            CMPLX(-0.4_DP,0.9_DP,kind=DP)]
  x(:,2) = [CMPLX(0.3_DP,0.8_DP,kind=DP), &
            CMPLX(1.1_DP,-0.5_DP,kind=DP)]
  y(:,1) = [CMPLX(-0.6_DP,0.4_DP,kind=DP), &
            CMPLX(0.2_DP,1.3_DP,kind=DP)]
  y(:,2) = [CMPLX(0.9_DP,0.1_DP,kind=DP), &
            CMPLX(-0.8_DP,-0.7_DP,kind=DP)]
  tx(:,1) = -CONJG(x(:,2))
  tx(:,2) =  CONJG(x(:,1))
  ty(:,1) = -CONJG(y(:,2))
  ty(:,2) =  CONJG(y(:,1))
  DO is1 = 1, 2
     DO is2 = 1, 2
        is = hub_spin_index(is1,is2)
        DO m1 = 1, ldim
           DO m2 = 1, ldim
              dns(m1,m2,is,1) = CONJG(tx(m1,is1)) * ty(m2,is2)
              rref(m1,m2,is) = CONJG(y(m1,is1)) * x(m2,is2)
           END DO
        END DO
     END DO
  END DO
  CALL hubbard_kramers_partner_nc(ldim,nat,dns,a)
  CALL check_block(a(:,:,:,1),rref,'auxiliary -B bilinear to physical bra',failures)
  DO is = 1, 4
     reff(:,:,is) = TRANSPOSE(dns(:,:,hub_spin_transpose(is),1))
  END DO
  IF (MAXVAL(ABS(reff-rref)) <= 1.0E-8_DP) THEN
     WRITE(*,'(A)') 'FAIL ordinary transpose accidentally matches physical bra'
     failures = failures + 1
  ENDIF

  ! Full magnetic bare-RHS identity on a small complex projector problem.
  ! The physical -q RHS is transformed explicitly with Theta=J K.  The
  ! auxiliary branch is then assembled independently using the two production
  ! operator maps: J V^* J^dagger for static V and J dV(q)^T J^dagger for a
  ! finite-q response.  All three terms in d(Phi V Phi^dagger) are retained.
  jspin = (0.0_DP,0.0_DP)
  DO m1 = 1, ldim
     jspin(m1,ldim+m1) = CMPLX(-1.0_DP,0.0_DP,kind=DP)
     jspin(ldim+m1,m1) = CMPLX(1.0_DP,0.0_DP,kind=DP)
  END DO
  DO j = 1, hdim
     DO i = 1, hdim
        phi_k(i,j) = CMPLX(0.07_DP*(2*i-j), &
             0.05_DP*(i+2*j),kind=DP)
        phi_q(i,j) = CMPLX(0.04_DP*(-i+3*j), &
             0.06_DP*(2*i+j),kind=DP)
        dphi_k0(i,j) = CMPLX(0.03_DP*(i+j), &
             0.02_DP*(-2*i+j),kind=DP)
        dphi_q0(i,j) = CMPLX(0.025_DP*(3*i-j), &
             0.035_DP*(i-2*j),kind=DP)
        v_static(i,j) = CMPLX(0.09_DP*(i-2*j), &
             0.055_DP*(2*i+j),kind=DP)
        dvq0(i,j) = CMPLX(0.045_DP*(i+3*j), &
             0.065_DP*(-i+2*j),kind=DP)
     END DO
     psi_minus(j) = CMPLX(0.13_DP*(2*j-3), &
          0.08_DP*(j+1),kind=DP)
  END DO
  v_static = 0.5_DP * (v_static + TRANSPOSE(CONJG(v_static)))
  DO i = 1, hdim
     v_static(i,i) = v_static(i,i) + CMPLX(0.17_DP*i,0.0_DP,kind=DP)
  END DO
  u_mode = CMPLX(0.6_DP,0.8_DP,kind=DP)
  dv_q = u_mode * dvq0
  CALL independent_full_adjoint(dv_q,dv_minus)
  dphi_q_minus = CONJG(u_mode) * dphi_q0
  ! The k-side derivative is a bra buffer: its dagger supplies u^* in the
  ! physical -q projected derivative, just as DOT_PRODUCT does in PH.
  dphi_k_minus = u_mode * dphi_k0
  proj_minus = MATMUL(TRANSPOSE(CONJG(phi_k)),psi_minus)
  dproj_minus = MATMUL(TRANSPOSE(CONJG(dphi_k_minus)),psi_minus)

  rhs_term_minus = MATMUL(dphi_q_minus,MATMUL(v_static,proj_minus))
  rhs_minus = rhs_term_minus
  rhs_term_minus = MATMUL(phi_q,MATMUL(v_static,dproj_minus))
  rhs_minus = rhs_minus + rhs_term_minus
  rhs_term_minus = MATMUL(phi_q,MATMUL(dv_minus,proj_minus))
  rhs_minus = rhs_minus + rhs_term_minus
  CALL theta_vector(rhs_minus,jspin,rhs_ref)

  CALL theta_matrix(phi_k,jspin,phi_k_aux)
  CALL theta_matrix(phi_q,jspin,phi_q_aux)
  CALL theta_matrix(dphi_k_minus,jspin,dphi_k_aux)
  CALL theta_matrix(dphi_q_minus,jspin,dphi_q_aux)
  CALL theta_vector(psi_minus,jspin,psi_aux)
  proj_aux = MATMUL(TRANSPOSE(CONJG(phi_k_aux)),psi_aux)
  dproj_aux = MATMUL(TRANSPOSE(CONJG(dphi_k_aux)),psi_aux)
  CALL theta_vector(proj_minus,jspin,rhs_wrong)
  CALL check_vector(proj_aux,rhs_wrong,'Theta projector equivalence',failures)
  CALL theta_vector(dproj_minus,jspin,rhs_wrong)
  CALL check_vector(dproj_aux,rhs_wrong,'Theta derivative-projector equivalence',failures)

  CALL full_to_blocks(v_static,v_block(:,:,:,1))
  CALL hubbard_time_reverse_nc(ldim,nat,v_block,work_block)
  CALL blocks_to_full(work_block(:,:,:,1),v_aux)
  CALL theta_matrix(v_static,jspin,v_double)
  CALL check_matrix(v_aux,v_double,'Theta static Hubbard operator',failures)
  CALL full_to_blocks(dv_q,v_block(:,:,:,1))
  work_block = v_block
  CALL hubbard_response_branch_inplace_nc(ldim,nat,work_block)
  CALL blocks_to_full(work_block(:,:,:,1),dv_aux)
  CALL theta_matrix(dv_minus,jspin,dv_conj_wrong)
  CALL check_matrix(dv_aux,dv_conj_wrong, &
       'Theta finite-q response operator',failures)

  rhs_term_aux = MATMUL(dphi_q_aux,MATMUL(v_aux,proj_aux))
  CALL theta_vector(MATMUL(dphi_q_minus, &
       MATMUL(v_static,proj_minus)),jspin,rhs_term_minus)
  CALL check_vector(rhs_term_aux,rhs_term_minus, &
       'Theta RHS projector-derivative term',failures)
  CALL check_nonzero(rhs_term_aux,'projector-derivative RHS term',failures)
  rhs_aux = rhs_term_aux
  rhs_term_aux = MATMUL(phi_q_aux,MATMUL(v_aux,dproj_aux))
  CALL theta_vector(MATMUL(phi_q, &
       MATMUL(v_static,dproj_minus)),jspin,rhs_term_minus)
  CALL check_vector(rhs_term_aux,rhs_term_minus, &
       'Theta RHS bra-projector derivative term',failures)
  CALL check_nonzero(rhs_term_aux,'bra-projector derivative RHS term',failures)
  rhs_aux = rhs_aux + rhs_term_aux
  rhs_term_aux = MATMUL(phi_q_aux,MATMUL(dv_aux,proj_aux))
  CALL theta_vector(MATMUL(phi_q, &
       MATMUL(dv_minus,proj_minus)),jspin,rhs_term_minus)
  CALL check_vector(rhs_term_aux,rhs_term_minus, &
       'Theta RHS response-potential term',failures)
  CALL check_nonzero(rhs_term_aux,'response-potential RHS term',failures)
  rhs_aux = rhs_aux + rhs_term_aux
  CALL check_vector(rhs_aux,rhs_ref,'complete three-term Theta RHS',failures)

  ! Counterexample 1: a second static Theta restores V and must not describe
  ! the auxiliary -B branch.
  CALL full_to_blocks(v_aux,work_block(:,:,:,1))
  CALL hubbard_time_reverse_nc(ldim,nat,work_block,v_block)
  CALL blocks_to_full(v_block(:,:,:,1),v_double)
  rhs_wrong = MATMUL(dphi_q_aux,MATMUL(v_double,proj_aux)) + &
       MATMUL(phi_q_aux,MATMUL(v_double,dproj_aux)) + &
       MATMUL(phi_q_aux,MATMUL(dv_aux,proj_aux))
  CALL check_vector_separated(rhs_wrong,rhs_ref, &
       'double-Theta RHS counterexample',failures)

  ! Counterexample 2: treating dV(q) as a fixed-q operator uses complex
  ! conjugation instead of the required transpose after dV(-q)=dV(q)^dagger.
  CALL full_to_blocks(dv_q,v_block(:,:,:,1))
  CALL hubbard_time_reverse_nc(ldim,nat,v_block,work_block)
  CALL blocks_to_full(work_block(:,:,:,1),dv_conj_wrong)
  rhs_wrong = MATMUL(dphi_q_aux,MATMUL(v_aux,proj_aux)) + &
       MATMUL(phi_q_aux,MATMUL(v_aux,dproj_aux)) + &
       MATMUL(phi_q_aux,MATMUL(dv_conj_wrong,proj_aux))
  CALL check_vector_separated(rhs_wrong,rhs_ref, &
       'finite-q conjugation-vs-transpose counterexample',failures)

  ! Counterexample 3: use u^* on the auxiliary branch.  The output derivative,
  ! daggered k-side derivative buffer, and response potential are all changed
  ! consistently so that an accidental real mode cannot hide this error.
  CALL theta_matrix(dphi_q0,jspin,dphi_q_wrong)
  dphi_q_wrong = CONJG(u_mode) * dphi_q_wrong
  CALL theta_matrix(dphi_k0,jspin,dphi_k_wrong)
  dphi_k_wrong = u_mode * dphi_k_wrong
  dproj_wrong = MATMUL(TRANSPOSE(CONJG(dphi_k_wrong)),psi_aux)
  dv_u_wrong = (CONJG(u_mode)/u_mode) * dv_aux
  rhs_wrong = MATMUL(dphi_q_wrong,MATMUL(v_aux,proj_aux)) + &
       MATMUL(phi_q_aux,MATMUL(v_aux,dproj_wrong)) + &
       MATMUL(phi_q_aux,MATMUL(dv_u_wrong,proj_aux))
  CALL check_vector_separated(rhs_wrong,rhs_ref, &
       'u-vs-u-star magnetic RHS counterexample',failures)

  ! Restore the deterministic matrix used by the remaining kernel tests.
  DO is = 1, 4
     DO m2 = 1, ldim
        DO m1 = 1, ldim
           dns(m1,m2,is,1) = CMPLX(100*is+10*m2+m1, &
                7*is-3*m1+m2,kind=DP)
        END DO
     END DO
  END DO

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

  CALL hubbard_dns_from_dv_nc(ldim,nat,u_atom,dv,a)
  CALL check_array(a,dns,'Dudarev response inverse',failures)
  u_atom = 0.0_DP
  CALL hubbard_dns_from_dv_nc(ldim,nat,u_atom,dv,a)
  b = (0.0_DP,0.0_DP)
  CALL check_array(a,b,'Dudarev zero-U inverse',failures)
  u_atom = 2.75_DP

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

  CALL hubbard_adjoint_nc(ldim,nat,dns,a)
  b = 0.5_DP * (dns + a)
  CALL hubbard_adjoint_nc(ldim,nat,b,a)
  CALL check_array(a,b,'electric-field Hermitian completion',failures)

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
  ! Explicit reference for hubbard_rotate_nc's direct combined rotation.
  ! sym_dns_nc obtains the response pullback by passing TRANSPOSE(D_l) and
  ! S^dagger to this kernel.
  rref = (0.0_DP,0.0_DP)
  DO is1 = 1, 2
     DO is2 = 1, 2
        is = hub_spin_index(is1,is2)
        DO is3 = 1, 2
           DO is4 = 1, 2
              js = hub_spin_index(is3,is4)
              DO m1 = 1, ldim
                 DO m2 = 1, ldim
                    DO m3 = 1, ldim
                       DO m4 = 1, ldim
                          rref(m1,m2,is) = rref(m1,m2,is) + &
                               CONJG(s1(is1,is3)) * d1(m1,m3) * &
                               r1(m3,m4,js) * s1(is2,is4) * d1(m2,m4)
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
  END DO
  explicit_err = MAXVAL(ABS(r2-rref))
  WRITE(*,'(A,1X,ES12.4)') 'explicit direct rotation maxerr=', explicit_err
  CALL check_block(r2,rref,'explicit direct orbital-spin rotation',failures)
  CALL hubbard_rotate_nc(ldim,TRANSPOSE(d1),sinv,r1,reff)
  inverse_difference = MAXVAL(ABS(reff-rref))
  WRITE(*,'(A,1X,ES12.4)') 'inverse pullback vs direct difference=', inverse_difference
  IF (inverse_difference <= 1.0E-8_DP) THEN
     WRITE(*,'(A)') 'FAIL rotation test does not distinguish pullback direction'
     failures = failures + 1
  ENDIF
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

  PURE SUBROUTINE full_to_blocks(full,blocks)
    COMPLEX(DP), INTENT(IN) :: full(hdim,hdim)
    COMPLEX(DP), INTENT(OUT) :: blocks(ldim,ldim,4)
    INTEGER :: row, col, io1, io2, iss1, iss2, iblock
    DO iss1 = 1, 2
       DO iss2 = 1, 2
          iblock = hub_spin_index(iss1,iss2)
          DO io2 = 1, ldim
             DO io1 = 1, ldim
                row = io1 + ldim*(iss1-1)
                col = io2 + ldim*(iss2-1)
                blocks(io1,io2,iblock) = full(row,col)
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE full_to_blocks

  PURE SUBROUTINE blocks_to_full(blocks,full)
    COMPLEX(DP), INTENT(IN) :: blocks(ldim,ldim,4)
    COMPLEX(DP), INTENT(OUT) :: full(hdim,hdim)
    INTEGER :: row, col, io1, io2, iss1, iss2, iblock
    DO iss1 = 1, 2
       DO iss2 = 1, 2
          iblock = hub_spin_index(iss1,iss2)
          DO io2 = 1, ldim
             DO io1 = 1, ldim
                row = io1 + ldim*(iss1-1)
                col = io2 + ldim*(iss2-1)
                full(row,col) = blocks(io1,io2,iblock)
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE blocks_to_full

  PURE COMPLEX(DP) FUNCTION dudarev_energy(u_value,matrix) RESULT(energy)
    COMPLEX(DP), INTENT(IN) :: u_value
    COMPLEX(DP), INTENT(IN) :: matrix(ldim,ldim,4)
    COMPLEX(DP) :: full(hdim,hdim), trace_n, trace_n2
    INTEGER :: ii, jj
    CALL blocks_to_full(matrix,full)
    trace_n = (0.0_DP,0.0_DP)
    trace_n2 = (0.0_DP,0.0_DP)
    DO ii = 1, hdim
       trace_n = trace_n + full(ii,ii)
       DO jj = 1, hdim
          trace_n2 = trace_n2 + full(ii,jj)*full(jj,ii)
       END DO
    END DO
    energy = 0.5_DP*u_value*(trace_n-trace_n2)
  END FUNCTION dudarev_energy

  PURE COMPLEX(DP) FUNCTION dudarev_mixed_difference( &
       u_value,base,left_direction,right_direction,step) RESULT(value)
    COMPLEX(DP), INTENT(IN) :: u_value
    COMPLEX(DP), INTENT(IN) :: base(ldim,ldim,4)
    COMPLEX(DP), INTENT(IN) :: left_direction(ldim,ldim,4)
    COMPLEX(DP), INTENT(IN) :: right_direction(ldim,ldim,4)
    REAL(DP), INTENT(IN) :: step
    value = (dudarev_energy(u_value,base+step*left_direction+ &
             step*right_direction) - &
         dudarev_energy(u_value,base+step*left_direction- &
             step*right_direction) - &
         dudarev_energy(u_value,base-step*left_direction+ &
             step*right_direction) + &
         dudarev_energy(u_value,base-step*left_direction- &
             step*right_direction)) / (4.0_DP*step*step)
  END FUNCTION dudarev_mixed_difference

  PURE SUBROUTINE independent_adjoint(matrix_in,matrix_out)
    COMPLEX(DP), INTENT(IN) :: matrix_in(ldim,ldim,4)
    COMPLEX(DP), INTENT(OUT) :: matrix_out(ldim,ldim,4)
    INTEGER :: io1, io2, iblock
    DO iblock = 1, 4
       DO io2 = 1, ldim
          DO io1 = 1, ldim
             matrix_out(io1,io2,iblock) = CONJG(matrix_in(io2,io1, &
                  hub_spin_transpose(iblock)))
          END DO
       END DO
    END DO
  END SUBROUTINE independent_adjoint

  PURE COMPLEX(DP) FUNCTION joint_transpose_pair( &
       left,right,conjugate_left) RESULT(value)
    COMPLEX(DP), INTENT(IN) :: left(ldim,ldim,4), right(ldim,ldim,4)
    LOGICAL, INTENT(IN) :: conjugate_left
    INTEGER :: io1, io2, iblock
    value = (0.0_DP,0.0_DP)
    DO iblock = 1, 4
       DO io2 = 1, ldim
          DO io1 = 1, ldim
             IF (conjugate_left) THEN
                value = value + CONJG(left(io1,io2,iblock)) * &
                     right(io2,io1,hub_spin_transpose(iblock))
             ELSE
                value = value + left(io1,io2,iblock) * &
                     right(io2,io1,hub_spin_transpose(iblock))
             ENDIF
          END DO
       END DO
    END DO
  END FUNCTION joint_transpose_pair

  PURE SUBROUTINE independent_full_adjoint(matrix_in,matrix_out)
    COMPLEX(DP), INTENT(IN) :: matrix_in(hdim,hdim)
    COMPLEX(DP), INTENT(OUT) :: matrix_out(hdim,hdim)
    matrix_out = TRANSPOSE(CONJG(matrix_in))
  END SUBROUTINE independent_full_adjoint

  PURE SUBROUTINE theta_matrix(matrix_in,jmatrix,matrix_out)
    COMPLEX(DP), INTENT(IN) :: matrix_in(hdim,hdim), jmatrix(hdim,hdim)
    COMPLEX(DP), INTENT(OUT) :: matrix_out(hdim,hdim)
    matrix_out = MATMUL(jmatrix,MATMUL(CONJG(matrix_in), &
         TRANSPOSE(CONJG(jmatrix))))
  END SUBROUTINE theta_matrix

  PURE SUBROUTINE theta_vector(vector_in,jmatrix,vector_out)
    COMPLEX(DP), INTENT(IN) :: vector_in(hdim), jmatrix(hdim,hdim)
    COMPLEX(DP), INTENT(OUT) :: vector_out(hdim)
    vector_out = MATMUL(jmatrix,CONJG(vector_in))
  END SUBROUTINE theta_vector

  SUBROUTINE check_separated(actual,reference,label,nfail)
    COMPLEX(DP), INTENT(IN) :: actual, reference
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nfail
    IF (ABS(actual-reference) <= 1.0E-8_DP) THEN
       WRITE(*,'(A,1X,ES12.4)') 'FAIL '//TRIM(label),ABS(actual-reference)
       nfail = nfail + 1
    ENDIF
  END SUBROUTINE check_separated

  SUBROUTINE check_matrix(actual,expected,label,nfail)
    COMPLEX(DP), INTENT(IN) :: actual(hdim,hdim), expected(hdim,hdim)
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nfail
    IF (MAXVAL(ABS(actual-expected)) > tol) THEN
       WRITE(*,'(A,1X,ES12.4)') 'FAIL '//TRIM(label), &
            MAXVAL(ABS(actual-expected))
       nfail = nfail + 1
    ENDIF
  END SUBROUTINE check_matrix

  SUBROUTINE check_vector(actual,expected,label,nfail)
    COMPLEX(DP), INTENT(IN) :: actual(hdim), expected(hdim)
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nfail
    IF (MAXVAL(ABS(actual-expected)) > tol) THEN
       WRITE(*,'(A,1X,ES12.4)') 'FAIL '//TRIM(label), &
            MAXVAL(ABS(actual-expected))
       nfail = nfail + 1
    ENDIF
  END SUBROUTINE check_vector

  SUBROUTINE check_vector_separated(actual,reference,label,nfail)
    COMPLEX(DP), INTENT(IN) :: actual(hdim), reference(hdim)
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nfail
    IF (MAXVAL(ABS(actual-reference)) <= 1.0E-8_DP) THEN
       WRITE(*,'(A,1X,ES12.4)') 'FAIL '//TRIM(label), &
            MAXVAL(ABS(actual-reference))
       nfail = nfail + 1
    ENDIF
  END SUBROUTINE check_vector_separated

  SUBROUTINE check_nonzero(actual,label,nfail)
    COMPLEX(DP), INTENT(IN) :: actual(hdim)
    CHARACTER(LEN=*), INTENT(IN) :: label
    INTEGER, INTENT(INOUT) :: nfail
    IF (MAXVAL(ABS(actual)) <= 1.0E-8_DP) THEN
       WRITE(*,'(A,1X,ES12.4)') 'FAIL '//TRIM(label),MAXVAL(ABS(actual))
       nfail = nfail + 1
    ENDIF
  END SUBROUTINE check_nonzero

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
