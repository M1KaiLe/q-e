!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dnsq_bare 
  !-----------------------------------------------------------------------
  !! DFPT+U : This routine calculates dnsbare, i.e. the bare variation 
  !! of the occupation matrix ns under the perturbation q, due to a 
  !! displacement of the atom L in the direction icart. 
  !
  !! $$ \text{dnsbare}(m1,m2,\text{ispin},I,\text{icart},L) = 
  !!    \sum_{k,n} [ \langle\psi(n,k,\text{ispin})| S_{k}\phi_(k,I,m1)\rangle \cdot 
  !!                 \langle\Delta^{L \text{icart}}_{-q}(S_{k+q}\phi(k+q,I,m2))
  !!                  |\psi(n,k,\text{ispin})\rangle  
  !!               + \langle \psi(n,k,\text{ispin})| S_{k}\phi(k,I,m2)> \cdot 
  !!                 \langle\Delta^{L \text{icart}}_{-q}(S_{k+q}\phi(k+q,I,m1))
  !!                  |\psi(n,k,\text{ispin})\rangle ] $$
  !
  !! Written  by A. Floris
  !! Modified by I. Timrov (01.10.2018)
  ! 
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : nwordwfcU
  USE units_lr,      ONLY : iuwfc, lrwfc
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE klist,         ONLY : xk, ngk, igk_k
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, offsetU, is_hubbard, nwfcU
  USE ldaU_ph,       ONLY : wfcatomk, wfcatomkpq, dwfcatomkpq, &
                            sdwfcatomk, sdwfcatomkpq, dvkb, vkbkpq, dvkbkpq, &
                            dnsbare, dnsbare_all_modes, proj1, proj2, read_dns_bare
  USE ldaU_lr,       ONLY : swfcatomk, swfcatomkpq
  USE wvfct,         ONLY : npwx, wg, nbnd 
  USE uspp,          ONLY : vkb, nkb, ofsbeta
  USE qpoint,        ONLY : nksq, ikks, ikqs
  USE control_lr,    ONLY : lgamma
  USE units_lr,      ONLY : iuatwfc, iuatswfc
  USE uspp_param,    ONLY : nh, nhm 
  USE lsda_mod,      ONLY : lsda, nspin, current_spin, isk
  USE noncollin_module, ONLY : noncolin, domag
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE mp_pools,      ONLY : inter_pool_comm
  USE mp_bands,      ONLY : intra_bgrp_comm 
  USE mp,            ONLY : mp_sum, mp_bcast
  USE mp_world,      ONLY : world_comm
  USE io_files,      ONLY : seqopn  
  USE buffers,       ONLY : get_buffer
  USE control_flags, ONLY : iverbosity
  USE wavefunctions, ONLY : evc  
  USE mp_images,     ONLY : intra_image_comm
  USE uspp_init,        ONLY : init_us_2
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  INTEGER :: i, j, k, icart, ios, is, na, nt, l, ih,        &
             ibeta, n, ihubst, ihubst1, ihubst2, nah, m, m1, m2,  &
             ibnd, ldim, npw, npwq, ik, ikk, ikq, iundnsbare
  COMPLEX(DP), ALLOCATABLE :: dqsphi(:,:), dmqsphi(:,:), dwfcatom_(:) 
  COMPLEX(DP), EXTERNAL :: ZDOTC
  LOGICAL :: exst
  ! 
  CALL start_clock( 'dnsq_bare' )
  !
  IF (noncolin) THEN
     CALL dnsq_bare_nc()
     CALL stop_clock( 'dnsq_bare' )
     RETURN
  ENDIF
  ! 
  ios = 0
  !
  ALLOCATE (dqsphi(npwx,nwfcU)) 
  ALLOCATE (dmqsphi(npwx,nwfcU))
  ALLOCATE (dwfcatom_(npwx))
  ALLOCATE (proj1(nbnd,nwfcU))
  ALLOCATE (proj2(nbnd,nwfcU))
  ! 
  proj1   = (0.d0, 0.d0) 
  proj2   = (0.d0, 0.d0) 
  dnsbare = (0.d0, 0.d0) 
  !
  ! Note: sdwfcatomkpq is not needed in this routine when we call delta_sphi,
  ! so we set it to zero and do not calculate.
  !
  sdwfcatomk   = (0.d0, 0.d0)
  sdwfcatomkpq = (0.d0, 0.d0)
  !
  ldim = 2 * Hubbard_lmax + 1
  !
  exst = .FALSE.
  !
  ! The unit number
  iundnsbare = 37
  !
  IF (ionode) CALL seqopn (iundnsbare, 'dnsbare', 'formatted', exst)
  !
  ! Read (unsymmetrized) dnsbare from file (if it was already computed)
  !
  IF (read_dns_bare) THEN 
     !   
     WRITE(stdout,*) 'Does the file dnsbare exist in tmp_dir and is it read correctly?'
     ! 
     IF (ionode) THEN 
        IF (exst) THEN
           READ(iundnsbare,*,iostat=ios) dnsbare
           IF (ios==0) THEN 
              WRITE(stdout,*) '...yes. the dnsbare matrix was read correctly from file'
           ELSE
              WRITE(stdout,*) '...no. the file exists but it seems to be corrupted'
           ENDIF
        ELSE
           WRITE(stdout,*) '...no. the file does not exist. '
        ENDIF
     ENDIF
     !
     CALL mp_bcast(exst, ionode_id, world_comm)           
     CALL mp_bcast(ios,  ionode_id, world_comm)
     IF (exst .AND. ios==0) CALL mp_bcast(dnsbare, ionode_id, world_comm) 
     !
     ! IT: Is it needed to broadcast for intra_image_comm?
     CALL mp_bcast(ios,  ionode_id, intra_image_comm)
     IF (exst .AND. ios==0) CALL mp_bcast(dnsbare, ionode_id, intra_image_comm)
     !
  ENDIF
  !
  ! Compute dnsbare (if it was not already done)
  !
  IF ((.NOT.exst) .OR. (ios/=0) .OR. (.NOT.read_dns_bare)) THEN
     !
     WRITE(stdout,'(/5x,"Calculating the dnsbare matrix...")')
     dnsbare = (0.d0, 0.d0)
     ! 
     DO ik = 1, nksq
        !
        ikk = ikks(ik)
        ikq = ikqs(ik)
        npw = ngk(ikk)
        npwq= ngk(ikq)
        !
        IF (lsda) current_spin = isk(ikk)
        !
        ! Read unperturbed KS wavefuctions psi(k)
        !
        IF (nksq.GT.1) CALL get_buffer (evc, lrwfc, iuwfc, ikk)
        !
        ! Compute beta functions at k (vkb) and at k+q (vkbkpq)
        ! 
        CALL init_us_2 (npw, igk_k(1,ikk), xk(:,ikk), vkb)
        IF (.NOT.lgamma) CALL init_us_2 (npwq, igk_k(1,ikq), xk(:,ikq), vkbkpq)
        ! 
        ! Compute the derivatives of beta functions d^{icart}beta at k and k+q
        ! for all bands and for the 3 cartesian directions
        !
        DO icart = 1, 3
           DO na = 1, nat 
              nt = ityp(na)
              DO ih = 1, nh(nt)
                 !
                 ibeta = ofsbeta(na) + ih
                 !
                 CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                            vkb(:,ibeta), dvkb(:,ibeta,icart))
                 !
                 IF (.NOT.lgamma) &
                 CALL dwfc (npwq, igk_k(1,ikq), ikq, icart, & 
                            vkbkpq(:,ibeta), dvkbkpq(:,ibeta,icart))
                 !
              ENDDO
           ENDDO
        ENDDO
        ! 
        ! Read the atomic orbitals \phi at k and k+q from file (unit iuatwfc)
        ! 
        CALL get_buffer (wfcatomk, nwordwfcU, iuatwfc, ikk)
        IF (.NOT.lgamma) CALL get_buffer (wfcatomkpq, nwordwfcU, iuatwfc, ikq)
        ! 
        ! Read S*\phi at k from file (unit iuatswfc)
        !
        CALL get_buffer (swfcatomk, nwordwfcU, iuatswfc, ikk)      
        !
        DO na = 1, nat
           ! 
           DO icart = 1, 3 
              ! 
              dqsphi  = (0.d0, 0.d0)  
              dmqsphi = (0.d0, 0.d0)
              !
              DO nah = 1, nat
                 !
                 nt = ityp(nah)
                 !
                 IF (is_hubbard(nt)) THEN
                    !    
                    DO m = 1, 2*Hubbard_l(nt)+1
                       !
                       ihubst = offsetU(nah) + m ! I m index
                       ! 
                       IF (nah==na) THEN 
                          !  
                          ! Calculate |d_icart\phi_(k,I,m)) >
                          !
                          CALL dwfc (npw, igk_k(1,ikk), ikk, icart, &
                                     wfcatomk(:,ihubst), dwfcatom_) 
                          !
                          ! Apply the S matrix: | S_{k} d_^(I,icart)\phi_(k,I,m) > 
                          !
                          CALL swfc (npw, 1, vkb, dwfcatom_, sdwfcatomk(:,ihubst))
                          !
                       ENDIF
                       !
                       ! Calculate dmqsphi = |\Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) >
                       !
                       CALL delta_sphi (ikk, ikq, na, icart, nah, ihubst, wfcatomk, wfcatomkpq, &
                                        sdwfcatomk, sdwfcatomkpq, vkb, vkbkpq, dvkb(:,:,icart), &
                                        dvkbkpq(:,:,icart), dqsphi, dmqsphi, 0)  
                       !
                       ! Calculate: 
                       ! proj1 (ibnd,ihubst) = < S_{k}\phi_(k,I,m) | psi(ibnd,k) >
                       ! proj2 (ibnd,ihubst) = < \Delta_{-q}(S_{k+q} \phi_(k+q,I,m)) | psi(ibnd,k)>
                       !
                       DO ibnd = 1, nbnd
                          proj1 (ibnd,ihubst) = ZDOTC (npw, swfcatomk(1,ihubst), 1, evc(1,ibnd), 1)
                          proj2 (ibnd,ihubst) = ZDOTC (npw, dmqsphi(1,ihubst ), 1, evc(1,ibnd), 1)
                       ENDDO
                       !  
                    ENDDO 
                    !
                 ENDIF
                 !
              ENDDO
              !
              CALL mp_sum (proj1, intra_bgrp_comm)
              CALL mp_sum (proj2, intra_bgrp_comm)
              !
              ! Finally calculate dnsbare.
              ! dnsbare is symmetric: dnsbare(m1,m2) = dnsbare(m2,m1)
              !
              DO nah = 1, nat
                 ! 
                 nt = ityp(nah)
                 !
                 IF (is_hubbard(nt)) THEN
                    !  
                    DO m1 = 1, 2*Hubbard_l(nt)+1
                       !
                       ihubst1 = offsetU(nah) + m1
                       !
                       DO m2 = m1, 2*Hubbard_l(nt)+1
                          !
                          ihubst2 = offsetU(nah) + m2
                          !
                          DO ibnd = 1, nbnd
                             dnsbare (m1,m2,current_spin,nah,icart,na) = &
                             dnsbare (m1,m2,current_spin,nah,icart,na) + wg(ibnd,ikk) * &
                                  ( CONJG(proj1(ibnd,ihubst1)) * proj2(ibnd,ihubst2)  + &
                                    CONJG(proj1(ibnd,ihubst2)) * proj2(ibnd,ihubst1) )
                          ENDDO
                          ! 
                       ENDDO
                       !
                    ENDDO
                    !
                 ENDIF
                 !
              ENDDO
              !
           ENDDO ! icart
           !
        ENDDO ! na
        !
     ENDDO ! ik
     !
     CALL mp_sum (dnsbare, inter_pool_comm)
     !
     ! Filling the m1 m2 lower triangular part 
     !
     DO m1 = 2, ldim
        DO m2 = 1, m1-1
           dnsbare(m1,m2,:,:,:,:) = dnsbare(m2,m1,:,:,:,:)
        ENDDO
     ENDDO
     !
     ! In nspin.eq.1 k-point weight is normalized to 2 el/band 
     ! in the whole BZ but we are interested in dns of one spin component
     !
     IF (nspin.EQ.1) dnsbare = 0.5d0 * dnsbare
     ! 
     ! Write (unsymmetrized) dnsbare to file
     !  
     IF (ionode) WRITE(iundnsbare,*) dnsbare
     !
  ENDIF
  !
  ! Close the unit iundnsbare
  !
  IF (ionode) CLOSE (unit=iundnsbare,status='keep')   
  !
  ! Symmetrize dnsbare, and also save the result
  ! in the pattern basis dnsbare_all_modes
  !
  CALL sym_dns_wrapper (ldim, dnsbare, dnsbare_all_modes) 
  !
  ! Write symmetryzed dnsbare in cartesian coordinates 
  ! to the standard output
  ! 
  IF (iverbosity==1) THEN 
     WRITE(stdout,*) 'DNS_BARE SYMMETRIZED IN CARTESIAN COORDINATES'
     DO na = 1, nat        
        DO icart = 1, 3
           WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') 'displaced atom L =', na, 'ipol=', icart
           DO nah = 1, nat  
              nt = ityp(nah)
              IF (is_hubbard(nt)) THEN
                 DO is = 1, nspin
                    WRITE(stdout,'(a,1x,i2,2x,a,1x,i2)') ' Hubbard atom', nah, 'spin', is
                    DO m1 = 1, ldim
                       WRITE(stdout,'(10(f15.10,1x))') dnsbare (m1,:,is,nah,icart,na)
                    ENDDO
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
     ENDDO
     WRITE(stdout,*)
  ENDIF
  ! 
  DEALLOCATE (dqsphi)
  DEALLOCATE (dmqsphi)
  DEALLOCATE (dwfcatom_)
  DEALLOCATE (proj1 )
  DEALLOCATE (proj2 )
  ! 
  CALL stop_clock ('dnsq_bare')
  ! 
  RETURN
  !
END SUBROUTINE dnsq_bare
!----------------------------------------------------------------------------------------

!----------------------------------------------------------------------------
SUBROUTINE dnsq_bare_nc()
  !----------------------------------------------------------------------------
  !! Bare projector contribution to the noncollinear Hubbard occupation
  !! response. This path is restricted to norm-conserving pseudopotentials.
  !
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : nwordwfcU, seqopn
  USE units_lr,      ONLY : iuwfc, lrwfc, iuatwfc, iuatswfc
  USE ions_base,     ONLY : nat, ityp
  USE klist,         ONLY : ngk, igk_k
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, offsetU, is_hubbard, nwfcU
  USE ldaU_ph,       ONLY : wfcatomk, dnsbare, dnsbare_all_modes, read_dns_bare
  USE ldaU_lr,       ONLY : swfcatomk
  USE wvfct,         ONLY : npwx, wg, nbnd
  USE qpoint,        ONLY : nksq, ikks
  USE qpoint_aux,    ONLY : ikmks
  USE noncollin_module, ONLY : npol, domag
  USE lsda_mod,      ONLY : nspin
  USE wavefunctions, ONLY : evc
  USE buffers,       ONLY : get_buffer
  USE mp,            ONLY : mp_sum, mp_bcast
  USE mp_pools,      ONLY : inter_pool_comm
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp_world,      ONLY : world_comm
  USE mp_images,     ONLY : intra_image_comm
  USE io_global,     ONLY : ionode, ionode_id, stdout
  USE control_ph,    ONLY : current_iq
  USE hubbard_nc_response, ONLY : hub_spin_index, hubbard_nc_format, &
                                  hubbard_branch_indices_nc, &
                                  hubbard_kramers_indices_nc
  !
  IMPLICIT NONE
  INTEGER :: isolv, nsolv, ik, ikk, ikmk, npw, na, nah, nt, icart, ibnd
  INTEGER :: is1, is2, js, is, m, m1, m2, ldim, ldim_nt
  INTEGER :: ihubst, ihubst1, ihubst2, iunit, ipattern, ios
  INTEGER :: mbra, mket, sbra, sket, bsign, ksign
  INTEGER :: mbra_k, mket_k, sbra_k, sket_k, ihubst1_k, ihubst2_k
  LOGICAL :: exst, valid_restart
  CHARACTER(LEN=80) :: header
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  COMPLEX(DP), ALLOCATABLE :: dphi(:,:), dtmp(:)
  COMPLEX(DP), ALLOCATABLE :: proj(:,:), dproj(:,:)
  COMPLEX(DP), ALLOCATABLE :: dns_branch(:,:,:,:,:,:)
  !
  IF (npol /= 2 .OR. nspin /= 4) CALL errore('dnsq_bare_nc', &
       'inconsistent noncollinear spin dimensions', 1)
  ldim = 2 * Hubbard_lmax + 1
  ALLOCATE(dphi(npwx*npol,nwfcU), dtmp(npwx))
  ALLOCATE(proj(nbnd,nwfcU), dproj(nbnd,nwfcU))
  ALLOCATE(dns_branch(ldim,ldim,4,nat,3,nat))
  dnsbare = (0.0_DP, 0.0_DP)
  nsolv = MERGE(2, 1, domag)
  !
  iunit = 37
  exst = .FALSE.
  valid_restart = .FALSE.
  IF (ionode) CALL seqopn(iunit, 'dnsbare_nc', 'formatted', exst)
  IF (read_dns_bare .AND. ionode .AND. exst) THEN
     READ(iunit,'(A)',IOSTAT=ios) header
     IF (ios == 0 .AND. TRIM(header) == hubbard_nc_format) THEN
        READ(iunit,*,IOSTAT=ios) dnsbare
        valid_restart = (ios == 0)
     ENDIF
  ENDIF
  CALL mp_bcast(valid_restart, ionode_id, world_comm)
  IF (valid_restart) THEN
     CALL mp_bcast(dnsbare, ionode_id, world_comm)
     CALL mp_bcast(dnsbare, ionode_id, intra_image_comm)
  ELSE
     WRITE(stdout,'(/5x,a)') 'Calculating noncollinear dnsbare matrix...'
      DO isolv = 1, nsolv
         dns_branch = (0.0_DP, 0.0_DP)
         DO ik = 1, nksq
            ikk = ikks(ik)
            IF (isolv == 1) THEN
               ikmk = ikk
            ELSE
               ikmk = ikmks(ik)
            ENDIF
           ! apply_trev stores the second branch in the direct-k G ordering.
           npw = ngk(ikk)
           CALL get_buffer(evc, lrwfc, iuwfc, ikmk)
           CALL get_buffer(wfcatomk, nwordwfcU, iuatwfc, ikk)
           CALL get_buffer(swfcatomk, nwordwfcU, iuatswfc, ikk)
           !
           proj = (0.0_DP, 0.0_DP)
           DO nah = 1, nat
              nt = ityp(nah)
              IF (.NOT. is_hubbard(nt)) CYCLE
              ldim_nt = 2 * Hubbard_l(nt) + 1
              DO is1 = 1, 2
                 DO m = 1, ldim_nt
                    ihubst = offsetU(nah) + m + ldim_nt*(is1-1)
                    DO ibnd = 1, nbnd
                       proj(ibnd,ihubst) = DOT_PRODUCT( &
                            swfcatomk(1:npw,ihubst),evc(1:npw,ibnd)) + &
                            DOT_PRODUCT(swfcatomk(npwx+1:npwx+npw,ihubst), &
                                        evc(npwx+1:npwx+npw,ibnd))
                    END DO
                 END DO
              END DO
           END DO
           CALL mp_sum(proj, intra_bgrp_comm)
           !
           DO na = 1, nat
              DO icart = 1, 3
                 dphi = (0.0_DP, 0.0_DP)
                 nt = ityp(na)
                 IF (is_hubbard(nt)) THEN
                    ldim_nt = 2 * Hubbard_l(nt) + 1
                    DO is1 = 1, 2
                       DO m = 1, ldim_nt
                          ihubst = offsetU(na) + m + ldim_nt*(is1-1)
                          DO js = 1, 2
                              CALL dwfc(npw,igk_k(1,ikk),ikk,icart, &
                                  wfcatomk(1+(js-1)*npwx,ihubst),dtmp)
                             dphi(1+(js-1)*npwx:npw+(js-1)*npwx,ihubst) = &
                                  dtmp(1:npw)
                          END DO
                       END DO
                    END DO
                 ENDIF
                 !
                 dproj = (0.0_DP, 0.0_DP)
                 DO nah = 1, nat
                    nt = ityp(nah)
                    IF (.NOT. is_hubbard(nt)) CYCLE
                    ldim_nt = 2 * Hubbard_l(nt) + 1
                    DO is1 = 1, 2
                       DO m = 1, ldim_nt
                          ihubst = offsetU(nah) + m + ldim_nt*(is1-1)
                          DO ibnd = 1, nbnd
                             dproj(ibnd,ihubst) = DOT_PRODUCT( &
                                  dphi(1:npw,ihubst),evc(1:npw,ibnd)) + &
                                  DOT_PRODUCT( &
                                  dphi(npwx+1:npwx+npw,ihubst), &
                                  evc(npwx+1:npwx+npw,ibnd))
                          END DO
                       END DO
                    END DO
                 END DO
                 CALL mp_sum(dproj, intra_bgrp_comm)
                 !
                 DO nah = 1, nat
                    nt = ityp(nah)
                    IF (.NOT. is_hubbard(nt)) CYCLE
                    ldim_nt = 2 * Hubbard_l(nt) + 1
                    DO is1 = 1, 2
                       DO is2 = 1, 2
                          is = hub_spin_index(is1,is2)
                          DO m1 = 1, ldim_nt
                              DO m2 = 1, ldim_nt
                                 CALL hubbard_branch_indices_nc(isolv,m1,m2,is1,is2, &
                                      mbra,mket,sbra,sket,bsign)
                                 ihubst1 = offsetU(nah) + mbra + ldim_nt*(sbra-1)
                                 ihubst2 = offsetU(nah) + mket + ldim_nt*(sket-1)
                                 DO ibnd = 1, nbnd
                                     dns_branch(m1,m2,is,nah,icart,na) = &
                                          dns_branch(m1,m2,is,nah,icart,na) + &
                                          REAL(bsign,DP) * wg(ibnd,ikk) * &
                                          CONJG(proj(ibnd,ihubst1)) * &
                                          dproj(ibnd,ihubst2)
                                      ! Complete the nonmagnetic spinor response
                                      ! with its fixed-q Kramers bra partner.
                                      IF (.NOT. domag) THEN
                                         CALL hubbard_kramers_indices_nc(m1,m2,is1,is2, &
                                              mbra_k,mket_k,sbra_k,sket_k,ksign)
                                         ihubst1_k = offsetU(nah) + mbra_k + &
                                              ldim_nt*(sbra_k-1)
                                         ihubst2_k = offsetU(nah) + mket_k + &
                                              ldim_nt*(sket_k-1)
                                         dns_branch(m1,m2,is,nah,icart,na) = &
                                              dns_branch(m1,m2,is,nah,icart,na) + &
                                              REAL(ksign,DP) * wg(ibnd,ikk) * &
                                              CONJG(proj(ibnd,ihubst1_k)) * &
                                              dproj(ibnd,ihubst2_k)
                                      ENDIF
                                END DO
                             END DO
                          END DO
                       END DO
                    END DO
                 END DO
              END DO
           END DO
         END DO
         CALL mp_sum(dns_branch, inter_pool_comm)
         dnsbare = dnsbare + dns_branch
      END DO
      IF (ionode) THEN
        REWIND(iunit)
        WRITE(iunit,'(A)') hubbard_nc_format
        WRITE(iunit,*) dnsbare
     ENDIF
  ENDIF
  IF (ionode) CLOSE(iunit,STATUS='keep')
  !
  CALL sym_dns_wrapper(ldim, dnsbare, dnsbare_all_modes)
  IF (ionode) THEN
     ipattern = 38
     CALL seqopn(ipattern,'dnsbare_nc_pattern_q' // &
          TRIM(int_to_char(current_iq)),'formatted',exst)
     WRITE(ipattern,'(A)') hubbard_nc_format
     WRITE(ipattern,'("ldim=",I0,1X,"nspin=4 nat=",I0,1X,"nmodes=",I0,1X,&
          &"spin_order=uu,ud,du,dd")') ldim, nat, 3*nat
     WRITE(ipattern,*) dnsbare_all_modes
     CLOSE(ipattern,STATUS='keep')
  ENDIF
  DEALLOCATE(dphi, dtmp, proj, dproj, dns_branch)
END SUBROUTINE dnsq_bare_nc
!----------------------------------------------------------------------------
