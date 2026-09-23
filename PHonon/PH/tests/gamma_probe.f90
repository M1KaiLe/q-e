program probe
 use kinds
 use noncollin_module
 use fft_base
 implicit none
 integer,parameter :: nr=3
 complex(dp),allocatable,target :: dvscfin(:,:,:)
 complex(dp),pointer :: dvscfins(:,:,:)
 complex(dp),allocatable :: base(:,:),good(:,:),psi(:,:),out(:,:),ref(:,:),orig(:,:),diff(:,:)
 complex(dp) :: def(1), c
 real(dp) :: olderr,newerr,formulaerr,griderr
 logical :: doublegrid,lmetq0
 integer :: itype,itg,idg,ish,ir,is,ipert,spin,ncase,imet
 ncase=0
 lmetq0=.true.
 ipert=1
 do imet=0,1
 lmetq0=imet==1
 do itype=1,4
  noncolin=itype>=3
  domag=itype==4
  npol=1
  if(noncolin) npol=2
  nspin_mag=1
  if(itype==2) nspin_mag=2
  if(domag) nspin_mag=4
  do itg=0,1
   ! Collinear task-group code expects a preselected spin channel;
   ! this probe uses the full channel array, so test LSDA without task groups.
   if(itype==2.and.itg==1) cycle
   dffts%has_task_groups=itg==1
   do idg=0,1
    doublegrid=idg==1
    do ish=1,3
     def(1)=cmplx(0.125_dp,0._dp,dp)
     if(ish==2) def(1)=cmplx(0.125_dp,0.0625_dp,dp)
     if(ish==3) def(1)=cmplx(0._dp,0._dp,dp)
     c=def(1)
     if(.not.lmetq0) c=cmplx(0._dp,0._dp,dp)
     allocate(dvscfin(nr,nspin_mag,1),base(nr,nspin_mag),good(nr,nspin_mag))
     allocate(psi(nr,npol),out(nr,npol),ref(nr,npol),orig(nr,npol),diff(nr,npol))
     do is=1,nspin_mag
      do ir=1,nr
       base(ir,is)=cmplx(0.23_dp*ir+0.17_dp*is,-0.11_dp*ir+0.08_dp*is,dp)
      enddo
     enddo
     do is=1,npol
      do ir=1,nr
       psi(ir,is)=cmplx(0.19_dp*ir-0.13_dp*is,0.07_dp*ir+0.29_dp*is,dp)
      enddo
     enddo
     dvscfin(:,:,1)=base
     if(doublegrid) then
      allocate(dvscfins(nr,nspin_mag,1))
      dvscfins(:,:,1)=base
     else
      dvscfins=>dvscfin
     endif
     ! Verbatim block extracted from the current production solve_linter.f90.
     include 'export_shift.inc'
     griderr=maxval(abs(dvscfins-dvscfin))
     good=base
     if(noncolin) then
      good(:,1)=good(:,1)-c
     else
      good=good-c
     endif
     olderr=0._dp
     newerr=0._dp
     formulaerr=0._dp
     do spin=1,merge(2,1,itype==2)
      orig=psi
      call apply_dpot(nr,orig,base,spin)
      ref=orig-c*psi
      out=psi
      call apply_dpot(nr,out,dvscfins(:,:,1),spin)
      olderr=max(olderr,maxval(abs(out-ref)))
      diff=cmplx(0._dp,0._dp,dp)
      formulaerr=max(formulaerr,maxval(abs(out-ref-diff)))
      out=psi
      call apply_dpot(nr,out,good,spin)
      newerr=max(newerr,maxval(abs(out-ref)))
     enddo
     if(max(newerr,formulaerr,griderr)>1.e-12_dp) stop 1
     if(olderr>1.e-12_dp) stop 3
     ncase=ncase+1
     write(*,'(a,5i3,4es14.5)') 'met,type,tg,dg,shift; actual,reference,formula,grid errors: ', &
         imet,itype,itg,idg,ish,olderr,newerr,formulaerr,griderr
     if(doublegrid) deallocate(dvscfins)
     nullify(dvscfins)
     deallocate(dvscfin,base,good,psi,out,ref,orig,diff)
    enddo
   enddo
  enddo
 enddo
 enddo
 write(*,'(a,i0)') 'PASS: fixed-source configurations = ',ncase
end program
