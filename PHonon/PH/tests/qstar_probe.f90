! Isolated audit: real QE numerical routines, captured output instead of disk I/O.
module io_global
 integer, parameter :: stdout=6
end module
module symm_base
 integer :: t_rev(48)=0
end module
module control_ph
 logical :: xmldyn=.true.
end module
module io_dyn_mat
 use kinds, only: dp
 complex(dp), allocatable :: saved(:,:,:,:,:)
 real(dp) :: saved_q(3,96)
 integer :: saved_n=0
contains
 subroutine write_dyn_mat(nat, counter, q, phi)
 integer :: nat,counter
 real(dp) :: q(3)
 complex(dp) :: phi(3,3,nat,nat)
 saved(:,:,:,:,counter)=phi
 saved_q(:,counter)=q
 saved_n=counter
 end subroutine
end module
subroutine errore(a,b,n)
 character(*) :: a,b
 integer :: n
 print *,a,b,n
 stop 1
end subroutine
subroutine zcopy(n,x,incx,y,incy)
 use kinds, only: dp
 integer :: n,incx,incy,i
 complex(dp) :: x(*),y(*)
 do i=1,n
 y(1+(i-1)*incy)=x(1+(i-1)*incx)
 enddo
end subroutine
subroutine write_dyn_on_file(q,phi,nat,unit)
 use kinds, only: dp
 integer :: nat,unit
 real(dp) :: q(3)
 complex(dp) :: phi(3,3,nat,nat)
 stop 'Unexpected non-XML output'
end subroutine
subroutine audit_run(nat,nsym,nsymq,tr,s,invs,irt,rtau,at,bg,q,dyn,project, &
                     nwritten,qq,output,before,after)
 use kinds, only: dp
 use symm_base, only: t_rev
 use io_dyn_mat
 integer :: nat,nsym,nsymq,tr(48),s(3,3,48),invs(48),irt(48,nat)
 integer :: project,nwritten,nq,isq(48),imq,na,nb
 real(dp) :: rtau(3,48,nat),at(3,3),bg(3,3),q(3),qq(3,96),sxq(3,48)
 complex(dp) :: dyn(3*nat,3*nat),output(3,3,nat,nat,96)
 complex(dp) :: before(3*nat,3*nat),after(3*nat,3*nat),phi(3,3,nat,nat)
 t_rev=tr
 allocate(saved(3,3,nat,nat,96))
 saved=0; saved_q=0; saved_n=0
 before=dyn
 call star_q1(q,at,bg,nsym,s,invs,nq,sxq,isq,imq,.false.,tr)
 call q2qstar_ph(dyn,at,bg,nat,nsym,s,invs,irt,rtau,nq,sxq,isq,imq,0,tr)
 after=dyn
 nwritten=saved_n; qq=saved_q; output=saved
 deallocate(saved)
end subroutine
subroutine get_constants(amu,cm)
 use constants, only: amu_ry,ry_to_cmm1
 use kinds, only: dp
 real(dp) :: amu,cm
 amu=amu_ry; cm=ry_to_cmm1
end subroutine
