module kinds
  implicit none
  integer, parameter :: dp=kind(1.d0)
end module
module noncollin_module
  logical :: noncolin=.true., domag=.true.
end module
module fft_types
  type fft_type_descriptor
    integer :: nnr, ngm
  end type
end module
module gc_lr
  use kinds
  real(dp), allocatable :: gmag(:,:,:), vsgga(:), segni(:)
end module
module probe_fft
  use kinds
contains
  subroutine derivative(n, input, q, output)
    integer, intent(in) :: n
    complex(dp), intent(in) :: input(n)
    real(dp), intent(in) :: q
    complex(dp), intent(out) :: output(n)
    complex(dp) :: coefficient, phase
    real(dp) :: angle, pi
    integer :: j,k,l,wave
    pi=acos(-1.d0)
    output=0
    do k=0,n-1
      wave=k
      if(k>n/2) wave=k-n
      coefficient=0
      do j=0,n-1
        angle=-2*pi*k*j/n
        coefficient=coefficient+input(j+1)*cmplx(cos(angle),sin(angle),dp)/n
      enddo
      do l=0,n-1
        angle=2*pi*k*l/n
        phase=cmplx(cos(angle),sin(angle),dp)
        output(l+1)=output(l+1)+cmplx(0.d0,wave+q,dp)*coefficient*phase
      enddo
    enddo
  end subroutine
end module
subroutine fft_qgradient(dfft, input, q, g, output)
  use kinds
  use fft_types
  use probe_fft
  implicit none
  type(fft_type_descriptor),intent(in) :: dfft
  complex(dp),intent(in) :: input(dfft%nnr)
  real(dp),intent(in) :: q(3),g(3,dfft%ngm)
  complex(dp),intent(out) :: output(3,dfft%nnr)
  call derivative(dfft%nnr,input,q(1),output(1,:))
  output(2,:)=cmplx(0.d0,q(2),dp)*input
  output(3,:)=cmplx(0.d0,q(3),dp)*input
end subroutine
subroutine fft_qgraddot(dfft, input, q, g, output)
  use kinds
  use fft_types
  use probe_fft
  implicit none
  type(fft_type_descriptor),intent(in) :: dfft
  complex(dp),intent(in) :: input(3,dfft%nnr)
  real(dp),intent(in) :: q(3),g(3,dfft%ngm)
  complex(dp),intent(out) :: output(dfft%nnr)
  call derivative(dfft%nnr,input(1,:),q(1),output)
  output=output+cmplx(0.d0,q(2),dp)*input(2,:)+cmplx(0.d0,q(3),dp)*input(3,:)
end subroutine
subroutine errore(where, message, code)
  character(*) :: where, message
  integer :: code
  write(*,*) where, message, code
  stop 1
end subroutine
program probe
  use kinds
  use fft_types
  use gc_lr
  implicit none
  type(fft_type_descriptor) :: dfft
  integer :: n,k,is
  real(dp) :: q(3),a,c
  real(dp),allocatable :: rho(:,:),grho(:,:,:),zero(:,:,:),ss(:,:,:),g(:,:)
  complex(dp),allocatable :: drho(:,:),dv(:,:)
  read(*,*) n,q,a,c
  dfft%nnr=n; dfft%ngm=n
  allocate(rho(n,4),grho(3,n,2),zero(n,2,2),ss(n,2,2),g(3,n))
  allocate(drho(n,4),dv(n,4),gmag(3,n,4),vsgga(n),segni(n))
  read(*,*) rho
  read(*,*) drho
  read(*,*) gmag
  read(*,*) vsgga
  read(*,*) segni
  grho=0; zero=0; g=0; dv=0; ss=c
  ss(:,1,1)=a;ss(:,2,2)=a
  call dgradcorr(dfft,rho,grho,zero,zero,zero,ss,q,drho,4,2,g,dv)
  do is=1,4
    do k=1,n
      write(*,'(2ES25.16)') real(dv(k,is)),aimag(dv(k,is))
    enddo
  enddo
end program
