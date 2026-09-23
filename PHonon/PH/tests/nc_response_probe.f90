program check_nc_response
  use kinds, only: dp
  use ph_nc_response_conv
  use, intrinsic :: ieee_arithmetic
  implicit none
  complex(dp) :: a(8,4,2), b(8,4,2), ar(8,4,2), br(8,4,2)
  real(dp) :: sums(3,2,2), left(3,2,2), right(3,2,2), changed(3,2,2)
  real(dp) :: err, referr, rms(2,2), reference(2,2), scaled(2,2), nan
  logical :: mask(8), half(8)
  integer :: i,c,j
  nan=ieee_value(0.0_dp,ieee_quiet_nan)
  do j=1,2
    do c=1,4
      do i=1,8
        a(i,c,j)=cmplx(0.03_dp*i*c*j,0.02_dp*(i+c),dp)
        b(i,c,j)=a(i,c,j)+cmplx(0.001_dp*c*j,-0.002_dp*i,dp)
      enddo
    enddo
  enddo
  mask=.true.; mask(5:8)=.false.
  call nc_response_sums(8,2,mask,a,b,sums)
  call nc_response_metric(4.0_dp,2,sums,referr,rms,reference,scaled)
  if (.not.(referr>0.0_dp)) stop 1
  ! Padding, including NaNs, must never enter the physical residual.
  a(5:8,:,:)=cmplx(nan,nan,dp); b(5:8,:,:)=cmplx(nan,nan,dp)
  call nc_response_sums(8,2,mask,a,b,changed)
  if(any(changed/=sums)) stop 2
  ! Splitting the physical grid between ranks must preserve the result.
  half=mask; half(3:4)=.false.
  call nc_response_sums(8,2,half,a,b,left)
  half=mask.and..not.half
  call nc_response_sums(8,2,half,a,b,right)
  if(maxval(abs(left+right-sums))>1e-14_dp) stop 3
  ! A rotation of the magnetic coordinate axes preserves the B-vector norm.
  ar=a; br=b
  ar(:,2,:)=a(:,3,:); ar(:,3,:)=-a(:,2,:)
  br(:,2,:)=b(:,3,:); br(:,3,:)=-b(:,2,:)
  call nc_response_sums(8,2,mask,ar,br,changed)
  if(maxval(abs(changed-sums))>1e-14_dp) stop 4
  ! Replicating samples or an irrep must not artificially improve convergence.
  call nc_response_metric(8.0_dp,2,2*sums,err,rms,reference,scaled)
  if(abs(err-referr)>1e-15_dp) stop 5
  changed=sums; changed(:,:,1)=sums(:,:,2)
  call nc_response_metric(4.0_dp,2,changed,err,rms,reference,scaled)
  if(abs(err-referr)>1e-15_dp) stop 6
  if(.not.nc_response_stop(referr,2*referr,.true.)) stop 7
  if(nc_response_stop(referr,2*referr,.false.)) stop 8
  if(nc_response_stop(referr,referr,.true.)) stop 9
  if(nc_response_stop(nan,1.0_dp,.true.)) stop 10
  if(nc_response_stop(-1.0_dp,1.0_dp,.true.)) stop 11
  changed=sums; changed(1,1,1)=nan
  call nc_response_metric(4.0_dp,2,changed,err,rms,reference,scaled)
  if(nc_response_stop(err,1.0_dp,.true.)) stop 12
  print *, 'PASS NC residual: padding, partition, magnetic rotation, grid/irrep scaling and stop gates'
end program
