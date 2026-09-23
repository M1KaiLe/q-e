module kinds
 integer, parameter :: dp=kind(1.d0)
end module
module noncollin_module
 logical :: noncolin,domag
 integer :: npol,nspin_mag
end module
module mp_bands
 integer :: me_bgrp=0
end module
module fft_base
 type fft_t
  logical :: has_task_groups=.false.
  integer :: nr1x=3,nr2x=1,my_nr3p=1
 end type
 type(fft_t) :: dffts
end module
