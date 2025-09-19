PROGRAM compute_vsedi
!----------------------------------------------------------------
!This scheme computes optical depth assuming it as twice the 
!surface area
!----------------------------------------------------------------

USE netcdf
USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
IMPLICIT NONE
INTEGER, PARAMETER       ::                                         &
       wp    = SELECTED_REAL_KIND (12,200),                       &
       iintegers = KIND  (1)

REAL(KIND=wp), PARAMETER :: &   !collection of parameters of size distributions
    !cloud= cloud_nue1mue1:
    a_CD = .124, & !m/kg^beta      D=ax^b
    b_CD = 1./3., &
    nuCD = 1., &                    !nu of f(x) cloud
    muCD = 1., &                    !mu  "  " cloud
    avelCD = 3.75e+05, &
    bvelCD = 2./3., &


    !rain= rainSBB:
    a_RD = 0.124 , &
    b_RD = 1./3., &
    nuRD = 1.0, &
    muRD = 1./3., &
    avelRD = 114.0137, &
    bvelRD = 0.234370, &

    !ice= ice_cosmo5:
    a_ID = 0.835, &
    b_ID = 0.39, &
    nuID = 0.0, &
    muID = 1./3., &
    avelID = 2.77e+01, &
    bvelID = 0.215790, &

    !snow= snowSBB:
    a_SD = 5.13, &
    b_SD = 0.5, &
    nuSD = 0.0, &
    muSD = 0.5, &
    avelSD = 2.77e+01, &
    bvelSD = 0.215790, &


    !graupel= graupelhail_cosmo5:
    a_GD = 0.142, &
    b_GD = 0.314, &
    nuGD = 1.0, &
    muGD = 1./3., &
    avelGD = 86.89371, &
    bvelGD = 0.268325, &

    rho_vel=0.4, &   !taking from ICON routine
    rho_vel_c=0.2, & !taking from ICON routine
    rho0=1.225       !kg/m3


!number of bins for dr bins in size distribution
integer(KIND=iintegers)::     nbins=1000

real(kind=wp), parameter :: pi = acos(-1.0_wp)

REAL (KIND = wp):: xxCD(2), xxRD(2), xxID(2), xxSD(2), xxGD(2), &
                   ddCD(2), ddRD(2), ddID(2), ddSD(2), ddGD(2)


INTEGER(KIND=iintegers)         :: ze,ncells,n,k,ierr
REAL(KIND=wp) :: qc1, qnc1, apr, qr1, qnr1, qg1, qng1, qs1, qns1, &
                 qi1, qni1, A, l

CHARACTER(:), ALLOCATABLE :: tmp
CHARACTER(:), ALLOCATABLE :: fname, fname2, fname3, fname4, fname5, fname6
REAL(KIND=wp), ALLOCATABLE :: rhocorr_c(:,:), rhocorr(:,:), rho(:,:), qr(:,:), qc(:,:), &
                              qi(:,:), qs(:,:), qg(:,:), qnr(:,:), qnc(:,:), &
                              qni(:,:), qns(:,:), qng(:,:), &
                              vm_c(:,:), & !mass-weighted sedimentation speed cloud
                              vm_r(:,:), & !mass-weighted sedimentation speed rain
                              vm_i(:,:), & !ice crystals
                              vm_s(:,:), & !snow
                              vm_g(:,:)    !graupel

LOGICAL :: l_ice

INTEGER :: i
INTEGER :: argc

! Define size ranges [x_min, x_max]
xxCD = [4.2e-15, 2.6e-10] ! kg
xxRD = [2.6e-10, 6.5e-5] ! kg
xxID = [1.0e-12, 1.0e-5] ! kg
xxSD = [1.0e-10, 2.0e-5] ! kg
xxGD = [4.19e-09, 5.3e-4] ! kg
! Calculate dd values
ddCD(1) = a_CD * xxCD(1)**b_CD
ddCD(2) = a_CD * xxCD(2)**b_CD
ddRD(1) = a_RD * xxRD(1)**b_RD
ddRD(2) = a_RD * xxRD(2)**b_RD
ddID(1) = a_ID * xxID(1)**b_ID
ddID(2) = a_ID * xxID(2)**b_ID
ddSD(1) = a_SD * xxSD(1)**b_SD
ddSD(2) = a_SD * xxSD(2)**b_SD
ddGD(1) = a_GD * xxGD(1)**b_GD
ddGD(2) = a_GD * xxGD(2)**b_GD

call get_command_argument(1, length=i)
allocate(character(i) :: fname)
call get_command_argument(1, fname)
call get_command_argument(2, length=i)
allocate(character(i) :: fname2)
call get_command_argument(2, fname2)
call get_command_argument(3, length=i)
allocate(character(i) :: fname3)
call get_command_argument(3, fname3)
call get_command_argument(4, length=i)
allocate(character(i) :: fname4)
call get_command_argument(4, fname4)
call get_command_argument(5, length=i)
allocate(character(i) :: fname5)
call get_command_argument(5, fname5)
call get_command_argument(6, length=i)
allocate(character(i) :: fname6)
call get_command_argument(6, fname6)
call get_command_argument(7, length=i)
allocate(character(i) :: tmp)
call get_command_argument(7, tmp)
read(tmp, *) ze
deallocate(tmp)
call get_command_argument(8, length=i)
allocate(character(i) :: tmp)
call get_command_argument(8, tmp)
read(tmp, *) ncells
deallocate(tmp)
call get_command_argument(9, length=i)
allocate(character(i) :: tmp)
call get_command_argument(9, tmp)
read(tmp, *) l_ice
deallocate(tmp)

write(*,*) "start allocating"
allocate(rhocorr(ncells,ze))
allocate(rhocorr_c(ncells,ze))
allocate(rho(ncells,ze))
allocate(qc(ncells,ze))
allocate(qnc(ncells,ze))
allocate(qr(ncells,ze))
allocate(qnr(ncells,ze))
if (l_ice) then
  allocate(qi(ncells,ze))
  allocate(qni(ncells,ze))
  allocate(qs(ncells,ze))
  allocate(qns(ncells,ze))
  allocate(qg(ncells,ze))
  allocate(qng(ncells,ze))
endif
allocate(vm_c(ncells,ze))
allocate(vm_r(ncells,ze))
if (l_ice) then
  allocate(vm_i(ncells,ze))
  allocate(vm_s(ncells,ze))
  allocate(vm_g(ncells,ze))
  vm_i(:,:) = 0.0_wp
  vm_s(:,:) = 0.0_wp
  vm_g(:,:) = 0.0_wp
endif

write(*,*) "finished allocation"
vm_c(:,:) = 0.0_wp
vm_r(:,:) = 0.0_wp

!------------------------------------------------------------------
!calculate mass-weighted sedimentation speed for cloud droplets
!------------------------------------------------------------------
!read in files
!---------------------------------------------------------------
write(*,*) "start qc"

call read_in_2d(rho,ncells,ze,"rho",fname)
write(*,*) "read in rho"
call read_in_2d(qc,ncells,ze,"qc",fname)
write(*,*) "read in qc"
call read_in_2d(qnc,ncells,ze,"qnc",fname)
write(*,*) "read in qnc"
!calculate density correction factors
rhocorr = exp(-rho_vel*log(rho/rho0))
rhocorr_c = exp(-rho_vel_c*log(rho/rho0))
do k=1,ze
  do n=1,ncells
    !cloud optical depth:
      qc1=qc(n,k)*rho(n,k)
      qnc1=qnc(n,k)*rho(n,k)
    if (qc1 > 1.0e-20_wp .and. qnc1 > 1.0e-30_wp) then
       call calc_dist_facs(a_CD, b_CD, nuCD, muCD, qc1, qnc1,l, A)   
       call mass_weighted_v_array(xxCD, A, nuCD, l, muCD, avelCD, bvelCD, rhocorr_c(n,k), 1000, vm_c(n,k))
    endif
  enddo
enddo

deallocate(qc)
deallocate(qnc)
write(*,*) "did qc"
!------------------------------------------------------------------
!calculate optical depth for rain
!------------------------------------------------------------------
!read in files
!---------------------------------------------------------------
call read_in_2d(qr,ncells,ze,"qr",fname)
call read_in_2d(qnr,ncells,ze,"qnr",fname)
do k=1,ze
  do n=1,ncells
    !cloud optical depth:
      qr1=qr(n,k)*rho(n,k)
      qnr1=qnr(n,k)*rho(n,k)
    if (qr1 > 1.0e-20_wp .and. qnr1 > 1.0e-30_wp) then
      call calc_dist_facs(a_RD, b_RD, nuRD, muRD, qr1, qnr1,l, A)
      call mass_weighted_v_array(xxRD, A, nuRD, l, muRD, avelRD, bvelRD, rhocorr(n,k), 1000, vm_r(n,k))
    endif
  enddo
enddo

deallocate(qr)
deallocate(qnr)

write(*,*) "did qr"
if (l_ice) then
  !------------------------------------------------------------------
  !calculate optical depth for cloud ice
  !------------------------------------------------------------------
  !read in files
  !---------------------------------------------------------------
  call read_in_2d(qi,ncells,ze,"qi",fname)
  call read_in_2d(qni,ncells,ze,"qni",fname)
  do k=1,ze
    do n=1,ncells
      !cloud optical depth:
        qi1=qi(n,k)*rho(n,k)
        qni1=qni(n,k)*rho(n,k)
      if (qi1 > 1.0e-20_wp .and. qni1 > 1.0e-30_wp) then
        call calc_dist_facs(a_ID, b_ID, nuID, muID, qi1, qni1,l, A)
        call mass_weighted_v_array(xxID, A, nuID, l, muID, avelID, bvelID, rhocorr(n,k), 1000, vm_i(n,k))
      endif
    enddo
  enddo

  deallocate(qi)
  deallocate(qni)
  write(*,*) "did qi"
  !------------------------------------------------------------------
  !calculate optical depth for snow
  !------------------------------------------------------------------
  !read in files
  !---------------------------------------------------------------
  call read_in_2d(qs,ncells,ze,"qs",fname)
  call read_in_2d(qns,ncells,ze,"qns",fname)
  do k=1,ze
    do n=1,ncells
      !cloud optical depth:
        qs1=qs(n,k)*rho(n,k)
        qns1=qns(n,k)*rho(n,k)
      if (qs1 > 1.0e-20_wp .and. qns1 > 1.0e-30_wp) then
        call calc_dist_facs(a_SD, b_SD, nuSD, muSD, qs1, qns1,l, A)
        call mass_weighted_v_array(xxSD, A, nuSD, l, muSD, avelSD, bvelSD, rhocorr(n,k), 1000, vm_s(n,k))
      endif
    enddo
  enddo

  deallocate(qs)
  deallocate(qns)
  write(*,*) "did qs"

  !------------------------------------------------------------------
  !calculate optical depth for graupel
  !------------------------------------------------------------------
  !read in files
  !---------------------------------------------------------------
  call read_in_2d(qg,ncells,ze,"qg",fname)
  call read_in_2d(qng,ncells,ze,"qng",fname)
  do k=1,ze
    do n=1,ncells
      !cloud optical depth:
        qg1=qg(n,k)*rho(n,k)
        qng1=qng(n,k)*rho(n,k)
      if (qg1 > 1.0e-20_wp .and. qng1 > 1.0e-30_wp) then
        call calc_dist_facs(a_GD, b_GD, nuGD, muGD, qg1, qng1,l, A)
        call mass_weighted_v_array(xxGD, A, nuGD, l, muGD, avelGD, bvelGD, rhocorr(n,k), 1000, vm_g(n,k))
      endif
    enddo
  enddo

  deallocate(qg)
  deallocate(qng)
  write(*,*) "did qg"
endif
!------------------------------------------------------------------
!write output
!------------------------------------------------------------------
write(*,*) "write output"
call write_out_2d(vm_c,ncells,ze,"vmlc",fname2)
call write_out_2d(vm_r,ncells,ze,"vmlr",fname3)
if (l_ice) then
  call write_out_2d(vm_i,ncells,ze,"vmli",fname4)
  call write_out_2d(vm_s,ncells,ze,"vmls",fname5)
  call write_out_2d(vm_g,ncells,ze,"vmlg",fname6)
endif

CONTAINS

  DOUBLE PRECISION FUNCTION gfct(x)
    !*******************************************************************************
    !                                                                              *
    !       gamma function from Numerical Recipes (F77)
    !       *
    !       (slightly modified, accounting for inlinig and verctorization
    !       *
    !*******************************************************************************
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(in) :: x

    DOUBLE PRECISION :: tmp, p

    DOUBLE PRECISION, PARAMETER :: c1 =  76.18009173d0
    DOUBLE PRECISION, PARAMETER :: c2 = -86.50532033d0
    DOUBLE PRECISION, PARAMETER :: c3 =  24.01409822d0
    DOUBLE PRECISION, PARAMETER :: c4 = -1.231739516d0
    DOUBLE PRECISION, PARAMETER :: c5 =  0.120858003d-2
    DOUBLE PRECISION, PARAMETER :: c6 = -0.536382d-5
    DOUBLE PRECISION, PARAMETER :: stp = 2.50662827465d0

    tmp = x + 4.5d0;
    p = stp * (1d0 + c1/x + c2/(x+1d0) + c3/(x+2d0) + c4/(x+3d0) + c5/(x+4d0) +c6/(x+5d0))
    gfct = p * EXP( (x-0.5d0) * LOG(tmp) - tmp )

    RETURN
  END FUNCTION gfct

!------------------------------------------------------------------------
!subroutines for computing functions for optical depth
!------------------------------------------------------------------------
subroutine gamma_dist(x, A, nu, ll, mu, gamma_val)
    INTEGER, PARAMETER       ::                         &
       wp    = SELECTED_REAL_KIND (12,200)
    real(KIND=wp), intent(in) :: x, A, nu, ll, mu
    real(KIND=wp), intent(out) :: gamma_val

    gamma_val = A * (x**nu) * exp(-ll * (x**mu))

end subroutine gamma_dist

subroutine calc_dist_facs(a_CD,b_CD,nuCD,muCD,qc,nc, lCD, Acd)
   INTEGER, PARAMETER       ::                         &
       wp    = SELECTED_REAL_KIND (12,200)
    real(KIND=wp), intent(in) :: a_CD,b_CD,nuCD,muCD,qc,nc
    real(KIND=wp), intent(out) :: lCD, Acd
    real(KIND=wp) :: xmCD
    xmCD = qc/nc
    !lambda and A of size distribution (Seifert & Beheng 2006)
    lCD = (xmCD*gfct((nuCD+1)/muCD)/gfct((nuCD+2)/muCD))**(-muCD)
    Acd = muCD*nc*(lCD**((nuCD+1)/muCD))/gfct((nuCD+1)/muCD)

end subroutine calc_dist_facs

subroutine mass_weighted_v_array(xlim, A, nu, ll, mu, avel, bvel, rhocorr, nps, vel_int)
   INTEGER, PARAMETER       ::                         &
       wp    = SELECTED_REAL_KIND (12,200)
   real(KIND=wp), intent(in) :: xlim(2), A, nu, ll, mu, avel, bvel, rhocorr
   integer, intent (in)      ::nps
   real(KIND=wp), intent(out) :: vel_int
   real(KIND=wp)            :: x_bounds(nps), x_diff(nps-1), x_bin(nps-1), p(nps-1), step, x_int
   integer                  :: i
   step = (xlim(2) - xlim(1)) / real(nps - 1)  

    do i = 1, nps
        x_bounds(i) = xlim(1) / 2.0 + real(i - 1) * step
    end do

    do i = 1, nps-1
        x_diff(i) = x_bounds(i+1) - x_bounds(i)
        x_bin(i) = (x_bounds(i+1) + x_bounds(i)) / 2.0
        call gamma_dist(x_bin(i), A, nu, ll, mu, p(i))
    end do
    
    vel_int = sum (x_bin *p *avel*x_bin**bvel/ &
                   rhocorr *x_diff) 
                  
    x_int = sum (x_bin * p * x_diff)

    if (vel_int < 1.0e-6_wp .or. x_int < 1.0e-6_wp) then
      vel_int =0.0
    else
      vel_int =-1.0 *vel_int/x_int
    end if
end subroutine mass_weighted_v_array

subroutine transform_x2r (a_CD, b_CD, nuCD, muCD, qc, nc, Kcd, dltCD, psiCD, xiCD)
    INTEGER, PARAMETER       ::                          &
       wp    = SELECTED_REAL_KIND (12,200)
    real(KIND=wp), intent(in) :: a_CD, b_CD, nuCD, muCD, qc, nc
    real(KIND=wp), intent(out) :: Kcd, dltCD, psiCD, xiCD
    real(KIND=wp) :: xmCD, dmCD, lCD, Acd

    xmCD = real(qc, kind=wp) / real(nc, kind=wp)
    dmCD = a_CD * xmCD**b_CD

    ! lambda and A of size distribution (Seifert & Beheng 2006)
    lCD = (xmCD * gamma((nuCD+1)/muCD) / gamma((nuCD+2)/muCD))**(-muCD)
    Acd = muCD * nc * (lCD**((nuCD+1)/muCD)) / gamma((nuCD+1)/muCD)

    ! conversion constants to get to f(r)=K*r**psi*exp(-delta*r**xi)
    Kcd = (Acd / b_CD) * (2.0/a_CD)**((nuCD+1)/b_CD)
    dltCD = lCD * (2.0/a_CD)**(muCD/b_CD)
    psiCD = (nuCD+1) / b_CD - 1
    xiCD = muCD / b_CD

end subroutine transform_x2r

subroutine Ap_r (Dlim, K, psi, dlt, xi, nps, Ap_val)
   INTEGER, PARAMETER       ::                          &
       wp    = SELECTED_REAL_KIND (12,200)
    real(KIND=wp), intent(in) :: Dlim(2), K, psi, dlt, xi
    integer, intent(in)      :: nps  !number of points
    real(KIND=wp), intent(out) :: Ap_val
    real(KIND=wp)            :: r_bounds(nps), r_diff(nps-1), r_bin(nps-1), p(nps-1), step
    integer                  :: i

    step = (Dlim(2) - Dlim(1)) / (2.0 * real(nps - 1))  ! D/2 = Radius

    do i = 1, nps
        r_bounds(i) = Dlim(1) / 2.0 + real(i - 1) * step
    end do

    do i = 1, nps-1
        r_diff(i) = r_bounds(i+1) - r_bounds(i)
        r_bin(i) = (r_bounds(i+1) + r_bounds(i)) / 2.0
        ! Assuming gamma_dist is defined elsewhere
        call gamma_dist(r_bin(i), K, psi, dlt, xi, p(i))
    end do

    Ap_val = 2*sum(p * pi * r_bin**2 * r_diff)
end subroutine Ap_r

! -----------------------------------------------------------------------
!subroutines for reading and writing netcdf
! -----------------------------------------------------------------------
 SUBROUTINE read_in_1d (poutput,ke, varname, filename)

  IMPLICIT NONE
  character(*), intent(in) :: varname, filename
  integer, intent(in) :: ke
  real(kind=wp), dimension(ke), intent(inout) :: poutput
  integer :: ncId, rhVarId, status, nz
  integer, dimension(nf90_max_var_dims) :: dimIDs
  real(kind=wp), allocatable, dimension(:) ::  zvar
 

  status = nf90_open(filename, nf90_NoWrite, ncid)
  status = nf90_inq_varid(ncid,varname, rhVarId)
  status = nf90_inquire_variable(ncid, rhVarId, dimids = dimIDs)
  status = nf90_inquire_dimension(ncid, dimIDs(1), len = nz)

  allocate(zvar(nz))
  status = nf90_get_var(ncid, rhVarId, zvar)
  status = nf90_close(ncid)
  poutput = zvar
  deallocate(zvar)
  
  END SUBROUTINE read_in_1d

 SUBROUTINE read_in_2d (poutput,ie,je, varname, filename)

  IMPLICIT NONE
  character(*), intent(in) :: varname, filename
  integer, intent(in) :: ie, je
  real(kind=wp), dimension(ie,je), intent(inout) :: poutput
  integer :: ncId, rhVarId, status, nx, ny
  integer, dimension(nf90_max_var_dims) :: dimIDs
  real(kind=wp), allocatable, dimension(:,:) ::  zvar
 

  status = nf90_open(filename, nf90_NoWrite, ncid)
  status = nf90_inq_varid(ncid,varname, rhVarId)
  status = nf90_inquire_variable(ncid, rhVarId, dimids = dimIDs)
  status = nf90_inquire_dimension(ncid, dimIDs(1), len = nx)
  status = nf90_inquire_dimension(ncid, dimIDs(2), len = ny)

  allocate(zvar(nx,ny))
  !status = nf90_get_var(ncid, rhVarId, zvar(:,:), start=(/1,1,1/), count=(/nx,ny,ti/))
  status = nf90_get_var(ncid, rhVarId, zvar)
  status = nf90_close(ncid)
  poutput = zvar
  deallocate(zvar)
  
  END SUBROUTINE read_in_2d
! -----------------------------------------------------------------------

  SUBROUTINE read_in_3d (poutput,ie,je,ke,varname, filename)

  IMPLICIT NONE
  character(*), intent(in) :: varname, filename
  integer, intent(in) :: ie, je, ke 
  real(kind=wp), dimension(ie,je,ke), intent(inout) :: poutput
  integer :: ncId, rhVarId, status, nx, ny, nz
  integer, dimension(nf90_max_var_dims) :: dimIDs
  real(kind=wp), allocatable, dimension(:,:,:) ::  zvar

  status = nf90_open(filename, nf90_NoWrite, ncid)
  status = nf90_inq_varid(ncid,varname, rhVarId)
  status = nf90_inquire_variable(ncid, rhVarId, dimids = dimIDs)
  status = nf90_inquire_dimension(ncid, dimIDs(1), len = nx)
  status = nf90_inquire_dimension(ncid, dimIDs(2), len = ny)
  status = nf90_inquire_dimension(ncid, dimIDs(3), len = nz)

  allocate(zvar(nx,ny, nz))

  status = nf90_get_var(ncid, rhVarId, zvar)
  status = nf90_close(ncid)

  poutput = zvar

  deallocate(zvar)
  
  END SUBROUTINE read_in_3d
! -----------------------------------------------------------------------

  SUBROUTINE write_out_1d (poutput, ie, varname, filename)

  IMPLICIT NONE
  character(*), intent(in) :: varname, filename
  integer, intent(in) :: ie
  real(kind=wp), dimension(ie), intent(in) :: poutput
  integer :: ncid, rhVarId, status, xDimID

  status = nf90_create(filename, nf90_clobber, ncid)
  status = nf90_def_dim(ncid, 'ncells', ie, xDimID)
  status = nf90_def_var(ncid, varname, nf90_double, &
                            (/ xDimId /), RhVarId)
  status = nf90_enddef(ncid)
  status = nf90_put_var(ncid, rhVarId, poutput)
  status = nf90_close(ncid)

  END SUBROUTINE write_out_1d
! -----------------------------------------------------------------------

  SUBROUTINE write_out_2d (poutput, ie, je, varname, filename)

  IMPLICIT NONE
  character(*), intent(in) :: varname, filename
  integer, intent(in) :: ie, je
  real(kind=wp), dimension(ie,je), intent(in) :: poutput
  integer :: ncid, rhVarId, status, xDimID, yDimID


  status = nf90_create(filename, nf90_clobber, ncid)
  status = nf90_def_dim(ncid, 'lon', ie, xDimID)
  status = nf90_def_dim(ncid, 'lat', je, yDimID)
  status = nf90_def_var(ncid, varname, nf90_double, &
                            (/ xDimId, yDimID /), RhVarId)
  status = nf90_enddef(ncid)
  status = nf90_put_var(ncid, rhVarId, poutput)
  status = nf90_close(ncid)

  END SUBROUTINE write_out_2d
! -----------------------------------------------------------------------

  SUBROUTINE write_out_3d (poutput, ie, je, ke, varname, filename)

  IMPLICIT NONE
  character(*), intent(in) :: varname, filename
  integer, intent(in) :: ie, je, ke
  real(kind=wp), dimension(ie,je,ke), intent(in) :: poutput
  integer :: ncId, rhVarId, status, xDimID, yDimID, zDimID

  status = nf90_create(filename, nf90_clobber, ncid)
  status = nf90_def_dim(ncid, 'lon', ie, xDimID)
  status = nf90_def_dim(ncid, 'lat', je, yDimID)
  status = nf90_def_dim(ncid, 'lev', ke, zDimID)
  status = nf90_def_var(ncid, varname, nf90_double, &
                            (/ xDimId, yDimID, zDimID /), RhVarId)
  status = nf90_enddef(ncid)
  status = nf90_put_var(ncid, rhVarId, poutput)
  status = nf90_close(ncid)

  END SUBROUTINE write_out_3d
! ---------------------------------------------------------------------- 


END PROGRAM compute_vsedi
