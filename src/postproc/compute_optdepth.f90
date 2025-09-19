! Compute optical depth.
PROGRAM compute_optdepth
  USE netcdf
  !USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE

  IMPLICIT NONE

  INTEGER, PARAMETER :: &
    wp = SELECTED_REAL_KIND(12, 200), &
    integers = KIND(1)

  ! Collection of PARAMETERs of size distributions.
  REAL(KIND=wp), PARAMETER :: &
    !cloud = cloud_nue1mue1:
    a_CD = .124, & ! m/kg^beta, D = ax^b
    b_CD = 1./3., &
    nuCD = 1., & ! nu of f(x) cloud
    muCD = 1., & ! mu "  " cloud

    ! rain = rainSBB:
    a_RD = 0.124 , &
    b_RD = 1./3., &
    nuRD = 1.0, &
    muRD = 1./3., &

    ! ice = ice_cosmo5:
    a_ID = 0.835, &
    b_ID = 0.39, &
    nuID = 0.0, &
    muID = 1./3., &

    ! snow= snowSBB:
    a_SD = 5.13, &
    b_SD = 0.5, &
    nuSD = 0.0, &
    muSD = 0.5, &

    ! graupel = graupelhail_cosmo5:
    a_GD = 0.142, &
    b_GD = 0.314, &
    nuGD = 1.0, &
    muGD = 1./3.

  ! Number of bins for dr bins in size distribution.
  INTEGER(KIND=integers) :: nbins = 1000

  REAL(kind=wp), PARAMETER :: pi = acos(-1.0_wp)

  REAL(KIND = wp) :: xxCD(2), xxRD(2), xxID(2), xxSD(2), xxGD(2), &
    ddCD(2), ddRD(2), ddID(2), ddSD(2), ddGD(2), &
    Kconst, dlt, psi, xi

  INTEGER(KIND=integers) :: ze, ncells, n, k, ierr

  REAL(KIND=wp) :: qc1, qnc1, apr, qr1, qnr1, qg1, qng1, qs1, qns1, &
    qi1, qni1, opt_tmp

  CHARACTER(:), ALLOCATABLE :: tmp
  CHARACTER(:), ALLOCATABLE :: fname, fname2, fname3, fname4
  REAL(KIND=wp), ALLOCATABLE :: h(:,:), dh(:,:), rho(:,:), qr(:,:), qc(:,:), &
    qi(:,:), qs(:,:), qg(:,:), qnr(:,:), qnc(:,:), &
    qni(:,:), qns(:,:), qng(:,:), &
    opt_cld(:,:), & ! Cloud optical depth
    opt_3d(:,:), & ! Optical depth 3D
    clt_2d(:), & ! Cloud fraction 2D
    opt_2d(:), & ! Optical depth 2D
    opt_2d_cld(:) ! Cloud optical depth 2D

  LOGICAL :: l_ice

  INTEGER :: i
  INTEGER :: argc

  ! Define size ranges [x_min, x_max].
  xxCD = [4.2e-15, 2.6e-10] ! kg
  xxRD = [2.6e-10, 6.5e-5] ! kg
  xxID = [1.0e-12, 1.0e-5] ! kg
  xxSD = [1.0e-10, 2.0e-5] ! kg
  xxGD = [4.19e-09, 5.3e-4] ! kg

  ! Calculate dd values.
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

  CALL get_command_argument(1, length=i)
  allocate(CHARACTER(i) :: fname)
  CALL get_command_argument(1, fname)
  CALL get_command_argument(2, length=i)
  allocate(CHARACTER(i) :: fname2)
  CALL get_command_argument(2, fname2)
  CALL get_command_argument(3, length=i)
  allocate(CHARACTER(i) :: fname3)
  CALL get_command_argument(3, fname3)
  CALL get_command_argument(4, length=i)
  allocate(CHARACTER(i) :: fname4)
  CALL get_command_argument(4, fname4)
  CALL get_command_argument(5, length=i)
  allocate(CHARACTER(i) :: tmp)
  CALL get_command_argument(5, tmp)
  READ(tmp, *) l_ice
  DEALLOCATE(tmp)

  !
  ! Calculate optical depth of cloud liquid.
  !
  WRITE(*,"(A)") "Computing cloud liquid optical depth..."
  CALL read_2d(h, "z_ifc", fname)
  CALL read_2d(rho, "rho", fname)
  CALL read_2d(qc, "qc", fname)
  CALL read_2d(qnc, "qnc", fname)

  ncells = SIZE(rho, 1)
  ze = SIZE(rho, 2)

  allocate(dh(ncells,ze))
  allocate(opt_3d(ncells,ze))
  allocate(opt_2d(ncells))
  allocate(clt_2d(ncells))
  allocate(opt_2d_cld(ncells))
  allocate(opt_cld(ncells,ze))

  opt_3d(:,:) = 0.0_wp
  opt_cld(:,:) = 0.0_wp
  clt_2d(:) = 0.0_wp
  dh = ABS(h(:,2:) - h(:,1:ze))

  DO k=1,ze
    DO n=1,ncells
      ! cloud optical depth:
      qc1 = qc(n,k)*rho(n,k)
      qnc1 = qnc(n,k)*rho(n,k)
      IF (qc1 > 1.0e-20_wp .AND. qnc1 > 1.0e-8_wp) THEN
        CALL transform_x2r(a_CD, b_CD, nuCD, muCD, qc1, qnc1, Kconst, dlt, psi, xi)
        CALL Ap_r(ddCD, Kconst, psi, dlt, xi, nbins, opt_cld(n,k))
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(qc)
  DEALLOCATE(qnc)

  !
  ! Calculate optical depth of rain.
  !
  WRITE(*,"(A)") "Computing rain optical depth..."
  CALL read_2d(qr, "qr", fname)
  CALL read_2d(qnr, "qnr", fname)
  DO k=1,ze
    DO n=1,ncells
      opt_tmp = 0.0
      ! cloud optical depth:
      qr1 = qr(n,k)*rho(n,k)
      qnr1 = qnr(n,k)*rho(n,k)
      IF (qr1 > 1.0e-20_wp .AND. qnr1 > 1.0e-10_wp) THEN
        CALL transform_x2r(a_RD, b_RD, nuRD, muRD, qr1, qnr1, Kconst, dlt, psi, xi)
        CALL Ap_r(ddRD, Kconst, psi, dlt, xi, nbins, opt_tmp)
        opt_3d(n,k) = opt_3d(n,k) + opt_tmp
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(qr)
  DEALLOCATE(qnr)

  IF (l_ice) THEN
    !
    ! Calculate optical depth of cloud ice.
    !
    WRITE(*,"(A)") "Computing cloud ice optical depth..."
    CALL read_2d(qi, "qi", fname)
    CALL read_2d(qni, "qni", fname)
    DO k=1,ze
      DO n=1,ncells
        ! Cloud optical depth:
        qi1 = qi(n,k)*rho(n,k)
        qni1 = qni(n,k)*rho(n,k)
        IF (qi1 > 1.0e-20_wp .AND. qni1 > 1.0e-10_wp) THEN
          opt_tmp = 0.0
          CALL transform_x2r(a_ID, b_ID, nuID, muID, qi1, qni1, Kconst, dlt, psi, xi)
          CALL Ap_r(ddID, Kconst, psi, dlt, xi, nbins, opt_tmp)
          opt_3d(n,k) = opt_3d(n,k) + opt_tmp
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(qi)
    DEALLOCATE(qni)

    !
    ! Calculate optical depth of snow.
    !
    WRITE(*,"(A)") "Computing snow optical depth..."
    CALL read_2d(qs, "qs", fname)
    CALL read_2d(qns, "qns", fname)
    DO k=1,ze
      DO n=1,ncells
        opt_tmp = 0.0
        ! cloud optical depth:
        qs1 = qs(n,k)*rho(n,k)
        qns1 = qns(n,k)*rho(n,k)
        IF (qs1 > 1.0e-20_wp .AND. qns1 > 1.0e-10_wp) THEN
          CALL transform_x2r(a_SD, b_SD, nuSD, muSD, qs1, qns1, Kconst, dlt, psi, xi)
          CALL Ap_r(ddSD, Kconst, psi, dlt, xi, nbins, opt_tmp)
          opt_3d(n,k) = opt_3d(n,k) + opt_tmp
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(qs)
    DEALLOCATE(qns)

    !
    ! Calculate optical depth of graupel
    !
    WRITE(*,"(A)") "Computing graupel optical depth..."
    CALL read_2d(qg, "qg", fname)
    CALL read_2d(qng, "qng", fname)
    DO k=1,ze
      DO n=1,ncells
        opt_tmp=0.0
        ! cloud optical depth:
        qg1 = qg(n,k)*rho(n,k)
        qng1 = qng(n,k)*rho(n,k)
        IF (qg1 > 1.0e-20_wp .AND. qng1 > 1.0e-10_wp) THEN
          CALL transform_x2r(a_GD, b_GD, nuGD, muGD, qg1, qng1, Kconst, dlt, psi, xi)
          CALL Ap_r(ddGD, Kconst, psi, dlt, xi, nbins, opt_tmp)
          opt_3d(n,k) = opt_3d(n,k) + opt_tmp
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(qg)
    DEALLOCATE(qng)
  ENDIF

  opt_2d_cld = sum(opt_cld*dh, dim=2)
  opt_2d = sum(opt_3d*dh, dim=2)
  opt_2d = opt_2d + opt_2d_cld
  DO n=1,ncells
    IF (opt_2d(n) .GE. 2.0_wp) THEN
      clt_2d(n) = 1.0_wp
    ENDIF
  END DO

  WRITE(*,"(A)") "Writing output..."
  CALL write_1d(opt_2d_cld, "odlc", fname2)
  CALL write_1d(opt_2d, "od", fname3)
  CALL write_1d(clt_2d, "clt", fname4)

CONTAINS

  !
  ! Subroutines for computing functions for optical depth.
  !

  SUBROUTINE gamma_dist(x, A, nu, ll, mu, gamma_val)
    INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12, 200)
    REAL(KIND=wp), INTENT(in) :: x, A, nu, ll, mu
    REAL(KIND=wp), INTENT(out) :: gamma_val

    gamma_val = A * (x**nu) * exp(-ll * (x**mu))
  END SUBROUTINE gamma_dist

  SUBROUTINE transform_x2r (a_CD, b_CD, nuCD, muCD, qc, nc, Kcd, dltCD, psiCD, xiCD)
    INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12, 200)
    REAL(KIND=wp), INTENT(in) :: a_CD, b_CD, nuCD, muCD, qc, nc
    REAL(KIND=wp), INTENT(out) :: Kcd, dltCD, psiCD, xiCD
    REAL(KIND=wp) :: xmCD, dmCD, lCD, Acd

    xmCD = REAL(qc, kind=wp) / REAL(nc, kind=wp)
    dmCD = a_CD * xmCD**b_CD

    ! lambda and A of size distribution (Seifert & Beheng 2006)
    lCD = (xmCD * gamma((nuCD+1)/muCD) / gamma((nuCD+2)/muCD))**(-muCD)
    Acd = muCD * nc * (lCD**((nuCD+1)/muCD)) / gamma((nuCD+1)/muCD)

    ! conversion constants to get to f(r)=K*r**psi*exp(-delta*r**xi)
    Kcd = (Acd / b_CD) * (2.0/a_CD)**((nuCD+1)/b_CD)
    dltCD = lCD * (2.0/a_CD)**(muCD/b_CD)
    psiCD = (nuCD+1) / b_CD - 1
    xiCD = muCD / b_CD

  END SUBROUTINE transform_x2r

  SUBROUTINE Ap_r(Dlim, K, psi, dlt, xi, nps, Ap_val)
    INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12, 200)
    REAL(KIND=wp), INTENT(in) :: Dlim(2), K, psi, dlt, xi
    INTEGER, INTENT(in) :: nps ! number of points
    REAL(KIND=wp), INTENT(out) :: Ap_val
    REAL(KIND=wp) :: r_bounds(nps), r_diff(nps-1), r_bin(nps-1), p(nps-1), step
    INTEGER :: i

    step = (Dlim(2) - Dlim(1)) / (2.0 * REAL(nps - 1)) ! D/2 = Radius

    DO i = 1, nps
      r_bounds(i) = Dlim(1) / 2.0 + REAL(i - 1) * step
    END DO

    DO i = 1, nps-1
      r_diff(i) = r_bounds(i+1) - r_bounds(i)
      r_bin(i) = (r_bounds(i+1) + r_bounds(i)) / 2.0
      ! Assuming gamma_dist is defined elsewhere
      CALL gamma_dist(r_bin(i), K, psi, dlt, xi, p(i))
    END DO

    Ap_val = 2*sum(p * pi * r_bin**2 * r_diff)
  END SUBROUTINE Ap_r

  !
  ! Subroutines for reading and writing NetCDF.
  !

  SUBROUTINE read_2d(poutput, varname, filename)
    IMPLICIT NONE
    REAL(kind=wp), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) ::  poutput
    CHARACTER(*), INTENT(in) :: varname, filename
    INTEGER, DIMENSION(nf90_max_var_dims) :: dimIDs
    INTEGER :: ncId, rhVarId, status, nx, ny

    status = nf90_open(filename, nf90_nowrite, ncid)
    status = nf90_inq_varid(ncid,varname, rhVarId)
    status = nf90_inquire_variable(ncid, rhVarId, dimids=dimIDs)
    status = nf90_inquire_dimension(ncid, dimIDs(1), len=nx)
    status = nf90_inquire_dimension(ncid, dimIDs(2), len=ny)

    allocate(poutput(nx,ny))
    !status = nf90_get_var(ncid, rhVarId, zvar(:,:), start=(/1,1,1/), count=(/nx,ny,ti/))
    status = nf90_get_var(ncid, rhVarId, poutput)
    status = nf90_close(ncid)
  END SUBROUTINE read_2d

  SUBROUTINE write_1d(poutput, varname, filename)
    IMPLICIT NONE
    CHARACTER(*), INTENT(in) :: varname, filename
    REAL(kind=wp), DIMENSION(:), INTENT(in) :: poutput
    INTEGER :: ncid, rhVarId, status, xDimID

    status = nf90_create(filename, nf90_clobber, ncid)
    status = nf90_def_dim(ncid, 'ncells', SIZE(poutput), xDimID)
    status = nf90_def_var(ncid, varname, nf90_double, (/ xDimId /), RhVarId)
    status = nf90_enddef(ncid)
    status = nf90_put_var(ncid, rhVarId, poutput)
    status = nf90_close(ncid)
  END SUBROUTINE write_1d
END PROGRAM compute_optdepth
