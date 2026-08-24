module module_bl_mynnedmf_diags
  use module_bl_mynnedmf_common, only: kind_phys,grav

  implicit none
  real(kind_phys),parameter::zero=0.0
  
  private
  public :: cloud_water_path, wspd_at_hgts !, pblh, cloud_ceiling

  contains

!===================================================================
! Subroutine to calculate LWP, IWP, and SWP
!===================================================================
  subroutine cloud_water_path (kts, kte, p1, qc1, qi1, qs1, qc_bl1, qi_bl1,      &
                cldfra_bl1, lwp1, iwp1, swp1)
    implicit none

    integer, intent(in) :: kts, kte
    real(kind_phys), dimension(kts:kte), intent(in) :: p1, qc1, qi1, qs1
    real(kind_phys), dimension(kts:kte), intent(in), optional :: qc_bl1, qi_bl1, &
                                                      cldfra_bl1 
    real(kind_phys), intent(inout), optional :: lwp1, iwp1, swp1

    ! local variables
    real(kind_phys) :: dp, sum1, sum2, sum3, qctotal, qitotal, qstotal
    integer :: k

       
    sum1=zero
    sum2=zero
    sum3=zero

    do k = kts, kte-1
      dp = p1(k) - p1(k+1)
      if (present(qc_bl1) .and. present(cldfra_bl1) .and. qc1(k)<1e-6 .and.      &
        cldfra_bl1(k)>0.01) then
        qctotal = qc_bl1(k)
      else
        qctotal = qc1(k)
      endif

      if (present(qi_bl1) .and. present(cldfra_bl1) .and. qi1(k)<1e-9 .and.      &
        cldfra_bl1(k)>0.01) then
        qitotal = qi_bl1(k)
      else
        qitotal = qi1(k)
      endif

      qstotal = qs1(k) !there currently is no sgs snow, it's ice+snow...

      sum1 = sum1 + max((dp/grav) * qctotal, zero)
      sum2 = sum2 + max((dp/grav) * (qitotal+qstotal), zero) !actually frozen water path
      sum3 = sum3 + max((dp/grav) * qstotal, zero)
    enddo

    lwp1 = sum1 * 1000._kind_phys ! kg m-2  --> g m-2
    iwp1 = sum2 * 1000._kind_phys
    swp1 = sum3 * 1000._kind_phys
    print*, 'lwp', lwp1

  end subroutine cloud_water_path

!===================================================================
! Subroutine to calculate wind speed at specific heights
!===================================================================
  subroutine wspd_at_hgts (kts, kte, dz1, u1, v1, wspd101, wspd801, wspd1601)
    implicit none

    integer, intent(in) :: kts, kte
    real(kind_phys), dimension(kts:kte), intent(in) :: dz1, u1, v1
    real(kind_phys), dimension(kts:kte) :: zw1
    real(kind_phys) :: depth, z_agl0, z_agl

    real(kind_phys), intent(out) :: wspd101, wspd801, wspd1601
    integer :: k

    zw1(kts)  = zero
    z_agl     = zero
    wspd101   = -99._kind_phys
    wspd801   = -99._kind_phys
    wspd1601  = -99._kind_phys
    do k = kts, kte-1
      z_agl0 = z_agl
      z_agl = zw1(k) + 0.5 * dz1(k)  ! at mass point AGL
      zw1(k+1) = zw1(k) + dz1(k)     ! at interface
      depth = z_agl - z_agl0
      !print*, 'depth', depth,'z_agl',z_agl, 'wspd', sqrt(u1(k)**2+v1(k)**2)

      ! 10-m wind speed
      if (z_agl > 10. .and. wspd101 == -99.) then
        wspd101 = interpolate_wind (k, kts, 10.0, z_agl, u1(k), v1(k),            &
                    u1(k-1), v1(k-1), depth)
        print*, 'k', k, 'wspd10', wspd101
      endif

      ! 80-m wind speed
      if (z_agl > 80. .and. wspd801 == -99.) then
        wspd801 = interpolate_wind (k, kts, 80.0, z_agl, u1(k), v1(k),            &
                    u1(k-1), v1(k-1), depth)
        print*, 'k', k, 'wspd80', wspd801
      endif

      ! 160-m wind speed
      if (z_agl > 160. .and. wspd1601 == -99.) then
        wspd1601 = interpolate_wind (k, kts, 160.0, z_agl, u1(k), v1(k),          &
                    u1(k-1), v1(k-1), depth)
        print*, 'k', k, 'wspd160', wspd1601
      endif

    enddo

  end subroutine wspd_at_hgts


  function interpolate_wind(k_idx, kts, target_z, z_agl,                      &
                            u1_curr,v1_curr, u1_prev, v1_prev, depth)         &
                              result(wspd1_hgt)

    integer, intent(in) :: k_idx, kts
    real, intent(in)    :: target_z, z_agl, u1_curr, v1_curr, u1_prev,        &
            v1_prev, depth
    real :: wspd1_hgt

    real :: wgt, wspd1_curr, wspd1_prev

    wspd1_curr = sqrt(u1_curr**2 + v1_curr**2)

    if (k_idx == kts) then
      wspd1_hgt = wspd1_curr
    else
      wgt = (z_agl - target_z) /depth
      wgt = min(max(zero, wgt), 1.0_kind_phys)
      wspd1_prev = sqrt(u1_prev**2 + v1_prev**2)
      wspd1_hgt = (1.0_kind_phys - wgt) * wspd1_curr + wgt * wspd1_prev
    endif

  end function interpolate_wind

!=================================================================================================================
 end module module_bl_mynnedmf_diags
!=================================================================================================================
