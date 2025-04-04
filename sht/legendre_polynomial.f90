subroutine up(l, jx, y, pm1, pm2, pm)
    implicit none

    integer :: j
    double precision :: fl
    integer, intent(in) :: l ! order of Legendre polynomial
    integer, intent(in) :: jx

    
    double precision, dimension(0:jx-1), intent(in) :: y ! colatitude
    double precision, dimension(0:jx-1), intent(in) :: pm1, pm2
    double precision, dimension(0:jx-1), intent(out) :: pm

    fl = dble(l)
    do j = 0, jx-1
        pm(j) = + sqrt((2.d0*fl+1.d0)*(2.d0*fl-1.d0))/fl*cos(y(j))*pm1(j) &
                - (fl-1.d0)/fl*sqrt((2.d0*fl+1.d0)/abs(2.d0*fl-3.d0))*pm2(j)
    enddo

    return
end subroutine up
!##################################################################

subroutine polynomial(y, jx, l_target, pm) bind(C)
    use iso_fortran_env
    implicit none

    integer :: j, l
    integer, intent(in) :: jx 
    integer, intent(in) :: l_target
    double precision, dimension(0:jx-1) :: pm1, pm2

    double precision, dimension(0:jx-1), intent(in) :: y
    double precision, dimension(0:jx-1), intent(out) :: pm

    ! l = 0
    do j = 0, jx-1
        pm(j) = 1.d0/sqrt(2.d0)
        pm1(j) = pm(j)
        pm2(j) = 0.d0
    enddo

    if (l_target == 0) then
        return
    end if

    do l = 1, l_target
        call up(l, jx, y, pm1, pm2, pm)

        do j = 0, jx-1
            pm2(j) = pm1(j)
            pm1(j) = pm(j)
        enddo
    enddo

    return
end subroutine polynomial

!##################################################################

subroutine forward(qq, y, jx, fqq) bind(C)
    use iso_fortran_env
    implicit none

    integer :: j, l
    double precision :: dy

    integer, intent(in) :: jx
    double precision, dimension(0:jx-1), intent(in) :: qq ! order of Legendre polynomial
    double precision, dimension(0:jx-1), intent(in) :: y ! colatitude
    double precision, dimension(0:jx-1), intent(out) :: fqq

    double precision, dimension(0:jx-1) :: pm, pm1, pm2

    ! l = 0
    do j = 0, jx-1
        pm(j) = 1.d0/sqrt(2.d0)
        pm1(j) = pm(j)
        pm2(j) = 0.d0
    enddo

    dy = y(1) - y(0)

    do j = 0, jx-1
        fqq(0) = fqq(0) + qq(j)*sin(y(j))*dy*pm(j)
    enddo

    do l = 1, jx - 1
        call up(l, jx, y, pm1, pm2, pm)
        do j = 0, jx - 1
            fqq(l) = fqq(l) + qq(j)*sin(y(j))*dy*pm(j)
        enddo

        do j = 0, jx-1
            pm2(j) = pm1(j)
            pm1(j) = pm(j)
        enddo
    enddo

    return
end subroutine forward

!##################################################################

subroutine backward(fqq, y, jx, qq) bind(C)
    use iso_fortran_env
    implicit none
    integer :: j, l

    integer, intent(in) :: jx
    double precision, dimension(0:jx-1), intent(in) :: fqq ! order of Legendre polynomial
    double precision, dimension(0:jx-1), intent(in) :: y ! colatitude
    double precision, dimension(0:jx-1), intent(out) :: qq

    double precision, dimension(0:jx-1) :: pm, pm1, pm2

    ! l = 0
    do j = 0, jx-1
        pm(j) = 1.d0/sqrt(2.d0)
        pm1(j) = pm(j)
        pm2(j) = 0.d0
    enddo

    do j = 0, jx - 1
        qq(j) = qq(j) + fqq(0)*pm(j)
    enddo

    do l = 1, jx - 1
        call up(l, jx, y, pm1, pm2, pm)
        do j = 0, jx - 1
            qq(j) = qq(j) + fqq(l)*pm(j)
        enddo

        do j = 0, jx-1
            pm2(j) = pm1(j)
            pm1(j) = pm(j)
        enddo
    enddo

    return
end subroutine backward
!##################################################################
