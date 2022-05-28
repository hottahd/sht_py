subroutine forward(qq,y,z,jx,kx,fqq)
  implicit none

  integer :: j,k,m,l
  real :: dy, dz
  real :: fm, fl ! real of m, l
  real :: faca, facb
  real, parameter :: pi = 3.141592653589e0
  integer, intent(in) :: jx,kx
  real, dimension(0:jx-1) :: pm, pn, pn0, pm0, pm1, pm2 ! legendre function
  real, dimension(0:jx-1), intent(in) :: y  
  real, dimension(0:kx-1), intent(in) :: z
  complex, dimension(0:jx-1,0:kx-1), intent(in)  :: qq
  complex, dimension(0:jx-1,0:kx-1), intent(out) :: fqq

  dy = y(1) - y(0)
  dz = z(1) - z(0)
  
  do k = 0,kx-1
  do j = 0,jx-1
     fqq(j,k) = (0.e0,0.e0)
  enddo
  enddo
  
  do m = 0,kx/2-1
     fm = real(m)
     ! l = m case
     if (m == 0) then
        do j = 0,jx-1
           pm(j) = 1.e0/sqrt(4.e0*pi)
        enddo
     else
        do j = 0,jx-1
           pm(j) = -sqrt((2.0*fm+1.0)/(2.0*fm))*sin(y(j))*pm(j)
        enddo
        do j = 0,jx-1
           if (abs(pm(j)) < 1.e-35) then
              pm(j) = 0.e0
           endif
        enddo
     endif

     ! negative m
     do j = 0,jx-1
        pn(j) = (-1.e0)**fm*pm(j)
     enddo

     ! Integration
     do j = 0,jx-1
        fqq(m,m) = fqq(m,m) + qq(j,m)*pm(j)*sin(y(j))*dy
     enddo
        
     if (m /= 0) then
        do j = 0,jx-1
           fqq(m,kx-m) = fqq(m,kx-m) + qq(j,kx-m)*pn(j)*sin(y(j))*dy
        enddo
     endif

     do j = 0,jx-1
        pm1(j) = pm(j)
        pm2(j) = 0.e0
     enddo

     do l = m+1,kx/2-1
        fl = real(l)

        faca = sqrt((2.0*fl-1.0)*(2.0*fl+1.0)/(fl+fm)/(fl-fm))
        facb = sqrt((2.0*fl+1.0)/(2.0*fl-3.0)*(fl+fm-1.0)*(fl-fm-1.0)/(fl+fm)/(fl-fm))

        do j = 0,jx-1
           pm0(j) = faca*cos(y(j))*pm1(j) - facb*pm2(j)
           ! Integration
           fqq(l,m) = fqq(l,m) + qq(j,m)*pm0(j)*sin(y(j))*dy
        enddo

        ! negative m        
        if (m /= 0) then
           do j = 0,jx-1
              pn0(j) = (-1.0)**fm*pm0(j)
              fqq(l,kx-m) = fqq(l,kx-m) + qq(j,kx-m)*pn0(j)*sin(y(j))*dy
           enddo
        endif
                
        do j = 0,jx-1        
           pm2(j) = pm1(j)
           pm1(j) = pm0(j)
        enddo
     enddo
  enddo
     
  return
end subroutine forward

subroutine backward(qq,y,z,jx,kx,fqq)
  implicit none

  integer :: j,k,m,l
  real :: dy, dz
  real :: fm, fl ! real of m, l
  real :: faca, facb
  real, parameter :: pi = 3.141592653589e0
  integer, intent(in) :: jx,kx
  real, dimension(0:jx-1) :: pm, pn, pn0, pm0, pm1, pm2 ! legendre function
  real, dimension(0:jx-1), intent(in) :: y  
  real, dimension(0:kx-1), intent(in) :: z
  complex, dimension(0:jx-1,0:kx-1), intent(in)  :: qq
  complex, dimension(0:jx-1,0:kx-1), intent(out) :: fqq

  dy = y(1) - y(0)
  dz = z(1) - z(0)
  
  do k = 0,kx-1
  do j = 0,jx-1
     fqq(j,k) = (0.e0,0.e0)
  enddo
  enddo
  
  do m = 0,kx/2-1
     fm = real(m)
     ! l = m case
     if (m == 0) then
        do j = 0,jx-1
           pm(j) = 1.e0/sqrt(4.e0*pi)
        enddo
     else
        do j = 0,jx-1
           pm(j) = -sqrt((2.0*fm+1.0)/(2.0*fm))*sin(y(j))*pm(j)
        enddo
        do j = 0,jx-1
           if (abs(pm(j)) < 1.e-35) then
              pm(j) = 0.e0
           endif
        enddo
     endif

     ! negative m
     do j = 0,jx-1
        pn(j) = (-1.e0)**fm*pm(j)
     enddo

     ! Integration
     do j = 0,jx-1
        fqq(j,m) = fqq(j,m) + qq(m,m)*pm(j)
     enddo
     if (m /= 0) then
        do j = 0,jx-1
           fqq(j,kx-m) = fqq(j,kx-m) + qq(m,kx-m)*pn(j)
        enddo
     endif
     
     do j = 0,jx-1
        pm1(j) = pm(j)
        pm2(j) = 0.e0
     enddo

     do l = m+1,kx/2-1
        fl = real(l)

        faca = sqrt((2.0*fl-1.0)*(2.0*fl+1.0)/(fl+fm)/(fl-fm))
        facb = sqrt((2.0*fl+1.0)/(2.0*fl-3.0)*(fl+fm-1.0)*(fl-fm-1.0)/(fl+fm)/(fl-fm))

        do j = 0,jx-1
           pm0(j) = faca*cos(y(j))*pm1(j) - facb*pm2(j)
           ! Integration
           fqq(j,m) = fqq(j,m) + qq(l,m)*pm0(j)
        enddo

        ! negative m        
        if (m /= 0) then
           do j = 0,jx-1
              pn0(j) = (-1.0)**fm*pm0(j)
              fqq(j,kx-m) = fqq(j,kx-m) + qq(l,kx-m)*pn0(j)
           enddo
        endif
                
        do j = 0,jx-1
           pm2(j) = pm1(j)
           pm1(j) = pm0(j)
        enddo        
     enddo
  enddo
     
  return
end subroutine backward
