!##################################################################
! calculate associated Legendre function P_m^m from P_{m-1}^{m-1}
subroutine legendre_m_up(m,y,pm,jx,pn)
  implicit none

  integer :: j
  integer, intent(in) :: m,jx
  double precision :: fm
  double precision, dimension(0:jx-1), intent(inout) :: pm
  double precision, dimension(0:jx-1), intent(out) :: pn
  double precision, dimension(0:jx-1), intent(in) :: y
  
  fm = dble(m)
  ! l = m case
  if (m == 0) then
     do j = 0,jx-1
        pm(j) = 1.d0
     enddo
  else
     do j = 0,jx-1
        pm(j) = -sqrt((2.d0*fm+1.d0)/(2.d0*fm))*sin(y(j))*pm(j)
     enddo
     do j = 0,jx-1
        if (abs(pm(j)) < 1.d-300) then
           pm(j) = 0.d0
        endif
     enddo
  endif

  ! negative m
  do j = 0,jx-1
     pn(j) = (-1.d0)**fm*pm(j)
  enddo

  return
end subroutine legendre_m_up

!##################################################################

subroutine legendre_l_up(l,m,y,pm1,pm2,jx,pm0,pn0)
  implicit none

  integer :: j
  integer, intent(in) :: l,m,jx
  double precision :: fm,fl,faca,facb
  double precision, dimension(0:jx-1), intent(in) :: y,pm1,pm2
  double precision, dimension(0:jx-1), intent(out) :: pm0,pn0

  fm = dble(m)
  fl = dble(l)

  faca = sqrt((2.d0*fl-1.d0)*(2.d0*fl+1.d0)/(fl+fm)/(fl-fm))
  facb = sqrt((2.d0*fl+1.d0)/(2.d0*fl-3.d0)*(fl+fm-1.d0)*(fl-fm-1.d0)/(fl+fm)/(fl-fm))

  do j = 0,jx-1
     pm0(j) = faca*cos(y(j))*pm1(j) - facb*pm2(j)
  enddo

  if (m /= 0) then
     do j = 0,jx-1
        pn0(j) = (-1.0)**fm*pm0(j)        
     enddo
  endif
  
  return
end subroutine legendre_l_up

!##################################################################
subroutine forward(N,qq,y,z,jx,kx,fqq)
  implicit none

  integer :: j,k,m,l
  integer, intent(in) :: N,jx,kx
  double precision :: dy, dz
  double precision, dimension(0:jx-1) :: pm, pn, pn0, pm0, pm1, pm2 ! legendre function
  double precision, dimension(0:jx-1), intent(in) :: y  
  double precision, dimension(0:kx-1), intent(in) :: z
  double complex, dimension(0:jx-1,0:kx-1), intent(in)  :: qq
  double complex, dimension(0:N*kx/2-1,0:kx-1), intent(out) :: fqq
  
  dy = y(1) - y(0)
  dz = z(1) - z(0)
  
  do k = 0,kx-1
  do j = 0,N*kx/2-1
     fqq(j,k) = (0.d0,0.d0)
  enddo
  enddo
  
  do m = 0,N*kx/2-1
     call legendre_m_up(m,y,pm,jx,pn)
     
     if (mod(m,n) == 0) then
        ! Integration
        do j = 0,jx-1
           fqq(m,m/N) = fqq(m,m/N) + qq(j,m/N)*pm(j)*sin(y(j))*dy
        enddo
        
        if (m /= 0) then
           do j = 0,jx-1
              fqq(m,kx-m/N) = fqq(m,kx-m/N) + qq(j,kx-m/N)*pn(j)*sin(y(j))*dy
           enddo
        endif

        do j = 0,jx-1
           pm1(j) = pm(j)
           pm2(j) = 0.d0
        enddo

        do l = m+1,N*kx/2-1
           call legendre_l_up(l,m,y,pm1,pm2,jx,pm0,pn0)

           do j = 0,jx-1
              ! Integration
              fqq(l,m/N) = fqq(l,m/N) + qq(j,m/N)*pm0(j)*sin(y(j))*dy
           enddo

           ! negative m        
           if (m /= 0) then
              do j = 0,jx-1
                 fqq(l,kx-m/N) = fqq(l,kx-m/N) + qq(j,kx-m/N)*pn0(j)*sin(y(j))*dy
              enddo
           endif
                
           do j = 0,jx-1        
              pm2(j) = pm1(j)
              pm1(j) = pm0(j)
           enddo
        enddo 
     endif ! mod
  enddo

  do k = 0,kx-1
  do j = 0,N*kx/2-1
      fqq(j,k) = 0.5d0*fqq(j,k)
   enddo
   enddo

     
  return
end subroutine forward

subroutine backward(N,qq,y,z,jx,kx,fqq)
  implicit none

  integer :: j,k,m,l
  integer, intent(in) :: N,jx,kx
  double precision :: dy, dz
  double precision, dimension(0:jx-1) :: pm, pn, pn0, pm0, pm1, pm2 ! legendre function
  double precision, dimension(0:jx-1), intent(in) :: y  
  double precision, dimension(0:kx-1), intent(in) :: z
  double complex, dimension(0:N*jx-1,0:kx-1), intent(in)  :: qq
  double complex, dimension(0:jx-1,0:kx-1), intent(out) :: fqq

  dy = y(1) - y(0)
  dz = z(1) - z(0)
  
  do k = 0,kx-1
  do j = 0,jx-1
     fqq(j,k) = (0.d0,0.d0)
  enddo
  enddo
  
  do m = 0,N*kx/2-1
     call legendre_m_up(m,y,pm,jx,pn)

     if(mod(m,n) == 0) then
        ! Integration
        do j = 0,jx-1
           fqq(j,m/N) = fqq(j,m/N) + qq(m,m/N)*pm(j)
        enddo
        if (m /= 0) then
           do j = 0,jx-1
              fqq(j,kx-m/N) = fqq(j,kx-m/N) + qq(m,kx-m/N)*pn(j)
           enddo
        endif
     
        do j = 0,jx-1
           pm1(j) = pm(j)
           pm2(j) = 0.d0
        enddo
        
        do l = m+1,N*kx/2-1
           call legendre_l_up(l,m,y,pm1,pm2,jx,pm0,pn0)

           do j = 0,jx-1
              ! Integration
              fqq(j,m/N) = fqq(j,m/N) + qq(l,m/N)*pm0(j)
           enddo

           ! negative m
           if (m /= 0) then
              do j = 0,jx-1
                 ! Integration
                 fqq(j,kx-m/N) = fqq(j,kx-m/N) + qq(l,kx-m/N)*pn0(j)
              enddo
           endif
                
           do j = 0,jx-1
              pm2(j) = pm1(j)
              pm1(j) = pm0(j)
           enddo
        enddo
     endif
  enddo
     
  return
end subroutine backward
