!##################################################################
! initialize the Legendre funnction
subroutine legendre_init(y,jx,kx,N,siny,cosy,sinydy,epm,faca,facb)
   implicit none
   
   integer :: l,m,j
   integer, intent(in) :: jx,kx,N
   double precision :: fm,fl,dy
   double precision, dimension(0:jx-1), intent(in) :: y
   double precision, dimension(0:jx-1), intent(out) :: siny,cosy,sinydy
   double precision, dimension(0:N*kx/2-1), intent(out) :: epm
   double precision, dimension(1:N*kx/2-1,0:N*kx/2-1), intent(out) :: faca,facb

   dy = y(1) - y(0)

   faca = 0.d0
   facb = 0.d0

   do m = 0,N*kx/2-1
      fm = dble(m)
      epm(m) = -sqrt((2.d0*fm+1.d0)/(2.d0*fm))
      do l = m+1,N*kx/2-1
         fl = dble(l)
         faca(l,m) = sqrt((2.d0*fl-1.d0)*(2.d0*fl+1.d0)/(fl+fm)/(fl-fm))
         facb(l,m) = sqrt((2.d0*fl+1.d0)/(2.d0*fl-3.d0)*(fl+fm-1.d0)*(fl-fm-1.d0)/(fl+fm)/(fl-fm))
      enddo
   enddo

   do j = 0,jx-1
      cosy(j) = cos(y(j))
      siny(j) = sin(y(j))
      sinydy(j) = sin(y(j))*dy
   enddo

   return
end subroutine legendre_init
!##################################################################
! calculate associated Legendre function P_m^m from P_{m-1}^{m-1}
subroutine legendre_m_up(m,siny,epm_m,pm,jx,kx,N,pn)
  implicit none

  integer :: j
  integer, intent(in) :: m,jx,kx,N
  
  double precision, dimension(0:jx-1), intent(in) :: siny
  double precision, intent(in) :: epm_m
  double precision, dimension(0:jx-1), intent(inout) :: pm
  double precision, dimension(0:jx-1), intent(out) :: pn
  
  ! l = m case
  do j = 0,jx-1
      pm(j) = epm_m*siny(j)*pm(j)
  enddo
  do j = 0,jx-1
      if (abs(pm(j)) < 1.d-300) then
         pm(j) = 0.d0
      endif
   enddo

  ! negative m
  do j = 0,jx-1
     pn(j) = (-1.d0)**m*pm(j)
  enddo

  return
end subroutine legendre_m_up

!##################################################################

subroutine legendre_l_up(l,m,cosy,faca_lm,facb_lm,pm1,pm2,jx,kx,N,pm0,pn0)
  implicit none

  integer :: j
  integer, intent(in) :: l,m,jx,kx,N
  
  double precision, dimension(0:jx-1), intent(in) :: cosy,pm1,pm2
  double precision, dimension(0:jx-1), intent(out) :: pm0,pn0
  double precision, intent(in) :: faca_lm,facb_lm
    

  !faca = sqrt((2.d0*fl-1.d0)*(2.d0*fl+1.d0)/(fl+fm)/(fl-fm))
  !facb = sqrt((2.d0*fl+1.d0)/(2.d0*fl-3.d0)*(fl+fm-1.d0)*(fl-fm-1.d0)/(fl+fm)/(fl-fm))

  do j = 0,jx-1
     pm0(j) = faca_lm*cosy(j)*pm1(j) - facb_lm*pm2(j)
  enddo

  do j = 0,jx-1
      pn0(j) = (-1.0)**m*pm0(j)        
  enddo
  
  return
end subroutine legendre_l_up

!##################################################################
subroutine forward(N,qq,y,jx,kx,fqq)
  implicit none

  integer :: j,k,m,l,m_N
  integer, intent(in) :: N,jx,kx
  double precision, dimension(0:jx-1) :: siny,cosy,sinydy
  double precision, dimension(0:N*kx/2-1) :: epm
  double precision, dimension(1:N*kx/2-1,0:N*kx/2-1) :: faca,facb
  double precision, dimension(0:jx-1) :: pm, pn, pn0, pm0, pm1, pm2 ! legendre function
  double precision, dimension(0:jx-1), intent(in) :: y
  double complex, dimension(0:jx-1,0:kx-1), intent(in)  :: qq
  double complex, dimension(0:N*kx/2-1,0:kx-1), intent(out) :: fqq
    
  call legendre_init(y,jx,kx,N,siny,cosy,sinydy,epm,faca,facb)

  do k = 0,kx-1
  do j = 0,N*kx/2-1
     fqq(j,k) = (0.d0,0.d0)
  enddo
  enddo
  
  m = 0
  do j = 0,jx-1
   pm(j) = 1.d0
  enddo

  do j = 0,jx-1
   pm1(j) = pm(j)
   pm2(j) = 0.d0
  enddo

  do l = m+1,N*kx/2-1
   call legendre_l_up(l,m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,kx,N,pm0,pn0)

   do j = 0,jx-1
      ! Integration
      fqq(l,m/N) = fqq(l,m/N) + qq(j,m/N)*pm0(j)*sinydy(j)
   enddo
        
   do j = 0,jx-1        
      pm2(j) = pm1(j)
      pm1(j) = pm0(j)
   enddo
  enddo   

  do m = 1,N*kx/2-1
   call legendre_m_up(m,siny,epm(m),pm,jx,kx,N,pn)
     
     if (mod(m,n) == 0) then
         m_N = m/N
        ! Integration
        do j = 0,jx-1
           fqq(m,m_N) = fqq(m,m_N) + qq(j,m_N)*pm(j)*sinydy(j)
        enddo
        
        do j = 0,jx-1
            fqq(m,kx-m_N) = fqq(m,kx-m_N) + qq(j,kx-m_N)*pn(j)*sinydy(j)
        enddo

        do j = 0,jx-1
         pm1(j) = pm(j)
         pm2(j) = 0.d0
        enddo

        do l = m+1,N*kx/2-1
           call legendre_l_up(l,m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,kx,N,pm0,pn0)

            do j = 0,jx-1
              ! Integration
              fqq(l,m_N) = fqq(l,m_N) + qq(j,m_N)*pm0(j)*sinydy(j)
           enddo

           do j = 0,jx-1
               fqq(l,kx-m_N) = fqq(l,kx-m_N) + qq(j,kx-m_N)*pn0(j)*sinydy(j)
           enddo
                
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

!##################################################################
subroutine backward(N,qq,y,z,jx,kx,fqq)
  implicit none

  integer :: j,k,m,l
  integer, intent(in) :: N,jx,kx
  double precision, dimension(0:jx-1) :: siny,cosy,sinydy
  double precision, dimension(0:N*kx/2-1) :: epm
  double precision, dimension(1:N*kx/2-1,0:N*kx/2-1) :: faca,facb
  double precision, dimension(0:jx-1) :: pm, pn, pn0, pm0, pm1, pm2 ! legendre function
  double precision, dimension(0:jx-1), intent(in) :: y  
  double precision, dimension(0:kx-1), intent(in) :: z
  double complex, dimension(0:N*kx/2-1,0:kx-1), intent(in)  :: qq
  double complex, dimension(0:jx-1,0:kx-1), intent(out) :: fqq
  
  call legendre_init(y,jx,kx,N,siny,cosy,sinydy,epm,faca,facb)

  do k = 0,kx-1
  do j = 0,jx-1
     fqq(j,k) = (0.d0,0.d0)
  enddo
  enddo

  m = 0
  do j = 0,jx-1
   pm(j) = 1.d0
  enddo

  do j = 0,jx-1
   pm1(j) = pm(j)
   pm2(j) = 0.d0
  enddo

  do l = m+1,N*kx/2-1
   call legendre_l_up(l,m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,kx,N,pm0,pn0)

   do j = 0,jx-1
      ! Integration
      fqq(j,m/N) = fqq(j,m/N) + qq(l,m/N)*pm0(j)
   enddo
        
   do j = 0,jx-1        
      pm2(j) = pm1(j)
      pm1(j) = pm0(j)
   enddo
  enddo 
  
  do m = 1,N*kx/2-1
     call legendre_m_up(m,siny,epm(m),pm,jx,kx,N,pn)   

     if(mod(m,n) == 0) then
        ! Integration
        do j = 0,jx-1
           fqq(j,m/N) = fqq(j,m/N) + qq(m,m/N)*pm(j)
        enddo
        do j = 0,jx-1
           fqq(j,kx-m/N) = fqq(j,kx-m/N) + qq(m,kx-m/N)*pn(j)
        enddo
     
        do j = 0,jx-1
           pm1(j) = pm(j)
           pm2(j) = 0.d0
        enddo
        
        do l = m+1,N*kx/2-1
         call legendre_l_up(l,m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,kx,N,pm0,pn0)

         do j = 0,jx-1
            ! Integration
            fqq(j,m/N) = fqq(j,m/N) + qq(l,m/N)*pm0(j)
         enddo

         ! negative m
         do j = 0,jx-1
            ! Integration
            fqq(j,kx-m/N) = fqq(j,kx-m/N) + qq(l,kx-m/N)*pn0(j)
         enddo
                
         do j = 0,jx-1
            pm2(j) = pm1(j)
            pm1(j) = pm0(j)
         enddo
        enddo
     endif
  enddo
     
  return
end subroutine backward
