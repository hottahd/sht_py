!##################################################################
! initialize the Legendre function
subroutine legendre_init(y,jx,kx,N,siny,cosy,sinydy,epm,faca,facb)
   implicit none
   
   integer :: l,m,j
   integer, intent(in) :: jx,kx,N
   double precision :: fm,fl,dy
   double precision, dimension(0:jx-1), intent(in) :: y
   double precision, dimension(0:jx-1), intent(out) :: siny,cosy,sinydy
   double precision, dimension(1:N*kx/2), intent(out) :: epm
   double precision, dimension(0:N*kx/2,0:N*kx/2), intent(out) :: faca,facb

   dy = y(1) - y(0)
   faca = 0.d0
   facb = 0.d0
   do m = 1,N*kx/2
      fm = dble(m)
      epm(m) = -sqrt((2.d0*fm+1.d0)/(2.d0*fm))
   enddo
   
   do m = 0,N*kx/2
      fm = dble(m)
      do l = m+1,N*kx/2
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
subroutine legendre_m_up(m,siny,epm_m,pm,jx)!,pn)
  implicit none

  integer :: j
  integer, intent(in) :: m,jx
  
  double precision, dimension(0:jx-1), intent(in) :: siny
  double precision, intent(in) :: epm_m
  double precision, dimension(0:jx-1), intent(inout) :: pm
  !double precision, dimension(0:jx-1), intent(out) :: pn
  
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
  !do j = 0,jx-1
  !   pn(j) = (-1.d0)**m*pm(j)
  !enddo

  return
end subroutine legendre_m_up

!##################################################################

subroutine legendre_l_up(m,cosy,faca_lm,facb_lm,pm1,pm2,jx,pm0)!,pn0)
  implicit none

  integer :: j
  integer, intent(in) :: m,jx
  
  double precision, dimension(0:jx-1), intent(in) :: cosy,pm1,pm2
  double precision, dimension(0:jx-1), intent(out) :: pm0!,pn0
  double precision, intent(in) :: faca_lm,facb_lm
  do j = 0,jx-1
     pm0(j) = faca_lm*cosy(j)*pm1(j) - facb_lm*pm2(j)
  enddo

  !do j = 0,jx-1
  !    pn0(j) = (-1.0)**m*pm0(j)
  !enddo
  
  return
end subroutine legendre_l_up

!##################################################################
! forword Spherical Harmonic Expansion
subroutine forward(N,qqg,yg,jxg,kx,fqqg)
   use omp_lib
   implicit none

  integer :: j,k,m,l,m_n,jx
  double precision, dimension(1:N*kx/2) :: epm
  double precision, dimension(0:N*kx/2,0:N*kx/2) :: faca,facb
  integer, intent(in) :: N,jxg,kx

  ! Global variables
  double precision, dimension(0:jxg-1), intent(in) :: yg
  double complex, dimension(0:jxg-1,0:kx/2), intent(in) :: qqg
  double complex, dimension(0:N*kx/2,0:kx/2), intent(out) :: fqqg

  ! OpenMP local variables
  integer :: jx0,jx1
  double precision, allocatable, dimension(:) :: siny,cosy,sinydy
  double precision, allocatable, dimension(:) :: pm, pm0, pm1, pm2 ! legendre function
  double precision, allocatable, dimension(:) :: y
  double complex, allocatable, dimension(:,:) :: qq
  double complex, allocatable, dimension(:,:) :: fqq
  double complex, allocatable, dimension(:,:,:) :: fqq_OMP

  integer :: OMP_N, OMP_ID
  
  do k = 0,kx/2
  do j = 0,N*kx/2
       fqqg(j,k) = (0.d0,0.d0)
  enddo
  enddo  

  !$OMP parallel private(OMP_ID,jx,jx0,jx1,y,siny,cosy,sinydy,qq &
  !$OMP ,pm,pm0,pm1,pm2,j,k,l,m,m_n,fqq)
  OMP_N  = OMP_GET_NUM_THREADS()
  OMP_ID = OMP_GET_THREAD_NUM()

  jx = jxg/OMP_N
  jx0 = jx*OMP_ID
  if (OMP_ID == OMP_N-1) then
   jx = jxg-jx*(OMP_N-1)
  endif
  jx1 = jx0 + jx - 1

  allocate(fqq(0:N*kx/2,0:kx/2))
  
  if (OMP_ID == 0) then
   allocate(fqq_OMP(0:N*kx/2,0:kx/2,0:OMP_N-1))
  endif

  allocate(y(0:jx-1))
  allocate(siny(0:jx-1))
  allocate(cosy(0:jx-1))
  allocate(sinydy(0:jx-1))
  allocate(qq(0:jx-1,0:kx/2))

  allocate(pm(0:jx-1))
  allocate(pm0(0:jx-1))
  allocate(pm1(0:jx-1))
  allocate(pm2(0:jx-1))

  y = yg(jx0:jx1)
  qq = qqg(jx0:jx1,:)

  call legendre_init(y,jx,kx,N,siny,cosy,sinydy,epm,faca,facb)

  do k = 0,kx/2
  do j = 0,N*kx/2
      fqq(j,k) = (0.d0,0.d0)
   enddo
   enddo
   
   m = 0
   do j = 0,jx-1
    pm(j) = 1.d0
   enddo

   ! integration
   do j = 0,jx-1
      fqq(m,m) = fqq(m,m) + qq(j,m)*pm(j)*sinydy(j)
   enddo

   do j = 0,jx-1
    pm1(j) = pm(j)
    pm2(j) = 0.d0
   enddo
 
   do l = m+1,N*kx/2
      call legendre_l_up(m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,pm0)
 
      do j = 0,jx-1
       ! Integration
       fqq(l,m/N) = fqq(l,m/N) + qq(j,m/N)*pm0(j)*sinydy(j)
       enddo
         
    do j = 0,jx-1
       pm2(j) = pm1(j)
       pm1(j) = pm0(j)
    enddo
   enddo

   do m = 1,N*kx/2
      call legendre_m_up(m,siny,epm(m),pm,jx)
      
      if (mod(m,n) == 0) then
          m_N = m/N
         ! Integration

          !$OMP SIMD
          do j = 0,jx-1             
            fqq(m,   m_N) = fqq(m,   m_N) + qq(j,   m_N)*pm(j)*sinydy(j)
        enddo
      
         do j = 0,jx-1
          pm1(j) = pm(j)
          pm2(j) = 0.d0
         enddo
 
         do l = m+1,N*kx/2
            call legendre_l_up(m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,pm0)
 
            !$OMP SIMD
            do j = 0,jx-1
               ! Integration
               fqq(l,   m_N) = fqq(l,   m_N) + qq(j,   m_N)*pm0(j)*sinydy(j)
            enddo

            do j = 0,jx-1        
               pm2(j) = pm1(j)
               pm1(j) = pm0(j)
            enddo
         enddo 
      endif ! mod
   enddo
 
   do k = 0,kx/2
   do j = 0,N*kx/2
      fqq_OMP(j,k,OMP_ID) = fqq(j,k)
    enddo
    enddo
  !$OMP end parallel
  
   do OMP_ID = 0,OMP_N-1
   do k = 0,kx/2
   do j = 0,N*kx/2
      fqqg(j,k) = fqqg(j,k) + 0.5d0*fqq_OMP(j,k,OMP_ID)
   enddo         
   enddo
   enddo
  return
end subroutine forward

!##################################################################
! backword Spherical Harmonic Expansion
subroutine backward(N,qq,yg,jxg,kx,fqqg)
   use omp_lib
   implicit none

   integer :: j,k,m,l,m_N,jx
   double precision, dimension(1:N*kx/2) :: epm
   double precision, dimension(0:N*kx/2,0:N*kx/2) :: faca,facb
   integer, intent(in) :: N,jxg,kx

   ! Global variables
   double precision, dimension(0:jxg-1), intent(in) :: yg
   double complex, dimension(0:N*kx/2,0:kx/2), intent(in) :: qq
   double complex, dimension(0:jxg-1,0:kx/2), intent(out) :: fqqg

   ! OpenMP local variables
   integer :: jx0,jx1
   double precision, allocatable, dimension(:) :: siny,cosy,sinydy
   double precision, allocatable, dimension(:) :: pm, pm0, pm1, pm2 ! legendre function
   double precision, allocatable, dimension(:) :: y  
   double complex, allocatable, dimension(:,:) :: fqq

   integer :: OMP_N, OMP_ID

   do k = 0,kx/2
   do j = 0,jxg-1
      fqqg(j,k) = (0.d0,0.d0)
   enddo
   enddo

   !$OMP parallel private(OMP_ID,jx,jx0,jx1,y,siny,cosy,sinydy &
   !$OMP ,pm,pm0,pm1,pm2,j,k,l,m,m_n,fqq)
   OMP_N  = OMP_GET_NUM_THREADS()
   OMP_ID = OMP_GET_THREAD_NUM()
      
   jx = jxg/OMP_N
   jx0 = jx*OMP_ID
   if (OMP_ID == OMP_N-1) then
    jx = jxg-jx*(OMP_N-1)
   endif  
   jx1 = jx0 + jx - 1
   
   allocate(fqq(0:jx-1,0:kx/2))
   
   allocate(y(0:jx-1))
   allocate(siny(0:jx-1))
   allocate(cosy(0:jx-1))
   allocate(sinydy(0:jx-1))
 
   allocate(pm(0:jx-1))
   allocate(pm0(0:jx-1))
   allocate(pm1(0:jx-1))
   allocate(pm2(0:jx-1))
 
   y = yg(jx0:jx1)
 
   call legendre_init(y,jx,kx,N,siny,cosy,sinydy,epm,faca,facb)

   m = 0
   do j = 0,jx-1
      pm(j) = 1.d0
   enddo

   do k = 0,kx/2
   do j = 0,jx-1
      fqq(j,k) = (0.d0,0.d0)
   enddo
   enddo

   ! Integration
   do j = 0,jx-1
      fqq(j,m/N) = fqq(j,m/N) + qq(m,m/N)*pm(j)
   enddo

   do j = 0,jx-1
   pm1(j) = pm(j)
   pm2(j) = 0.d0
  enddo

  do l = m+1,N*kx/2
   call legendre_l_up(m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,pm0)

   do j = 0,jx-1
      ! Integration
      fqq(j,m/N) = fqq(j,m/N) + qq(l,m/N)*pm0(j)
   enddo
        
   do j = 0,jx-1        
      pm2(j) = pm1(j)
      pm1(j) = pm0(j)
   enddo
  enddo
  
  do m = 1,N*kx/2
     call legendre_m_up(m,siny,epm(m),pm,jx)

     if(mod(m,n) == 0) then
      m_n = m/N
      ! Integration
      do j = 0,jx-1
         fqq(j,   m_N) = fqq(j,   m_N) + qq(m,   m_N)*pm(j)
      enddo

      do j = 0,jx-1
           pm1(j) = pm(j)
           pm2(j) = 0.d0
      enddo
        
        do l = m+1,N*kx/2
         call legendre_l_up(m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,pm0)

         do j = 0,jx-1
            ! Integration
            fqq(j,   m_N) = fqq(j,   m_N) + qq(l,   m_N)*pm0(j)
         enddo
                
         do j = 0,jx-1
            pm2(j) = pm1(j)
            pm1(j) = pm0(j)
         enddo
        enddo
     endif
  enddo

  do k = 0,kx/2
  do j = 0,jx-1
      fqqg(jx0+j,k) = fqq(j,k)
  enddo
  enddo

  !$OMP end parallel
     
  return
end subroutine backward
