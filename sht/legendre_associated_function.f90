!##################################################################
! initialize the Legendre function factors
subroutine init_factor(kx,N,epm,faca,facb)
   implicit none
   
   integer :: l,m,j
   integer, intent(in) :: kx,N
   double precision :: fm,fl
   double precision, dimension(1:N*kx/2), intent(out) :: epm
   double precision, dimension(0:N*kx/2,0:N*kx/2), intent(out) :: faca,facb

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

   return
end subroutine init_factor
!##################################################################
! initialize the Legendre function
subroutine init_geometry(y,jx,cosy,siny,sinydy)
   implicit none
   
   integer :: j
   integer, intent(in) :: jx
   double precision :: dy
   double precision, dimension(0:jx-1), intent(in) :: y
   double precision, dimension(0:jx-1), intent(out) :: cosy,siny,sinydy

   dy = y(1) - y(0)

   do j = 0,jx-1
      cosy(j) = cos(y(j))
      siny(j) = sin(y(j))
      sinydy(j) = sin(y(j))*dy
   enddo

   return
end subroutine init_geometry
!##################################################################
! calculate associated Legendre function P_m^m from P_{m-1}^{m-1}
subroutine m_up(m,siny,epm_m,pm,jx)
  implicit none

  integer :: j
  integer, intent(in) :: m,jx
  
  double precision, dimension(0:jx-1), intent(in) :: siny
  double precision, intent(in) :: epm_m
  double precision, dimension(0:jx-1), intent(inout) :: pm
  
  ! l = m case
  do j = 0,jx-1
      pm(j) = epm_m*siny(j)*pm(j)
  enddo
  do j = 0,jx-1
      if (abs(pm(j)) < 1.d-100) then
         pm(j) = 0.d0
      endif
   enddo

  return
end subroutine m_up

!##################################################################

subroutine l_up(m,cosy,faca_lm,facb_lm,pm1,pm2,jx,pm0) !,pn0)
  implicit none

  integer :: j
  integer, intent(in) :: m,jx
  
  double precision, dimension(0:jx-1), intent(in) :: cosy,pm1,pm2
  double precision, dimension(0:jx-1), intent(out) :: pm0
  double precision, intent(in) :: faca_lm,facb_lm
  do j = 0,jx-1
     pm0(j) = faca_lm*cosy(j)*pm1(j) - facb_lm*pm2(j)
  enddo

  return
end subroutine l_up

!##################################################################
! forword Legendre associated function transform
subroutine forward(N,qqg,yg,jxg,kx,fqq) bind(C)
  use omp_lib
  use iso_fortran_env
   implicit none

  integer :: j,k,m,l,m_n,jx
  double precision, dimension(1:N*kx/2) :: epm
  double precision, dimension(0:N*kx/2,0:N*kx/2) :: faca,facb
  integer, intent(in) :: N,jxg,kx

  ! Global variables
  double precision, dimension(0:jxg-1), intent(in) :: yg
  double complex, dimension(0:jxg-1,0:kx/2), intent(in) :: qqg
  double complex, dimension(0:N*kx/2,0:kx/2), intent(out) :: fqq
  double complex, allocatable :: fqq_local(:,:)

  ! OpenMP local variables
  integer :: jx0,jx1
  double precision, allocatable, dimension(:) :: siny,cosy,sinydy
  double precision, allocatable, dimension(:) :: pm, pm0, pm1, pm2 ! legendre function
  double precision, allocatable, dimension(:) :: y
  double complex, allocatable, dimension(:,:) :: qq

  integer :: OMP_N, OMP_ID
  
  fqq = (0.d0,0.d0)

  call init_factor(kx,N,epm,faca,facb)

  !$OMP parallel private(OMP_N,OMP_ID,jx,jx0,jx1,y,siny,cosy,sinydy,qq &
  !$OMP ,pm,pm0,pm1,pm2,j,k,l,m,m_n,fqq_local) shared(qqg,fqq)
  
  allocate(fqq_local(0:N*kx/2,0:kx/2))
  fqq_local = (0.d0,0.d0)

  OMP_N  = OMP_GET_NUM_THREADS()
  OMP_ID = OMP_GET_THREAD_NUM()

  jx = jxg/OMP_N
  jx0 = jx*OMP_ID
  if (OMP_ID == OMP_N-1) then
   jx = jxg-jx*(OMP_N-1)
  endif
  jx1 = jx0 + jx - 1

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

  call init_geometry(y,jx,cosy,siny,sinydy)
   
   m = 0
   do j = 0,jx-1
    pm(j) = sinydy(j)/sqrt(2.d0)
   enddo

   ! integration
   do j = 0,jx-1
      fqq_local(m,m) = fqq_local(m,m) + qq(j,m)*pm(j)
   enddo

   do j = 0,jx-1
    pm1(j) = pm(j)
    pm2(j) = 0.d0
   enddo
 
   do l = m+1,N*kx/2
      call l_up(m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,pm0)
 
      do j = 0,jx-1
       ! Integration
       fqq_local(l,m/N) = fqq_local(l,m/N) + qq(j,m/N)*pm0(j)
      enddo
         
    do j = 0,jx-1
       pm2(j) = pm1(j)
       pm1(j) = pm0(j)
    enddo
   enddo

   do m = 1,N*kx/2
      call m_up(m,siny,epm(m),pm,jx)
      
      if (mod(m,n) == 0) then
          m_N = m/N
         ! Integration

          do j = 0,jx-1             
            fqq_local(m,   m_N) = fqq_local(m,   m_N) + qq(j,   m_N)*pm(j)
        enddo
      
         do j = 0,jx-1
          pm1(j) = pm(j)
          pm2(j) = 0.d0
         enddo
 
         do l = m+1,N*kx/2
            call l_up(m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,pm0)
 
            do j = 0,jx-1
               ! Integration
               fqq_local(l,   m_N) = fqq_local(l,   m_N) + qq(j,   m_N)*pm0(j)
            enddo

            do j = 0,jx-1        
               pm2(j) = pm1(j)
               pm1(j) = pm0(j)
            enddo
         enddo 
      endif ! mod
   enddo

  !$omp critical
   do k = 0,kx/2
   do j = 0,N*kx/2
       fqq(j,k) = fqq(j,k) + fqq_local(j,k)
   enddo
   enddo   
  !$omp end critical
  !$OMP end parallel
 
  return
end subroutine forward

!##################################################################
! backward Legendre associated function transform
subroutine backward(N,qq,yg,jxg,kx,fqqg) bind(C)
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

  call init_factor(kx,N,epm,faca,facb)

  !$OMP parallel private(OMP_N,OMP_ID,jx,jx0,jx1,y,siny,cosy,sinydy &
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
 
   call init_geometry(y,jx,cosy,siny,sinydy)

   m = 0
   do j = 0,jx-1
      pm(j) = 1.d0/sqrt(2.d0)
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
   call l_up(m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,pm0)

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
     call m_up(m,siny,epm(m),pm,jx)

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
         call l_up(m,cosy,faca(l,m),facb(l,m),pm1,pm2,jx,pm0)

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

  do j = 0,jx-1
  do k = 0,kx/2
      fqqg(jx0+j,k) = fqq(j,k)
  enddo
  enddo

 !$OMP end parallel
     
  return
end subroutine backward
