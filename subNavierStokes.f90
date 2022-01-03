!    navier-stokes equations for constant properties
      subroutine ns_momentum(fluidCellCount,fluidIndexPtr,deltax, deltay, deltat, &
							 deltx2, delty2, u, v, p, ut, vt, alpha, re)
      INTEGER:: i,j, n
	  INTEGER ,intent(in)    :: fluidCellCount
	  INTEGER (KIND = 8),intent(in)  :: fluidIndexPtr(:, :)
	  REAL (KIND = 8), intent(in)        ::  deltax(:), deltay(:)
	  REAL (KIND = 8), intent(in)        ::  deltat, alpha, re 
	  REAL (KIND = 8), intent(out)        ::  deltx2, delty2
	  REAL (KIND = 8),intent(in):: u(:, :), v(:, :),  &
                                   p(:, :)
	  REAL (KIND = 8),intent(inout)::  ut(:, :), vt(:, :)
                                       
      REAL(KIND=8):: u1a,u22,u3,u4,u5,u6,u7,u8,u13,u14,&
       v1a,v22,v3,v4,v5,v6,v7,v8,v9,v10,d2udx2,d2udy2, &
       duudx,dvudy,dpdx,d2vdx2,d2vdy2,duvdx,dvvdy,dpdy
      
      DO 70 n = 1, fluidCellCount
         i = fluidIndexPtr(n, 1)
         j = fluidIndexPtr(n, 2) 
!     duu/dx         
        deltx2 = deltax(i)**2
        delty2 = deltay(j)**2
        u1a = u(i-1,j)  +  u(i,j)
        u22 = u(i-1,j)  -  u(i,j)
        u3  = u(i,j)    +  u(i+1,j)
        u4  = u(i,j)    -  u(i+1,j)
        

!     duv/dy
        u5 = u(i,j-1)   +  u(i,j)
        u6 = u(i,j-1)   -  u(i,j)
        u7 = u(i,j)     +  u(i,j+1)
        u8 = u(i,j)     -  u(i,j+1)

!     dwu/dz
       ! u9 = u(i,j-1)   + u(i,j)
        !u10= u(i,j-1)   - u(i,j)
        !u11= u(i,j)     + u(i,j+1)
        !u12= u(i,j)     - u(i,j+1)

!     duv/dx
        u13 = u(i-1,j)  + u(i-1,j+1)
        u14 = u7

!     duw / dx
        !u15 = u(i-1,j)  + u(i-1,j+1)
        !u16 = u11


!cccccccccccccccccccccc    diff -v    ccccccccccccccccccccccc
!     dvu / dx
        v1a = v(i,j-1)   + v(i+1,j-1)     
        v22 = v(i,j)     + v(i+1,j) 

!     duv / dx
        v3 = v(i-1,j)    + v(i,j)
        v4 = v(i-1,j)    - v(i,j)
        v5 = v(i,j)      + v(i+1,j)
        v6 = v(i,j)      - v(i+1,j)
!c
!c     dvv / dy
        v7 = v(i,j-1)    + v(i,j)
        v8 = v(i,j-1)    - v(i,j)
        v9 = v(i,j)      + v(i,j+1)
        v10= v(i,j)      - v(i,j+1)  
!c
!c     dwu / dz
       ! v11 = v(i,j-1)   + v(i,j)
       ! v12 = v(i,j-1)   - v(i,j) 
        !v13 = v(i,j)     + v(i,j+1)
       ! v14 = v(i,j)     - v(i,j+1)    
!c
!c     dvw / dy
       ! v15 = v(i,j-1)   + v(i,j-1+1)
       ! v16 = v13
!c
!c
!c
!c
!cccccccccccccccccc  diff - w  cccccccccccccccccccccc
!c     dwu / dz

!cccccccccccccccccccccccc u - deq cccccccccccccccccccccccccccccccc
!c
      d2udx2 = (u22 - u4 ) / deltx2
      d2udy2 = (u6  - u8 ) / delty2
!c
      duudx = 0.25*( u1a*u1a + alpha*dabs(u1a)*u22 - &
     &              u3*u3 -   alpha*dabs(u3)*u4)/deltax(i)
!c
      dvudy = 0.25*( v1a*u5 + alpha*dabs(v1a)*u6 - &
     &               v22*u7 - alpha*dabs(v22)*u8)/deltay(j)
!c
      dpdx=( p(i+1,j) - p(i,j) )/deltax(i)  
      
       
      ut(i,j) = deltat*(+duudx+dvudy+1./re*(d2udx2+d2udy2))+u(i,j) &
     &          -deltat*dpdx
      d2vdx2 = (v4 - v6)/deltx2
      d2vdy2 = (v8 - v10)/delty2
      !d2vdz2 = (v12 - v14)/deltz2
      duvdx = 0.25*(u13*v3 + alpha*dabs(u13)*v4- &
     &              u14*v5 - alpha*dabs(u14)*v6)/deltax(i)
!c
      dvvdy = 0.25*(v7*v7 + alpha*dabs(v7)*v8 - &
     &              v9*v9 - alpha*dabs(v9)*v10)/deltay(j)
!c
      dpdy= ( p(i,j+1) - p(i,j))/deltay(j)
      vt(i,j) = deltat*(+duvdx+dvvdy+1./re*(d2vdx2+d2vdy2 &
     &              ))+v(i,j)-deltat*dpdy
      !print*, u(i,j), ut(i,j), v(i,j), vt(i,j)  
      !print*,i, j,ut(i,j), u(i,j-1)

70   continue
!c    
      !stop
      !write(6,*) 'nsMomentum complete'
      
      end subroutine ns_momentum

!cssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
