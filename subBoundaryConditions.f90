

      
!c********************************************************************            
SUBROUTINE velocity_bc (nx, ny, cell, ut, vt)
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      REAL (KIND = 8) :: t1, t2, yl, yushift, uc
      INTEGER (KIND = 8):: i, j, n, vbcOption
	  INTEGER ( KIND = 8), intent(in)   :: nx, ny, cell(:)
	  REAL (KIND = 8), intent(inout) :: ut(:,:), vt(:,:)
      yushift = 0._rk
      yl = 0._rk
      uc = 1._rk
      DO 30 j = 2, ny+1 
        !yl = (j-2)*deltay + deltay/2._rk 
        !ut(1, j)     =  0.
        !IF ((j-1).GT.ny/3.AND.(j-1).LE.2*ny/3) ut(1, j) = -6*(yl-yushift)*(yl-yushift-1._rk)    !sudden expansion - parabolic inlet     
        !ut(1, j) = -6/ly**2*(yl-yushift)*(yu-ylshift-ly)                        !parabolic inlet
        !t1 = (1.0_rk - yl)
        !t2 = (yl - 0.5_rk)
        !ut(1, j) = dmax1(24._rk*t1*t2 ,0._rk)                                  !back-step
        !print*, j-1, ut(1,j)
        !ut(1, j) = ut(2, j)      
		ut(1, j) = ut(2, j)                                          
        !t1 = (1.0_rk - (deltay*(j-1)-deltay/2._rk))
        !t2 = (deltay*(j-1)-deltay/2._rk - 0.5_rk)
        !ut(1, j) = dmax1(24._rk*t1*t2 ,0._rk)
        !ut(1, j)     =  1._rk                                                 !uniform inlet
        vt(1, j)     =  vt(2, j)  
        ut(nx+1, j)  =  ut(nx, j)     
        !ut(nx+2, j)  =  ut(nx+1, j)                                            !Neumann - low Re convective flow
        vt(nx+2, j)  =  vt(nx+1,j)  
        !ut(nx+1, j)  =  u(nx+1,j)-(deltat/deltax(i))*uc*(u(nx+1,j)-u(nx,j))       !Orlanski - vortex shedding Re Convective flow
		  !ut(nx+2, j)  =  u(nx+2,j)-(deltat/deltax)*uc*(u(nx+2,j)-u(nx+1,j)) 
        !vt(nx+1, j)  =  v(nx+1,j)-(deltat/deltax(i))*uc*(v(nx+1,j)-v(nx,j))
        !vt(nx+2, j)  = -vt(nx+1, j) + v(nx+2,j) + v(nx+1,j) - (2._rk*deltat/deltax)*uc*(v(nx+2,j)-v(nx+1,j))
       ! ut(nx+2, j)  =  2._rk*ut(nx+1, j ) - ut(nx, j)                                       
 30   CONTINUE   
      !GOTO 2                                                                 !uncomment for internal flows
      DO 10 i = 2, nx+1
		  !ut(i, 1)        = -ut(i, 2)   
         ut(i, 1)        =  -ut(i, 2)     
         vt(i, 1)        =  0.!vt(i, 2)                        
         !ut(i, ny+2)     = 2._rk-ut(i, ny+1)                                      !wall no slip - closed channel
         ut(i, ny+2)     =  -ut(i, ny+1)                                       !symmetric bc - open channel
         vt(i, ny+1)     =  0.!vt(i, ny)  
         !vt(i, ny+2)     =  2._rk*vt(i, ny+1) - vt(i, ny)  
         !print*,  ut(i, ny+2), i 
 10   CONTINUE
      goto 100
      DO 40 j = 2, ny+1
      DO 40 i = 2, nx+1
         n   = i-1  + (nx)*(j-2)
         IF (cell(n).ne.0) THEN
           ! ut(i, j)  = 0
            !ut(i-1, j)  = 0
            vt(i, j)  = 0
            vt(i, j-1)  = 0
         ENDIF 
 40   CONTINUE
 100  CONTINUE          
      END SUBROUTINE velocity_bc
      
      SUBROUTINE solid_cell_bc(solidCellCount, solidIndexPtr, ut, vt, p)
         
         INTEGER, PARAMETER :: rk = selected_real_kind(8)
         INTEGER (KIND = 8):: i, j, n, vbcOption  
		 INTEGER, intent(in):: solidCellCount
		 INTEGER (KIND = 8), intent(in):: solidIndexPtr(:,:)
		 REAL (KIND = 8), intent(inout) :: ut(:,:), vt(:,:), p(:,:)
         DO n = 1, solidCellCount
            i = solidIndexPtr(n,1)
            j = solidIndexPtr(n,2)
            ut(i,j) = 0._rk
            vt(i,j) = 0._rk
            p(i,j) = 0._rk
         END DO
      END SUBROUTINE solid_cell_bc


