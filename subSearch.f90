     SUBROUTINE shift_surface_nodes(ibNodes, xShift, yShift, xnode1, ynode1, &
								  xnode, ynode, xt, xdot, xddot, yt, ydot, yddot)
             
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  i
		REAL (KIND = 8), intent(inout) :: xnode(:), ynode(:)
		REAL (KIND = 8), intent(in)    :: xShift, yShift
		INTEGER (KIND=8), intent(in)  ::ibNodes
		REAL (KIND = 8), intent(out)    ::xt, xdot, xddot, yt, ydot, yddot
		REAL (KIND = 8), intent(out) :: xnode1(ibNodes), ynode1(ibNodes) 
        
        xnode(:) = xnode(:) + xShift
        ynode(:) = ynode(:) + yShift   
        xnode1 = xnode
        ynode1 = ynode
        xt    = 0.!a0*(1. - cos(ang_theta*totime))
        xdot  = 0.!a0*ang_theta*sin(ang_theta*totime)
        xddot = 0.!a0*ang_theta**2*cos(ang_theta*totime)
        
        yt    = 0._rk
        ydot  = 0._rk
        yddot = 0._rk
                Print*, xt, xdot, xddot
      END SUBROUTINE shift_surface_nodes
!***********************************************************      
     SUBROUTINE compute_surface_variables(pi, freq, xt_temp, a0, totime, tempTime,&
										  xdot_temp, xddot_temp, yt_temp, ydot_temp, yddot_temp,&
										  xnode1, ynode1,xnode, ynode)

        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  i
        REAL (KIND=8)      :: ang_theta
		REAL (KIND =8), intent(in)    ::pi, freq, a0, totime, temptime
		REAL (KIND =8), intent(out)    ::xt_temp, yt_temp, xdot_temp, xddot_temp, ydot_temp, yddot_temp
		REAL (KIND = 8),intent(inout) :: xnode(:), ynode(:),  xnode1(:), ynode1(:)
        
        ang_theta =  2._rk*pi*freq     
                       
        xt_temp    = a0*(1. - cos(ang_theta*(totime-tempTime)))
        xdot_temp  = a0*ang_theta*sin(ang_theta*(totime-tempTime))
        xddot_temp = a0*ang_theta**2*cos(ang_theta*(totime-tempTime))
        
        yt_temp    = 0._rk
        ydot_temp  = 0._rk
        yddot_temp = 0._rk
        
        xnode(:) = xnode1(:) + xt_temp
        ynode(:) = ynode1(:) + yt_temp   

        Print*, xt_temp, xdot_temp, xddot_temp, totime_temp
		!Print*,'---', xnode(1)
        !Do i = 1, ibnodes
           Print*, 'var computed'
        !enddo
       ! pause
      END SUBROUTINE compute_surface_variables
! !***********************************************************   
      SUBROUTINE compute_surface_norm(ibElems,xnode, ynode, ibElP1, ibElP2, xcent, ycent, cosAlpha, cosBeta, area)
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n!, c1, c2, c3, c4
        REAL (KIND=8)      :: p1x, p1y, p2x, p2y, lenEL, inor
		INTEGER (KIND=8),intent(in)    :: ibElP1(:), ibElP2(:), ibElems
		REAL (KIND = 8), intent(in) :: xnode(:), ynode(:)
		REAL (KIND = 8), intent(out) :: xcent(ibElems), ycent(ibElems), cosAlpha(ibElems), cosBeta(ibElems), area(ibElems)
                             
        
        inor = -1._rk
        !compute centroid and direction cosines           
        DO n = 1, ibElems
           p1x = xnode(ibElP1(n))                       !x coordinate element node 1
           p1y = ynode(ibElP1(n))                       !y coordinate element node 1
           p2x = xnode(ibElP2(n))                       !x coordinate element node 2
           p2y = ynode(ibElP2(n))                       !y coordinate element node 2
           lenEL = dsqrt((p2y-p1y)**2 + (p2x-p1x)**2)   !length of element
           area(n) = lenEl
           xcent(n) =  (p2x+p1x)/2._rk                  !centroid x coordinate element
           ycent(n) =  (p2y+p1y)/2._rk                  !centroid y coordinate element          
           cosAlpha(n) = -(p2y-p1y)/lenEl*inor               !direction cosine unit normal along x
           cosBeta(n)  = (p2x-p1x)/lenEl*inor               !direction cosine unit normal along y
           !WRITE(*,*) xcent(n), ycent(n), cosAlpha(n)
           WRITE(196,*) xcent(n), ycent(n), cosBeta(n)
        ENDDO
        !print*, 'ibel done'     
       !stop
     END SUBROUTINE compute_surface_norm
! !******************************************************************

     SUBROUTINE tagging(nx, ny, x1, y1, ibElems, xcent, ycent, cosAlpha, cosBeta, cell, xp, yp, &
						ibCellCount, fluidCellCount, solidCellCount, nodeId )
               
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n, m, i, j, n1, n2, n3, n4, iNew, jNew, &
                               nel2n, sumId, ar_intcptd
        INTEGER            :: iPt, iPt1
        REAL (KIND=8)      :: n1x, n1y, n2x, n2y, n3x, n3y, dis , minDis, n1dotn,  &
                              cent_x, cent_y, sStart, sStop, ic1x, ic2x, ic3x, ic1y, ic2y, ic3y 
		INTEGER (KIND=8),intent(in)    :: ibElems, nx, ny
		INTEGER (KIND=8),intent(inout)    :: cell(:)
		REAL (KIND = 8), intent(in) :: xcent(:), ycent(:), cosAlpha(:), cosBeta(:), x1(:), y1(:) , xp(:), yp(:) 
		INTEGER ,intent(out) :: ibCellCount, fluidCellCount, solidCellCount
		INTEGER (KIND = 8), intent(out)  :: nodeId(nx+3,ny+3)

    
 
        !CALL cpu_time(sStart) 
        ibCellCount = 0  
        fluidCellCount = 0
        solidCellCount = 0
	     nodeId = 0    
        DO 10 j = 1, ny+3
        DO 10 i = 1, nx+3
           m = 0
           minDis = 1e14
           n1x = x1(i)
           n1y = y1(j)      
           n1dotn = 0
	        DO m = 1, ibElems
              cent_x = xcent(m)
              cent_y = ycent(m)
              dis  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2 )     
		        IF (dis.LT.minDis) THEN
                 minDis   = dis                 
                 nel2n    = m
              ENDIF
           ENDDO
           !print*, i, j, minDis
	        n1dotn = (n1x - xcent(nel2n))*cosAlpha(nel2n) + (n1y - ycent(nel2n))*cosBeta(nel2n)
           IF (n1dotn.LE.-1e-14) THEN
              nodeId(i,j) = 1  
           !ELSEIF (abs(n1dotn).LT.1e-16.AND.ABS(n1dotn).GT.0) THEN
              !nodeId((j-1)*(nx+3) + i) = 1  
           !ELSEIF (abs(n1dotn).EQ.0) THEN
              !nodeId((j-1)*(nx+3) + i) = 1  
           ENDIF
 10      CONTINUE
	      cell = 0
         DO 20 j = 2, ny+1
         DO 20 i = 2, nx+1
            n   = i-1  + (nx)*(j-2)
            sumId = 0
            sumId = nodeId(i,j) + nodeId(i+1,j) + nodeId(i,j+1) + nodeId(i+1,j+1)

            IF (sumId.eq.4) THEN
		         cell(n) = 1
		         solidCellCount = solidCellCount + 1
            ELSEIF (sumId.eq.0) THEN
               cell(n) = 0  
               fluidCellCount  = fluidCellCount + 1 
            ELSEIF (sumId.NE.0.and.sumId.NE.4) THEN 
               cell(n) = 2               
               ibCellCount = ibCellCount + 1 
            ENDIF 
 20      CONTINUE  
    GOTO 100
         print*, 'search done', ibCellCount, fluidCellCount, solidCellCount
         ibCellCount = 0
         DO 30 j = 2, ny+1
         DO 30 i = 2, nx+1
            n   = i-1  + (nx)*(j-2)  
            IF (Cell(n).eq.2) THEN
              IF (cell(n+1).eq.0.OR.cell(n-1).eq.0.OR.cell(n+nx).eq.0.OR.cell(n-nx).eq.0) THEN      
                 ibCellCount = ibCellCount + 1   
              ELSE
                 cell(n) = 1  
                 solidCellCount = solidCellCount + 1 
              ENDIF
            ENDIF
 30      CONTINUE
 100 CONTINUE
 
        open(82,file='inner.dat',status='unknown')
        write(82,*)'variables = "x", "y","z"'!, "dnorm","ibmloc"'
       !do i=1,iim
       do j=2,ny+1
       do i=2,nx+1
       n   = i-1  + (nx)*(j-2) 
       if(cell(n).eq.2)then
       write(82,*)xp(i),yp(j), 0
        !+ ,ibmloc(i,j,k)
       endif
       end do
       end do
      ! end do
      !stop
        print*, 'search done', ibCellCount, fluidCellCount, solidCellCount
        !stop             
     END SUBROUTINE tagging

! !******************************************************************

     SUBROUTINE selective_retagging(ibCellCount, interceptedIndexPtr, nodeId, x1, y1, &
									ibElems, xcent, ycent, cosAlpha, cosBeta, solidCellCount,&
									fluidCellCount, nx, ny, cell, xp, yp, ibCellCount_temp,&
									fluidCellCount_temp, solidCellCount_temp)
            
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n, m, i, j, n1, n2, n3, n4, i1, j1, nn, &
                               nel2n, sumId
        INTEGER            :: iPt, iPt1
        REAL (KIND=8)      :: n1x, n1y, dis , minDis, n1dotn,  &
                              cent_x, cent_y, sStart, sStop, ic1x, ic2x, ic3x, ic1y, ic2y, ic3y 
        INTEGER (KIND=8),intent(in)    :: ibElems, nx, ny
		INTEGER (KIND=8),intent(inout)    :: cell(:)
		REAL (KIND = 8), intent(in) :: xcent(:), ycent(:), cosAlpha(:), cosBeta(:), x1(:), y1(:) , xp(:), yp(:) 
		INTEGER ,intent(inout) :: ibCellCount, fluidCellCount, solidCellCount
		INTEGER (KIND = 8), intent(inout)  :: nodeId(:,:)
		INTEGER (KIND = 8), intent(in) :: interceptedIndexPtr(:,:)
		INTEGER ,intent(out) :: ibCellCount_temp, fluidCellCount_temp, solidCellCount_temp
		ibCellCount_temp=ibCellCount
		fluidCellCount_temp=fluidCellCount
		solidCellCount_temp= solidCellCount
		
        DO nn = 1, ibCellCount_temp  
           i1 = interceptedIndexPtr(nn, 1)
           j1 = interceptedIndexPtr(nn, 2)
           DO 10 j = j1-2, j1+3
           DO 10 i = i1-2, i1+3
              nodeId(i,j) = 0
              m = 0
              minDis = 1e14
              n1x = x1(i)
              n1y = y1(j)      
              n1dotn = 0
	           DO m = 1, ibElems
                 cent_x = xcent(m)
                 cent_y = ycent(m)
                 dis  = dsqrt( (n1y-cent_y)**2 + (n1x-cent_x)**2 )     
		           IF (dis.LT.minDis) THEN
                    minDis   = dis                 
                    nel2n    = m
                 ENDIF
              ENDDO
              !print*, i, j, minDis
	           n1dotn = (n1x - xcent(nel2n))*cosAlpha(nel2n) + (n1y - ycent(nel2n))*cosBeta(nel2n)
              IF (n1dotn.LE.1e-14) THEN
                 nodeId(i,j) = 1  
              !ELSEIF (abs(n1dotn).LT.1e-16.AND.ABS(n1dotn).GT.0) THEN
                 !nodeId((j-1)*(nx+3) + i) = 1  
              !ELSEIF (abs(n1dotn).EQ.0) THEN
                 !nodeId((j-1)*(nx+3) + i) = 1  
              ENDIF
 10        CONTINUE
	      ENDDO
	      ibCellCount_temp = 0  
	      solidCellCount_temp = 0
	      fluidCellCount_temp = 0
         DO 20 j = 2, ny+1
         DO 20 i = 2, nx+1
            n   = i-1  + (nx)*(j-2)
			
            sumId = 0
            sumId = nodeId(i,j) + nodeId(i+1,j) + nodeId(i,j+1) + nodeId(i+1,j+1)
            IF (sumId.eq.4) THEN
		         cell(n) = 1
		         solidCellCount_temp = solidCellCount_temp + 1
            ELSEIF (sumId.eq.0) THEN
               cell(n) = 0  
               fluidCellCount_temp  = fluidCellCount_temp + 1 
            ELSEIF (sumId.NE.0.and.sumId.NE.4) THEN 
               cell(n) = 2               
               ibCellCount_temp = ibCellCount_temp + 1 
			   
            ENDIF 
 20      CONTINUE  
    GOTO 100
         print*, 'search done', ibCellCount_temp, fluidCellCount_temp, solidCellCount_temp
         ibCellCount_temp = 0
         DO 30 j = 2, ny+1
         DO 30 i = 2, nx+1
            n   = i-1  + (nx)*(j-2)  
            IF (Cell(n).eq.2) THEN
              IF (cell(n+1).eq.0.OR.cell(n-1).eq.0.OR.cell(n+nx).eq.0.OR.cell(n-nx).eq.0) THEN      
                 ibCellCount_temp = ibCellCount_temp + 1   
              ELSE
                 cell(n) = 1  
                 solidCellCount_temp = solidCellCount_temp + 1 
              ENDIF
            ENDIF
 30      CONTINUE

 
        open(82,file='inner1.dat',status='unknown')
        write(82,*)'variables = "x", "y","z"'!, "dnorm","ibmloc"'
       !do i=1,iim
       do j=2,ny+1
       do i=2,nx+1
       n   = i-1  + (nx)*(j-2) 
	   
       if(cell(n).eq.2)then
       write(82,*)xp(i),yp(j), 0
        !+ ,ibmloc(i,j,k)
       endif
       end do
       end do
      ! end do
      !stop
        print*, 'search done', ibCellCount_temp, fluidCellCount_temp, solidCellCount_temp
        100 CONTINUE         
        print*, 'selective retagging', ibCellCount_temp, fluidCellCount_temp, solidCellCount_temp
     END SUBROUTINE selective_retagging
     
! !*************************************************
 SUBROUTINE cell_count(nx, ny, ibCellCount, fluidCellCount, solidCellCount, cell, interceptedIndexPtr, pNormDis, nel2p, &
					  u1NormDis, u2NormDis , v1NormDis, v2NormDis, fluidIndexPtr, solidIndexPtr)
                
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n, iPt, iPt1, iPt2, i, j
		INTEGER (KIND=8),intent(in)    :: nx, ny, ibCellCount, fluidCellCount, solidCellCount, cell(:)
		INTEGER (KIND = 8), intent(out) :: interceptedIndexPtr(ibCellCount,2), fluidIndexPtr(fluidCellCount, 2), &
										   solidIndexPtr(solidCellCount, 2), nel2p(ibCellCount)
		REAL(KIND = 8), intent(out)::pNormDis(ibCellCount), u1NormDis(ibCellCount), u2NormDis(ibCellCount) , &
									 v1NormDis(ibCellCount), v2NormDis(ibCellCount)
      
         WRITE(199,*)'variables= "x","y","z","cellid"'
         WRITE(199,*) 'zone, ', 'i = ', nx,' j = ', ny, ' k = ', 1
         iPt  = 0
         iPt1 = 0
         iPt2 = 0
        
         DO 30 j = 2, ny+1
         DO 30 i = 2, nx+1
            n   = i-1  + (nx)*(j-2)
            !print*, i, j
            IF (cell(n).eq.0) THEN
               iPt1 = iPt1 + 1
               fluidIndexPtr(iPt1, 1) = i
               fluidIndexPtr(iPt1, 2) = j
            ELSEIF (cell(n).eq.1) THEN
               iPt2 = iPt2 + 1
               solidIndexPtr(iPt2, 1) = i
               solidIndexPtr(iPt2, 2) = j  
            ELSEIF (cell(n).eq.2) THEN
               iPt = iPt + 1
               interceptedIndexPtr(iPt, 1) = i
               interceptedIndexPtr(iPt, 2) = j
            ENDIF 
            
            !write(199,*) xp(i), yp(j), 0, cell(n)
 30      CONTINUE     
         !stop
 END SUBROUTINE cell_count


! !**************************************************************************

     ! SUBROUTINE compute_norm_distance(ibCellCount,interceptedIndexPtr, xcent, ycent, cosAlpha, &
									  ! cosBeta,v2NormDis,v1NormDis, u1NormDis, u2NormDis, pNormDis,nel2p)      
        ! INTEGER (KIND=8),intent(in)    :: ibCellCount
		! INTEGER (KIND = 8), intent(in) :: interceptedIndexPtr(:,:)										   
		! INTEGER (KIND = 8), intent(inout) ::  nel2p(:)
		! REAL(KIND = 8), intent(inout)::pNormDis(:), u1NormDis(:), u2NormDis(:), v2NormDis(:),v1NormDis(:)     
        ! INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        ! INTEGER (kind = 8) ::  n, m, i, j, n1, n2, n3, n4, iNew, jNew, &
                               ! nel2n, sumId, ar_intcptd
        ! INTEGER            :: k, iPt
        ! REAL (KIND=8)      :: n3x1, n3x2, n3y1, n3y2, n3x, n3y, dis , minDis, n1dotn,  &
                              ! cent_x, cent_y, sStart, sStop, ic1x, ic2x, ic3x, ic1y, ic2y, ic3y 
       ! CHARACTER*70  filename1
    
        ! nel2p = 0
            
        ! DO 10 k = 1, ibCellCount
           ! m = 0
           ! minDis = 1e14    
           ! n3x = xp(interceptedIndexPtr(k, 1))
           ! n3y = yp(interceptedIndexPtr(k, 2))
           ! n1dotn = 0
           ! nel2n = 0
           ! !print*,  normDisPtr(k, 2), n3x, normDisPtr(k, 2), n3y
	        ! DO m = 1, ibElems
              ! cent_x = xcent(m)
              ! cent_y = ycent(m)
              ! dis  = dsqrt( (n3y-cent_y)**2 + (n3x-cent_x)**2 )     
	           ! IF (dis.LT.minDis) THEN
                 ! minDis   = dis                 
                 ! nel2n    = m
              ! ENDIF
           ! ENDDO
           
           ! nel2p(k) = nel2n
           ! pNormDis(k) = (n3x - xcent(nel2n))*cosAlpha(nel2n) + (n3y - ycent(nel2n))*cosBeta(nel2n)
           
           ! n3y  = yu(interceptedIndexPtr(k, 2))
           ! n3x1 = xu(interceptedIndexPtr(k, 1))
           ! n3x2 = xu(interceptedIndexPtr(k, 1)+1)
           ! u1NormDis(k) = (n3x1 - xcent(nel2n))*cosAlpha(nel2n) + (n3y - ycent(nel2n))*cosBeta(nel2n)          
           ! u2NormDis(k) = (n3x2 - xcent(nel2n))*cosAlpha(nel2n) + (n3y - ycent(nel2n))*cosBeta(nel2n)
           ! n3x  = xv(interceptedIndexPtr(k, 1))
           ! n3y1 = yv(interceptedIndexPtr(k, 2))
           ! n3y2 = yv(interceptedIndexPtr(k, 2)+1)                      
           ! v1NormDis(k) = (n3x - xcent(nel2n))*cosAlpha(nel2n) + (n3y1 - ycent(nel2n))*cosBeta(nel2n)           
           ! v2NormDis(k) = (n3x - xcent(nel2n))*cosAlpha(nel2n) + (n3y2 - ycent(nel2n))*cosBeta(nel2n)
           ! !print*, pNormDis(k), u1NormDis(k),  u2NormDis(k), v1NormDis(k), v2NormDis(k)
 ! 10      CONTINUE
	
         ! print*, 'computeNormDistance done'
        ! !stop
             
     ! END SUBROUTINE compute_norm_distance




     SUBROUTINE compute_norm_distance(ibElems, ibCellCount, xp, yp, nel2p, xcent, ycent, pNormDis, cosAlpha, cosBeta, yu, xu,&
									interceptedIndexPtr, u1NormDis, u2NormDis, xv, yv, v1NormDis, v2NormDis)
              
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n, m, i, j, n1, n2, n3, n4, iNew, jNew, &
                               nel2n, sumId, ar_intcptd
        INTEGER            :: k, iPt
        REAL (KIND=8)      :: n3x1, n3x2, n3y1, n3y2, n3x, n3y, dis , minDis, n1dotn,  &
                              cent_x, cent_y, sStart, sStop, ic1x, ic2x, ic3x, ic1y, ic2y, ic3y 
		INTEGER (KIND=8),intent(in)    :: ibCellCount, ibElems
		INTEGER (KIND = 8), intent(in) :: interceptedIndexPtr(:,:)										   
		INTEGER (KIND = 8), intent(inout) ::  nel2p(:)
		REAL(KIND = 8), intent(in)::xp(:), yp(:),xcent(:), ycent(:)
		REAL(KIND = 8), intent(inout)::pNormDis(:), u1NormDis(:), u2NormDis(:), v2NormDis(:),v1NormDis(:),&
									   cosAlpha(:), cosBeta(:),xu(:), yu(:), xv(:), yv(:)
       CHARACTER*70  filename1
    
        nel2p = 0
            
        DO 10 k = 1, ibCellCount
           m = 0
           minDis = 1e14    
           n3x = xp(interceptedIndexPtr(k, 1))
           n3y = yp(interceptedIndexPtr(k, 2))
           n1dotn = 0
           nel2n = 0
           !print*,  normDisPtr(k, 2), n3x, normDisPtr(k, 2), n3y
	        DO m = 1, ibElems
              cent_x = xcent(m)
              cent_y = ycent(m)
              dis  = dsqrt( (n3y-cent_y)**2 + (n3x-cent_x)**2 )     
	           IF (dis.LT.minDis) THEN
                 minDis   = dis                 
                 nel2n    = m
              ENDIF
           ENDDO
           
           nel2p(k) = nel2n
           pNormDis(k) = (n3x - xcent(nel2n))*cosAlpha(nel2n) + (n3y - ycent(nel2n))*cosBeta(nel2n)
           
           n3y  = yu(interceptedIndexPtr(k, 2))
           n3x1 = xu(interceptedIndexPtr(k, 1))
           n3x2 = xu(interceptedIndexPtr(k, 1)+1)
           u1NormDis(k) = (n3x1 - xcent(nel2n))*cosAlpha(nel2n) + (n3y - ycent(nel2n))*cosBeta(nel2n)          
           u2NormDis(k) = (n3x2 - xcent(nel2n))*cosAlpha(nel2n) + (n3y - ycent(nel2n))*cosBeta(nel2n)
           n3x  = xv(interceptedIndexPtr(k, 1))
           n3y1 = yv(interceptedIndexPtr(k, 2))
           n3y2 = yv(interceptedIndexPtr(k, 2)+1)                      
           v1NormDis(k) = (n3x - xcent(nel2n))*cosAlpha(nel2n) + (n3y1 - ycent(nel2n))*cosBeta(nel2n)           
           v2NormDis(k) = (n3x - xcent(nel2n))*cosAlpha(nel2n) + (n3y2 - ycent(nel2n))*cosBeta(nel2n)
           !print*, pNormDis(k), u1NormDis(k),  u2NormDis(k), v1NormDis(k), v2NormDis(k)
 10      CONTINUE
	
         print*, 'computeNormDistance done'
		 print*, 'n3x1',n3x1
		 print*, 'nel2n',nel2n
		 print*, 'n3y',n3y
        !stop
             
     END SUBROUTINE compute_norm_distance
	 
	 
	 ! !*************************************************
 SUBROUTINE cell_count_2(nx, ny, cell, interceptedIndexPtr, &
					   fluidIndexPtr, solidIndexPtr)
                
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n, iPt, iPt1, iPt2, i, j
		INTEGER (KIND=8),intent(in)    :: nx, ny, cell(:)
		INTEGER (KIND = 8), intent(inout) :: interceptedIndexPtr(:,:), fluidIndexPtr(:,:), &
										     solidIndexPtr(:,:)
		
      
         WRITE(199,*)'variables= "x","y","z","cellid"'
         WRITE(199,*) 'zone, ', 'i = ', nx,' j = ', ny, ' k = ', 1
         iPt  = 0
         iPt1 = 0
         iPt2 = 0
        
         DO 30 j = 2, ny+1
         DO 30 i = 2, nx+1
            n   = i-1  + (nx)*(j-2)
            !print*, i, j
            IF (cell(n).eq.0) THEN
               iPt1 = iPt1 + 1
               fluidIndexPtr(iPt1, 1) = i
               fluidIndexPtr(iPt1, 2) = j
            ELSEIF (cell(n).eq.1) THEN
               iPt2 = iPt2 + 1
               solidIndexPtr(iPt2, 1) = i
               solidIndexPtr(iPt2, 2) = j  
            ELSEIF (cell(n).eq.2) THEN
               iPt = iPt + 1
               interceptedIndexPtr(iPt, 1) = i
               interceptedIndexPtr(iPt, 2) = j
            ENDIF 
            
            !write(199,*) xp(i), yp(j), 0, cell(n)
 30      CONTINUE     
         !stop
 END SUBROUTINE cell_count_2
! !***********************************************************   
      SUBROUTINE compute_surface_norm_2(ibElems,xnode, ynode, ibElP1, ibElP2, xcent, ycent, cosAlpha, cosBeta, area)
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) ::  n!, c1, c2, c3, c4
        REAL (KIND=8)      :: p1x, p1y, p2x, p2y, lenEL, inor
		INTEGER (KIND=8),intent(in)    :: ibElP1(:), ibElP2(:), ibElems
		REAL (KIND = 8), intent(in) :: xnode(:), ynode(:)
		REAL (KIND = 8), intent(inout) :: xcent(:), ycent(:), cosAlpha(:), cosBeta(:), area(:)
                             
        
        inor = -1._rk
        !compute centroid and direction cosines           
        DO n = 1, ibElems
           p1x = xnode(ibElP1(n))                       !x coordinate element node 1
           p1y = ynode(ibElP1(n))                       !y coordinate element node 1
           p2x = xnode(ibElP2(n))                       !x coordinate element node 2
           p2y = ynode(ibElP2(n))                       !y coordinate element node 2
           lenEL = dsqrt((p2y-p1y)**2 + (p2x-p1x)**2)   !length of element
           area(n) = lenEl
           xcent(n) =  (p2x+p1x)/2._rk                  !centroid x coordinate element
           ycent(n) =  (p2y+p1y)/2._rk                  !centroid y coordinate element          
           cosAlpha(n) = -(p2y-p1y)/lenEl*inor               !direction cosine unit normal along x
           cosBeta(n)  = (p2x-p1x)/lenEl*inor               !direction cosine unit normal along y
           !WRITE(*,*) xcent(n), ycent(n), cosAlpha(n)
           WRITE(196,*) xcent(n), ycent(n), cosBeta(n)
		   !print*,'-if',xcent(n)
        ENDDO
        !print*, 'ibel done'     
       !stop
     END SUBROUTINE compute_surface_norm_2

! !**************************************************************************