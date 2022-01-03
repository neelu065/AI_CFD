subroutine read_input(input_filename, gridx_filename, gridy_filename,nx,ny,&
					  dt_order, a0, freq,	u0, v0, re,fx, fy, xShift, yShift,&
					  itamax, pcItaMax,itaSola, surGeoPoints,istart )
	   CHARACTER (LEN = 128), intent(in) ::input_filename, gridx_filename, gridy_filename
	   REAL (KIND =8),intent(out):: dt_order, a0, freq,	u0, v0, re,fx, fy, xShift, yShift
	   INTEGER ( KIND = 8), intent(out) :: itamax, pcItaMax,itaSola, surGeoPoints,istart
	   INTEGER ( KIND = 8), intent(out) :: nx,ny
         
	   
        OPEN(60, FILE = input_filename, FORM = 'formatted')
		
        READ(60,*) nx, ny,  &
                   fx, fy,  &
                   re,      &
                   itamax, itaSola, pcItaMax, &                                    
                   u0, v0,  &
                   surGeoPoints, a0, freq,  &
                   xShift, yShift,            &
                   istart, dt_order
        CLOSE(60) 
		!print*, 'nx',nx
	   
	 
		
end subroutine read_input
subroutine create_grid(input_filename, gridx_filename, gridy_filename,nx,ny,&
					  dt_order, a0, freq,	u0, v0, re,fx, fy, xShift, yShift,&
					  itamax, pcItaMax,itaSola, surGeoPoints,istart,&
					  deltax, deltay, x1, y1, xu, yu, xv, yv, xp, yp, alpha, deltat, pi )
	   CHARACTER (LEN = 128), intent(in) ::input_filename, gridx_filename, gridy_filename
	   REAL (KIND =8),intent(in):: dt_order, a0, freq,	u0, v0, re,fx, fy, xShift, yShift
	   INTEGER ( KIND = 8), intent(in) :: itamax, pcItaMax,itaSola, surGeoPoints,istart
	   INTEGER ( KIND = 8), intent(in) :: nx,ny
	   REAL (KIND =8),intent(out)::alpha, deltat, pi  
	   REAL (KIND = 8),intent(out) :: x1(nx+3), y1(ny+3), deltax(nx+2), deltay(ny+2), xu(nx+3),&
									  yu(ny+2), xv(nx+2), yv(ny+3), xp(nx+2), yp(ny+2)
	   INTEGER, PARAMETER :: rk = selected_real_kind(8)
       INTEGER (KIND=8) :: i, j
		OPEN(61, FILE = gridx_filename, FORM = 'formatted')
        DO i = 2, nx+2
           READ(61, *) x1(i)
        END DO
        CLOSE(61)  

        DO i = 2, nx+1
           deltax(i) = x1(i+1) - x1(i)
        END DO  
        
        deltax(1)    = deltax(2)
        deltax(nx+2) = deltax(nx+1)
        x1(1)        = x1(2) - deltax(1)
        x1(nx+3)     = x1(nx+2) + deltax(nx+2)

        OPEN(62, FILE = gridy_filename, FORM = 'formatted')
        DO i = 2, ny+2
           READ(62, *) y1(i)
        END DO
        CLOSE(62)  
        
        DO i = 2, ny+1
           deltay(i) = y1(i+1) - y1(i)
        END DO  
        
        deltay(1)    = deltay(2)
        deltay(ny+2) = deltay(ny+1)
        y1(1)        = y1(2) - deltay(1)
        y1(ny+3)     = y1(ny+2) + deltay(ny+2)
        
        DO i = 1, nx+3
           xu(i) = x1(i)
        ENDDO
        DO i = 1, ny+2
           yu(i) = 0.5_rk*(y1(i)+y1(i+1))
           yp(i) = yu(i)
        END DO
        
        DO i = 1, ny+3
           yv(i) = y1(i)
        ENDDO
        DO i = 1, nx+2
           xv(i) = 0.5_rk*(x1(i)+x1(i+1))
           xp(i) = xv(i)
        END DO
        DO i = 1, nx+2
        DO j = 1, ny+2
           !Print*, deltax(i),deltay(j)
        ENDDO   
        ENDDO
        pi = 4.D0*ATAN(1.D0)
        deltat = dt_order*0.002!5._rk/1200!pi/1200._rk!0.0041666666666_rk! dt_order*5e-4
        alpha  = 1._rk      
         print*, 'dt =',  deltat
end subroutine create_grid
!************************************************************************************************
subroutine read_geometry_nodes(geometry_filename,surGeoPoints,ibNodes, ibElems)
                  
       INTEGER (KIND = 8) :: n, i1, i2, i3, i4, i5, i6, i7, gPoints
       CHARACTER (LEN = 72) :: cLine
	   CHARACTER (LEN = 128) :: line
	   CHARACTER (LEN = 128), intent(in) ::geometry_filename
	   INTEGER ( KIND = 8), intent(in) :: surGeoPoints
	   INTEGER (KIND=8), intent(out)   :: ibElems, ibNodes
       !print*, 'hi'
       !OPEN(121, FILE = 'geometries/NACA0012.msh', form = 'formatted')               !READ SURFACE MESH FILE  
       OPEN(121, FILE = geometry_filename, form = 'formatted')               !READ SURFACE MESH FILE  
        DO n = 1, 4
           READ (121,*) cLine
        END DO
        READ (121,*) ibNodes  !nsurf=total no. of points in file  
        !ALLOCATE ( xnode(ibNodes), ynode(ibNodes) )
        DO n = 1, ibNodes
           READ (121,*) 
           !print*, i1, xnode(n), ynode(n)
        END DO
        DO n = 1, 2
          READ (121,*) line
        END DO        
        READ (121,*) ibElems   !no. of elements 
        DO n = 1, surGeoPoints
          READ (121,*) cLine
        END DO
        ibElems = ibElems-surGeoPoints
        
        
        
       !do n = 1, ibNodes
         !print*, n, xnode(n), ynode(n)
       !end do
       
       CLOSE(121)     
       PRINT *, 'SURFACE MESH READING COMPLETE'
       PRINT *, 'ibNodes =', ibNodes, 'ibElems =', ibElems
       
end subroutine read_geometry_nodes
   
subroutine read_surface_mesh(geometry_filename,surGeoPoints,ibNodes, ibElems, xnode, ynode,ibElP1, ibElP2)
                  
       INTEGER (KIND = 8) :: n, i1, i2, i3, i4, i5, i6, i7, gPoints
       CHARACTER (LEN = 72) :: cLine
	   CHARACTER (LEN = 128) :: line
	   CHARACTER (LEN = 128), intent(in) ::geometry_filename
	   INTEGER (KIND=8), intent(in)   :: ibElems, ibNodes, surGeoPoints
	   REAL (KIND = 8),intent(out) :: xnode(ibNodes), ynode(ibNodes),ibElP1(ibElems), ibElP2(ibElems)
       !print*, 'hi'
       !OPEN(121, FILE = 'geometries/NACA0012.msh', form = 'formatted')               !READ SURFACE MESH FILE  
       OPEN(121, FILE = geometry_filename, form = 'formatted')               !READ SURFACE MESH FILE  
        DO n = 1, 4
           READ (121,*) cLine
        END DO
        READ (121,*)   !nsurf=total no. of points in file  
        
        DO n = 1, ibNodes
           READ (121,*) i1, xnode(n), ynode(n), i2
           !print*, i1, xnode(n), ynode(n)
        END DO
        DO n = 1, 2
          READ (121,*) line
        END DO        
        READ (121,*)    !no. of elements 
        DO n = 1, surGeoPoints
          READ (121,*) cLine
        END DO
        !ibElems = ibElems-surGeoPoints
        
        ibElP1 = 0
        ibElP2 = 0
        DO n = 1, ibElems
          READ (121,*) i1, i2, i3, i4, i5, ibElP1(n), ibElP2(n)
          !Print*, n, ibElP1(n), ibElP2(n)
        END DO
       !do n = 1, ibNodes
         !print*, n, xnode(n), ynode(n)
       !end do
       
       CLOSE(121)     
       PRINT *, 'SURFACE MESH READING COMPLETE'
       PRINT *, 'ibNodes =', ibNodes, 'ibElems =', ibElems
       
end subroutine read_surface_mesh
      
!*******************************************************************

