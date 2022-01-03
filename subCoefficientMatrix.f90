      SUBROUTINE coefficient_matrix(nx, ny, Acx, Acy, xp, yp )     
                
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        INTEGER (kind = 8) :: n, i, j, k, iType
		INTEGER (kind = 8), intent(in) ::nx, ny
		REAL (KIND = 8), intent(inout)    ::Acx(:, :),&
											Acy(:, :), xp(:), yp(:)   
        REAL (KIND = 8)    ::  rx1, rx2, rxsum, ry1, ry2, rysum
         
         DO j = 2, ny+1
            ry1   = yp(j)   - yp(j-1)       
            ry2   = yp(j+1) - yp(j)     
            rysum = ry1 + ry2       
            Acy(j-1, 1) =   2._rk/(ry1*rysum)
            Acy(j-1, 2) =  -2._rk/(ry1*ry2)
            Acy(j-1, 3) =   2._rk/(ry2*rysum) 
            !print*, Acy(j-1,:)           
         END DO 

         DO i = 2, nx+1   
            rx1   = xp(i)   - xp(i-1)       
            rx2   = xp(i+1) - xp(i)     
            rxsum = rx1 + rx2                  
            Acx(i-1, 1) =   2._rk/(rx1*rxsum)
            Acx(i-1, 2) =  -2._rk/(rx1*rx2)
            Acx(i-1, 3) =   2._rk/(rx2*rxsum)    
            !print*, Acy(i-1,:)          
         END DO         
         !inlet, i = 1
         Acx(1, 2)   =  Acx(1, 2)  + Acx(1, 1)
         Acx(1, 1)   =  0._rk
         !outlet, i = nx
         Acx(nx, 2)  =  Acx(nx, 2) + Acx(nx, 3) 
         Acx(nx, 3)  =  0._rk 
         !bottom, j = 1
         Acy(1, 2)   =  Acy(1, 2)  + Acy(1, 1)
         Acy(1, 1)   =  0._rk
         !top, i = ny
         Acy(ny, 2)  =  Acy(ny, 2) + Acy(ny, 3) 
         Acy(ny, 3)  =  0._rk 
         !!$acc loop collapse(2) 
         print*, "Coefficient Matrix generated"

      END SUBROUTINE 
