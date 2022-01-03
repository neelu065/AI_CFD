!csssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE write_output(ita, nx, ny, xp, yp, u, v, p, totime, cell)
       
       INTEGER, PARAMETER :: rk = selected_real_kind(8) 
       CHARACTER*70  filename1
       CHARACTER*70  filename2
       INTEGER  :: n, k, i, j  
	   REAL (KIND = 8), intent(in)   ::p(:,:), v(:,:), u(:,:), totime,&
									   xp(:), yp(:)
	   INTEGER (KIND = 8), intent(in) :: ita, nx, ny, cell(:)	   
         IF(ita.GT.0.AND.mod(ita,100).EQ.0) THEN        
        
            WRITE(filename1,1) ita
 1          FORMAT('out/fielddata.',i9.9,".dat")
            OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')
            WRITE(786,*)'variables= "x","y","z","u","v","p","totime","cellid"'
            WRITE(786,*) 'zone, ', 'i = ', nx,' j = ', ny, ' k = ', 1
            k = 1
            DO 30 j = 2, ny+1
            DO 30 i = 2, nx+1
               n    = i-1 + nx*(j-2)
               WRITE(786,*) xp(i), yp(j) , 0, &
               0.5*(u(i,j)+u(i-1,j)),  0.5*(v(i,j)+v(i,j-1)),  p(i,j),   totime, cell(n)           
 30         CONTINUE 
            CLOSE(786)      
         ENDIF        
      END SUBROUTINE write_output
 ! !csssssssssssssssssssssssssssssssssssssssssssssssssssss
      ! SUBROUTINE writeInlineData
       ! USE global
       ! INTEGER, PARAMETER :: rk = selected_real_kind(8) 
       ! CHARACTER*70  filename1
       ! CHARACTER*70  filename2
       ! INTEGER  :: n, k, i, j       
         ! IF(mod(ita,200).eq.0.OR.mod(ita,201).EQ.0) THEN        
        
            ! WRITE(filename1,1) ita
 ! 1          FORMAT('inlineData/linedata_u1_.',i9.9,".dat")
            ! OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')            
            ! WRITE(filename2,2) ita
 ! 2          FORMAT('inlineData/linedata_u2_.',i9.9,".dat")
            ! OPEN(UNIT = 787, FILE = filename2, STATUS = 'unknown')
            ! DO j = 1, ny+2
               ! WRITE(786,*) 0.5_rk*(u(304,j)+u(378,j)), 0.5_rk*(u(401,j)+u(402,j))            
               ! WRITE(787,*) 0.5_rk*(u(425,j)+u(426,j)), 0.5_rk*(u(449,j)+u(450,j))            
            ! ENDDO
            ! CLOSE(786)      
            ! CLOSE(787)
            ! WRITE(filename1,3) ita
 ! 3          FORMAT('inlineData/linedata_v1_.',i9.9,".dat")
            ! OPEN(UNIT = 786, FILE = filename1, STATUS = 'unknown')            
            ! WRITE(filename2,4) ita
 ! 4          FORMAT('inlineData/linedata_v2_.',i9.9,".dat")
            ! OPEN(UNIT = 787, FILE = filename2, STATUS = 'unknown')
            ! DO j = 1, ny+2
               ! WRITE(786,*) 0.5_rk*(v(377,j)+v(378,j)), 0.5_rk*(v(401,j)+v(402,j))            
               ! WRITE(787,*) 0.5_rk*(v(425,j)+v(426,j)), 0.5_rk*(v(449,j)+v(450,j))            
            ! ENDDO
            ! CLOSE(786)      
            ! CLOSE(787)
         ! ENDIF        
      ! END SUBROUTINE writeInlineData
    
! !csssssssssssssssssssssssssssss
    
! !csssssssssssssssssssssssssssssssssssssssssssssssssssss
      ! SUBROUTINE checkPoint
         ! USE global
         ! INTEGER, PARAMETER :: rk = selected_real_kind(8)  
         ! CHARACTER*70  filename3
         ! IF (mod(ita,50000).eq.0) THEN 
	     ! WRITE(filename3,3) ita          
 ! 3    	     FORMAT('checkPoint/checkData.',i9.9,".dat")
            ! open(UNIT = 790, FILE = filename3, STATUS = 'unknown')
            ! DO j = 1, ny+2
            ! DO i = 1, nx+2
               ! write(790, *) u(i,j), v(i,j), p(i, j), totime
            ! END DO
            ! ENDDO
            ! CLOSE(790)      
        ! ENDIF
     ! END SUBROUTINE
     
! !csssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE fft_data(u, totime, ita, v, p)
         
         INTEGER, PARAMETER :: rk = selected_real_kind(8)  
         CHARACTER*70  filename3
         INTEGER  :: il_x1, il_x2, il_x3, jl_y1, jl_y2, jl_y3
		 REAL (KIND = 8), intent(in)   ::p(:,:), v(:,:), u(:,:), totime
		 INTEGER (KIND = 8), intent(in) :: ita
         IF (ita.GT.50000) THEN 
            il_x1 = 212
            il_x2 = 222
            il_x3 = 232
            jl_y1 = 192
            jl_y2 = 202
            jl_y3 = 212
	     WRITE(filename3,2)       
 2    	     FORMAT('FFT/fftData_u'".dat")
            OPEN(UNIT = 791, FILE = filename3, ACCESS = 'Append', STATUS = 'unknown')
            WRITE(791,*) u(il_x1, jl_y1), u(il_x1, jl_y2), u(il_x1, jl_y3), &
                         u(il_x2, jl_y1), u(il_x2, jl_y2), u(il_x2, jl_y3), &
                         u(il_x3, jl_y1), u(il_x3, jl_y2), u(il_x3, jl_y3), totime             
            CLOSE(792)
            WRITE(filename3,3)          
 3    	     FORMAT('FFT/fftData_v'".dat")
            OPEN(UNIT = 792, FILE = filename3, ACCESS = 'Append', STATUS = 'unknown')
            WRITE(792,*) v(il_x1, jl_y1), v(il_x1, jl_y2), v(il_x1, jl_y3), &
                         v(il_x2, jl_y1), v(il_x2, jl_y2), v(il_x2, jl_y3), &
                         v(il_x3, jl_y1), v(il_x3, jl_y2), v(il_x3, jl_y3), totime             
            CLOSE(793)
            GOTO 100
            WRITE(filename3,4)           
 4    	     FORMAT('FFT/fftData_p'".dat")
            OPEN(UNIT = 793, FILE = filename3, ACCESS = 'Append', STATUS = 'unknown')
            WRITE(793,*) p(il_x1, jl_y1), p(il_x1, jl_y2), p(il_x1, jl_y3), &
                         p(il_x2, jl_y1), p(il_x2, jl_y2), p(il_x2, jl_y3), &
                         p(il_x3, jl_y1), p(il_x3, jl_y2), p(il_x3, jl_y3), totime             
            CLOSE(793)
 100        CONTINUE      
        ENDIF
     END SUBROUTINE
!cssssssssssssssssssssssssssssssssssssssssssssssssssssssss
       SUBROUTINE stress_2d(ibElems, nx, ny, xcent, ycent, x1, y1, cosAlpha,&
						   cosBeta, xu, yu, xv, yv, xp, yp, ut, vt, re, area, p, ita, totime)
       
       INTEGER:: i, j, ielem, i_x, j_y
       REAL (KIND = 8), intent(in) :: xcent(:), ycent(:), cosAlpha(:), cosBeta(:), area(:), totime, &
									  re, x1(:), y1(:), xu(:), yu(:), xv(:), yv(:), xp(:), yp(:), vt(:, :),&
									  ut(:, :), p(:, :)
	   INTEGER (KIND = 8), intent(in) :: ita, ibElems, nx, ny
	   
       INTEGER, DIMENSION(ibElems+1):: ielem_i_cell, jelem_j_cell

       REAL(KIND = 8):: lxx, dx1pt, dy1pt, p1x, p1y, p11, C1_P

       REAL(KIND = 8), DIMENSION(ibElems+1):: xpnt1, xpnt2, xpnt3, & 
     & xpnt4, ypnt1, ypnt2, ypnt3, ypnt4, xpnt2new, ypnt2new 

       REAL(KIND = 8), DIMENSION(ibElems+1):: xvelpnt2, yvelpnt2
    ! & xvelpnt3, yvelpnt3
       
       REAL(KIND = 8), DIMENSION(ibElems+1):: xprpnt2, yprpnt2, &
     & xprpnt3, yprpnt3
    
       REAL(KIND = 8), DIMENSION(ibElems+1):: u_pnt1, v_pnt1, u_pnt2, & 
     & v_pnt2, b1, b2, P_pnt2, P_pnt3

       !REAL(KIND = 8), DIMENSION(ibElems+1):: ustar_pnt2, ustar_pnt3, 
     !& vstar_pnt2, vstar_pnt3, b3, u_pnt3, v_pnt3
       
       REAL(KIND = 8), DIMENSION(ibElems+1):: DN_pnt12, del_X, del_Y

       REAL(KIND = 8), DIMENSION(ibElems+1):: DNpr_pnt2, DNpr_pnt3

       !REAL(KIND = 8), DIMENSION(ibElems+1):: DNv_pnt2, DNv_pnt3

       REAL(KIND = 8), DIMENSION(ibElems+1)::  shear_t_force, &
     & shear_x_force, shear_y_force, p_surf, f_surf, f_surf_x, f_surf_y

       REAL(KIND = 8):: area_Sx, area_Px, cx1, cy1, xlast, ylast 

       !REAL(KIND = 8), DIMENSION(nx+1)::  x1, xp, xu, xv

       !REAL(KIND = 8), DIMENSION(ny+1)::  y1, yp, yu, yv

       REAL(KIND = 8):: pressure_drag,viscous_drag,viscousLift,PressureLift,surf_area, &
     & viscousdragcoefficient,pressuredragcoefficient,ViscousLiftcoefficient, &
     & PressureLiftcoefficient,TotalLiftcoefficient,alternativeLiftcofficient
     
       CHARACTER(LEN=70):: filename1, filename2, filename3, filename4,filename5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!       OPEN(unit=787,file='velocity.dat',status='unknown')
!       DO j = 2, ny+1-1
!       DO i = 2, nx+1-1
!       READ(787,*) ut(i,j), vt(i,j), p(i,j)
!       END DO
!       END DO
!       CLOSE(787)         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
       
!declare global deltay, deltax, nx+1, ny+1, ibElems, area(ibElems), cosAlpha(ibElems), cosBeta(ibElems)
!!!!!!!!!!!!!!!!!!!!!!!!!!!normal grid!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !x1(1) = 0.0
        !xlast = x1(1)
        !DO i = 2, nx+1 
       ! x1(i) = xlast + deltax
        !xlast = x1(i)
        !END DO
        
        !y1(1) = 0.0
       ! ylast = y1(1)
        !DO j = 2, ny+1
        !y1(j)= ylast + deltay
        !ylast = y1(j)
        !END DO        
!!!!!!!!!!!!!!!!!!!!!!!!!!!normal grid part end!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!pressure grid start!!!!!!!!!!!!!!!!!!!!!!!!!
        !DO j = 1, ny
        !DO i = 1, nx
       ! xp(i) = 0.5*(x1(i)+x1(i+1))
        !yp(j) = 0.5*(y1(j)+y1(j+1))
        !END DO
        !END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pressure grid end!!!!!!!!!!!!!!!!!!!!!!!!!!
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!u-velocity grid start!!!!!!!!!!!!!!!!!!!!!!
	     !DO i = 1, nx+1
	      !xu(i) = x1(i)
	    ! END DO

	     !DO j = 1, ny
        !yu(j) = 0.5*(y1(j)+y1(j+1))
	     !END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!u-velocity grid end!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!v-velocity grid start!!!!!!!!!!!!!!!!!!!!!!!!!
	     !DO j = 1, ny+1
	      !yv(j) = y1(j) 
	    ! END DO
	
	    ! DO i = 1, nx
       ! xv(i)=0.5*(x1(i)+x1(i+1))
	    ! END DO
!!!!!!!!!!!!!!!!!!!!!!!!v-velocity grid end!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !IF (ita.GT.20000) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
       DO ielem = 1, ibElems

       do i = 2, nx+1
       if((xcent(ielem).ge.x1(i)).and.(xcent(ielem).lt.x1(i+1)))then
       ielem_i_cell(ielem) = i
       end if
       end do

       do j = 2, ny+1
       if((ycent(ielem).ge.y1(j)).and.(ycent(ielem).lt.y1(j+1)))then
       jelem_j_cell(ielem) = j
       end if
       end do

       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!FIND POINT 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO ielem = 1, ibElems

       xpnt1(ielem) = xcent(ielem)
       ypnt1(ielem) = ycent(ielem)

       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF POINT 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!find velocity point 2!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO ielem = 1, ibElems

       del_X(ielem) = x1(ielem_i_cell(ielem)+1)-x1(ielem_i_cell(ielem))
       del_Y(ielem) = y1(jelem_j_cell(ielem)+1)-y1(jelem_j_cell(ielem))

       DN_pnt12(ielem) = dsqrt(del_X(ielem)**2 + del_Y(ielem)**2)


       
       lxx = DN_pnt12(ielem)
       
       !lxx = Sqrt(deltax**2 + deltay**2)
       
       xvelpnt2(ielem) = xpnt1(ielem)+lxx*cosAlpha(ielem)
       yvelpnt2(ielem) = ypnt1(ielem)+lxx*cosBeta(ielem)
       
       xprpnt2(ielem) = xpnt1(ielem)+lxx*cosAlpha(ielem)
       yprpnt2(ielem) = ypnt1(ielem)+lxx*cosBeta(ielem)
       
       lxx = 1.5*lxx
       
       xprpnt3(ielem) = xpnt1(ielem)+lxx*cosAlpha(ielem)
       yprpnt3(ielem) = ypnt1(ielem)+lxx*cosBeta(ielem)
    
   
       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!End velocity point 2!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!find pressure point 2, 3!!!!!!!!!!!!!!!!!!!!!!!!!
       !DO ielem = 1, ibElems

       

       !END DO
!!!!!!!!!!!!!!!!!!!!!!!end pressure point 2, 3!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!velocity at point 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO ielem = 1, ibElems

       u_pnt1(ielem) = 0.0
       v_pnt1(ielem) = 0.0

       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!velocity interpolation at point 2!!!!!!!!!!!!!!!!!
       DO ielem = 1, ibElems

!!!!!!!!!!!!!!!!!!!!!u velocity interpolation in point 2!!!!!!!!!!!!!!!!
       DO i = 2, nx
       if(xvelpnt2(ielem).ge.xu(i).and.xvelpnt2(ielem).lt.xu(i+1)) then
       i_x = i
       end if
       END DO
	
       DO j = 2, ny
       if(yvelpnt2(ielem).ge.yu(j).and.yvelpnt2(ielem).lt.yu(j+1)) then
       j_y = j
       end if
       END DO
     
       dx1pt = xvelpnt2(ielem)-xu(i_x)
       dy1pt = yvelpnt2(ielem)-yu(j_y)
	 
       p1x = ut(i_x,j_y) + (ut(i_x+1,j_y) &
     & -ut(i_x,j_y))*(dx1pt/(xu(i_x+1)-xu(i_x)))
     
       p1y = ut(i_x,j_y+1) + (ut(i_x+1,j_y+1) &
     & -ut(i_x,j_y+1))*(dx1pt/(xu(i_x+1)-xu(i_x)))

       p11 = p1x + (p1y-p1x)*(dy1pt/(yu(j_y+1)-yu(j_y)))
	 
       u_pnt2(ielem) = p11
!!!!!!!!!!!!!!!!!!end u velocity interpolation at point 2!!!!!!!!!!!!!!!
         
!!!!!!!!!!!!!!!!!!v velocity interpolation in point 2!!!!!!!!!!!!!!!!!!!	 
       DO i = 2, nx
       if(xvelpnt2(ielem).ge.xv(i).and.xvelpnt2(ielem).lt.xv(i+1)) then
       i_x = i
       end if
       END DO
	
       DO j = 2, ny
       if(yvelpnt2(ielem).ge.yv(j).and.yvelpnt2(ielem).lt.yv(j+1)) then
       j_y = j
       end if
       END DO
        
       dx1pt = xvelpnt2(ielem)-xv(i_x)
       dy1pt = yvelpnt2(ielem)-yv(j_y)
	 
       p1x = vt(i_x,j_y) + (vt(i_x+1,j_y) &
     & -vt(i_x,j_y))*(dx1pt/(xv(i_x+1)-xv(i_x)))
           
       p1y = vt(i_x,j_y+1) + (vt(i_x+1,j_y+1) &
     & -vt(i_x,j_y+1))*(dx1pt/(xv(i_x+1)-xv(i_x)))
       
       p11 = p1x + (p1y-p1x)*(dy1pt/(yv(j_y+1)-yv(j_y)))
	
       v_pnt2(ielem) = p11
!!!!!!!!!!!!!!!!!!!!!end v velocity interpolation at point 2!!!!!!!!!!!!
       
       END DO
!!!!!!!!!!!!!!!!!!!end velocity interpolation at point 2!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!velocity interpolation at point 3!!!!!!!!!!!!!!!!!!
!       DO ielem = 1, ibElems

!!!!!!!!!!!!!!!!!!!!!u velocity interpolation at point 3!!!!!!!!!!!!!!!!
!       DO i = 2, nx+1-1
!       if(xvelpnt3(ielem).ge.xu(i).and.xvelpnt3(ielem).lt.xu(i+1)) then
!       i_x = i
!       end if
!       END DO
	
!       DO j = 2, ny+1-1
!       if(yvelpnt3(ielem).ge.yu(j).and.yvelpnt3(ielem).lt.yu(j+1)) then
!       j_y = j
!       end if
!       END DO
     
!       dx1pt = abs(xvelpnt3(ielem)-xu(i_x))
!       dy1pt = abs(yvelpnt3(ielem)-yu(j_y))
	 
!       p1x = ut(i_x-1,j_y) + (ut(i_x,j_y) &
!     & -ut(i_x-1,j_y))*(dx1pt/(xu(i_x+1)-xu(i_x)))
     
!       p1y = ut(i_x-1,j_y+1) + (ut(i_x,j_y+1) &
!     & -ut(i_x-1,j_y+1))*(dx1pt/(xu(i_x+1)-xu(i_x)))

!       p11 = p1x + (p1y-p1x)*(dy1pt/(yu(j_y+1)-yu(j_y)))
	 
!       u_pnt3(ielem) = p11
!!!!!!!!!!!!!!!!!!end u velocity interpolation at point 3!!!!!!!!!!!!!!!
         
!!!!!!!!!!!!!!!!!!v velocity interpolation at point 3!!!!!!!!!!!!!!!!!!!	 
!       DO i = 2, nx+1-1
!       if(xvelpnt3(ielem).ge.xv(i).and.xvelpnt3(ielem).lt.xv(i+1)) then
!       i_x = i
!       end if
!       END DO
	
!       DO j = 2, ny+1-1
!       if(yvelpnt3(ielem).ge.yv(j).and.yvelpnt3(ielem).lt.yv(j+1)) then
!       j_y = j
!       end if
!       END DO
        
!       dx1pt = abs(xvelpnt3(ielem)-xv(i_x))
!       dy1pt = abs(yvelpnt3(ielem)-yv(j_y))
	 
!       p1x = vt(i_x,j_y-1) + (vt(i_x+1,j_y-1) &
!     & -vt(i_x,j_y-1))*(dx1pt/(xv(i_x+1)-xv(i_x)))
           
!       p1y = vt(i_x,j_y) + (vt(i_x+1,j_y) &
!     & -vt(i_x,j_y))*(dx1pt/(xv(i_x+1)-xv(i_x)))
       
!       p11 = p1x + (p1y-p1x)*(dy1pt/(yv(j_y+1)-yv(j_y)))
	
!       v_pnt3(ielem) = p11
!!!!!!!!!!!!!!!!!!!!!end v velocity interpolation at point 3!!!!!!!!!!!!
       
!       END DO
!!!!!!!!!!!!!!!!!!!end velocity interpolation at point 3!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       DO ielem = 1, ibElems

!       DNv_pnt2(ielem) = sqrt((xvelpnt2(ielem)-xpnt1(ielem))**2 + &
!     & (yvelpnt2(ielem)-ypnt1(ielem))**2)
   
!       DNv_pnt3(ielem) = sqrt((xvelpnt3(ielem)-xpnt1(ielem))**2 + &
!     & (yvelpnt3(ielem)-ypnt1(ielem))**2) 

!       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!calculate dudn, dvdn!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!       DO ielem = 1, ibElems

!       ustar_pnt2(ielem) = u_pnt2(ielem)-u_pnt1(ielem)
!       ustar_pnt3(ielem) = u_pnt3(ielem)-u_pnt1(ielem)
!       vstar_pnt2(ielem) = v_pnt2(ielem)-v_pnt1(ielem)
!       vstar_pnt3(ielem) = v_pnt3(ielem)-v_pnt1(ielem)

!       b1(ielem) = (DNv_pnt3(ielem)**2*ustar_pnt2(ielem) - &
!     & DNv_pnt2(ielem)**2*ustar_pnt3(ielem))/ &
!     & (DNv_pnt3(ielem)*DNv_pnt2(ielem)*(DNv_pnt3(ielem)-DNv_pnt2(ielem)))
   
!       b2(ielem) = (DNv_pnt3(ielem)**2*vstar_pnt2(ielem) - &
!     & DNv_pnt2(ielem)**2*vstar_pnt3(ielem))/ &
!     & (DNv_pnt3(ielem)*DNv_pnt2(ielem)*(DNv_pnt3(ielem)-DNv_pnt2(ielem)))

!       END DO
!!!!!!!!!!!!!!!!!!!!!!!end of calculation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO ielem = 1, ibElems

       DN_pnt12(ielem) = dsqrt((xvelpnt2(ielem)-xpnt1(ielem))**2 + &
     & (yvelpnt2(ielem)-ypnt1(ielem))**2)
   
       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!calculate dudn, dvdn!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
       DO ielem = 1, ibElems

       b1(ielem) = (u_pnt2(ielem)-u_pnt1(ielem))/DN_pnt12(ielem)
   
       b2(ielem) = (v_pnt2(ielem)-v_pnt1(ielem))/DN_pnt12(ielem)
         
       END DO
!!!!!!!!!!!!!!!!!!!!!!!end of calculation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!calculation of shear forces!!!!!!!!!!!!!!!!!!!!!!
       area_Sx = 0.0 
       DO ielem = 1, ibElems   

       cx1 = b1(ielem) - &
     & (b1(ielem)*cosAlpha(ielem)+b2(ielem)*cosBeta(ielem))*cosAlpha(ielem)

       cy1 = b2(ielem) - &
     & (b1(ielem)*cosAlpha(ielem)+b2(ielem)*cosBeta(ielem))*cosBeta(ielem)

       !shear_t_force(ielem)= (2./re)*sqrt(cx1**2 + cy1**2)*area(ielem)

       shear_x_force(ielem) = (2.0/re)*cx1*area(ielem)
       
       shear_y_force(ielem) = (2.0/re)*cy1*area(ielem)

       area_Sx = area_Sx + abs(cosBeta(ielem)*area(ielem))
     
       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!end of shear force calculation!!!!!!!!!!!!!!!!
       area_Sx = area_Sx/2.
!!!!!!!!!!!!!!!!!!!!!!!!pressure interpolation!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!Pressure interpolation at point 2!!!!!!!!!!!!!!!!!!
       DO ielem=1,ibElems

       DO i = 2, nx
       if(xprpnt2(ielem).ge.xp(i).and.xprpnt2(ielem).lt.xp(i+1)) i_x = i
       END DO
	
       DO j = 2, ny
       if(yprpnt2(ielem).ge.yp(j).and.yprpnt2(ielem).lt.yp(j+1)) j_y = j
       END DO
	
       dx1pt = abs(xprpnt2(ielem)-xp(i_x))
       dy1pt = abs(yprpnt2(ielem)-yp(j_y))
	 
       p1x = p(i_x,j_y) + &
     & (p(i_x+1,j_y)-p(i_x,j_y))*(dx1pt/(xp(i_x+1)-xp(i_x)))   
     
       p1y = p(i_x,j_y+1) + &
     & (p(i_x+1,j_y+1)-p(i_x,j_y+1))*(dx1pt/(xp(i_x+1)-xp(i_x)))
            
       p11 = p1x + (p1y-p1x)*(dy1pt/(yp(j_y+1)-yp(j_y)))
       	 
       P_pnt2(ielem) = p11
 
       END DO
!!!!!!!!!!!!!!!!!!!!!end of pressure interpolation at point 2!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!pressure interpolation at point 3!!!!!!!!!!!!!!!!!!!
       DO ielem=1,ibElems

       DO i = 2, nx
       if(xprpnt3(ielem).ge.xp(i).and.xprpnt3(ielem).lt.xp(i+1)) i_x = i
       END DO
	
       DO j = 2, ny
       if(yprpnt3(ielem).ge.yp(j).and.yprpnt3(ielem).lt.yp(j+1)) j_y = j
       END DO
		
       dx1pt = abs(xprpnt3(ielem)-xp(i_x))
       dy1pt = abs(yprpnt3(ielem)-yp(j_y))
      
       p1x = p(i_x,j_y) + &
     & (p(i_x+1,j_y)-p(i_x,j_y))*(dx1pt/(xp(i_x+1)-xp(i_x)))   
     
       p1y = p(i_x,j_y+1) + &
     & (p(i_x+1,j_y+1)-p(i_x,j_y+1))*(dx1pt/(xp(i_x+1)-xp(i_x)))
            
       p11 = p1x + (p1y-p1x)*(dy1pt/(yp(j_y+1)-yp(j_y)))
         	 
       P_pnt3(ielem) = p11

       END DO
!!!!!!!!!!!!!!!!!end of pressure interpolation at point 3!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO ielem = 1, ibElems

       DNpr_pnt2(ielem) = sqrt((xprpnt2(ielem)-xpnt1(ielem))**2 + &
     & (yprpnt2(ielem)-ypnt1(ielem))**2)
   
       DNpr_pnt3(ielem) = sqrt((xprpnt3(ielem)-xpnt1(ielem))**2 + &
     & (yprpnt3(ielem)-ypnt1(ielem))**2) 

       END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!pressure force calculation!!!!!!!!!!!!!!!!!!!!!!!!!!!
       area_Px = 0.0
       DO ielem = 1, ibElems

       !C1_P =  (P_pnt2(ielem)*DNpr_pnt3(ielem)**2 - &
     !& P_pnt3(ielem)*DNpr_pnt2(ielem)**2)/ &
     !& (DNpr_pnt3(ielem)**2-DNpr_pnt2(ielem)**2)
       C1_P =  (P_pnt2(ielem)*DNpr_pnt3(ielem) - &
     & P_pnt3(ielem)*DNpr_pnt2(ielem))/ &
     & (DNpr_pnt3(ielem)-DNpr_pnt2(ielem))

       p_surf(ielem) = C1_P

       f_surf(ielem) = C1_P*abs(area(ielem)) 

       f_surf_x(ielem) = -2.*f_surf(ielem)*cosAlpha(ielem)
       
       f_surf_y(ielem) = -2.*f_surf(ielem)*cosBeta(ielem)

       area_Px = area_Px + abs(cosAlpha(ielem)*area(ielem))

       END DO
!!!!!!!!!!!!!!!!!!!!!!end of pressure force calculation!!!!!!!!!!!!!!!!!
       area_Px = area_Px/2.

       

       !DO ielem = 1, ibElems
      

       !write(794,*) xpnt1(ielem),p_surf(ielem)
       
       !write(795,*) xpnt1(ielem),b1(ielem),b2(ielem)

       !END DO

       !close(794)
       !close(795)

!!!!!!!!!!!!!!!!!!!!!!!drag calculation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       pressure_drag = 0.
       viscous_drag = 0.
       surf_area = 0.
       viscousLift = 0.
       PressureLift = 0.
       
       !ViscousLiftcoefficient = 0.
       !PressureLiftcoefficient = 0.


       DO ielem = 1, ibElems
       
       viscous_drag = viscous_drag + shear_x_force(ielem)
       pressure_drag = pressure_drag + f_surf_x(ielem)
       viscousLift = viscousLift + shear_y_force(ielem)
       PressureLift = PressureLift + f_surf_y(ielem)
       surf_area = surf_area + abs(area(ielem))

       END DO       
      
       viscousdragcoefficient = viscous_drag/area_Sx
       pressuredragcoefficient = pressure_drag/area_Sx
       ViscousLiftcoefficient = viscousLift/area_Sx 
       PressureLiftcoefficient = PressureLift/area_Sx
       TotalLiftcoefficient = ViscousLiftcoefficient + &
     & PressureLiftcoefficient
       alternativeLiftcofficient = (viscousLift+PressureLift)/surf_area
!!!!!!!!!!!!!!!!!!!!!end of drag calculation!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !IF(mod(ita,1000).EQ.0) THEN
       !PRINT*, 'IN STRESS'
       !write(filename1,108) ita
 !108   format('Stress/pressuredata.',i9.9,".dat")
 
       !write(filename2,109) ita
 !109   format('Stress/veldata.',i9.9,".dat")
 
 
       !open(unit=794,file=filename1,status='unknown')
       !open(unit=795,file=filename2,status='unknown')
       write(filename3,110) ita
 110   format('Stress/stress',i9.9,".dat")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       open(unit=796,file=filename3,Access='Append',status='unknown')!!!!!!!!!!!!!!!!!!!!!!!!!       
       WRITE(796,*)  viscousdragcoefficient, pressuredragcoefficient, &
     & ViscousLiftcoefficient,  PressureLiftcoefficient, TotalLiftcoefficient, &
     & alternativeLiftcofficient, surf_area, area_Sx, totime, ita
     
       CLOSE(796)
       !PRINT*, 'IN STRESS1'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       write(filename4,111) 
 111   format('meanLift'".dat")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       open(unit=797,file=filename4,Access='Append',status='unknown')!!!!!!!!!!!!!!!!!!!!!!!!!       
       WRITE(797,*)   pressuredragcoefficient,totime, ita
     
       CLOSE(797)
       write(filename5,112) 
 112   format('drag'".dat")
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       open(unit=798,file=filename5,Access='Append',status='unknown')!!!!!!!!!!!!!!!!!!!!!!!!!       
       WRITE(798,*)  viscousdragcoefficient+pressuredragcoefficient, ita
     
       CLOSE(798)
       !PRINT*, 'IN STRESS2'
       END SUBROUTINE stress_2d
