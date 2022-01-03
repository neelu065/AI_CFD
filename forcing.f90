SUBROUTINE pressure_forcing_1(ibCellCount, interceptedIndexPtr, ita, nel2p,&
							  xcent, xp, pNormDis, xddot, p, deltax, deltay,&
							  cosAlpha, cosBeta, yp)
      
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: k, j, i, ii, il, jl, i_x1, i_x2, i_x3, i_y1, i_y2, i_y3
      REAL (KIND = 8) :: diagCell, n1, n2, n3, ns, pos1_x, pos1_y, pos2_x, pos2_y, pos3_x, pos3_y, pt1, pt2, pt3, &
                         aval, bval, cval, p_pos1, p_pos2, p_pos3, sur2nodeDis, dpdn, p_x1, p_x2
      REAL (KIND = 8), DIMENSION(4,4) :: invA, outA
      REAL (KIND = 8), DIMENSION(3,3) :: inv3A, out3A 
      REAL (KIND = 8), DIMENSION(4) :: ca0
      REAL (KIND = 8), DIMENSION(3) :: ca1
	  INTEGER , intent(in)        ::ibCellCount
	  INTEGER (KIND = 8), intent(in) :: interceptedIndexPtr(:, :), ita, nel2p(:)
	  REAL (KIND = 8), intent(in)   ::xddot, pNormDis(:), xp(:), xcent(:), deltax(:), deltay(:),&
									  cosAlpha(:), cosBeta(:), yp(:)
	  REAL (KIND = 8), intent(inout)   ::p(:,:)
             CHARACTER*70  filename2
                    CHARACTER*70  filename3
      !IF (mod(ita,10000).eq.0) THEN 
	   WRITE(filename2,2) ita
	   WRITE(filename3,3) ita
           !WRITE(filename4,4) ita
	   !WRITE(filename5,5) ita
           !WRITE(filename6,6) ita
	   !WRITE(filename7,7) ita
 2    	   FORMAT('pressure_intercepted/p_int.',i9.9,".dat")
 3    	   FORMAT('pressure_intercepted/p_sur.',i9.9,".dat")
 !4    	   FORMAT('pressure_intercepted/du1dn.',i9.9,".dat")
 !5    	   FORMAT('pressure_intercepted/du2dn.',i9.9,".dat")
 !6    	   FORMAT('pressure_intercepted/dv1dn.',i9.9,".dat")
 !7    	   FORMAT('pressure_intercepted/dv2dn.',i9.9,".dat")
        		open(UNIT = 787, FILE = filename2, STATUS = 'unknown')
        		open(UNIT = 788, FILE = filename3, STATUS = 'unknown')
                        !open(UNIT = 789, FILE = filename4, STATUS = 'unknown')
        		!open(UNIT = 790, FILE = filename5, STATUS = 'unknown')
        		!open(UNIT = 791, FILE = filename6, STATUS = 'unknown')
        		!open(UNIT = 792, FILE = filename7, STATUS = 'unknown')
        !ENDIF
        !print*, 1
      DO k = 1, ibCellCount
         dpdn = -xddot
         i = interceptedIndexPtr(k, 1) 
         j = interceptedIndexPtr(k, 2) 
          sur2nodeDis = pNormDis(k)
         pt1 = 1._rk*dsqrt(deltax(i)**2 + deltay(j)**2)
         !pos1_x = xp(i) + pt1*cosAlpha(nel2p(k))
         !pos1_y = yp(j) + pt1*cosBeta(nel2p(k))
        ! print*,  i, j , xp(i), yp(j)
        ! pause
         !if(dabs(pos1_x-xp(i)).gt.(2.*deltax(i))) pt1 = 1.5_rk*deltax(i)
	  !if(dabs(pos1_y-yp(j)).gt.(2.*deltay(j))) pt1 = 1.5_rk*deltay(j)
	  !normal distance from intercepted cell pressure node - pt1, pt2, pt3
	      pt2 = 2._rk*dsqrt(deltax(i)**2 + deltay(j)**2)    
	      !pt3 = 3_rk*pt1  

	  !coordinates of three points from interceptd cell pressure node
	  pos1_x = xp(i) + pt1*cosAlpha(nel2p(k))
	  pos2_x = xp(i) + pt2*cosAlpha(nel2p(k))
	  !pos3_x = xp(i) + pt3*cosAlpha(nel2p(k))
	  pos1_y = yp(j) + pt1*cosBeta(nel2p(k))
         pos2_y = yp(j) + pt2*cosBeta(nel2p(k))
         !pos3_y = yp(j) + pt3*cosBeta(nel2p(k))
         !Print*, pt1, pt2, pt3, cosAlpha(nel2p(k)), cosBeta(nel2p(k))
         !Print*, pos1_x, pos2_x, pos3_x, xp(i), i, cosAlpha(nel2p(k))
         !Print*, pos1_y, pos2_y, pos3_y, yp(j), j, cosBeta(nel2p(k))
         DO 10 il = i-7, i+7         
            if(pos1_x.ge.xp(il).and.pos1_x.lt.xp(il+1)) i_x1 = il
	     if(pos2_x.ge.xp(il).and.pos2_x.lt.xp(il+1)) i_x2 = il
	     !if(pos3_x.ge.xp(il).and.pos3_x.lt.xp(il+1)) i_x3 = il
 10      CONTINUE        
         DO 20 jl = j-7, j+7       
            if(pos1_y.ge.yp(jl).and.pos1_y.lt.yp(jl+1)) i_y1 = jl
	     if(pos2_y.ge.yp(jl).and.pos2_y.lt.yp(jl+1)) i_y2 = jl
	     !if(pos3_y.ge.yp(jl).and.pos3_y.lt.yp(jl+1)) i_y3 = jl
 20      CONTINUE
         	                           !print*, i, j, i_y1, i_y2, i_y3
         !pressure at pos1
         !interpolation along x  
         p_x1 = p(i_x1, i_y1)   + (p(i_x1+1, i_y1)   - p(i_x1, i_y1))  *(pos1_x - xp(i_x1))/(xp(i_x1+1) - xp(i_x1))
         p_x2 = p(i_x1, i_y1+1) + (p(i_x1+1, i_y1+1) - p(i_x1, i_y1+1))*(pos1_x - xp(i_x1))/(xp(i_x1+1) - xp(i_x1))
         !interpolation along y
         p_pos1 = p_x1 + (p_x2 - p_x1)*(pos1_y - yp(i_y1))/(yp(i_y1+1)-yp(i_y1))
         !pressure at pos2
         !interpolation along x  
         p_x1 = p(i_x2, i_y2)   + (p(i_x2+1, i_y2)   - p(i_x2, i_y2))  *(pos2_x - xp(i_x2))/(xp(i_x2+1) - xp(i_x2))
         p_x2 = p(i_x2, i_y2+1) + (p(i_x2+1, i_y2+1) - p(i_x2, i_y2+1))*(pos2_x - xp(i_x2))/(xp(i_x2+1) - xp(i_x2))
         !interpolation along y
         p_pos2 = p_x1 + (p_x2 - p_x1)*(pos2_y - yp(i_y2))/(yp(i_y2+1)-yp(i_y2))
         !pressure at pos3
         !interpolation along x  
         !p_x1 = p(i_x3, i_y3)   + (p(i_x3+1, i_y3)   - p(i_x3, i_y3))  *(pos3_x - xp(i_x3))/(xp(i_x3+1) - xp(i_x3))
         !p_x2 = p(i_x3, i_y3+1) + (p(i_x3+1, i_y3+1) - p(i_x3, i_y3+1))*(pos3_x - xp(i_x3))/(xp(i_x3+1) - xp(i_x3))
         !interpolation along y
         !p_pos3 = p_x1 + (p_x2 - p_x1)*(pos3_y - yp(i_y3))/(yp(i_y3+1)-yp(i_y3))
         !quadratic interpolation for correct pressure at intercepted cell
                  	 
		
		 
         !IF (sur2nodeDis.GE.0) THEN
             n1 = pt1 + sur2nodeDis
             n2 = pt2 + sur2nodeDis 
         !ELSE
            ! n1 = pt2 + sur2nodeDis
            ! n2 = pt3 + sur2nodeDis
            ! p_pos1 = p_pos2
            ! p_pos2 = p_pos3                        
         !ENDIF
         bval = dpdn
         aval = ((p_pos2-p_pos1) - bval*(n2-n1))/((n2+n1)*(n2-n1))
         cvaL = p_pos2 - aval*n2**2 - bval*n2
         p(i,j) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval 
                  	                   !        print*,k, 1
         !IF (mod(ita,10000).eq.0) THEN
		   WRITE(787,*) xp(i),  p(i,j)
	          WRITE(788,*) xcent(nel2p(k)), cval
         !ENDIF
                       !  print*,k
      ENDDO
      
END SUBROUTINE pressure_forcing_1


SUBROUTINE velocity_forcing_1(ibCellCount, xdot, ydot, interceptedIndexPtr, u2NormDis, &
							deltax, deltay, xu, cosAlpha, nel2p, yu, cosBeta, ut, u1NormDis,& 
							vt, v2NormDis, v1NormDis, xv, yv) 
      
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: k, j, i, ii, il, jl, i_x1, i_x2, i_x3, i_y1, i_y2, i_y3
      REAL (KIND = 8) :: diagCell, n1, n2, n3, ns, pos1_x, pos1_y, pos2_x, pos2_y, pos3_x, pos3_y, pt1, pt2, pt3, &
                         aval, bval, cval, usurf, u_pos1, u_pos2, u_pos3, sur2nodeDis,  u_x1, u_x2, &
                         vsurf, v_pos1, v_pos2, v_pos3, v_x1, v_x2
	  REAL (KIND = 8), intent(in)   ::xdot, ydot, u1NormDis(:), u2NormDis(:) , v1NormDis(:), &
									  v2NormDis(:), cosAlpha(:), cosBeta(:), deltax(:), deltay(:),&
									  xu(:), yu(:), xv(:), yv(:)
	  REAL (KIND = 8), intent(inout)   :: ut(:, :), vt(:, :)
	  INTEGER , intent(in)        ::ibCellCount
	  INTEGER (KIND = 8), intent(in) :: interceptedIndexPtr(:, :),nel2p(:)
	  

      DO k = 1, ibCellCount
         usurf = xdot
         vsurf = ydot
         i = interceptedIndexPtr(k, 1) 
         j = interceptedIndexPtr(k, 2) 

         sur2nodeDis = u2NormDis(k)

          pt1 = 1._rk*dsqrt(deltax(i)**2 + deltay(j)**2)
         
         
         !if(dabs(pos1_x-xu(i+1)).gt.(2.*deltax(i))) pt1 = 1.5_rk*deltax(i)
	  !if(dabs(pos1_y-yu(j)).gt.(2.*deltay(j)))   pt1 = 1.5_rk*deltay(j)
	  
	  !normal distance from intercepted cell pressure node - pt1, pt2, pt3
	 pt2 =2._rk*dsqrt(deltax(i)**2 + deltay(j)**2)
	  pt3 = 3_rk*pt1  
	  !coordinates of three points from interceptd cell pressure node
	  pos1_x = xu(i+1) + pt1*cosAlpha(nel2p(k))
	  pos2_x = xu(i+1) + pt2*cosAlpha(nel2p(k))
	  pos3_x = xu(i+1) + pt3*cosAlpha(nel2p(k))
	  pos1_y = yu(j)   + pt1*cosBeta(nel2p(k))
         pos2_y = yu(j)   + pt2*cosBeta(nel2p(k))
         pos3_y = yu(j)   + pt3*cosBeta(nel2p(k))
         !Print*, pt1, pt2, pt3, cosAlpha(nel2p(k)), cosBeta(nel2p(k))
         !Print*, pos1_x, pos2_x, pos3_x, xp(i), i, cosAlpha(nel2p(k))
         !Print*, pos1_y, pos2_y, pos3_y, yp(j), j, cosBeta(nel2p(k))
         DO 10 il = i-7, i+7         
            if(pos1_x.ge.xu(il).and.pos1_x.lt.xu(il+1)) i_x1 = il
	     if(pos2_x.ge.xu(il).and.pos2_x.lt.xu(il+1)) i_x2 = il
	     if(pos3_x.ge.xu(il).and.pos3_x.lt.xu(il+1)) i_x3 = il
 10      CONTINUE        
         DO 20 jl = j-7, j+7         
            if(pos1_y.ge.yu(jl).and.pos1_y.lt.yu(jl+1)) i_y1 = jl
	     if(pos2_y.ge.yu(jl).and.pos2_y.lt.yu(jl+1)) i_y2 = jl
	     if(pos3_y.ge.yu(jl).and.pos3_y.lt.yu(jl+1)) i_y3 = jl
 20      CONTINUE
         !velocity u right side at pos1
         !interpolation along x  
         u_x1 = ut(i_x1-1, i_y1)   + (ut(i_x1, i_y1)   - ut(i_x1-1, i_y1))  *(pos1_x - xu(i_x1))/(xu(i_x1) - xu(i_x1-1))
         u_x2 = ut(i_x1-1, i_y1+1) + (ut(i_x1, i_y1+1) - ut(i_x1-1, i_y1+1))*(pos1_x - xu(i_x1))/(xu(i_x1) - xu(i_x1-1))
         !interpolation along y
         u_pos1 = u_x1 + (u_x2 - u_x1)*(pos1_y - yu(i_y1))/(yu(i_y1+1)-yu(i_y1))
         !velocity u right side at pos2
         !interpolation along x  
         u_x1 = ut(i_x2-1, i_y2)   + (ut(i_x2, i_y2)   - ut(i_x2-1, i_y2))  *(pos2_x - xu(i_x2))/(xu(i_x2) - xu(i_x2-1))
         u_x2 = ut(i_x2-1, i_y2+1) + (ut(i_x2, i_y2+1) - ut(i_x2-1, i_y2+1))*(pos2_x - xu(i_x2))/(xu(i_x2) - xu(i_x2-1))
         !interpolation along y
         u_pos2 = u_x1 + (u_x2 - u_x1)*(pos2_y - yu(i_y2))/(yu(i_y2+1)-yu(i_y2))
         !velocity u right side at pos3
         !interpolation along x  
         u_x1 = ut(i_x3-1, i_y3)   + (ut(i_x3, i_y3)   - ut(i_x3-1, i_y3))  *(pos3_x - xu(i_x3))/(xu(i_x3) - xu(i_x3-1))
         u_x2 = ut(i_x3-1, i_y3+1) + (ut(i_x3, i_y3+1) - ut(i_x3-1, i_y3+1))*(pos3_x - xu(i_x3))/(xu(i_x3) - xu(i_x3-1))
         !interpolation along y
         u_pos3 = u_x1 + (u_x2 - u_x1)*(pos3_y - yu(i_y3))/(yu(i_y3+1)-yu(i_y3))
         !quadratic interpolation for correct pressure at intercepted cell
		 
         IF (sur2nodeDis.GT.0) THEN
             n1 = pt1 + sur2nodeDis
             n2 = pt2 + sur2nodeDis 
         ELSE
             n1 = pt2 + sur2nodeDis
             n2 = pt3 + sur2nodeDis
             u_pos1 = u_pos2
             u_pos2 = u_pos3                        
         ENDIF
         cval = usurf
         bval = u_pos2-cval - (u_pos1-cval)*(n2**2/n1**2)
         bval = bval/(n2 - (n2**2)/n1)
         avaL = (u_pos1 - cval - bval*n1)/n1**2
         ut(i,j) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval     
         
         
 !******************************************************************   
         
         sur2nodeDis = u1NormDis(k)

        pt1 = 1._rk*dsqrt(deltax(i)**2 + deltay(j)**2)
         !pos1_x = xu(i)   + pt1*cosAlpha(nel2p(k))
         !pos1_y = yu(j)   + pt1*cosBeta(nel2p(k))
         
         !if(dabs(pos1_x-xu(i)).gt.(2.*deltax(i))) pt1 = 1.5_rk*deltax(i)
	  !if(dabs(pos1_y-yu(j)).gt.(2.*deltay(j))) pt1 = 1.5_rk*deltay(j)
	  !normal distance from intercepted cell pressure node - pt1, pt2, pt3
	   pt2 = 2._rk*dsqrt(deltax(i)**2 + deltay(j)**2)   
	  pt3 = 3._rk*pt1  
	  !coordinates of three points from interceptd cell pressure node
	  pos1_x = xu(i) + pt1*cosAlpha(nel2p(k))
	  pos2_x = xu(i) + pt2*cosAlpha(nel2p(k))
	  pos3_x = xu(i) + pt3*cosAlpha(nel2p(k))
	  pos1_y = yu(j) + pt1*cosBeta(nel2p(k))
         pos2_y = yu(j) + pt2*cosBeta(nel2p(k))
         pos3_y = yu(j) + pt3*cosBeta(nel2p(k))
         !Print*, pt1, pt2, pt3, cosAlpha(nel2p(k)), cosBeta(nel2p(k))
         !Print*, pos1_x, pos2_x, pos3_x, xp(i), i, cosAlpha(nel2p(k))
         !Print*, pos1_y, pos2_y, pos3_y, yp(j), j, cosBeta(nel2p(k))
         DO 30 il = i-7, i+7         
            if(pos1_x.ge.xu(il).and.pos1_x.lt.xu(il+1)) i_x1 = il
	     if(pos2_x.ge.xu(il).and.pos2_x.lt.xu(il+1)) i_x2 = il
	     if(pos3_x.ge.xu(il).and.pos3_x.lt.xu(il+1)) i_x3 = il
 30      CONTINUE        
         DO 40 jl = j-7, j+7         
            if(pos1_y.ge.yu(jl).and.pos1_y.lt.yu(jl+1)) i_y1 = jl
	     if(pos2_y.ge.yu(jl).and.pos2_y.lt.yu(jl+1)) i_y2 = jl
	     if(pos3_y.ge.yu(jl).and.pos3_y.lt.yu(jl+1)) i_y3 = jl
 40      CONTINUE
         !velocity u left side at pos1
         !interpolation along x  
         u_x1 = ut(i_x1-1, i_y1)   + (ut(i_x1, i_y1)   - ut(i_x1-1, i_y1))  *(pos1_x - xu(i_x1))/(xu(i_x1) - xu(i_x1-1))
         u_x2 = ut(i_x1-1, i_y1+1) + (ut(i_x1, i_y1+1) - ut(i_x1-1, i_y1+1))*(pos1_x - xu(i_x1))/(xu(i_x1) - xu(i_x1-1))
         !interpolation along y
         u_pos1 = u_x1 + (u_x2 - u_x1)*(pos1_y - yu(i_y1))/(yu(i_y1+1)-yu(i_y1))
         !velocity u left side at pos2
         !interpolation along x  
         u_x1 = ut(i_x2-1, i_y2)   + (ut(i_x2, i_y2)   - ut(i_x2-1, i_y2))  *(pos2_x - xu(i_x2))/(xu(i_x2) - xu(i_x2-1))
         u_x2 = ut(i_x2-1, i_y2+1) + (ut(i_x2, i_y2+1) - ut(i_x2-1, i_y2+1))*(pos2_x - xu(i_x2))/(xu(i_x2) - xu(i_x2-1))
         !interpolation along y
         u_pos2 = u_x1 + (u_x2 - u_x1)*(pos2_y - yu(i_y2))/(yu(i_y2+1)-yu(i_y2))
         !velocity u left side at pos3
         !interpolation along x  
         u_x1 = ut(i_x3-1, i_y3)   + (ut(i_x3, i_y3)   - ut(i_x3-1, i_y3))  *(pos3_x - xu(i_x3))/(xu(i_x3) - xu(i_x3-1))
         u_x2 = ut(i_x3-1, i_y3+1) + (ut(i_x3, i_y3+1) - ut(i_x3-1, i_y3+1))*(pos3_x - xu(i_x3))/(xu(i_x3) - xu(i_x3-1))
         !interpolation along y
         u_pos3 = u_x1 + (u_x2 - u_x1)*(pos3_y - yu(i_y3))/(yu(i_y3+1)-yu(i_y3))
         !quadratic interpolation for correct velocity at intercepted cell
	
         IF (sur2nodeDis.GT.0) THEN
             n1 = pt1 + sur2nodeDis
             n2 = pt2 + sur2nodeDis 
         ELSE
             n1 = pt2 + sur2nodeDis
             n2 = pt3 + sur2nodeDis
             u_pos1 = u_pos2
             u_pos2 = u_pos3                        
         ENDIF
         cval = usurf
         bval = u_pos2-cval - (u_pos1-cval)*(n2**2/n1**2)
         bval = bval/(n2 - (n2**2)/n1)
         avaL = (u_pos1 - cval - bval*n1)/n1**2
         ut(i-1,j) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval 
                  
  !******************************************************************   
         
         sur2nodeDis = v2NormDis(k)

          pt1 = 1._rk*dsqrt(deltax(i)**2 + deltay(j)**2)
         
         !if(dabs(pos1_x-xv(i)).gt.(2.*deltax(i)))     pt1 = 1.5_rk*deltax(i)
	  !if(dabs(pos1_y-yv(j+1)).gt.(2.*deltay(j)))   pt1 = 1.5_rk*deltay(j)
	  !normal distance from intercepted cell pressure node - pt1, pt2, pt3
 pt2 = 2._rk*dsqrt(deltax(i)**2 + deltay(j)**2)
	  pt3 = 3_rk*pt1  
	  !coordinates of three points from interceptd cell pressure node
	  pos1_x = xv(i) + pt1*cosAlpha(nel2p(k))
	  pos2_x = xv(i) + pt2*cosAlpha(nel2p(k))
	  pos3_x = xv(i) + pt3*cosAlpha(nel2p(k))
	  pos1_y = yv(j+1) + pt1*cosBeta(nel2p(k))
         pos2_y = yv(j+1) + pt2*cosBeta(nel2p(k))
         pos3_y = yv(j+1) + pt3*cosBeta(nel2p(k))
         !Print*, pt1, pt2, pt3, cosAlpha(nel2p(k)), cosBeta(nel2p(k))
         !Print*, pos1_x, pos2_x, pos3_x, xp(i), i, cosAlpha(nel2p(k))
         !Print*, pos1_y, pos2_y, pos3_y, yp(j), j, cosBeta(nel2p(k))
         DO 50 il = i-7, i+7         
            if(pos1_x.ge.xv(il).and.pos1_x.lt.xv(il+1)) i_x1 = il
	     if(pos2_x.ge.xv(il).and.pos2_x.lt.xv(il+1)) i_x2 = il
	     if(pos3_x.ge.xv(il).and.pos3_x.lt.xv(il+1)) i_x3 = il
 50      CONTINUE        
         DO 60 jl = j-7, j+7         
            if(pos1_y.ge.yv(jl).and.pos1_y.lt.yv(jl+1)) i_y1 = jl
	     if(pos2_y.ge.yv(jl).and.pos2_y.lt.yv(jl+1)) i_y2 = jl
	     if(pos3_y.ge.yv(jl).and.pos3_y.lt.yv(jl+1)) i_y3 = jl
 60      CONTINUE
         !interpolation along x  
         v_x1 = vt(i_x1, i_y1-1)   + (vt(i_x1, i_y1)   - vt(i_x1, i_y1-1))  *(pos1_y - yv(i_y1))/(yv(i_y1+1) - yv(i_y1))
         v_x2 = vt(i_x1+1, i_y1-1) + (vt(i_x1+1, i_y1) - vt(i_x1+1, i_y1-1))*(pos1_y - yv(i_y1))/(yv(i_y1+1) - yv(i_y1))
         !interpolation along y
         v_pos1 = v_x1 + (v_x2 - v_x1)*(pos1_x - xv(i_x1))/(xv(i_x1+1)-xv(i_x1))
         !velocity u right side at pos2
         !interpolation along x  
         v_x1 = vt(i_x2, i_y2-1)   + (vt(i_x2, i_y2)   - vt(i_x2, i_y2-1))  *(pos2_y - yv(i_y2))/(yv(i_y2+1) - yv(i_y2))
         v_x2 = vt(i_x2+1, i_y2-1) + (vt(i_x2+1, i_y2) - vt(i_x2+1, i_y2-1))*(pos2_y - yv(i_y2))/(yv(i_y2+1) - yv(i_y2))
         !interpolation along y
         v_pos2 = v_x1 + (v_x2 - v_x1)*(pos2_x - xv(i_x2))/(xv(i_x2+1)-xv(i_x2))
         !velocity u right side at pos3
         !interpolation along x  
         v_x1 = vt(i_x3, i_y3-1)   + (vt(i_x3, i_y3)   - vt(i_x3, i_y3-1))  *(pos3_y - yv(i_y3))/(yv(i_y3+1) - yv(i_y3))
         v_x2 = vt(i_x3+1, i_y3-1) + (vt(i_x3+1, i_y3) - vt(i_x3+1, i_y3-1))*(pos3_y - yv(i_y3))/(yv(i_y3+1) - yv(i_y3))
         !interpolation along y
         v_pos3 = v_x1 + (v_x2 - v_x1)*(pos3_x - xv(i_x3))/(xv(i_x3+1)-xv(i_x3))
         !quadratic interpolation for correct pressure at intercepted cell
         !quadratic interpolation for correct pressure at intercepted cell
		 
         IF (sur2nodeDis.GT.0) THEN
             n1 = pt1 + sur2nodeDis
             n2 = pt2 + sur2nodeDis 
         ELSE
             n1 = pt2 + sur2nodeDis
             n2 = pt3 + sur2nodeDis
             v_pos1 = v_pos2
             v_pos2 = v_pos3                        
         ENDIF
         cval = vsurf
         bval = v_pos2-cval - (v_pos1-cval)*(n2**2/n1**2)
         bval = bval/(n2 - (n2**2)/n1)
         avaL = (v_pos1 - cval - bval*n1)/n1**2
         vt(i,j) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval 

 !******************************************************************   
         
         sur2nodeDis = v1NormDis(k)

          pt1 = 1._rk*dsqrt(deltax(i)**2 + deltay(j)**2)
         !pos1_x = xv(i)   + pt1*cosAlpha(nel2p(k))
        ! pos1_y = yv(j) + pt1*cosBeta(nel2p(k))
         !if(dabs(pos1_x-xv(i)).gt.(2.*deltax(i)))     pt1 = 1.5_rk*deltax(i)
	  !if(dabs(pos1_y-yv(j)).gt.(2.*deltay(j)))     pt1 = 1.5_rk*deltay(j)
	  !normal distance from intercepted cell pressure node - pt1, pt2, pt3
 pt2 = 2._rk*dsqrt(deltax(i)**2 + deltay(j)**2)
	  pt3 = 3._rk*pt1  
	  !coordinates of three points from interceptd cell pressure node
	  pos1_x = xv(i) + pt1*cosAlpha(nel2p(k))
	  pos2_x = xv(i) + pt2*cosAlpha(nel2p(k))
	  pos3_x = xv(i) + pt3*cosAlpha(nel2p(k))
	  pos1_y = yv(j) + pt1*cosBeta(nel2p(k))
         pos2_y = yv(j) + pt2*cosBeta(nel2p(k))
         pos3_y = yv(j) + pt3*cosBeta(nel2p(k))
         !Print*, pt1, pt2, pt3, cosAlpha(nel2p(k)), cosBeta(nel2p(k))
         !Print*, pos1_x, pos2_x, pos3_x, xp(i), i, cosAlpha(nel2p(k))
         !Print*, pos1_y, pos2_y, pos3_y, yp(j), j, cosBeta(nel2p(k))
         DO 70 il = i-7, i+7         
            if(pos1_x.ge.xv(il).and.pos1_x.lt.xv(il+1)) i_x1 = il
	     if(pos2_x.ge.xv(il).and.pos2_x.lt.xv(il+1)) i_x2 = il
	     if(pos3_x.ge.xv(il).and.pos3_x.lt.xv(il+1)) i_x3 = il
 70      CONTINUE        
         DO 80 jl = j-7, j+7         
            if(pos1_y.ge.yv(jl).and.pos1_y.lt.yv(jl+1)) i_y1 = jl
	     if(pos2_y.ge.yv(jl).and.pos2_y.lt.yv(jl+1)) i_y2 = jl
	     if(pos3_y.ge.yv(jl).and.pos3_y.lt.yv(jl+1)) i_y3 = jl
 80      CONTINUE
         !velocity u right side at pos1
         !interpolation along x  
         v_x1 = vt(i_x1, i_y1-1)   + (vt(i_x1, i_y1)   - vt(i_x1, i_y1-1))  *(pos1_y - yv(i_y1))/(yv(i_y1+1) - yv(i_y1))
         v_x2 = vt(i_x1+1, i_y1-1) + (vt(i_x1+1, i_y1) - vt(i_x1+1, i_y1-1))*(pos1_y - yv(i_y1))/(yv(i_y1+1) - yv(i_y1))
         !interpolation along y
         v_pos1 = v_x1 + (v_x2 - v_x1)*(pos1_x - xv(i_x1))/(xv(i_x1+1)-xv(i_x1))
         !velocity u right side at pos2
         !interpolation along x  
         v_x1 = vt(i_x2, i_y2-1)   + (vt(i_x2, i_y2)   - vt(i_x2, i_y2-1))  *(pos2_y - yv(i_y2))/(yv(i_y2+1) - yv(i_y2))
         v_x2 = vt(i_x2+1, i_y2-1) + (vt(i_x2+1, i_y2) - vt(i_x2+1, i_y2-1))*(pos2_y - yv(i_y2))/(yv(i_y2+1) - yv(i_y2))
         !interpolation along y
         v_pos2 = v_x1 + (v_x2 - v_x1)*(pos2_x - xv(i_x2))/(xv(i_x2+1)-xv(i_x2))
         !velocity u right side at pos3
         !interpolation along x  
         v_x1 = vt(i_x3, i_y3-1)   + (vt(i_x3, i_y3)   - vt(i_x3, i_y3-1))  *(pos3_y - yv(i_y3))/(yv(i_y3+1) - yv(i_y3))
         v_x2 = vt(i_x3+1, i_y3-1) + (vt(i_x3+1, i_y3) - vt(i_x3+1, i_y3-1))*(pos3_y - yv(i_y3))/(yv(i_y3+1) - yv(i_y3))
         !interpolation along y
         v_pos3 = v_x1 + (v_x2 - v_x1)*(pos3_x - xv(i_x3))/(xv(i_x3+1)-xv(i_x3))
         !quadratic interpolation for correct pressure at intercepted cell
	
         IF (sur2nodeDis.GT.0) THEN
             n1 = pt1 + sur2nodeDis
             n2 = pt2 + sur2nodeDis 
         ELSE
             n1 = pt2 + sur2nodeDis
             n2 = pt3 + sur2nodeDis
             v_pos1 = v_pos2
             v_pos2 = v_pos3                        
         ENDIF
         cval = vsurf
         bval = v_pos2-cval - (v_pos1-cval)*(n2**2/n1**2)
         bval = bval/(n2 - (n2**2)/n1)
         avaL = (v_pos1 - cval - bval*n1)/n1**2
         vt(i,j-1) = aval*sur2nodeDis**2 + bval*sur2nodeDis + cval
         !goto 100


        !print*,'divergence', ( ut(i,j) - ut(i-1,j) )/deltax(i) +  ( vt(i,j) - vt(i,j-1) )/deltay(j), i, j, ut(i,j)
          
 100 continue         
      ENDDO
      
END SUBROUTINE velocity_forcing_1

