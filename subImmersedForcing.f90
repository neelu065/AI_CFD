SUBROUTINE pressureForcing
      USE global
      IMPLICIT NONE
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: k, j, i, ii, il, jl, i_x1, i_x2, i_x3, i_y1, i_y2, i_y3
      REAL (KIND = 8) :: diagCell, n1, n2, n3, ns, pos1_x, pos1_y, pos2_x, pos2_y, pos3_x, pos3_y, pt1, pt2, pt3, &
                         xp1, xp2, xp3, xp4, yp1, yp2, yp3, yp4, p_pos1, p_pos2, p_pos3, sur2nodeDis, dpdn
      REAL (KIND = 8), DIMENSION(4,4) :: invA, outA
      REAL (KIND = 8), DIMENSION(3,3) :: inv3A, out3A 
      REAL (KIND = 8), DIMENSION(4) :: ca0
      REAL (KIND = 8), DIMENSION(3) :: ca1
      DO k = 1, ibCellCount
         invA = 0._rk
         outA = 0._rk
         inv3A = 0._rk
         out3A = 0._rk 
         dpdn = 0._rk
         i = interceptedIndexPtr(k, 1) 
         j = interceptedIndexPtr(k, 2) 
         sur2nodeDis = pNormDis(k)
         pt1 = dsqrt(deltax(i)**2 + deltay(j)**2)
         pos1_x = xp(i) + pt1*cosAlpha(nel2p(k))
         pos1_y = yp(j) + pt1*cosBeta(nel2p(k))
         if(dabs(pos1_x-xp(i)).gt.(2.*deltax(i))) pt1 = 1.5_rk*deltax(i)
	  if(dabs(pos1_y-yp(j)).gt.(2.*deltay(j))) pt1 = 1.5_rk*deltay(j)
	  !normal distance from intercepted cell pressure node - pt1, pt2, pt3
	  pt2 = 2_rk*pt1    
	  pt3 = 3_rk*pt1  
	  !coordinates of three points from interceptd cell pressure node
	  pos1_x = xp(i) + pt1*cosAlpha(nel2p(k))
	  pos2_x = xp(i) + pt2*cosAlpha(nel2p(k))
	  pos3_x = xp(i) + pt3*cosAlpha(nel2p(k))
	  pos1_y = yp(j) + pt1*cosBeta(nel2p(k))
         pos2_y = yp(j) + pt2*cosBeta(nel2p(k))
         pos3_y = yp(j) + pt3*cosBeta(nel2p(k))
         !Print*, pt1, pt2, pt3, cosAlpha(nel2p(k)), cosBeta(nel2p(k))
         !Print*, pos1_x, pos2_x, pos3_x, xp(i), i, cosAlpha(nel2p(k))
         !Print*, pos1_y, pos2_y, pos3_y, yp(j), j, cosBeta(nel2p(k))
         DO 10 il = i-7, i+7         
            if(pos1_x.ge.xp(il).and.pos1_x.lt.xp(il+1)) i_x1 = il
	     if(pos2_x.ge.xp(il).and.pos2_x.lt.xp(il+1)) i_x2 = il
	     if(pos3_x.ge.xp(il).and.pos3_x.lt.xp(il+1)) i_x3 = il
 10      CONTINUE        
         DO 20 jl = j-7, j+7         
            if(pos1_y.ge.yp(jl).and.pos1_y.lt.yp(jl+1)) i_y1 = jl
	     if(pos2_y.ge.yp(jl).and.pos2_y.lt.yp(jl+1)) i_y2 = jl
	     if(pos3_y.ge.yp(jl).and.pos3_y.lt.yp(jl+1)) i_y3 = jl
 20      CONTINUE  
        ! Print*, i_x1, i_x2, i_x3, 'x'
        ! Print*, i_y1, i_y2, i_y3, 'y'
         xp1 = xp(i_x1) 
         yp1 = yp(i_y1) 
         
         xp2 = xp(i_x1)
         yp2 = yp(i_y1+1)  
                
         xp3 = xp(i_x1+1)
         yp3 = yp(i_y1+1)    
              
         xp4 = xp(i_x1+1)
         yp4 = yp(i_y1)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*p(i_x1, i_y1) + outA(i,2)*p(i_x1, i_y1+1) + outA(i,3)*p(i_x1+1, i_y1+1) + outA(i,4)*p(i_x1+1, i_y1)
         ENDDO
         
         !*******************************
         p_pos1 = ca0(1) + ca0(2)*pos1_x + ca0(3)*pos1_y + ca0(4)*pos1_x*pos1_y  
         !********************************
         
         xp1 = xp(i_x2) 
         yp1 = yp(i_y2) 
         
         xp2 = xp(i_x2)
         yp2 = yp(i_y2+1)  
                
         xp3 = xp(i_x2+1)
         yp3 = yp(i_y2+1)    
              
         xp4 = xp(i_x2+1)
         yp4 = yp(i_y2)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*p(i_x2, i_y2) + outA(i,2)*p(i_x2, i_y2+1) + outA(i,3)*p(i_x2+1, i_y2+1) + outA(i,4)*p(i_x2+1, i_y2)
         ENDDO
         
         !**************************
         p_pos2 = ca0(1) + ca0(2)*pos2_x + ca0(3)*pos2_y + ca0(4)*pos2_x*pos2_y         
         !***************************
         
         xp1 = xp(i_x3) 
         yp1 = yp(i_y3) 
         
         xp2 = xp(i_x3)
         yp2 = yp(i_y3+1)  
                
         xp3 = xp(i_x3+1)
         yp3 = yp(i_y3+1)    
              
         xp4 = xp(i_x3+1)
         yp4 = yp(i_y3)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*p(i_x3, i_y3) + outA(i,2)*p(i_x3, i_y3+1) + outA(i,3)*p(i_x3+1, i_y3+1) + outA(i,4)*p(i_x3+1, i_y3)
         ENDDO
         
         !*************************************
         p_pos3 = ca0(1) + ca0(2)*pos3_x + ca0(3)*pos3_y + ca0(4)*pos3_x*pos3_y    
         !*************************************
         
         
         !quadratic interpolation for correct pressure at intercepted cell
         IF (sur2nodeDis.GE.0._rk) THEN
             n1 = pt1 + sur2nodeDis
             n2 = pt2 + sur2nodeDis 
         ELSE
             n1 = pt2 - sur2nodeDis
             n2 = pt3 - sur2nodeDis
             p_pos1 = p_pos2
             p_pos2 = p_pos3                        
         ENDIF
         ns = 0._rk
         inv3A(1,1) =  n1*n1
         inv3A(1,2) =  n1
         inv3A(1,3) =  1._rk
         inv3A(2,1) =  n2*n2
         inv3A(2,2) =  n2
         inv3A(2,3) =  1._rk 
         inv3A(3,1) =  2._rk*ns
         inv3A(3,2) =  1._rk
         inv3A(3,3) =  0._rk    
         Call matinv3(inv3A, out3A)
         DO ii = 1, 3
            ca1(ii) = out3A(i,1)*p_pos1 + out3A(i,2)*p_pos2 + out3A(i,3)*dpdn
         ENDDO
         
         !*************************************************************
         p(i,j) = ca1(1)*sur2nodeDis**2 + ca1(2)*sur2nodeDis + ca1(3)
         !*************************************************************
                  !print*, ca1(1), ca1(2), ca1(3), 'p'
      ENDDO
      
END SUBROUTINE pressureForcing

!*********************************************************************
SUBROUTINE velocityForcing
      USE global
      IMPLICIT NONE
      INTEGER, PARAMETER :: rk = selected_real_kind(8)
      INTEGER :: k, j, i, ii, il, jl, i_x1, i_x2, i_x3, i_y1, i_y2, i_y3
      REAL (KIND = 8) :: diagCell, n1, n2, n3, ns, pos1_x, pos1_y, pos2_x, pos2_y, pos3_x, pos3_y, pt1, pt2, pt3, &
                         xp1, xp2, xp3, xp4, yp1, yp2, yp3, yp4, u_pos1, u_pos2, u_pos3, v_pos1, v_pos2, v_pos3, &
                         sur2nodeDis, usurf, vsurf
      REAL (KIND = 8), DIMENSION(4,4) :: invA, outA
      REAL (KIND = 8), DIMENSION(3,3) :: inv3A, out3A 
      REAL (KIND = 8), DIMENSION(4) :: ca0
      REAL (KIND = 8), DIMENSION(3) :: ca1

      DO k = 1, ibCellCount
      
         i = interceptedIndexPtr(k, 1) 
         j = interceptedIndexPtr(k, 2) 
         !u velocity right
         invA = 0._rk
         outA = 0._rk
         inv3A = 0._rk
         out3A = 0._rk 
         usurf = 0._rk
         !print*, i,j         
         !right side u velocity
         sur2nodeDis = u2NormDis(k)
         pt1 = dsqrt(deltax(i)**2 + deltay(j)**2)
         pos1_x = xu(i+1) + pt1*cosAlpha(nel2p(k))
         pos1_y = yu(j)   + pt1*cosBeta(nel2p(k))
         if(dabs(pos1_x-xu(i+1)).gt.(2.*deltax(i))) pt1 = 1.5_rk*deltax(i)
	  if(dabs(pos1_y-yu(j)).gt.(2.*deltay(j)))   pt1 = 1.5_rk*deltay(j)
	  !normal distance from intercepted cell pressure node - pt1, pt2, pt3
	  pt2 = 2_rk*pt1    
	  pt3 = 3_rk*pt1  
	  !coordinates of three points from interceptd cell pressure node
	  pos1_x = xu(i+1) + pt1*cosAlpha(nel2p(k))
	  pos2_x = xu(i+1) + pt2*cosAlpha(nel2p(k))
	  pos3_x = xu(i+1) + pt3*cosAlpha(nel2p(k))
	  pos1_y = yu(j)   + pt1*cosBeta(nel2p(k))
         pos2_y = yu(j)   + pt2*cosBeta(nel2p(k))
         pos3_y = yu(j)   + pt3*cosBeta(nel2p(k))

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
 
         !Print*, i_x1, i_x2, i_x3, 'x'
         !Print*, i_y1, i_y2, i_y3, 'y'
         !if (ita.gt.0) pause
         xp1 = xu(i_x1) 
         yp1 = yu(i_y1) 
         
         xp2 = xu(i_x1)
         yp2 = yu(i_y1+1)  
                
         xp3 = xu(i_x1+1)
         yp3 = yu(i_y1+1)    
              
         xp4 = xu(i_x1+1)
         yp4 = yu(i_y1)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*ut(i_x1-1, i_y1) + outA(i,2)*ut(i_x1-1, i_y1+1) + &
                      outA(i,3)*ut(i_x1, i_y1+1) + outA(i,4)*ut(i_x1, i_y1)
         ENDDO
          !*******************************
         u_pos1 = ca0(1) + ca0(2)*pos1_x + ca0(3)*pos1_y + ca0(4)*pos1_x*pos1_y  
         !********************************
         
         xp1 = xu(i_x2) 
         yp1 = yu(i_y2) 
         
         xp2 = xu(i_x2)
         yp2 = yu(i_y2+1)  
                
         xp3 = xu(i_x2+1)
         yp3 = yu(i_y2+1)    
              
         xp4 = xu(i_x2+1)
         yp4 = yu(i_y2)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*ut(i_x2-1, i_y2) + outA(i,2)*ut(i_x2-1, i_y2+1) + &
                      outA(i,3)*ut(i_x2, i_y2+1) + outA(i,4)*ut(i_x2, i_y2)
         ENDDO
         
         !**************************
         u_pos2 = ca0(1) + ca0(2)*pos2_x + ca0(3)*pos2_y + ca0(4)*pos2_x*pos2_y         
         !***************************
         
         xp1 = xu(i_x3) 
         yp1 = yu(i_y3) 
         
         xp2 = xu(i_x3)
         yp2 = yu(i_y3+1)  
                
         xp3 = xu(i_x3+1)
         yp3 = yu(i_y3+1)    
              
         xp4 = xu(i_x3+1)
         yp4 = yu(i_y3)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*ut(i_x3-1, i_y3) + outA(i,2)*ut(i_x3-1, i_y3+1) + &
                      outA(i,3)*ut(i_x3, i_y3+1) + outA(i,4)*ut(i_x3, i_y3)
         ENDDO
         
         !*************************************
         u_pos3 = ca0(1) + ca0(2)*pos3_x + ca0(3)*pos3_y + ca0(4)*pos3_x*pos3_y    
         !*************************************
         !print*, u_pos1, u_pos2, u_pos3
        ! print*, 
         !if (ita.gt.1) pause
         !quadratic interpolation for correct pressure at intercepted cell
         IF (sur2nodeDis.GE.0._rk) THEN
             n1 = pt1 + sur2nodeDis
             n2 = pt2 + sur2nodeDis 
         ELSE
             n1 = pt2 - sur2nodeDis
             n2 = pt3 - sur2nodeDis
             u_pos1 = u_pos2
             u_pos2 = u_pos3                        
         ENDIF
         ns = 0._rk
         inv3A(1,1) =  n1*n1
         inv3A(1,2) =  n1
         inv3A(1,3) =  1._rk
         inv3A(2,1) =  n2*n2
         inv3A(2,2) =  n2
         inv3A(2,3) =  1._rk 
         inv3A(3,1) =  ns*ns
         inv3A(3,2) =  ns
         inv3A(3,3) =  1._rk    
         Call matinv3(inv3A, out3A)
         DO ii = 1, 3
            ca1(ii) = out3A(i,1)*u_pos1 + out3A(i,2)*u_pos2 + out3A(i,3)*usurf
         ENDDO
         
         !*****************************
         ut(i,j) = ca1(1)*sur2nodeDis**2 + ca1(2)*sur2nodeDis + ca1(3)
         !*****************************
         
         !u velocity left
         invA = 0._rk
         outA = 0._rk
         inv3A = 0._rk
         out3A = 0._rk 
         usurf = 0._rk
                  
         !right side u velocity
         sur2nodeDis = u1NormDis(k)
         pt1 = dsqrt(deltax(i)**2 + deltay(j)**2)
         pos1_x = xu(i)   + pt1*cosAlpha(nel2p(k))
         pos1_y = yu(j)   + pt1*cosBeta(nel2p(k))
         if(dabs(pos1_x-xu(i)).gt.(2.*deltax(i)))   pt1 = 1.5_rk*deltax(i)
	  if(dabs(pos1_y-yu(j)).gt.(2.*deltay(j)))   pt1 = 1.5_rk*deltay(j)
	  !normal distance from intercepted cell pressure node - pt1, pt2, pt3
	  pt2 = 2_rk*pt1    
	  pt3 = 3_rk*pt1  
	  !coordinates of three points from interceptd cell pressure node
	  pos1_x = xu(i)   + pt1*cosAlpha(nel2p(k))
	  pos2_x = xu(i)   + pt2*cosAlpha(nel2p(k))
	  pos3_x = xu(i)   + pt3*cosAlpha(nel2p(k))
	  pos1_y = yu(j)   + pt1*cosBeta(nel2p(k))
         pos2_y = yu(j)   + pt2*cosBeta(nel2p(k))
         pos3_y = yu(j)   + pt3*cosBeta(nel2p(k))
         
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
        ! Print*, i_x1, i_x2, i_x3, 'x'
        ! Print*, i_y1, i_y2, i_y3, 'y'
         xp1 = xu(i_x1) 
         yp1 = yu(i_y1) 
         
         xp2 = xu(i_x1)
         yp2 = yu(i_y1+1)  
                
         xp3 = xu(i_x1+1)
         yp3 = yu(i_y1+1)    
              
         xp4 = xu(i_x1+1)
         yp4 = yu(i_y1)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*ut(i_x1-1, i_y1) + outA(i,2)*ut(i_x1-1, i_y1+1) + &
                      outA(i,3)*ut(i_x1, i_y1+1) + outA(i,4)*ut(i_x1, i_y1)
         ENDDO
          !*******************************
         u_pos1 = ca0(1) + ca0(2)*pos1_x + ca0(3)*pos1_y + ca0(4)*pos1_x*pos1_y  
         !********************************
         
         xp1 = xu(i_x2) 
         yp1 = yu(i_y2) 
         
         xp2 = xu(i_x2)
         yp2 = yu(i_y2+1)  
                
         xp3 = xu(i_x2+1)
         yp3 = yu(i_y2+1)    
              
         xp4 = xu(i_x2+1)
         yp4 = yu(i_y2)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*ut(i_x2-1, i_y2) + outA(i,2)*ut(i_x2-1, i_y2+1) + &
                      outA(i,3)*ut(i_x2, i_y2+1) + outA(i,4)*ut(i_x2, i_y2)
         ENDDO
         
         !**************************
         u_pos2 = ca0(1) + ca0(2)*pos2_x + ca0(3)*pos2_y + ca0(4)*pos2_x*pos2_y         
         !***************************
         
         xp1 = xu(i_x3) 
         yp1 = yu(i_y3) 
         
         xp2 = xu(i_x3)
         yp2 = yu(i_y3+1)  
                
         xp3 = xu(i_x3+1)
         yp3 = yu(i_y3+1)    
              
         xp4 = xu(i_x3+1)
         yp4 = yu(i_y3)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*ut(i_x3-1, i_y3) + outA(i,2)*ut(i_x3-1, i_y3+1) + &
                      outA(i,3)*ut(i_x3, i_y3+1) + outA(i,4)*ut(i_x3, i_y3)
         ENDDO
         
         !*************************************
         u_pos3 = ca0(1) + ca0(2)*pos3_x + ca0(3)*pos3_y + ca0(4)*pos3_x*pos3_y    
         !*************************************
         
         
         !quadratic interpolation for correct pressure at intercepted cell
         IF (sur2nodeDis.GE.0._rk) THEN
             n1 = pt1 + sur2nodeDis
             n2 = pt2 + sur2nodeDis 
         ELSE
             n1 = pt2 - sur2nodeDis
             n2 = pt3 - sur2nodeDis
             u_pos1 = u_pos2
             u_pos2 = u_pos3                        
         ENDIF
         ns = 0._rk
         inv3A(1,1) =  n1*n1
         inv3A(1,2) =  n1
         inv3A(1,3) =  1._rk
         inv3A(2,1) =  n2*n2
         inv3A(2,2) =  n2
         inv3A(2,3) =  1._rk 
         inv3A(3,1) =  ns*ns
         inv3A(3,2) =  ns
         inv3A(3,3) =  1._rk    
         Call matinv3(inv3A, out3A)
         DO ii = 1, 3
            !print*, inv3A(ii, :)
            ca1(ii) = out3A(i,1)*u_pos1 + out3A(i,2)*u_pos2 + out3A(i,3)*usurf
                       ! print*, out3A(ii, :), 'o'
         ENDDO
         !print*, ca1(1), ca1(2), ca1(3), 'u'
         
         !**************************
         ut(i-1,j) = ca1(1)*sur2nodeDis**2 + ca1(2)*sur2nodeDis + ca1(3)
         !***************************
         goto 1
         
         !v velocity top
         invA = 0._rk
         outA = 0._rk
         inv3A = 0._rk
         out3A = 0._rk 
         vsurf = 0._rk
                  
         !v velocity top
         sur2nodeDis = v2NormDis(k)
         pt1 = dsqrt(deltax(i)**2 + deltay(j)**2)
         pos1_x = xv(i)     + pt1*cosAlpha(nel2p(k))
         pos1_y = yv(j+1)   + pt1*cosBeta(nel2p(k))
         if(dabs(pos1_x-xv(i)).gt.(2.*deltax(i)))     pt1 = 1.5_rk*deltax(i)
	  if(dabs(pos1_y-yv(j+1)).gt.(2.*deltay(j)))   pt1 = 1.5_rk*deltay(j)
	  !normal distance from intercepted cell pressure node - pt1, pt2, pt3
	  pt2 = 2_rk*pt1    
	  pt3 = 3_rk*pt1  
	  !coordinates of three points from interceptd cell pressure node
	  pos1_x = xv(i)     + pt1*cosAlpha(nel2p(k))
	  pos2_x = xv(i)     + pt2*cosAlpha(nel2p(k))
	  pos3_x = xv(i)     + pt3*cosAlpha(nel2p(k))
	  pos1_y = yv(j+1)   + pt1*cosBeta(nel2p(k))
         pos2_y = yv(j+1)   + pt2*cosBeta(nel2p(k))
         pos3_y = yv(j+1)   + pt3*cosBeta(nel2p(k))
         
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
        ! Print*, i_x1, i_x2, i_x3, 'x'
        ! Print*, i_y1, i_y2, i_y3, 'y'
         xp1 = xv(i_x1) 
         yp1 = yv(i_y1) 
         
         xp2 = xv(i_x1)
         yp2 = yv(i_y1+1)  
                
         xp3 = xv(i_x1+1)
         yp3 = yv(i_y1+1)    
              
         xp4 = xv(i_x1+1)
         yp4 = yv(i_y1)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*vt(i_x1, i_y1-1) + outA(i,2)*vt(i_x1, i_y1) + &
                      outA(i,3)*vt(i_x1+1, i_y1) + outA(i,4)*vt(i_x1+1, i_y1-1)
         ENDDO
          !*******************************
         v_pos1 = ca0(1) + ca0(2)*pos1_x + ca0(3)*pos1_y + ca0(4)*pos1_x*pos1_y  
         !********************************
         
         xp1 = xv(i_x2) 
         yp1 = yv(i_y2) 
         
         xp2 = xv(i_x2)
         yp2 = yv(i_y2+1)  
                
         xp3 = xv(i_x2+1)
         yp3 = yv(i_y2+1)    
              
         xp4 = xv(i_x2+1)
         yp4 = yv(i_y2)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*vt(i_x2, i_y2-1) + outA(i,2)*vt(i_x2, i_y2) + &
                      outA(i,3)*vt(i_x2+1, i_y2) + outA(i,4)*vt(i_x2+1, i_y2-1)
         ENDDO
         
         !**************************
         v_pos2 = ca0(1) + ca0(2)*pos2_x + ca0(3)*pos2_y + ca0(4)*pos2_x*pos2_y         
         !***************************
         
         xp1 = xv(i_x3) 
         yp1 = yv(i_y3) 
         
         xp2 = xv(i_x3)
         yp2 = yv(i_y3+1)  
                
         xp3 = xv(i_x3+1)
         yp3 = yv(i_y3+1)    
              
         xp4 = xv(i_x3+1)
         yp4 = yv(i_y3)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*vt(i_x3, i_y3-1) + outA(i,2)*vt(i_x3, i_y3) + &
                      outA(i,3)*vt(i_x3+1, i_y3) + outA(i,4)*vt(i_x3+1, i_y3-1)
         ENDDO
         
         !*************************************
         v_pos3 = ca0(1) + ca0(2)*pos3_x + ca0(3)*pos3_y + ca0(4)*pos3_x*pos3_y    
         !*************************************
         
         
         !quadratic interpolation for correct pressure at intercepted cell
         IF (sur2nodeDis.GE.0._rk) THEN
             n1 = pt1 + sur2nodeDis
             n2 = pt2 + sur2nodeDis 
         ELSE
             n1 = pt2 - sur2nodeDis
             n2 = pt3 - sur2nodeDis
             v_pos1 = v_pos2
             v_pos2 = v_pos3                        
         ENDIF
         ns = 0._rk
         inv3A(1,1) =  n1*n1
         inv3A(1,2) =  n1
         inv3A(1,3) =  1._rk
         inv3A(2,1) =  n2*n2
         inv3A(2,2) =  n2
         inv3A(2,3) =  1._rk 
         inv3A(3,1) =  ns*ns
         inv3A(3,2) =  ns
         inv3A(3,3) =  1._rk    
         Call matinv3(inv3A, out3A)
         DO ii = 1, 3
            !print*, inv3A(ii, :)
            ca1(ii) = out3A(i,1)*v_pos1 + out3A(i,2)*v_pos2 + out3A(i,3)*vsurf
                        !print*, out3A(ii, :), 'o'
         ENDDO
         !print*, ca1(1), ca1(2), ca1(3), 'u'
         vt(i,j) = ca1(1)*sur2nodeDis**2 + ca1(2)*sur2nodeDis + ca1(3)
         
                  !v velocity top
         invA = 0._rk
         outA = 0._rk
         inv3A = 0._rk
         out3A = 0._rk 
         vsurf = 0._rk
                  
         !v velocity top
         sur2nodeDis = v1NormDis(k)
         pt1 = dsqrt(deltax(i)**2 + deltay(j)**2)
         pos1_x = xv(i)   + pt1*cosAlpha(nel2p(k))
         pos1_y = yv(j)   + pt1*cosBeta(nel2p(k))
         if(dabs(pos1_x-xv(i)).gt.(2.*deltax(i)))   pt1 = 1.5_rk*deltax(i)
	  if(dabs(pos1_y-yv(j)).gt.(2.*deltay(j)))   pt1 = 1.5_rk*deltay(j)
	  !normal distance from intercepted cell pressure node - pt1, pt2, pt3
	  pt2 = 2_rk*pt1    
	  pt3 = 3_rk*pt1  
	  !coordinates of three points from interceptd cell pressure node
	  pos1_x = xv(i)   + pt1*cosAlpha(nel2p(k))
	  pos2_x = xv(i)   + pt2*cosAlpha(nel2p(k))
	  pos3_x = xv(i)   + pt3*cosAlpha(nel2p(k))
	  pos1_y = yv(j)   + pt1*cosBeta(nel2p(k))
         pos2_y = yv(j)   + pt2*cosBeta(nel2p(k))
         pos3_y = yv(j)   + pt3*cosBeta(nel2p(k))
         
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
        ! Print*, i_x1, i_x2, i_x3, 'x'
        ! Print*, i_y1, i_y2, i_y3, 'y'
         xp1 = xv(i_x1) 
         yp1 = yv(i_y1) 
         
         xp2 = xv(i_x1)
         yp2 = yv(i_y1+1)  
                
         xp3 = xv(i_x1+1)
         yp3 = yv(i_y1+1)    
              
         xp4 = xv(i_x1+1)
         yp4 = yv(i_y1)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*vt(i_x1, i_y1-1) + outA(i,2)*vt(i_x1, i_y1) + &
                      outA(i,3)*vt(i_x1+1, i_y1) + outA(i,4)*vt(i_x1+1, i_y1-1)
         ENDDO
          !*******************************
         v_pos1 = ca0(1) + ca0(2)*pos1_x + ca0(3)*pos1_y + ca0(4)*pos1_x*pos1_y  
         !********************************
         
         xp1 = xv(i_x2) 
         yp1 = yv(i_y2) 
         
         xp2 = xv(i_x2)
         yp2 = yv(i_y2+1)  
                
         xp3 = xv(i_x2+1)
         yp3 = yv(i_y2+1)    
              
         xp4 = xv(i_x2+1)
         yp4 = yv(i_y2)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*vt(i_x2, i_y2-1) + outA(i,2)*vt(i_x2, i_y2) + &
                      outA(i,3)*vt(i_x2+1, i_y2) + outA(i,4)*vt(i_x2+1, i_y2-1)
         ENDDO
         
         !**************************
         v_pos2 = ca0(1) + ca0(2)*pos2_x + ca0(3)*pos2_y + ca0(4)*pos2_x*pos2_y         
         !***************************
         
         xp1 = xv(i_x3) 
         yp1 = yv(i_y3) 
         
         xp2 = xv(i_x3)
         yp2 = yv(i_y3+1)  
                
         xp3 = xv(i_x3+1)
         yp3 = yv(i_y3+1)    
              
         xp4 = xv(i_x3+1)
         yp4 = yv(i_y3)
         
         invA(1,1) =  1._rk
         invA(1,2) = xp1
         invA(1,3) = yp1
         invA(1,4) = xp1*yp1
         
         invA(2,1) = 1._rk
         invA(2,2) = xp2
         invA(2,3) = yp2
         invA(2,4) = xp2*yp2
         
         invA(3,1) = 1._rk
         invA(3,2) = xp3
         invA(3,3) = yp3
         invA(3,4) = xp3*yp3
         
         invA(4,1) = 1._rk
         invA(4,2) = xp4
         invA(4,3) = yp4
         invA(4,4) = xp4*yp4
         
         Call matinv4(invA, outA)
         
         DO ii = 1, 4
            ca0(ii) = outA(i,1)*vt(i_x3, i_y3-1) + outA(i,2)*vt(i_x3, i_y3) + &
                      outA(i,3)*vt(i_x3+1, i_y3) + outA(i,4)*vt(i_x3+1, i_y3-1)
         ENDDO
         
         !*************************************
         v_pos3 = ca0(1) + ca0(2)*pos3_x + ca0(3)*pos3_y + ca0(4)*pos3_x*pos3_y    
         !*************************************
         
         
         !quadratic interpolation for correct pressure at intercepted cell
         IF (sur2nodeDis.GE.0._rk) THEN
             n1 = pt1 + sur2nodeDis
             n2 = pt2 + sur2nodeDis 
         ELSE
             n1 = pt2 - sur2nodeDis
             n2 = pt3 - sur2nodeDis
             v_pos1 = v_pos2
             v_pos2 = v_pos3                        
         ENDIF
         ns = 0._rk
         inv3A(1,1) =  n1*n1
         inv3A(1,2) =  n1
         inv3A(1,3) =  1._rk
         inv3A(2,1) =  n2*n2
         inv3A(2,2) =  n2
         inv3A(2,3) =  1._rk 
         inv3A(3,1) =  ns*ns
         inv3A(3,2) =  ns
         inv3A(3,3) =  1._rk    
         Call matinv3(inv3A, out3A)
         DO ii = 1, 3
            !print*, inv3A(ii, :)
            ca1(ii) = out3A(i,1)*v_pos1 + out3A(i,2)*v_pos2 + out3A(i,3)*vsurf
                        !print*, out3A(ii, :), 'o'
         ENDDO
         !print*, ca1(1), ca1(2), ca1(3), 'u'
         vt(i,j-1) = ca1(1)*sur2nodeDis**2 + ca1(2)*sur2nodeDis + ca1(3)
 1 continue
      ENDDO
      
END SUBROUTINE velocityForcing

!*********************************************************************
SUBROUTINE matinv4(A, B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    REAL (kind = 8), intent(in) :: A(4,4)   !! Matrix
    REAL (kind =8), intent(out)            :: B(4,4)   !! Inverse matrix
    REAL  (kind = 8)           :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
      1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
  end subroutine
!*********************************************************************  
  subroutine matinv3(A, B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    REAL(KIND=8), intent(in) :: A(3,3)   !! Matrix
    REAL(KIND=8)             :: B(3,3)   !! Inverse matrix
    REAL(KIND=8)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
  end subroutine
