
      SUBROUTINE poisson_solver_init(fluidCellCount, pc, pco, b, derr1,&
									  eps1,  dudt, dvdt, derrStdSt, nIterPcor, divmax)
        
        INTEGER(KIND=8) :: i, j, n
        INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        REAL (KIND = 8)    :: dalt,  div, dab, sStart, sStop, sStart1, sStop1, solaTime 
		INTEGER ( KIND = 8), intent(out)   ::nIterPcor
		INTEGER,  intent(in)   ::fluidCellCount
        REAL (KIND = 8), intent(out)    :: b(fluidCellCount), derr1, eps1,  dudt, dvdt, derrStdSt, divmax
		REAL (KIND = 8), intent(inout)    :: pc(:,:), pco(:,:)        
         !ALLOCATE ( b(fluidCellCount) )

        ! ALLOCATE (r(nx*ny), rs(nx*ny), p_((nx+2)*(ny+2)), ap(nx*ny), s((nx+2)*(ny+2)), as(nx*ny), &
        !          u_((nx+2)*(ny+2)), q((nx+2)*(ny+2)))          
         !pct = pct1
         !sumItaSola   = 0
         !solaTime  = 0._rk
         !sumIterPc = 0
         !pcItaMax  = 32
         eps1 = 1e-6
         derrStdSt = 0._rk
         !print*, eps1
 1       CONTINUE

         !initialize variables
         nIterPcor = 0  
         divmax = 0._rk
         dalt   = 0._rk         
         derr1  = 0._rk
         b(:)   = 0._rk 
         pc(:,:) = 0.0_rk
         pco(:,:) = 0.0_rk
       END SUBROUTINE poisson_solver_init      
	   
         ! CALL computeDiv    !divergence vector  
         ! !$acc data copyin(Acx, Acy, b, pc, pco, fluidIndexPtr )
         ! !CALL jacobi(eps1, nIterPcor, derr1) 
         ! !CALL SOR(eps1, nIterPcor, derr1)  
        ! CALL REDBLACKSOR(eps1, nIterPcor, derr1) 
                  ! !CALL BICGSTAB(0, eps1, nIterPcor, derr1) 
         ! !$acc update host(pc) 
         ! !$acc end data 
         ! !print*, ita, nIterPcor, derr1
         ! CALL correctPressure!pressure correction            
         ! CALL correctVelocity !velocity correction   
         ! CALL velocityBC      !correct velocity at boundaries
		 
		 
	  SUBROUTINE poisson_solver(fluidCellCount, fluidIndexPtr, dudt, dvdt,&
								ut, vt, deltat, derrStdSt, u, v, nIterPcor, derr1, ita)
		INTEGER,  intent(in)   ::fluidCellCount
		INTEGER (KIND = 8), intent(in):: fluidIndexPtr(:,:), nIterPcor, ita
		REAL (KIND = 8), intent(in)    :: derr1,   deltat
		REAL (KIND = 8), intent(inout)    :: ut(:,:), vt(:,:), u(:,:), v(:,:), dudt, dvdt, derrStdSt
         DO n = 1, fluidCellCount
           i = fluidIndexPtr(n, 1)
           j = fluidIndexPtr(n, 2)  
           dudt = dabs((ut(i,j) - u(i,j)))/deltat
           dvdt = dabs((vt(i,j) - v(i,j)))/deltat
           derrStdSt = dmax1(derrStdSt, dudt, dvdt)                           
           
         ENDDO
         WRITE(*,16) ita, nIterPcor, derr1, derrStdSt   
 16      FORMAT(' ',I8, I10, 2E15.6)
         PRINT*, '********************************'
         u(:,:) = ut(:,:)
         v(:,:) = vt(:,:)
         !DEALLOCATE (b)
      END SUBROUTINE poisson_solver
! !*******************************************************************
      SUBROUTINE compute_div(fluidCellCount, fluidIndexPtr, b, ut, deltax, vt, deltay)
        

         INTEGER :: n, i, j
		 INTEGER,  intent(in)   ::fluidCellCount
		 INTEGER (KIND = 8), intent(in):: fluidIndexPtr(:,:)
		 REAL (KIND = 8), intent(in)    :: ut(:,:), vt(:,:)
		 REAL (KIND = 8), intent(in)      :: deltax(:), deltay(:)
		 REAL (KIND = 8), intent(inout)   :: b(:)
         DO 70 n = 1, fluidCellCount
           i = fluidIndexPtr(n, 1)
           j = fluidIndexPtr(n, 2)     
           b(n) =  ( ut(i,j) - ut(i-1,j) )/deltax(i) +  ( vt(i,j) - vt(i,j-1) )/deltay(j) 
           !print*, b(n), i, j, n, ut(i,j) - ut(i-1,j) , vt(i,j) - vt(i,j-1)
 70      CONTINUE
      END SUBROUTINE compute_div
! !*******************************************************************

! !*******************************************************************

      SUBROUTINE correct_pressure(fluidCellCount, fluidIndexPtr, p, pc)
         INTEGER,  intent(in)   ::fluidCellCount
		 INTEGER (KIND = 8), intent(in):: fluidIndexPtr(:,:)
		 REAL (KIND = 8), intent(inout)    :: p(:,:)
		 REAL (KIND = 8), intent(in)    :: pc(:,:)
         INTEGER(KIND=8) :: n, i, j
         DO 30 n = 1, fluidCellCount 
            i = fluidIndexPtr(n, 1)
            j = fluidIndexPtr(n, 2)              
            p(i,j)  = p(i,j) + pc(i,j) !- (0.5d0*( p(2, ny/2) + p(2, ny/2+1)))
 30      CONTINUE   
      END SUBROUTINE correct_pressure

! !c *********************************************************************

      SUBROUTINE correct_velocity(fluidCellCount, fluidIndexPtr, &
								  deltat, deltax, deltay, ut, vt, pc)
         INTEGER,  intent(in)   ::fluidCellCount
		 INTEGER (KIND = 8), intent(in):: fluidIndexPtr(:,:)
		 REAL (KIND = 8), intent(in)      :: deltax(:), deltay(:), deltat
		 REAL (KIND = 8), intent(inout)    :: ut(:,:), vt(:,:)
		 REAL (KIND = 8), intent(in)    :: pc(:,:)
         INTEGER(KIND=8) :: n, i, j
         DO 30 n = 1, fluidCellCount 
            i = fluidIndexPtr(n, 1)
            j = fluidIndexPtr(n, 2)      
            ut(i,j) = ut(i,j)-deltat/deltax(i)*(pc(i+1,j)-pc(i,j))
            vt(i,j) = vt(i,j)-deltat/deltay(j)*(pc(i,j+1)-pc(i,j))                  
 30      CONTINUE
      END SUBROUTINE correct_velocity
      
! !c ******************************************************************************

      ! SUBROUTINE JACOBI(epsi, isum, derr)      
         ! USE global
         ! IMPLICIT NONE
         ! INTEGER, PARAMETER :: rk = selected_real_kind(8)   
         ! INTEGER(KIND=8) :: n, i, j
         ! REAL (KIND = 8), INTENT(IN)     :: epsi
         ! REAL (KIND = 8), INTENT(OUT)    :: derr
         ! INTEGER (KIND = 8), INTENT(OUT) :: isum
         ! REAL (KIND = 8)                 :: dres 
         ! isum = 0   
         
 ! 3       isum = isum + 1  
          ! !CALL PBC  
         ! derr = 0._rk  
         ! dres = 0._rk
        ! !$acc parallel loop present(pc,b, Acx, Acy ,pco, fluidIndexPtr) firstprivate(deltat,fluidCellCount)   
         ! DO 10 n = 1, fluidCellCount   
           ! i = fluidIndexPtr(n, 1)
           ! j = fluidIndexPtr(n, 2)            
           ! pc(i,j) = (b(n)/deltat &
     ! &         - Acy(j-1,1)*pco(i,j-1) - Acx(i-1,1)*pco(i-1,j) &
     ! &         - Acx(i-1,3)*pco(i+1,j) - Acy(j-1,3)*pco(i,j+1))  / (Acx(i-1,2)+ Acy(j-1,2))  
           ! !print*, i, j, pc(i,j), isum, n, b(n) 
 ! 10      CONTINUE
         ! !$acc end parallel      
         

         ! !$acc parallel loop reduction(max:derr) present(pc, pco, fluidIndexPtr) firstprivate(fluidCellCount)  
          ! DO 20 n = 1, fluidCellCount   
             ! i = fluidIndexPtr(n, 1)
             ! j = fluidIndexPtr(n, 2) 
             ! derr = dmax1(derr,abs(pc(i,j)-pco(i,j)))                           

             ! !print*, pco(i,j), pc(i,j)
 ! 20      CONTINUE
         ! !$acc end parallel 
         ! !$acc parallel loop present(pc, pco, fluidIndexPtr) firstprivate(fluidCellCount)      
          ! DO 30 n = 1, fluidCellCount   
             ! i = fluidIndexPtr(n, 1)
             ! j = fluidIndexPtr(n, 2) 
             ! pco(i,j) = pc(i,j)
 ! 30      CONTINUE
         ! !$acc end parallel 
         ! IF (mod(isum,50000).EQ.0)WRITE(*,*) isum, derr
         ! !IF (derr.gt.1e-5.AND.isum.lt.pcItaMax) GOTO 3
         ! IF (derr.gt.epsi) GOTO 3
         ! !if(ita.ge.2) stop     
         ! !IF (ita.lt.5.AND.isum.lt.1500) GOTO 3         
      ! END SUBROUTINE JACOBI

! !*********************************************************************
! !*********************************************************************
          
      SUBROUTINE redblacksor(epsi, isum, derr, fluidCellCount, fluidIndexPtr, pc, pco, b, deltat, Acy, Acx)      
         
         INTEGER, PARAMETER :: rk = selected_real_kind(8)   
         INTEGER(KIND=8) :: n, i, j
         REAL (KIND = 8) :: omega, dres, derr2
		 REAL (KIND = 8), intent(in)    :: deltat
		 INTEGER,  intent(in)   ::fluidCellCount
		 INTEGER (KIND = 8), intent(in):: fluidIndexPtr(:,:)
         REAL (KIND = 8), intent(in)     :: Acy(:,:), Acx(:,:), b(:), epsi
         REAL (KIND = 8), intent(inout)    ::  pc(:,:), pco(:,:)
         INTEGER (KIND = 8), intent(out) :: isum
		 REAL (KIND = 8), intent(out)    :: derr
         !epsi=1e-6
		 isum = 0   
         omega = 1.667_rk      
 3       isum = isum + 1  
         derr = 0._rk 
         derr2 = 0._rk
         dres = 0._rk
         !print*, 1111
        !$acc parallel loop present(pc,b, Acx, Acy ,pco, fluidIndexPtr) firstprivate(deltat,fluidCellCount)   
         DO 10 n = 1, fluidCellCount   
             i = fluidIndexPtr(n, 1)
             j = fluidIndexPtr(n, 2)  
            IF (MOD(i+j,2).eq.1) THEN       
            pc(i,j) = (b(n)/deltat &
     &         - Acy(j-1,1)*pco(i,j-1) - Acx(i-1,1)*pco(i-1,j) &
     &         - Acx(i-1,3)*pco(i+1,j) - Acy(j-1,3)*pco(i,j+1))  / (Acx(i-1,2)+ Acy(j-1,2)) 
            pc(i,j) = (1._rk-omega)*pco(i,j) + omega*pc(i,j)
            ENDIF
 10      CONTINUE 
         !$acc end parallel  
         !print*, 222   
        !$acc parallel loop present(pc,b, Acx, Acy ,pco, fluidIndexPtr) firstprivate(deltat,fluidCellCount)   
         DO 20 n = 1, fluidCellCount   
             i = fluidIndexPtr(n, 1)
             j = fluidIndexPtr(n, 2)   
            IF (MOD(i+j,2).eq.0) THEN
                pc(i,j) = (b(n)/deltat &
     &        - Acy(j-1,1)*pc(i,j-1) - Acx(i-1,1)*pc(i-1,j) &
     &        - Acx(i-1,3)*pc(i+1,j) - Acy(j-1,3)*pc(i,j+1) ) / (Acx(i-1,2)+ Acy(j-1,2)) 
                pc(i,j) = omega*pc(i,j) + (1._rk-omega)*pco(i,j)
            ENDIF
 20      CONTINUE          
         !$acc end parallel  
              
         !$acc parallel loop reduction(max:derr) present(pc, pco, fluidIndexPtr) firstprivate(fluidCellCount)  
         DO 30 n = 1, fluidCellCount   
             i = fluidIndexPtr(n, 1)
             j = fluidIndexPtr(n, 2)   
            derr = dmax1(derr,abs(pc(i,j)-pco(i,j)))            
            pco(i,j) = pc(i,j)  
 30      CONTINUE 
         !$acc end parallel
           
         !IF (mod(isum,25000).EQ.0)WRITE(*,*) isum, derr
         !IF (ita.LE.2.AND.isum.LT.50000) GOTO 3
         !IF (derr.gt.epsi) GOTO 3
		 !print*, 'inside fortran derr',derr,'epsi', epsi,'isum', isum
         IF (derr.GE.epsi) GOTO 3
         !IF (ita.lt.15.AND.isum.lt.100) GOTO 3         
      END SUBROUTINE redblacksor
! !c ******************************************************************************

      ! SUBROUTINE SOR(epsi, isum, derr)      
         ! USE global
         ! IMPLICIT NONE
         ! INTEGER, PARAMETER :: rk = selected_real_kind(8)   
         ! INTEGER(KIND=8) :: n, i, j
         ! REAL (KIND = 8), INTENT(IN)     :: epsi
         ! REAL (KIND = 8), INTENT(OUT)    :: derr
         ! INTEGER (KIND = 8), INTENT(OUT) :: isum
         ! REAL (KIND = 8)                 :: dres, omega
         ! isum = 0   
         ! omega = 1.667_rk
 ! 3       isum = isum + 1  
          ! !CALL PBC  
         ! derr = 0._rk  
         ! dres = 0._rk
        ! !!$acc parallel loop collapse(2) present(pc,b,Acx, Acy ,pco) firstprivate(deltat)   
         ! DO 10 n = 1, fluidCellCount   
           ! i = fluidIndexPtr(n, 1)
           ! j = fluidIndexPtr(n, 2)            
           ! pc(i,j) = (b(n)/deltat &
     ! &         - Acy(j-1,1)*pc(i,j-1) - Acx(i-1,1)*pc(i-1,j) &
     ! &         - Acx(i-1,3)*pco(i+1,j) - Acy(j-1,3)*pco(i,j+1))  / (Acx(i-1,2)+ Acy(j-1,2))  
           ! pc(i,j) = omega*pc(i,j) + (1._rk-omega)*pco(i,j) 
           ! !print*, i, j, pc(i,j), isum, n, b(n) 
 ! 10      CONTINUE
         ! !!$acc end parallel      
         

         ! !!$acc parallel loop collapse(2) reduction(max:derr) present(pc, pco)    
          ! DO 20 n = 1, fluidCellCount   
             ! i = fluidIndexPtr(n, 1)
             ! j = fluidIndexPtr(n, 2) 
             ! derr = dmax1(derr,abs(pc(i,j)-pco(i,j)))                           
             ! pco(i,j) = pc(i,j) 
             ! !print*, pco(i,j), pc(i,j)
 ! 20      CONTINUE
         ! !!$acc end parallel 
         ! IF (mod(isum,1000).EQ.0)WRITE(*,*) isum, derr
         ! !IF (derr.gt.1e-5.AND.isum.lt.pcItaMax) GOTO 3
         ! IF (derr.gt.epsi) GOTO 3
         ! !if(ita.ge.2) stop     
         ! !IF (ita.lt.5.AND.isum.lt.1500) GOTO 3         
      ! END SUBROUTINE SOR

! !*********************************************************************


! !c ******************************************************************************
      ! SUBROUTINE BICGSTAB(id, epsi, isum, norm)
        ! USE global        
        ! IMPLICIT NONE   
        ! INTEGER, PARAMETER :: rk = selected_real_kind(8) 
        ! INTEGER , INTENT(IN) :: id
        ! REAL (KIND = 8), INTENT(IN)     :: epsi
        ! REAL (KIND = 8), INTENT(OUT)    :: norm
        ! INTEGER (KIND = 8), INTENT(OUT) :: isum
        ! REAL (KIND = 8)  :: temp, aprs, rrs, alpha1, ass, asas, r1rs, omega, beta
        
        ! INTEGER (KIND = 8) :: i, j, n, n1
        ! n = 0
        ! n1 = 0
        ! p_ = 0._rk
        ! r  = 0._rk
        ! rs = 0._rk
        ! ap = 0._rk
        ! s  = 0._rk
        ! temp = 0._rk 
        ! aprs = 0._rk
        ! rrs  = 0._rk
        ! alpha1 = 0._rk
        ! ass = 0._rk
        ! asas = 0._rk
        ! r1rs = 0._rk
        ! omega = 0._rk
        ! beta = 0._rk
    ! !$acc parallel loop collapse(2) present(cell, b, Acx, Acy, pc) copyout(r, p_, rs) firstprivate(deltat)    
	     ! DO 10 j = 2, ny+1
        ! DO 10 i = 2, nx+1
         ! n1    = i + (nx+2)*(j-1)
         ! n     = i-1 + nx*(j-2)
         ! !IF (cell(n).NE.1) THEN
         ! IF (cell(n).EQ.id) THEN
	        ! r(n)  = b(n)/deltat - ( Acx(i-1,3)*pc(i+1,j) + Acx(i-1,1)*pc(i-1,j) &
     ! &                           + (Acx(i-1,2)+ Acy(j-1,2))*pc(i,j)                         &
     ! &                           + Acy(j-1,3)*pc(i,j+1) + Acy(j-1,1)*pc(i,j-1) )
           ! p_(n1) = r(n)
           ! rs(n)  = r(n)
           ! !print*, r(n)
        ! END IF
! 10      CONTINUE
    ! !$acc end parallel  
    ! !!$acc update host(r, p_, rs)    
         
       ! isum = 0
       ! !if(ita.ge.1) stop     
 ! 1     CONTINUE 
        ! norm  = -1000._rk
    ! !$acc parallel loop collapse(2) present(Acx, Acy, cell) copyin(p_) copyout(ap)   
        ! DO 20 j = 2, ny+1
        ! DO 20 i = 2, nx+1

         ! n1    = i + (nx+2)*(j-1)
         ! n     = i-1 + nx*(j-2)
         ! !IF (cell(n).NE.1) &
         ! IF (cell(n).EQ.id) &

           ! ap(n) =  Acx(i-1,1)*p_(n1-(nx+2)) + Acy(j-1,1)*p_(n1-1)    &
     ! &           + (Acx(i-1,2)+ Acy(j-1,2))*p_(n1)                              &
     ! &           +  Acx(i-1,3)*p_(n1+1)      + Acy(j-1,3)*p_(n1+(nx+2))
 ! 20     CONTINUE  
    ! !$acc end parallel  
    ! !!$acc update host(ap)             
        ! aprs = 0.
        ! rrs  = 0.
    ! !$acc parallel loop copyin(ap, r, rs) present(cell) reduction(+:aprs, rrs) 
         ! DO n = 1, (nx)*(ny)
         ! !IF (cell(n).NE.1) THEN
         ! IF (cell(n).EQ.id) THEN

	         ! aprs = aprs + ap(n)*rs(n)
	         ! rrs  = rrs  + r(n)*rs(n)
         ! END IF
        ! END DO
    ! !$acc end parallel
          ! !print*, aprs, rrs
	     ! alpha1 = 0._rk
	     ! alpha1 = rrs/aprs

    ! !$acc parallel loop collapse(2) copyin(r, ap) present(cell) firstprivate(alpha1) copyout(s)
        ! DO 12 j = 2, ny+1
        ! DO 12 i = 2, nx+1 
           ! n1    = i + (nx+2)*(j-1)
           ! n     = i-1 + nx*(j-2)
          ! ! IF (cell(n).NE.1) &
         ! IF (cell(n).EQ.id) &
            ! s(n1) = r(n) - alpha1*ap(n)
 ! 12     CONTINUE      
    ! !$acc end parallel  
       
    ! !$acc parallel loop collapse(2) present(cell, Acx, Acy) copyin(s) copyout(as)  
	     ! DO 30 j = 2, ny+1
        ! DO 30 i = 2, nx+1
           ! n   = i-1 + nx*(j-2)
           ! n1  = i + (nx+2)*(j-1)          
           ! !IF (cell(n).NE.1) &
         ! IF (cell(n).EQ.id) &
           ! as(n)  =   Acx(i-1,1)*s(n1-(nx+2)) + Acy(j-1,1)*s(n1-1)     &
     ! &              + (Acx(i-1,2)+ Acy(j-1,2))*s(n1)                              &
     ! &              + Acx(i-1,3)*s(n1+1)      + Acy(j-1,3)*s(n1+nx+2) 

 ! 30     CONTINUE  
    ! !$acc end parallel     
	     ! ass   = 0._rk
	     ! asas  = 0._rk
    ! !$acc parallel loop collapse(2) copyin(as, s) present(cell) reduction(+:ass, asas) 
	     ! DO 31 j = 2, ny+1
        ! DO 31 i = 2, nx+1
         ! n1    = i + (nx+2)*(j-1)
         ! n     = i-1 + nx*(j-2)
         ! !IF (cell(n).NE.1) THEN
         ! IF (cell(n).EQ.id) THEN 
	         ! ass   = ass  + as(n)*s(n1)
	         ! asas  = asas + as(n)*as(n)
         ! END IF
 ! 31     CONTINUE   
    ! !$acc end parallel      
 
        ! omega = ass/asas      
    ! !$acc parallel loop collapse(2) present(cell) copyin(s, as, p_) firstprivate(alpha1,omega) copyout(r) copy(pc)
	     ! DO 40 j = 2, ny+1
        ! DO 40 i = 2, nx+1
         ! n1   = i + (nx+2)*(j-1)
         ! n    = i-1 + nx*(j-2)  
         ! !IF (cell(n).NE.1) THEN 
         ! IF (cell(n).EQ.id) THEN
  
	        ! r(n)    = s(n1) - omega*as(n)
           ! pc(i,j) = pc(i,j) + alpha1*p_(n1)  + omega*s(n1)
         ! END IF
 ! 40     CONTINUE  
    ! !$acc end parallel   

    ! !$acc parallel loop copyin(r) present(cell) reduction(max: norm) 
        ! DO 50 n = 1, nx*ny
          ! !IF (cell(n).NE.1) &
         ! IF (cell(n).EQ.id) &

          ! norm = dmax1(norm ,abs(r(n)))
 ! 50     CONTINUE 
    ! !$acc end parallel  
 
	     ! r1rs  = 0._rk
   ! !$acc parallel loop copyin(r, rs) present(cell) reduction(+:r1rs)
	     ! DO n = 1, (nx)*(ny)
             ! !  IF (cell(n).NE.1) &
         ! IF (cell(n).EQ.id) &

	         ! r1rs  = r1rs  + r(n)*rs(n)
	     ! END DO
    ! !$acc end parallel     
  
	     ! beta = (r1rs/rrs)*(alpha1/omega)
    ! !$acc parallel loop collapse(2) present(cell) firstprivate(beta,omega) copyin(ap, r) copy(p_)
	     ! DO 60 j = 2, ny+1
        ! DO 60 i = 2, nx+1
           ! n    = i-1 + nx*(j-2)
           ! n1   = i + (nx+2)*(j-1)
! !	        IF (cell(n).NE.1) &
         ! IF (cell(n).EQ.id) &
  
	        ! p_(n1)  = r(n)  + beta*(p_(n1) - omega*ap(n))
 ! 60     CONTINUE	
    ! !$acc end parallel 
	     ! isum = isum + 1
               ! ! WRITE(*,*) isum, norm      
	     ! If (norm.gt.epsi) GOTO 1

      ! END SUBROUTINE BICGSTAB
