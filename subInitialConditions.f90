
      SUBROUTINE initial_conditions( fluidCellCount, fluidIndexPtr, u, ut, v, vt, p)     
      
       INTEGER, PARAMETER :: rk = selected_real_kind(8)
	   INTEGER, intent(in)::fluidCellCount
	   INTEGER(KIND = 8), intent(in)::fluidIndexPtr(:, :)
	   
	   REAL (KIND = 8), intent(inout) :: u(:, :), ut(:, :),  &
                                         v(:, :), vt(:, :),  &
                                         p(:, :)
       INTEGER::  i, j, k, n
       REAL (KIND = 8) :: t1, t2
        k = 1 
        ! cell variables 
        u = 0._rk
        ut = 0._rk
        v = 0._rk
        vt = 0._rk
        p = 0._rk 
        DO k = 1, fluidCellCount
           i = fluidIndexPtr(k, 1)
             j = fluidIndexPtr(k, 2) 
           u(i,j) = 0._rk
          v(i,j) = 0._rk
          ut(i,j) =0._rk
          vt(i,j) =0._rk
         
        END DO
        !pc = 0.0_rk
         !pco(:,:) = 0.0_rk
        print*, 'initial'
 10     CONTINUE            
      END SUBROUTINE initial_conditions
