
      PROGRAM main
        USE global
        IMPLICIT NONE
        CALL cpu_time(dStart)
        PRINT*, 'Time start =', dStart 
        
        CALL readInput
        CALL allocateArrays
        CALL readSurfaceMesh
        CALL shiftSurfaceNodes
        !CALL comp22uteSurfaceVariables
        CALL computeSurfaceNorm
        CALL tagging
        CALL cellCount
        CALL computeNormDistance
        !stop
        IF (iStart.eq.0) CALL initialConditions        
        CALL coefficientMatrix
        ita = 0
        totime = 0
        tempTime = 0
 1      CONTINUE
        !GOTO 2

        CALL nsMomentum
        CALL velocityBC
        CALL solidCellBC
		DO 755 n = 1, fluidCellCount
         i = fluidIndexPtr(n, 1)
         j = fluidIndexPtr(n, 2) 
		 print*,'ut', ut(i, j)
755		 continue
        CALL velocityForcing1

        CALL poissonSolver
        PRINT*, 'CONTI DONE'
        CALL pressureForcing1
        PRINT*, 'P FORCING DONE'
        CALL solidCellBC
        PRINT*, 'SOLID BC DONE'

        CALL fftData
        PRINT*, 'FFT DONE'
        CALL stress2D
        PRINT*, 'STRESS DONE'        
        !CALL writeInlineData
 !2      CONTINUE
        CALL writeOutput
        PRINT*, 'OUTPUT DONE'                
        ita = ita + 1
        totime = ita*deltat
        !IF(ita.EQ.3000) tempTime = totime 
        !IF(ita.GT.3000) THEN
           !DEALLOCATE(xcent, ycent, cosAlpha, cosBeta, area)
           ! CALL computeSurfaceVariables   
           ! PRINT*, 'SURFACE BC DONE'                   
           ! CALL computeSurfaceNorm  
           ! PRINT*, 'SURFACE NORM DONE'                   
           ! !CALL tagging      
           ! CALL selectiveRetagging
           ! DEALLOCATE(interceptedIndexPtr,pNormDis,nel2p,u1NormDis, & 
           ! u2NormDis,v1NormDis,v2NormDis,fluidIndexPtr,solidIndexPtr)        
           ! CALL cellCount        
           ! CALL computeNormDistance
        !ENDIF        
        !stop
        !pause
        GOTO 1
      END PROGRAM main