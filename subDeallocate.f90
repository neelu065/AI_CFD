      SUBROUTINE deallocate_1(xcent, ycent, cosAlpha, cosBeta, area)
                                                            
        REAL (KIND = 8), intent(inout) :: xcent(:), ycent(:), cosAlpha(:), cosBeta(:), area(:)
        DEALLOCATE(xcent, ycent, cosAlpha, cosBeta, area)
        
      END SUBROUTINE deallocate_1
	  
	  SUBROUTINE deallocate_2(interceptedIndexPtr, pNormDis, nel2p, u1NormDis,&
					          u2NormDis , v1NormDis, v2NormDis, fluidIndexPtr, solidIndexPtr)     
                                                            
        REAL(KIND = 8), intent(inout)::pNormDis(:), u1NormDis(:), u2NormDis(:), v2NormDis(:),v1NormDis(:)
									   
		INTEGER (KIND = 8), intent(inout) :: interceptedIndexPtr(:,:), fluidIndexPtr(:,:), &
										     solidIndexPtr(:,:), nel2p(:)
		DEALLOCATE(interceptedIndexPtr, pNormDis, nel2p, u1NormDis,&
					  u2NormDis , v1NormDis, v2NormDis, fluidIndexPtr, solidIndexPtr)      
        
      END SUBROUTINE deallocate_2

