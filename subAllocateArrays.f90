      SUBROUTINE allocatearrays(nx, ny, u, ut, v, vt, p, cell, Acx, Acy, r, rs, p_, ap, as_, s, u_, q, bcSurf, pc, pco)
        INTEGER ( KIND = 8), intent(in) :: nx,ny
		REAL (KIND = 8), intent(out) :: u(nx+2,ny+2), ut(nx+2,ny+2), v(nx+2,ny+2),Acx(nx,3), Acy(ny,3),  &
                                        vt(nx+2,ny+2), p(nx+2,ny+2),  pc(nx+2,ny+2), pco(nx+2,ny+2), bcSurf(3),&
										r(nx*ny), rs(nx*ny), p_((nx+2)*(ny+2)), ap(nx*ny), s((nx+2)*(ny+2)), as_(nx*ny), &
										u_((nx+2)*(ny+2)), q((nx+2)*(ny+2))
										
        INTEGER (KIND = 8), intent(out) :: cell(nx*ny)                                                    
       
           
        
      END SUBROUTINE allocatearrays
