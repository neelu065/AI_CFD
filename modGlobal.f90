       MODULE global
       IMPLICIT NONE
       CHARACTER (LEN = 128) :: line
       INTEGER               :: istart, id1, ibCellCount, fluidCellCount, solidCellCount
       INTEGER ( KIND = 8)   :: nx, ny, itamax, pcItaMax, &
                                ita, nIterPcor, nc,     &
                                sumIterPc,itaSola, totIterPc
       REAL (KIND =8)        :: lx, ly,  dt_order,   & 
                                deltx2, delty2,  &
                                a0, freq, &
                                u0, v0, p0,         &  
                                eps, epsDiv, re, divmax,   &
                                fx, fy,  alpha, pct, pct1, &
                                deltat, totime, temptime, dfinish, dstart, solverTime, pi , msTime, mindx
       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:)    :: cell, cell1
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:)       :: b, x, y, deltax, deltay, x1, y1, xu, yu, xv, yv, xp, yp   
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:, :)    :: Ac, Acx, Acy
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:, :) :: u, ut,  &
                                                            v, vt,  &
                                                            p, pc, pco
       INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:, :) :: interceptedIndexPtr, fluidIndexPtr, solidIndexPtr, nodeId
       ! ibm variables
       INTEGER  :: bcType
       REAL (KIND = 8)    :: xShift, yShift, aoa, piv_pt, var_surf, &
       xt, xdot, xddot, yt, ydot, yddot
       INTEGER (KIND=8)   :: surGeoPoints, ibElems, ibNodes 
       INTEGER (KIND=8), ALLOCATABLE, DIMENSION (:)   :: ibElP1, ibElP2, nel2p
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:) :: xnode, ynode,  xnode1, ynode1, xcent, ycent, &
                                                      cosAlpha, cosBeta, bcSurf, area, &
                                                      pNormDis, u1NormDis, u2NormDis , v1NormDis, v2NormDis
       !REAL (KIND = 8), ALLOCATABLE, DIMENSION (:)   :: r, rs, p_, ap, s, as, q, u_
      ! INTEGER (KIND = 8), ALLOCATABLE, DIMENSION (:) :: nn2el, nodeId
       REAL (KIND = 8), ALLOCATABLE, DIMENSION (:)   :: r, rs, p_, ap, s, as, q, u_               
       END MODULE global
                
