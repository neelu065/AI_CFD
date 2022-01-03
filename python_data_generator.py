import numpy as np
import subReadInput
import subAllocateArrays
import subSearch
import subInitialConditions
import subCoefficientMatrix
import subNavierStokes
import subBoundaryConditions
import forcing
import subPcorVcor
import subWriteOutput
import cmath

input_filename = 'inputdata'
gridx_filename = 'gridx.dat'
gridy_filename = 'gridy.dat'
geometry_filename = 'geometries/cyl1.msh'

### following lines including loop represents the list generator
input_data = "X"
parenthesis = "[]"
input_size = 19

for i in range(1, input_size + 1):
    command_variable = ""
    command_variable = input_data + str(i) + " = " + parenthesis
    exec(command_variable)

# X1 = []
# X2 = []
# X3 = []
# X4 = []
# X5 = []
# X6 = []
# X7 = []
# X8 = []
# X9 = []
# X10 = []
# X11 = []
# X12 = []
# X13 = []
# X14 = []
# X15 = []
# X16 = []
# X17 = []
# X18 = []
# X19 = []
Y = []

ita = 0
totime = 0
tempTime = 0
atStart = 0

nx, ny, dt_order, a0, freq, u0, v0, re, fx, fy, xShift, yShift, itamax, pcItaMax, itaSola, surGeoPoints, istart = subReadInput.read_input(
    input_filename,
    gridx_filename,
    gridy_filename)
min_sur_dis = np.zeros((nx + 2, ny + 2))

deltax, deltay, x1, y1, xu, yu, xv, yv, xp, yp, alpha, deltat, pi = subReadInput.create_grid(input_filename,
                                                                                             gridx_filename,
                                                                                             gridy_filename, nx, ny,
                                                                                             dt_order, a0, freq, u0, v0,
                                                                                             re, fx,
                                                                                             fy, xShift, yShift, itamax,
                                                                                             pcItaMax, itaSola,
                                                                                             surGeoPoints, istart)

u, ut, v, vt, p, cell, Acx, Acy, r, rs, p_, ap, as_, s, u_, q, bcSurf, pc, pco = subAllocateArrays.allocatearrays(nx,
                                                                                                                  ny)

ibNodes, ibElems = subReadInput.read_geometry_nodes(geometry_filename, surGeoPoints)

xnode, ynode, ibElP1, ibElP2 = subReadInput.read_surface_mesh(geometry_filename, surGeoPoints, ibNodes, ibElems)

xnode1, ynode1, xt, xdot, xddot, yt, ydot, yddot = subSearch.shift_surface_nodes(ibNodes, xShift, yShift, xnode, ynode)

xcent, ycent, cosAlpha, cosBeta, area = subSearch.compute_surface_norm(ibElems, xnode, ynode, ibElP1, ibElP2)

ibCellCount, fluidCellCount, solidCellCount, nodeId = subSearch.tagging(nx, ny, x1, y1, ibElems, xcent, ycent, cosAlpha,
                                                                        cosBeta, cell, xp, yp)

interceptedIndexPtr, pNormDis, nel2p, u1NormDis, u2NormDis, v1NormDis, v2NormDis, fluidIndexPtr, solidIndexPtr = subSearch.cell_count(
    nx, ny, ibCellCount,
    fluidCellCount, solidCellCount,
    cell)

subSearch.compute_norm_distance(ibElems, ibCellCount, xp, yp, nel2p, xcent, ycent, pNormDis, cosAlpha, cosBeta, yu, xu,
                                interceptedIndexPtr, u1NormDis, u2NormDis,
                                xv, yv, v1NormDis, v2NormDis)

if atStart == 0:
    subInitialConditions.initial_conditions(fluidCellCount, fluidIndexPtr, u, ut, v, vt, p)

subCoefficientMatrix.coefficient_matrix(nx, ny, Acx, Acy, xp, yp)

for iterationnn in range(1001):

    deltx2, delty2 = subNavierStokes.ns_momentum(fluidCellCount, fluidIndexPtr, deltax, deltay, deltat, u, v, p, ut, vt,
                                                 alpha, re)

    subBoundaryConditions.velocity_bc(nx, ny, cell, ut, vt)

    subBoundaryConditions.solid_cell_bc(solidCellCount, solidIndexPtr, ut, vt, p)

    forcing.velocity_forcing_1(ibCellCount, xdot, ydot, interceptedIndexPtr, u2NormDis,
                               deltax, deltay, xu, cosAlpha, nel2p, yu, cosBeta, ut, u1NormDis,
                               vt, v2NormDis, v1NormDis, xv, yv)

    if iterationnn > 0:
        xy_present = []

        for k in range(ibCellCount):
            i = interceptedIndexPtr[k, 0] - 1
            j = interceptedIndexPtr[k, 1] - 1

            for jj in range(j - 4, j + 5):
                for ii in range(i - 4, i + 5):
                    n = ii + (nx) * (jj - 1)

                    if (cell[n - 1] == 0 or cell[n - 1] == 2):
                        if (n not in xy_present):
                            m = 0
                            minDis = 1e14
                            n1x = xp[ii]
                            n1y = yp[jj]

                            for m in range(ibElems):
                                cent_x = xcent[m]
                                cent_y = ycent[m]
                                dis = (
                                    cmath.sqrt((n1y - cent_y) * (n1y - cent_y) + (n1x - cent_x) * (n1x - cent_x))).real
                                if (dis < minDis):
                                    minDis = dis
                                    nel2n = m

                            X1 = np.append(X1, p[ii, jj])
                            X2 = np.append(X2, (
                                    (ut[ii, jj] - ut[ii - 1, jj]) / deltax[ii] + (vt[ii, jj] - vt[ii, jj - 1]) /
                                    deltay[jj]))
                            X3 = np.append(X3, minDis)
                            X4 = np.append(X4, -xddot)
                            X5 = np.append(X5, 0.002)
                            X6 = np.append(X6, xp[ii])
                            X7 = np.append(X7, yp[jj])
                            X8 = np.append(X8, deltax[ii])
                            X9 = np.append(X9, deltay[jj])
                            X10 = np.append(X10, re)
                            X11 = np.append(X11, ut[ii - 1, jj])
                            X12 = np.append(X12, ut[ii, jj])
                            X13 = np.append(X13, vt[ii, jj])
                            X14 = np.append(X14, vt[ii, jj - 1])
                            X15 = np.append(X15, a0)
                            X16 = np.append(X16, freq)
                            X17 = np.append(X17, xcent[nel2n])
                            X18 = np.append(X18, ycent[nel2n])
                            X19 = np.append(X19, totime - tempTime)
                            xy_present = np.append(xy_present, n)
        #    print('xy_present',len(xy_present))
        del xy_present
    ##

    #
    b, derr1, eps1, dudt, dvdt, derrStdSt, nIterPcor, divmax = subPcorVcor.poisson_solver_init(fluidCellCount, pc, pco)

    subPcorVcor.compute_div(fluidCellCount, fluidIndexPtr, b, ut, deltax, vt, deltay)

    epsi = 0.000001
    isum, derr = subPcorVcor.redblacksor(epsi, fluidCellCount, fluidIndexPtr, pc, pco, b, deltat, Acy, Acx)

    subPcorVcor.correct_pressure(fluidCellCount, fluidIndexPtr, p, pc)

    subPcorVcor.correct_velocity(fluidCellCount, fluidIndexPtr, deltat, deltax, deltay, ut, vt, pc)

    subBoundaryConditions.velocity_bc(nx, ny, cell, ut, vt)

    nIterPcor = isum
    derr1 = derr

    subPcorVcor.poisson_solver(fluidCellCount, fluidIndexPtr, dudt, dvdt, ut, vt, deltat, derrStdSt, u, v, nIterPcor,
                               derr1, ita)
    del b

    forcing.pressure_forcing_1(ibCellCount, interceptedIndexPtr, ita, nel2p, xcent, xp, pNormDis, xddot, p, deltax,
                               deltay, cosAlpha, cosBeta, yp)

    subBoundaryConditions.solid_cell_bc(solidCellCount, solidIndexPtr, ut, vt, p)
    if (iterationnn > 0):
        xy_present = []

        for k in range(ibCellCount):
            i = interceptedIndexPtr[k, 0] - 1
            j = interceptedIndexPtr[k, 1] - 1

            for jj in range(j - 4, j + 5):
                for ii in range(i - 4, i + 5):
                    n = ii + (nx) * (jj - 1)

                    if cell[n - 1] == 0 or cell[n - 1] == 2:
                        if (n not in xy_present):
                            xy_present = np.append(xy_present, n)
                            Y = np.append(Y, p[ii, jj])
        del xy_present

    subWriteOutput.fft_data(u, totime, ita, v, p)
    subWriteOutput.stress_2d(ibElems, nx, ny, xcent, ycent, x1, y1, cosAlpha, cosBeta, xu, yu, xv, yv, xp, yp, ut, vt,
                             re, area, p, ita, totime)
    subWriteOutput.write_output(ita, nx, ny, xp, yp, u, v, p, totime, cell)
    ita = ita + 1
    atStart = atStart + 1
    totime = ita * deltat

    del xcent, ycent, cosAlpha, cosBeta, area
    xt_temp, xdot_temp, xddot_temp, yt_temp, ydot_temp, yddot_temp = subSearch.compute_surface_variables(pi, freq, a0,
                                                                                                         totime,
                                                                                                         tempTime,
                                                                                                         xnode1, ynode1,
                                                                                                         xnode, ynode)
    xt = xt_temp
    xdot = xdot_temp
    xddot = xddot_temp
    yt = yt_temp
    ydot = ydot_temp
    yddot = yddot_temp

    xcent, ycent, cosAlpha, cosBeta, area = subSearch.compute_surface_norm(ibElems, xnode, ynode, ibElP1, ibElP2)

    ibCellCount_temp, fluidCellCount_temp, solidCellCount_temp = subSearch.selective_retagging(ibCellCount,
                                                                                               interceptedIndexPtr,
                                                                                               nodeId, x1, y1, ibElems,
                                                                                               xcent, ycent,
                                                                                               cosAlpha, cosBeta,
                                                                                               solidCellCount,
                                                                                               fluidCellCount, nx, ny,
                                                                                               cell, xp, yp)
    ibCellCount = ibCellCount_temp
    fluidCellCount = fluidCellCount_temp
    solidCellCount = solidCellCount_temp

    del interceptedIndexPtr, pNormDis, nel2p, u1NormDis, u2NormDis, v1NormDis, v2NormDis, fluidIndexPtr, solidIndexPtr
    interceptedIndexPtr, pNormDis, nel2p, u1NormDis, u2NormDis, v1NormDis, v2NormDis, fluidIndexPtr, solidIndexPtr = subSearch.cell_count(
        nx, ny, ibCellCount,
        fluidCellCount, solidCellCount,
        cell)

    subSearch.compute_norm_distance(ibElems, ibCellCount, xp, yp, nel2p, xcent, ycent, pNormDis, cosAlpha, cosBeta, yu,
                                    xu, interceptedIndexPtr, u1NormDis,
                                    u2NormDis, xv, yv, v1NormDis, v2NormDis)
print('len(X1)', len(X1), len(Y))
from xlwt import Workbook

wb = Workbook()

sheet1 = wb.add_sheet('Sheet 1')
for kkkk in range(1):
    sheet1.write(kkkk + 1, 0, X1[kkkk])
    sheet1.write(kkkk + 1, 1, X2[kkkk])
    sheet1.write(kkkk + 1, 2, X3[kkkk])
    sheet1.write(kkkk + 1, 3, X4[kkkk])
    sheet1.write(kkkk + 1, 4, X5[kkkk])
    sheet1.write(kkkk + 1, 5, X6[kkkk])
    sheet1.write(kkkk + 1, 6, X7[kkkk])
    sheet1.write(kkkk + 1, 7, X8[kkkk])
    sheet1.write(kkkk + 1, 8, X9[kkkk])
    sheet1.write(kkkk + 1, 9, X10[kkkk])
    sheet1.write(kkkk + 1, 10, X11[kkkk])
    sheet1.write(kkkk + 1, 11, X12[kkkk])
    sheet1.write(kkkk + 1, 12, X13[kkkk])
    sheet1.write(kkkk + 1, 13, X14[kkkk])
    sheet1.write(kkkk + 1, 14, X15[kkkk])
    sheet1.write(kkkk + 1, 15, X16[kkkk])
    sheet1.write(kkkk + 1, 16, X17[kkkk])
    sheet1.write(kkkk + 1, 17, X18[kkkk])
    sheet1.write(kkkk + 1, 18, X19[kkkk])
    sheet1.write(kkkk + 1, 19, Y[kkkk])
wb.save('training1.xls')
from xlwt import Workbook

wb = Workbook()

sheet1 = wb.add_sheet('Sheet 2')
for kkkk in range(64995, 64995 * 2):
    sheet1.write(kkkk + 1 - 64995, 0, X1[kkkk])
    sheet1.write(kkkk + 1 - 64995, 1, X2[kkkk])
    sheet1.write(kkkk + 1 - 64995, 2, X3[kkkk])
    sheet1.write(kkkk - 64995 + 1, 3, X4[kkkk])
    sheet1.write(kkkk - 64995 + 1, 4, X5[kkkk])
    sheet1.write(kkkk - 64995 + 1, 5, X6[kkkk])
    sheet1.write(kkkk - 64995 + 1, 6, X7[kkkk])
    sheet1.write(kkkk - 64995 + 1, 7, X8[kkkk])
    sheet1.write(kkkk - 64995 + 1, 8, X9[kkkk])
    sheet1.write(kkkk - 64995 + 1, 9, X10[kkkk])
    sheet1.write(kkkk - 64995 + 1, 10, X11[kkkk])
    sheet1.write(kkkk - 64995 + 1, 11, X12[kkkk])
    sheet1.write(kkkk - 64995 + 1, 12, X13[kkkk])
    sheet1.write(kkkk - 64995 + 1, 13, X14[kkkk])
    sheet1.write(kkkk - 64995 + 1, 14, X15[kkkk])
    sheet1.write(kkkk - 64995 + 1, 15, X16[kkkk])
    sheet1.write(kkkk - 64995 + 1, 16, X17[kkkk])
    sheet1.write(kkkk - 64995 + 1, 17, X18[kkkk])
    sheet1.write(kkkk - 64995 + 1, 18, X19[kkkk])
    sheet1.write(kkkk - 64995 + 1, 19, Y[kkkk])
wb.save('training2.xls')

from xlwt import Workbook

wb = Workbook()

sheet1 = wb.add_sheet('Sheet 3')
for kkkk in range(64995 * 2, 64995 * 3):
    sheet1.write(kkkk + 1 - 64995 * 2, 0, X1[kkkk])
    sheet1.write(kkkk + 1 - 64995 * 2, 1, X2[kkkk])
    sheet1.write(kkkk + 1 - 64995 * 2, 2, X3[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 3, X4[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 4, X5[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 5, X6[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 6, X7[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 7, X8[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 8, X9[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 9, X10[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 10, X11[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 11, X12[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 12, X13[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 13, X14[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 14, X15[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 15, X16[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 16, X17[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 17, X18[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 18, X19[kkkk])
    sheet1.write(kkkk - 64995 * 2 + 1, 19, Y[kkkk])
wb.save('training3.xls')

from xlwt import Workbook

wb = Workbook()

sheet1 = wb.add_sheet('Sheet 4')
for kkkk in range(64995 * 3, 64995 * 4):
    sheet1.write(kkkk + 1 - 64995 * 3, 0, X1[kkkk])
    sheet1.write(kkkk + 1 - 64995 * 3, 1, X2[kkkk])
    sheet1.write(kkkk + 1 - 64995 * 3, 2, X3[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 3, X4[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 4, X5[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 5, X6[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 6, X7[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 7, X8[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 8, X9[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 9, X10[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 10, X11[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 11, X12[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 12, X13[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 13, X14[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 14, X15[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 15, X16[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 16, X17[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 17, X18[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 18, X19[kkkk])
    sheet1.write(kkkk - 64995 * 3 + 1, 19, Y[kkkk])
wb.save('training4.xls')

from xlwt import Workbook

wb = Workbook()

sheet1 = wb.add_sheet('Sheet 5')
for kkkk in range(64995 * 4, 64995 * 5):
    sheet1.write(kkkk + 1 - 64995 * 4, 0, X1[kkkk])
    sheet1.write(kkkk + 1 - 64995 * 4, 1, X2[kkkk])
    sheet1.write(kkkk + 1 - 64995 * 4, 2, X3[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 3, X4[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 4, X5[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 5, X6[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 6, X7[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 7, X8[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 8, X9[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 9, X10[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 10, X11[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 11, X12[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 12, X13[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 13, X14[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 14, X15[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 15, X16[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 16, X17[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 17, X18[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 18, X19[kkkk])
    sheet1.write(kkkk - 64995 * 4 + 1, 19, Y[kkkk])
wb.save('training5.xls')

from xlwt import Workbook

wb = Workbook()

sheet1 = wb.add_sheet('Sheet 6')
for kkkk in range(64995 * 5, len(X1)):
    sheet1.write(kkkk + 1 - 64995 * 5, 0, X1[kkkk])
    sheet1.write(kkkk + 1 - 64995 * 5, 1, X2[kkkk])
    sheet1.write(kkkk + 1 - 64995 * 5, 2, X3[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 3, X4[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 4, X5[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 5, X6[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 6, X7[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 7, X8[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 8, X9[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 9, X10[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 10, X11[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 11, X12[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 12, X13[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 13, X14[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 14, X15[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 15, X16[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 16, X17[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 17, X18[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 18, X19[kkkk])
    sheet1.write(kkkk - 64995 * 5 + 1, 19, Y[kkkk])
wb.save('training6.xls')
