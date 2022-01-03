#making of make file
compiler=gfortran
cc=-c
name = 'sudden'
all:name
name:modGlobal.o subReadInput.o subAllocateArrays.o subInitialConditions.o subCheckPoint.o subCoefficientMatrix.o subSearch.o subNavierStokes.o subBoundaryConditions.o subPcorVcor.o subWriteOutput.o subStressCallinear.o subDeallocateArrays.o main.o
	$(compiler) modGlobal.o subReadInput.o subAllocateArrays.o subInitialConditions.o subCheckPoint.o subCoefficientMatrix.o subSearch.o subNavierStokes.o subBoundaryConditions.o subPcorVcor.o subWriteOutput.o subStressCal.o subDeallocateArrays.o main.o -o name
modGlobal.o:modGlobal.f90 
	$(compiler) $(cc) modGlobal.f90
subReadInput.o:subReadInput.f90 
	$(compiler) $(cc) subReadInput.f90
subAllocateArrays.o:subAllocateArrays.f90 
	$(compiler) $(cc) subAllocateArrays.f90
subInitialConditions.o:subInitialConditions.f90
	$(compiler) $(cc) subInitialConditions.f90
subCheckPoint.o:subCheckPoint.f90
	$(compiler) $(cc) subCheckPoint.f90
subCoefficientMatrix.o:subCoefficientMatrix.f90
	$(compiler) $(cc) subCoefficientMatrix.f90
subSearch.o:subSearch.f90
	$(compiler) $(cc) subSearch.f90
subNavierStokes.o:subNavierStokes.f90
	$(compiler) $(cc) subNavierStokes.f90
subBoundaryConditions.o:subBoundaryConditions.f90
	$(compiler) $(cc) subBoundaryConditions.f90
subPcorVcor.o:subPcorVcor.f90
	$(compiler) $(cc) subPcorVcor.f90
subWriteOutput.o:subWriteOutput.f90
	$(compiler) $(cc) subWriteOutput.f90
subStressCallinear.o:subStressCal.f90
	$(compiler) $(cc) subStressCal.f90
subDeallocateArrays.o:subDeallocateArrays.f90
	$(compiler) $(cc) subDeallocateArrays.f90
main.o:main.f90
	$(compiler) $(cc) main.f90
clean:
	rm *.o *.mod name	
