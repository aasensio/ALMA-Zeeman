COMPILER = gfortran
OPTIONS = -c -ffree-line-length-none -fPIC -O3

pyformal: formal.o
# 	$(COMPILER) $(OPTIONS) formal.f90
	
clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 2 -name "*.a" -delete ;
	find . -maxdepth 1 -name "pyformal.so" -delete ; find . -maxdepth 1 -name "*.pyf" -delete
	find ../ -maxdepth 2 -name "pyformal.so" -delete ; find . -maxdepth 1 -name "*.c" -delete
	
formal.o: formal.f90
	$(COMPILER) $(OPTIONS)  formal.f90
	