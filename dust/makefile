COMPILER = gfortran
OPTIONS = -c -ffree-line-length-none -fPIC -O3
  
pydust: kappa_dust.o

clean: 
	find . -maxdepth 2 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete ;
	find ../ -maxdepth 2 -name "pydust.so" -delete ; find . -maxdepth 1 -name "*.c" -delete
	
kappa_dust.o: kappa_dust.f90
	$(COMPILER) $(OPTIONS)  kappa_dust.f90
