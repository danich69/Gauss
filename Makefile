calculate: integral.o integrator.o func.o
	gfortran -fcheck=all $^ -o $@ -g
mesh: main.o gauss.o linear_system.o func.o
	gfortran -fcheck=all $^ -o $@ -g
main.o: main.f90 gauss.mod linear_system.mod func.o
	gfortran -fcheck=all main.f90 -c -g
gauss.mod gauss.o: gauss.f90 linear_system.mod func.mod
	gfortran -fcheck=all $^ -c -g
linear_system.mod linear_system.o: linear_system.f90
	gfortran -fcheck=all $^ -c -g
func.mod func.o: func.f90
	gfortran -fcheck=all $^ -c -g
integrator.mod integrator.o: integrator.f90
	gfortran -fcheck=all $^ -c -g
integral.o: integral.f90 func.mod integrator.mod
	gfortran -fcheck=all $^ -c -g
clean: 
	rm -f *.o *mod
