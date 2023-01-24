# Compiler
MPIFC = mpifort   
FC = gfortran  
IC = ifort

# The CFLAGS variable sets compile flags:
#  -g                            compile with debug information
#  -Wall                         give verbose compiler warnings
#  -O0                           do not optimize generated code
#  -O2 -ftree-vectorize          apply level 2 optimization with auto-vectorization only
#  -O3                           apply level 3 optimizations
#  -std=gnu99                    use the GNU99 standard language definition
#  -fopt-info-all=opt_report.all save optimization report to file
# -fopenmp (-qopenmp)            enable openmp for gnu (intel)
# -fcheck=all                    enable array bounds checking
# -CB                            intel compiler array bounds checkingg
# -ffpe-trap=invalid             stops program at first NaN (gnu compiler)
# -fpe0                          stops program at first NaN (intel compiler)

CFLAGS = -O3 -qopenmp 

OBJS = constants_mod.o grid_arrays_mod.o mover_mod.o file_io_mod.o particles_init_mod.o generate_bfield_mod.o test_particles_main.o 

test_particles_run : $(OBJS)
	$(IC) $(CFLAGS) -o ../test_particles_run $(OBJS)
	
constants_mod.o : ../src/constants_mod.f90
	$(IC) $(CFLAGS) -c ../src/constants_mod.f90

grid_arrays_mod.o : ../src/grid_arrays_mod.f90 constants_mod.o
	$(IC) $(CFLAGS) -c ../src/grid_arrays_mod.f90

generate_bfield_mod.o : ../src/generate_bfield_mod.f90 constants_mod.o grid_arrays_mod.o
	$(IC) $(CFLAGS) -c ../src/generate_bfield_mod.f90
    
mover_mod.o : ../src/mover_mod.f90 constants_mod.o grid_arrays_mod.o generate_bfield_mod.o
	$(IC) $(CFLAGS) -c ../src/mover_mod.f90				

file_io_mod.o : ../src/file_io_mod.f90 constants_mod.o grid_arrays_mod.o
	$(IC) $(CFLAGS) -c ../src/file_io_mod.f90		
    
particles_init_mod.o : ../src/particles_init_mod.f90 constants_mod.o grid_arrays_mod.o
	$(IC) $(CFLAGS) -c ../src/particles_init_mod.f90   
   
test_particles_main.o : ../src/test_particles_main.f90 constants_mod.o grid_arrays_mod.o mover_mod.o file_io_mod.o particles_init_mod.o generate_bfield_mod.o
	$(IC) $(CFLAGS) -c ../src/test_particles_main.f90
    
clean:
	rm  -f *.o *.mod 