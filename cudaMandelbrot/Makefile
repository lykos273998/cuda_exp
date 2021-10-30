f1=cuM
inc_path=pgm_res
f1: $(f1).cu
	nvcc -Xcompiler -fopenmp $(f1).cu -o $(f1) -arch sm_50 -I ${inc_path} -O3

prof: $(f1).cu
	nvcc -Xcompiler -fopenmp $(f1).cu -o $(f1) -arch sm_50 -I ${inc_path}  -DPROF -O3

r1:
	./$(f1)