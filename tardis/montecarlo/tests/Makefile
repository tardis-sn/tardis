all: test_cmontecarlo

rk_isaac.o: ../src/randomkit/rk_isaac.c
	g++ -fpermissive -c ../src/randomkit/rk_isaac.c

rk_mt.o: ../src/randomkit/rk_mt.c
	g++ -fpermissive -c ../src/randomkit/rk_mt.c

rk_primitive.o: ../src/randomkit/rk_primitive.c
	g++ -fpermissive -c ../src/randomkit/rk_primitive.c

rk_sobol.o: ../src/randomkit/rk_sobol.c
	g++ -fpermissive -c ../src/randomkit/rk_sobol.c

cmontecarlo.o: ../src/cmontecarlo.c
	g++ -c ../src/cmontecarlo.c -Irandomkit

test_cmontecarlo.o: test_cmontecarlo.c
	g++ -c test_cmontecarlo.c -Irandomkit

test_cmontecarlo: rk_isaac.o rk_mt.o rk_primitive.o rk_sobol.o cmontecarlo.o test_cmontecarlo.o
	g++ -o test_cmontecarlo rk_isaac.o rk_mt.o rk_primitive.o rk_sobol.o cmontecarlo.o test_cmontecarlo.o
