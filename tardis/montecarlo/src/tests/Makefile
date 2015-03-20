all: test_cmontecarlo

rk_isaac.o: ../randomkit/rk_isaac.c
	g++ -fpermissive -c ../randomkit/rk_isaac.c

rk_mt.o: ../randomkit/rk_mt.c
	g++ -fpermissive -c ../randomkit/rk_mt.c

rk_primitive.o: ../randomkit/rk_primitive.c
	g++ -fpermissive -c ../randomkit/rk_primitive.c

rk_sobol.o: ../randomkit/rk_sobol.c
	g++ -fpermissive -c ../randomkit/rk_sobol.c

cmontecarlo.o: ../cmontecarlo.c
	g++ -c ../cmontecarlo.c -Irandomkit

test_cmontecarlo.o: test_cmontecarlo.c
	g++ -c test_cmontecarlo.c -Irandomkit

test_cmontecarlo: rk_isaac.o rk_mt.o rk_primitive.o rk_sobol.o cmontecarlo.o test_cmontecarlo.o
	g++ -o test_cmontecarlo rk_isaac.o rk_mt.o rk_primitive.o rk_sobol.o cmontecarlo.o test_cmontecarlo.o
