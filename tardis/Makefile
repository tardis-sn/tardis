all: cmontecarlo_test

rk_isaac.o: randomkit/rk_isaac.c
	g++ -fpermissive -c randomkit/rk_isaac.c

rk_mt.o: randomkit/rk_mt.c
	g++ -fpermissive -c randomkit/rk_mt.c

rk_primitive.o: randomkit/rk_primitive.c
	g++ -fpermissive -c randomkit/rk_primitive.c

rk_sobol.o: randomkit/rk_sobol.c
	g++ -fpermissive -c randomkit/rk_sobol.c

cmontecarlo.o: cmontecarlo.c
	g++ -c cmontecarlo.c -Irandomkit

cmontecarlo_test.o: cmontecarlo_test.c
	g++ -c cmontecarlo_test.c -Irandomkit

CuTest.o: CuTest.c
	g++ -c CuTest.c

cmontecarlo_test: rk_isaac.o rk_mt.o rk_primitive.o rk_sobol.o cmontecarlo.o CuTest.o cmontecarlo_test.o
	g++ -o cmontecarlo_test rk_isaac.o rk_mt.o rk_primitive.o rk_sobol.o cmontecarlo.o CuTest.o cmontecarlo_test.o
