n_steps=20000
n_burn=1000

all: Probabilidades_poblaciones.pdf

Probabilidades_poblaciones.pdf: poblaciones.dat poblaciones.py
	python poblaciones.py

poblaciones.dat: mcmc_loktavolterra.x archivo.dat 
	./mcmc_loktavolterra.x $(n_steps) $(n_burn) > poblaciones.dat

archivo.dat: lotka_volterra_obs.dat
	tail -n +2 lotka_volterra_obs.dat > archivo.dat

mcmc_loktavolterra.x : poblaciones.c
	cc poblaciones.c  -lgsl -lgslcblas -lm -o mcmc_loktavolterra.x

clean:
	rm -f mcmc_loktavolterra.x poblaciones.dat Probabilidades_poblaciones.pdf archivo.dat

