n_steps=20000
n_burn=1000

all: Probabilidades_solar.pdf

Probabilidades_solar.pdf: solar.py solar.dat
	python solar.py

solar.dat: mcmc_solar.x 
	./mcmc_solar.x $(n_steps) $(n_burn) > solar.dat

mcmc_solar.x : solar.c
	cc solar.c  -lgsl -lgslcblas -lm -o mcmc_solar.x

clean:
	rm -f mcmc_solar.x solar.dat Probabilidades_solar.pdf
