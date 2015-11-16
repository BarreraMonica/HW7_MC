#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#define PI 3.1416
#define USAGE "./mcmc_solar.x n_steps n_burn"

void likelihood(float *chi, float *y_obs, float *y_model, int n_puntos);
void my_model(float *mine, float *anio ,float a, float b, float c, float d, int n_puntos);
void mcmc(int n_steps, float *a, float *b, float *c, float *d, float *l, float *x_obs, float *y_obs,int n_puntos, float *mine);
void min(float *mine, float *x_obs, float *a, float *b, float *c, float *d, float *l, int n_steps, int n_puntos);

void leer(float *anio, float *manc, int n_puntos);
float *reserva(int n_puntos);
void print_array(float *a, float *b, float *c, float *d , float *l, int n_steps, int n_burn);

int main(int argc, char **argv){
	
	srand48(time(NULL));
	
	//Valores que entran
	int n_steps=atof(argv[1]);
	int n_burn=atof(argv[2]);
	//Inicializacion de los arrays que guardan los datos
	int n_puntos= 4141;
	float *manc;
	float *anio;
	anio=reserva(n_puntos);
	manc=reserva(n_puntos);
	
	leer(anio, manc, n_puntos);
	
	//Array con los resultados del ajuste
	float *mine;
	mine=reserva(n_puntos);
	
	//Arrays que guardan los datos de las variables
	float *a;
	float *b;
	float *c;
	float *d;
	float *l;
	a=reserva(n_steps);
	b=reserva(n_steps);
	c=reserva(n_steps);
	d=reserva(n_steps);
	l=reserva(n_steps);
	//Primer paso en la cadena
	a[0]=drand48();
	b[0]=drand48();
	c[0]=drand48();
	d[0]=drand48();
	float chi;
	my_model(mine,anio,a[0],b[0],c[0],d[0],n_puntos);
	likelihood(&chi,manc,mine,n_puntos);
	l[0]=chi;
	
	//Iteraciones
	mcmc(n_steps,a,b,c,d,l,anio,manc,n_puntos,mine);
	min(mine,anio,a,b,c,d,l,n_steps,n_puntos);
	
	print_array(a,b,c,d,l,n_steps,n_burn);
	
	free(mine);
	free(anio);
	free(manc);
	free(a);
	free(b);
	free(c);
	free(d);
	
	return 0;
}

//OBTENCION DEL MEJOR AJUSTE
void min(float *mine, float *x_obs, float *a, float *b, float *c, float *d, float *l, int n_steps, int n_puntos){
	int i;
	float temp, a2, b2, c2, d2;
	temp=l[0];
	a2=a[0];
	b2=b[0];
	c2=c[0];
	d2=d[0];
	for(i=1;i<n_steps;i++){
		if(temp > l[i] ){
			temp=l[i];
			a2=a[i];
			b2=b[i];
			c2=c[i];
			d2=d[i];
		}
	}
	my_model(mine,x_obs,a2,b2,c2,d2,n_puntos);
}

//ITERACIONES
void mcmc(int n_steps, float *a, float *b, float *c, float *d, float *l, float *x_obs, float *y_obs,int n_puntos, float *mine){
	int i;
	float a2,b2,c2,d2, chi, chi2, alpha, beta;
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	
    for(i=0;i<(n_steps-1);i++){
		a2 = a[i] + gsl_ran_gaussian(r, 0.1);
		b2 = b[i] + gsl_ran_gaussian(r, 0.1);
		c2 = c[i] + gsl_ran_gaussian(r, 0.1);
		d2 = c[i] + gsl_ran_gaussian(r, 0.1);
		
		float *y_init;
		float *y_prime;
		y_init=reserva(n_puntos);
		y_prime=reserva(n_puntos);
		
		my_model(y_init,x_obs,a[i],b[i],c[i],d[i],n_puntos);
		my_model(y_prime,x_obs,a2,b2,c2,d2,n_puntos);
	
		likelihood(&chi,y_obs,y_init,n_puntos);
		likelihood(&chi2,y_obs,y_prime,n_puntos);
		
		alpha = -chi2 + chi;
		if(alpha >= 0.0){
			a[i+1]=a2;
			b[i+1]=b2;
			c[i+1]=c2;
			d[i+1]=d2;
			l[i+1]=chi2;
		}
		else{
			beta = drand48();
			if(alpha < -6000.0 ){
				a[i+1]=a[i];
				b[i+1]=b[i];
				c[i+1]=c[i];
				d[i+1]=d[i];
				l[i+1]=chi;
			}
			else if(beta <= exp(alpha)){
				a[i+1]=a2;
				b[i+1]=b2;
				c[i+1]=c2;
				d[i+1]=d2;
				l[i+1]=chi2;
			}
			else{
				a[i+1]=a[i];
				b[i+1]=b[i];
				c[i+1]=c[i];
				d[i+1]=d[i];
				l[i+1]=chi;
			}
		}
		
		free(y_init);
		free(y_prime);
	}

  gsl_rng_free (r);
}

//LECTURA DEL ARCHIVO
void leer(float *anio, float *manc, int n_puntos){
	
	float *mes;
	mes=reserva(n_puntos);
	
	FILE *in;
	int i;
	in=fopen("monthrg.dat","r");
	if(!in){
    printf("problems opening the file %s\n","monthrg.dat");
    exit(1);
	}
	float temp1;
	float temp2;
	for(i=0;i<n_puntos;i++){
	  do{
		  fscanf(in, "%f %f %f %f %f\n", &anio[i], &mes[i], &temp1, &manc[i] , &temp2);
	  }while(manc[i] == -99.0);
	}
	fclose(in);
	
	//Calculo del anio como decimal
	for(i=0;i<n_puntos;i++){
		anio[i] = anio[i] + (mes[i]-1.0)/12.0;
	}
	free(mes);
}

//IMPRESION DE LOS ARRAYS
void print_array(float *a, float *b, float *c, float *d , float *l, int n_steps, int n_burn){
  int i;
  for(i=n_burn;i<n_steps;i++){
    printf("%f %f %f %f %f\n", a[i], b[i], c[i] , d[i], l[i]);
  }
}

//EVALUAR LOS PUNTOS CON EL AJUSTE
void my_model(float *mine, float *anio ,float a, float b, float c, float d, int n_puntos){
	int i;
	for(i=0;i<n_puntos;i++){
		mine[i] = a * cos((2*PI/d)*anio[i] + b) + c;
	}
}

//RESERVA DE MEMORIA
float *reserva(int n_puntos){
  float *array;
  int i;
  if(!(array = malloc(n_puntos * sizeof(float)))){
    printf("Problema en reserva\n");
    exit(1);
  }
  for(i=0;i<n_puntos;i++){
    array[i] = 0.0;
  }
  return array;
}

//CALCULO DE LA VARIACION DEL AJUSTE
void likelihood(float *chi, float *y_obs, float *y_model, int n_puntos){
	int i;
	float temp=0;
	for(i=0;i<n_puntos;i++){
		temp=temp+pow(y_obs[i]-y_model[i],2.0);
	}
	*chi=(1.0/2.0)*temp;
}


