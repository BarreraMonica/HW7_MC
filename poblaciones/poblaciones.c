#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#define USAGE "./mcmc_loktavolterra.x n_steps n_burn"

void rg4th(float *resultado, float x, float y, float a, float b, float c, float d, float h);
void x_prime(float *valor, float x, float y, float alpha, float beta);
void y_prime(float *valor, float x, float y, float gamma, float delta);

void my_model(float *mine_presas,float *mine_depredadores,float *tiempo,float ini_presas, float ini_depredadores,float a, float b, float c, float d, int n_puntos);
void likelihood(float *chi,float *presas,float *depredadores,float *mine_presas,float *mine_depredadores,int n_puntos);
void mcmc(int n_steps,float *a,float *b,float *c,float *d,float *l,float *tiempo,float *presas,float *depredadores,int n_puntos,float *mine_presas,float *mine_depredadores, float ini_presas, float ini_depredadores);
void min(float *mine_presas,float *mine_depredadores,float ini_presas,float ini_depredadores,float *tiempo,float *a,float *b,float *c,float *d,float *l,int n_steps,int n_puntos);

void leer(float *tiempo, float *presas, float *depredadores, int n_puntos);
float *reserva(int n_puntos);
void print_array(float *a, float *b, float *c, float *d , float *l, int n_steps, int n_burn);

int main(int argc, char **argv){
	
	srand48(time(NULL));
	
	//Valores que entran
	int n_steps=atof(argv[1]);
	int n_burn=atof(argv[2]);
	//Inicializacion de los arrays que guardan los datos
	int n_puntos= 96;
	float *tiempo;
	float *presas;
	float *depredadores;
	tiempo=reserva(n_puntos);
	presas=reserva(n_puntos);
	depredadores=reserva(n_puntos);
	
	leer(tiempo, presas, depredadores, n_puntos);
	//Condiciones iniciales para resolver las ecuaciones diferenciales
	float ini_presas, ini_depredadores;
	ini_presas=presas[0];
	ini_depredadores=depredadores[0];
	
	//Array con los resultados del ajuste
	float *mine_presas;
	mine_presas=reserva(n_puntos);
	float *mine_depredadores;
	mine_depredadores=reserva(n_puntos);
	
	//Arrays que guardan los datos de las variables
	float *a; //alpha
	float *b; //beta
	float *c; //gamma
	float *d; //delta
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
	my_model(mine_presas,mine_depredadores,tiempo,ini_presas,ini_depredadores,a[0],b[0],c[0],d[0],n_puntos);
	likelihood(&chi,presas,depredadores,mine_presas,mine_depredadores,n_puntos);
	l[0]=chi;
	
	mcmc(n_steps,a,b,c,d,l,tiempo,presas,depredadores,n_puntos,mine_presas,mine_depredadores,ini_presas,ini_depredadores);
	min(mine_presas,mine_depredadores,ini_presas,ini_depredadores,tiempo,a,b,c,d,l,n_steps,n_puntos);

	print_array(a,b,c,d,l,n_steps,n_burn);
	
	free(tiempo);
	free(presas);
	free(mine_presas);
	free(depredadores);
	free(mine_depredadores);
	free(a);
	free(b);
	free(c);
	free(d);
	
	return 0;
}

//OBTENCION DEL MEJOR AJUSTE
void min(float *mine_presas,float *mine_depredadores,float ini_presas,float ini_depredadores,float *tiempo,float *a,float *b,float *c,float *d,float *l,int n_steps,int n_puntos){
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
	my_model(mine_presas,mine_depredadores,tiempo,ini_presas,ini_depredadores,a2,b2,c2,d2,n_puntos);	
}

//ITERACIONES
void mcmc(int n_steps,float *a,float *b,float *c,float *d,float *l,float *tiempo,float *presas,float *depredadores,int n_puntos,float *mine_presas,float *mine_depredadores, float ini_presas, float ini_depredadores){
	int i;
	float a2,b2,c2,d2,chi,chi2,alpha,beta;
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
		
		float *presas_init;
		float *presas_prime;
		presas_init=reserva(n_puntos);
		presas_prime=reserva(n_puntos);
		float *depredadores_init;
		float *depredadores_prime;
		depredadores_init=reserva(n_puntos);
		depredadores_prime=reserva(n_puntos);
		
		my_model(presas_init,depredadores_init,tiempo,ini_presas,ini_depredadores,a[i],b[i],c[i],d[i],n_puntos);
		likelihood(&chi,presas,depredadores,presas_init,depredadores_init,n_puntos);
		
		my_model(presas_prime,depredadores_prime,tiempo,ini_presas,ini_depredadores,a2,b2,c2,d2,n_puntos);
		likelihood(&chi2,presas,depredadores,presas_prime,depredadores_prime,n_puntos);
		
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
		
		free(presas_init);
		free(presas_prime);
		free(depredadores_init);
		free(depredadores_prime);
	}

  gsl_rng_free (r);
}

//IMPRESION DE LOS ARRAYS
void print_array(float *a, float *b, float *c, float *d , float *l, int n_steps, int n_burn){
  int i;
  for(i=n_burn;i<n_steps;i++){
    printf("%f %f %f %f %f\n", a[i], b[i], c[i] , d[i], l[i]);
  }
}

//DEFINICION DEL CALCULO DEL AJUSTE
void my_model(float *mine_presas,float *mine_depredadores,float *tiempo,float ini_presas, float ini_depredadores,float a, float b, float c, float d, int n_puntos){	
	
	//Condiciones iniciales
	mine_presas[0] = ini_presas;
	mine_depredadores[0] = ini_depredadores;  
	
	//Recorrido R-K
	int i;
	for(i=1;i<n_puntos;i++){
		float *r;
		float h;
		h = tiempo[i]-tiempo[i-1];
		r=reserva(2);
		rg4th(r, mine_presas[i-1], mine_depredadores[i-1], a, b, c, d, h);
		mine_presas[i]=r[0];
		mine_depredadores[i]=r[1];
		free(r);
	}
}

//RUNGE-KUTTA DE CUARTO ORDEN
void rg4th(float *resultado, float x, float y, float a, float b, float c, float d, float h){
    float k1_x,k1_y;
    x_prime(&k1_x,x,y,a,b);
    y_prime(&k1_y,x,y,c,d);
    
    //First step
    float x1,y1;
    x1 = x + (h/2.0) * k1_x;
    y1 = y + (h/2.0) * k1_y;
	
	float k2_x,k2_y;
    x_prime(&k2_x,x1,y1,a,b);
    y_prime(&k2_y,x1,y1,c,d);
    
    //Second step
    float x2,y2;
    x2 = x + (h/2.0) * k2_x;
    y2 = y + (h/2.0) * k2_y;
	
	float k3_x,k3_y;
	x_prime(&k3_x,x2,y2,a,b);
    y_prime(&k3_y,x2,y2,c,d);
     
    //Third step
    float x3,y3;
    x3 = x + h * k3_x;
    y3 = y + h * k3_y;
	
	float k4_x,k4_y;
	x_prime(&k4_x,x3,y3,a,b);
    y_prime(&k4_y,x3,y3,c,d);
     
    //Fourth step
    float average_x, average_y;
    average_x = (1.0/6.0)*(k1_x + 2.0*k2_x + 2.0*k3_x + k4_x);
    average_y = (1.0/6.0)*(k1_y + 2.0*k2_y + 2.0*k3_y + k4_y);
    
    resultado[0] = x + h * average_x;
    resultado[1] = y + h * average_y;
}

//ECUACION DIFERENCIAL DE PRESAS
void x_prime(float *valor, float x, float y, float alpha, float beta){
	*valor= x * (alpha - beta * y);
}

//ECUACION DIFERENCIAL DE DEPREDADORES
void y_prime(float *valor, float x, float y, float gamma, float delta){
	*valor = -1 * y * (gamma - delta * x);
}


//RESERVA DEL AJUSTE
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
void likelihood(float *chi,float *presas,float *depredadores,float *mine_presas,float *mine_depredadores,int n_puntos){
	int i;
	float temp1=0;
	float temp2=0;
	for(i=0;i<n_puntos;i++){
		temp1=temp1+pow(presas[i]-mine_presas[i],2.0);
		temp2=temp2+pow(depredadores[i]-mine_depredadores[i],2.0);
	}
	*chi=(1.0/2.0)*(temp1+temp2);
}

//LECTURA DEL ARCHIVO
void leer(float *tiempo, float *presas, float *depredadores,int n_puntos){
	
	FILE *in;
	int i;
	in=fopen("archivo.dat","r");
	if(!in){
    printf("problems opening the file %s\n","archivo.dat");
    exit(1);
	}
	for(i=0;i<(n_puntos);i++){

		  fscanf(in, "%f %f %f\n", &tiempo[i], &presas[i], &depredadores[i]);
		
	}
	fclose(in);
}


