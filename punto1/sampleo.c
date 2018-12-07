#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <omp.h>

#define N 1000
//int i;


double f(double x)
{
 return pow(2*3.1415926535,-1/2)*exp(-x*x/2.0);   
}

void print(int rank, double *cadenaAImprimir)
{
    char *satan = "s";
    sprintf(satan, "gauss%d.txt", rank);
	FILE *output = fopen(satan, "w+");
	for(int j=0;j<N;j+=1) {
            fprintf(output, "%f\n",cadenaAImprimir[j]);
			}
	fclose(output);
}

int main(int argc, char ** argv)
{
    
double *cadena;
double candidato;
double r;
int thread_id;
double paso = 0.3;
#pragma omp parallel //private(cadena,r,i,thread_id) //private(cadena, candidato,r,i,thread_id)
{
 int thread_id= omp_get_thread_num();
int thread_count = omp_get_num_threads();
srand48(thread_id);




double *cadena = (double*) calloc(N,sizeof(double));

cadena[0] = paso*(drand48()-0.5);
//printf("sirve\n");
//#pragma omp for
int i = 0;
for(i=1;i<N;i+=1){
 candidato = cadena[i-1]+(paso)*(drand48()-0.5);
 r = f(candidato)/f(cadena[i-1]);
 r = fmin(1,r);
 if(drand48() < r){
  cadena[i] = candidato;   
 }
 else{
 cadena[i] = cadena[i-1];   
 }
}

//print(thread_id, cadena);
char *satan = (char *) malloc(sizeof(char)*1000);
thread_id= omp_get_thread_num();
sprintf(satan, "gauss%d.txt", thread_id);
FILE *output = fopen(satan, "w+");
for(int j=0;j<N;j+=1) {
          fprintf(output, "%f\n",cadena[j]);
			}
	fclose(output);

}
free(cadena);
return 0;
}




