#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h> 
#include <omp.h>

#define G 667408E-16
// CONSTANTE DE ELASTICIDADE
#define E 0.5     
typedef struct galaxy{
   double x;
   double y;
   double z;
   double m;
   double ax;
   double ay; 
   double az;
   double oldvx;
   double oldvy;
   double oldvz;
   double fx;
   double fy; 
   double fz;
   double vx;
   double vy;
   double vz; 
   double newx;
   double newy;
   double newz;

} galaxy;
int iteracao, count, count2;
double positive(double in){
	if(in<0)
		return(-in);
	else
		return(in);
}
double normaV(float x, float y, float z){
	double aux;
	aux =  pow(x,2);
	aux += pow(y,2);
	aux += pow(z,2);
	aux = sqrt(aux);
	return(aux);
}

void Fijstruc ( galaxy *planeta, int size, int inctempo){

	
	
	FILE * fp;
	fp = fopen ("particles.dat", "a");
	fprintf(fp,"\n");
	int i ,j;
	//printf("A comecar a simulacao com %i processadores. Hello P%i.\n", omp_get_num_threads(), omp_get_thread_num());
	
		#pragma omp parallel for num_threads(omp_get_max_threads()) private(i,j) shared(planeta)
		for (i=0; i<size; i++){
		
		planeta[i].fx=0;	//reset vetor Forca
		planeta[i].fy=0;
		planeta[i].fz=0;
		//printf("A comecar a simulacao com %i processadores. Hello P%i.\n", omp_get_num_threads(), omp_get_thread_num());
			for(j=0; j<size; j++){
			if(j!=i){
			double Px, Py, Pz;
			double GMM_norm3;
			double Modulo;
			Px = planeta[j].x-planeta[i].x;					//Pjx - Pix
			Py = planeta[j].y-planeta[i].y;					//Pjy - Piy
			Pz = planeta[j].z-planeta[i].z;					//Pjz - Piz

			Modulo = normaV(Px,Py,Pz);					// ||Pj-Pi||
			GMM_norm3 = planeta[i].m * planeta[j].m * G / pow(Modulo,3);   // G * mi *mj / ||Pj-Pi||^3     6.87 * (1/pow(10, 11))

			//printf("\nPlaneta %d Planeta %d Pjx - Pix: %f   Pjy -Piy: %f  Pjz-Piz: %f\n",i, j, Px, Py, Pz);
			//printf("MODULO: %f de ||Pj - Pi||\n", Modulo);
			
			//#pragma omp critical
			//{	
				//printf("%f PRocessor :%d\n",planeta[i].fx, omp_get_thread_num());
				planeta[i].fx += Px * GMM_norm3;				//acumulador da componente Forca Resultante de i componente x
				//printf("%lf * %lf =%lf  PRocessor :%d\n\n", Px, GMM_norm3, planeta[i].fx, omp_get_thread_num());
				planeta[i].fy += Py * GMM_norm3;				// para y
				planeta[i].fz += Pz * GMM_norm3;				// para z
				
			//}
				}

			}
			

		planeta[i].ax = planeta[i].fx / planeta[i].m;				// componente x de acelaracao do Corpo i
		//printf("P:%d A:%lf = %lf / %lf\n",i, planeta[i].ax, planeta[i].fx, planeta[i].m);
		planeta[i].ay = planeta[i].fy / planeta[i].m;				//ay i
		planeta[i].az = planeta[i].fz / planeta[i].m;				//az i

		planeta[i].oldvx = planeta[i].vx;					//update old to new point
		planeta[i].oldvy = planeta[i].vy;
		planeta[i].oldvz = planeta[i].vz;	
	
		planeta[i].vx = planeta[i].oldvx + planeta[i].ax * inctempo;		// Vx(i)=Vx(i-1) + ax * stepTempo
		planeta[i].vy = planeta[i].oldvy + planeta[i].ay * inctempo;		// Vy(i)=Vy(i-1) + ay * stepTempo
		planeta[i].vz = planeta[i].oldvz + planeta[i].az * inctempo;		// Vz(i)=Vz(i-1) + az * stepTempo

		
		//printf( "\n\nPlaneta %d Velocidade X:%lf  Y:%lf  Z:%lf \nVold: X:%lf  Y:%lf  Z:%lf\n\n\n",i, planeta[i].vx,planeta[i].vy,planeta[i].vz,planeta[i].oldvx,planeta[i].oldvy,planeta[i].oldvz);
		planeta[i].newx= planeta[i].x + planeta[i].vx * inctempo;
		planeta[i].newy= planeta[i].y + planeta[i].vy * inctempo;
		planeta[i].newz= planeta[i].z + planeta[i].vz * inctempo;
		
		//printf("Planeta %d Iteracao:%d Novas coordenadas: %f %f %f %f\n", i,iteracao, planeta[i].newx, planeta[i].newy, planeta[i].newz, planeta[i].m);
			
		
		}	
	
	#pragma omp parallel for collapse(2) private(i,j) shared(planeta)	
	for (i=0; i<size; i++){
		for(j=size; 0<=j; j--){
			if(j<i){
			double Px = planeta[j].newx-planeta[i].newx;   //calcular
			double Py = planeta[j].newy-planeta[i].newy;
			double Pz = planeta[j].newz-planeta[i].newz;
      				if(normaV(Px,Py,Pz)<100){
					float mi, mj;
					count++;
				
					//if(positive(planeta[i].m-planeta[j].m)<500){
							
						mi=planeta[i].m;
						mj=planeta[j].m;
						if(planeta[i].vx / planeta[j].vx<0)
						planeta[i].vx=-(mi*planeta[i].vx + mj * planeta[j].oldvx + mj * E *(planeta[j].oldvx-planeta[i].oldvx)) / (mi + mj);
						else 
						planeta[i].vx=(mi*planeta[i].vx + mj * planeta[j].oldvx + mj * E *(planeta[j].oldvx-planeta[i].oldvx)) / (mi + mj);

						if(planeta[i].vy / planeta[j].vy<0)
						planeta[i].vy=-(mi*planeta[i].vy + mj * planeta[j].oldvy + mj * E *(planeta[j].oldvy-planeta[i].oldvy)) / (mi + mj);
						else planeta[i].vy=(mi*planeta[i].vy + mj * planeta[j].oldvy + mj * E *(planeta[j].oldvy-planeta[i].oldvy)) / (mi + mj);
						
						if(planeta[i].vy / planeta[j].vy<0)
						planeta[i].vz=-(mi*planeta[i].vz + mj * planeta[j].oldvz + mj * E *(planeta[j].oldvz-planeta[i].oldvz)) / (mi + mj);
						else
						planeta[i].vz=(mi*planeta[i].vz + mj * planeta[j].oldvz + mj * E *(planeta[j].oldvz-planeta[i].oldvz)) / (mi + mj);

						if(planeta[i].vx / planeta[j].vx<0)
						planeta[j].vx=-(mi*planeta[i].vx + mj * planeta[j].oldvx + mi * E *(planeta[i].oldvx-planeta[j].oldvx)) / (mi + mj);
						else
						planeta[j].vx=-(mi*planeta[i].vx + mj * planeta[j].oldvx + mi * E *(planeta[i].oldvx-planeta[j].oldvx)) / (mi + mj);

						if(planeta[i].vy / planeta[j].vy<0)
						planeta[j].vy=-(mi*planeta[i].vy + mj * planeta[j].oldvy + mi * E *(planeta[i].oldvy-planeta[j].oldvy)) / (mi + mj);
						else
						planeta[j].vy=(mi*planeta[i].vy + mj * planeta[j].oldvy + mi * E *(planeta[i].oldvy-planeta[j].oldvy)) / (mi + mj);

						if(planeta[i].vy / planeta[j].vy<0)
						planeta[j].vz=-(mi*planeta[i].vz + mj * planeta[j].oldvz + mi * E *(planeta[i].oldvz-planeta[j].oldvz)) / (mi + mj);
						else planeta[j].vz=(mi*planeta[i].vz + mj * planeta[j].oldvz + mi * E *(planeta[i].oldvz-planeta[j].oldvz)) / (mi + mj);
						
						//printf( "Planeta %d Velocidade X:%f  Y:%f  Z:%f \nVold: X:%f  Y:%f  Z:%f\n",i, planeta[i].vx,planeta[i].vy,planeta[i].vz,planeta[i].oldvx,planeta[i].oldvy,planeta[i].oldvz);
						planeta[i].newx= planeta[i].x + planeta[i].vx * inctempo;
						planeta[i].newy= planeta[i].y + planeta[i].vy * inctempo;
						planeta[i].newz= planeta[i].z + planeta[i].vz * inctempo;
						planeta[j].newx= planeta[j].x + planeta[j].vx * inctempo;
						planeta[j].newy= planeta[j].y + planeta[j].vy * inctempo;
						planeta[j].newz= planeta[j].z + planeta[j].vz * inctempo;
						//}

					}
				}
			}
		}

		for (i=0; i<size; i++){
			fprintf(fp,"%0.1f %0.1f %0.1f %0.1f \n", planeta[i].newx, planeta[i].newy, planeta[i].newz, planeta[i].m);
			planeta[i].x= planeta[i].newx;
			planeta[i].y= planeta[i].newy;	
			planeta[i].z= planeta[i].newz;
			//planeta[i].m= planeta[i].m/sqrt(1-(pow(norma_V, 2)/pow(300000000, 2)));	//chage mass due to velocity
			//fprintf(fp,"%4.1f %4.1f %4.1f %5.1f\n", planeta[i].x, planeta[i].y, planeta[i].z, planeta[i].m);
		}
		fclose(fp);
	
}			

int main(){

  galaxy *location = (galaxy *)malloc(1000 * sizeof (galaxy));
  int tempototal, inctempo;
  FILE * fp;
  fp = fopen ("particles.dat", "r");
  int size = 0;
  fscanf(fp, "%d %d", &tempototal, &inctempo);
  while(fscanf(fp, "%lf %lf %lf %lf",  &location[size].x, &location[size].y,  &location[size].z , &location[size].m) == 4){
	  printf("Planeta %d : X: %f  Y: %f  Z: %f M: %f\n", size, location[size].x, location[size].y, location[size].z, location[size].m);
	  location[size].newx=0;
	  location[size].newy=0;
	  location[size].newz=0;
	  location[size].vx=0;
	  location[size].vy=0;
	  location[size].vz=0;
	  location[size].oldvx=0;
	  location[size].oldvy=0;
	  location[size].oldvz=0;
      size +=1;
  }

fclose(fp);


struct timespec start_seq, end_seq;
clock_gettime(CLOCK_MONOTONIC, &start_seq);
	for( int t=0; t<=tempototal; t+=inctempo){
		Fijstruc(location, size, inctempo);
		iteracao++;
		printf("Iteracao: %d Countof crashes: %d\n\n\n",iteracao, count);	

	}
 clock_gettime(CLOCK_MONOTONIC, &end_seq);

    double initialTimeSeq=(start_seq.tv_sec*1e3)+(start_seq.tv_nsec*1e-6);
    double finalTimeSeq=(end_seq.tv_sec*1e3)+(end_seq.tv_nsec*1e-6);

    printf("Simulacao: %f s\n", (finalTimeSeq-initialTimeSeq) / 1000.0);
free(location);


return 0;
}
 

			

	
