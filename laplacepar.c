#define CONV 1e-6
#define INI 0.2 
#define TAG_DOWN 666
#define TAG_UP 999
#define _USE_MATH_DEFINES // for pi 
#include <stdio.h>      /* printf */
#include <math.h>       /* sin */
#include "mpi.h"
#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))// define the max function
#endif 

int main(int argc, char **argv) {
	FILE *out = fopen("serial.dat", "a"); // write only 
	double **T, **Told;
	double  *y;
	int   i,j, iter,rank,size;
	int down,top;
	int rowt, rowd;
	int NROWS, NCOLS;
	double time1[8],time2[8];
	double error,maxerror; // value for store the max of the error 
	NROWS=atoi(argv[1]);
	NCOLS=atoi(argv[1]);
		T = (double *) malloc((NROWS+1)*sizeof(double));
	if(T == NULL)
		printf(stderr, "out of memory\n");

	for(i = 0; i < NROWS+1; i++)
	{
		T[i] = (double *)malloc((NCOLS+1)* sizeof(double));
		if(T[i] == NULL)
		{
			printf(stderr, "out of memory\n");
			return 1;
		}
	}

	Told =(double *) malloc((NROWS+1)*sizeof(double));
	if(Told == NULL)
		fprintf(stderr, "out of memory\n");

	for(i = 0; i < NROWS+1; i++)
	{
		Told[i] =(double *) malloc((NCOLS+1)* sizeof(double));
		if(Told[i] == NULL)
		{
			fprintf(stderr, "out of memory\n");
			return 1;
		}
	}

	y =(double*) malloc((NROWS+1)*sizeof(double));
	if(Told == NULL)
		fprintf(stderr, "out of memory\n");
	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	down = rank-1; if ( rank == 0 ) down=MPI_PROC_NULL;
	top = rank+1; if (rank == size-1 ) top=MPI_PROC_NULL;
	rowd=(int)(rank*NROWS/size);
	rowt=(int)(rank+1)*NROWS/size;
	time1[rank] = MPI_Wtime();
iter=0;
maxerror=0.0;
	/*    Initial and Boundary Values     */
	if(rank==0){
		for(i=0;i<(rowt+1);i++){
			y[i]=(double)i/(NROWS);
		}
		for( i=1; i<rowt+1; i++ ){
			for ( j=1; j<NCOLS; j++ ){
				T[i][j] = INI;
				T[i][NCOLS]=0.0;
			}
			T[i][0]= (double)(sin(M_PI*y[i])*sin(M_PI*y[i])); //left boundarie
		}
		for(j=0;j<NCOLS+1;j++){
			T[0][j]=0.0;
		}
		for( i=0; i<rowt+1; i++ )
			for( j=0; j<NCOLS+1; j++ ){
				Told[i][j] = T[i][j];
			}

			/*    Do Computation on Sub-grid for Niter iterations     */

			do{ 
error=0.0;
				/* exchange stripe with down neighbour */
				for(j=1; j<NCOLS; j+=2){
					MPI_Sendrecv(&T[rowt][j], 1, MPI_DOUBLE, top,TAG_UP,&T[rowt+1][j], 1, MPI_DOUBLE, top,TAG_UP,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				for( i=1; i<rowt+1; i++ )
					for( j=1+(i)%2; j<NCOLS; j+=2){
						T[i][j] = 0.25 * ( T[i+1][j] + T[i-1][j] +
							T[i][j+1] + T[i][j-1] );
						error=max(error,fabs((double)(Told[i][j]-T[i][j])));
						Told[i][j] = T[i][j];
					}

				for(j=0; j<NCOLS; j+=2){
						MPI_Sendrecv(&T[rowt][j], 1, MPI_DOUBLE, top,TAG_UP,&T[rowt+1][j], 1, MPI_DOUBLE, top,TAG_UP,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
					for( i=1; i<rowt+1; i++ ){
						for( j=1+(i+1)%2; j<NCOLS; j+=2){
							T[i][j] = 0.25 * ( T[i+1][j] + T[i-1][j] +
								T[i][j+1] + T[i][j-1] );
							error=max(error,fabs((double)(Told[i][j]-T[i][j])));
							Told[i][j] = T[i][j];
						}}
					MPI_Allreduce(&error, &maxerror, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
					iter++;
			}while(maxerror>CONV&&iter<5000);/* End of iteration */
  for( i=0; i<rowt+1; i++ )
                for( j=0; j<NCOLS+1; j++ ){
                        Told[i][j] = T[i][j];
                        printf("%d %d %f \n",i,j,T[i][j]);
                }

	}
	else if (rank==size -1){
	for(i=rowd+1;i<NROWS+1;i++){
		y[i]=(double)i/(NROWS);
	}
	for( i=rowd+1; i<rowt; i++ ){
		for ( j=1; j<NCOLS; j++ ){
			T[i][j] = INI;
			T[NROWS][j]=0.0;
		}
		T[i][0]= (double)(sin(M_PI*y[i])*sin(M_PI*y[i])); //left boundarie
		T[i][NCOLS]=0.0;
	
}
T[NROWS][NCOLS]=0.0;
T[NROWS][0]= (double)(sin(M_PI*y[i])*sin(M_PI*y[i]));
 //left boundarie

        for( i=rowd+1; i<rowt+1; i++ )
                for( j=0; j<NCOLS+1; j++ ){
                        Told[i][j] = T[i][j];
                }


		/*    Do Computation on Sub-grid for Niter iterations     */

		do {
error=0.0;
			/* exchange stripe with down neighbour */
			for(j=1; j<NCOLS; j+=2){
				MPI_Sendrecv(&T[rowd+1][j], 1, MPI_DOUBLE, down,TAG_UP,&T[rowd][j], 1 , MPI_DOUBLE, down, TAG_UP,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			for( i=rowd+1; i<NROWS; i++ )
				for( j=1+(i)%2; j<NCOLS; j+=2){
					T[i][j] = 0.25 * ( T[i+1][j] + T[i-1][j] +
						T[i][j+1] + T[i][j-1] );
					error=max(error,fabs((double)(Told[i][j]-T[i][j])));
					Told[i][j] = T[i][j];
				}

				for(j=0; j<NCOLS; j+=2){
					MPI_Sendrecv(&T[rowd+1][j], 1, MPI_DOUBLE, down,TAG_UP,&T[rowd][j], 1 , MPI_DOUBLE, down, TAG_UP,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}
				for( i=rowd+1; i<NROWS; i++ ){
					for( j=1+(i+1)%2; j<NCOLS; j+=2){
						T[i][j] = 0.25 * ( T[i+1][j] + T[i-1][j] +
							T[i][j+1] + T[i][j-1] );
						error=max(error,fabs((double)(Told[i][j]-T[i][j])));
						Told[i][j] = T[i][j];
					}}
				MPI_Allreduce(&error, &maxerror, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				iter++;
		}while(maxerror>CONV&&iter<5000);/* End of iteration */
			
		for( i=rowd+1; i<NROWS+1; i++ )
		for( j=0; j<NCOLS+1; j++ ){
			Told[i][j] = T[i][j];
			printf("%d %d %f \n",i,j,T[i][j]);
		}
	}
	else{
	for(i=rowd+1;i<(rowt+1);i++){
		y[i]=(double)i/(NROWS);
	}
	for( i=rowd+1; i<rowt+1; i++ ){
		for ( j=1; j<NCOLS; j++ ){
			T[i][j] = INI;
		}
		T[i][0]= (double)(sin(M_PI*y[i])*sin(M_PI*y[i])); //left boundarie
		T[i][NCOLS]=0.0;
	}
 for( i=rowd+1; i<rowt+1; i++ )
                for( j=0; j<NCOLS+1; j++ ){
                        Told[i][j] = T[i][j];
                }

		/*    Do Computation on Sub-grid for Niter iterations     */

		do {
error=0.0;
			/* exchange stripe with down neighbour */
			for(j=1; j<NCOLS; j+=2){
				MPI_Sendrecv(&T[rowd+1][j], 1, MPI_DOUBLE, down,TAG_UP,&T[rowd][j], 1 , MPI_DOUBLE, down, TAG_UP,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Sendrecv(&T[rowt][j], 1, MPI_DOUBLE, top,TAG_UP,&T[rowt+1][j], 1, MPI_DOUBLE, top,TAG_UP,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			for( i=rowd+1; i<rowt+1; i++ )
				for( j=1+(i)%2; j<NCOLS; j+=2){
					T[i][j] = 0.25 * ( T[i+1][j] + T[i-1][j] +
						T[i][j+1] + T[i][j-1] );
					error=max(error,fabs((double)(Told[i][j]-T[i][j])));
					Told[i][j] = T[i][j];

				}

				for(j=0; j<NCOLS; j+=2){
					MPI_Sendrecv(&T[rowd+1][j], 1, MPI_DOUBLE, down,TAG_UP,&T[rowd][j], 1 , MPI_DOUBLE, down, TAG_UP,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Sendrecv(&T[rowt][j], 1, MPI_DOUBLE, top,TAG_UP,&T[rowt+1][j], 1, MPI_DOUBLE, top,TAG_UP,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
				for( i=rowd+1; i<rowt+1; i++ ){
					for( j=1+(i+1)%2; j<NCOLS; j+=2){
						T[i][j] = 0.25 * ( T[i+1][j] + T[i-1][j] +
							T[i][j+1] + T[i][j-1] );
						error=max(error,fabs((double)(Told[i][j]-T[i][j])));
						Told[i][j] = T[i][j];

					}}
				MPI_Allreduce(&error, &maxerror, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
				iter++;

		}while(maxerror>CONV&&iter<5000);/* End of iteration */
		for( i=rowd+1; i<rowt+1; i++ )
			for( j=0; j<NCOLS+1; j++ ){
				Told[i][j] = T[i][j];
				printf("%d %d %f \n",i,j,T[i][j]);
			}
	}
	time2[rank] = MPI_Wtime();
	fprintf(out, " %d %d %d %lf",size,NROWS,iter,time2[rank]-time1[rank]);
	//fprintf(out,"Time elapsed for processor %d: %lf  iter %d \n", rank, time2[rank]-time1[rank],iter);
	fclose(out);	
MPI_Finalize();
	return 0 ;
}    /* End of Program */
