#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

double bench_t_start, bench_t_end;

static
double rtclock()
{
struct timeval Tp;
int stat;
stat = gettimeofday (&Tp, NULL);
if (stat != 0)
printf ("Error return from gettimeofday: %d", stat);
return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start()
{
bench_t_start = rtclock ();
}

void bench_timer_stop()
{
bench_t_end = rtclock ();
}

void bench_timer_print()
{
printf ("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start);
}

#define MAXSIZE 3000
float Matrix[MAXSIZE][MAXSIZE];

void FillMatrix(float Matrix[MAXSIZE][MAXSIZE] , int start, int stop )
{
	for (int i = 0; i < MAXSIZE; ++i)
	{
		for (int j = 0; j < MAXSIZE; ++j)
		{
			Matrix[i][j]=0;
		}
	}


	for (int i = start; i <= stop; ++i)
	{
		for (int j = 0; j < MAXSIZE; ++j)
		{
			//Matrix[i][j] = i+j+i*j+(i+1)/(j+1) - (i+1)%(j+1);// RESULT TEST [3,4,5,6,7] Matrix
			if (i==j) Matrix[i][j]=1; else Matrix[i][j]=0;      // TIME TEST 
			//Matrix[i][j] = MAXSIZE*i+j+1;					 // Another Result test
			//Matrix[i][j]=0;							    // Funny TEst					
		}
	}
}

void PrintMatrix( float Matrix[][MAXSIZE] , int start, int stop )
{
	for (int i = start; i <=stop; ++i)
	{
		printf("\n");
		printf("%d ::: ",i );
		for (int j = 0; j < MAXSIZE; ++j)
		{
			printf("%f    ",Matrix[i][j]);
		}
	}
	printf("\n");
}

int Sort(int minus[MAXSIZE])
{
	int count=0;
	for (int i = 0; i < MAXSIZE; ++i)
	{
		for (int j = 0; j < MAXSIZE-1; ++j)
		{
			if (minus[j] > minus[j+1]) {int a = minus[j];minus[j]=minus[j+1];minus[j+1]=a;count++;}
		}
	}
	return count;
}

int Round(int i, int rank,int count )
{
	int key=rank;
	while( key<i)
	{
		key+=count;
	}
	return key;	
}

int main ()
{
	MPI_Init(NULL, NULL);
	int rank;
    int count;
    float RESULT=1;
    MPI_Status stat;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &count);
    int minus[MAXSIZE]={0};

    bench_timer_start();
   	
	FillMatrix( Matrix,0,MAXSIZE-1);

	for (int i = 0; i < MAXSIZE; ++i)
	{
		// Congig OUR matrix
		if( rank == i%count)
		{
			for (int k = 0; k < count; ++k)
			{
				if(k != i%count)
				{	
					MPI_Send(&Matrix[i],MAXSIZE,MPI_FLOAT,k,0,MPI_COMM_WORLD);
				}	
			}  // можно попытаться сделать неблокирующую но проблема будет при синхронизации
				// хотя можно и попробовать будет интересно
			int   maxj=0;
			float maxel = 0;
			float cand;
			for (int j = 0; j < MAXSIZE; j++)
			{
				cand = Matrix[i][j];
				if(Matrix[i][j] < 0) cand=cand*(-1);
				if(cand > maxel) {maxel=cand;maxj=j;}
			}
			RESULT=RESULT*Matrix[i][maxj];
			//printf("I AM %d ; AND I MULTIPLAY %f ;\n",rank,Matrix[i][maxj] );
			minus[maxj]=i+1;

			float kf=0;
			for (int k = i+count; k < MAXSIZE; k = k + count)
			{
				kf=Matrix[k][maxj]/Matrix[i][maxj];
				for (int j = 0; j < MAXSIZE; ++j)
				{
					Matrix[k][j] = Matrix[k][j] - (Matrix[i][j]*kf); 
				}
			}
		}
		else
		{

			MPI_Recv(&Matrix[i],MAXSIZE,MPI_FLOAT,i%count,0,MPI_COMM_WORLD,&stat);
			
			int   maxj=0;
			float maxel = 0;
			float cand;
			for (int j = 0; j < MAXSIZE; j++)
			{
				cand = Matrix[i][j];
				if(Matrix[i][j] < 0) cand=cand*(-1);
				if(cand > maxel) {maxel=cand;maxj=j;}
			}

			float kf=0;	
			for (int k = Round(i,rank,count); k < MAXSIZE; k=k+count)
			{
				kf=Matrix[k][maxj]/Matrix[i][maxj];
				for (int j = 0; j < MAXSIZE; ++j)
				{
					Matrix[k][j] = Matrix[k][j] - (Matrix[i][j]*kf); 
				}
			}
		}

		//if(rank==3){PrintMatrix(Matrix,0,MAXSIZE-1);}

	}	

	if(rank==0)
	{
		float payback=1;
		int mes_sent[MAXSIZE];
		for (int i = 1; i < count; ++i)
		{
			MPI_Recv(&payback,1,MPI_FLOAT,i,0,MPI_COMM_WORLD,&stat);
			RESULT=RESULT*payback;
			MPI_Recv(&mes_sent,MAXSIZE,MPI_INT,i,0,MPI_COMM_WORLD,&stat);
			for (int i = 0; i < MAXSIZE; ++i)
			{
				minus[i]=minus[i]+mes_sent[i];
			}
		}
	}
	else
	{
		MPI_Send(&RESULT,1,MPI_FLOAT,0,0,MPI_COMM_WORLD);
		MPI_Send(&minus,MAXSIZE,MPI_INT,0,0,MPI_COMM_WORLD);
	}	

	if(rank==0){
		int c =Sort(minus);
		printf("SORT=%d\n",c);
		if(c%2) {printf("RESULT: %f\n",RESULT*-1);}
		else {printf("RESULT: %f\n",RESULT);}
	bench_timer_stop();
	bench_timer_print();
	}

	if(rank!=0){
	bench_timer_stop();}
	MPI_Finalize();
}



