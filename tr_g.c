#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


#define MAXSIZE 5
#define NUM_THREADS 8
//#define OTLADKA_ARRAY 1
//#define OTLADKA_MATRIX 1
//#define OTLADKA_CALC 1
//#define OTLADKA_THREADS 1

float Matrix[MAXSIZE][MAXSIZE];
int minus[MAXSIZE]={0};

void FillMatrix()
{
	for (int i = 0; i < MAXSIZE; ++i)
	{
		for (int j = 0; j < MAXSIZE; ++j)
		{
			//Matrix[i][j] = i+j+i*j+(i+1)/(j+1) - (i+1)%(j+1);// RESULT TEST [3,4,5,6,7] Matrix
			//if (i==j) Matrix[i][j]=1; else Matrix[i][j]=0;  // TIME TEST 
			//Matrix[i][j] = MAXSIZE*i+j+1;					 // Another Result test
			//Matrix[i][j]=0;							    // Funny TEst					
		}
	}

	// Matrix[0][0]=2;Matrix[0][1]=2;Matrix[0][2]=1;Matrix[0][3]=1; // Izlom test
	// Matrix[1][0]=2;Matrix[1][1]=2;Matrix[1][2]=45;Matrix[1][3]=24;
	// Matrix[2][0]=0;Matrix[2][1]=34;Matrix[2][2]=129;Matrix[2][3]=254;
	// Matrix[3][0]=934;Matrix[3][1]=113;Matrix[3][2]=0;Matrix[3][3]=132;

	// Matrix[0][0]=1;Matrix[0][1]=2;Matrix[0][2]=3;Matrix[0][3]=4;Matrix[0][4]=5; //wolfram test
	// Matrix[1][0]=0;Matrix[1][1]=2;Matrix[1][2]=3;Matrix[1][3]=4;Matrix[1][4]=5;
	// Matrix[2][0]=1;Matrix[1][1]=0;Matrix[2][2]=3;Matrix[2][3]=4;Matrix[2][4]=5;
	// Matrix[3][0]=1;Matrix[3][1]=2;Matrix[3][2]=0;Matrix[3][3]=4;Matrix[3][4]=5;
	// Matrix[4][0]=1;Matrix[4][1]=2;Matrix[4][2]=3;Matrix[4][3]=0;Matrix[4][4]=5;


}

void PrintMatrix()
{
	for (int i = 0; i < MAXSIZE; ++i)
	{
		printf("\n");
		for (int j = 0; j < MAXSIZE; ++j)
		{
			printf("%f    ",Matrix[i][j]);
		}
	}
	printf("\n");
}

void PrintArray()
{
	for (int i = 0; i < MAXSIZE; ++i)
	{
		printf("%d,  ",minus[i]);
	}
	printf("\n");
}

int Sort()
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

void Gaus()
{
	float result=1;
	for (int i = 0; i < MAXSIZE; ++i)
	{
		int threads;
		if (MAXSIZE < NUM_THREADS) threads=MAXSIZE; else threads=NUM_THREADS;
		float elem[NUM_THREADS]={0};
		int   oper[NUM_THREADS]={0};

		#pragma omp parallel num_threads(threads) 
		{	
			int maxj=0;
			float maxel = 0;
			float cand;

			int id = omp_get_thread_num();
			int num = omp_get_num_threads();
			#ifdef OTLADKA_THREADS
				printf("id=%d;num=%d;\n",id,num);
			#endif
			for (int j = id; j < MAXSIZE; j=j+num)
			{
				cand = Matrix[i][j];
				if(Matrix[i][j] < 0) cand=cand*(-1);
				if(cand > maxel) {maxel=cand;maxj=j;}
			}
			elem[id]=maxel;
			oper[id]=maxj;
		}
		float maxel=elem[0];
		int   maxj=oper[0];
		for (int j = 0; j < threads; ++j)
			{
				if (elem[j] > maxel) {maxel=elem[j];maxj=oper[j];} 
				#ifdef OTLADKA_THREADS	
					printf("-maxel=%f;maxj=%d;\n",elem[j],oper[j]);
				#endif
			}	
		#ifdef OTLADKA_THREADS
			printf("---->maxel=%f;maxj=%d;\n",maxel,maxj);
		#endif
	
		

		#ifdef OTLADKA_CALC
			printf("Строка: %d; Столбец: %d; Число: %f\n",i+1,maxj+1,maxel);
		#endif

		if (maxel == 0 ) {result=0; break;}
		result=result* (Matrix[i][maxj]);	
		minus[maxj]=i+1;
		#ifdef OTLADKA_ARRAY
			PrintArray();
		#endif

		float koef = 0;	
		#pragma omp parallel for private(koef)
		for (int k = i+1; k < MAXSIZE; ++k)
		{
			koef=Matrix[k][maxj]/Matrix[i][maxj];
			#ifdef OTLADKA_CALC
			   printf("KF = %f\n",koef);
			#endif
			//#pragma omp parallel for private(koef)	
			for (int j = 0; j < MAXSIZE; ++j)
			{
				#ifdef OTLADKA_CALC
					printf("Decree: %f\n",(Matrix[i][j]*koef));
				#endif
				Matrix[k][j] = Matrix[k][j] - (Matrix[i][j]*koef); 
			}
		}
		#ifdef OTLADKA_MATRIX
			PrintMatrix();
		#endif	
	}
	printf("\n");
	printf("*************************---*************************\n");
	printf("\n");
	int c =Sort();
	printf("SORT:%d\n",c);
	#ifdef OTLADKA_ARRAY
		PrintArray();
	#endif
	if(c%2) {printf("RESULT:%f\n",result*-1);}
	else {printf("RESULT:%f\n",result);}
	printf("\n");
	printf("*************************---*************************\n");
	//printf("\n");
}

int main ()
{
	FillMatrix();
	#ifdef OTLADKA_MATRIX
		PrintMatrix();
	#endif
	Gaus();
}