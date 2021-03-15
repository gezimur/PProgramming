#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

const unsigned int A = 576;

void InsertionSort(double* Arr, unsigned long N)
{
	unsigned int k, i;
	for (i = 0; i < N; i++)
	{
		k = i;
		while(k > 0 && Arr[k] <= Arr[k - 1])
		{
			double E  = Arr[k];
			Arr[k] = Arr[k - 1];
			Arr[k - 1] = E;
			
			k--;
		}
	}
}

void showArr(double* Arr, unsigned long N)
{
	printf("\n");
	int i;
	for(i = 0; i < N; i++){
		printf("%lf ", Arr[i]);
	}
	printf("\n");
}

int main(int argc, char* argv[])
{
	unsigned long i = 0, N;
	unsigned int ExpCount = 50;
	struct timeval StartTime, FinishTime, StartCalcTime, FinishCalcTime;
	long  TimeDelta_ms = 0;
	
	N = atoi(argv[1]);
	
	gettimeofday(&StartTime, NULL);
	
	//double* Arr1 = malloc(sizeof(double) * N * ExpCount);
	//double* Arr2 = malloc(sizeof(double) * N / 2 * ExpCount);
	unsigned long* RandV = malloc(sizeof(long) * (N + N / 2) * ExpCount);
	for (i = 0; i < ExpCount; i++)
	{
		//double* Arr1 = malloc(sizeof(double) * N);
		//double* Arr2 = malloc(sizeof(double) * N / 2);
		//unsigned long* RandV = malloc(sizeof(double) * (N + N / 2) );
		srand(i);
		
		unsigned long j;
		
		unsigned int s = i;
		for (j = 0; j < (N + N / 2); j++)
		{
			RandV[ (N + N / 2) * i + j] = rand_r(&s);
		}
		
	}//*/
	double* Arr1 = malloc(sizeof(double) * N);
	double* Arr2 = malloc(sizeof(double) * N / 2);
	gettimeofday(&StartCalcTime, NULL);
	
	for (i = 0; i < ExpCount; i++)
	{
		unsigned long j;
		unsigned long* RandVPtr = (RandV + (N + N / 2) * i);
		
		for (j = 0; j < N; j++)
		{
			Arr1[j] = ( RandVPtr[j] % A ) + (double)( RandVPtr[j] % 1000 ) / 1000;
		}
		RandVPtr += N;
		for (j = 0; j < N / 2; j++)
		{
			Arr2[j] = ( (RandVPtr[j] % (A * 9) ) + A) + (double)( RandVPtr[j] % 1000 ) / 1000;
		}
		
		for(j = 0; j < N; j++)
		{
			Arr1[ j] = pow((Arr1[ j] / exp(1.0)), 1.0 / 3.0);
		}
		
		/*
		for(j = N / 2 - 1; j > 0; j--)
		{
			Arr2[(N / 2) * i + j] = log(fabs(tan(Arr2[(N / 2) * i + j]+Arr2[ (N / 2) * i + j - 1])));
		}
		Arr2[(N / 2) * i] = log(fabs(tan(Arr2[(N / 2) * i])));//*/
		
		for (j  = 0; j < N / 2; j++)
		{
			Arr2[j] = pow(Arr1[ j], Arr2[ j]);
		}
		
		/*
		unsigned long j;
		//showArr((Arr2 + (N / 2 * i)), N / 2);
		InsertionSort( (Arr2 + (N / 2 * i)), N / 2);
		//showArr((Arr2 + (N / 2 * i)), N / 2);
		
		long double Sum = 0;
		for (j = 0; j < N / 2; j++){
			if ( (long)(round(Arr2[j]) / round(Arr2[0])) % 2 == 0 )
					Sum += sin(Arr2[j] ) ;
		}//*/
	}	
	gettimeofday(&FinishCalcTime, NULL);
	TimeDelta_ms = 1000 * (FinishCalcTime.tv_sec - StartCalcTime.tv_sec) + (FinishCalcTime.tv_usec - StartCalcTime.tv_usec) / 1000;
	printf("\nCalculated time = %ld", TimeDelta_ms);
	
	free(Arr1);
	free(Arr2);

	free(RandV);
	
	gettimeofday(&FinishTime, NULL);
	
	TimeDelta_ms = 1000 * (FinishTime.tv_sec - StartTime.tv_sec) + (FinishTime.tv_usec - StartTime.tv_usec) / 1000;
	
	printf("\nN=%ld. Milliseconds passsed: %ld\n", N, TimeDelta_ms);
	
	return 0;
}