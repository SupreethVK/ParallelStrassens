#include <stdio.h>
#include <stdlib.h>
#include <time.h>

float ** initMatrix(int n)
{
	float ** m = (float **) malloc(sizeof(float *) * n);
	for(int i=0; i<n; ++i)
	{
		m[i] = (float *) malloc(sizeof(float) * n);
	}
	return m;
}

float ** genMatrix(float ** m, int n)
{
	float ** backup = m;
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			m[i][j] = (((float)rand()/(float)RAND_MAX) * (rand()%100));
		}
	}
	return backup;
}

void displayMatrix(float ** m, int n)
{
	for(int i=0; i<n; ++i)
	{
		for(int j=0; j<n; ++j)
		{
			printf("%f ", m[i][j]);
		}
		printf("\n");
	}
}

int main() 
{
	int n;
	scanf("%d", &n);
	srand(time(NULL));
	float ** M1 = initMatrix(n);
	float ** M2 = initMatrix(n); 
	printf("%d\n", n);	
	M1 = genMatrix(M1, n);	
	displayMatrix(M1, n);
	M2 = genMatrix(M2, n);	
	displayMatrix(M2, n);
	return 0;
}
