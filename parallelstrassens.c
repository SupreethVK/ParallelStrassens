/* Parallel Strassen's Matrix multiplication */

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<pthread.h>

//Structure to hold two matrices and their dimension
typedef struct 
{
	float ** M1;
	float ** M2;
	int n; 
} mat;

//Structure to hold a single matrix and its dimension
typedef struct
{
	float ** m;
	int n;
} matS;

//function to calculate execution time
double calc_time(struct timespec start, struct timespec stop) 
{
	double t;
	t = (stop.tv_sec - start.tv_sec) * 1000; 
	t += (stop.tv_nsec - start.tv_nsec) * 0.000001; 
	return t;
}

//function to allocate memory to the matrix types
float ** initMatrix(int n)
{
	float ** m = (float **) calloc(n, sizeof(float *));
	for(int i=1; i<=n; ++i)
	{
		m[i] = (float *) calloc(n, sizeof(float));
	}
	return m;
}

//function to initialize the matrices with values from the input
float ** initMatrixV(int n)
{
	float ** m = (float **) calloc(n, sizeof(float *));
	for(int i=1; i<=n; ++i)
	{
		m[i] = (float *) calloc(n, sizeof(float));
	}
	for(int i=1; i<=n; ++i)
	{
		for(int j=1; j<=n; ++j)
		{
			scanf("%f", &m[i][j]);
		}
	}
	return m;
}

void displayMatrix(float ** m, int n)
{
	printf("\n-----------------------------------------------------------------\n");
	for(int i=1; i<=n; ++i)
	{
		for(int j=1; j<=n; ++j)
		{
			printf("%f ", m[i][j]);
		}
		printf("\n");
	}
}

//function to allocate memory to the two-matrix structure defined above
mat * initStruct2(int N)
{
	mat * m = (mat *) malloc(sizeof(mat));
	m->M1 = initMatrix(N);
	m->M2 = initMatrix(N);
	m->n = N;
	return m;
}

//function to allocate memory to the one-matrix structure defined above
matS * initStruct(int N)
{
	matS * M = (matS *) malloc(sizeof(matS));
	M->m = initMatrix(N);
	M->n = N;
	return M;
}

//function to free allocated memory from matrix type
void freeMat(float ** m, int N)
{
	for(int i=1; i<N; ++i)
	{
		free(m[i]);
	}
	free(m);
}

//function to free allocated memory from mat type structures
void freeStruct2(mat * M)
{
	freeMat(M->M1, M->n);
	freeMat(M->M2, M->n);
	free(M);
}

//function to free allocated memory from matS type structures
void freeStruct(matS * M)
{
	freeMat(M->m, M->n);
	free(M);
}

//function to add two matrices
float ** addMatrix(float ** a, float ** b, int n)
{
	float ** c = initMatrix(n);
	for(int i=1; i<=n; ++i)
	{
		for(int j=1; j<=n; ++j)
		{
			c[i][j] = a[i][j] + b[i][j];
		}
	}
	return c;
}

float ** subMatrix(float ** a, float ** b, int n)
{
	float ** c = initMatrix(n);
	for(int i=1; i<=n; ++i)
	{
		for(int j=1; j<=n; ++j)
		{
			c[i][j] = a[i][j] - b[i][j];
		}
	}
	return c;
}

//Thread Routine that is called recursively
void * remul(void * M)
{
	mat * MM = (mat *) M;
	int ord = MM->n;
	int N = ord/2;
	
	matS * c = initStruct(ord); //allocating memory to the resultant matrix
	
	if(ord==1)
	{
		c->m[1][1] = MM->M1[1][1] * MM->M2[1][1] ;
	}
	else
	{	
		float ** temp;
		//allocating memory to the submatrices
		float ** A11 = initMatrix(N);
		float ** A12 = initMatrix(N);
		float ** A21 = initMatrix(N);
		float ** A22 = initMatrix(N);
		float ** B11 = initMatrix(N);
		float ** B12 = initMatrix(N);
		float ** B21 = initMatrix(N);
		float ** B22 = initMatrix(N);
		//initialising the submatrices
		for(int i=1; i<=(N); ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				A11[i][j] = MM->M1[i][j];
				B11[i][j] = MM->M2[i][j];
			}
		}
		int k=1;
		for(int i=(N)+1; i<=ord; ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				A21[k][j] = MM->M1[i][j];
				B21[k][j] = MM->M2[i][j];
			}
			k++;
		}
		int l = 1;
		for(int i=1; i<=(N); ++i)
		{
			for(int j=(N)+1; j<=ord; ++j)
			{
				A12[i][l] = MM->M1[i][j];
				B12[i][l] = MM->M2[i][j];
				l++;
			}
			l = 1;
		}
		k = 1;
		l = 1;
		for(int i=(N)+1; i<=ord; ++i)
		{
			for(int j=(N)+1; j<=ord; ++j)
			{
				A22[k][l] = MM->M1[i][j];
				B22[k][l] = MM->M2[i][j];
				l++;
			}
			l = 1;
			k++;
		}
		//Creating the structures which will be used to store the matrix pairs following the algorithm
		mat * m1 = initStruct2(N);//A11, S1
		m1->M1 = A11;
		m1->M2 = subMatrix(B12, B22, N);
		mat * m2 = initStruct2(N);//S2, B22
		m2->M1 = addMatrix(A11, A12, N);
		m2->M2 = B22;
		mat * m3 = initStruct2(N);//S3, B11
		m3->M1 = addMatrix(A21, A22, N);
		m3->M2 = B11;
		mat * m4 = initStruct2(N);//A22, S4	
		m4->M1 = A22;
		m4->M2 = subMatrix(B21, B11, N);
		mat * m5 = initStruct2(N);//S5, S6	
		m5->M1 = addMatrix(A11, A22, N);
		m5->M2 = addMatrix(B11, B22, N);
		mat * m6 = initStruct2(N);//S7, S8
		m6->M1 = subMatrix(A12, A22, N);
		m6->M2 = addMatrix(B21, B22, N);
		mat * m7 = initStruct2(N);//S9, S10	
		m7->M1 = subMatrix(A11, A21, N);
		m7->M2 = addMatrix(B11, B12, N);
		
		//Recursive calls 
		matS * p1 = remul((void *) m1);
		matS * p2 = remul((void *) m2);
		matS * p3 = remul((void *) m3);
		matS * p4 = remul((void *) m4);
		matS * p5 = remul((void *) m5);
		matS * p6 = remul((void *) m6);
		matS * p7 = remul((void *) m7);
		
		freeStruct2(m1);
		freeStruct2(m2);
		freeStruct2(m3);
		freeStruct2(m4);
		freeStruct2(m5);
		freeStruct2(m6);
		freeStruct2(m7);
		
		//After getting the resultant matrices from the recursive calls, 
		//pair matrices for further operations specified by the algorithm
		float ** t1;
		t1 = addMatrix(p5->m, p4->m, N);
		float ** t2;
		t2 = addMatrix(t1, p6->m, N);
		temp = subMatrix(t2, p2->m, N);
		for(int i=1; i<=(N); ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				c->m[i][j] = temp[i][j];
			}
		}	
		freeMat(t1, N);
		freeMat(t2, N);
		freeMat(temp, N);	
		temp = addMatrix(p1->m, p2->m, N);
		l = 1;
		for(int i=1; i<=(N); ++i)
		{
			for(int j=(N)+1; j<=ord; ++j)
			{
				c->m[i][j] = temp[i][l];
				l++;
			}
			l = 1;
		}
		freeMat(temp, N);
	
	
		temp = addMatrix(p3->m, p4->m, N);
		k = 1;
		for(int i=(N)+1; i<=ord; ++i)
		{
			for(int j=1; j<=(N); ++j)
			{
				c->m[i][j] = temp[k][j];		
			}
			k++;
		}
		freeMat(temp, N);	
		
		t1 = addMatrix(p5->m, p1->m, N);
		t2 = subMatrix(t1, p3->m, N);
		temp = subMatrix(t2, p7->m, N);
		k = 1;
		l = 1;
		for(int i=(N)+1; i<=ord; ++i)
		{
			for(int j=(N)+1; j<=ord; ++j)
			{
				c->m[i][j] = temp[k][l];
				l++;
			}
			l = 1;
			k++;
		}
		freeMat(t1, N);
		freeMat(t2, N);
		freeMat(temp, N);
		freeStruct(p1);
		freeStruct(p2);
		freeStruct(p3);
		freeStruct(p4);
		freeStruct(p5);
		freeStruct(p6);
		freeStruct(p7);		
	}
	return (void *) c;
}
		
//Wrapper function for the routine function,
//Here is where the threads are created and initialised		
matS * multiply(mat * M)
{
	int n = M->n;
	int N = n/2;
	pthread_t tid[8];
	
	void * stat1;
	void * stat2;
	void * stat3;
	void * stat4;
	void * stat5;
	void * stat6;
	void * stat7;
	
	matS * p1;
	matS * p2;
	matS * p3;
	matS * p4;
	matS * p5;
	matS * p6;
	matS * p7;
	
	float ** temp;
	matS * c = initStruct(n);
	
	float ** A11 = initMatrix(N);
	float ** A12 = initMatrix(N);
	float ** A21 = initMatrix(N);
	float ** A22 = initMatrix(N);
	float ** B11 = initMatrix(N);
	float ** B12 = initMatrix(N);
	float ** B21 = initMatrix(N);
	float ** B22 = initMatrix(N);
	for(int i=1; i<=(N); ++i)
	{
		for(int j=1; j<=(N); ++j)
		{
			A11[i][j] = M->M1[i][j];
			B11[i][j] = M->M2[i][j];
		}
	}
	int k=1;
	for(int i=(N)+1; i<=n; ++i)
	{
		for(int j=1; j<=(N); ++j)
		{
			A21[k][j] = M->M1[i][j];
			B21[k][j] = M->M2[i][j];
		}
		k++;
	}
	int l = 1;
	for(int i=1; i<=(N); ++i)
	{
		for(int j=(N)+1; j<=n; ++j)
		{
			A12[i][l] = M->M1[i][j];
			B12[i][l] = M->M2[i][j];
			l++;
		}
		l = 1;
	}
	k = 1;
	l = 1;
	for(int i=(N)+1; i<=n; ++i)
	{
		for(int j=(N)+1; j<=n; ++j)
		{
			A22[k][l] = M->M1[i][j];
			B22[k][l] = M->M2[i][j];
			l++;
		}
		l = 1;
		k++;
	}
	
	//Splitting the input matrices into their submatrices and pairing them accordingly
	mat * m1 = initStruct2(N);//A11, S1
	m1->M1 = A11;
	m1->M2 = subMatrix(B12, B22, N);
	
	mat * m2 = initStruct2(N);//S2, B22
	m2->M1 = addMatrix(A11, A12, N);
	m2->M2 = B22;
	
	mat * m3 = initStruct2(N);//S3, B11
	m3->M1 = addMatrix(A21, A22, N);
	m3->M2 = B11;
	
	mat * m4 = initStruct2(N);//A22, S4	
	m4->M1 = A22;
	m4->M2 = subMatrix(B21, B11, N);
	
	mat * m5 = initStruct2(N);//S5, S6	
	m5->M1 = addMatrix(A11, A22, N);
	m5->M2 = addMatrix(B11, B22, N);
	
	mat * m6 = initStruct2(N);//S7, S8
	m6->M1 = subMatrix(A12, A22, N);
	m6->M2 = addMatrix(B21, B22, N);
	
	mat * m7 = initStruct2(N);//S9, S10	
	m7->M1 = subMatrix(A11, A21, N);
	m7->M2 = addMatrix(B11, B12, N);
	
	//Creating the threads and passing the matrix pairs to the thread routine
	pthread_create(&tid[0], NULL, remul, (void*) m1);
	pthread_create(&tid[1], NULL, remul, (void*) m2);
	pthread_create(&tid[2], NULL, remul, (void*) m3);
	pthread_create(&tid[3], NULL, remul, (void*) m4);
	pthread_create(&tid[4], NULL, remul, (void*) m5);
	pthread_create(&tid[5], NULL, remul, (void*) m6);
	pthread_create(&tid[6], NULL, remul, (void*) m7);
	
	//Collecting the threads and storing the return values in the second parameter passed
	pthread_join(tid[0], &stat1);
	
	p1 = (matS *) stat1; //A11, B11
	pthread_join(tid[1], &stat2);
	
	p2 = (matS *) stat2;//A12, B21
	pthread_join(tid[2], &stat3);
	
	p3 = (matS *) stat3;//A11, B12
	pthread_join(tid[3], &stat4);
	
	p4 = (matS *) stat4;//A12, B22
	pthread_join(tid[4], &stat5);
	
	p5 = (matS *) stat5;//A21, B11
	pthread_join(tid[5], &stat6);
	
	p6 = (matS *) stat6;//A22, B21
	pthread_join(tid[6], &stat7);
	
	p7 = (matS *) stat7;//A21, B12
	
	freeStruct2(m1);
	freeStruct2(m2);
	freeStruct2(m3);
	freeStruct2(m4);
	freeStruct2(m5);
	freeStruct2(m6);
	freeStruct2(m7);
	
	//Resultant matrices are paired accordingly for further operations as specified by the algorithm
	float ** t1;
	t1 = addMatrix(p5->m, p4->m, N);
	float ** t2;
	t2 = addMatrix(t1, p6->m, N);
	temp = subMatrix(t2, p2->m, N);
	for(int i=1; i<=(N); ++i)
	{
		for(int j=1; j<=(N); ++j)
		{
			c->m[i][j] = temp[i][j];
		}
	}	
	freeMat(t1, N);
	freeMat(t2, N);
	freeMat(temp, N);	
	temp = addMatrix(p1->m, p2->m, N);
	l = 1;
	for(int i=1; i<=(N); ++i)
	{
		for(int j=(N)+1; j<=n; ++j)
		{
			c->m[i][j] = temp[i][l];
			l++;
		}
		l = 1;
	}
	freeMat(temp, N);
	
	
	temp = addMatrix(p3->m, p4->m, N);
	k = 1;
	for(int i=(N)+1; i<=n; ++i)
	{
		for(int j=1; j<=(N); ++j)
		{
			c->m[i][j] = temp[k][j];		
		}
		k++;
	}
	freeMat(temp, N);	
		
	t1 = addMatrix(p5->m, p1->m, N);
	t2 = subMatrix(t1, p3->m, N);
	temp = subMatrix(t2, p7->m, N);
	k = 1;
	l = 1;
	for(int i=(N)+1; i<=n; ++i)
	{
		for(int j=(N)+1; j<=n; ++j)
		{
			c->m[i][j] = temp[k][l];
			l++;
		}
		l = 1;
		k++;
	}
	freeMat(t1, N);
	freeMat(t2, N);
	freeMat(temp, N);
	freeStruct(p1);
	freeStruct(p2);
	freeStruct(p3);
	freeStruct(p4);
	freeStruct(p5);
	freeStruct(p6);
	freeStruct(p7);	
	return c;
}


int main()
{
	struct timespec start, stop;	
	clock_t start_t, end_t;
	int N;
	scanf("%d", &N);
	mat * Matrix = (mat *) malloc(sizeof(mat));
	Matrix->n = N;
	Matrix->M1 = initMatrixV(N);
	Matrix->M2 = initMatrixV(N);
	//displayMatrix(Matrix->M1, Matrix->n);
	//displayMatrix(Matrix->M2, Matrix->n);
	matS * res;
	clock_gettime(CLOCK_REALTIME,&start);
	start_t = clock();
	res = multiply(Matrix);
	//displayMatrix(res->m, res->n);
	end_t = clock();
	clock_gettime(CLOCK_REALTIME,&stop);
	
	printf("%lf milliseconds\n",calc_time(start,stop));
	printf("Clock ticks: %ld\n", end_t - start_t);
	freeStruct2(Matrix);
	freeStruct(res);
	return 0;
}
