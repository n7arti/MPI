#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include "mpi.h"
#include <iostream>
#include <time.h>
#include <math.h>
#include <string>
#include <vector>
#include <complex>

using namespace std;
typedef complex<double> base;

//Алгоритм быстрого преобразования Фурье
void fft(vector<base>& a, bool invert) {
	int n = (int)a.size();
	if (n == 1)  return;

	vector<base> a0(n / 2), a1(n / 2);
	for (int i = 0, j = 0; i < n; i += 2, ++j) {
		a0[j] = a[i];
		a1[j] = a[i + 1];
	}
	fft(a0, invert);
	fft(a1, invert);

	double ang = 2 * M_PI / n * (invert ? -1 : 1);
	base w(1), wn(cos(ang), sin(ang));
	for (int i = 0; i < n / 2; ++i) {
		a[i] = a0[i] + w * a1[i];
		a[i + n / 2] = a0[i] - w * a1[i];
		if (invert)
			a[i] /= 2, a[i + n / 2] /= 2;
		w *= wn;
	}
}

void multiply(const vector<int>& a, const vector<int>& b, vector<int>& res) {
	vector<base> fa(a.begin(), a.end()), fb(b.begin(), b.end());
	size_t n = 1;
	//Получаем новую длину для числа
	while (n < max(a.size(), b.size()))
		n <<= 1;
	n <<= 1;
	fa.resize(n), fb.resize(n);

	fft(fa, false), fft(fb, false);
	for (size_t i = 0; i < n; ++i)
		fa[i] *= fb[i];
	fft(fa, true);

	res.resize(n);
	for (size_t i = 0; i < n; ++i)
		res[i] = int(fa[i].real() + 0.5);
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	int length = 10;
	int lengthResult = (length - 1) * 2 + 1;
	double dt = 0;
	int ProcNum, ProcRank, RecvRank, CalcRank;

	vector<int> A(lengthResult);
	vector<int> B(lengthResult);
	int* result = new int[lengthResult];
	int* num1 = new int[length];
	int* num2 = new int[length];
	int* num3 = new int[length];
	for (int i = 0; i < length; i++) {
		num1[i] = rand() % 10;
		num2[i] = rand() % 10;
	}
	MPI_Status Status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		printf("Number 1: ");
		for (int i = length - 1; i >= 0; i--)
			printf("%d", num1[i]);
		printf("\n");
		printf("Number 2: ");
		for (int i = length - 1; i >= 0; i--)
			printf("%d", num2[i]);
		printf("\n");
		dt -= MPI_Wtime();
	}

	//Создаем новый тип числа
	MPI_Datatype number;
	MPI_Type_contiguous(length, MPI_INT, &number);
	MPI_Type_commit(&number);

	MPI_Datatype numberResult;
	MPI_Type_contiguous(lengthResult, MPI_INT, &numberResult);
	MPI_Type_commit(&numberResult);

	MPI_Comm GridComm; // коммуникатор для декартовых решеток
	//Создаем двумерную решетку 4х4
	int GridNum, GridRank;
	int dims[2], periods[2], reorder = 1;
	dims[0] = dims[1] = 2;
	periods[0] = periods[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &GridComm);

	//Коммуникаторы для строк и столбцов
	MPI_Comm RowComm, ColComm;
	MPI_Comm_size(GridComm, &GridNum);
	MPI_Comm_rank(GridComm, &GridRank);

	//Коммуникаторы для строк
	int subdims[2];
	subdims[0] = 0;
	subdims[1] = 1;
	MPI_Cart_sub(GridComm, subdims, &RowComm);

	//Коммуникаторы для столбцов
	subdims[0] = 1;
	subdims[1] = 0;
	MPI_Cart_sub(GridComm, subdims, &ColComm);
	int coord[2];
	MPI_Cart_coords(GridComm, GridRank, dims[0], coord);
	int coordsum = coord[0] + coord[1];
	int RowRank, RowNum, ColNum, ColRank;
	MPI_Comm_size(RowComm, &RowNum);
	MPI_Comm_rank(RowComm, &RowRank);
	MPI_Comm_size(ColComm, &ColNum);
	MPI_Comm_rank(ColComm, &ColRank);


	if (coordsum % 2 == 0)
	{
		//(0,0)
		if (GridRank == 0)
		{
			MPI_Send(num1, 1, number, RowRank + 1, 0, RowComm);
			MPI_Send(num2, 1, number, ColRank + 1, 0, ColComm);
		}
		//(1,1)
		else
		{
			for (int i = 1; i < GridNum - 1; i++)
			{
				if (i % 2 == 1)
				{
					MPI_Recv(result, 1, numberResult, 0, 0, ColComm, &Status);
				}
				else
				{
					MPI_Recv(result, 1, numberResult, 0, 0, RowComm, &Status);
				}
				if (i == 1)
				{
					for (int j = 0; j < lengthResult; j++)
					{
						A[j] = result[j];
					}
				}
				else
				{
					for (int j = 0; j < lengthResult; j++)
					{
						B[j] = result[j];
					}
					multiply(A, B, A);
				}
			}
			for (int i = 0; i < lengthResult; i++)
			{
				if (i < A.size())
					result[i] = A[i];
				else
					result[i] = 0;
			}
			for (int i = 0; i < lengthResult - 1; i++)
			{
				if (result[i] / 10 > 0)
				{
					result[i + 1] += result[i] / 10;
					result[i] %= 10;
				}
			}
			dt += MPI_Wtime();
			printf("Result: ");
			for (int i = lengthResult - 1; i >= 0; i--)
				printf("%d", result[i]);
			printf("\n");
			printf("Time: %d\n", dt);
		}
	}
	// (0,1) и (1,0)
	else
	{
		A.resize(length);
		B.resize(length);
		int recvcount = 1;
		for (int i = 0; i < recvcount; i++)
		{
			if (coord[0] == 1)
			{
				MPI_Recv(num3, 1, number, 0, 0, ColComm, &Status);
			}
			else
			{
				MPI_Recv(num3, 1, number, 0, 0, RowComm, &Status);
			}
			if (i == 0)
			{
				for (int j = 0; j < length; j++)
				{
					A[j] = num3[j];
				}
			}
			else
			{
				for (int j = 0; j < length; j++)
				{
					B[j] = num3[j];
				}
				multiply(A, B, A);
			}
		}
		for (int i = 0; i < lengthResult; i++)
		{
			if (i < A.size())
				result[i] = A[i];
			else  
				result[i] = 0;
		}
		// отправляем (1,1)
		if (RowRank == 1)
		{
			MPI_Send(result, 1, numberResult, ColRank + 1, 0, ColComm);
		}
		if (ColRank == 1)
		{
			MPI_Send(result, 1, numberResult, RowRank + 1, 0, RowComm);
		}
	}
	MPI_Type_free(&number);
	MPI_Type_free(&numberResult);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;

}
