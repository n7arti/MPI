#include <stdio.h>
#include <mpi.h>
#include <iostream>
#include <windows.h>
#include <time.h>
#include <memory>
#include <vector>
#include <string>
using namespace std;

void multiplyToomCook_3(int16_t*, int16_t*, int, int16_t*);

void carry(int16_t* a, int size) {
	int16_t cr = 0;
	for (int i = 0; i < size; i++) {
		a[i] += cr;
		if (a[i] < 0) {
			cr = -(-(a[i] + 1) / 10 + 1);
		}
		else {
			cr = a[i] / 10;
		}
		a[i] -= cr * 10;
	}
}

void ToomCook(int16_t* A, int16_t* B, int16_t* Z, int size)
{
	int16_t* a = new int16_t[size];
	int16_t* b = new int16_t[size];
	int i;

	for (i = 0; i < size; i++)
	{
		a[i] = A[i];
		b[i] = B[i];
	}
	
	multiplyToomCook_3(a, b, size, Z);
	carry(Z, size * 2);
}

void multiply(int16_t* a, int16_t* b, int size, int16_t* z)
{
	for (int i = 0; i < size * 2; i++)
		z[i] = 0;

	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
				z[i+j] += a[i] * b[j];
}

void multiplyToomCook_3(int16_t* m, int16_t* n, int size, int16_t* z)
{
	if (size <= 9)
	{
		multiply(m, n, size, z);
		return;
	}

	int16_t** a = new int16_t * [3];
	int16_t** b = new int16_t * [3];
	int16_t** r = new int16_t * [5];
	int16_t** p = new int16_t * [5];
	int16_t** q = new int16_t * [5];
	int16_t** r_k = new int16_t * [5];
	for (int i = 0; i < 3; i++) {
		a[i] = &m[i * size / 3];
		b[i] = &n[i * size / 3];
	}	
	for (int i = 0; i < 5; i++) {
		if (i % 2 == 0)
			r[i] = &z[(size / 3) * i];
		else 
			r[i] = new int16_t[(size / 3) * 2];
		p[i] = new int16_t[size / 3];
		q[i] = new int16_t[size / 3];
		r_k[i] = new int16_t[(size / 3)*2];

	} 
	//Вычисление в точках
	//p(0) = a0
	//p(1) = a2 + a1 + a0
	//p(-1) = a2 - a1 + a0
	//p(-2) = 4 * a2 - 2 * a1 + a0
	//p(inf) = a2
	int number = 0; 
	for (int i = 0; i < size / 3; i++)
	{
		for (int j = 0; j < 5; j++) {
			if (j == 2)
				number = -1;
			else if (j == 3)
				number = -2;
			else number = j;
			if (j == 4) {
				p[j][i] = a[2][i];
				q[j][i] = b[2][i];
			}
			else {
				p[j][i] = a[2][i] * pow(number, 2) + a[1][i] * number + a[0][i];
				q[j][i] = b[2][i] * pow(number, 2) + b[1][i] * number + b[0][i];
			}
		}
	}

	// поточечное умножение
	multiplyToomCook_3(p[0], q[0], size / 3, r[0]);
	multiplyToomCook_3(p[1], q[1], size / 3, r[1]);
	multiplyToomCook_3(p[2], q[2], size / 3, r[2]);
	multiplyToomCook_3(p[3], q[3], size / 3, r[3]);
	multiplyToomCook_3(p[4], q[4], size / 3, r[4]);

	//Интерполяция
	//r0 = r(0)
	//r1 = 3 * r(0) + 2 * r(1) - 6 * r(-1) + r(-2) - 12 * r(inf) / 6
	//r2 = - 6 * r(0) + 3 * r(1) + 3 * r(-1)   - 6 * r(inf) / 6
	//r3 = - 3 * r(0) + r(1) + 3 * r(-1) - r(-2) + 12 * r(inf) / 6
	//r4 = r(inf)
	for (int i = 0; i < (size / 3) * 2; i++) {
		r_k[0][i] = r[0][i];
		r_k[1][i] = (3 * r[0][i] + 2 * r[1][i] - 6 * r[2][i] + r[3][i] - 12 * r[4][i]) / 6;
		r_k[2][i] = (-6 * r[0][i] + 3 * r[1][i] + 3 * r[2][i] - 6 * r[4][i]) / 6;
		r_k[3][i] = (-3 * r[0][i] +  r[1][i] + 3 * r[2][i] - r[3][i] + 12 * r[4][i])/6;
		r_k[4][i] = r[4][i];
	}
	//z = r4 * x^4 + r3 * x^3 + r2 * x^2 + r1 * x + r0
	for (int i = 0; i < (size / 3) * 2; i++) {
			z[i + size / 3] += r_k[1][i];
	}
	for (int i = 0; i < (size / 3) * 2; i++) {
		z[i + (size / 3) * 3] += r_k[3][i];
	}
	delete[] a, b, r, r_k, p, q;
}

string printRes(int16_t* a, int16_t* b, int16_t* z, int size)
{
	string result = "";
	int aLen = size, bLen = size, zLen = 18;
	while (a[aLen - 1] == 0) 
		if (a[aLen - 1] == 0) 
			aLen--;
	while (b[bLen - 1] == 0) 
		if (b[bLen - 1] == 0) 
			bLen--;
	while (z[zLen - 1] == 0) 
		if (z[zLen - 1] == 0) 
			zLen--;
	result += "a = ";
	for (int i = aLen - 1; i >= 0; i--) {
		result += to_string(a[i]);
	}
	result += "\n";

	result += "b = ";
	for (int i = bLen - 1; i >= 0; i--) {
		result += to_string(b[i]);
	}
	result += "\n";

	result += "z = ";
	for (int i = zLen - 1; i >= 0; i--) {
		result += to_string(z[i]);
	}
	result += "\n";

	return result;
}

void removeRanksNumbers(vector<int>& v, int numtasks)
{
	for (int i = v.size() - 1; i > 0; i--)
		v.erase(v.begin() + i);
}

string getNumber(int16_t* number, int size)
{
	string result = "";
	int i;
	int aLen = 72;
	while (number[aLen - 1] == 0) 
		if (number[aLen - 1] == 0) 
			aLen--;
	for (i = aLen - 1; i >= 0; i--) {
		result += to_string(number[i]);
	}
	result += "\n";
	return result;
}


int main(int* argc, char** argv)
{
	int ProcNum, ProcRank;
	int size = 9;
	double dt = 0;
	
	MPI_Status status;
	MPI_Init(argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	int resultSize = 3 * size * pow(3, ProcNum);
	MPI_Datatype MPI_LONG_NUMBER;
	MPI_Type_contiguous(resultSize, MPI_SHORT, &MPI_LONG_NUMBER);
	MPI_Type_commit(&MPI_LONG_NUMBER);

	MPI_Group MY_GROUP;
	MPI_Comm MY_COMM;
	MPI_Comm_group(MPI_COMM_WORLD, &MY_GROUP);
	MPI_Comm_create(MPI_COMM_WORLD, MY_GROUP, &MY_COMM);

	srand(time(NULL) + ProcRank);
	if (ProcRank == 0)
		dt -= MPI_Wtime();

	int16_t* A = new int16_t[resultSize];
	int16_t* B = new int16_t[resultSize];
	int16_t* Z = new int16_t[resultSize];

	for (int i = 0; i < size; i++)
	{
		A[i] = rand() % 10;
		B[i] = rand() % 10;
		Z[i] = 0;
	}
	for (int i = size; i < resultSize; i++)
	{
		A[i] = 0;
		B[i] = 0;
		Z[i] = 0;
	}
	
	ToomCook(A, B, Z, resultSize);
	
	cout << to_string(ProcRank) + "\n" + printRes(A, B, Z, resultSize) << endl;

	vector<int> newGroup;
	for (int i = 0; i < ProcNum; i++)
	{
		newGroup.push_back(i);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	int iteration = 1;

	while (newGroup.size() != 1)
	{
		MPI_Group_incl(MY_GROUP, newGroup.size(), newGroup.data(), &MY_GROUP);

		MPI_Comm_create(MY_COMM, MY_GROUP, &MY_COMM);

		for (int i = 0; i < newGroup.size(); i++)
		{
			if (newGroup[i] == ProcRank)
			{
				if (i % 2 == 0)
				{
					for (int i = 0; i < resultSize; i++)
						A[i] = Z[i];
				}
				else
				{
					MPI_Send(Z, 1, MPI_LONG_NUMBER, newGroup[i - 1] / iteration, 0, MY_COMM);
				}
				break;
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		for (int i = 0; i < newGroup.size(); i++)
		{
			if (newGroup[i] == ProcRank)
			{
				if (i % 2 == 0)
				{
					MPI_Recv(B, 1, MPI_LONG_NUMBER, newGroup[i + 1] / iteration, 0, MY_COMM, &status);
					ToomCook(A, B, Z, resultSize);
				}
			}
		}
		removeRanksNumbers(newGroup, ProcNum);

		iteration *= 2;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	dt += MPI_Wtime();
	if (ProcRank == 0)
	{
		cout << "Result:\n" + getNumber(Z, resultSize) << endl;
		cout << "Time: "<<  dt << endl;
	}


	MPI_Group_free(&MY_GROUP);
	MPI_Comm_free(&MY_COMM);

	delete[] A, B, Z;
	MPI_Finalize();

	return 0;
}