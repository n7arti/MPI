////#include <iostream>
//#include "mpi.h"
//#include <ctime>
//
//int main(int argc, char* argv[])
//{
//	srand(time(0));
//	int length = 3;
//	int ProcNum, ProcRank;
//	MPI_Status Status;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//	int resLenght = length*ProcNum;
//	printf("%d", resLenght);
//	long int* sendarray = new long int[resLenght];
//	long int* temp = new long int[resLenght];
//	long int* result= new long int[resLenght];
//	MPI_Datatype long_number;
//	int count = 0;
//
//	
//	// Каждый процесс генерирует массив полиномов
//	for (int i = 0; i < resLenght; i++)
//		sendarray[i] = 0;
//	for (int i = 0; i < length; i++)
//		sendarray[i] = rand() % 10;
//	MPI_Barrier(MPI_COMM_WORLD);
//	printf("Number %d: ", ProcRank);
//	for (int i = 0; i < length; i++)
//		printf("%d", sendarray[i]);
//	printf("\n");
//	
//	
//	MPI_Type_contiguous(resLenght, MPI_INT, &long_number);
//	MPI_Type_commit(&long_number);
//	int k = 2, a = 0, b = 1;
//	int max = length;
//	// 0->1 2->3; 1->3
//	
//	while (k <= ProcNum) {
//		if ((ProcRank % k) == a)
//			MPI_Send(sendarray, 1, long_number, ProcRank + k / 2, 0, MPI_COMM_WORLD);
//
//		else if ((ProcRank % k) == b)
//		{
//			MPI_Recv(temp, 1, long_number, ProcRank - k / 2, 0, MPI_COMM_WORLD, &Status);
//			for (int i = 0; i < resLenght; i++)
//				result[i] = 0;
//			for (int i = 0; i < length; i++)
//				for (int j = 0; j < length; j++) {
//					if (i + j > max)
//						max = i + j;
//					result[i + j] += (sendarray[i] * temp[j]);
//				}
//			for (int i = resLenght - 1; i >= 0; i--) {
//				if (result[i] > 10) {
//					if (i != 0) {
//						sendarray[i] = result[i] % 10;
//						result[i - 1] += result[i] / 10;
//					}
//					else {
//						sendarray[i] = result[i];
//					}
//				}
//				else {
//					sendarray[i] = result[i];
//				}
//			}
//				
//		}
//		MPI_Barrier(MPI_COMM_WORLD);
//		length = max+1;
//		k *= 2;
//		a = b;
//		b = k - 1;
//	}
//	
//	if (ProcRank == 3)
//	{
//		printf("Result\n");
//		for (int i = 0; i < length; i++)
//			printf("%d", sendarray[i]);
//		printf("\n");
//	}
//	delete[] result;
//	delete[] sendarray;
//	delete[] temp;
//
//	MPI_Type_free(&long_number);
//	MPI_Finalize();
//
//	return 0;
//}