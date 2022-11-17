//#include <iostream>
//#include "mpi.h"
//#include <ctime>
//
//int main(int argc, char* argv[])
//{
//	srand(time(0));
//	int pow = 5;
//	int ProcNum, ProcRank;
//	MPI_Status Status;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//	int respow = pow*ProcNum;
//	long int* sendarray = new long int[respow];
//	long int* temp = new long int[respow];
//	long int* result= new long int[respow];
//	MPI_Datatype polinom_type;
//	int count = 0;
//
//	
//	// Каждый процесс генерирует массив полиномов
//	for (int i = 0; i < respow; i++)
//		sendarray[i] = 0;
//	for (int i = 0; i < pow; i++)
//		sendarray[i] = (ProcRank + rand()) % 10 - 5;
//	MPI_Barrier(MPI_COMM_WORLD);
//	printf("Polinom %d ", ProcRank);
//	for (int i = 0; i < respow; i++)
//		printf("%d ", sendarray[i]);
//	printf("\n");
//	
//	
//	MPI_Type_contiguous(respow, MPI_INT, &polinom_type);
//	MPI_Type_commit(&polinom_type);
//	int k = 2, a = 0, b =1;
//	int max = pow;
//	// 0->1 2->3; 1->3
//	
//	while (k <= ProcNum) {
//		if ((ProcRank % k) == a)
//			MPI_Send(sendarray, 1, polinom_type, ProcRank + k / 2, 0, MPI_COMM_WORLD);
//
//		else if ((ProcRank % k) == b)
//		{
//			MPI_Recv(temp, 1, polinom_type, ProcRank - k/2, 0, MPI_COMM_WORLD, &Status);
//			for (int i = 0; i < respow; i++)
//				result[i] = 0;
//			for (int i = 0; i < pow; i++)
//				for (int j = 0; j < pow; j++) {
//					if (i + j > max)
//						max = i + j;
//					result[i + j] += (sendarray[i] * temp[j]);
//				}
//			for (int i = 0; i < respow; i++) 
//				sendarray[i] = result[i];
//		}
//		MPI_Barrier(MPI_COMM_WORLD);
//		pow = max+1;
//		k *= 2;
//		a = b;
//		b = k - 1;
//	}
//	
//	if (ProcRank == 3)
//	{
//		printf("Result\n");
//		for (int i = 0; i < respow; i++)
//			printf("%d ", sendarray[i]);
//		printf("\n");
//	}
//	delete[] result;
//	delete[] sendarray;
//	delete[] temp;
//
//	MPI_Type_free(&polinom_type);
//	MPI_Finalize();
//
//	return 0;
//}