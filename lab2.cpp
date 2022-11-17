////#include <stdio.h>
//#include <iostream>
//#include <ctime>
//#include "mpi.h"
//#include "windows.h"
//
//int main(int argc, char* argv[])
//{
//	int ProcNum, ProcRank, recv = 0, count = 0, send;
//	bool flag = true;
//	MPI_Status Status;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//
//	while (true) {
//
//
//		//printf("Chislo %d", send);
//
//		MPI_Reduce(&send, &recv, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
//
//		if (recv == -1 || send == -1) {
//			break;
//		}
//		else {
//			count += ProcNum - 1;
//		}
//	}
//
//	if (ProcRank == 0) {
//		printf("Count = %d\n", count);
//	}
//	MPI_Finalize();
//	return 0;
//}