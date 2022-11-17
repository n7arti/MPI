//#include <stdio.h>
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
//	while (flag) {
//		if (ProcRank == 0)
//		{
//			// Действия, выполняемые только процессом с рангом 0
//			printf("\n Hello from process %3d", ProcRank);
//			for (int i = 1; i < ProcNum; i++) {
//				MPI_Recv(&recv, 1, MPI_INT, MPI_ANY_SOURCE,
//					MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
//				printf("\n Recv %3d\n", recv);
//				if ((recv != -1) && flag) {
//					count++;
//					printf("\n Count++");
//				}
//				if ((recv == -1) && flag) {
//					flag = false;
//					printf("\n Count = %3d\n", count);
//				}
//				MPI_Send(&flag, 1, MPI_C_BOOL, i, 0, MPI_COMM_WORLD);
//			}
//
//		}
//		else // Сообщение, отправляемое всеми процессами, кроме процесса с рангом 0
//		{
//			static time_t tval = time(0);
//			tval += 10;
//			srand(tval);
//			send = rand() % 5 - 1;
//			//printf("\n Process %3d send %3d\n", ProcRank, send);
//			MPI_Send(&send, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
//			MPI_Recv(&flag, 1, MPI_C_BOOL, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &Status);
//		}
//	}
//	MPI_Finalize();
//	return 0;
//}