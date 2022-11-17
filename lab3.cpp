//#include <iostream>
//#include "mpi.h"
//#include <windows.h>
//#include <ctime>
//
//using namespace std;
//
//int main(int argc, char* argv[])
//{
//    int Procnum, ProcRank, count = 0, send;
//
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//    MPI_Comm_size(MPI_COMM_WORLD, &Procnum);
//
//    srand((unsigned int)time(0) + ProcRank);
//    bool flag = true;
//    int* mes = new int[Procnum];
//    while (flag)
//    {
//        static time_t tval = time(0);
//        tval += 10;
//        srand(tval);
//        send = rand() % 5 - 1;
//        MPI_Gather(&send, 1, MPI_INT,
//            mes, 1, MPI_INT, 0, MPI_COMM_WORLD);
//        if (ProcRank == 0)
//        {
//            for (int i = 0; i < Procnum; i++)
//            {
//                if (mes[i] == -1)
//                    flag = false;
//                if (flag) count++;
//            }
//        }
//        MPI_Bcast(&flag, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
//    }
//    if (ProcRank == 0)
//        printf("Count = %d\n", count);
//
//    MPI_Finalize();
//
//    delete[] mes;
//
//    return 0;
//}
