#include <iostream>
#include <vector>
#include <mpi.h>
#include <time.h> 

using namespace std;

struct matrix
{
	int length;
	int high;
	vector<int> value;
};

void MatrixWithVector(matrix* mA, matrix* mB, int place, matrix* mC) {
	int Csize = mA->value.size();
	for (int i = 0; i < Csize; i++) {
		if (i % mA->length == place) {
			mC->value.push_back(mB->value[i / mA->length]);
		}
		else {
			mC->value.push_back(mA->value[i]);
		}
	}
}

void MatrixDelColRow(matrix* mA, int col, int row) {
	for (int i = col; i < mA->value.size(); i += mA->length - 1) {
		mA->value.erase(mA->value.begin() + i);
	}
	mA->length--;
	for (int j = row * (mA->length) + mA->length - 1; j >= row * mA->length; j--) {
		mA->value.erase(mA->value.begin() + j);
	}
	mA->high--;
}

int matrixDet(matrix* mA) {
	int det = 0;
	int degree = 1;
	if (mA->value.size() == 1) {
		return mA->value[0];
	}
	for (int j = 0; j < mA->high; j++) {
		matrix mB;
		mB = *mA;
		MatrixDelColRow(&mB, 0, j);
		det += (degree * mA->value[j * mA->length] * matrixDet(&mB));
		degree *= -1;
	}
	return det;
}

int main(int argc, char** argv)
{
	clock_t start;
	int size, rank;
	int n;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
		cout << "Введите число n уравнений и неизвестных:";
		cin >> n;

		cout << "Введите последовательно коэффициенты при неизвестных:";
		int* mA = new int[n * n];
		for (int i = 0; i < n * n; i++) {
			int val;
			cin >> val;
			mA[i] = val;
		}

		cout << "Введите последовательно свободные члены";
		int* mB = new int[n];
		for (int i = 0; i < n; i++) {
			int val;
			cin >> val;
			mB[i] = val;
		}
		start = clock();
		for (int to_thread = 1; to_thread < size; to_thread++) {
			MPI_Send(&n, 1, MPI_INT, to_thread, 0, MPI_COMM_WORLD);
			MPI_Send(mA, n * n, MPI_INT, to_thread, 1, MPI_COMM_WORLD);
			MPI_Send(mB, n, MPI_INT, to_thread, 2, MPI_COMM_WORLD);
		}
	}
	else
	{
		matrix mA, mB;
		MPI_Status status1, status2, status3;
		int count, count2;
		MPI_Recv(&n, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status1);

		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status2);
		MPI_Get_count(&status2, MPI_INT, &count);
		int* buf = new int[count];
		MPI_Recv(buf, count, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status2);
		for (int i = 0; i < n * n; i++) {
			mA.value.push_back(buf[i]);
		}
		mA.length = n;
		mA.high = n;

		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status3);
		MPI_Get_count(&status3, MPI_INT, &count2);
		int* buf2 = new int[count2];
		MPI_Recv(buf2, count2, MPI_INT, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status3);
		for (int i = 0; i < n; i++) {
			mB.value.push_back(buf2[i]);
		}
		mB.length = 1;
		mB.high = n;

		int range = (n - 1) / (size - 1) + 1;
		int ibeg = (rank - 1) * range;
		int iend = rank * range;

		if (rank == 1) {
			int detA = matrixDet(&mA);
			if (detA == 0) {
				cout << "Система имеет бесконечно много решений или ни одного";
				return 0;
			}
			MPI_Send(&detA, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
		}

		for (int i = ibeg; i < n && i < iend; i++) {
			matrix mC;
			mC.length = n;
			mC.high = n;
			MatrixWithVector(&mA, &mB, i, &mC);
			int det = matrixDet(&mC);
			MPI_Send(&det, 1, MPI_INT, 0, i + 3, MPI_COMM_WORLD);
		}

	}
	if (rank == 0)
	{
		int detA;
		MPI_Status status1;
		MPI_Recv(&detA, 1, MPI_INT, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, &status1);

		int* arrayDet = new int[n];
		for (int i = 0; i < n; i++) {
			MPI_Status status;
			int det;
			MPI_Recv(&det, 1, MPI_INT, MPI_ANY_SOURCE, i + 3, MPI_COMM_WORLD, &status);
			arrayDet[i] = det;
		}

		for (int i = 0; i < n; i++) {
			cout << i + 1 << " Неизвестная:" << arrayDet[i] * 1.0 / detA << "\n";
		}

		clock_t end = clock();
		double seconds = (double)(end - start) / CLOCKS_PER_SEC;
		cout << "Время выполнения: " << seconds;
	}

	MPI_Finalize();
	return 0;
}