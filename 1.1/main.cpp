#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "Header.h"
#include <fstream>
#include<time.h>
using namespace std;
ifstream fin;
ofstream fout;
int main()
{
	srand(time(0));
	setlocale(LC_ALL, "ru");
	Matr M(10), N(M), L(N);
	Sist A(M), B(A);
	//1
	cout << "������� � �������: " << endl << M << endl;
	A.Gauss();
	cout << "������� �X = f ����� ������� ������� ������: " << endl;
	cout << A << endl;
	cout << endl;
	//2
	cout << "������� �:" << endl;
	cout << M << endl;
	N = M.Obrat();
	cout << "������� N, �������� �:" << endl;
	cout << N << endl;
	cout << "������� M*N:" << endl;
	L = M * N;
	cout << L << endl;
	cout << endl;
	cout << "������� �������:           ������������� �����������:" << endl;
	for (int i = 6; i <= 106; i += 10)
	{
		Sist A1(i);
		A1.Gauss();
		cout << i << "                          ";
		A1.ShowPogr();
		cout << endl;
	}
	system("pause");
	return 0;
}
