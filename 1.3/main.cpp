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
	Matr M(10);
	Sist A1(M), A2(M);
	cout << A1 << endl;
	A1.MMN();
	cout << endl;
	cout << "  W\t\t Количество итераций q\t\t Норма\t\t\t Абсолютная погрешность" << endl;
	for (double w = 0.2; w <= 1.8; w += 0.2)
	{
		A2.Relax(w);
		cout << w << "\t\t\t";
		A2.ShowQ();
		cout << "\t\t\t";
		A2.ShowNorm();
		cout << "\t\t\t";
		A2.ShowPogr();
		cout << endl;
	}
	system("pause");
	return 0;
}
