#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <iostream>
#include "Header.h"
#include <fstream>
#include<time.h>
using namespace std;
int main()
{
	setlocale(LC_ALL, "ru");
	cout << "Явный метод Эйлера:" << endl;
	double* Y_me = JaME();
	cout << endl << "Метод Рунге-Кутты:" << endl;
	double* Y_mrk = MRCute();
	cout << endl << "Метод последовательного повышения порядка точности:" << endl;
	double* Y_mpppt = MPPPT(Y_me);
	cout << endl << "Метод Адамса:" << endl;
	double* Y_ma = NMA(Y_mrk[1]);
	cout << endl << "Невязка для явного метода Эйлера:" << endl;
	for (int j = 0; j < N; j++) {
		cout << "Невязка в точке " << fixed << showpoint << setprecision(2) << a + (j + 1) * tau << ":  " << scientific << abs(Y_me[j+1]-Y_ma[j+1]) << endl;
	}

	cout << endl << "Невязка для метода Рунге-Кутты:" << endl;
	for (int j = 0; j < N; j++) {
		cout << "Невязка в точке " << fixed << showpoint << setprecision(2) << a + (j + 1) * tau << ":  " << scientific << abs(Y_mrk[j + 1] - Y_ma[j + 1]) << endl;
	}

	cout << endl << "Невязка для метода ПППТ:" << endl;
	for (int j = 0; j < N; j++) {
		cout << "Невязка в точке " << fixed << showpoint << setprecision(2) << a + (j + 1) * tau << ":  " << scientific << abs(Y_mpppt[j + 1] - Y_ma[j + 1]) << endl;
	}
	system("pause");
	return 0;
}
