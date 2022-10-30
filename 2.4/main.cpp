#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "Header.h"
#include <fstream>
#include<time.h>
using namespace std;
int main()
{
	setlocale(LC_ALL, "ru");
	cout << "»нтерпал€ционный многочлен в узлах „ебышЄва:" << endl;
	Cheb1(5);
	Cheb1(10);
	Cheb1(15);
	Cheb1(20);
	cout << endl;
	Cheb2(5);
	Cheb2(10);
	Cheb2(15);
	Cheb2(20);

	cout << endl << "»нтерпал€ционный многочлен в равносто€щих узлах:" << endl << endl;
	Ravn1(5);
	Ravn1(10);
	Ravn1(15);
	Ravn1(20);
	cout << endl;
	Ravn2(5);
	Ravn2(10);
	Ravn2(15);
	Ravn2(20);
	system("pause");
	return 0;
}
