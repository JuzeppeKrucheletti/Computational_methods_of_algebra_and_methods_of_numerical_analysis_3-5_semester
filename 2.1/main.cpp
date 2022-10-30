#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <string>
#include "Header.h"
using namespace std;
ifstream fin;
ofstream fout;
long main()
{
	setlocale(LC_ALL, "ru");
	cout << "Дихотомия:" << endl;
	double* mass = Dihotomia(-1, 0);
	cout << endl << "Проверка условий теоремы о сходимости МПИ:" << endl;
	T1(mass[0], mass[1], mass[2]);
	cout << endl;
	cout << "МПИ:" << endl;
	MPI(mass[0], mass[1], mass[2]);
	cout << endl << "Проверка условий теоремы о сходимости метода Ньютона:" << endl;
	T2(mass[0], mass[1], mass[2]);
	cout << endl << "Метод Ньютона:" << endl;
	MN(mass[0], mass[1], mass[2]);
	cout << endl << "Метод Ньютона c постоянной производной:" << endl;
	MNpost(mass[0], mass[1], mass[2]);
	return 0;
}
