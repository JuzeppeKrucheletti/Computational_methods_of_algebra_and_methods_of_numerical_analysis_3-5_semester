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
	cout << "���������:" << endl;
	double* mass = Dihotomia(-1, 0);
	cout << endl << "�������� ������� ������� � ���������� ���:" << endl;
	T1(mass[0], mass[1], mass[2]);
	cout << endl;
	cout << "���:" << endl;
	MPI(mass[0], mass[1], mass[2]);
	cout << endl << "�������� ������� ������� � ���������� ������ �������:" << endl;
	T2(mass[0], mass[1], mass[2]);
	cout << endl << "����� �������:" << endl;
	MN(mass[0], mass[1], mass[2]);
	cout << endl << "����� ������� c ���������� �����������:" << endl;
	MNpost(mass[0], mass[1], mass[2]);
	return 0;
}
