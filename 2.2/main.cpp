#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include "Header.h"
#include <fstream>
#include<time.h>
using namespace std;
int main()
{
	setlocale(LC_ALL, "ru");
	cout << "����� �������:" << endl;
	Newton(1, 2);
	cout << endl << "����� �������:" << endl;
	Sek(1, 2, 1.2, 2.4);
	cout << endl << "����� ������-�������:" << endl;
	Gauss_Zeidel(1, 2);
	system("pause");
	return 0;
}
