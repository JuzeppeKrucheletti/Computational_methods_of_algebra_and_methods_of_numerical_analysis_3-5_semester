#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include<cmath>
#include <iostream>
#include "Header.h"
#include <fstream>
#include<time.h>
using namespace std;
int main()
{
	setlocale(LC_ALL, "ru");
	cout << "1.\n�� ��:" << endl;
	LP_Runge(-0.5, 0.7);
	cout << "2.\n�� �:" << endl;
	Trap(-0.5, 0.7);
	cout << "�� C�:" << endl;
	SP(-0.5, 0.7);
	cout << "2.\n�� ���� ��� n=3:" << endl;
	NAST(-0.5, 0.7);
	system("pause");
	return 0;
}
