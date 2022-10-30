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
	TD_Sist A1(10);
	cout << A1 << endl;
	A1.Progonka();
	cout << A1 << endl;
	system("pause");
	return 0;
}
