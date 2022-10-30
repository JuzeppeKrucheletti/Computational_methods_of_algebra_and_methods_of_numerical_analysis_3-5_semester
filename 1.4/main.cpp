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
	Matr M(5);
	cout << M << endl;
	cout << endl << "Степенной метод:" << endl;
	M.StepennoiMetod();
	cout << endl << "Итерационный метод вращений:" << endl;
	M.MetodVrashchenij();
	system("pause");
	return 0;
}