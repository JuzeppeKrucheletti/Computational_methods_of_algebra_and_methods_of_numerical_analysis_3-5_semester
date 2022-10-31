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
	cout << "Метод механических квадратур для ИУФ:" << endl;
	double mmk_f = MMK_F(0, 1);
	cout << endl << "Метод механических квадратур для ИУВ:" << endl;
	double mmk_v = MMK_V(0, 1);
	cout << endl << "Метод последовательных приближений для ИУФ:" << endl;
	double mpp_f = MPP_F(0, 1);
	cout << endl << "Метод последовательных приближений для ИУВ:" << endl;
	double mpp_v = MPP_V(0, 1);

	cout << endl << "Значения в точке для ИУФ:" << endl;
	cout << "ММК: " << mmk_f << endl;
	cout << "МПП: " << mpp_f << endl;
	cout << endl << "Значения в точке для ИУВ:" << endl;
	cout << "ММК: " << mmk_v << endl;
	cout << "МПП: " << mpp_v << endl;
	system("pause");
	return 0;
}
