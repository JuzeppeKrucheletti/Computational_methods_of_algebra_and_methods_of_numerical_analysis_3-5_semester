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
	cout << "����� ������������ ��������� ��� ���:" << endl;
	double mmk_f = MMK_F(0, 1);
	cout << endl << "����� ������������ ��������� ��� ���:" << endl;
	double mmk_v = MMK_V(0, 1);
	cout << endl << "����� ���������������� ����������� ��� ���:" << endl;
	double mpp_f = MPP_F(0, 1);
	cout << endl << "����� ���������������� ����������� ��� ���:" << endl;
	double mpp_v = MPP_V(0, 1);

	cout << endl << "�������� � ����� ��� ���:" << endl;
	cout << "���: " << mmk_f << endl;
	cout << "���: " << mpp_f << endl;
	cout << endl << "�������� � ����� ��� ���:" << endl;
	cout << "���: " << mmk_v << endl;
	cout << "���: " << mpp_v << endl;
	system("pause");
	return 0;
}
