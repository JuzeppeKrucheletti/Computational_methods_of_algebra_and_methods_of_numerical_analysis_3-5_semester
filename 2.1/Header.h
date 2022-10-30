#include <iostream>
#include <stdio.h>
#include <Math.h>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;
double E = pow(10, -7);
double F(double x) {
	double F = sin(x) - 2 * x * x + 0.5;
	return F;
}
double F_(double x) {
	double F = cos(x) - 4 * x;
	return F;
}
double F__(double x) {
	double F = (sin(x) + 4) * (-1);
	return F;
}

double f(double x) {
	return (-1) * sqrt((sin(x) + 0.5) / 2);
}
double f_(double x) {
	//return pow((sin(x) + 0.5) / 2,-0.5)*cos(x)*0.25;
	return (-1) * cos(x) / (2 * sqrt((2 * sin(x) + 1)));
}
double* Dihotomia(double a, double b) {
	double x = (a + b) / 2;
	cout << "k" << "\t    " << "a_k   " << "\t" << "b_k   " << "\t" << "   F(a_k)   " << "\t" << "F(b_k)" << "\t         " << "(a_k + b_k) / 2" << "\t  " << " | a_k - b_k|" << endl;
	for (int k = 0; abs(a - b) > 0.2; k++) {
		cout << k << "\t   " << a << "\t       " << b << "\t   " << F(a) << "\t   " << F(b) << "\t        " << (a + b) / 2 << "\t                     " << abs(b - a) << endl;
		if (F(x) * F(a) < 0) b = x;
		else a = x;
		x = (a + b) / 2;
	}
	double* mass = new double[3];
	mass[0] = x;
	mass[1] = a;
	mass[2] = b;
	return mass;
}
bool T1(double x, double a, double b) {
	//1
	double delta = abs(a - b) / 2;
	cout << "1. Дельта = " << delta << endl;
	//2
	//найдём максимум производной фуекции мпи
	//заметим, что производная функции монотонна, значит можем найти максимум,
	//подставив значение на концах
	double q;
	if (abs(f_(a)) > abs(f_(b))) q = abs(f_(a));
	else q = abs(f_(b));
	cout << "2. q = " << q << endl;
	if (q >= 1) {
		cout << "Условия теоремы не выполнены: q>=1" << endl;
		return false;
	}
	//3
	double m = abs(x - f(x));
	cout << "3. m = " << m << endl;
	if ((m / (1 - q)) > delta) {
		cout << "Условия теоремы не выполнены: m / (1 - q) > delta" << endl;
		return false;
	}
	cout << "Условия теоремы выполнены." << endl;
	return true;
}
double MPI(double x, double a, double b) {

	cout << "k" << "\t   " << "Xk" << "\t          " << "|Xk+1-Xk|" << endl;
	double x_ = x;
	int k = 0;
	do
	{
		x_ = x;
		x = f(x);
		cout << dec << k << "\t " << scientific << x_ << "\t " << abs(x - x_) << endl;
		k++;
	} while (abs(x - x_) > E);
	return x;
}
bool T2(double x, double a, double b) {
	//1
	double h = (-1) * F(x) / F_(x);
	cout << "1. h0 = " << h << endl;
	cout << "Значения F(x)*F_(x) на концах отрезка: " << F(x) * F_(x) << " и " << F(x + 2 * h) * F_(x + 2 * h) << endl;
	if (abs(F(x) * F_(x)) < 0.00001 || abs(F(x + 2 * h) * F_(x + 2 * h)) < 0.00001) {
		cout << "Условия теоремы не выполнены: значения на концах отрезка близки к 0." << endl;
		return false;
	}
	//2
	double M;
	if (abs(F__(x)) > abs(F__(x + 2 * h))) M = abs(F__(x));
	else M = abs(F__(x + 2 * h));
	cout << "2. M = " << M << endl;
	if (2 * abs(h) * M > abs(F_(x))) {
		cout << "Условия теоремы не выполнены: 2*|h|*M > |F_(x)|" << endl;
		return false;
	}
	//3
	cout << "Условия теоремы выполнены." << endl;
	return true;
}
double MNpost(double x, double a, double b) {


	cout << "k" << "\t   " << "Xk" << "\t          " << "|Xk+1-Xk|" << endl;
	double x_ = x, F_post = F_(x);
	int k = 0;
	do
	{
		x_ = x;
		x = x - (F(x) / F_post);
		cout << dec << k << "\t " << scientific << x_ << "\t " << abs(x - x_) << endl;
		k++;
	} while (abs(x - x_) > E);
	return x;
}
double MN(double x, double a, double b) {


	cout << "k" << "\t   " << "Xk" << "\t          " << "|Xk+1-Xk|" << endl;
	double x_ = x;
	int k = 0;
	do
	{
		x_ = x;
		x = x - (F(x) / F_(x));
		cout << dec << k << "\t " << scientific << x_ << "\t " << abs(x - x_) << endl;
		k++;
	} while (abs(x - x_) > E);
	return x;
}
