#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include<cmath>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;
const double E = pow(10, -5);
const double e = M_E;
const double M1 = 97.12;
const double M2 = 1.2 * pow(10, 9);
const double I_t = 2.084776558;

double F(double x) {
	return (pow(e, 2 * x) / sqrt(1 - x * x));
}
void LP_Runge(double a, double b) {
	int N = 2;
	double h = (b - a) / N;
	double Q_h, Q_h_2;
	double x_i;
	do {
		Q_h = 0;
		x_i = a;
		for (int i = 0; i < N; i++) {
			x_i = a + h * i;
			Q_h += F(x_i);
		}
		Q_h *= h;

		h /= 2;
		N *= 2;
		Q_h_2 = 0;
		x_i = a;
		for (int i = 0; i < N; i++) {
			x_i = a + h * i;
			Q_h_2 += F(x_i);
		}
		Q_h_2 *= h;
	} while (abs(Q_h_2 - Q_h) > E);

	cout << "Полученное приближенное значение:\nI = " << fixed << showpoint << setprecision(9) << Q_h_2 << endl;
	cout << "Значение шага h: " << scientific << h << endl;
	cout << "Количество разбиений N: " << N << endl;
	cout << "Невязка: " << scientific << abs(Q_h_2 - I_t) << endl;
}
void Trap(double a, double b) {
	int N = 1184;
	double h = (b - a) / N;
	double I = 0;
	for (int i = 1; i < N; i++) {

		I += F(a + i * h);
	}
	I += (F(a) + F(b)) / 2;
	I *= h;
	cout << "Полученное приближенное значение:\nI = " << fixed << showpoint << setprecision(9) << I << endl;
	cout << "Значение шага h: " << scientific << h << endl;
	cout << "Количество разбиений N: " << N << endl;
	cout << "Невязка: " << scientific << abs(I_t - I) << endl;
}
void SP(double a, double b) {
	int N = 837;
	double h = (b - a) / N;
	double I = 0;
	for (int i = 0; i < N; i++) {

		I += F(a + i * h + h / 2);
	}
	I *= h;
	cout << "Полученное приближенное значение:\nI = " << fixed << showpoint << setprecision(9) << I << endl;
	cout << "Значение шага h: " << scientific << h << endl;
	cout << "Количество разбиений N: " << N << endl;
	cout << "Невязка: " << scientific << abs(I_t - I) << endl;
}
void NAST(double a, double b) {
	int N = 4;
	double* X = new double[4];
	double* A = new double[4];
	X[0] = -0.8611363115;
	X[1] = -0.3399810435;
	X[2] = -X[1];
	X[3] = -X[0];
	A[0] = 0.3468548451;
	A[1] = 0.6521451548;
	A[2] = A[1];
	A[3] = A[0];
	double I = 0;
	for (int i = 0; i < N; i++) {
		double x_i = (a + b + (b - a) * X[i]) / 2;
		I += A[i] * F(x_i);
	}
	I *= (b - a) / 2;
	cout << "Полученное приближенное значение:\nI = " << fixed << showpoint << setprecision(9) << I << endl;
	cout << "Количество разбиений N: " << N << endl;
	cout << "Невязка: " << scientific << abs(I_t - I) << endl;
}
