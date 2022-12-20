#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;
const double E = pow(10, -6);





double F(double t, double u) {
	return ((2*u*u*u*t)/(1-t*t*u*u));
}
const double u_0 = 1;
const double a = 2;
const double b = 2.4;
const double N = 10;
const double tau = (b - a) / N;

double* JaME() {
	double* Y = new double[N+1];
	Y[0] = u_0;
	for (int j = 0; j < N; j++) {
		Y[j + 1] = Y[j] + tau * F(a + tau * j, Y[j]);
		cout << "Значение в точке " << fixed << showpoint << setprecision(2) <<  a + (j + 1) * tau << ":  " << fixed << showpoint << setprecision(8) << Y[j + 1] << endl;
	}

	return Y;
}

double K_2(double t, double u) {
	return F(t + tau / 2, u + (tau / 2) * F(t, u));
}

double* MRCute() {
	double* Y = new double[N + 1];
	Y[0] = u_0;
	for (int j = 0; j < N; j++) {
		Y[j + 1] = Y[j] + tau * K_2(a + tau * j, Y[j]);
		cout << "Значение в точке " << fixed << showpoint << setprecision(2) << a + (j + 1) * tau << ":  " << fixed << showpoint << setprecision(8) << Y[j + 1] << endl;
	}

	return Y;
}
double* MPPPT(double* Y_JaME) {
	double* Y = new double[N + 1];
	Y[0] = u_0;
	for (int j = 0; j < N; j++) {
		Y[j + 1] = Y[j] + (tau / 2) * (F(a + tau * j, Y[j]) + F(a + tau * (j + 1), Y_JaME[j + 1]));
		cout << "Значение в точке " << fixed << showpoint << setprecision(2) << a + (j + 1) * tau << ":  " << fixed << showpoint << setprecision(8) << Y[j + 1] << endl;
	}

	return Y;
}

double* NMA(double Y_MRK) {
	double* Y = new double[N + 1];
	Y[0] = u_0;
	Y[1] = Y_MRK;
	cout << "Значение в точке " << fixed << showpoint << setprecision(2) << a + 1 * tau << ":  " << fixed << showpoint << setprecision(8) << Y[1] << endl;
	for (int j = 1; j < N; j++) {
		int k = 0;
		Y[j + 1] = Y[j];
		double pred_y;
		do {
			pred_y = Y[j + 1];
			Y[j + 1] = Y[j] + (tau / 12) * (5 * F(a + tau * (j + 1), pred_y) + 8 * F(a + tau * j, Y[j]) - F(a + tau * (j - 1), Y[j - 1]));
		} while (abs(Y[j + 1] - pred_y) > pow(tau, 5));
		cout << "Значение в точке " << fixed << showpoint << setprecision(2) << a + (j + 1) * tau << ":  " << fixed << showpoint << setprecision(8) << Y[j + 1] << endl;
	}

	return Y;
}
