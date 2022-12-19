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



// Гаусс
class Matr
{
	friend class Sist;

	double** M;
	int n;
public:
	Matr()
	{
		this->n = 1;
		M = new double* [n];
		for (int i = 0; i < n; i++)
		{
			M[i] = new double[n];
		}
		srand(time(0));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				M[i][j] = (((double)rand() / RAND_MAX) * 2 - 1) * 100;
			}
		}
	}
	Matr(int n, double** M1)
	{
		this->n = n;
		M = new double* [n];
		for (int i = 0; i < n; i++)
		{
			M[i] = new double[n];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				M[i][j] = M1[i][j];
			}
		}
	}
	Matr(const Matr& B)
	{
		this->n = B.n;
		M = new double* [n];
		for (int i = 0; i < n; i++)
		{
			M[i] = new double[n];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				M[i][j] = B.M[i][j];
			}
		}
	}
	~Matr()
	{
		for (int i = 0; i < n; i++)
		{
			delete[] M[i];
		}
		delete[] M;
	}
	friend ostream& operator <<(ostream& os, Matr& A)
	{
		os << fixed << showpoint << setprecision(3);
		for (int i = 0; i < A.n; i++)
		{

			for (int j = 0; j < A.n; j++)
			{
				os << A.M[i][j] << " ";
			}
			os << endl;
		}
		return os;
	}

	Matr& operator = (const Matr& B)
	{
		if (this == &B)
			return *this;
		for (int i = 0; i < n; i++)
		{
			delete[] M[i];
		}
		n = B.n;
		delete[] M;
		M = new double* [n];
		for (int i = 0; i < n; i++)
		{
			M[i] = new double[n];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				M[i][j] = B.M[i][j];
			}
		}
		return *this;
	}

};
class Sist
{
	friend class Matr;
private:
	Matr A;
	double* X;
	double* F;
	double* X_;
	double pogr;
public:
	int size()
	{
		return A.n;
	}
	void ShowPogr()
	{
		cout << scientific << pogr << fixed << showpoint << setprecision(3);
	}
	Sist(const Matr& B, double* F1)
	{
		A = B;
		X = new double[B.n];
		for (int i = 0; i < B.n; i++)
		{
			X[i] = 1 + i;
		}
		F = F1;
		X_ = new double[B.n];
		for (int i = 0; i < size(); i++)
		{
			X_[i] = 0;
		}
	}

	~Sist()
	{
		delete[] X;
		delete[] X_;
		delete[] F;
	}
	double* Gauss()
	{
		//выбор главного элемента по столбцу
		for (int i = 0; i < size(); i++)
		{
			double max = abs(A.M[i][i]);
			int gl = i;
			for (int j = i; j < size(); j++)
			{
				if (abs(A.M[j][i]) > max)
				{
					max = A.M[j][i];
					gl = j;
				}
			}
			if (max == 0) throw "Все главные элементы нулевые";
			//исключение в случае, когда все главные элементы нулевые
			for (int j = i; j < size(); j++)
			{
				swap(A.M[i][j], A.M[gl][j]);
			}
			swap(F[i], F[gl]);
			//прямой ход метода Гаусса
			for (int j = i + 1; j < size(); j++)
			{
				A.M[i][j] /= A.M[i][i];
			}
			F[i] /= A.M[i][i];
			A.M[i][i] /= A.M[i][i];
			for (int j = i + 1; j < size(); j++)
			{
				for (int k = i + 1; k < size(); k++)
				{
					A.M[j][k] -= A.M[i][k] * A.M[j][i];
				}
				F[j] -= F[i] * A.M[j][i];
				A.M[j][i] -= A.M[i][i] * A.M[j][i];
			}
		}
		//обратный ход метода Гаусса
		for (int i = (size() - 1); i >= 0; i--)
		{
			X_[i] = F[i];
			for (int j = size() - 1; j > i; j--)
			{
				X_[i] -= A.M[i][j] * X_[j];
			}
		}
		//вычисление относительной погрешности
		double maxX = abs(X[0]), max_X = abs(X[0] - X_[0]);
		for (int i = 0; i < size(); i++)
		{
			if (max_X < abs(X[i] - X_[i]))
			{
				max_X = abs(X[i] - X_[i]);
			}
			if (maxX < abs(X[i]))
			{
				maxX = abs(X[i]);
			}
		}
		pogr = max_X / maxX;
		return X_;
	}

};


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
