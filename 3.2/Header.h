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


double F(double x) {
	return (2 - sin(M_PI * x));
}
double K(double x, double s) {
	return 1 / (2 + sin(M_PI * (x + s)));
}

double MMK_F(double a, double b) {
	int N = 11;
	double h = (b - a) / 11;
	double** B = new double* [N - 1];
	for (int i = 0; i < N - 1; i++) {
		B[i] = new double[N - 1];
	}
	double x_i, x_j;
	for (int i = 0; i < N - 1; i++) {
		for (int j = 0; j < N - 1; j++) {
			x_i = a + i * h;
			x_j = a + j * h;
			if (i == j) {
				B[i][j] = 1 - 0.5 * h * K(x_i, x_j);
			}
			else {
				B[i][j] = -0.5 * h * K(x_i, x_j);
			}
		}
	}
	double* F_ = new double[N - 1];
	for (int i = 0; i < N - 1; i++) {
		x_i = a + i * h;
		F_[i] = F(x_i);
	}
	Matr M(N - 1, B);
	Sist S(M, F_);
	double* Y = S.Gauss();
	cout << "Значения функции u(x) в точках x_i:" << fixed << showpoint << setprecision(7) << endl;
	for (int i = 0; i < N - 1; i++) {
		cout << Y[i] << endl;
	}
	double x_ = (a + b) / 2.2;
	double U_x_ = 0;
	for (int i = 0; i < N - 1; i++) {
		x_i = a + i * h;
		U_x_ += K(x_, x_i) * Y[i];
	}
	U_x_ *= 0.5 * h;
	U_x_ += F(x_);
	return U_x_;
}
double MMK_V(double a, double b) {
	int N = 11;
	double h = (b - a) / 11;
	double* Y = new double[N];
	double x_i, x_k, A_i;
	for (int i = 0; i < N; i++) {
		Y[i] = 0;
		x_i = a + i * h;
		for (int k = 0; k < i; k++) {
			x_k = a + k * h;
			A_i = h;
			if (k == 0 || k == i) {
				A_i = h / 2;
			}
			Y[i] += 0.5 * A_i * K(x_i, x_k) * Y[k];
		}
		Y[i] += F(x_i);
		Y[i] /= 1 - (h / 2) * 0.5 * K(x_i, x_i);
	}
	cout << "Значения функции u(x) в точках x_i:" << endl;
	for (int i = 0; i < N; i++) {
		cout << Y[i] << endl;
	}
	double x_ = (a + b) / 2.2;
	double U_x_ = 0, S = 1;

	for (int i = 0; i < N; i++) {
		double x_i = a + i * h;
		S = Y[i];
		for (int j = 0; j < N; j++) {
			double x_j = a + j * h;
			if (i != j) {
				S *= (x_ - x_j);
				S /= (x_i - x_j);
			}
		}
		U_x_ += S;
	}

	return U_x_;
}

double MPP_F(double a, double b) {
	int N = 11, n = 5;
	double h = (b - a) / 11;
	double* Y_0 = new double[N - 1];
	double* Y_1 = new double[N - 1];
	double x_j, x_k, A_i;
	for (int i = 0; i < N - 1; i++) {
		Y_0[i] = F(a + i * h);
	}
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < N - 1; j++) {
			x_j = a + j * h;
			Y_1[j] = F(x_j);
			for (int k = 0; k < N - 1; k++) {
				x_k = a + k * h;
				Y_1[j] += 0.5 * h * K(x_j, x_k) * Y_0[k];
			}
		}
		for (int j = 0; j < N - 1; j++) {
			Y_0[j] = Y_1[j];
		}
	}
	cout << "Значения функции u(x) в точках x_i:" << endl;
	for (int i = 0; i < N - 1; i++) {
		cout << Y_0[i] << endl;
	}

	double x_ = (a + b) / 2.2;
	double U_x_ = 0;
	for (int i = 0; i < N - 1; i++) {
		double x_i = a + i * h;
		U_x_ += K(x_, x_i) * Y_0[i];
	}
	U_x_ *= 0.5 * h;
	U_x_ += F(x_);
	return U_x_;
}

double MPP_V(double a, double b) {
	int N = 11, n = 5;
	double h = (b - a) / 11;
	double* Y_0 = new double[N];
	double* Y_1 = new double[N];
	double x_j, x_k, A_i;
	for (int i = 0; i < N; i++) {
		Y_0[i] = F(a + i * h);
	}
	for (int i = 0; i < n - 1; i++) {
		for (int j = 0; j < N; j++) {
			x_j = a + j * h;
			Y_1[j] = F(x_j);
			for (int k = 0; k <= j; k++) {
				x_k = a + k * h;
				A_i = h;
				if (k == 0 || k == j) {
					A_i = h / 2;
				}
				Y_1[j] += 0.5 * A_i * K(x_j, x_k) * Y_0[k];
			}
		}
		for (int j = 0; j < N; j++) {
			Y_0[j] = Y_1[j];
		}
	}
	cout << "Значения функции u(x) в точках x_i:" << endl;
	for (int i = 0; i < N; i++) {
		cout << Y_0[i] << endl;
	}

	double x_ = (a + b) / 2.2;
	double U_x_ = 0, S = 1;

	for (int i = 0; i < N; i++) {
		double x_i = a + i * h;
		S = Y_0[i];
		for (int j = 0; j < N; j++) {
			x_j = a + j * h;
			if (i != j) {
				S /= (x_i - x_j);
				S *= (x_ - x_j);
			}
		}
		U_x_ += S;
	}

	return U_x_;
}