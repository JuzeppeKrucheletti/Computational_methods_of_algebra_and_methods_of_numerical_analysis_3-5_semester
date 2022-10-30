#pragma once
#define _USE_MATH_DEFINES
#include<cmath>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;
double E = pow(10, -6);
double PI = acos(-1.0);


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





double F1(double x) {
	return (x * cos(x + 5));
}
double F2(double x) {
	return (1 / (1 + 25 * x * x));
}
//double Cheb(double x) {
//	return (x/5);
//}
void Cheb1(int N) {
	double* X_k = new double[N + 1];
	double* F_k = new double[N + 1];
	double a = -5, b = 5;
	//cout << "Таблица значений:" << endl << "X_k" << "           " << "F_k" << endl;
	for (int i = 0; i < N + 1; i++) {
		X_k[i] = (a + b) / 2 + ((b - a) / 2) * cos(((2 * i + 1) * PI) / (2 * (N + 1)));
		F_k[i] = F1(X_k[i]);
		//cout <<X_k[i]<<"      "<< endl;
	}






	double* Rez1 = new double[N + 1];
	cout << "Интерполяционный многочлен " << N << " степени для функции f1:" << endl;

	for (int i = 0; i < N + 1; i++) {
		Rez1[i] = F_k[i];
		for (int j = 0; j < N + 1; j++) {
			if (i != j) {
				Rez1[i] /= (X_k[i] - X_k[j]);
			}
		}
	}
	/*for (int i = 0; i < N + 1; i++) {
		if (Rez1[i] < 0) cout << fixed << showpoint << setprecision(20) << Rez1[i];
		else cout << "+" << Rez1[i];
		for (int j = 0; j < N + 1; j++) {
			if (i != j) {
				if (X_k[j] < 0) cout << "(x+" << abs(X_k[j]) << ")";
				else cout << "(x-" << X_k[j] << ")";
			}
		}
	}*/

	/*cout <<endl<< "Оценка погрешности:" <<endl;
	double max = 0;
	for (int i = 0; i < 100; i++) {
		int x_k = a + i * (b - a) / 100;
		double Pn_xk = 0;
		for (int j = 0; j < N + 1; j++) {
			Pn_xk += Rez1[j] * pow(x_k, j);
		}
		if (abs(F1(x_k) - Pn_xk) > max) max = abs(F1(x_k) - Pn_xk);
	}

	cout << max << endl;
	cout << endl;*/
	cout << "Оценка погрешности:" << endl;
	double max = 0;
	for (int k = 0; k < 100; k++) {
		double x_k = a + k * (b - a) / 100;
		double summ = 0;
		for (int i = 0; i < N + 1; i++) {
			double el = Rez1[i];
			for (int j = 0; j < N + 1; j++) {
				if (i != j) {
					el *= (x_k - X_k[j]);
				}
			}
			summ += el;
		}
		if (abs(summ - F1(x_k)) > max) max = abs(summ - F1(x_k));
	}
	cout << max << endl;
}
void Cheb2(int N) {
	double* X_k = new double[N + 1];
	double* F_k = new double[N + 1];
	double a = -5, b = 5;
	//cout << "Таблица значений:" << endl << "X_k" << "           " << "F_k" << endl;
	for (int i = 0; i < N + 1; i++) {
		X_k[i] = (a + b) / 2 + ((b - a) / 2) * cos(((2 * i + 1) * PI) / (2 * (N + 1)));
		F_k[i] = F2(X_k[i]);
		//cout <<X_k[i]<<"      "<< endl;
	}






	double* Rez1 = new double[N + 1];
	cout << "Интерполяционный многочлен " << N << " степени для функции f2:" << endl;

	for (int i = 0; i < N + 1; i++) {
		Rez1[i] = F_k[i];
		for (int j = 0; j < N + 1; j++) {
			if (i != j) {
				Rez1[i] /= (X_k[i] - X_k[j]);
			}
		}
	}
	/*for (int i = 0; i < N + 1; i++) {
		if (Rez1[i] < 0) cout << fixed << showpoint << setprecision(20) << Rez1[i];
		else cout << "+" << Rez1[i];
		for (int j = 0; j < N + 1; j++) {
			if (i != j) {
				if (X_k[j] < 0) cout << "(x+" << abs(X_k[j]) << ")";
				else cout << "(x-" << X_k[j] << ")";
			}
		}
	}*/

	/*cout <<endl<< "Оценка погрешности:" <<endl;
	double max = 0;
	for (int i = 0; i < 100; i++) {
		int x_k = a + i * (b - a) / 100;
		double Pn_xk = 0;
		for (int j = 0; j < N + 1; j++) {
			Pn_xk += Rez1[j] * pow(x_k, j);
		}
		if (abs(F1(x_k) - Pn_xk) > max) max = abs(F1(x_k) - Pn_xk);
	}

	cout << max << endl;
	cout << endl;*/
	cout << "Оценка погрешности:" << endl;
	double max = 0;
	for (int k = 0; k < 100; k++) {
		double x_k = a + k * (b - a) / 100;
		double summ = 0;
		for (int i = 0; i < N + 1; i++) {
			double el = Rez1[i];
			for (int j = 0; j < N + 1; j++) {
				if (i != j) {
					el *= (x_k - X_k[j]);
				}
			}
			summ += el;
		}
		if (abs(summ - F2(x_k)) > max) max = abs(summ - F2(x_k));
	}
	cout << max << endl;
}
void Ravn1(int N) {
	double* X_k = new double[N + 1];
	double* F_k = new double[N + 1];
	double a = -5, b = 5;
	//cout << "Таблица значений:" << endl << "X_k" << "           " << "F_k" << endl;
	for (int i = 0; i < N + 1; i++) {
		X_k[i] = a + i * (b - a) / N;
		F_k[i] = F1(X_k[i]);
		//cout <<X_k[i]<<"      "<< endl;
	}






	double* Rez1 = new double[N + 1];
	cout << "Интерполяционный многочлен " << N << " степени для функции f1:" << endl;

	for (int i = 0; i < N + 1; i++) {
		Rez1[i] = F_k[i];
		for (int j = 0; j < N + 1; j++) {
			if (i != j) {
				Rez1[i] /= (X_k[i] - X_k[j]);
			}
		}
	}
	/*for (int i = 0; i < N + 1; i++) {
		if (Rez1[i] < 0) cout << fixed << showpoint << setprecision(20) << Rez1[i];
		else cout << "+" << Rez1[i];
		for (int j = 0; j < N + 1; j++) {
			if (i != j) {
				if (X_k[j] < 0) cout << "(x+" << abs(X_k[j]) << ")";
				else cout << "(x-" << X_k[j] << ")";
			}
		}
	}*/

	/*cout <<endl<< "Оценка погрешности:" <<endl;
	double max = 0;
	for (int i = 0; i < 100; i++) {
		int x_k = a + i * (b - a) / 100;
		double Pn_xk = 0;
		for (int j = 0; j < N + 1; j++) {
			Pn_xk += Rez1[j] * pow(x_k, j);
		}
		if (abs(F1(x_k) - Pn_xk) > max) max = abs(F1(x_k) - Pn_xk);
	}

	cout << max << endl;
	cout << endl;*/
	cout << "Оценка погрешности:" << endl;
	double max = 0;
	for (int k = 0; k < 100; k++) {
		double x_k = a + k * (b - a) / 100;
		double summ = 0;
		for (int i = 0; i < N + 1; i++) {
			double el = Rez1[i];
			for (int j = 0; j < N + 1; j++) {
				if (i != j) {
					el *= (x_k - X_k[j]);
				}
			}
			summ += el;
		}
		if (abs(summ - F1(x_k)) > max) max = abs(summ - F1(x_k));
	}
	cout << max << endl;
}



void Ravn2(int N) {
	double* X_k = new double[N + 1];
	double* F_k = new double[N + 1];
	double a = -5, b = 5;
	//cout << "Таблица значений:" << endl << "X_k" << "           " << "F_k" << endl;
	for (int i = 0; i < N + 1; i++) {
		X_k[i] = a + i * (b - a) / N;
		F_k[i] = F2(X_k[i]);
		//cout <<X_k[i]<<"      "<< endl;
	}






	double* Rez1 = new double[N + 1];
	cout << "Интерполяционный многочлен " << N << " степени для функции f2:" << endl;

	for (int i = 0; i < N + 1; i++) {
		Rez1[i] = F_k[i];
		for (int j = 0; j < N + 1; j++) {
			if (i != j) {
				Rez1[i] /= (X_k[i] - X_k[j]);
			}
		}
	}
	/*for (int i = 0; i < N + 1; i++) {
		if (Rez1[i] < 0) cout << fixed << showpoint << setprecision(20) << Rez1[i];
		else cout << "+" << Rez1[i];
		for (int j = 0; j < N + 1; j++) {
			if (i != j) {
				if (X_k[j] < 0) cout << "(x+" << abs(X_k[j]) << ")";
				else cout << "(x-" << X_k[j] << ")";
			}
		}
	}*/
	/*cout <<endl<< "Оценка погрешности:" <<endl;
	double max = 0;
	for (int i = 0; i < 100; i++) {
		int x_k = a + i * (b - a) / 100;
		double Pn_xk = 0;
		for (int j = 0; j < N + 1; j++) {
			Pn_xk += Rez1[j] * pow(x_k, j);
		}
		if (abs(F1(x_k) - Pn_xk) > max) max = abs(F1(x_k) - Pn_xk);
	}

	cout << max << endl;
	cout << endl;*/
	cout << "Оценка погрешности:" << endl;
	double max = 0;
	for (int k = 0; k < 100; k++) {
		double x_k = a + k * (b - a) / 100;
		double summ = 0;
		for (int i = 0; i < N + 1; i++) {
			double el = Rez1[i];
			for (int j = 0; j < N + 1; j++) {
				if (i != j) {
					el *= (x_k - X_k[j]);
				}
			}
			summ += el;
		}
		if (abs(summ - F2(x_k)) > max) max = abs(summ - F2(x_k));
	}
	cout << max << endl;
}

