#pragma once
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;
double E = pow(10, -6);



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





double F1(double x, double y) {
	return (sin(x + 0.5 * y) - 0.5 * (x - y) - 1.3);
}
double F2(double x, double y) {
	return (x - y + x * x + y * y - 6.35);
}
double DF1_Dx(double x, double y) {
	return (cos(x + 0.5 * y) - 0.5);
}
double DF1_Dy(double x, double y) {
	return (0.5 * cos(x + 0.5 * y) + 0.5);
}
double DF2_Dx(double x, double y) {
	return (1 + 2 * x);
}
double DF2_Dy(double x, double y) {
	return (2 * y - 1);
}
void Newton(double x_0, double y_0) {
	double x_k, y_k, x_k_1, y_k_1, delta_x_k, delta_y_k;
	x_k = x_0;
	y_k = y_0;
	cout << "k" << "     " << "  x_k" << "                " << "y_k" << "              " << " || delta_X || " << endl;
	cout << "0       " << x_0 << "                     " << y_0 << "                  " << "-" << endl;
	int k = 1;
	while (true) {
		double** M = new double* [2];
		M[0] = new double[2];
		M[1] = new double[2];
		M[0][0] = DF1_Dx(x_k, y_k);
		M[0][1] = DF1_Dy(x_k, y_k);
		M[1][0] = DF2_Dx(x_k, y_k);
		M[1][1] = DF2_Dy(x_k, y_k);
		double* F = new double[2];
		F[0] = (-1) * F1(x_k, y_k);
		F[1] = (-1) * F2(x_k, y_k);
		Matr A(2, M);
		Sist S(A, F);
		double* delta_X = S.Gauss();
		x_k_1 = x_k + delta_X[0];
		y_k_1 = y_k + delta_X[1];
		double norm = abs(delta_X[0]);
		if (abs(delta_X[1]) > norm) norm = abs(delta_X[1]);
		cout << dec << k << "    " << scientific << x_k_1 << "        " << y_k_1 << "        " << norm << dec << endl;
		x_k = x_k_1;
		y_k = y_k_1;
		k++;
		if (norm <= E) break;
	}
	double norm_f = F1(x_k_1, y_k_1);
	if (F2(x_k_1, y_k_1) > norm_f) norm_f = F2(x_k_1, y_k_1);
	cout << "||F(X_n)|| = " << norm_f << endl;
}
void Sek(double x_0, double y_0, double x_1, double y_1) {
	double x_k, y_k, x_k_1, y_k_1, x_k_2, y_k_2, delta_x_k, delta_y_k;
	x_k = x_0;
	y_k = y_0;
	x_k_1 = x_1;
	y_k_1 = y_1;
	cout << dec << "k" << "     " << "  x_k" << "                " << "y_k" << "              " << " || delta_X || " << endl;
	cout << "0    " << x_0 << "        " << y_0 << "          " << "-" << endl;
	double norm1 = abs(x_1 - x_0);
	if (abs(y_1 - y_0) > norm1) norm1 = abs(y_1 - y_0);
	cout << "1    " << x_1 << "        " << y_1 << "        " << norm1 << endl;
	int k = 2;
	while (true) {
		double** M = new double* [2];
		M[0] = new double[2];
		M[1] = new double[2];
		M[0][0] = (F1(x_k_1, y_k_1) - F1(x_k, y_k_1)) / (x_k_1 - x_k);
		M[0][1] = (F1(x_k_1, y_k_1) - F1(x_k_1, y_k)) / (y_k_1 - y_k);
		M[1][0] = (F2(x_k_1, y_k_1) - F2(x_k, y_k_1)) / (x_k_1 - x_k);
		M[1][1] = (F2(x_k_1, y_k_1) - F2(x_k_1, y_k)) / (y_k_1 - y_k);
		double* F = new double[2];
		F[0] = (-1) * F1(x_k_1, y_k_1);
		F[1] = (-1) * F2(x_k_1, y_k_1);
		Matr A(2, M);
		Sist S(A, F);
		double* delta_X = S.Gauss();
		x_k_2 = x_k_1 + delta_X[0];
		y_k_2 = y_k_1 + delta_X[1];
		double norm = abs(delta_X[0]);
		if (abs(delta_X[1]) > norm) norm = abs(delta_X[1]);
		cout << dec << k << "    " << scientific << x_k_2 << "        " << y_k_2 << "        " << norm << dec << endl;
		x_k = x_k_1;
		y_k = y_k_1;
		x_k_1 = x_k_2;
		y_k_1 = y_k_2;
		k++;
		if (norm <= E) break;
	}
	double norm_f = F1(x_k_2, y_k_2);
	if (F2(x_k_2, y_k_2) > norm_f) norm_f = F2(x_k_2, y_k_2);
	cout << "||F(X_n)|| = " << norm_f << endl;
}



void Gauss_Zeidel(double x_0, double y_0) {
	double x_k, y_k, x_k_1, y_k_1, norm_f;
	x_k = x_0;
	y_k = y_0;
	cout << dec << "k" << "     " << "  x_k" << "                " << "y_k" << "              " << " ||F(X_k)|| " << endl;
	double norm_f1 = F1(x_k, y_k);
	if (F2(x_k, y_k) > norm_f1) norm_f1 = F2(x_k, y_k);
	cout << "0    " << x_0 << "        " << y_0 << "          " << norm_f1 << endl;
	int k = 1;
	while (true) {
		//1
		double x_s, x_s_1;
		x_s = x_k;
		while (true)
		{
			x_s_1 = x_s - (F1(x_s, y_k) / DF1_Dx(x_s, y_k));
			if (abs(x_s - x_s_1) < E) break;
			x_s = x_s_1;
		}
		x_k_1 = x_s_1;
		//
		// 
		//2
		double y_s, y_s_1;
		y_s = y_k;
		while (true)
		{
			y_s_1 = y_s - (F2(x_k_1, y_s) / DF2_Dy(x_k_1, y_s));
			if (abs(y_s - y_s_1) < E) break;
			y_s = y_s_1;
		}
		y_k_1 = y_s_1;
		//
		x_k = x_k_1;
		y_k = y_k_1;
		norm_f = abs(F1(x_k_1, y_k_1));
		if (abs(F2(x_k_1, y_k_1)) > norm_f) norm_f = abs(F2(x_k_1, y_k_1));
		cout << dec << k << "    " << scientific << x_k_1 << "        " << y_k_1 << "        " << norm_f << dec << endl;
		k++;

		if (norm_f <= E) break;
	}
}
