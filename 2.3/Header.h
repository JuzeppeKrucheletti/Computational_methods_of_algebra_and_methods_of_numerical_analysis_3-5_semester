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





double F(double x) {
	return (x * cos(x + 5));
}
void MNK() {
	int N = 60, N1 = 5, N2 = 7;
	double* X_k = new double[N + 1];
	double* F_k = new double[N + 1];
	double x = -3;
	//cout << "Таблица значений:" << endl << "X_k" << "           " << "F_k" << endl;
	for (int i = 0; i < N + 1; x = x + 0.1, i++) {
		X_k[i] = x;
		F_k[i] = F(x);
		//cout <<X_k[i]<<"      "<< endl;
	}

	double* S1 = new double[N1 * 2];
	double* M1 = new double[N1];
	double* S2 = new double[N2 * 2];
	double* M2 = new double[N2];
	//4
	for (int i = 0; i < N1; i++) {
		M1[i] = 0;
		for (int j = 0; j < N + 1; j++) {
			M1[i] += F_k[j] * pow(X_k[j], i);
		}

	}
	for (int i = 0; i < N1 * 2; i++) {
		S1[i] = 0;
		for (int j = 0; j < N + 1; j++) {
			S1[i] += pow(X_k[j], i);
		}
	}
	//6
	for (int i = 0; i < N2; i++) {
		M2[i] = 0;
		for (int j = 0; j < N + 1; j++) {
			M2[i] += F_k[j] * pow(X_k[j], i);
		}
	}
	for (int i = 0; i < N2 * 2; i++) {
		S2[i] = 0;
		for (int j = 0; j < N + 1; j++) {
			S2[i] += pow(X_k[j], i);
		}
	}
	double** Matr1 = new double* [N1];
	for (int i = 0; i < N1; i++) {
		Matr1[i] = new double[N1];
	}
	for (int i = 0; i < N1; i++) {
		for (int j = 0; j < N1; j++) {
			Matr1[i][j] = S1[i + j];
			//cout << i + j << " ";

		}
		//cout << endl;
	}

	double** Matr2 = new double* [N2];
	for (int i = 0; i < N2; i++) {
		Matr2[i] = new double[N2];
	}
	for (int i = 0; i < N2; i++) {
		for (int j = 0; j < N2; j++) {
			Matr2[i][j] = S2[i + j];
			//cout << i + j << " ";

		}
		//cout << endl;
	}

	Matr Mt1(N1, Matr1);
	Matr Mt2(N2, Matr2);
	Sist Sist1(Mt1, M1);
	Sist Sist2(Mt2, M2);
	double* Rez1 = Sist1.Gauss();
	double* Rez2 = Sist2.Gauss();
	cout << "Коэффициенкты аппроксимирующего многочлена 4 степени:" << endl;
	for (int i = 0; i < N1; i++) {
		cout << "c" << i << "  " << Rez1[i] << endl;
	}
	cout << "Наилучшее среднеквадратичное приближение для аппроксимирующего многочлена 4 степени:" << endl;
	double delta1 = 0;
	for (int i = 0; i < N + 1; i++) {
		double fi_k = 0;
		for (int j = 0; j < N1; j++) {
			fi_k += Rez1[j] * pow(X_k[i], j);
		}
		delta1 += pow((F_k[i] - fi_k), 2);
	}
	cout << delta1 << endl;
	cout << endl;


	cout << "Коэффициенкты аппроксимирующего многочлена 6 степени:" << endl;
	for (int i = 0; i < N2; i++) {
		cout << "c" << i << "  " << Rez2[i] << endl;;
	}
	cout << "Наилучшее среднеквадратичное приближение для аппроксимирующего многочлена 6 степени:" << endl;
	double delta2 = 0;
	for (int i = 0; i < N + 1; i++) {
		double fi_k = 0;
		for (int j = 0; j < N2; j++) {
			fi_k += Rez2[j] * pow(X_k[i], j);
		}
		delta2 += pow((F_k[i] - fi_k), 2);
	}
	cout << delta2 << endl;
	cout << endl;
}
