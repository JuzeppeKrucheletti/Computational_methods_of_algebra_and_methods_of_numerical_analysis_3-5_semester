#pragma once
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;
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
	Matr(int n)
	{
		this->n = n;
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
	Matr operator * (const Matr& B)
	{
		if (n != B.n) throw "Матрицы имеют разные размеры";
		Matr P(n);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				P.M[i][j] = 0;
				for (int k = 0; k < n; k++)
				{
					P.M[i][j] += M[i][k] * B.M[k][j];
				}
			}
		}
		return P;
	}

	Matr Obrat()
	{
		Matr E(n), Obr(n), O(*this);
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j) E.M[i][j] = 1;
				else E.M[i][j] = 0;
			}
		}
		for (int i = 0; i < n; i++)
		{
			//выбор главного элемента по столбцу
			double max = abs(O.M[i][i]);
			int gl = i;
			for (int j = i; j < n; j++)
			{
				if (abs(O.M[j][i]) > max)
				{
					max = O.M[j][i];
					gl = j;
				}
			}
			if (max == 0) throw "Все главные элементы нулевые";
			//исключение в случае, когда все главные элементы нулевые
			for (int j = i; j < n; j++)
			{
				swap(O.M[i][j], O.M[gl][j]);
			}
			for (int j = 0; j < n; j++)
			{
				swap(E.M[i][j], E.M[gl][j]);
			}
			//прямой ход метода Гаусса
			for (int j = i + 1; j < n; j++)
			{
				O.M[i][j] /= O.M[i][i];
			}
			for (int j = 0; j < n; j++)
			{
				E.M[i][j] /= O.M[i][i];
			}
			O.M[i][i] /= O.M[i][i];
			for (int j = i + 1; j < n; j++)
			{
				for (int k = i + 1; k < n; k++)
				{
					O.M[j][k] -= O.M[i][k] * O.M[j][i];
				}
				for (int k = 0; k < n; k++)
				{
					E.M[j][k] -= E.M[i][k] * O.M[j][i];
				}
				O.M[j][i] -= O.M[i][i] * O.M[j][i];
			}
		}
		//обратный ход метода Гаусса
		for (int k = 0; k < n; k++)
		{
			for (int i = (n - 1); i >= 0; i--)
			{
				Obr.M[i][k] = E.M[i][k];
				for (int j = n - 1; j > i; j--)
				{
					Obr.M[i][k] -= O.M[i][j] * Obr.M[j][k];
				}
			}
		}
		return Obr;
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
	Sist(const Matr& B)
	{
		A = B;
		X = new double[B.n];
		for (int i = 0; i < B.n; i++)
		{
			X[i] = 1 + i;
		}
		F = new double[B.n];
		for (int i = 0; i < size(); i++)
		{
			F[i] = 0;
			for (int j = 0; j < size(); j++)
			{
				F[i] += A.M[i][j] * X[j];
			}
		}
		X_ = new double[B.n];
		for (int i = 0; i < size(); i++)
		{
			X_[i] = 0;
		}
	}
	Sist(const Sist& B)
	{
		A = B.A;
		X = new double[size()];
		F = new double[size()];
		X_ = new double[size()];
		for (int i = 0; i < size(); i++)
		{
			for (int j = 0; j < size(); j++)
			{
				A.M[i][j] = B.A.M[i][j];
			}
		}
		for (int i = 0; i < size(); i++)
		{
			X[i] = B.X[i];
		}
		for (int i = 0; i < size(); i++)
		{
			F[i] = B.F[i];
		}
		for (int i = 0; i < size(); i++)
		{
			X_[i] = B.X_[i];
		}
	}
	~Sist()
	{
		delete[] X;
		delete[] X_;
		delete[] F;
	}
	void Gauss()
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
	}
	friend ostream& operator <<(ostream& os, Sist& B)
	{
		os << "Матрица системы:" << endl << fixed << showpoint << setprecision(3);
		os << B.A << endl;
		os << "Столбец свободных членов:" << endl;
		for (int i = 0; i < B.size(); i++)
		{
			os << B.F[i] << " ";
			os << endl;
		}
		os << "Точное решение X:" << endl;
		for (int i = 0; i < B.size(); i++)
		{
			os << B.X[i] << " ";
			os << endl;
		}
		os << "Приближенное решение X_:" << endl;
		for (int i = 0; i < B.size(); i++)
		{
			os << fixed << showpoint << setprecision(20) << B.X_[i] << " ";
			os << endl;
		}
		os << "Относительная погрешность: " << fixed << showpoint << setprecision(20) << B.pogr << endl << fixed << showpoint << setprecision(3);
		return os;
	}

	Sist& operator = (const Sist& B)
	{
		if (this == &B)
			return *this;
		A = B.A;
		delete[] X;
		delete[] X_;
		delete[] F;
		X = new double[A.n];
		F = new double[A.n];
		X_ = new double[A.n];
		for (int i = 0; i < A.n; i++)
		{
			X[i] = B.X[i];
		}
		for (int i = 0; i < A.n; i++)
		{
			F[i] = B.F[i];
		}
		for (int i = 0; i < A.n; i++)
		{
			X_[i] = B.X_[i];
		}

		return *this;
	}

};
