#pragma once
#include <iostream>
#include <cstdio>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>
#include <cmath>
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

		M[0][0] = rand() % 200 - 100;

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
			for (int j = 0; j <= i; j++)
			{
				if (i != j)
					M[i][j] = M[j][i] = rand() % 200 - 100;
				else
					M[i][j] = 0;
			}
		}
		for (int i = 0; i < n; i++)
		{
			int summ = 0;
			for (int j = 0; j < n; j++)
			{
				summ += abs(M[i][j]);
			}
			M[i][i] = summ + 5 + rand() % 45;
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
};
class Sist
{
	friend class Matr;
private:
	Matr A;
	double* X;
	double* F;
	double* X_;
	const double E = 1 * pow(10, -7);
	const  int k_max = 5000;
	double norm = 0;
	int q = 0;
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
	void ShowNorm()
	{
		cout << scientific << norm << fixed << showpoint << setprecision(3);
	}
	void ShowQ()
	{
		cout << dec << q << fixed << showpoint << setprecision(3);
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
			X_[i] = F[i];
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
		q = B.q;
	}
	~Sist()
	{
		delete[] X;
		delete[] X_;
		delete[] F;
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
		os << "Начальное приближение:" << endl;
		for (int i = 0; i < B.size(); i++)
		{
			os << B.F[i] << " ";
			os << endl;
		}
		return os;
	}
	void Show()
	{
		cout << "Номер итерации q, при которой достигнута требуемая точность:" << q << endl;

		cout << "Приближенное решение X_:" << endl;
		for (int i = 0; i < size(); i++)
		{
			cout << fixed << showpoint << setprecision(20) << X_[i] << " ";
			cout << endl;
		}
		cout << "Норма AX-F: " << scientific << norm << endl << fixed << showpoint << setprecision(3);
		cout << "Абсолютная погрешность: " << scientific << pogr << endl << fixed << showpoint << setprecision(3);



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
	void MMN()
	{
		for (int i = 0; i < size(); i++)
		{
			X_[i] = F[i];
		}
		double* V = new double[size()];
		double* Ax = new double[size()];
		double* Av = new double[size()];
		double T;
		norm = E;
		for (q = 0; q < k_max; q++)
		{
			//1
			norm = 0;
			for (int i = 0; i < size(); i++)
			{
				Ax[i] = 0;
				for (int j = 0; j < size(); j++)
				{
					Ax[i] += A.M[i][j] * X_[j];
				}
				V[i] = Ax[i] - F[i];
				norm += abs(V[i]);
			}
			if (norm < E) break;
			//2
			for (int i = 0; i < size(); i++)
			{
				Av[i] = 0;
				for (int j = 0; j < size(); j++)
				{
					Av[i] += A.M[i][j] * V[j];
				}
			}
			double Av_v = 0, Av_Av = 0;
			for (int i = 0; i < size(); i++)
			{
				Av_v += Av[i] * V[i];
				Av_Av += Av[i] * Av[i];
			}
			T = Av_v / Av_Av;
			//3
			for (int i = 0; i < size(); i++)
			{
				X_[i] = X_[i] - T * V[i];
			}
		}
		if (q == k_max)
			cout << "Итерационный процесс остановлен при достижении максимального числа итераций" << endl;
		delete[] Av;
		delete[] Ax;
		delete[] V;
		pogr = 0;
		for (int i = 0; i < size(); i++)
		{
			pogr += abs(X[i] - X_[i]);
		}
		this->Show();
	}

	void Relax(double w)
	{
		for (int i = 0; i < size(); i++)
		{
			X_[i] = F[i];
		}
		if (w < 0 || w > 2)
			cout << "w должно принадлежать промежутку (0;2)" << endl;
		double* Ax = new double[size()];
		for (q = 0; q < k_max; q++)
		{
			for (int i = 0; i < size(); i++)
			{
				double summ1 = 0, summ2 = 0;
				for (int j = 0; j < i; j++)
				{
					summ1 += A.M[i][j] * X_[j];
				}
				for (int j = i + 1; j < size(); j++)
				{
					summ2 += A.M[i][j] * X_[j];
				}
				X_[i] = (1 - w) * X_[i] + (w / A.M[i][i]) * (F[i] - (summ1 + summ2));
			}
			norm = 0;
			for (int i = 0; i < size(); i++)
			{
				Ax[i] = 0;
				for (int j = 0; j < size(); j++)
				{
					Ax[i] += A.M[i][j] * X_[j];
				}
				norm += abs(Ax[i] - F[i]);
			}
			if (norm < E) break;
		}
		if (q == k_max)
			cout << "Итерационный процесс остановлен при достижении максимального числа итераций" << endl;
		pogr = 0;
		for (int i = 0; i < size(); i++)
		{
			pogr += abs(X[i] - X_[i]);
		}
		delete[] Ax;
	}
};

