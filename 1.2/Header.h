#pragma once
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;
class TD_Sist
{
private:
	int n;
	double* A;
	double* B;
	double* C;
	double* Y;
	double* F;
	double* Y_;
	double pogr;
public:
	int size()
	{
		return n;
	}
	void ShowPogr()
	{
		cout << scientific << pogr << fixed << showpoint << setprecision(3);
	}
	TD_Sist(int k)
	{
		this->n = k - 1;
		A = new double[n + 1];
		B = new double[n + 1];
		C = new double[n + 1];
		Y = new double[n + 1];
		Y_ = new double[n + 1];
		F = new double[n + 1];
		srand(time(0));
		for (int i = 1; i <= n; i++)
		{
			A[i] = (((double)rand() / RAND_MAX) * 2 - 1) * 100;
		}
		for (int i = 0; i <= n - 1; i++)
		{
			B[i] = (((double)rand() / RAND_MAX) * 2 - 1) * 100;
		}
		for (int i = 0; i <= n; i++)
		{
			Y[i] = i + 1;
			Y_[i] = 0;
			if (i == 0)
			{
				C[0] = abs(B[0]) + 5 + (((double)rand() / RAND_MAX) * 5);
			}
			else if (i == n)
			{
				C[n] = abs(A[n]) + 5 + (((double)rand() / RAND_MAX) * 5);
			}
			else
			{
				C[i] = abs(A[i]) + abs(B[i]) + 5 + (((double)rand() / RAND_MAX) * 5);
			}
		}
		F[0] = C[0] * Y[0] - B[0] * Y[1];
		F[n] = C[n] * Y[n] - A[n] * Y[n - 1];
		for (int i = 1; i < n; i++)
		{
			F[i] = C[i] * Y[i] - A[i] * Y[i - 1] - B[i] * Y[i + 1];
		}
	}

	~TD_Sist()
	{
		delete[] A;
		delete[] B;
		delete[] C;
		delete[] Y;
		delete[] Y_;
		delete[] F;
	}

	void Progonka()
	{
		//Прямая прогонка
		double* a, * b;
		a = new double[n + 1];
		b = new double[n + 2];
		a[1] = B[0] / C[0];
		b[1] = F[0] / C[0];
		for (int i = 1; i <= n - 1; i++)
		{
			double zn = C[i] - A[i] * a[i];
			//if (abs(zn) < 0,00000001) throw "Знаменатель мал";
			a[i + 1] = B[i] / zn;
			b[i + 1] = (F[i] + A[i] * b[i]) / zn;
		}
		b[n + 1] = (F[n] + A[n] * b[n]) / (C[n] - A[n] * a[n]);
		//Обратная прогонка
		Y_[n] = b[n + 1];
		for (int i = n - 1; i >= 0; i--)
		{
			Y_[i] = a[i + 1] * Y_[i + 1] + b[i + 1];
		}


		//pogr
		double maxY = abs(Y[0]), max_Y = abs(Y[0] - Y_[0]);
		for (int i = 0; i <= n; i++)
		{
			if (max_Y < abs(Y[i] - Y_[i]))
			{
				max_Y = abs(Y[i] - Y_[i]);
			}
			if (maxY < abs(Y[i]))
			{
				maxY = abs(Y[i]);
			}
		}
		pogr = max_Y / maxY;
	}
	friend ostream& operator <<(ostream& os, TD_Sist& M)
	{
		os << "Вектор А:" << endl << fixed << showpoint << setprecision(5);
		for (int i = 1; i <= M.size(); i++)
		{
			os << M.A[i] << " ";
		}
		os << endl;
		os << "Вектор B:" << endl;
		for (int i = 0; i <= M.size() - 1; i++)
		{
			os << M.B[i] << " ";
		}
		os << endl;
		os << "Вектор C:" << endl;
		for (int i = 0; i <= M.size(); i++)
		{
			os << M.C[i] << " ";
		}
		os << endl;
		os << "Точное решение Y:" << endl << dec;
		for (int i = 0; i <= M.size(); i++)
		{
			os << M.Y[i] << " ";
		}
		os << endl;
		os << "Вектор F:" << endl << fixed << showpoint << setprecision(5);
		for (int i = 0; i <= M.size(); i++)
		{
			os << M.F[i] << " ";
		}
		os << endl;
		os << "Полученное решение Y_:" << endl << fixed << showpoint << setprecision(20);
		for (int i = 0; i <= M.size(); i++)
		{
			os << M.Y_[i] << " " << endl;
		}
		os << "Относительная погрешность: " << scientific << M.pogr << endl << dec;
		return os;
	}
};
