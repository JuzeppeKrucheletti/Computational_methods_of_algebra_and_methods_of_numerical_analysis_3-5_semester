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
			M[i] = new double[i + 1];
		}
		srand(time(0));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				M[i][j] = rand() % 200 - 100;
			}
		}
	}
	Matr(const Matr& B)
	{
		this->n = B.n;
		M = new double* [n];
		for (int i = 0; i < n; i++)
		{
			M[i] = new double[i + 1];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j <= i; j++)
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
			for (int j = 0; j < i; j++)
			{
				os << A.M[i][j] << " ";
			}
			for (int j = i; j < A.n; j++)
			{
				os << A.M[j][i] << " ";
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
			M[i] = new double[n - i];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n - i; j++)
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
	void StepennoiMetod()
	{
		double E = pow(10, -7);
		double* Y = new double[n];
		double* U = new double[n];
		double norm_y;
		double norm;
		double norm_pred = 0;
		double L;
		int k_max = 0;
		int k = 0;
		for (int i = 0; i < n; i++)
		{
			Y[i] = 1;
		}
		while (true)
		{
			norm_y = 0;
			for (int i = 0; i < n; i++)
			{
				norm_y += Y[i] * Y[i];
			}
			norm_y = sqrt(norm_y);
			for (int i = 0; i < n; i++)
			{
				U[i] = Y[i] / norm_y;
			}
			//


			for (int i = 0; i < n; i++)
			{
				Y[i] = 0;
				for (int j = 0; j < i; j++)
				{
					Y[i] += M[i][j] * U[j];
				}
				for (int j = i; j < n; j++)
				{
					Y[i] += M[j][i] * U[j];
				}

			}
			L = 0;
			for (int i = 0; i < n; i++)
			{
				L += Y[i] * U[i];
			}
			//
			norm = 0;
			for (int i = 0; i < n; i++)
			{
				norm += (Y[i] - L * U[i]) * (Y[i] - L * U[i]);
			}
			norm = sqrt(norm);
			if (norm < E) break;
			else
			{
				if (norm > norm_pred)
					k_max++;
				else
					k_max = 0;
			}
			if (k_max > 50)
			{
				cout << "Алгоритм расходится. Количество итераций: " << k << endl;
				break;
			}
			//
			norm_pred = norm;
			k++;
		}
		cout << "Начальный вектор Y0:" << endl;
		for (int i = 0; i < n; i++)
		{
			cout << "1" << endl;
		}
		cout << "Номер итерации k, при котором была достигнута требуемая точность:" << k << endl;
		cout << "Приближенное наибольшее по модулю собственное значение:" << L << endl;
		cout << "Соответствующий ему вектор U:" << endl;
		for (int i = 0; i < n; i++)
		{
			cout << U[i] << endl;
		}
		cout << "Вектор AU-LU:" << scientific << endl;
		for (int i = 0; i < n; i++)
		{
			cout << Y[i] - L * U[i] << endl;
		}
		cout << "Норма AU-LU:" << scientific << norm << endl;
	}
	void MetodVrashchenij()
	{
		bool t_cr = true;
		double E = pow(10, -7);
		double cos_f, sin_f;
		int max_j, max_i;
		double norm;
		int k = 0, k_max = 0;
		double norm_pred = 0;
		double** A = new double* [n];
		double** A1 = new double* [n];
		for (int i = 0; i < n; i++)
		{
			A[i] = new double[n];
			A1[i] = new double[n];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < i; j++)
			{
				A1[i][j] = A[i][j] = M[i][j];
			}
			for (int j = i; j < n; j++)
			{
				A1[i][j] = A[i][j] = M[j][i];
			}
		}
		double** T = new double* [n];
		for (int i = 0; i < n; i++)
		{
			T[i] = new double[n];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)
					T[i][j] = 1;
				else
					T[i][j] = 0;
			}
		}
		while (true)
		{
			double max = 0;
			max_j = 1, max_i = 0;
			for (int i = 0; i < n; i++)
			{
				for (int j = i + 1; j < n; j++)
				{
					if (abs(A[i][j]) > max)
					{
						max = abs(A[i][j]);
						max_i = i;
						max_j = j;
					}
				}
			}
			if (A[max_i][max_i] == A[max_j][max_j])
			{
				cos_f = 1 / sqrt(2);
				sin_f = -1 / sqrt(2);
			}
			else
			{
				double m = (2 * A[max_i][max_j]) / (A[max_i][max_i] - A[max_j][max_j]);
				cos_f = sqrt((0.5 + (0.5 / sqrt(1 + m * m))));
				sin_f = m * sqrt((0.5 - (0.5 / sqrt(1 + m * m)))) / abs(m);
			}

			double* str1, * str2;
			str1 = new double[n];
			str2 = new double[n];
			if (t_cr)
			{
				t_cr = false;
				T[max_i][max_i] = cos_f;
				T[max_j][max_j] = cos_f;
				T[max_i][max_j] = -sin_f;
				T[max_j][max_i] = sin_f;
			}
			else
			{
				double* mass1, * mass2;
				mass1 = new double[n];
				mass2 = new double[n];
				for (int i = 0; i < n; i++)
				{
					mass1[i] = (cos_f * T[i][max_i] + (sin_f)*T[i][max_j]);
					mass2[i] = (-sin_f) * T[i][max_i] + (cos_f)*T[i][max_j];
				}
				for (int i = 0; i < n; i++)
				{
					T[i][max_i] = mass1[i];
					T[i][max_j] = mass2[i];
				}
			}


			for (int i = 0; i < n; i++)
			{
				str1[i] = cos_f * A[max_i][i] + (sin_f)*A[max_j][i];
				str2[i] = (-sin_f) * A[max_i][i] + (cos_f)*A[max_j][i];
			}
			for (int i = 0; i < n; i++)
			{
				A[max_i][i] = str1[i];
				A[max_j][i] = str2[i];
			}
			for (int i = 0; i < n; i++)
			{
				str1[i] = cos_f * A[i][max_i] + (sin_f)*A[i][max_j];
				str2[i] = (-sin_f) * A[i][max_i] + (cos_f)*A[i][max_j];
			}
			for (int i = 0; i < n; i++)
			{
				A[i][max_i] = str1[i];
				A[i][max_j] = str2[i];
			}
			norm = 0;
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					if (i != j)
						norm += A[i][j] * A[i][j];
					//cout << A[i][j] << " ";
				}
				//cout << endl;
			}

			if (norm < E) break;
			else
			{
				if (norm > norm_pred)
					k_max++;
				else
					k_max = 0;
			}
			if (k_max > 50)
			{
				cout << "Алгоритм расходится. Количество итераций: " << k << endl;
				break;
			}
			//
			norm_pred = norm;
			k++;
		}


		cout << "Номер итерации k, при котором была достигнута требуемая точность:" << k << endl;
		cout << "Приближенные собственные значения:" << endl;
		for (int i = 0; i < n; i++)
		{
			cout << fixed << showpoint << setprecision(3) << A[i][i] << endl;
		}
		cout << "Соответствующиe им собственные векторa:" << endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << scientific << T[i][j] << "    ";
			}
			cout << endl;
		}
		cout << "Векторы AX-LX:" << endl;
		double** R = new double* [n];
		for (int i = 0; i < n; i++)
		{
			R[i] = new double[n];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				R[i][j] = 0;
				for (int k = 0; k < n; k++)
				{
					R[i][j] += A1[i][k] * T[k][j];
				}
				R[i][j] -= T[i][j] * A[j][j];
				cout << scientific << R[i][j] << "   ";
			}
			cout << endl;
		}
		cout << fixed << showpoint << setprecision(3);
	}
};

