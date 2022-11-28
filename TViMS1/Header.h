#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <list>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <time.h>
using namespace std;
class MarkovChain {
private:
	int M;
	bool H0;
	double* p0, * p1;
	double** P0, ** P1;
public:
	MarkovChain() {
		/*cout << "Введите расмерность дискретного пространства X:" << endl;
		cin >> M;*/
		M = 2;
		p0 = new double[M];
		p1 = new double[M];
		P0 = new double* [M];
		P1 = new double* [M];
		for (int i = 0; i < M; i++) {
			P0[i] = new double[M];
			P1[i] = new double[M];
		}
		p0[0] = 0.6;
		p0[1] = 0.4;
		p1[0] = 0.3;
		p1[1] = 0.7;
		P0[0][0] = 0.6;
		P0[0][1] = 0.4;
		P0[1][0] = 0.7;
		P0[1][1] = 0.3;

		P1[0][0] = 0.2;
		P1[0][1] = 0.8;
		P1[1][0] = 0.4;
		P1[1][1] = 0.6;
		

		/*p0 = new double[M];
		p1 = new double[M];
		P0 = new double* [M];
		P1 = new double* [M];
		for (int i = 0; i < M; i++) {
			P0[i] = new double[M];
			P1[i] = new double[M];
		}
	     
		double summ = 0;
		while (summ != 1) {
			cout << "Задайте начальные вероятности p0:" << endl;
			for (int i = 0; i < M; i++) {
				cin >> p0[i];
				summ += p0[i];
				if (summ > 1) {
					cout << "Сумма начальных вероятностей дожна равняться 1....." << endl;
					summ = 0;
					break;
				}
				if (p0[i] < 0) {
					cout << "Вероятности не могут быть отрицательными." << endl;
					summ = 0;
					break;
				}
			}
		}

		summ = 0;
		while (summ != 1) {
			cout << "Задайте начальные вероятности p1:" << endl;
			for (int i = 0; i < M; i++) {
				cin >> p1[i];
				summ += p1[i];
				if (summ > 1) {
					cout << "Сумма начальных вероятностей дожна равняться 1." << endl;
					summ = 0;
					break;
				}
				if (p1[i] < 0) {
					cout << "Вероятности не могут быть отрицательными." << endl;
					summ = 0;
					break;
				}
			}
		}

		cout << "Задайте матрицу состояний P0:" << endl;
		for (int i = 0; i < M; i++) {
			summ = 0;
			while (summ != 1) {
				cout << "Вероятности перехода из состояния " << i <<":" << endl;
				for (int j = 0; j < M; j++) {
					cin >> P0[i][j];
					summ += P0[i][j];
					if (summ > 1) {
						cout << "Сумма начальных вероятностей дожна равняться 1." << endl;
						summ = 0;
						break;
					}
					if (P0[i][j] < 0) {
						cout << "Вероятности не могут быть отрицательными." << endl;
						summ = 0;
						break;
					}
				}
			}
		}

		cout << "Задайте матрицу состояний P1:" << endl;
		for (int i = 0; i < M; i++) {
			summ = 0;
			while (summ != 1) {
				cout << "Вероятности перехода из состояния " << i << ":" << endl;
				for (int j = 0; j < M; j++) {
					cin >> P1[i][j];
					summ += P1[i][j];
					if (summ > 1) {
						cout << "Сумма начальных вероятностей дожна равняться 1." << endl;
						summ = 0;
						break;
					}
					if (P1[i][j] < 0) {
						cout << "Вероятности не могут быть отрицательными." << endl;
						summ = 0;
						break;
					}
				}
			}
		}*/

	};

	~MarkovChain() {
		delete[] p0, p1;
		for (int i = 0; i < M; i++) {
			delete[] P0[i];
			delete[] P1[i];
		}
		delete[] P0, P1;
	}
	void MC_Print() {
		cout << "Начальные вероятности p0:" << endl;
		for (int i = 0; i < M; i++) {
			cout << " " << p0[i] << endl;
		}
		cout << "Начальные вероятности p1:" << endl;
		for (int i = 0; i < M; i++) {
			cout << " " << p1[i] << endl;
		}
		cout << "Матрица состояний P0:" << endl;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				cout << P0[i][j] << " ";
			}
			cout << endl;
		}
		cout << "Матрица состояний P1:" << endl;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++) {
				cout << P1[i][j] << " ";
			}
			cout << endl;
		}
	}
	void P_Test(int N, double alpha, double beta) {
		cout << "/////////\n////////\n" << "Последовательная проверка при alpha = " << alpha << " и beta = " << beta << endl<<endl;
		double E_n_0 = 0, E_n_1 = 0;
		double re_alpha = 0, re_beta = 0;
		double A = (1 - beta) / alpha, B = alpha/(1-beta);
		for (int nabl = 0; nabl < N; nabl++) {
			double res;
			double p0_m, p1_m;
			int X_i_ = 0, X_i = 0;
			double random = double(rand()) / RAND_MAX;
			double summ = 0;
			for (int i = 0; i < M; i++) {
				summ += p0[i];
				if (random < summ) {
					X_i = i;
					break;
				}
			}
			p0_m = p0[X_i];
			p1_m = p1[X_i];
			int j = 1;
			while (true) {
				//cout << X_i << "  ";
				if (p1_m / p0_m >= A) {
					res = 1;
					break;
				}
				else if (p1_m / p0_m <= B) {
					res = 0;
					break;
				}
				else {
					random = double(rand()) / RAND_MAX;
					summ = 0;
					for (int i = 0; i < M; i++) {
						summ += P0[X_i][i];
						if (random < summ) {
							X_i_ = i;
							break;
						}
					}
					p0_m *= P0[X_i][X_i_];
					p1_m *= P1[X_i][X_i_];
					j++;
					X_i = X_i_;
				}
			}
			//cout << endl;
			re_alpha += double(res);
			E_n_0 += j;
		}
		re_alpha /= N;
		E_n_0 /= N;
		cout << "Реальная вероятность ошибки первого рода: " << re_alpha << endl;
		cout << "Среднее число наблюдений при истинной гипотезе H0: " << E_n_0 << endl;


		for (int nabl = 0; nabl < N; nabl++) {
			double res;
			double p0_m, p1_m;
			int X_i_ = 0, X_i = 0;
			double random = double(rand()) / RAND_MAX;
			double summ = 0;
			for (int i = 0; i < M; i++) {
				summ += p1[i];
				if (random < summ) {
					X_i = i;
					break;
				}
			}
			p0_m = p0[X_i];
			p1_m = p1[X_i];
			int j = 1;
			while (true) {
				//cout << X_i << "  ";
				if (p1_m / p0_m >= A) {
					res = 0;
					break;
				}
				else if (p1_m / p0_m <= B) {
					res = 1;
					break;
				}
				else {
					random = double(rand()) / RAND_MAX;
					summ = 0;
					for (int i = 0; i < M; i++) {
						summ += P1[X_i][i];
						if (random < summ) {
							X_i_ = i;
							break;
						}
					}
					p0_m *= P0[X_i][X_i_];
					p1_m *= P1[X_i][X_i_];
					j++;
					X_i = X_i_;
				}
			}
			//cout << endl;
			re_beta += double(res);
			E_n_1 += j;
		}
		re_beta /= N;
		E_n_1 /= N;
		cout <<"Реальная вероятность ошибки второго рода: " << re_beta << endl;
		cout << "Среднее число наблюдений при истинной гипотезе H1: " << E_n_1 << endl << "/////////\n////////\n";
	}

	


};
