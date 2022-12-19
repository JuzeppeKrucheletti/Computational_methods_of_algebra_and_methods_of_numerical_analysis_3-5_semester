#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <iostream>
#include "Header.h"
#include <fstream>
#include<time.h>
using namespace std;
int main()
{
	setlocale(LC_ALL, "ru");
	srand(time(NULL));
	MarkovChain MC;
	double alpha[] = { 0.1, 0.05, 0.01 }, beta[] = { 0.1, 0.05, 0.01 };
	MC.MC_Print();
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			cout << "ѕоследовательный тест при ошибках первого и второго рода соответственно: " << alpha[i] << "  и  " << beta[j] << endl;
			MC.P_Test(10000, alpha[i], beta[j]);
		}
	}
	system("pause");
	return 0;
}