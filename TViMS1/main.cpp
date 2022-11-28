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
	MC.MC_Print();
	MC.P_Test(10000, 0.4, 0.3);
	system("pause");
	return 0;
}