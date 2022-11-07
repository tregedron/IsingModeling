#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <list>
#include <time.h>
#include <algorithm>

using namespace std;

const double kT = 10;

//вывод матриц
void print_matrix(double** matrix, int N)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout.width(5);
			cout.precision(3);
			cout << matrix[i][j]<<" ";
		}
		cout << endl;
	}
	cout << endl;
}

//функция, считающая матрицу взаимодействия и дистанций раз и навсегда.
void calc_Jmatrix(double** Jmatrix, int N, double a, double b, double** DistMatrix)
{
	double r = 0;
	double full_lap = (a + b) * N / 2;
	int chered = 1;
	for (int i = 0; i < N; i++)
	{
		double temp_for = 0;
		double temp_rev = full_lap;
		for (int j = 0; j < N; j++)
		{
			if (j >= i + 1)
			{

				if (chered % 2 == 0)
				{
					temp_for += a;
					temp_rev -= a;
					DistMatrix[i][j] = min(temp_for, temp_rev); 
					DistMatrix[j][i] = DistMatrix[i][j];
				}
				else
				{
					temp_for += b;
					temp_rev -= b;
					DistMatrix[i][j] = min(temp_for, temp_rev);
					DistMatrix[j][i] = DistMatrix[i][j];
				}
				Jmatrix[i][j] = 1.0 / (DistMatrix[i][j] * DistMatrix[i][j]);
				Jmatrix[j][i] = Jmatrix[i][j];
			}
			chered++;
		}
	}

}



//Функция, расчитывающая энергию системы "все со всеми, J=1". Параметры: указатель на массив спинов, его размер, указатель на переменную для хранения изменения энергии.
void calc_energy(double* spins_array, int N, double* Energy)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{

			*Energy += spins_array[i] * spins_array[j];
		}
	}
	cout << "Start Energy: " << *Energy << endl;
}

//Функция, расчитывающая энергию системы с матрицей взаимодействия. Параметры: указатель на массив спинов, его размер, указатель на переменную для хранения изменения энергии, матрица взаимодействия.
void calc_energy_Jmatrix(double* spins_array, int N, double* Energy, double** Jmatrix)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{

			*Energy += Jmatrix[i][j] * spins_array[i] * spins_array[j];
		}
	}
	cout << "Start Energy: " << *Energy << endl;
}

//Функция, меняющая случайно один из спинов, расчитывающая изменение энергии при таком изменении и принимающая шаг. Параметры: указатель на массив спинов, его размер, указатель на переменную для хранения изменения энергии,
void change_spin(double* spins_array, int N, double* Energy, double* dEnergy, int* success)
{
	//выбрали спин, перещёлкнули его
	int i = rand() % N;
	spins_array[i] = -1 * spins_array[i];

	//обнулили изменение энергии с прошлого шага и посчитали её
	*dEnergy = 0;
	for (int j = 0; j < N; j++)
	{
		if (i == j)
		{
			continue;
		}
		else
		{
			*dEnergy += spins_array[i] * spins_array[j];
		}
	}

	//принимаем решение о принятии/отклонении шага. Шаг принимаем пропорционально e^(-dE).
	if (1.0 * rand() / RAND_MAX < exp(-2.0 * (*dEnergy) / (kT)))
	{
		*Energy += 2.0 * (*dEnergy);
		*success += 1;

	}
	else
	{
		spins_array[i] = -1 * spins_array[i];
	}
}

//Функция, меняющая случайно один из спинов, расчитывающая изменение энергии при таком изменении и принимающая шаг. Параметры: указатель на массив спинов, его размер, указатель на переменную для хранения изменения энергии
//указатель на переменную подсчёта числа успешных принятий шагов и матрицу взаимодействий.
void change_spin_Jmatrix(double* spins_array, int N, double* Energy, double* dEnergy, int* success, double** Jmatrix)
{
	//выбрали спин, перещёлкнули его
	int i = rand() % N;
	spins_array[i] = -1 * spins_array[i];

	//обнулили изменение энергии с прошлого шага и посчитали её
	*dEnergy = 0;
	for (int j = 0; j < N; j++)
	{
			*dEnergy += Jmatrix[i][j] * spins_array[i] * spins_array[j];
			//cout <<"!!!"<< Jmatrix[i][j] <<" "<< spins_array[i] * spins_array[j] <<" "<<"!!!";
	}
	//принимаем решение о принятии/отклонении шага. Шаг принимаем пропорционально e^(-dE).
	double p = 1.0 * rand() / RAND_MAX;
	if (p < exp(-2.0 * (*dEnergy) / (kT)))
	{
		*Energy += 2.0 * (*dEnergy);
		*success += 1;
		//cout <<"accepted "<< *dEnergy << " "<< p<<" "<< exp(-2.0 * (*dEnergy) / (kT)) <<endl;
	}
	else
	{
		spins_array[i] = -1 * spins_array[i];
		//cout <<"rejected "<< *dEnergy << " " << p << " " << exp(-2.0 * (*dEnergy) / (kT)) << endl;
	}
}

//функция вывода состояния системы на данном шаге моделирования.
void print_state(double* spins_array, int N, double* Energy, double *dEnergy, int Step)
{
	double minus = 0; double plus = 0;
	cout << endl;
	cout << Step << " Energy:" << *Energy<<" dEnergy: "<<*dEnergy << endl;
	cout << "System: ";
	for (int i = 0; i < N; i++)
	{
		cout << spins_array[i] << " ";
		if (spins_array[i] == 1)
		{
			plus++;
		}
		else
		{
			minus++;
		}
	}
	cout << endl;
	cout << "minus: " << minus << " plus: " << plus << " delta:" << plus - minus << endl;
	cout << endl;
}

//так писать очень нехорошо, нужно быть аккуратнее с округлениями. 
//Идея кода: взять массив спинов, пройти по всем парам спинов, каждый раз считая их произведение. Результат произведения класть в соответствующую ячейку массива с функцией распределения.
//Соответствующая ячейка берётся из таблицы расстояний. Итог: распределение произведения спинов от расстояния.
void calc_destribution(double* spins_array, int N, double* Spin_destribution, double** DistMatrix)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = i; j < N; j++)
		{
			Spin_destribution[(int)(DistMatrix[i][j])] += spins_array[i] * spins_array[j];
		}
	}
}

int main()
{
	srand(time(NULL));
	double Energy = 0;
	double dEnergy = 0;
	int N = 0; int Avg = 0;
	cout << "Number of spins:" << endl;
	cin >> N;
	if (N % 2 == 1)
	{
		N++; cout << "Your N wasn't even, I've changed it :)" << endl;
	}
	double* spins_array = new double[N];
	for (int i = 0; i < N; i++)
	{
		spins_array[i] = 1.0;
	}
	//Здесть заводим матрицу взаимодействия. Организуем J_{i,j} ~ 1/r^2, r = aphpa*a + beta*b, это позволит задать спины на сетке с шагом a,b.
	double a, b = 0;
	a = 1; b = 3;
	double** Jmatrix = new double* [N] {};
	for (int i = 0; i < N; i++)
	{
		Jmatrix[i] = new double[N];
		for (int j = 0; j < N; j++)
		{
			Jmatrix[i][j] = 0;
		}
	}
	double** DistMatrix = new double* [N] {};
	for (int i = 0; i < N; i++)
	{
		DistMatrix[i] = new double[N];
		for (int j = 0; j < N; j++)
		{
			DistMatrix[i][j] = 0;
		}
	}
	calc_Jmatrix(Jmatrix, N, a, b, DistMatrix);
	print_matrix(Jmatrix, N); 
	print_matrix(DistMatrix, N);

	//calc_energy_Jmatrix(spins_array, N, &Energy);
	calc_energy_Jmatrix(spins_array, N, &Energy, Jmatrix);

	int Steps = 0; int success = 0;
	cout << "Number of spin changes:" << endl;
	cin >> Steps;

	double* spin_destribution = new double[(int)((a + b) * N / 2)] {0};
	
	cout << "START";
	print_state(spins_array, N, &Energy, &dEnergy, 0);

	//выход на равновесие.
	for (int i = 0; i < 1000; i++)
	{
		change_spin_Jmatrix(spins_array, N, &Energy, &dEnergy, &success, Jmatrix);
	}
	cout << "AVERAGING...";
	print_state(spins_array, N, &Energy, &dEnergy, 0);

	for (int i = 0; i < Steps; ++i)
	{
		//change_spin(spins_array, N, &Energy, &dEnergy, &success);
		change_spin_Jmatrix(spins_array, N, &Energy, &dEnergy, &success, Jmatrix);
		//вывод на экран (впоследствии в файл) результатов.
	
		if (i % 25 == 0)
		{
			calc_destribution(spins_array, N, spin_destribution, DistMatrix); Avg++;
		}
		if (i % 10000 == 0)
		{
			print_state(spins_array, N, &Energy, &dEnergy, i);
		}
	}

	ofstream file;
	file.open("destribution.txt");
	for (int i = 0; i < (int)((a + b) * N / 2); i++)
	{
		file << i << "   " << (double)( spin_destribution[i] / Avg) << endl;
	}


	cout << "rate: " << success << endl;
	delete [] Jmatrix;
	delete [] spins_array;
	delete [] spin_destribution;
	return 0;
}

