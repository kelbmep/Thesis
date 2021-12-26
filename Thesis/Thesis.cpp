#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <map>
#include <boost/multiprecision/cpp_dec_float.hpp>

using boost::multiprecision::cpp_dec_float_50;

int main()
{
	std::ifstream fin("src/speedtest-500MB.bin", std::ios::binary);

	cpp_dec_float_50 SH[40], TS[40], RE[40];
	fin.seekg(0, std::ios::end);

	//Длина последовательности
	int len = (int)fin.tellg();
	len /= 8;
	int tem = len;

	unsigned int start, end, global, globalend, time;
	global = clock();
	fin.seekg(0, std::ios::beg);
	std::cout << len << std::endl;

	//Параметр, необходимый для вычисления оценок энтропии Тсаллиса и Реньи
	const int r = 2;

	//Фиксированное значение такое, чтобы массив помещался в ОЗУ
	const int s0 = 26;

	//Количество s-грамм фиксированной длины s0
	const int N1 = (int)pow(2.0, s0);

	//Значение длины s-граммы
	int s = 42;

	//Значение разницы необходимого количества слагаемых для достижения заданной точности и параметра распределения Пуассона
	int dif = 0;

	while (s >= 3)
	{
		//Отрезание лишних битов на конце так, чтобы длина последовательности была кратна длине s-граммы
		len = tem - tem % s;
		start = clock();

		std::cout << "s = " << s << std::endl;
		//Количество различных всевозможных s-грамм заданной длины
		long long N = (long long)pow(2.0, s);

		//Общее количество всех s-грамм заданной длины для последовательности
		int n = len / s * 8;

		//Значение параметра распределения Пуассона
		cpp_dec_float_50 l = 1.0 * n / N;

		//Максимальное число итераций для вычисления значений математического ожидания и дисперсии оценки энтропии Шеннона
		cpp_dec_float_50 K;

		if (2 * l < 500)
			K = 500;
		else
			K = l * 2;
		//Вспомогательные переменные для вычисления значений математического ожидания и дисперсии оценки энтропии Шеннона
		cpp_dec_float_50 S = 0, b = 1, a, A = 0, B = 0;

		//Значение математического ожидания на предыдущем шаге для оценки энтропии Шеннона
		cpp_dec_float_50 mew1 = log(n);

		//Значение дисперсии на предыдущем шаге для оценки энтропии Шеннона
		cpp_dec_float_50 sig1 = 0;

		//Цикл для вычисления математического ожидания и дисперсии оценки энтропии Шеннона
		for (int k = 1; k <= K; k++)
		{
			b *= l / k;
			if (k > l - dif * (s + 1) * 2 / s || dif == 0)
			{
				a = b * log(k + 1);
				S += a;
				A += a * (k + 1) * log(k + 1);
				B += a * (k - l + 1);
				if (k > l)
				{
					//Значение математического ожидания на данном шаге для оценки энтропии Шеннона
					cpp_dec_float_50 mew = log(n) - exp(-l) * S;

					//Значение дисперсии на данном шаге для оценки энтропии Шеннона
					cpp_dec_float_50 sig = exp(-l) / n * (A - exp(-l) * (pow(S, 2) * l + pow(B, 2)));

					//Проверка достижения заданной точности значения дисперсии для прерывания цикла
					if (fabs(sig - sig1) < sig1 * pow(10, -17))
					{
						dif = k - (int)l;
						break;
					}
					mew1 = mew;
					sig1 = sig;
				}
			}
		}
		//Значение математического ожидания для оценки энтропии Шеннона
		cpp_dec_float_50 mew = log(n) - exp(-l) * S;

		//Значение дисперсии для оценки энтропии Шеннона
		cpp_dec_float_50 sig = exp(-l) / n * (A - exp(-l) * (pow(S, 2) * l + pow(B, 2)));

		std::cout << std::setprecision(15) << mew << " " << sig << std::endl;

		//Строка, в которой записана последовательность
		char* buf = new char[len];
		fin.seekg(0, std::ios::beg);
		fin.read(buf, len);

		//Массив частот s-грамм
		int* counter;

		//Вычисление размерности массива частот, в зависимости от длины s-граммы
		if (s > s0)
			counter = new int[N1];
		else
			counter = new int[N];

		//Вспомогательные переменные для вычиления частот s-грамм 
		int t = 0;
		long long cur = 0, temp;

		//Частная сумма для оценок энтропии Тсаллиса и Реньи 
		long long Z1 = 0;

		//Частотные оценки распределения вероятности
		double* p;

		//Значение оценки энтропии Шеннона
		double h = 0;

		//Частная сумма для оценок энтропии Тсаллиса и Реньи 
		cpp_dec_float_50 Z = 0;

		//Использование алгоритма, если массив частот помещается в ОЗУ
		if (s <= s0)
		{
			for (int i = 0; i < N; i++)
				counter[i] = 0;
			for (int i = 0; i < len; i++)
				for (int j = 7; j >= 0; j--)
				{
					cur <<= 1;
					cur += (buf[i] >> j) & 1;
					t++;
					if (t == s)
					{
						counter[cur]++;
						cur = 0;
						t = 0;
					}
				}
			p = new double[N];
			if (s < 10)
			{
				for (int i = 0; i < N; i++)
				{
					p[i] = counter[i] * 1.0 / n;
					Z += pow(counter[i], 2) - counter[i];
				}
			}
			else
			{
				for (int i = 0; i < N; i++)
				{
					p[i] = counter[i] * 1.0 / n;
					Z1 += (long long)pow(counter[i], 2) - counter[i];
				}
				Z = Z1;
			}
			for (int i = 0; i < N; i++)
				if (p[i] != 0)
					h -= p[i] * log(p[i]);

			delete[] p;
			delete[] counter;
		}
		//Использование алгоритма, если массив частот не помещается в ОЗУ
		else
		{
			//Алгоритм, без использования контейнера std::map
			if (s <= 30)
			{
				for (int i = 0; i < pow(2, s - s0); i++)
				{
					for (int g = 0; g < N1; g++)
						counter[g] = 0;
					for (int j = 0; j < len; j++)
						for (int k = 7; k >= 0; k--)
						{
							cur <<= 1;
							cur += (buf[j] >> k) & 1;
							t++;
							if (t == s - s0)
							{
								temp = cur;
								cur = 0;
							}
							else if (t == s)
							{
								if (temp == i)
									counter[cur]++;
								cur = 0;
								t = 0;
							}
						}
					p = new double[N1];
					for (int m = 0; m < N1; m++)
						p[m] = counter[m] * 1.0 / n;

					for (int m = 0; m < N1; m++)
						if (p[m] != 0)
							h -= p[m] * log(p[m]);

					for (int m = 0; m < N1; m++)
						Z1 += (long long)pow(counter[m], 2) - counter[m];

					delete[] p;
				}
				Z = Z1;
			}
			//Алгоритм, с использованием контейнера std::map
			else
			{
				std::map<long long, int> countermap;
				for (int i = 0; i < len; i++)
					for (int j = 7; j >= 0; j--)
					{
						cur <<= 1;
						cur += (buf[i] >> j) & 1;
						t++;
						if (t == s)
						{
							countermap[cur]++;
							cur = 0;
							t = 0;
						}
					}
				for (auto it = countermap.begin(); it != countermap.end(); it++)
				{
					h -= it->second * 1.0 / n * log(it->second * 1.0 / n);
					Z1 += (long long)pow(it->second, 2) - it->second;
				}
				Z = Z1;
				countermap.clear();
			}
			delete[] counter;
		}
		//Значение оценок энтропии Тсаллиса и Реньи
		cpp_dec_float_50 Hr, Sr;

		Sr = 1 / (r - 1) * (1 - Z / pow(n, r));
		Hr = log(n) + 1 / (r - 1) * (log(n) - log(Z));

		SH[s - 3] = (h - mew) / (sqrt(sig) * 1.96);
		TS[s - 3] = (Sr - (1 - 1.0 / N)) / ((sqrt(2.0 / N) / n) * 1.96);
		RE[s - 3] = (Hr - log(N)) / (sqrt(2 / (n * l)) * 1.96);

		std::cout << std::setprecision(15) << "Shennon - " << (h - mew) / (sqrt(sig) * 1.96) << std::endl;
		std::cout << "Tsallis - " << (Sr - (1 - 1.0 / N)) / ((sqrt(2.0 / N) / n) * 1.96) << std::endl;
		std::cout << "Renia - " << (Hr - log(N)) / (sqrt(2 / (n * l)) * 1.96) << std::endl;

		end = clock();
		time = end - start;
		std::cout << "Time - " << (double)time / 1000 << "s" << std::endl;

		delete[] buf;

		s--;
	}
	std::ofstream fout11("res/Shennon.txt");
	std::ofstream fout12("res/Tsallis.txt");
	std::ofstream fout13("res/Renia.txt");

	for (int i = 0; i <= 39; i++)
		fout11 << std::setprecision(15) << i + 3 << " " << SH[i] << std::endl;

	for (int i = 0; i <= 39; i++)
		fout12 << std::setprecision(15) << i + 3 << " " << TS[i] << std::endl;

	for (int i = 0; i <= 39; i++)
		fout13 << std::setprecision(15) << i + 3 << " " << RE[i] << std::endl;

	globalend = clock();
	time = globalend - global;
	std::cout << "AllTime - " << (double)time / 1000 << "s" << std::endl;

	return 0;
}