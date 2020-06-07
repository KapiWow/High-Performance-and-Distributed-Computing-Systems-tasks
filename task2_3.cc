#include <stdlib.h>
#include <omp.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <cfloat>
#include <climits>
#include <sstream>
#include <random>
#include <functional>

typedef float T;

inline T f(T a, T b) {
    return sin(a) + sin(b);
}

int main(int argc, char** argv)
{
	using namespace std::chrono;
    using func = std::function<T(T,T)>;
    std::stringstream ss;
	
	int countOfElement;
	ss << argv[1];
	ss >> countOfElement;
	int threads_count = 8;
	
	std::mt19937 mt_rand(time(0));

	int int_max = INT_MAX;
	srand(clock());
	T *arr = new T[countOfElement];
	for (int i = 0; i < countOfElement; i++)
		arr[i] = (T)(mt_rand()/(T)int_max)/2;
	T *arr2 = new T[countOfElement];
	for (int i = 0; i < countOfElement; i++)
		arr2[i] = (T)(mt_rand()/(T)int_max)/2;

	std::ofstream fout("test");

	int time;
	
	for (int count = 1; count <= threads_count; count++) {
		T sum = 0;
		omp_set_num_threads(count);
		high_resolution_clock::time_point _t0{high_resolution_clock::now()};
		if (count > 1) {
#pragma omp parallel shared(arr, arr2, sum)
			{
				T pSum = 0;
#pragma omp for 
				for (int i = 0; i < countOfElement; i++) {
					pSum += f(arr[i], arr2[i]);
                }
#pragma omp critical
				sum += pSum;
			}
		} else {
			for (int i = 0; i < countOfElement; i++)
			    sum += f(arr[i], arr2[i]);
		}
		high_resolution_clock::time_point _t1{high_resolution_clock::now()};
        std::cout << sum << std::endl;
		auto dt = duration_cast<microseconds>(_t1 - _t0);
		time = dt.count();

		fout << time << std::endl;
	}

}
