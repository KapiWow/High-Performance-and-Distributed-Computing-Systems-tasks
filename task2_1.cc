#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <limits>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include <random>

typedef float T;
typedef std::vector<T> vec_type;
typedef std::vector<int> vec_int_type;

#define COUNT 100000

void f(vec_type &v, float &dest, int countOfElement)
{
    dest = 0;
    int ind = 0;
    int tick_count = 8;
    int count = countOfElement/tick_count;
    for (size_t i = 0; i < COUNT; i++) {
        ind = (ind + int(dest)) % count;
        for (size_t j = 0; j < tick_count; j++)
            dest += v[ind +count*j];
        dest *= 0.9;
    }
}

int main(int argc, char** argv) {
	using namespace std::chrono;
    std::stringstream ss;

	int countOfElement;
	ss << argv[1];
	ss >> countOfElement;
	if (countOfElement == 0) 
        countOfElement = 10;

	vec_type vec(countOfElement, 0);
	vec_int_type ind(COUNT, 0);
	int int_max = std::numeric_limits<int>::max();
	std::mt19937 mt_rand(time(0));
	for (int i = 0; i < countOfElement; i++) {
		vec[i] = (T)(countOfElement*(mt_rand()/(T)int_max));
	}
	for (int i = 0; i < COUNT; i++) {
		ind[i] = mt_rand()%countOfElement;
	}
	float res;
	high_resolution_clock::time_point _t0{high_resolution_clock::now()};
    //=================
	f(vec, res, countOfElement);
	//=================
	high_resolution_clock::time_point _t1{high_resolution_clock::now()};
	std::clog << res << std::endl;
	auto dt = duration_cast<microseconds>(_t1 - _t0);
	std::ofstream fout("test");
	fout << dt.count();
	std::clog << "time: " << dt.count() << std::endl;
	return 0;
}
