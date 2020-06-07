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

void sequence(vec_type &v, float &dest)
{
    dest = 0;
    for (size_t i = 0; i < v.size(); i++) {
        dest += v[i];
    }
}

void shift(vec_type &v, float &dest)
{
    dest = 0;
    int shift = 111;
    int shift_count = 10;
    for (size_t i = 0; i < (v.size()-shift)/shift_count; i++) {
        for (size_t j = 0; j < shift_count; j++)
            dest += v[i + shift*j];
    }
}

void rand(vec_type &v, vec_int_type &ind, float &dest)
{
    dest = 0;
    for (size_t i = 0; i < v.size(); i++) {
        dest += v[ind[i]];
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
	vec_int_type ind(countOfElement, 0);
	int int_max = std::numeric_limits<int>::max();
	std::mt19937 mt_rand(time(0));
	for (int i = 0; i < countOfElement; i++) {
		vec[i] = (T)(mt_rand()/(T)int_max);
	}
	for (int i = 0; i < countOfElement; i++) {
		ind[i] = mt_rand()%countOfElement;
	}
	float res;
	std::ofstream fout("test");
    {
	    high_resolution_clock::time_point _t0{high_resolution_clock::now()};
        //=================
	    sequence(vec, res);
	    //=================
	    high_resolution_clock::time_point _t1{high_resolution_clock::now()};
	    std::clog << res << std::endl;
	    auto dt = duration_cast<microseconds>(_t1 - _t0);
	    fout << dt.count() << std::endl;
    }
    {
	    high_resolution_clock::time_point _t0{high_resolution_clock::now()};
        //=================
	    shift(vec, res);
	    //=================
	    high_resolution_clock::time_point _t1{high_resolution_clock::now()};
	    std::clog << res << std::endl;
	    auto dt = duration_cast<microseconds>(_t1 - _t0);
	    fout << dt.count() << std::endl;
    }
    {
	    high_resolution_clock::time_point _t0{high_resolution_clock::now()};
        //=================
	    rand(vec, ind, res);
	    //=================
	    high_resolution_clock::time_point _t1{high_resolution_clock::now()};
	    std::clog << res << std::endl;
	    auto dt = duration_cast<microseconds>(_t1 - _t0);
	    fout << dt.count() << std::endl;
    }
	return 0;
}
