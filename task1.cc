// Compile with 'g++ -march=native task1.cc -O4'

#include <iostream>
#include <thread>
#include <cmath>
#include <immintrin.h>

static const uint32_t xorshift_seed = 2463534242;

static const double eps = 1e-10;

inline double xorshift_rand(uint32_t &x) {
	x ^= x >> 13;
	x ^= x << 17;
	x ^= x >> 5;
	return 1. * x / ((uint32_t)(-1));
}

template <class T>
void do_not_opt_out(T &&x) {
	static auto ttid = std::this_thread::get_id();
	if (ttid == std::thread::id()) {
		const auto* p = &x;
		putchar(*reinterpret_cast<const char*>(p));
		std::abort();
	}
}

template <typename T>
double benchmark(T fn) {
	double voxel_size = 3;

	double corner[3];
	corner[0] = 10;
	corner[1] = 11;
	corner[2] = 12;

	double values[2][2][2][3];
	double dd = 100;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				values[i][j][k][0] = dd++;
				values[i][j][k][1] = dd++;
				values[i][j][k][2] = dd++;
			}
		}
	}

	size_t nr_runs = 1e8;
	double result = 0;

	auto t1 = std::chrono::high_resolution_clock::now();

	uint32_t rand_state = xorshift_seed;

	for (size_t i = 0; i < nr_runs; ++i) {
		double point[3];
		point[0] = corner[0] + voxel_size * xorshift_rand(rand_state);
		point[1] = corner[1] + voxel_size * xorshift_rand(rand_state);
		point[2] = corner[2] + voxel_size * xorshift_rand(rand_state);

		double r[3];
		fn(r, point, corner, voxel_size, values);

		result += r[0];
		result += r[1];
		result += r[2];
	}

	auto t2 = std::chrono::high_resolution_clock::now();
	auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

	do_not_opt_out(result);

	return dt * 1e-3;
}

template <typename T>
bool compare(T fn1, T fn2) {
	double voxel_size = 3;

	double corner[3];
	corner[0] = 10;
	corner[1] = 11;
	corner[2] = 12;

	double values[2][2][2][3];
	double dd = 100;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			for (int k = 0; k < 2; ++k) {
				values[i][j][k][0] = dd++;
				values[i][j][k][1] = dd++;
				values[i][j][k][2] = dd++;
			}
		}
	}

	double input[][3] = {
		{0, 0, 0},
		{0, 0, 1},
		{0, 1, 0},
		{0, 1, 1},
		{1, 0, 0},
		{1, 0, 1},
		{1, 1, 0},
		{1, 1, 1},
		{0.1, 0.1, 0.1},
		{0.1, 0.1, 0.8},
		{0.1, 0.8, 0.1},
		{0.1, 0.8, 0.8},
		{0.8, 0.1, 0.1},
		{0.8, 0.1, 0.8},
		{0.8, 0.8, 0.1},
		{0.8, 0.8, 0.8},
	};

	for (size_t i = 0; i < sizeof(input) / sizeof(input[0]); ++i) {
		double point[3];
		for (int j = 0; j < 3; ++j)
			point[j] = corner[j] + voxel_size * input[i][j];

		double r1[3];
		fn1(r1, point, corner, voxel_size, values);

		double r2[3];
		fn2(r2, point, corner, voxel_size, values);

		if (fabs(r1[0] - r2[0]) < eps && fabs(r1[1] - r2[1]) < eps && fabs(r1[2] - r2[2]) < eps)
			continue;

		std::cerr << "\nDifferent values\n";
		std::cerr << "at:     " << input[i][0] << " " << input[i][1] << "  " << input[i][2] << "\n";
		std::cerr << "point:  " << point[0] << " " << point[1] << "  " << point[2] << "\n";
		std::cerr << "r1:     " << r1[0] << " " << r1[1] << "  " << r1[2] << "\n";
		std::cerr << "r2:     " << r2[0] << " " << r2[1] << "  " << r2[2] << "\n";

		return false;
	}

	return true;
}

void tri_interp_original(
	double result[3],
	const double point[3],
	const double corner[3],
	double voxel_size,
	const double values[2][2][2][3]
) {
	double p[3];
	for (int i = 0; i < 3; ++i)
		p[i] = (point[i] - corner[i]) / voxel_size;

	for (int i = 0; i < 3; ++i) {
		result[i] = 0;
		result[i] += values[0][0][0][i] * (1 - p[0]) * (1 - p[1]) * (1 - p[2]);
		result[i] += values[0][0][1][i] * (1 - p[0]) * (1 - p[1]) * (    p[2]);
		result[i] += values[0][1][0][i] * (1 - p[0]) * (    p[1]) * (1 - p[2]);
		result[i] += values[0][1][1][i] * (1 - p[0]) * (    p[1]) * (    p[2]);
		result[i] += values[1][0][0][i] * (    p[0]) * (1 - p[1]) * (1 - p[2]);
		result[i] += values[1][0][1][i] * (    p[0]) * (1 - p[1]) * (    p[2]);
		result[i] += values[1][1][0][i] * (    p[0]) * (    p[1]) * (1 - p[2]);
		result[i] += values[1][1][1][i] * (    p[0]) * (    p[1]) * (    p[2]);
	}
}



#pragma GCC push_options
#pragma GCC optimize ("-ffast-math")

inline void my_tri_interp(
	double result[3],
	const double point[3],
	const double corner[3],
	double voxel_size,
	const double values[2][2][2][3]
) {
	double p[3];
	for (int i = 0; i < 3; ++i)
		p[i] = (point[i] - corner[i]) / voxel_size;
	//==============================================================
	// ready 2.26 sec
	//__m256d p1_0 = _mm256_set1_pd(1-p[0]);
	//__m256d p0_0 = _mm256_set1_pd(p[0]);
	//__m256d p1_1 = _mm256_set_pd(1-p[1], 1-p[1], p[1], p[1]) * _mm256_set_pd(1-p[2], p[2], 1-p[2], p[2]);
	//__m256d coef1 = p1_0*p1_1;
	//__m256d coef2 = p0_0*p1_1;
	//__m256i mask = _mm256_setr_epi64x(-1, -1, -1, 1);
	//__m256d result_256 = _mm256_setzero_pd();
   	//__m256d v1 = _mm256_set_pd(values[0][0][0][0], values[0][0][1][0], values[0][1][0][0], values[0][1][1][0]);
   	//__m256d v2 = _mm256_set_pd(values[1][0][0][0], values[1][0][1][0], values[1][1][0][0], values[1][1][1][0]);
   	//__m256d v3 = _mm256_set_pd(values[0][0][0][1], values[0][0][1][1], values[0][1][0][1], values[0][1][1][1]);
   	//__m256d v4 = _mm256_set_pd(values[1][0][0][1], values[1][0][1][1], values[1][1][0][1], values[1][1][1][1]);
   	//__m256d v5 = _mm256_set_pd(values[0][0][0][2], values[0][0][1][2], values[0][1][0][2], values[0][1][1][2]);
   	//__m256d v6 = _mm256_set_pd(values[1][0][0][2], values[1][0][1][2], values[1][1][0][2], values[1][1][1][2]);
	//v1 = v1*coef1 + v2*coef2;
	//v3 = v3*coef1 + v4*coef2;
	//v5 = v5*coef1 + v6*coef2;
	//result[0] = v1[0] + v1[1] + v1[2] + v1[3];
	//result[1] = v3[0] + v3[1] + v3[2] + v3[3];
	//result[2] = v5[0] + v5[1] + v5[2] + v5[3];
	//==============================================================
	// ready 2.552 sec
	//__m256d p1_0 = _mm256_set1_pd(1-p[0]);
	//__m256d p0_0 = _mm256_set1_pd(p[0]);
	//__m256d p1_1 = _mm256_set_pd(p[1], p[1], 1-p[1], 1-p[1]) * _mm256_set_pd(p[2], 1-p[2], p[2], 1-p[2]);
	//__m256d coef1 = p1_0*p1_1;
	//__m256d coef2 = p0_0*p1_1;
	//__m256i mask = _mm256_setr_epi64x(-1, -1, -1, 1);
	//__m256d result_256 = _mm256_setzero_pd();
   	//__m256d v000 = _mm256_maskload_pd((double*)values , mask);
   	//__m256d v001 = _mm256_maskload_pd((double*)values + 3, mask);
   	//__m256d v010 = _mm256_maskload_pd((double*)values + 6, mask);
   	//__m256d v011 = _mm256_maskload_pd((double*)values + 9, mask);
   	//__m256d v100 = _mm256_maskload_pd((double*)values + 12, mask);
   	//__m256d v101 = _mm256_maskload_pd((double*)values + 15, mask);
   	//__m256d v110 = _mm256_maskload_pd((double*)values + 18, mask);
   	//__m256d v111 = _mm256_maskload_pd((double*)values + 21, mask);
	//result_256 = 
	//	v000*coef1[0] + 
	//	v001*coef1[1] + 
	//	v010*coef1[2] + 
	//	v011*coef1[3] + 
	//	v100*coef2[0] + 
	//	v101*coef2[1] + 
	//	v110*coef2[2] + 
	//	v111*coef2[3];
	//_mm256_maskstore_pd((double*)result ,mask, result_256);
	//==============================================================
	//ready 2.408 sec
	//double p1 = p[0];
	//double p2 = p[1];
	//double p3 = p[2];
	//double p12 = p1*p2;
	//double p13 = p1*p3;
	//double p23 = p2*p3;
	//double p123 = p12*p3;
	//double coef[8];
	//coef[0] = (1 - p[0]) * (1 - p[1]) * (1 - p[2]);
	//coef[1] = (1 - p[0]) * (1 - p[1]) * (    p[2]);
	//coef[2] = (1 - p[0]) * (    p[1]) * (1 - p[2]);
	//coef[3] = (1 - p[0]) * (    p[1]) * (    p[2]);
	//coef[4] = (    p[0]) * (1 - p[1]) * (1 - p[2]);
	//coef[5] = (    p[0]) * (1 - p[1]) * (    p[2]);
	//coef[6] = (    p[0]) * (    p[1]) * (1 - p[2]);
	//coef[7] = (    p[0]) * (    p[1]) * (    p[2]);
	//__m256i mask = _mm256_setr_epi64x(-1, -1, -1, 1);
	//__m256d result_256 = _mm256_setzero_pd();
	//for (int i = 0; i < 8; i++) {
   	//	__m256d val = _mm256_maskload_pd((double*)values + i*3, mask);
	//	result_256 = (val * coef[i]) + result_256;
	//}
	//_mm256_maskstore_pd((double*)result ,mask, result_256);
	//=====================================================
	//// ready 1.994 sec
	//__m256i mask = _mm256_setr_epi64x(-1, -1, -1, 1);
	//__m256d result_256 = _mm256_setzero_pd();
   	//__m256d v000 = _mm256_maskload_pd((double*)values , mask);
   	//__m256d v001 = _mm256_maskload_pd((double*)values + 3, mask);
   	//__m256d v010 = _mm256_maskload_pd((double*)values + 6, mask);
   	//__m256d v011 = _mm256_maskload_pd((double*)values + 9, mask);
   	//__m256d v100 = _mm256_maskload_pd((double*)values + 12, mask);
   	//__m256d v101 = _mm256_maskload_pd((double*)values + 15, mask);
   	//__m256d v110 = _mm256_maskload_pd((double*)values + 18, mask);
   	//__m256d v111 = _mm256_maskload_pd((double*)values + 21, mask);
	//result_256 = ((v000*(1-p[2]) + v001*(p[2]))*(1-p[1]) + (v010*(1-p[2]) + v011*(p[2]))*(p[1]))*(1-p[0]) +
	//		((v100*(1-p[2]) + v101*(p[2]))*(1-p[1]) + (v110*(1-p[2]) + v111*(p[2]))*(p[1]))*(p[0]);
	//_mm256_maskstore_pd((double*)result ,mask, result_256);
	
	//=====================================================
	//ready 2.001 sec
	//__m256i mask = _mm256_setr_epi64x(-1, -1, -1, 1);
	//__m256d result_256 = _mm256_setzero_pd();
   	//__m256d v000 = _mm256_maskload_pd((double*)values , mask);
   	//__m256d v001 = _mm256_maskload_pd((double*)values + 3, mask);
   	//__m256d v010 = _mm256_maskload_pd((double*)values + 6, mask);
   	//__m256d v011 = _mm256_maskload_pd((double*)values + 9, mask);
   	//__m256d v100 = _mm256_maskload_pd((double*)values + 12, mask);
   	//__m256d v101 = _mm256_maskload_pd((double*)values + 15, mask);
   	//__m256d v110 = _mm256_maskload_pd((double*)values + 18, mask);
   	//__m256d v111 = _mm256_maskload_pd((double*)values + 21, mask);
	//__m256d v00 = v000 + p[0]*(v100 - v000);
	//__m256d v01 = v001 + p[0]*(v101 - v001);
	//__m256d v10 = v010 + p[0]*(v110 - v010);
	//__m256d v11 = v011 + p[0]*(v111 - v011);
	//__m256d v0 = v00 + p[1]*(v10 - v00);
	//__m256d v1 = v01 + p[1]*(v11 - v01);
	//result_256 = v0 + p[2]*(v1 - v0);
	//_mm256_maskstore_pd((double*)result ,mask, result_256);
	//=====================================================
	//ready 1.831 sec
	
	__m256d p0 = _mm256_set1_pd(p[0]);
	__m256d p1 = _mm256_set1_pd(p[1]);
	__m256d p2 = _mm256_set1_pd(p[2]);

	__m256i mask = _mm256_setr_epi64x(-1, -1, -1, 1);
	__m256d result_256 = _mm256_setzero_pd();
	
   	__m256d v000 = _mm256_maskload_pd((double*)values , mask);
   	__m256d v001 = _mm256_maskload_pd((double*)values + 3, mask);
    __m256d v00 = v000 + p2*(v001 - v000);
   	
	__m256d v010 = _mm256_maskload_pd((double*)values + 6, mask);
   	__m256d v011 = _mm256_maskload_pd((double*)values + 9, mask);
    __m256d v01 = v010 + p2*(v011 - v010);
   	
	__m256d v100 = _mm256_maskload_pd((double*)values + 12, mask);
   	__m256d v101 = _mm256_maskload_pd((double*)values + 15, mask);
    __m256d v10 = v100 + p2*(v101 - v100);
   	
	__m256d v110 = _mm256_maskload_pd((double*)values + 18, mask);
   	__m256d v111 = _mm256_maskload_pd((double*)values + 21, mask);
    __m256d v11 = v110 + p2*(v111 - v110);
   
	__m256d v0 = v00 + p1*(v01 - v00);
	__m256d v1 = v10 + p1*(v11 - v10);
	
	result_256 = v0 + p0*(v1 - v0);
	
	_mm256_maskstore_pd((double*)result ,mask, result_256);
	//=====================================================
	//ready 2.594 sec
	//double p1 = p[0];
	//double p2 = p[1];
	//double p3 = p[2];
	//double p12 = p1*p2;
	//double p13 = p1*p3;
	//double p23 = p2*p3;
	//double p123 = p12*p3;
	//double coef[8];

	//coef[7] = p123;
	//coef[6] = p12 - p123;
	//coef[5] = p13 - p123;
	//coef[4] = p1 - p12  - coef[5];
	//coef[3] = p23 - p123;
	//coef[2] = p2 - p12 - coef[3];
	//coef[1] = p3 - p13 - coef[3];
	//coef[0] = 1 - p1 -p2 -p3 + p12 + p13 + p23 - p123;
	//
	//for (int i = 0; i < 3; ++i) {
	//	result[i] = 0;
	//	result[i] += values[0][0][0][i] * coef[0];
	//	result[i] += values[0][0][1][i] * coef[1];
	//	result[i] += values[0][1][0][i] * coef[2];
	//	result[i] += values[0][1][1][i] * coef[3];
	//	result[i] += values[1][0][0][i] * coef[4];
	//	result[i] += values[1][0][1][i] * coef[5];
	//	result[i] += values[1][1][0][i] * coef[6];
	//	result[i] += values[1][1][1][i] * coef[7];
	//}
	//=====================================================
	// ready 2.344 sec
	//for (int i = 0; i < 3; ++i) {
	//	result[i] = 0;
	//	result[i] += values[0][0][0][i] * (1 - p[0]) * (1 - p[1]) * (1 - p[2]);
	//	result[i] += values[0][0][1][i] * (1 - p[0]) * (1 - p[1]) * (    p[2]);
	//	result[i] += values[0][1][0][i] * (1 - p[0]) * (    p[1]) * (1 - p[2]);
	//	result[i] += values[0][1][1][i] * (1 - p[0]) * (    p[1]) * (    p[2]);
	//	result[i] += values[1][0][0][i] * (    p[0]) * (1 - p[1]) * (1 - p[2]);
	//	result[i] += values[1][0][1][i] * (    p[0]) * (1 - p[1]) * (    p[2]);
	//	result[i] += values[1][1][0][i] * (    p[0]) * (    p[1]) * (1 - p[2]);
	//	result[i] += values[1][1][1][i] * (    p[0]) * (    p[1]) * (    p[2]);
	//}
}

#pragma GCC pop_options

int main()
{
	std::cerr << "tri_interp_original: ";
	std::cerr << benchmark(tri_interp_original); 
	std::cerr << " sec\n";

	std::cerr << "my_tri_interp: ";
	std::cerr << benchmark(my_tri_interp); 
	std::cerr << " sec\n";

	if (!compare(tri_interp_original, my_tri_interp))
		std::cerr << "not valid\n";
	else
		std::cerr << "valid\n";
}

