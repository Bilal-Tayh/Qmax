#ifndef QMAX_H
#define QMAX_H
#include <fstream>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include <iomanip> 
#include "RngFast.hpp"

class QMax
{
	int *_A;
	int _curIdx;
	int _q;
	int _actualsize;
	int _actualsizeMinusOne;
	int _qMinusOne;
	float _gamma;
	int _nminusq;
	int _phi;
    std::default_random_engine _generator;
	void maintenance();
	int PartitionAroundPivot(int left, int right, int pivot_idx, int* nums);
    double _delta;
    double _alpha;
    double _psi;
    int _K;
    int _Z;
    unsigned int * RandByteArray;
    rng::rng128 gen_arr;
    int counter;
    int bytecounter;
    int randArrLen;
    int bitsNum;
    int bytesNum;
    int mask;
    __m256i rand_bits;
public:
	void reset();
	int findKthLargestAndPivot();
	QMax(int q, float gamma);
	void insert(int id);
	int* largestQ();
	void print();
    int checkPivot(int value);
    int findValueIndex(int value);
    int GenerateRandom(int max);
    int GenerateRandom1(int max);
    uint32_t mm256_extract_epi32_var_indx(int i);
};
#endif

