#include "QmaxNew.hpp"
#include <iostream>
#include <random>
#include <math.h> 
#include <queue>
#include <vector>

#include <bitset>
#include <unistd.h>
#include <cstring>
#include <algorithm>
#include <memory.h>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include <iomanip> 
#include "RngFast.hpp"



void QMaxNew::print(){
	for (int i = 0; i < _actualsize; ++i)
		std::cout << _A[i] << " ";
	std::cout << std::endl;
	std::cout << "_phi  = " << _phi << std::endl;
}

QMaxNew::QMaxNew(int q, float gamma){
	_actualsize = q * (1+gamma);
	_actualsizeMinusOne = _actualsize - 1;
	_curIdx = 0;
	_A = (int*) malloc(sizeof(int) * _actualsize);
	if (!_A)
		exit(1);
	_gamma = gamma;
	_q = q;
	_qMinusOne = q - 1;
	_nminusq = _actualsize - q;
	_phi = -1;
    _delta = 1.0-0.999;
    _alpha=0.35;
    _psi = 2.0/3.0;
    _K=ceil( ((_alpha*_gamma*(2+_gamma - _alpha*_gamma)) / (pow(_gamma -_alpha*_gamma,2))) * log(1/_delta));
    _Z = (int)  ( (_K*(1+_gamma)) / (_alpha*_gamma));
    gen_arr();
    rand_bits = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
    counter=0;
}




void QMaxNew::insert(int v){
	if (v < _phi){
		return;
	}
	else{
        while(_A[_curIdx]>=_phi){
            _curIdx++;
            if(_curIdx==_actualsize){
                    maintenance();
            }
        }
        _A[_curIdx] = v;
	}
}







void QMaxNew::maintenance(){
    int ii=_phi;
	_phi = findKthLargestAndPivot();
 	_curIdx = 0;
}

int* QMaxNew::largestQ(){
	maintenance();
	return _A + _nminusq;
}

inline void swap(int &x, int &y){
	int z = x;
	x = y;
	y = z;
}
int QMaxNew::PartitionAroundPivot(int left, int right, int pivot_idx, int* nums) {
	int pivot_value = nums[pivot_idx];
	int new_pivot_idx = right;
	swap(nums[pivot_idx], nums[right]);
	for (int i = right-1; i >= left; --i) {
		if (nums[i] > pivot_value) {
			swap(nums[i], nums[--new_pivot_idx]);
		}
	}
	swap(nums[right], nums[new_pivot_idx]);
	return new_pivot_idx;
}






uint32_t QMaxNew::mm256_extract_epi32_var_indx(int i){   
        __m128i indx = _mm_cvtsi32_si128(i);
        __m256i val  = _mm256_permutevar8x32_epi32(rand_bits, _mm256_castsi128_si256(indx));
        return         _mm_cvtsi128_si32(_mm256_castsi256_si128(val));
}  






int QMaxNew::GenerateRandom(int max){
    int indx = 0;
    if(counter >= 8){
        rand_bits =_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
        counter=0;
    }
    indx = mm256_extract_epi32_var_indx(counter)%max;
    counter++;
    return indx;
}



// if value dont exist return -1
int QMaxNew::findValueIndex(int value){
    for(int i=0;i<_actualsize;i++){
        if(_A[i]==value){
            return i;
        }
    }
    return -1;
}







// check if the conditions holds for the possible pivot "value"
// return the pivot index in _A if it hold otherwise return -1
int QMaxNew::checkPivot(int value){
//     int index=0;
//     int bigger=0;
//     int smaller=0;
//     for(int i=0;i<_actualsize;i++){
//         if(_V[i]==value){
//             index=i;
//         }
//         else if(_V[i]>value){
//             bigger++;
//         }
//         else smaller++;
//     } 
//     
    
    
    
    
    int left = 0, right = _actualsizeMinusOne;
    int pivot_idx = findValueIndex(value);
    
    return pivot_idx;
    
//      if(pivot_idx==-1){
//         return -1;
//     }
//     
// //     int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _A);
//     if (_actualsize - pivot_idx <= _nminusq) {
//         if(pivot_idx >= _gamma*_q*_psi) {  // pivot_idx < _q - 1.
//             return pivot_idx;
//         }
//     }
//     return -1;
}











int QMaxNew::findKthLargestAndPivot(){
    
    _delta=_delta/2.0;
    _K=ceil( ((_alpha*_gamma*(2+_gamma - _alpha*_gamma)) / (pow(_gamma -_alpha*_gamma,2))) * log(1/_delta));
    _Z = (int)  ( (_K*(1+_gamma)) / (_alpha*_gamma) );

    int tries=2;
    while(tries!=0){
        // B should contain Z random values from _A
        std::priority_queue <int, std::vector<int>, std::greater<int> > p;
        for(int i=0;i<_K;i++){
            int j=GenerateRandom(_actualsize);
            p.push(_A[j]);
        }
        int top = p.top();
        
        for(int i=_K;i<_Z;i++){
            int j=GenerateRandom(_actualsize);
            if(top<_A[j]){
                p.pop();
                p.push(_A[j]);  
                top = p.top();
            }
            
        }
    
        
        int idx = checkPivot(top);
        if(idx!=-1){
            return _A[idx];   
        }
        // if the conditions dont hold try sample Z elemnts from _A again...
        tries--;
    }
    int left = 0, right = _actualsizeMinusOne;
	while (left <= right) {
		int pivot_idx = left;
		int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _A);
		if (new_pivot_idx == _nminusq) {
			return _A[new_pivot_idx];
		} else if (new_pivot_idx > _nminusq) {
			right = new_pivot_idx - 1;
		} else {  // new_pivot_idx < _q - 1.
			left = new_pivot_idx + 1;
		}
	}
	
}

void QMaxNew::reset(){
	_phi = -1;
	_curIdx = _actualsize;
}

