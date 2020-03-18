#include "Qmax.hpp"
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



void QMax::print(){
	for (int i = 0; i < _actualsize; ++i)
		std::cout << _A[i] << " ";
	std::cout << std::endl;
	std::cout << "_phi  = " << _phi << std::endl;
}

QMax::QMax(int q, float gamma){
	_actualsize = q * (1+gamma);
	_actualsizeMinusOne = _actualsize - 1;
	_curIdx = _actualsize;
	_A = (int*) malloc(sizeof(int) * _actualsize);
	if (!_A)
		exit(1);
	_gamma = gamma;
	_q = q;
	_qMinusOne = q - 1;
	_nminusq = _actualsize - q;
	_phi = -1;
    _delta = 1.0-0.999;
    _alpha=0.83;
    _psi = 2.0/3.0;
    _K=ceil( ((_alpha*_gamma*(2+_gamma - _alpha*_gamma)) / (pow(_gamma -_alpha*_gamma,2))) * log(1/_delta));
    _Z = (int)  ( (_K*(1+_gamma)) / (_alpha*_gamma) );
    gen_arr();
    rand_bits =_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
    RandByteArray=(unsigned char *)&rand_bits;
    counter=0;
    bytecounter=0;
    bitsNum = std::floor(log2(_actualsize))+1;
    bytesNum = std::max((int)ceil((bitsNum-1)/sizeof(int))-1,1);
    mask = std::pow(2,bitsNum)-1;
}
void QMax::insert(int v){
	if (v < _phi){
		return;
	}
	else{
		_A[--_curIdx] = v;
		if (_curIdx){
			return;
		}
		else {
			maintenance();
		}
	}
}

void QMax::maintenance(){
    int ii=_phi;
	_phi = findKthLargestAndPivot();
    if(ii>_phi){
        std::cout<<"wrong"<<std::endl;
    }
	_curIdx = _nminusq;
}

int* QMax::largestQ(){
	maintenance();
	return _A + _nminusq;
}

inline void swap(int &x, int &y){
	int z = x;
	x = y;
	y = z;
}
int QMax::PartitionAroundPivot(int left, int right, int pivot_idx, int* nums) {
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





/*
int QMax::GenerateRandom(int min,int max){
  // uniform_real_distribution documentation
  // http://www.cplusplus.com/reference/random/uniform_real_distribution/
  std::uniform_real_distribution<double> distribution(0.0,1.0);
  double u = distribution(_generator);
  int A = max-min+1;
  for(int i=0;i<A;i++){
      if(u<(double)(i+1)/(double)A){
         return min+i;
      }
  }
  
  return 0;
}*/




int QMax::GenerateRandom(int max){
  int indx = 0;
  do{
    
    if((bytecounter+bytesNum) > 31){
        rand_bits =_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
        RandByteArray=(unsigned char *)&rand_bits;
        bytecounter=0;
    }
    unsigned int l=0;
    std::memcpy(&l,&(RandByteArray[bytecounter]),bytesNum);
    bytecounter += bytesNum;
    indx = l&mask;

  }while(indx>max);

   return indx;
}



/*
int QMax::GenerateRandom(int max){
  // uniform_real_distribution documentation
  // http://www.cplusplus.com/reference/random/uniform_real_distribution/
//   std::uniform_real_distribution<double> distribution(min,max+1);
//   double u = distribution(_generator);
  
//   return int(u);
  

  
  int indx = 0;
  do{
    indx = 0;
    for(int i=0;i<bitsNum;i++){
        if(((RandByteArray[counter]>>bytecounter)&1)==1){
            indx*=2;
            indx+=1;
        }
        else{
            indx*=2;
        }
        bytecounter++;
        if(bytecounter == 8){
            bytecounter=0;
            counter++;
            if(counter==31){
                __m256i rand_bits =_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
                RandByteArray=(unsigned char *)&rand_bits;
                counter=0;
            }
        }
    }
  }while(indx>max);

   return indx;
}
*/





// if value dont exist return -1
int QMax::findValueIndex(int value){
    for(int i=0;i<_actualsize;i++){
        if(_A[i]==value){
            return i;
        }
    }
    return -1;
}







// check if the conditions holds for the possible pivot "value"
// return the pivot index in _A if it hold otherwise return -1
int QMax::checkPivot(int value){
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
    
    
     if(pivot_idx==-1){
        return -1;
    }
    
    int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _A);
    
    if (_actualsize - new_pivot_idx <= _nminusq) {
        if(new_pivot_idx >= _gamma*_q*_psi) {  // new_pivot_idx < _q - 1.
            return new_pivot_idx;
        }
    }
    return -1;
}











int QMax::findKthLargestAndPivot(){
    
//     for(int i=0;i<100;i++){
//         std::cout<<gen_arr()<<std::endl;
//     }
    
    
    int tries=2;
    while(tries!=0){
        // B should contain Z random values from _A
        std::priority_queue <int, std::vector<int>, std::greater<int> > p;
        int top=0;
        for(int i=0;i<_Z;i++){
            int j=GenerateRandom(_actualsize);
            
            if(p.size()<_K-1){
                p.push(_A[j]);
            }
            else if(p.size()==_K-1){
                p.push(_A[j]);
                top = p.top();
            }
            else if(top<_A[j]){
                p.pop();
                p.push(_A[j]);  
                top = p.top();
            }
            
        }
    
        
        // check if the conditions holds for the possible pivot B[Kth_minimal_idx]
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

void QMax::reset(){
	_phi = -1;
	_curIdx = _actualsize;
}

