#include "QmaxKV.hpp"
#include <iostream>
#include <random>
#include <math.h> 
#include <queue>
#include <vector>

#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <immintrin.h>
#include <iomanip> 
#include "RngFast.hpp"

#include <tgmath.h>


#define READFROM(data,index,size) (((1<<size))-1 << (index))

int bitExtracted(int arr, int k ,int p){
    return (((1<<k)-1)&(arr>>(p-1)));
}


QMaxKV::QMaxKV(int q, float gamma) {
	_actualsize = q * (1+gamma);
	_actualsizeMinusOne = _actualsize - 1;
	_curIdx = _actualsize;
	_K = (key*) malloc(sizeof(key) * _actualsize);
	_V = (val*) malloc(sizeof(val) * _actualsize);
	if (!_V || !_K) {
		exit(1);
	}
	_gamma = gamma;
	_q = q;
	_qMinusOne = q - 1;
	_nminusq = _actualsize - q;
	_phi = 0;
	_k_phi = 0;
    _delta = 1.0-0.999;
    _alpha=0.83;
    _psi = 2.0/3.0;
    _k=ceil( ((_alpha*_gamma*(2+_gamma - _alpha*_gamma)) / (pow(_gamma -_alpha*_gamma,2))) * log(1/_delta));
    _Z = (int)  ( (_k*(1+_gamma)) / (_alpha*_gamma) );
    gen_arr();
    __m256i rand_bits =_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
    RandByteArray=(char *)&rand_bits;
    counter=0;
    bitcounter=0;
}

QMaxKV::~QMaxKV() {
	free(_K);
	free(_V);
}

val QMaxKV::getMinimalVal() {
	return _phi;
}

key QMaxKV::getMinimalKey() {
	return _k_phi;
}

void QMaxKV::update(key k, val v) {
	for(int i=0; i < _actualsize; ++i) {
		if(_K[i] == k) {
			_V[i] = v;
			return;
		}
	}
}

void QMaxKV::insert(key k, val v) {
	if (v < _phi) { // DIRECT COMPARSION: compare if new val is smaller than smallest val
		return;
	} else {
		_K[--_curIdx] = k;
		_V[_curIdx] = v;
		if (_curIdx){
			return;
		} else {
			maintenance();
		}
	}
}

void QMaxKV::maintenance() {
	_phi = findKthLargestAndPivot();
	_curIdx = _nminusq;
}

outputkv QMaxKV::largestQ() {
	outputkv out;
	maintenance();
	out.keyArr = _K + _nminusq;
	out.valArr = _V + _nminusq;
	return out;
}

inline void QMaxKV::swap(int a, int b) {
	if (a==b) return;
	key k = _K[a];
	_K[a] = _K[b];
	_K[b] = k;
	val v = _V[a];
	_V[a] = _V[b];
	_V[b] = v;
}



inline void QMaxKV::swap1(int a, int b, val* num) {
	if (a==b) return;
	val v = num[a];
	num[a] = num[b];
	num[b] = v;
}



int QMaxKV::PartitionAroundPivot(int left, int right, int pivot_idx, val* nums) {
	val pivot_value = nums[pivot_idx];
	int new_pivot_idx = right;
	swap(pivot_idx, right);
	for (int i = right-1; i >= left; --i) {
		if (nums[i] > pivot_value) { // DIRECT COMPARSION: compare if nums[i] is bigger than pivot_value
			swap(i, --new_pivot_idx);
		}
	}
	swap(right, new_pivot_idx);
//     std::cout<<_V[new_pivot_idx] <<" " << _V[right]<<" " << _V[left]<<std::endl;
	return new_pivot_idx;
}
    
    
    
    
    
    
    
    
int QMaxKV::PartitionAroundPivot1(int left, int right, int pivot_idx, val* nums) {
	val pivot_value = nums[pivot_idx];
	int new_pivot_idx = right;
	swap1(pivot_idx, right,nums);
	for (int i = right-1; i >= left; --i) {
		if (nums[i] > pivot_value) { // DIRECT COMPARSION: compare if nums[i] is bigger than pivot_value
			swap1(i, --new_pivot_idx,nums);
		}
	}
	swap1(right, new_pivot_idx,nums);
	return new_pivot_idx;
}




int QMaxKV::GenerateRandom(int max){
  // uniform_real_distribution documentation
  // http://www.cplusplus.com/reference/random/uniform_real_distribution/
//   std::uniform_real_distribution<double> distribution(min,max+1);
//   double u = distribution(_generator);
  
//   return int(u);
  
  
  int bitsNum = std::floor(log2(max))+1;
  
  int indx = 0;
  
//   std::cout<<"*********"<< bitsNum<<"**********" <<std::endl;
  do{
//        std::cout<<"do"<<std::endl;
    indx = 0;
    for(int i=0;i<bitsNum;i++){
//         std::cout<< ((RandByteArray[counter]>>bitcounter)&1)<<std::endl;
        if(((RandByteArray[counter]>>bitcounter)&1)==1){
            indx*=2;
            indx+=1;
        }
        else{
            indx*=2;
        }
        bitcounter++;
        if(bitcounter == 8){
            bitcounter=0;
            counter++;
            if(counter==32){
                __m256i rand_bits =_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
                RandByteArray=(char *)&rand_bits;
                counter=0;
            }
        }
    }
  }while(indx>max);
//    std::cout<< indx <<std::endl;
  
    

//    std::cout<< RandByteArray[rand_pos >> 3] & (1 << ((rand_pos & 0x7) - 1)) <<std::endl;
   return indx;
}





int QMaxKV::findValueIndex(val value){
    for(int i=0;i<_actualsize;i++){
        if(_V[i]==value){
            return i;
        }
    }
    return -1;
}





// check if the conditions holds for the possible pivot "value"
// return the pivot index in _V if it hold otherwise return -1
int QMaxKV::checkPivot(val value){
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
    
    int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _V);
    
    
//         std::cout<<"*****" <<new_pivot_idx <<"*******"<<std::endl;
//         std::cout<<"bigger= "<< bigger << ",   n-q= "<< _nminusq<<std::endl;
//         std::cout<<"smaller= "<< smaller << ",   should = "<< _gamma*_q*_psi<<std::endl;
//     
    
    if (_actualsize - new_pivot_idx <= _nminusq) {
        if(new_pivot_idx >= _gamma*_q*_psi) {  // new_pivot_idx < _q - 1.
            return new_pivot_idx;
        }
    }
    return -1;
}







/*
val QMaxKV::findKthLargestAndPivot() {
	int left = 0, right = _actualsizeMinusOne;
	while (left <= right) {
		int pivot_idx = _q;
		int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _V);
		if (new_pivot_idx == _nminusq) {
			_k_phi = _K[new_pivot_idx];
			return _V[new_pivot_idx];
		} else if (new_pivot_idx > _nminusq) {
			right = new_pivot_idx - 1;
		} else {  // new_pivot_idx < _q - 1.
			left = new_pivot_idx + 1;
		}
	}
}*/





val QMaxKV::findKthLargestAndPivot() {

    int tries=2;
    while(tries!=0){
        
        // B should contain Z random values from _V
        
    std::priority_queue <val, std::vector<val>, std::greater<val> > p;
        for(int i=0;i<_Z;i++){
            int j=GenerateRandom(_actualsize);
            
            if(p.size()<_k){
                p.push(_V[j]);
            }
            else if(p.top()<_V[j]){
                p.pop();
                p.push(_V[j]);   
            }
            
        }
    
        
        // check if the conditions holds for the possible pivot B[Kth_minimal_idx]
        int idx = checkPivot(p.top());
        if(idx!=-1){
            return _V[idx];   
        }
        // if the conditions dont hold try sample Z elemnts from _V again...
        tries--;
    }
//     std::cout<<"Miss"<<std::endl;

    
 int left = 0, right = _actualsizeMinusOne;
	while (left <= right) {
		int pivot_idx = left;
		int new_pivot_idx = PartitionAroundPivot(left, right, pivot_idx, _V);
		if (new_pivot_idx == _nminusq) {
			return _V[new_pivot_idx];
		} else if (new_pivot_idx > _nminusq) {
			right = new_pivot_idx - 1;
		} else {  // new_pivot_idx < _q - 1.
			left = new_pivot_idx + 1;
		}
	}
	
}







/*
void QMaxKV::checkCorrectness(){
    std::cout<<_q<<std::endl;
    std::cout<<_p.size()<<std::endl;
    int i=0;
    while(!_p.empty()){
     if(findValueIndex(_p.top())==-1){
        i++;
     }
        _p.pop();
    }
     std::cout<<i<<std::endl;
      std::cout<<"Miss"<<std::endl;
    
}*/



