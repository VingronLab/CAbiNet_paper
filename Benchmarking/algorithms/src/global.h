#ifndef GLOBAL_H_
#define GLOBAL_H_

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>

#include<vector>
#include<list>

#include<algorithm>

#include<cstdlib>
#include<cmath>
#include<ctime>

using namespace std; 

#define L_BUF_SIZE 1000000
#define S_BUF_SIZE 600

// Original values. Modified by Clemens Kohl
//#define L_BUF_SIZE 100000
//#define S_BUF_SIZE 200

extern int cell_num, gene_num;
extern double sim_thresh;

class ElemInt {
public:
	int id;
	int degree;

public:
	ElemInt(){}

	ElemInt(const ElemInt &elem) {
		id = elem.id;
		degree = elem.degree;
	}

	ElemInt(int i, int d) {
		id = i;
		degree = d;
	}

	const ElemInt operator=(const ElemInt &elem) {
		this->id = elem.id;
		this->degree = elem.degree;
		return *this;
	}

	ElemInt set_val(int i, int d) {
		this->id = i;
		this->degree = d;
		return *this;
	}
};


class ElemDbl {
public:
	int id;
	double val;

public:
	ElemDbl() {}

	ElemDbl(const ElemDbl &elem) {
		id = elem.id;
		val = elem.val;
	}

	ElemDbl(int i, double v) {
		id = i;
		val = v;
	}

	friend bool operator< (const ElemDbl &elem1, const ElemDbl &elem2) {
		return (elem1.val > elem2.val || (elem1.val == elem2.val && elem1.id < elem2.id));
	}
};

#endif
