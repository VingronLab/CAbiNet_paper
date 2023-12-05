#ifndef DATA_H_
#define DATA_H_

#include "global.h"

class Data {
private:
	int row_sz, col_sz; //the number of rows & columns in the matrix
public:
	double ** matrix; //cell-by-gene data matrix
	vector<vector<ElemDbl>> sorted_matrix; //store nonzero values and sort each row in decreasing order of the values;

public:
	Data();
	~Data();

	void init(int row, int col);
	void read_data(string in_fn, int & row_num, int & col_num);
	int get_row_sz();
	int get_col_sz();
};


#endif // !DATA_H_
