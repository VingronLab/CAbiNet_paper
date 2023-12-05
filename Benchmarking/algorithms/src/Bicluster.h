#ifndef BICLUSTER_H_
#define BICLUSTER_H_

#include "global.h"
#include "Data.h"

using namespace std;

class Bicluster{
public:

	int seed_row;

	int * rows; //bitmap of the rows in the pattern
	int num_rows;

	vector<int> cols;
	int col_to_check;
	
public:
	Bicluster();
	void clear();

	int set_seed_row(int seed);
	int pick_seed_cols(Data &in_data, int seed_col_sz);

	int expand_rows(Data &in_data, double max_diff, double do_rate);
	void expand_cols(Data &in_data, int max_col_sz, double max_diff, double do_rate);

	void output_biclust(fstream& outfile);

};

#endif