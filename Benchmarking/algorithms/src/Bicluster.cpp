#include "Bicluster.h"

//int comparison
bool compare_int_inc(const int & a, const int & b) {
	return (a < b);
}



//Bicluster class
Bicluster::Bicluster() {
	rows = new int[cell_num];
	for (int i = 0; i < cell_num; i++)
		rows[i] = 0;
	num_rows = 0;
	return;
}

void Bicluster::clear() {
	for (int i = 0; i < cell_num; i++)
		rows[i] = 0;
	num_rows = 0;

	cols.clear();
	col_to_check = 0;

	seed_row = -1;

	return;
}

int Bicluster::set_seed_row(int seed) {
	seed_row = seed;
	rows[seed_row] = 1;
	num_rows += 1;

	return seed_row;
}


//pick seed columns which have the largest values in the seed row 
int Bicluster::pick_seed_cols(Data &in_data, int seed_col_sz) {

	int i;
	vector<ElemDbl> & seed_row_ptr = in_data.sorted_matrix[seed_row];
	vector<ElemDbl>::iterator it = seed_row_ptr.begin();
	
	for (i = 0; i < seed_col_sz&&it != seed_row_ptr.end(); i++) {
		cols.push_back(it->id);
		it++;
	}
	if (it != seed_row_ptr.end())
		col_to_check = i;
	else
		col_to_check = -1;
		
	return (int)cols.size();
} 



//expand the bicluster with the rows 
//which are similar with the seed row
int Bicluster::expand_rows(Data &in_data, double max_diff, double do_rate) {
	int i, j;
	double diff;

	int row_sz = in_data.get_row_sz();
	int num_cols = (int)cols.size();
	int max_bad_entries = (int)(num_cols*do_rate);
	int num_good_entries, num_bad_entries;

	for (i = 0; i < row_sz; i++) {
		if (rows[i] == 1) //row i is already in the bicluster.
			continue;

		num_good_entries = 0;
		num_bad_entries = 0;
		bool is_added = true;
		for (j = 0; j < num_cols; j++) {
			if (in_data.matrix[i][cols[j]] == 0) {
				num_bad_entries += 1;
				if (num_bad_entries > max_bad_entries) {
					is_added = false;
					break;
				}
				else
					continue;
			}
			diff = fabs(in_data.matrix[i][cols[j]] - in_data.matrix[seed_row][cols[j]]);
			diff = diff / in_data.matrix[seed_row][cols[j]];
			if (diff <= max_diff)
				num_good_entries += 1;
			else {
				num_bad_entries += 1;
				if (num_bad_entries > max_bad_entries) {
					is_added = false;
					break;
				}
			}
		} //end of checking current row

		if (is_added == true && num_good_entries >= num_cols*(1 - do_rate)) {
			rows[i] = 1;
			num_rows += 1;
		}
	} //end of expanding the current bicluster

	return (num_rows - 1);
}


//expand the bicluster with the columns 
//which satisfy both the similarity and dropout requirements
void Bicluster::expand_cols(Data& in_data, int max_col_sz, double max_diff, double do_rate) {
	int col, k;

	double diff;
	bool is_added;
	int num_zero_entries;

	if (col_to_check == -1) //no more columns for expansion
		return;

	vector<ElemDbl>::iterator it = in_data.sorted_matrix[seed_row].begin() + col_to_check;
	while(it != in_data.sorted_matrix[seed_row].end() && (int)cols.size()<max_col_sz){
		col = it->id; 
		it++;

		is_added = true;
		num_zero_entries = 0;
		for (k = 0; k < cell_num; k++) {
			if (rows[k] == 0)
				continue;
			if (in_data.matrix[k][col] == 0.0) {//ignore entries with value 0
				num_zero_entries += 1;
				if (num_zero_entries > num_rows*do_rate)
					break;
				else
					continue;
			}

			diff = fabs(in_data.matrix[k][col] - in_data.matrix[seed_row][col]);
			diff = diff / in_data.matrix[seed_row][col];
			if (diff > max_diff) {
				is_added = false;
				break;
			}
		}
		if (is_added == true && num_zero_entries <= num_rows*do_rate)
			cols.push_back(col);

	}
	sort(cols.begin(), cols.end(), compare_int_inc);
	return;
}

//output a bicluster to file
void Bicluster::output_biclust(fstream& outfile) {
	int i;
	outfile << "seed: " << seed_row << ' ' << num_rows << ' ' << cols.size()<< endl;
	outfile << "rows: ";
	for (i = 0; i < cell_num; i++) {
		if (rows[i] != 0 && rows[i] != -1)
			outfile << i << ' ';
	}
	outfile << "\ncolumns: ";
	for (i = 0; i < (int)cols.size(); i++)
		outfile << cols[i]<<' ';
	outfile << endl;
}