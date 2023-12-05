#include "Data.h"
#include "Pattern.h"

using namespace std;

Data::Data() {
	row_sz = -1;
	col_sz = -1;
	matrix = NULL;
	return;
}

Data::~Data() { 
	return;
}


/* Allocate space for matrix */
void Data::init(int row_num, int col_num) {
	row_sz = row_num;
	col_sz = col_num;

	matrix = new double *[row_sz];
	for (int i = 0; i < row_sz; i++)
		matrix[i] = new double[col_sz];

	return;
}


//read the input gene-by-cell data and store as a cell-by-gene matrix
void Data::read_data(string in_fn, int &cell_num, int &gene_num) {
	char char_buf[L_BUF_SIZE], name_buf[S_BUF_SIZE];
	string str_buf;

	fstream in_file;

	in_file.open(in_fn, ios::in | ios::out);
	if (!in_file.is_open()) {
		cout << "Failed to open data file: " << in_fn << endl;
		exit(1);
	}

	//get the column and row number;
	in_file.getline(char_buf, L_BUF_SIZE);
	istringstream in_str(char_buf);
	in_str >> gene_num >> cell_num;


	init(cell_num, gene_num);

	if (cell_num < 50) //use adjacency matrix to store patterns
		Pattern::ptn_type = 0;
	else //use adjacency list to store patterns
		Pattern::ptn_type = 1;

	int row_count = 0;
	//read the matrix
	while (in_file.getline(char_buf, L_BUF_SIZE)) {
		in_str.clear();
		in_str.str("");
		in_str.str(char_buf);
		in_str.get(name_buf, S_BUF_SIZE, ',');
		in_str.ignore(1);

		int col_count = 0;
		while (in_str.get(name_buf, S_BUF_SIZE, ',')) {
			in_str.ignore(1);
			matrix[col_count][row_count] = atof(name_buf);
			col_count++;
		}
		row_count++;
	}

	// Print the final values
	cout << "Final gene_num: " << gene_num << endl;
	cout << "row_count: " << row_count << endl;

	if (row_count != gene_num) {
		cout << "row counts mismatch" << endl;
		exit(1);
	}
	in_file.close();

	//store a sorted version of the matrix, used for biclustering
	for (int i = 0; i < row_sz; i++) {
		vector<ElemDbl> row_elems;
		ElemDbl elem;
		for (int j = 0; j < col_sz; j++) {
			if (matrix[i][j] != 0) {
				elem.id = j;
				elem.val = matrix[i][j];
				row_elems.push_back(elem);
			}
		}
		sort(row_elems.begin(), row_elems.end());
		sorted_matrix.push_back(row_elems);
	}
	return;
}


int Data::get_row_sz() {
	return row_sz;
}

int Data::get_col_sz() {
	return col_sz;
}
