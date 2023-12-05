//Convert mined biclusters to patterns which are ready for pattern merging
#include "global.h"
#include "Pattern.h"

using namespace std;

int Pattern::ptn_type;

Pattern::Pattern() {
	
	cells = new int[cell_num];
	for(int i=0; i<cell_num; i++)
		cells[i] = 0;
	num_cells = 0;

	if (ptn_type == 0) {
		adj_matrix = new int*[cell_num];
		for (int i = 0; i < cell_num; i++) {
			adj_matrix[i] = new int[cell_num];
			for (int j = 0; j < cell_num; j++)
				adj_matrix[i][j] = 0;
		}
	}
	else
		adj_matrix = NULL;

	is_avail = true;
	return;
}

//create a pattern based on a bicluster
Pattern::Pattern(Bicluster biclust) {
	int i, j;

	num_cells = biclust.num_rows;
	cells = new int[cell_num];
	for (i = 0; i < cell_num; i++)
		cells[i] = biclust.rows[i];

	if (ptn_type == 0) { //adjacency matrix
		adj_matrix = new int*[cell_num];
		for (i = 0; i < cell_num; i++) {
			adj_matrix[i] = new int[cell_num];
			for (j = 0; j < cell_num; j++)
				adj_matrix[i][j] = 0;
		}
		for (i = 0; i < cell_num; i++) {
			if (cells[i] == 0)
				continue;
			adj_matrix[i][i] += 1;
			for (j = i + 1; j < cell_num; j++) {
				if (cells[j] == 0)
					continue;
				adj_matrix[i][j] += 1;
				adj_matrix[j][i] += 1;
			}
		}
	}
	else { //adjacency list
		adj_matrix = NULL;

		vector<ElemInt> elem_list;
		for (i = 0; i < cell_num; i++)
			if (cells[i] == 1){
				ElemInt elem(i, 1);
				elem_list.push_back(elem);
			}
		for (i = 0; i < cell_num; i++)
			if (cells[i] == 1) {
				vector<ElemInt> tmp_list(elem_list);
				adj_list[i] = tmp_list;
			}
	}

	genes.assign(biclust.cols.begin(), biclust.cols.end()); //copy the gene set
	
	is_avail = true;
	return;
}

Pattern::~Pattern() {
	return;
}

const Pattern Pattern::operator=(const Pattern& ptn)
{
	int i, j;

	this->num_cells = ptn.num_cells;
	for (i = 0; i < cell_num; i++)
		this->cells[i] = ptn.cells[i];

	if (ptn_type == 0) {
		for (i = 0; i < cell_num; i++)
			for (j = 0; j < cell_num; j++)
				this->adj_matrix[i][j] = ptn.adj_matrix[i][j];
	}
	else
		this->adj_list = ptn.adj_list;

	for (i = 0; i<(int)ptn.genes.size(); i++)
		this->genes.push_back(ptn.genes[i]);

	return *this;
}

//sort patterns in non-decreasing order of cell ids
bool operator<(const Pattern &ptn1, const Pattern &ptn2) {

	int cell_count = 0;
	for (int i = 0; i < cell_num; i++) {
		if (ptn1.cells[i] == 0 && ptn2.cells[i] == 0)
			continue;
		else if (ptn1.cells[i] == 1 && ptn2.cells[i] == 0)
			return true;
		else if (ptn1.cells[i] == 0 && ptn2.cells[i] == 1)
			return false;
		else //if (first.cells[i] == 1 && second.cells[i] == 1)
			cell_count += 1;

		if (cell_count == ptn1.num_cells && cell_count < ptn2.num_cells)
			return true;
		else if (cell_count < ptn1.num_cells && cell_count == ptn2.num_cells)
			return false;
	}
	return false;
}

//compute the similarity between two patterns
//for pattern merging
double compute_p2p_sim(Pattern &ptn1, Pattern &ptn2) {

	int intersect_num = 0;
	for (int i = 0; i < cell_num; i++) {
		if (ptn1.cells[i] == 1 && ptn2.cells[i] == 1)
			intersect_num += 1;
	}

	if (ptn1.num_cells <= ptn2.num_cells)
		return (double)intersect_num / ptn1.num_cells;
	else
		return (double)intersect_num / ptn2.num_cells;
}

//merge two patterns
Pattern combine_ptn(Pattern &ptn1, Pattern &ptn2) {

	int i, j;
	Pattern combined_ptn;

	//merge the cell set
	combined_ptn.num_cells = 0;
	for (i = 0; i < cell_num; i++) {
		if (ptn1.cells[i] == 1 || ptn2.cells[i] == 1) {
			combined_ptn.cells[i] = 1;
			combined_ptn.num_cells += 1;
		}
		else
			combined_ptn.cells[i] = 0;
	}

	//merge the adjacency matrix
	if (Pattern::ptn_type == 0) {
		for (i = 0; i < cell_num; i++) {
			for (j = 0; j < cell_num; j++)
				combined_ptn.adj_matrix[i][j] = ptn1.adj_matrix[i][j] + ptn2.adj_matrix[i][j];
		}
	}
	//merge the adjacency list
	else {
		map<int, vector<ElemInt>>::const_iterator it1 = ptn1.adj_list.begin();
		map<int, vector<ElemInt>>::const_iterator it2 = ptn2.adj_list.begin();
		while (it1 != ptn1.adj_list.end() && it2 != ptn2.adj_list.end()) {
			if (it1->first < it2->first) {
				vector<ElemInt> tmp_list(it1->second);
				combined_ptn.adj_list[it1->first] = tmp_list;
				it1++;
			}
			else if (it1->first > it2->first) {
				vector<ElemInt> tmp_list(it2->second);
				combined_ptn.adj_list[it2->first] = tmp_list;
				it2++;
			}
			else {
				vector<ElemInt> tmp_list;
				vector<ElemInt>::const_iterator it11 = it1->second.begin();
				vector<ElemInt>::const_iterator it22 = it2->second.begin();
				while (it11 != it1->second.end() && it22 != it2->second.end()) {
					if (it11->id < it22->id) {
						ElemInt elem(it11->id, it11->degree);
						tmp_list.push_back(elem);
						it11++;
					}
					else if (it11->id > it22->id) {
						ElemInt elem(it22->id, it22->degree);
						tmp_list.push_back(elem);
						it22++;
					}
					else {//it11->id==it22->id
						ElemInt elem(it11->id, it11->degree + it22->degree);
						tmp_list.push_back(elem);
						it11++; it22++;
					}
				}
				while(it11!=it1->second.end()){
					ElemInt elem(it11->id, it11->degree);
					tmp_list.push_back(elem);
					it11++;
				}
				while (it22 != it2->second.end()) {
					ElemInt elem(it22->id, it22->degree);
					tmp_list.push_back(elem);
					it22++;
				}
				combined_ptn.adj_list[it1->first] = tmp_list;
				it1++; it2++;
			}
		}
		while (it1 != ptn1.adj_list.end()) {
			vector<ElemInt> tmp_list(it1->second);
			combined_ptn.adj_list[it1->first] = tmp_list;
			it1++;
		}
		while (it2 != ptn2.adj_list.end()) {
			vector<ElemInt> tmp_list(it2->second);
			combined_ptn.adj_list[it2->first] = tmp_list;
			it2++;
		}
	}

	//merge two gene sets by taking their intersection
	for (i = 0, j = 0; i < (int)ptn1.genes.size() && j < (int)ptn2.genes.size();) {
		if (ptn1.genes[i] < ptn2.genes[j])
			i++;
		else if (ptn1.genes[i] > ptn2.genes[j])
			j++;
		else {
			combined_ptn.genes.push_back(ptn1.genes[i]);
			i++; j++;
		}
	}
	return combined_ptn;
}


//remove a cell from a pattern
bool Pattern::delete_cell(int cid) {
	if (cells[cid] == 0)
		return false;
	cells[cid] = 0;
	num_cells -= 1;
	if (num_cells > 1) { //update adjacency info
		if (ptn_type == 0) {
			for (int i = 0; i < cell_num; i++) {
				adj_matrix[cid][i] = 0;
				adj_matrix[i][cid] = 0;
			}
		}
		else {
			adj_list.erase(cid);
			map<int, vector<ElemInt>>::iterator it1 = adj_list.begin();
			for (; it1 != adj_list.end();it1++) {
				vector<ElemInt>::iterator it2 = it1->second.begin();
				for (; it2 != it1->second.end();it2++) {
					if (it2->id == cid) {
						it1->second.erase(it2);
						break;
					}
				}
			}
		}
	}
	return true;
}

//delete a pattern and release the
//space allocated for it
void Pattern::delete_ptn() {

	if (ptn_type == 0) {
		for (int i = 0; i < cell_num; i++) {
			delete adj_matrix[i];
			adj_matrix[i] = NULL;
		}
		delete adj_matrix;
	}
	else {
		map<int, vector<ElemInt>>::iterator it = adj_list.begin();
		for (;it != adj_list.end(); it++)
			it->second.clear();
		adj_list.clear();
	}

	genes.clear();
	delete cells;

	return;
}

//return the largest degree between cid and any other cell in the pattern 
int Pattern::get_max_degree(int cid) {
	int curr_degree;
	if (ptn_type == 0) {
		curr_degree = 0;
		for (int i = 0; i < cell_num; i++) {
			if (cells[i] == 0)
				continue;
			if (curr_degree < adj_matrix[cid][i])
				curr_degree = adj_matrix[cid][i];
		}
	}
	else {
		curr_degree = 0;
		vector<ElemInt>::iterator it = adj_list[cid].begin();
		for (; it != adj_list[cid].end(); it++)
			if (curr_degree < it->degree)
				curr_degree = it->degree;
	}
	return curr_degree;
}

//compute the mean vector of a pattern
vector<double> Pattern::compute_mean_vec(Data &in_data) {
	int i, j;
	vector<double> mean_vec;
	for (i = 0; i < (int)genes.size(); i++) {
		double avg_val = 0;
		int num_nonzeros = 0;
		for (j = 0; j < cell_num; j++) {
			if (cells[j] == 1 && in_data.matrix[j][genes[i]] != 0){
				avg_val += in_data.matrix[j][genes[i]];
				num_nonzeros += 1;
			}
		}
		if (num_nonzeros != 0)
			avg_val = avg_val / num_nonzeros;
		else
			avg_val = 0.0;
		mean_vec.push_back(avg_val);
	}
	return mean_vec;
}


//compute the distance between a cell and a pattern
double compute_c2p_diff(Data &in_data, int cid, Pattern ptn, vector<double> mean_vec, double do_rate) {

	int i;
	int num_nonzeros = 0;
	double diff, diff_sum=0.0;

	vector<double>::iterator it = mean_vec.begin();
	for (i = 0; i < (int)ptn.genes.size(); i++) {
		if (in_data.matrix[cid][ptn.genes[i]] != 0) {
			num_nonzeros += 1;
			diff = (fabs(in_data.matrix[cid][ptn.genes[i]] - it[i]))/it[i];
			diff_sum += diff;
		}
	}
	if (num_nonzeros > 0) {
		diff_sum = diff_sum / (num_nonzeros*num_nonzeros);
		return diff_sum;
	}
	else
		return 1000.0;
}


//reset the bitmap cells
void Pattern::reset_cells() {
	for (int i = 0; i < cell_num; i++)
		cells[i] = 0;
	num_cells = 0;
	return;
}

//reset cells, genes, and the adjacency matrix/list to be 0
void Pattern::reset_all() {
	for (int i = 0; i < cell_num; i++)
		cells[i] = 0;
	num_cells = 0;
	if (ptn_type == 0) {
		for (int i = 0; i < cell_num; i++)
			for (int j = 0; j < cell_num; j++)
				adj_matrix[i][j] = 0;
	}
	else {
		map<int, vector<ElemInt>>::iterator it = adj_list.begin();
		while (it != adj_list.end())
			it->second.clear();
		adj_list.clear();
	}
	
	genes.clear();
	return;
}

void Pattern::output_ptn(fstream &outfile) const {

	outfile << "(#cells = " << num_cells << ") ";
	for (int i = 0; i < cell_num; i++)
		if (cells[i] == 1)
			outfile << i << ' ';
	outfile << endl;

	return;
}

void Pattern::print_ptn() const {
	cout << "(#cells = " << num_cells << ") ";
	for (int i = 0; i < cell_num; i++)
		if (cells[i] == 1)
			cout << i << ' ';
	cout << endl;
}