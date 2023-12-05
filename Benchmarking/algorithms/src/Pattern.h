//Convert mined biclusters to patterns which are ready for pattern merging
#ifndef PATTERN_H_
#define PATTERN_H_

using namespace std;

#include <map>
#include "global.h"
#include "Bicluster.h"

class Pattern {
public:
	static int ptn_type; // 0-adjacency matrix, 1-adjacency list

	int * cells;
	int num_cells;

	int ** adj_matrix;	
	map<int, vector<ElemInt>> adj_list;

	vector<int> genes;

	bool is_avail;


public:
	Pattern(); 
	Pattern(Bicluster biclust); //convert a bicluster to a pattern for pattern merging
	~Pattern();

	static void set_ptn_type(int type) {
		ptn_type = type;
	}

	const Pattern operator=(const Pattern& ptn);

	friend bool operator<(const Pattern &ptn1, const Pattern &ptn2);

	friend Pattern combine_ptn(Pattern &ptn1, Pattern &ptn2);
	friend double compute_p2p_sim(Pattern &ptn1, Pattern &ptn2);
	friend double compute_c2p_diff(Data &in_data, int cid, Pattern ptn, vector<double> mean_vec, double do_rate);
	
	int get_max_degree(int cid);
	vector<double> compute_mean_vec(Data &in_data);

	bool delete_cell(int cid);
	void delete_ptn();

	void reset_all();
	void reset_cells();

	void print_ptn() const;
	void output_ptn(fstream &outfile)  const;
};


#endif
