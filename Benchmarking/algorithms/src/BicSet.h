#ifndef BICSET_H_
#define BICSET_H_

#include<queue>

#include "global.h"
#include "Bicluster.h"
#include "Pattern.h"

class PtnSim {
public:
	double sim;
	int sz_diff;

	list<Pattern>::iterator ptn1;
	list<Pattern>::iterator ptn2;

public:
	PtnSim(double s, list<Pattern>::iterator p1, list<Pattern>::iterator p2) {
		sim = s;
		sz_diff = abs(p1->num_cells - p2->num_cells);
		ptn1 = p1;
		ptn2 = p2;
	}

	friend bool operator< (const PtnSim &psim1, const PtnSim & psim2){
		if (psim1.sim != psim2.sim) {
			if (psim1.sim < psim2.sim)
				return true;
			else
				return false;
		}
		else {
			if (psim1.sz_diff > psim2.sz_diff)
				return true;
			else
				return false;
		}
	}
};


class ElemList {
public:
	list<ElemInt> recs;

public:
	void get_recs(list<Pattern> ptn_set, ElemInt *max_rec, int cid);
	ElemInt get_front();
};


class BicSet {
public:
	list<Pattern> ptn_set;

	//used to construct a priority queue for pattern merging
	vector<PtnSim> ptn_sim_list;

	//cells covered by at least one pattern
	char * covered_cells;
	int num_covered_cells;

public:
	BicSet();
	//mine biclusters
	void mine_biclusts(Data & in_data, int & seed_col_sz, int & max_col_sz, double & max_diff, double & do_rate, fstream & outfile);
	//merge similar patterns
	void merge_biclusts(fstream & outfile);
	//assign unclustered cells to the nearest cluster
	int assign(Data &matrix, fstream & outfile, double & do_rate);

private:
	bool insert_biclust(Bicluster ptn); //insert mined biclusters to ptn_set

	void merge_patterns(); //merge pairs of simiar patterns
	void merge_patterns_staged(); // for data with large cell set & excessive number of biclusters are mined in the 1st phase
	void filter_patterns(); //make each item appear in exactly one pattern

	//for merge_patterns
	void compute_ptn_sims(double min_sim);
	void insert_ptn(Pattern ptn, double min_sim);
	void remove_ptn(list<Pattern>::iterator &ptr);

	//for filter_patterns
	// keep the item in the pattern where its degree is the largest
	void keep_max_only(int cid, ElemList rlist, ElemInt max_rec);

public:
	void print_ptn_set() const;
	void output_ptn_set(fstream &out_file) const;
};

#endif

