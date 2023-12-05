#include "global.h"
#include "BicSet.h"

using namespace std;

BicSet::BicSet() {
	num_covered_cells = 0;
	covered_cells = new char[cell_num];
	for (int i = 0; i < cell_num; i++)
		covered_cells[i] = 0;
	return;
}

//mine biclusters
void BicSet::mine_biclusts(Data & in_data, int & seed_col_sz, int & max_col_sz, double & max_diff, double &do_rate, fstream & outfile) {

	int num_rows = in_data.get_row_sz();

	int i, seed_row, num_expd_rows;

	Bicluster biclust;

	cout << "Phase 1: Pattern Mining" << endl;

	for (i = 0; i < num_rows; i++) {
		biclust.clear();

		//pick seed row(cell)
		seed_row = biclust.set_seed_row(i);

		//pick seed column(gene) set
		if (biclust.pick_seed_cols(in_data, seed_col_sz) == 0)
			continue;

        //search for rows similar with seed row
		num_expd_rows = biclust.expand_rows(in_data, max_diff, do_rate);
		if (num_expd_rows == 0)
			continue;
        //search for more columns with similar values in the selected rows
		biclust.expand_cols(in_data, max_col_sz, max_diff, do_rate);

		if (biclust.num_rows > 1 && insert_biclust(biclust) != true)
			cout << "Bicluster insertion fault!";
	}

	cout << "num_biclusters = " << ptn_set.size() << endl;
	cout << "num_covered_cells = " << num_covered_cells << endl;
	outfile << "Phase 1: Pattern Mining.\nnum_biclusters = " << ptn_set.size() << endl;
	outfile << "num_covered_cells = " << num_covered_cells << endl;
	if (ptn_set.size() == 0) {
		cout << "No pattern is mined" << endl;
		// EDITED
		//	exit(0);
	}
	return;
}


bool BicSet::insert_biclust(Bicluster biclust) {
	Pattern pattern(biclust);
	ptn_set.push_back(pattern);

	for (int i = 0; i < cell_num; i++) {
		if (biclust.rows[i] == 1 && covered_cells[i] == 0) {
			num_covered_cells += 1;
			covered_cells[i] = 1;
		}
	}
	return true;
}

//merge similar biclusters to generate cell clusters
void BicSet::merge_biclusts(fstream & outfile) {

	cout << "\nPhase 2: Pattern Merging\n";
	merge_patterns();
	/* When the dataset contains a large cell set and 
	* excessive number of biclusters are mined in the first phase
	*/
	//	merge_patterns_staged();
	filter_patterns();
	
	cout << "num_patterns = " << ptn_set.size() << endl;
	cout << "num_covered_cells = " << num_covered_cells << endl;

	outfile << "\nPhase 2: Pattern Merging\nnum_patterns = " << ptn_set.size() << endl;
	outfile << "num_covered_cells = " << num_covered_cells << endl;

	return;
}


void BicSet::merge_patterns() {

	compute_ptn_sims(sim_thresh);
	while (!ptn_sim_list.empty()) {
		PtnSim psim = ptn_sim_list.front();
		pop_heap(ptn_sim_list.begin(), ptn_sim_list.end());
		ptn_sim_list.pop_back();
		if (psim.ptn1->is_avail && psim.ptn2->is_avail) {
			Pattern new_ptn = combine_ptn(*psim.ptn1, *psim.ptn2);
			remove_ptn(psim.ptn1);
			remove_ptn(psim.ptn2);
			insert_ptn(new_ptn, sim_thresh);
		}
	}

	list<Pattern>::iterator pset_it = ptn_set.begin();
	while (pset_it != ptn_set.end()) {
		if (pset_it->is_avail == false)
			pset_it = ptn_set.erase(pset_it);
		else
			pset_it++;
	}
	return;
}

//construct a priority queue to store the indices of
//pattern pairs eligible for merging
void BicSet::compute_ptn_sims(double min_sim) {

	list<Pattern>::iterator pset_it1 = ptn_set.begin(), pset_it2;
	for (; pset_it1 != ptn_set.end();) {
		if (pset_it1->is_avail == false)
			pset_it1 = ptn_set.erase(pset_it1);
		else {
			pset_it2 = pset_it1; pset_it2++;
			for (; pset_it2 != ptn_set.end();) {
				if (pset_it2->is_avail == false)
					pset_it2 = ptn_set.erase(pset_it2);
				else {
					double sim = compute_p2p_sim(*pset_it1, *pset_it2);
					if (sim >= min_sim) {
						PtnSim psim(sim, pset_it1, pset_it2);
						ptn_sim_list.push_back(psim);
					}
					pset_it2++;
				}
			}
			pset_it1++;
		}
	}
	make_heap(ptn_sim_list.begin(), ptn_sim_list.end());
	return;
}


//insert a pattern to ptn_set
//update the priority queue ptn_sim_list
void BicSet::insert_ptn(Pattern ptn, double min_sim) {

	list<Pattern>::iterator new_ptr = ptn_set.end();
	new_ptr = ptn_set.insert(new_ptr, ptn);

	list<Pattern>::iterator pset_it = ptn_set.begin();
	for (int i = 0; i < (int)ptn_set.size() - 1; i++) {
		if (pset_it->is_avail) {
			double sim = compute_p2p_sim(*pset_it, *new_ptr);
			if (sim >= min_sim) {
				PtnSim psim(sim, pset_it, new_ptr);
				ptn_sim_list.push_back(psim);
				push_heap(ptn_sim_list.begin(), ptn_sim_list.end());
			}
		}
		pset_it++;
	}

	return;
}

//when a large number of biclusters are mined in the 1st phase
void BicSet::merge_patterns_staged() {
	double sim = 1.0, sim_steps[3] = { 0.3, 0.15, 0.5 }; //the sum of steps equals to (1-sim_thresh)
	for (int i = 0; i < 3; i++) {
		sim = sim - sim_steps[i];
		compute_ptn_sims(sim);
		while (!ptn_sim_list.empty()) {
			PtnSim psim = ptn_sim_list.front();
			pop_heap(ptn_sim_list.begin(), ptn_sim_list.end());
			ptn_sim_list.pop_back();
			if (psim.ptn1->is_avail && psim.ptn2->is_avail) {
				Pattern new_ptn = combine_ptn(*psim.ptn1, *psim.ptn2);
				remove_ptn(psim.ptn1);
				remove_ptn(psim.ptn2);
				insert_ptn(new_ptn, sim);
			}
		}
	}

	list<Pattern>::iterator pset_it = ptn_set.begin();
	while (pset_it != ptn_set.end()) {
		if (pset_it->is_avail == false)
			pset_it = ptn_set.erase(pset_it);
		else
			pset_it++;
	}
	return;
}

//if a cell appears in multiple clusters
//keep it in the cluster with which its adhesive score is maximal
void BicSet::filter_patterns() {

	ElemInt max_rec;
	for (int cid = 0; cid < cell_num; cid++) {
		ElemList reclist;
		reclist.get_recs(ptn_set, &max_rec, cid);
		if (reclist.recs.size() > 1) { // find the max degree
			keep_max_only(cid, reclist, max_rec);
		}
	}

	//remove patterns with <=1 elements
	list<Pattern>::iterator p_it = ptn_set.begin();
	for (p_it = ptn_set.begin(); p_it != ptn_set.end();) {
		if (p_it->num_cells <= 1) {
			for (int cid = 0; cid < cell_num; cid++) {
				if (p_it->cells[cid] == 1 && covered_cells[cid]==1) {
					covered_cells[cid] = 0;
					num_covered_cells -= 1;
					break;
				}
			}
			p_it = ptn_set.erase(p_it);
		}
		else 
			p_it++;
	}

	return;
}


//if cell cid appears in more than one pattern
//keep it in the pattern where it has the largest degree
void BicSet::keep_max_only(int cid, ElemList reclist, ElemInt max_rec) {

	ElemInt curr_rec;
	list<Pattern>::iterator ptr = ptn_set.begin();
	int pid = 0;
	while (reclist.recs.size() != 0) {
		curr_rec = reclist.get_front();
		if (curr_rec.id == max_rec.id)  //keep the item in the pattern
			continue;
		else { //remove item from the pattern with pid = curr_rec.pid
			while (pid < curr_rec.id && ptr != ptn_set.end()) {
				pid++;
				ptr++;
			}
			if (pid == curr_rec.id)
				ptr->delete_cell(cid);
		}
	}
	return;
}


void BicSet::remove_ptn(list<Pattern>::iterator & pset_it) {
	pset_it->delete_ptn();
	pset_it->is_avail = false;
}


//assign the unclustered cells to the nearest cluster
// bool BicSet::assign(Data & in_data, fstream & outfile, double & do_rate) {
int BicSet::assign(Data & in_data, fstream & outfile, double & do_rate) {

	if (num_covered_cells == cell_num)
		return false;

	cout << "\nPhase 3: Post-hoc Assignment" << endl;

	int i, j;
	vector<ElemInt> c2p_list;

	//compute the center of each pattern
	list<vector<double>> ptn_mean_vecs;
	list<Pattern>::iterator it1 = ptn_set.begin();
	for (; it1 != ptn_set.end(); it1++) {
		vector<double> mean_vec = it1->compute_mean_vec(in_data);
		ptn_mean_vecs.push_back(mean_vec);
	}

	for (int cid = 0; cid < cell_num; cid++) {
		if (covered_cells[cid] == 1) //skip covered cell
			continue;

		double diff, min_diff = 1000;
		int near_ptn = -1;
		list<vector<double>>::iterator it2 = ptn_mean_vecs.begin();
		it1 = ptn_set.begin();
		for (j = 0; it2 != ptn_mean_vecs.end(); it2++, it1++, j++) {
			diff = compute_c2p_diff(in_data, cid, *it1, *it2, do_rate);
			if (diff < min_diff) {
				min_diff = diff;
				near_ptn = j;
			}
		}
		//assign cell c_id to near_ptn;
		if (near_ptn != -1) {
			ElemInt elem(cid, near_ptn);
			c2p_list.push_back(elem);
		}
	}

	for (i = 0, it1 = ptn_set.begin(); i < (int)ptn_set.size(); i++, it1++) {
		vector<ElemInt>::iterator it3 = c2p_list.begin();
		while (it3 != c2p_list.end()) {
			if (it3->degree == i) {
				it1->cells[it3->id] = 1;
				it1->num_cells += 1;
				covered_cells[it3->id] = 1;
				num_covered_cells += 1;
				it3 = c2p_list.erase(it3);
			}
			else
				it3++;
		}
	}

	cout << "num_patterns = " << ptn_set.size() << endl;
	cout << "num_covered_cells = " << num_covered_cells << endl;
	outfile << "\nPhase 3: Post-hoc Assignment\nnum_patterns = " << ptn_set.size() << endl;
	outfile << "num_covered_cells = " << num_covered_cells << endl;
	return ptn_set.size();
}


void BicSet::print_ptn_set() const {
	cout << ptn_set.size() << " patterns " << endl;

	list<Pattern>::const_iterator p_it = ptn_set.begin();
	for (int i = 0; p_it != ptn_set.end(); p_it++, i++) {
		cout << "Pattern " << i << ": ";
		p_it->print_ptn();
	}
	return;
}


//output patterns to file
void BicSet::output_ptn_set(fstream &outfile) const {

	outfile << "num_patterns = " << ptn_set.size() << endl;
	list<Pattern>::const_iterator p_it = ptn_set.begin();
	for (int i = 0; p_it != ptn_set.end(); p_it++, i++) {
		outfile << "Pattern " << i << ": ";
		p_it->output_ptn(outfile);
	}
	return;
}

//return the front element in the list
ElemInt ElemList::get_front() {
	ElemInt curr_rec;
	if (!recs.empty()) {
		curr_rec = recs.front();
		recs.pop_front();
	}
	return curr_rec;
}


void ElemList::get_recs(list<Pattern> ptn_set, ElemInt * max_rec, int cid) {
	ElemInt curr_rec;
	max_rec->set_val(-1, -1);

	list<Pattern>::iterator ptr = ptn_set.begin();
	for (int pid=0; ptr != ptn_set.end(); ptr++, pid++) {
		if (ptr->cells[cid] == 0) continue; //cell cid does not in current pattern

		curr_rec.id = pid;
		curr_rec.degree = ptr->get_max_degree(cid);

		if (curr_rec.degree > max_rec->degree) {
			(*max_rec) = curr_rec;
		}
		recs.push_back(curr_rec);
	}

	return;
}
