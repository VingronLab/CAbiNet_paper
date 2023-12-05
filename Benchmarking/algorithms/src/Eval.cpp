//compute the ARI score
#include "Eval.h"

using namespace std;

double eval_ARI(string gt_fn, BicSet clust_set) {
	char char_buf[L_BUF_SIZE];
	istringstream in_str;
	int cell_id;

	fstream gt_file;
	double ARI_score;

	list<Bicluster> gt_set;

	gt_file.open(gt_fn, ios::in);
	if (!gt_file.is_open()) {
		cout << "Fail to open ground truth file" << endl;
		exit(0);
	}

	while (gt_file.getline(char_buf, L_BUF_SIZE)) {
		Bicluster gt_biclust;
		in_str.clear();
		in_str.str("");
		in_str.str(char_buf);
		while (!in_str.eof()) {
			in_str >> cell_id;
			gt_biclust.rows[cell_id] = 1;
			gt_biclust.num_rows += 1;
		}
		gt_set.push_back(gt_biclust);
	}
	gt_file.close();

	ARI_score = compute_ARI(clust_set, gt_set);

	return ARI_score;

}


double compute_ARI(BicSet clust_set, list<Bicluster> gt_set) {
	int i, j, ptn_id, gt_id;
	int ** cont_table; //contigency table with ptn_pool for X dim; gt_set for Y dim
	int * row_sum; // the row sum of contigency table
	int * col_sum; // the col_sum of contigency table
	int num_rows, num_cols, total_sum;

	double ari_score, idx, max_idx, exp_idx;

	num_rows = (int)clust_set.ptn_set.size();
	num_cols = (int)gt_set.size();
	cont_table = new int*[num_rows];
	for (i = 0; i < num_rows; i++) {
		cont_table[i] = new int[num_cols];
		for (j = 0; j < num_cols; j++)
			cont_table[i][j] = 0;
	}
	row_sum = new int[num_rows];
	for (i = 0; i < num_rows; i++)
		row_sum[i] = 0;
	col_sum = new int[num_cols];
	for (j = 0; j < num_cols; j++)
		col_sum[j] = 0;

	//construct the contigency table, row_sum and col_sum
	total_sum = 0;
	list<Pattern>::iterator ptn_it = clust_set.ptn_set.begin();
	for (ptn_id = 0; ptn_it != clust_set.ptn_set.end(); ptn_it++, ptn_id++) {
		list<Bicluster>::iterator gt_it = gt_set.begin();
		for (gt_id=0; gt_it != gt_set.end(); gt_it++, gt_id++) {
			for (i = 0; i < cell_num; i++) {
				if (ptn_it->cells[i] == 1 && gt_it->rows[i] == 1)
					cont_table[ptn_id][gt_id] += 1;
			}
			row_sum[ptn_id] += cont_table[ptn_id][gt_id];
			col_sum[gt_id] += cont_table[ptn_id][gt_id];
			total_sum += cont_table[ptn_id][gt_id];
		}
	}
	
	//compute adjusted rand index (ari)
	idx = 0.0;
	max_idx = 0.0;
	exp_idx = 0.0;
	for (i = 0; i < num_rows; i++) {
		for (j = 0; j < num_cols; j++)
			idx += (double)cont_table[i][j] * (cont_table[i][j] - 1) / 2;
	}

	int row_sq_sum = 0, col_sq_sum = 0; // the sum of (row_sum choose 2) and (col_sum choose 2)
	for (i = 0; i < num_rows; i++)
		row_sq_sum += row_sum[i] * (row_sum[i] - 1) / 2;
	for (j = 0; j < num_cols; j++)
		col_sq_sum += col_sum[j] * (col_sum[j] - 1) / 2;

	max_idx = (double)(row_sq_sum + col_sq_sum) / 2;
	exp_idx = (double)(row_sq_sum*col_sq_sum) / (total_sum*(total_sum - 1) / 2);

	ari_score = (idx - exp_idx) / (max_idx - exp_idx);
	cout << "ARI = " << ari_score << endl;

	for (i = 0; i < num_rows; i++) {
		delete cont_table[i];
		cont_table[i] = NULL;
	}
	delete cont_table;
	delete row_sum;
	delete col_sum;
	return ari_score;
}