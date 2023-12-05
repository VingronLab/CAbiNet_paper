/*
* DivBiclust
* A biclustering-based framwork for identifying
* cell subpopulations based on scRNA-Seq data
*
* Programming Language: C++ 
* Platform: Visual Studio 2015
*
* Author: Qiong Fang, Dewei Su 
* Contact: fang.qiong@gmail.com, 
*
* Updated on 5 November 2019
*
*/

#include "global.h"
#include "Data.h"
#include "BicSet.h"
#include "Eval.h"
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;


// // difference threshold for biclustering
// double max_diff = 0.15; 

// //fraction of missing values in a bicluster
// double do_rate = 0.1;

// // size of seed gene set
// int seed_col_sz = 50;  

// // maximum size of gene set, fixed to 100
// int max_col_sz = 100; 

// //similarity threshold for pattern merging, fixed to 0.5
double sim_thresh = 0.5; 

    // string ds_type = "../data-hcc/HPIC";// HEC, HPIC, HUM1, HUM2, HUM3
string in_fn; //name the input gene-by-cell matrix file as "HCC_matrix.txt";
string gt_fn; //name the ground truth file as "HCC_gt.txt";
string out_fn; // output file "HCC_output.txt";
int gene_num, cell_num;
Data in_data; //input cell-by-gene data matrix

// double set_global_var(double sim_thresh){
//     return sim_thresh;
// }

// sim_thresh = set_global_var()

struct MyVariables {
    double ari;
    int nclust;
};

// ' @export
// [[Rcpp::export]]
std::string DivBiclust(string ds_type,
		 string in_file,
         double max_diff = 0.15, 
         double do_rate = 0.1, 
         int seed_col_sz = 50, 
         int max_col_sz = 100,
         double simthresh = 0.5){
         // double sim_thresh = 0.5
        // ) {
	clock_t start, end;
	double ari_score;
    
    sim_thresh = simthresh;
    
    // cout << "The value of sim_thresh is " << sim_thresh << "\n";

	fstream out_file;

	// if (argc != 1 && argc != 5) {
	// 	cout << "Usage: DivBi data_type max_diff seed_col_sz do_rate" << endl;
	// 	cout << "data_type: HCC, HPIC, HUM1, ...\n";
	// 	cout << "max_diff: maximum expression difference allowed for biclustering\n";
	// 	cout << "seed_col_sz: size of seed column set for biclustering\n";
	// 	cout << "do_rate: dropout rate \n";
	// 	exit(0);
	// }

	// if (argc == 5) {
	// 	ds_type = string(argv[1]);
	// 	max_diff = atof(argv[2]);
	// 	seed_col_sz = atoi(argv[3]);
	// 	do_rate = atof(argv[4]);
	// }

	in_fn = in_file + "_matrix.txt";
	gt_fn = in_file + "_gt.txt";
	out_fn = ds_type + "_output.txt";

	out_file.open(out_fn, ios::out);
	if (!out_file.is_open())
		cout << "Fail to open output file: " << out_fn << endl;

	cout << "data_file = " << in_fn << endl;
	cout << "gt_file = " << gt_fn << endl;
	cout << "output_file = " << out_fn << endl;
	cout << "max_diff = " << max_diff << endl;
	cout << "seed_col_sz = " << seed_col_sz << endl;
	cout << "do_rate = " << do_rate << endl;

	out_file << "data_file = " << in_fn << endl;
	out_file << "gt_file = " << gt_fn << endl;
	out_file << "max_diff = " << max_diff << "\tseed_col_sz = " << seed_col_sz
		<< "\tdo_rate = " << do_rate << endl << endl;

	start = clock();

	in_data.read_data(in_fn, cell_num, gene_num);


	BicSet bic_set;
    int nclust;
	//mine biclusters
	bic_set.mine_biclusts(in_data, seed_col_sz, max_col_sz, max_diff, do_rate, out_file);

	// EDITED
	if (bic_set.ptn_set.size() == 0) {
		
    	std::string concatenated_str = "NA_0";
		return concatenated_str;
	}


	//merge biclusters to generate cell clusters 
	bic_set.merge_biclusts(out_file);

	//post-hoc assignment of unclustered cells
	nclust = bic_set.assign(in_data, out_file, do_rate);

	//compute ARI score
	ari_score = eval_ARI(gt_fn, bic_set);
	out_file << "ARI = " << ari_score << endl;

	end = clock();
	cout << "runtime = " << (float)(end - start) / CLOCKS_PER_SEC << " sec" << endl << endl;
	out_file << "runtime = " << (float)(end - start) / CLOCKS_PER_SEC << " sec" << endl << endl;

	//output the clustering results to file
	bic_set.output_ptn_set(out_file);
	bic_set.print_ptn_set();

	out_file.close();
    
    // MyVariables result;
    // result.ari = ari_score;
    // result.nclust = nclust;
    
    std::string concatenated_str = std::to_string(ari_score) + "_" + std::to_string(nclust);

	return concatenated_str;
}
