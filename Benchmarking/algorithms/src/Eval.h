#ifndef EVAL_H_
#define EVAL_H_

#include "global.h"
#include "Bicluster.h"
#include "BicSet.h"

using namespace std;

double eval_ARI(string gt_fn, BicSet clust_set);
double compute_ARI(BicSet clust_set, list<Bicluster> gt_set);

#endif

