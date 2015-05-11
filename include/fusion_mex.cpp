#include "PointCloud.h"

// Calls main_function with correct template arguments
#include "mex_wrapper.h"
#include "QPBO-v1.3.src/QPBO.h"

// Catch errors
static void erfunc(char *err) {
	mexErrMsgTxt(err);
}

template<typename Data_cost, typename Regularization_cost>
void main_function(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// Parse input
	ASSERT(nrhs == 12);
	ASSERT(nlhs == 4);

	int curarg = 1;
	unsigned int dimensions = mxGetScalar(prhs[curarg++]);
	matrix<double> assignments(prhs[curarg++]);
	matrix<double> proposal(prhs[curarg++]);
	matrix<double> points(prhs[curarg++]);
	matrix<unsigned int> connectivity(prhs[curarg++]);
	matrix<double> connectivity_weights(prhs[curarg++]);
	double data_weight = mxGetScalar(prhs[curarg++]);
	double tol = mxGetScalar(prhs[curarg++]);
	double eps = mxGetScalar(prhs[curarg++]);
	bool verbose = mxGetScalar(prhs[curarg++]);
	bool improve = mxGetScalar(prhs[curarg++]);

	int num_points = points.N;
	int num_pairwise = connectivity.N;

	// Create problem instance
	QPBO<double> graph(num_points, num_pairwise, erfunc);
	graph.AddNode(num_points);


	// Data term
	for (int i = 0; i < num_points; ++i) {
		Data_cost data_cost(dimensions, data_weight, &points(0,i));
		graph.AddUnaryTerm(i, data_cost(&assignments(0,i)), data_cost(&proposal(0,i)));
	}

	// Regularization term
	for (int i = 0; i < num_pairwise; ++i) {
		int source = connectivity(0,i);
		int target = connectivity(1,i);

		Regularization_cost regularization_cost(dimensions, connectivity_weights(i), tol, &points(0,source), &points(0,target), eps);

		graph.AddPairwiseTerm(source, target,
				regularization_cost(&assignments(0,source)	, &assignments(0,target)),
				regularization_cost(&assignments(0,source)	, &proposal(0,target)),
				regularization_cost(&proposal(0,source)	, &assignments(0,target)),
				regularization_cost(&proposal(0,source)	, &proposal(0,target))
				);
	}

	// Merge edges
	graph.MergeParallelEdges();

	// Solve for optimimum and label weak persitencies
	graph.Solve();
	graph.ComputeWeakPersistencies();

	// Defining outout
	matrix<double> labelling(num_points);
	matrix<double> energy(1);
	matrix<double> lower_bound(1);
	matrix<double> num_unlabelled(1);

	plhs[0] = labelling;
	plhs[1] = energy;
	plhs[2] = lower_bound;
	plhs[3] = num_unlabelled;

	// Count unlabelled before improve
	num_unlabelled(0) = 0;
	for (int i = 0; i < num_points; ++i)
	{
		if ((graph.GetLabel(i) < 0)) {
			num_unlabelled(0)++;
		}
	}

	// if any labels are unlabelled run improve
	if (improve && (num_unlabelled(0) > 0)) {
		graph.Improve();
	}

	for (int i = 0; i < num_points; ++i) {
			labelling(i) = graph.GetLabel(i);
	}

	energy(0) = graph.ComputeTwiceEnergy()/2;
	lower_bound(0) = graph.ComputeTwiceLowerBound()/2;
}