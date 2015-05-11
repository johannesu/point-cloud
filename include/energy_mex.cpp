#include "PointCloud.h"

// Calls main_function with correct template arguments
#include "mex_wrapper.h"

template<typename Data_cost, typename Regularization_cost>
void main_function(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// Parse input
	ASSERT(nrhs == 10);
	ASSERT(nlhs <= 5);


	int curarg = 1;
	unsigned int dimensions = mxGetScalar(prhs[curarg++]);
	matrix<double> assignments(prhs[curarg++]);
	matrix<double> points(prhs[curarg++]);
	matrix<unsigned int> connectivity(prhs[curarg++]);
	matrix<double> connectivity_weights(prhs[curarg++]);
	double data_weight = mxGetScalar(prhs[curarg++]);
	double tol = mxGetScalar(prhs[curarg++]);
	double eps = mxGetScalar(prhs[curarg++]);
	bool verbose = mxGetScalar(prhs[curarg++]);

	int num_points = points.N;
	int num_pairwise = connectivity.N;

	matrix<double> U(1);
	matrix<double> B(1);
	matrix<double> E(1);

	U(0) = 0;
	B(0) = 0;

	// Data term
	for (int i = 0; i < num_points; ++i) {
		Data_cost data_cost(dimensions, data_weight, &points(0,i));
		U(0) += data_cost(&assignments(0,i));		
	}

	// Regularization term
	for (int i = 0; i < num_pairwise; ++i) {
		int source = connectivity(0,i);
		int target = connectivity(1,i);

		Regularization_cost regularization_cost(dimensions, connectivity_weights(i), tol, &points(0,source), &points(0,target), eps);
		B(0) += regularization_cost(&assignments(0,source), &assignments(0,target));
	}

	E(0) = U(0) + B(0);
	plhs[0] = E;
	plhs[1] = U;
	plhs[2] = B;

	// Return the cost for each term
	if (nlhs > 3) {
		matrix<double> UT(num_points,1);
		matrix<double> PT(num_pairwise,1);

		// Data term
		for (int i = 0; i < num_points; ++i) {
			Data_cost data_cost(dimensions, data_weight, &points(0,i));
			UT(i) = data_cost(&assignments(0,i));		
		}

		// Regularization term
		for (int i = 0; i < num_pairwise; ++i) {
			int source = connectivity(0,i);
			int target = connectivity(1,i);

			Regularization_cost regularization_cost(dimensions, connectivity_weights(i), tol, &points(0,source), &points(0,target), eps);
			PT(i) = regularization_cost(&assignments(0,source), &assignments(0,target));
		}

		plhs[3] = UT;
		plhs[4] = PT;
	}
}