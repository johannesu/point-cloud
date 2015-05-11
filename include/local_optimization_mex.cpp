#include <vector>
#include "PointCloud.h"


void mex_log_function(const std::string& str)
{
	mexPrintf("%s\n", str.c_str());
}

// Calls main_function with correct template arguments
#include "mex_wrapper.h"

template<typename Data_cost, typename Regularization_cost>
void main_function(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// Parse input
	ASSERT(nrhs == 11);
	ASSERT(nlhs == 1);

	int curarg = 1;
	unsigned int dimensions = mxGetScalar(prhs[curarg++]);
	matrix<double> in_assignments(prhs[curarg++]);
	matrix<double> points(prhs[curarg++]);
	matrix<unsigned int> connectivity(prhs[curarg++]);
	matrix<double> connectivity_weights(prhs[curarg++]);
	double data_weight = mxGetScalar(prhs[curarg++]);
	double tol = mxGetScalar(prhs[curarg++]);
	double eps = mxGetScalar(prhs[curarg++]);
	bool verbose = mxGetScalar(prhs[curarg++]);
	int max_iterations = mxGetScalar(prhs[curarg++]);

	int num_points = points.N;
	int num_pairwise = connectivity.N;

	// Copy values and then optimize on the copied values
	matrix<double> assignments(dimensions+1,num_points);

	for (int i = 0; i < in_assignments.numel(); ++i) {
		assignments(i) = in_assignments(i);
	}

	//

	// Function to be optimized
	spii::Function f;

	if (dimensions == 2) {
		// Setup the data term and variables
		for (int i = 0; i < num_points; ++i) {
			auto data_term = std::make_shared<spii::AutoDiffTerm<Data_cost, 3>>(dimensions, data_weight, &points(0,i));
			f.add_variable(&assignments(0,i), 3);
			f.add_term(data_term, &assignments(0,i));
		}

		// Regularization term
		for (int i = 0; i < num_pairwise; ++i) {
			int source = connectivity(0,i);
			int target = connectivity(1,i);

			auto regularization_term =
				std::make_shared<spii::AutoDiffTerm<Regularization_cost, 3, 3>>
				(dimensions, connectivity_weights(i), tol, &points(0,source), &points(0,target), eps);

			f.add_term(regularization_term,&assignments(0,source), &assignments(0,target));
		}
	} else {

		// Setup the data term and variables
		for (int i = 0; i < num_points; ++i) {
			auto data_term = std::make_shared<spii::AutoDiffTerm<Data_cost, 4>>(dimensions, data_weight, &points(0,i));
			f.add_variable(&assignments(0,i), 4);
			f.add_term(data_term, &assignments(0,i));
		}

		// Regularization term
		for (int i = 0; i < num_pairwise; ++i) {
			int source = connectivity(0,i);
			int target = connectivity(1,i);

			auto regularization_term =
				std::make_shared<spii::AutoDiffTerm<Regularization_cost, 4, 4>>
				(dimensions, connectivity_weights(i), tol, &points(0,source), &points(0,target), eps);

			f.add_term(regularization_term,&assignments(0,source), &assignments(0,target));
		}
	}

	double intial_function_value = f.evaluate();

	// Optimization
	std::unique_ptr<spii::Solver> solver;
	solver.reset(new spii::LBFGSSolver);
	solver->function_improvement_tolerance = 1e-8;
	solver->argument_improvement_tolerance = 1e-8;
	solver->maximum_iterations = max_iterations;

	if (verbose) {
		solver->log_function = mex_log_function;
	} 	else {
		solver->log_function = nullptr;
	}

	spii::SolverResults results;
	solver->solve(f, &results);

	double optimized_function_value = f.evaluate();

	if (verbose)
	{
		std::stringstream sout;
		sout << results << std::endl;
		f.print_timing_information(sout);
		mexPrintf("%s\n", sout.str().c_str());

		mexPrintf("Initial function value: %e\n", intial_function_value);
		mexPrintf("Final function value:   %e\n", optimized_function_value);
	}


	plhs[0] = assignments;

	return;
}