#pragma once


template<typename Data_cost, typename Regularization_cost>
void main_function(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction(int            nlhs,     /* number of expected outputs */
                 mxArray        *plhs[],  /* mxArray output pointer array */
                 int            nrhs,     /* number of inputs */
                 const mxArray  *prhs[]   /* mxArray input pointer array */)
{
 char problem_type[1024];
 if (mxGetString(prhs[0], problem_type, 1024)) {
   throw runtime_error("First argument must be a string.");
 }

  if (!strcmp(problem_type,"default")) {
    main_function<Quadratic_data, Default_regularization>(nlhs, plhs, nrhs, prhs);
  } else if (!strcmp(problem_type,"linear")) {
    main_function<Quadratic_data, Linear_regularization>(nlhs, plhs, nrhs, prhs);
  } else  if (!strcmp(problem_type,"quadratic")){
    main_function<Quadratic_data, Qudratic_regularization>(nlhs, plhs, nrhs, prhs);
  } else  if (!strcmp(problem_type,"length")){
    main_function<Quadratic_data, Length_regularization>(nlhs, plhs, nrhs, prhs);
  } else {
    mexErrMsgTxt("Unknown problem_type.");
   }
}