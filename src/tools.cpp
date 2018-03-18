#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /*
   * Calculating the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  bool validData = false;

  if ((estimations.size()!=0) && (estimations.size() == ground_truth.size()))
  {
      validData = true;
  }
  else
  {
      cout<<"CalculateRMSE(): Error - Incorrect vector length(s) \n";
  }

  if (validData == true)
  {
      VectorXd residual(4);
      VectorXd residualSum(4);
      residualSum << 0,0,0,0;
      /*
       * accumulate squared residuals
       */
      for(int i=0; i < estimations.size(); ++i){
            residual = (estimations[i] - ground_truth[i]);
            residual = residual.array().square();
            residualSum += residual;
      }
      /*
       * calculating mean of the residuals
       */
      residualSum = residualSum/estimations.size();

      /*
       * calculate the squared root
       */
      rmse = residualSum.array().sqrt();
  }

  return rmse;
}
