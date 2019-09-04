#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  if(estimations.size() != 0)
  {
  	if(estimations.size() == ground_truth.size())
  	{
        VectorXd sum(4);
        sum << 0, 0, 0, 0;
    	for(int count = 0; count < estimations.size(); count++)
    	{
            VectorXd diff = (estimations[count] - ground_truth[count]);
    		VectorXd mul = (diff.array() * diff.array());
            sum += mul;
     	 }
         rmse = (sum.array())/estimations.size();
         rmse = rmse.array().sqrt();
  	}
    else
    {
      std::cout << "RMSE: estimations size and gt size not equal " << std::endl;
    }
  }
  else
  {
    std::cout << "RMSE: estimations size is zero " << std::endl;
  }
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3, 4);
  Hj << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  double p_sqr = (px*px + py*py);
  if(0 < p_sqr)
  {
    double p_sqroot = sqrt(p_sqr);
    double p_sq_3_by_2 = p_sqr * p_sqroot;
    Hj << px/p_sqroot, py/p_sqroot, 0, 0, 
          -1 * py/p_sqr, px/p_sqr, 0, 0,
          py*(vx*py - vy*px)/p_sq_3_by_2, px*(vy*px - vx*py)/p_sq_3_by_2, px/p_sqroot, py/p_sqroot;
  }
  else
  {
    std::cout << "Jacobian: Position square is zero " << std::endl;
  }
  
  return Hj;
  
/*   MatrixXd Hj(3, 4);
  
  const double APPROX_ZERO = 0.0001;

    // Unpack the state vector
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    // Calculate frequently used calculations
    double px2 = px * px;
    double py2 = py * py;

    // If px2 is zero set it to a small value
    if (fabs(px2) < APPROX_ZERO) {
      px2 = APPROX_ZERO;
    }

    // If py2 is zero set it to a small value
    if (fabs(py2) < APPROX_ZERO) {
      py2 = APPROX_ZERO;
    }

    double ss = px2 + py2;
    double srss = sqrt(ss);

    // Create the Jacobian
    Hj(0, 0) = px / srss;
    Hj(0, 1) = py / srss;
    Hj(0, 2) = 0;
    Hj(0, 3) = 0;

    Hj(1, 0) = -py / ss;
    Hj(1, 1) = px / ss;
    Hj(1, 2) = 0;
    Hj(1, 3) = 0;

    Hj(2, 0) = (py * (vx * py - px * vy)) / (ss * srss);
    Hj(2, 1) = (px * (vy * px - py * vx)) / (ss * srss);
    Hj(2, 2) = px / srss;
    Hj(2, 3) = py / srss;

return Hj;*/
}
