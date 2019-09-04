#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;


/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() 
{
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  
 // std::cout << "After prediction " << x_ << std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
 
  VectorXd y = z - H_ * x_;
  
  //std::cout << "Radar y: " << y <<std::endl;
    
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;
  
  x_ = x_ + K * y;
  
  MatrixXd I = MatrixXd::Identity(4, 4);
  P_ = (I - K * H_) * P_;   
  
   //std::cout << "Radar z: " << z <<std::endl<< "After update " << x_ << std::endl;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd h_of_x = VectorXd(3);
  h_of_x(0) = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
  
  // Denominator should be greater than zero.
  if(h_of_x(0) > 0)
  {
    h_of_x(1) = atan2(x_(1), x_(0));
    h_of_x(2) = (x_(0) * x_(2) + x_(1) * x_(3))/h_of_x(0);

    // Normalize the calculated angle to be within -phi and +phi
    if(h_of_x(1) < -M_PI)
      h_of_x(1) += 2*M_PI;

    if(h_of_x(1) > M_PI)
      h_of_x(1) -= 2*M_PI;
    
   // h_of_x = cartesian_to_polar(x_);
   // std::cout << "After cartesian to polar: " << h_of_x << std::endl;

    VectorXd y = z - h_of_x;
    
    // Normalize the calculated angle to be within -phi and +phi
    if(y(1) < -M_PI)
      y(1) += 2*M_PI;

    if(y(1) > M_PI)
      y(1) -= 2*M_PI;
    
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * Ht * Si;

    x_ = x_ + K * y;

    MatrixXd I = MatrixXd::Identity(4, 4);
    P_ = (I - K * H_) * P_;  
  }
}
