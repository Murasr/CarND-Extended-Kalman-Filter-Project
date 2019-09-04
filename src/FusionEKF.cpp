#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  // measurement matrix - laser
  H_laser_ << 1, 0, 0, 0, 
              0, 1, 0, 0;

  // Initialize the state covariance matrix
  // Provide higher covariance for velocity
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {
 
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Initialize the first measurement
      // Convert from polar to cartesian coordinates
      // Assign zero for velocity
      ekf_.x_ << measurement_pack.raw_measurements_(0) * cos(measurement_pack.raw_measurements_(1)), measurement_pack.raw_measurements_(0)* sin(measurement_pack.raw_measurements_(1)), measurement_pack.raw_measurements_(2) * cos(measurement_pack.raw_measurements_(1)), measurement_pack.raw_measurements_(2)* sin(measurement_pack.raw_measurements_(1));
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize the first measurement
      // Assign zero for velocity
      ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
    }
    
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  // Do the normalization of time difference
  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
  //std::cout << "time diff: " << dt << std::endl;
  /**
   * Prediction
   */
  
  // State matrix calculation
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  
  // Set the process covariance matrix
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;
  
  double sig_ax_2 = 9;
  double sig_ay_2 = 9;
  
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4/4 * sig_ax_2, 0, dt_3/2 * sig_ax_2, 0,
             0, dt_4/4 * sig_ay_2, 0, dt_3/2 * sig_ay_2,
             dt_3/2 * sig_ax_2, 0, dt_2 * sig_ax_2, 0,
             0, dt_3/2 * sig_ay_2, 0, dt_2 * sig_ay_2;

  // Do prediction
  ekf_.Predict();

  /**
   * Update
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {
    // Calculate jacobian for radar to linearize the non linear measurements
    Tools toolInst;
    Hj_ = toolInst.CalculateJacobian(ekf_.x_);
    
    // Copy the H and R matrix of radar
	ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    
    // Do Update EKF method
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    //  Copy the H and R matrix of laser
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;

    // Do measurement update
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
 cout << "x_ = " << ekf_.x_ << endl;
 cout << "P_ = " << ekf_.P_ << endl;
}
