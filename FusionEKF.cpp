#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;




/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices and vectors
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  ekf_.H_ = MatrixXd(2, 4);

  

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

//initialize H_laser_
  ekf_.H_ << 1,0,0,0,
             0,1,0,0;



//initialize P_
  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1,0,0,0,
             0,1,0,0,
             0,0,1000,0,
             0,0,0,1000;




  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &meas_package) {

//was &measurement_pack
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    
    
    
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 0, 0;
    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /** sign? on py
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //set the state with the initial location and zero velocity

      float ro = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float ro_dot = meas_package.raw_measurements_[2];
      float px = ro*cos(theta);
      float py = ro*sin(theta);
      double p_x = px;
      double p_y = py;
      ekf_.x_<< p_x, p_y, 0, 0;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      float px = meas_package.raw_measurements_[0];

      float py = meas_package.raw_measurements_[1];

      double p_x = px;
      double p_y = py;
      ekf_.x_<< p_x, p_y, 0, 0;

      }
      previous_timestamp_ = meas_package.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
*/
// Radar and laser updates
  float dt = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;

  float dt2 = dt * dt;
  float dt3 = dt2 * dt;
  float dt4 = dt3 * dt;
  previous_timestamp_ = meas_package.timestamp_;
  
    //update state matrix 
  ekf_.F_ = MatrixXd(4,4);
    
  ekf_.F_ << 1,0,dt,0,
             0,1,0,dt,
             0,0,1,0,
             0,0,0,1;

    //update covariance matrix
  float noise_ax = 9.0;
  float noise_ay = 9.0;

  ekf_.Q_ = MatrixXd(4,4);
  ekf_.Q_ << (dt4/4)* noise_ax, 0, (dt3/2)*noise_ax, 0,
         0, (dt4/4)*noise_ay, 0, (dt3/2)*noise_ay,
         (dt3/2)*noise_ax, 0, dt2*noise_ax, 0,
         0, (dt3/2)*noise_ay, 0, dt2*noise_ay;

  
  ekf_.Predict();
  
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */


  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    
    ekf_.Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(meas_package.raw_measurements_);
  } 
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.Update(meas_package.raw_measurements_);
  }
  
}