#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"
#include "measurement_package.h"
class KalmanFilter {
public:

  // state vector
  Eigen::VectorXd x_;
  
  //add this state vector
  Eigen::VectorXd hx_;
  
  // state covariance matrix
  Eigen::MatrixXd P_;

  // state transistion matrix
  Eigen::MatrixXd F_;
  
  // measurement matrix
  Eigen::MatrixXd H_;

  // measurement covariance matrix
  Eigen::MatrixXd R_;

  // process covariance matrix
  Eigen::MatrixXd Q_;

  // add Hj_
  Eigen::MatrixXd Hj_;
  

  
  
  

  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Init Initializes Kalman filter
   * @param x_in Initial state
   * @param P_in Initial state covariance
   * @param F_in Transition matrix
   * @param H_in Measurement matrix
   * @param R_in Measurement covariance matrix
   * @param Q_in Process covariance matrix
   * add state vector hx_ in void
   * add matrix Hj_ in void
   */
  void Init(Eigen::VectorXd &x_in, Eigen::VectorXd &hx_in, Eigen::MatrixXd &P_in, Eigen::MatrixXd &F_in,
      Eigen::MatrixXd &H_in, Eigen::MatrixXd &R_in, Eigen::MatrixXd &Q_in, Eigen::MatrixXd &Hj_in);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param delta_T Time between k and k+1 in s
   */
  void Predict();

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z_);


  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z_);

  
  
  
  /**
  * Run the whole flow of the Kalman Filter from here.
  */

  /**
  * Kalman Filter update and prediction math lives in here.
  */

  
};

#endif /* KALMAN_FILTER_H_ */
