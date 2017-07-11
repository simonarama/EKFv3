#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;


KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, VectorXd &hx_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in, MatrixXd &Hj_in) {
  x_ = x_in;
  hx_ = hx_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  Hj_ = Hj_in;


}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  
  MatrixXd Ft = F_.transpose();

  MatrixXd FP = F_ * P_;
 
 // MatrixXd P_ = ((F_ * P_) * Ft) + Q_;
  MatrixXd P_ = FP * Ft + Q_;
} 

void KalmanFilter::Update(const VectorXd &z_) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  
  VectorXd Hx = H_ * x_;
 
  VectorXd y = z_ - Hx;
   
  MatrixXd HP = H_ * P_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = HP * Ht + R_;
  //MatrixXd S = (H_ * P_ * Ht) + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  x_ = x_ + K * y;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
  
}

void KalmanFilter::UpdateEKF(const VectorXd &z_) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  double p_x = x_(0);
  double p_y = x_(1);
  double v1 = x_(2);
  double v2 = x_(3);
  double c1 = p_x*p_x + p_y*p_y;
  double c2 = sqrt(c1);
  
  hx_ = VectorXd(3);
  hx_ << 1,1,1; //initialize
  hx_ << c2, atan2(p_y,p_x), (p_x*v1 + p_y*v2)/c2;
  
  VectorXd y = z_ - hx_;
  
  //keep theta between +pi and -pi
  while (abs(y[1] > M_PI)){ 
      y[1] -= 2.* M_PI;}

  MatrixXd HjP = Hj_ * P_;
  MatrixXd Hjt = Hj_.transpose();

  MatrixXd S = HjP * Hjt + R_;
  //MatrixXd S = (Hj_ * P_ * Hjt) + R_;
  MatrixXd Si = S.inverse();

  MatrixXd PHjt = P_ * Hjt;
  MatrixXd K = PHjt * Si;
  //MatrixXd K = P_ * Hjt * Si;
  
  VectorXd Ky = K * y;
 
  x_ = Ky + x_;
  
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  MatrixXd KHj = K * Hj_;

  MatrixXd IsubKHj = I - KHj;

  P_ = IsubKHj * P_;
  
}