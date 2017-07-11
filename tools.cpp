#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO://code from lecture
    * Calculate the RMSE here. calculate residual and sum in loop
  */
  VectorXd RMSE(4);
  RMSE << 0,0,0,0;
  
  for(unsigned int i=0; i<estimations.size(); ++i){
      VectorXd residual = estimations[i] - ground_truth[i];
      residual = residual.array()*residual.array();
      RMSE+=residual;
  
}
//mean of residual sum
RMSE=RMSE/estimations.size();

//square root of mean
RMSE=RMSE.array().sqrt();

//return rmse
return RMSE;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_) {
  /**
  TODO: code from lecture
    * Calculate a Jacobian here. 1st recover state parameters
  */
  MatrixXd Hj(3,4);
  double p_x = x_(0);
  double p_y = x_(1);
  double v1 = x_(2);
  double v2 = x_(3);
  
  //simplify calcs, c3 follows exponent addition rule
  double c1 = p_x*p_x + p_y*p_y;
  double c2 = sqrt(c1);
  double c3 = (c1*c2);
  
  //division by 0 check
  if(fabs(c1)<0.0001){
      cout<<"Jacobian Division by Zero Error"<<endl;
      return Hj;
  }
  //Jacobian calc
  Hj<<(p_x/c2),(p_y/c2),0,0,
      -(p_y/c1),(p_x/c1),0,0,
      p_y*(v1*p_y-v2*p_x)/c3, p_x*(p_x*v2-p_y*v1)/c3, p_x/c2, p_y/c2;
      return Hj;
}
