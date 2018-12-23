#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    if (estimations.size() == 0) {
        cout << "size of estimations is zero";
        return rmse;
    }

    if (estimations.size() != ground_truth.size()) {
        cout << "size of estimations is not equal to size of ground truth";
        return rmse;
    }

    // accumulate squared residuals
    for (int i=0; i < estimations.size(); ++i) {
        VectorXd diff = estimations.at(i)  - ground_truth.at(i);
        rmse = rmse.array() + diff.array() * diff.array() ;
    }

    // calculate the mean
    rmse /= estimations.size();

    // calculate the squared root
    rmse = rmse.array().sqrt();

    // return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
    MatrixXd Hj(3,4);
    Hj << 0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0;
    // recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    // check division by zero
    float c1 = px*px+py*py;
    if (fabs(c1) < 0.0001) {
        cout << "both px and py are zero";
        return Hj;
    }
    // compute the Jacobian matrix
    float px2_py2 = px*px + py*py;
    float px2_py2_sqrt = sqrt(px2_py2);
    float px2_py2_3_2 = pow(px2_py2, 3.0/2.0);

    Hj << px / px2_py2_sqrt, py / px2_py2_sqrt, 0, 0,
    - py / px2_py2, px / px2_py2, 0, 0,
    py * (vx * py - vy * px) / px2_py2_3_2, px *(vy*px - vx*py) / px2_py2_3_2, px / px2_py2_sqrt, py / px2_py2_sqrt ;


    return Hj;
}
