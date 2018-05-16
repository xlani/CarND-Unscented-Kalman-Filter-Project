#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  //set state dimension
  n_x_ = 5;

  //set augmented state dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/2;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  //set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for(int i = 1; i < 2 * n_aug_ + 1; i++) {
      weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //Time when the state is true, in us
  time_us_ = 0;

  //initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

}

UKF::~UKF() {
    output_.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if(!is_initialized_) {

      //init cout
      cout << "UKF initialized!" << endl;

      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
        float rho = meas_package.raw_measurements_[0];
        float phi = meas_package.raw_measurements_[1];
        x_ << sin(phi)*rho, cos(phi)*rho, 5, 0, 0;
      }

      else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        /**
        Initialize state.
        */
        x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 5, 0, 0;
      }

        //init covariance matrix
        P_ << 0.15, 0, 0, 0, 0,
              0, 0.15, 0, 0, 0,
              0, 0,  0.3, 0, 0,
              0, 0, 0, 0.03, 0,
              0, 0, 0, 0, 0.03;

      //Save timestamp of first measurement
      time_us_ = meas_package.timestamp_;

      //Done initializing, no need to predict or update
      is_initialized_ = true;

      //files to store nis parameters
      output_.open("./data.csv");
      output_ << "time, NIS_laser_, threshold_laser, NIS_radar_, threshold_radar" << endl;

      return;
  }

  /*****************************************************************************
  *  Skip measurement if sensor is disabled
  ****************************************************************************/

  if(!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      cout << endl << "Radar measurement skipped!" << endl;
      return;
  }
  else if(!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
      cout << endl << "Lidar measurement skipped!" << endl;
      return;
  }

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/

  /*** Update the state transition matrix F according to the new elapsed time.
       - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //compute the time elapsed between the current and previous measurements

  dt_ = (meas_package.timestamp_ - time_us_) / 1000000.0;	//dt - expressed in seconds
  time_us_ = meas_package.timestamp_;

  Prediction(dt_);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /*** Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    //radar updates
    UpdateRadar(meas_package);
  }

  else {
    //laser updates
    UpdateLidar(meas_package);
  }

  //print the output
  //cout << "time_us_ = " << time_us_ << endl;
  //cout << "dt = " << dt_ << endl;
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.setZero();
  x_aug.head(n_x_) = x_;

  //create augmented covariance matrix
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  //create augmented sigma points
  //set sigma points as columns of matrix Xsig
  Xsig_aug.col(0) = x_aug;

  //set remaining sigma points
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

      //extract values out of matrix
      double px       = Xsig_aug(0,i);
      double py       = Xsig_aug(1,i);
      double v        = Xsig_aug(2,i);
      double yaw      = Xsig_aug(3,i);
      double yawd     = Xsig_aug(4,i);
      double nu_a     = Xsig_aug(5,i);
      double nu_yawdd = Xsig_aug(6,i);

      //predicted state for positions
      double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = px + v/yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
        py_p = py + v/yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
    }
    else {
        px_p = px + v*delta_t*cos(yaw);
        py_p = py + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //predict state mean
  x_.setZero();
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {
      x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.setZero();
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {

      VectorXd diff = Xsig_pred_.col(i) - x_;

      //angle normalization
      while (diff(3)> M_PI) diff(3)-=2.*M_PI;
      while (diff(3)<-M_PI) diff(3)+=2.*M_PI;

      P_ = P_ + weights_(i) * diff * diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  //set measurement dimension, lidar can measure px and py
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //create vector for incoming lidar measurement (px, py)
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],   //px in m
       meas_package.raw_measurements_[1];   //py in m

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //transform sigma points into measurement space
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {

      //set values into right column
      Zsig(0,i) = Xsig_pred_(0,i);
      Zsig(1,i) = Xsig_pred_(1,i);

  }
  //calculate mean predicted measurement
  z_pred.setZero();
  for(int i = 0; i< 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //calculate innovation covariance matrix S
  S.setZero();
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {

      VectorXd diff = Zsig.col(i) - z_pred;
      S = S + weights_(i) * diff * diff.transpose();

  }

  //add noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;

  S = S + R;

  //calculate cross correlation matrix
  Tc.setZero();
  for(int i = 0; i< 2 * n_aug_ + 1; i++) {

      //vectors for residuals
      VectorXd diffx = Xsig_pred_.col(i) - x_;
      VectorXd diffz = Zsig.col(i) - z_pred;

      //cross correlation matrix
      Tc = Tc + weights_(i) * diffx * diffz.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc * Si;

  //residual vector
  VectorXd diffz = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * diffz;
  P_ = P_ - K * S * K.transpose();

  //update NIS for laser measurement
  NIS_laser_ = diffz.transpose()*Si*diffz;
  output_ << time_us_ << "," << NIS_laser_ << ",5.991,0.,7.815" << endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);

  //create vector for incoming radar measurement (rho, phi, rhod)
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0],   //rho in m
       meas_package.raw_measurements_[1],   //phi in rad
       meas_package.raw_measurements_[2];   //rho_dot in m/s

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //transform sigma points into measurement space
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {

      //extract values for better readability
      double px = Xsig_pred_(0,i);
      double py = Xsig_pred_(1,i);
      double v = Xsig_pred_(2,i);
      double yaw = Xsig_pred_(3,i);

      //calculate rho, phi and rhod
      double rho = sqrt(px*px + py*py);
      double phi = atan2(py,px);
      double rhod;
      if(rho > 0.001) {
          rhod = v * (px*cos(yaw) + py*sin(yaw)) / rho;
      }
      else {
          cout << "Division by zero avoided!" << endl;
          rhod = 0;
      }

      //set values into right column
      Zsig(0,i) = rho;
      Zsig(1,i) = phi;
      Zsig(2,i) = rhod;

  }
  //calculate mean predicted measurement
  z_pred.setZero();
  for(int i = 0; i< 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //calculate innovation covariance matrix S
  S.setZero();
  for(int i = 0; i < 2 * n_aug_ + 1; i++) {

      VectorXd diff = Zsig.col(i) - z_pred;

      //angle normalization
      while (diff(1)> M_PI) diff(1)-=2.*M_PI;
      while (diff(1)<-M_PI) diff(1)+=2.*M_PI;

      S = S + weights_(i) * diff * diff.transpose();

  }

  //add noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0, std_radrd_*std_radrd_;

  S = S + R;

  //calculate cross correlation matrix
  Tc.setZero();
  for(int i = 0; i< 2 * n_aug_ + 1; i++) {

      //vectors for residuals
      VectorXd diffx = Xsig_pred_.col(i) - x_;
      VectorXd diffz = Zsig.col(i) - z_pred;

      //angle normalization
      while(diffx(3)> M_PI) diffx(3) -= 2.*M_PI;
      while(diffx(3)<-M_PI) diffx(3) += 2.*M_PI;
      while(diffz(1)> M_PI) diffz(1) -= 2.*M_PI;
      while(diffz(1)<-M_PI) diffz(1) += 2.*M_PI;

      //cross correlation matrix
      Tc = Tc + weights_(i) * diffx * diffz.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd Si = S.inverse();
  MatrixXd K = Tc * Si;

  //residual vector and angle normalization
  VectorXd diffz = z - z_pred;
  while(diffz(1)> M_PI) diffz(1) -= 2.*M_PI;
  while(diffz(1)<-M_PI) diffz(1) += 2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * diffz;
  P_ = P_ - K * S * K.transpose();

  //update NIS for laser measurement
  NIS_radar_ = diffz.transpose()*Si*diffz;
  // cout << "NIS radar = " << NIS_radar_ << endl;
  // cout << "diffz = " << diffz << endl;
  output_ << time_us_ << ",0.,5.991," << NIS_radar_ << ",7.815" << endl;

}
