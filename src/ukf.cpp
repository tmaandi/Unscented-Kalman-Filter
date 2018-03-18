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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

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

  is_initialized_ = false;

  // Number of states
  n_x_ = 5;

  // Number of states after process noise augumentation
  n_aug_ = 7;

  time_us_ = 0;

  lambda_ = 3 - n_x_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd( 2 * n_aug_ + 1);


  ///* Lidar measurement noise covariance matrix
  R_Lidar_ = MatrixXd(2,2);
  R_Lidar_ << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;

  ///* Radar measurement noise covariance matrix
  R_Radar_ = MatrixXd(3,3);
  R_Radar_ << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (is_initialized_ == false)
  {
    /*****************
     * Set weights
     *****************/
    for (int i = 0; i < (2 * n_aug_ + 1); ++i)
    {
      if (i == 0)
      {
        weights_(i) = lambda_/(lambda_ + n_aug_);
      }
      else
      {
        weights_(i) = 1/(2 * (lambda_ + n_aug_));
      }
    }

    /*
     * Initialize timestamp
     */
    time_us_ = meas_package.timestamp_;

    /***************************************************
     * Initializing states and process covariance matrix
     ***************************************************/

    if (meas_package.sensor_type_== MeasurementPackage::LASER)
    {
      /*
       * Initializing state vector using first measurement readings
       */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.0, 0.0, 0.0;
      /*
       * Initializing process covariance matrix, assuming variance in the x and y
       * position states 'similar' to the ones in the sensor readings, and 1 for
       * the rest
       */
      P_ = MatrixXd::Identity(n_x_, n_x_);
      P_(0,0) = std_laspx_ * std_laspx_;
      P_(1,1) = std_laspy_ * std_laspy_;
    }
    else
    {
      /*
       * Initializing if the first measurement is from Radar
       */
      float rho = meas_package.raw_measurements_[0];
      float theta = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];

      /*
       * Convert radar from polar to Cartesian coordinates and initialize state.
       */
      float px = rho*cos(theta);
      float py = rho*sin(theta);
      float sigma_px = std_radr_ * cos(std_radphi_);
      float sigma_py = std_radr_ * sin(std_radphi_);

      x_ << px, py, rho_dot/3, rho_dot/3, 0.0;

      /*
       * Initializing process covariance matrix
       */
      P_ = MatrixXd::Identity(n_x_, n_x_);
      /*
       * std.dev of px, py < = std. dev of rho
       */
      P_(0,0) = sigma_px * sigma_px;
      P_(1,1) = sigma_py * sigma_py;
    }

    /*
     * done initializing, no need to predict or update
     */
    is_initialized_ = true;
    return;
  }

  /*
   * Calculating delta time between two measurements
   */
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  /*
   * Prediction Step
   */
  Prediction(dt);

  /*
   * Update Step
   */
  if ((meas_package.sensor_type_== MeasurementPackage::LASER)\
      && (use_laser_ == true))
  {
    UpdateLidar(meas_package);
  }
  else if ((meas_package.sensor_type_== MeasurementPackage::RADAR)\
      && (use_radar_ == true))
  {
    UpdateRadar(meas_package);
  }
  else
  {
    //cout<<"Unknown Sensor"<<endl;
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  double delta_t2 = delta_t * delta_t;

  /****************************************
   * Creating augmented sigma points
   ****************************************/

  /*
   * Create augmented mean vector
   */
  VectorXd x_aug = VectorXd(n_aug_);

  /*
   * Create augmented state covariance
   */
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  /*
   * Create sigma point matrix
   */
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  /*
   * Create augmented mean state
   */
  x_aug.setZero();
  x_aug.head(n_x_) = x_;
  x_aug.tail(n_aug_ - n_x_) << 0,0;

  /*
   * Create augmented covariance matrix
   */
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_ - n_x_, n_aug_ - n_x_) << std_a_ * std_a_, 0, 0, std_yawdd_ * std_yawdd_;

  /*
   * Create square root matrix
   */
  MatrixXd A_aug = MatrixXd(n_aug_,n_aug_);
  A_aug = P_aug.llt().matrixL();

  /*
   * Create augmented sigma points
   */
  Xsig_aug.col(0) = x_aug;

  for (int col = 0; col < n_aug_; ++col)
  {
    Xsig_aug.col(col + 1) = x_aug + sqrt(lambda_ + n_aug_) * A_aug.col(col);
    Xsig_aug.col(col + n_aug_ + 1) = x_aug - sqrt(lambda_ + n_aug_) * A_aug.col(col);
  }

  /*********************************************
   * Create matrix with predicted sigma points
   *********************************************/

  for (int col = 0; col < (2 * n_aug_ + 1); ++col)
  {
    double p_x, p_y, v_k, psi_k, psi_dot_k, nu_a_k, nu_psi_ddot_k;

    p_x = Xsig_aug(0, col);
    p_y = Xsig_aug(1, col);
    v_k = Xsig_aug(2, col);
    psi_k = Xsig_aug(3, col);
    psi_dot_k = Xsig_aug(4, col);
    nu_a_k = Xsig_aug(5, col);
    nu_psi_ddot_k =Xsig_aug(6, col);

    if (psi_dot_k != 0)
    {
      Xsig_pred_(0, col) = p_x + (v_k / psi_dot_k) * \
          (sin(psi_k + psi_dot_k * delta_t) - sin(psi_k)) + \
          0.5 * (delta_t2) * cos(psi_k) * nu_a_k;

      Xsig_pred_(1, col) = p_y + (v_k / psi_dot_k) * \
          (-cos(psi_k + psi_dot_k * delta_t) + cos(psi_k)) + \
          0.5 * delta_t2 * sin(psi_k) * nu_a_k;

      Xsig_pred_(2, col) = v_k + delta_t * nu_a_k;

      Xsig_pred_(3, col) = psi_k + psi_dot_k * delta_t + \
          0.5 * delta_t2 * nu_psi_ddot_k;

      Xsig_pred_(4, col) = psi_dot_k + delta_t * nu_psi_ddot_k;
    }
    else
    {
      Xsig_pred_(0, col) = p_x + v_k * cos(psi_k) * delta_t +\
         0.5 * delta_t2 * cos(psi_k) * nu_a_k;

      Xsig_pred_(1, col) = p_y + v_k * sin(psi_k) * delta_t + \
         0.5 * delta_t2 * sin(psi_k) * nu_a_k;

      Xsig_pred_(2, col) = v_k + delta_t * nu_a_k;

      Xsig_pred_(3, col) = psi_k + psi_dot_k * delta_t + \
          0.5 * delta_t2 * nu_psi_ddot_k;

      Xsig_pred_(4, col) = psi_dot_k + delta_t * nu_psi_ddot_k;
    }
  }

  /**********************
   * Predict state mean
   **********************/
  /*
   * Zeroing x
   */
  x_.setZero();

  for (int i = 0; i < (2 * n_aug_ + 1); ++i)
  {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  /*************************************
   * Predicting state covariance matrix
   *************************************/
  /*
   * Zeroing P
   */
  P_.setZero();

  for (int i = 0; i < (2 * n_aug_ + 1); ++i)
  {
    /*
     * State difference
     */
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    /*
     * Angle normalization
     */
    while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) <-M_PI) x_diff(3) += 2. * M_PI;
    P_ += weights_(i) * ((x_diff) * (x_diff.transpose()));
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  int n_z = 2;

  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;

  MatrixXd H = MatrixXd(n_z, n_x_);

  H.setZero();

  H << 1,0,0,0,0,
       0,1,0,0,0;

  /*
   * update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H * x_;
  VectorXd z_diff = z - z_pred;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ *Ht + R_Lidar_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  /*
   * new estimate
   */
  x_ = x_ + K * z_diff;
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;

  /*
   * Calcualte Lidar NIS
   * It checks for the consistency of the filter.
   * If too high, it means we underestimated the
   * variance in the process noise, if too low,
   * it means we overestimated the variance in
   * the process noise
   */

  double nisLidar;

  nisLidar = (z_diff.transpose()) * S.inverse() * z_diff;

//  cout<<"NIS = "<<nisLidar<<endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  /*
   * Length of radar measurement vector
   */
  int n_z = 3;
  /*
   * Create matrix for sigma points in measurement space
   */

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  /*
   * Mean predicted measurement
   */
  VectorXd z_pred = VectorXd(n_z);

  /*
   * Measurement covariance matrix S
   */
  MatrixXd S = MatrixXd(n_z, n_z);

  /*
   * Zeroing Zsigs
   */
  Zsig.setZero();

  /*
   * Transform sigma points into measurement space
   */

  for (int col = 0; col < 2 * n_aug_ + 1; ++col)
  {
    double px, py, v, psi, psi_dot, rho, phi, rho_dot;
    px = Xsig_pred_(0, col);
    py = Xsig_pred_(1, col);
    v = Xsig_pred_(2, col);
    psi = Xsig_pred_(3, col);
    psi_dot = Xsig_pred_(4, col);

    rho = sqrt(px * px + py * py);
    phi = atan2(py, px);

    /*
     * Zero division protection
     */
    if(rho != 0)
    {
      rho_dot = v * (px * cos(psi) + py * sin(psi))/rho;
    }
    else
    {
      rho_dot = 0;
    }

    Zsig.col(col) << rho, phi, rho_dot;
  }

  /*
   * Calculate mean predicted measurement
   */

  /*
   * Initialize z_pred
   */
  z_pred.setZero();

  for (int col = 0; col < (2 * n_aug_ + 1); ++col)
  {
    z_pred += weights_(col) * Zsig.col(col);
  }

  /*
   * Calculate innovation covariance matrix S
   */
  /*
   * Initialize S
   */
  S.setZero();

  /*
   * Difference between measurement and estimation
   */
  VectorXd z_diff = VectorXd(n_z);

  /*
   * Zeroing meas. diff
   */
  z_diff.setZero();

  for(int col = 0; col < (2 * n_aug_ + 1); ++col)
  {
    z_diff = Zsig.col(col) - z_pred;
    /*
     * Angle normalization
     */
    while (z_diff(1) >  M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;
    S += weights_(col) * z_diff * z_diff.transpose();
  }

  /*
   * Adding measurement noise covariance matrix
   */
  S += R_Radar_;

  /*
   * Create vector for incoming measurement
   */
  VectorXd z = VectorXd(n_z);
  z = meas_package.raw_measurements_;

  /*
   * Create matrix for cross correlation Tc
   */
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  /*
   * Zeroing Tc
   */
  Tc.setZero();

  /*
   * Calculate cross correlation matrix
   */
  for (int col = 0; col < (2 * n_aug_ + 1); ++col)
  {
    /*
     * Measurement Residual
     */
    VectorXd z_diff = Zsig.col(col) - z_pred;

    /*
     * Angle normalization
     */
    while(z_diff(1) >  M_PI) z_diff(1) -= 2. * M_PI;
    while(z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    /*
     * State Residual
     */
    VectorXd x_diff = Xsig_pred_.col(col) - x_;
    /*
     * Angle normalization
     */
    while (x_diff(3) >  M_PI) x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

    Tc += (weights_(col))*(x_diff)*((z_diff).transpose());

  }
  /*
   * Calculate Kalman Gain
   */
  MatrixXd K  = MatrixXd(n_x_, n_z);

  K = Tc * S.inverse();

  /*
   * Update state mean and covariance matrix
   */
  x_ = x_ + K * (z - z_pred);

  P_ = P_ - K * S * K.transpose();

  /*
   * Calcualte Lidar NIS
   * It checks for the consistency of the filter.
   * If too high, it means we underestimated the
   * variance in the process noise, if too low,
   * it means we overestimated the variance in
   * the process noise
   */

  double nis;

  z_diff = z - z_pred;

  nis = (z_diff.transpose()) * S.inverse() * z_diff;

//  cout<<"NIS = "<<nis<<endl;

}








