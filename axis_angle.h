#ifndef AXIS_ANGLE_H_
#define AXIS_ANGLE_H_

#include <iostream>

#include <Eigen/Dense>

//------------------------------------------------------------
// use Eigen library
//------------------------------------------------------------


//------------------------------------------------------------
void axisAngleToQuaternion(const Eigen::Vector3f &axis,
                           float angle,
                           Eigen::Quaternionf &quaternion);

//------------------------------------------------------------
void quaternionToAxisAngle(const Eigen::Quaternionf &quaternion,
                           Eigen::Vector3f &axis,
                           float &angle);

//------------------------------------------------------------
void axisAngleToRotationMatrix(const Eigen::Vector3f &axis,
                               float angle,
                               Eigen::Matrix3f &rotmat);

//------------------------------------------------------------
void rotationMatrixToAxisAngle(const Eigen::Matrix3f &rotmat,
                               Eigen::Vector3f &axis,
                               float &angle);

//------------------------------------------------------------
void rotationMatrixToQuaternion(const Eigen::Matrix3f &rotmat, Eigen::Quaternionf &quaternion);

//------------------------------------------------------------
void quaternionToRotationMatrix(const Eigen::Quaternionf &quaternion, Eigen::Matrix3f &rotmat);

//------------------------------------------------------------
void eulerAngleToQuaternion_XYZ(const Eigen::Vector3f &eulerAngles, Eigen::Quaternionf &quaternion);

//------------------------------------------------------------
void quaternionToEulerAngle_XYZ(const Eigen::Quaternionf &quaternion, Eigen::Vector3f &eulerAngles);

//------------------------------------------------------------
void rotationMatrixToEulerAngle_XYZ(const Eigen::Matrix3f &rotmat, Eigen::Vector3f &eulerAngles);

//------------------------------------------------------------
void eulerAnglesToRotationMatrix_XYZ(const Eigen::Vector3f &eulerAngles, Eigen::Matrix3f &rotmat);




//------------------------------------------------------------
// double array representation
//------------------------------------------------------------


//------------------------------------------------------------
void axisAngleToRotationMatrix(const double *axis,
                               double angle,
                               double *rotmat,
                               bool normalized_p = true);

//------------------------------------------------------------
void rotationMatrixToAxisAngle(const double *mat,
                               double *axis,
                               double &angle);

//------------------------------------------------------------
void rotationMatrixToAxisAngle_simple(const double *rotmat,
                                      double *axis,
                                      double &angle);

//------------------------------------------------------------
void axisAngleToQuaternion(const double *axis,
                           double angle,
                           double *quaternion);

//------------------------------------------------------------
void quaternionToAxisAngle(const double *quaternion,
                           double *axis,
                           double &angle);

//------------------------------------------------------------
void rotationMatrixToQuaternion(const double *rotmat, double *quaterntion);

//------------------------------------------------------------
void quaternionToRotationMatrix(const double *quaternion, double *rotmat);

//------------------------------------------------------------
void eulerAngleToQuaternion_XYZ(const double *eulerAngles, double *quaternion);

//------------------------------------------------------------
void quaternionToEulerAngle_XYZ(const double *quaternion, double *eulerAngles);

//------------------------------------------------------------
void rotationMatrixToEulerAngle_XYZ(const double *rotmat, double *eulerAngles);

//------------------------------------------------------------
void eulerAnglesToRotationMatrix_XYZ(const double *eulerAngles, double *rotmat);




#endif // AXIS_ANGLE_H_
