#ifndef AXIS_ANGLE_H_
#define AXIS_ANGLE_H_

#include <iostream>

#include "axis_angle.h"
#include "helpers.h"

#include <cmath>
#include <cfloat>
#include <Eigen/Dense>

////////////////////////////////////////////////////////////
/// \brief axisAngleRotationMatrix:
/// compose a rotation matrix from given rotation axis and angle,
/// make sure axis is of length 3, rotmat is of length 9.
/// \param axis: input, length 3
/// \param angle: input
/// \param rotmat: output, rotation matrix, length 9
/// \param normalized_p: whether the axis is normalized
///
void axisAngleToRotationMatrix(const double *axis, double angle, double *rotmat, bool normalized_p=true);

////////////////////////////////////////////////////////////
/// \brief rotationMatrixToAxisAngle:
/// decomposite the rotation matrix to extract the rotating axis and angle.
/// \param rotmat: input, length 9
/// \param axis: output, rotaion axis, length 3;
/// \param angle: output, rotation angle, [0, pi]
///
void rotationMatrixToAxisAngle(const double *mat, double *axis, double &angle);


////////////////////////////////////////////////////////////
/// \brief rotationMatrixToAxisAngle_simple:
/// decomposite the rotation matrix to extract the rotating axis and angle.
/// simpler approach
/// \param rotmat: input, length 9
/// \param axis: output, rotaion axis, length 3;
/// \param angle: output, rotation angle, [0, pi]
///
void rotationMatrixToAxisAngle_simple(const double *rotmat, double *axis, double &angle);


////////////////////////////////////////////////////////////
/// \brief rotationMatrixToQuaternion
/// extract Quaternion from rotation matrix
/// \param rotmat: input, rotation matrix
/// \param quaterntion, output
///
void rotationMatrixToQuaternion(const double *rotmat, double *quaterntion);

////////////////////////////////////////////////////////////
/// \brief quaternionToRotationMatrix_indirect
/// \param quaternion, input
/// \param rotmat, output, rotation matrix
///
void quaternionToRotationMatrix_indirect(const double *quaternion, double *rotmat);

////////////////////////////////////////////////////////////
/// \brief quaternionToRotationMatrix
/// compute rotation matrix directly from quaternion
/// \param quaternion
/// \param rotmat
///
void quaternionToRotationMatrix(const double *quaternion, double *rotmat);

////////////////////////////////////////////////////////////
/// \brief eulerAngleToQuaternion
/// convert euler angles to quaternion, rotatioin order: z->y->x
/// \param eulerAngles: input, euler angles, [angle_x, angle_y, angle_z]
/// \param quaternion: output
///
void eulerAngleToQuaternion_ZYX(const double *eulerAngles, double *quaternion);

////////////////////////////////////////////////////////////
/// \brief quaternionToEulerAngle
/// convert quaternion to euler angles
/// \param quaternion: input
/// \param eulerAngles: output, [angle_x, angle_y, angle_z]
///
void quaternionToEulerAngle(const double *quaternion, double *eulerAngles);

////////////////////////////////////////////////////////////
/// \brief rotationMatrixToEulerAngle_viaQuaternion
/// extract Euler rotation angles from rotation matrix, via quaternion convertion
/// \param rotmat: input, rotation matrix
/// \param eulerAngles: output, [angle_x, angle_y, angle_z];
///
void rotationMatrixToEulerAngle_viaQuaternion(const double *rotmat, double *eulerAngles);

////////////////////////////////////////////////////////////
/// \brief eulerAnglesToRotationMatrix_indirect
/// compute rotation matrix indiectly, by means of quaternion
/// \param eulerAngles: input, [angle_x, angle_y, angle_z]
/// \param rotmat: output, rotation matrix
///
void eulerAnglesToRotationMatrix_indirect(const double *eulerAngles, double *rotmat);

#endif // AXIS_ANGLE_H_
