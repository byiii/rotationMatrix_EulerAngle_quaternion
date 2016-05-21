#ifndef AXIS_ANGLE_H_
#define AXIS_ANGLE_H_

#include <iostream>

////////////////////////////////////////////////////////////
/// \brief axisAngleRotationMatrix:
/// compose a rotation matrix from given rotation axis and angle,
/// make sure axis is of length 3, rotmat is of length 9.
/// \param axis: input, length 3
/// \param angle: input
/// \param rotmat: output, length 9
///
void axisAngleToRotationMatrix(const double *axis, double angle, double *rotmat);


////////////////////////////////////////////////////////////
/// \brief rotationMatrixToAxisAngle:
/// decomposite the rotation matrix to extract the rotating axis and angle.
/// \param rotmat: input, length 9
/// \param axis: output, length 3;
/// \param angle: output
///
void rotationMatrixToAxisAngle(double *rotmat, double *axis, double &angle);

#endif // AXIS_ANGLE_H_
