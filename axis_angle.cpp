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
void axisAngleToRotationMatrix(const double *axis, double angle, double *rotmat, bool normalized_p)
{
    double ux = axis[0];
    double uy = axis[1];
    double uz = axis[2];

    // check if the axis is an unit vector, if not, convert to unit vector
    if(!normalized_p)
    {
        double axis_length = std::sqrt(ux*ux+uy*uy+uz*uz);
        ux = ux/axis_length;
        uy = uy/axis_length;
        uz = uz/axis_length;
    }

    // use c and s to represent cos angle and sin angle
    double c = cos(angle);
    double s = sin(angle);

    // compute the rotation matrix element
    rotmat[0] = ux*ux*(1-c)+c;    rotmat[1] = ux*uy*(1-c)-uz*s; rotmat[2] = ux*uz*(1-c)+uy*s;
    rotmat[3] = ux*uy*(1-c)+uz*s; rotmat[4] = uy*uy*(1-c)+c;    rotmat[5] = uy*uz*(1-c)-ux*s;
    rotmat[6] = ux*uz*(1-c)-uy*s; rotmat[7] = uy*uz*(1-c)+ux*s; rotmat[8] = uz*uz*(1-c)+c;

    //    std::cout << "\n------------------------------------------------------------"
    //              << std::endl;
    //    printMatrix(rotmat, 3, 3);
    //    std::cout << "\n------------------------------------------------------------"
    //              << std::endl;
}


////////////////////////////////////////////////////////////
/// \brief rotationMatrixToAxisAngle:
/// decomposite the rotation matrix to extract the rotating axis and angle.
/// \param rotmat: input, length 9
/// \param axis: output, rotaion axis, length 3;
/// \param angle: output, rotation angle, [0, pi]
///
void rotationMatrixToAxisAngle(const double *mat, double *axis, double &angle)
{
    // copy input mat to rotmat, because i do not know how to use [Eigen::Map]
    // with [const double *]
    double rotmat[9] = {0.0};
    for(size_t i=0; i<9; ++i)
        rotmat[i] = mat[i];

    // compute the matrix trace to get rotation angle.
    // here rotation angle is in the interval [0, pi],
    // so attention should paid to the rotation angle value outside this
    // interval.
    double trace = rotmat[0]+rotmat[4]+rotmat[8];
    double cos_theta = 0.5*(trace-1);
    if(fabs(cos_theta+1)<2*DBL_EPSILON)
        angle = M_PI;
    else
        angle = acos(cos_theta);

    // convert [double *] to [Eigen::Matrix3d],
    // take care that [Eigen::Matrix] is default to use the column major, so here
    // use the transpose.
    Eigen::Matrix3d rot = Eigen::Map<Eigen::MatrixXd>(rotmat, 3, 3).transpose();

    // compute [R-I], to get the eigen vector corresponding to the eigen velue 1
    Eigen::Matrix3d rot_i = rot - Eigen::Matrix3d::Identity();
    // do svd, the ideal eigen vector is the right eigen vector corresponding to
    // the smallest singular value 0
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(rot_i, Eigen::ComputeFullV);
    Eigen::Vector3d rot_axis = svd.matrixV().col(2);

    // to handle the direction of the rotation axis. because the calculated
    // angle is always within [0, pi], for other true angle value, the axis
    // direction should flipped so that the [angle] can be fit to [0, pi].
    double ux = fabs(rot_axis(0));
    double uy = fabs(rot_axis(1));
    double uz = fabs(rot_axis(2));

    // [sin angle] is always >= 0, so the sign of uz_sin is up to uz, so is the
    // uy_sin anf ux_sin.
    double uz_sin = rot(1,0)-rot(0,1);
    double uy_sin = rot(0,2)-rot(2,0);
    double ux_sin = rot(2,1)-rot(1,2);

    axis[0] = ux*getSign(ux_sin);
    axis[1] = uy*getSign(uy_sin);
    axis[2] = uz*getSign(uz_sin);
}



////////////////////////////////////////////////////////////
/// \brief rotationMatrixToAxisAngle_simple:
/// decomposite the rotation matrix to extract the rotating axis and angle.
/// simpler approach
/// \param rotmat: input, length 9
/// \param axis: output, rotaion axis, length 3;
/// \param angle: output, rotation angle, [0, pi]
///
void rotationMatrixToAxisAngle_simple(const double *rotmat, double *axis, double &angle)
{
    // compute the matrix trace to get rotation angle.
    // here rotation angle is in the interval [0, pi],
    // so attention should paid to the rotation angle value outside this
    // interval.
    double trace = rotmat[0]+rotmat[4]+rotmat[8];
    double cos_theta = 0.5*(trace-1.0);

    // three conditions for rotation angle
    if(fabs(cos_theta-1)<DBL_EPSILON)
    {
        // if angle = 0, the there is no rotation at all,
        // we can assign any value to the axis. here simply choose the x axis
        axis[0] = 1.0;
        axis[1] = 0.0;
        axis[2] = 0.0;

        angle =0.0;
    }
    else if((fabs(cos_theta+1)<2*DBL_EPSILON))
    {
        // if the rotation angle is equal to PI
        double ux_sq = (rotmat[0]+1)/2.0;
        double uy_sq = (rotmat[4]+1)/2.0;
        double uz_sq = (rotmat[8]+1)/2.0;
        double ux = sqrt(ux_sq);
        double uy = sqrt(uy_sq);
        double uz = sqrt(uz_sq);
        axis[0] = ux;
        axis[1] = uy*getSign(rotmat[1]+rotmat[3]);
        axis[2] = uz*getSign(rotmat[2]+rotmat[6]);

        angle = M_PI;
    }
    else
    {
        angle = acos(cos_theta);
        // other condition
        // notice that sin_theta is always positive.
        double sin_theta = sqrt(1-cos_theta*cos_theta);
        axis[0] = (rotmat[7]-rotmat[5])/2.0/sin_theta;
        axis[1] = (rotmat[2]-rotmat[6])/2.0/sin_theta;
        axis[2] = (rotmat[3]-rotmat[1])/2.0/sin_theta;
    }
}


////////////////////////////////////////////////////////////
/// \brief rotationMatrixToQuaternion
/// extract Quaternion from rotation matrix
/// \param rotmat: input, rotation matrix
/// \param quaterntion, output
///
void rotationMatrixToQuaternion(const double *rotmat, double *quaterntion)
{
    // first to calsulate the rotation axis and rotation angle
    double axis[3] = {0.0};
    double angle = 0.0;
    rotationMatrixToAxisAngle(rotmat, axis, angle);

    // then construct the quaternion with the axis and angle
    double c = cos(angle/2);
    double s = sin(angle/2);
    quaterntion[0] = c;
    quaterntion[1] = axis[0]*s;
    quaterntion[2] = axis[1]*s;
    quaterntion[3] = axis[2]*s;
}


////////////////////////////////////////////////////////////
/// \brief quaternionToRotationMatrix_indirect
/// \param quaternion, input
/// \param rotmat, output, rotation matrix
///
void quaternionToRotationMatrix_indirect(const double *quaternion, double *rotmat)
{
    // convert the quaternion to axis and angle first
    double theta_h = acos(quaternion[0]);
    double sin_theta_h = sin(theta_h);
    double angle = theta_h*2;
    double axis[3] = { quaternion[1]/sin_theta_h,
                       quaternion[2]/sin_theta_h,
                       quaternion[3]/sin_theta_h};
    // use the axis and angle to construct the rotation matrix
    axisAngleToRotationMatrix(axis, angle, rotmat);
}


////////////////////////////////////////////////////////////
/// \brief quaternionToRotationMatrix
/// compute rotation matrix directly from quaternion
/// \param quaternion
/// \param rotmat
///
void quaternionToRotationMatrix(const double *quaternion, double *rotmat)
{
    double q0 = quaternion[0];
    double q1 = quaternion[1];
    double q2 = quaternion[2];
    double q3 = quaternion[3];

    rotmat[0] = 1-2*q2*q2-2*q3*q3; rotmat[1] = 2*(q1*q2-q3*q0);   rotmat[2] = 2*(q1*q3+q2*q0);
    rotmat[3] = 2*(q1*q2+q3*q0);   rotmat[4] = 1-2*q1*q1-2*q3*q3; rotmat[5] = 2*(q2*q3-q1*q0);
    rotmat[6] = 2*(q1*q3-q2*q0);   rotmat[7] = 2*(q2*q3+q1*q0);   rotmat[8] = 1-2*q1*q1-2*q2*q2;
}


////////////////////////////////////////////////////////////
/// \brief eulerAngleToQuaternion
/// convert euler angles to quaternion, rotatioin order: z->y->x
/// \param eulerAngles: input, euler angles, [angle_x, angle_y, angle_z]
/// \param quaternion: output
///
void eulerAngleToQuaternion_ZYX(const double *eulerAngles, double *quaternion)
{
    // euler angle
    double psi = eulerAngles[2]; // z axis
    double theta = eulerAngles[1]; // y axis
    double phi = eulerAngles[0]; // x axis

    double c_psi_h = cos(psi/2.0);
    double s_psi_h = sin(psi/2.0);
    double c_theta_h = cos(theta/2.0);
    double s_theta_h = sin(theta/2.0);
    double c_phi_h = cos(phi/2.0);
    double s_phi_h = sin(phi/2.0);

    quaternion[0] = c_phi_h*c_theta_h*c_psi_h + s_phi_h*s_theta_h*s_psi_h;
    quaternion[1] = s_phi_h*c_theta_h*c_psi_h - c_phi_h*s_theta_h*s_psi_h;
    quaternion[2] = c_phi_h*s_theta_h*c_psi_h + s_phi_h*c_theta_h*s_psi_h;
    quaternion[3] = c_phi_h*c_theta_h*s_phi_h - s_phi_h*s_theta_h*c_psi_h;
}


////////////////////////////////////////////////////////////
/// \brief quaternionToEulerAngle
/// convert quaternion to euler angles
/// \param quaternion: input
/// \param eulerAngles: output, [angle_x, angle_y, angle_z]
///
void quaternionToEulerAngle(const double *quaternion, double *eulerAngles)
{
    double q0 = quaternion[0];
    double q1 = quaternion[1];
    double q2 = quaternion[2];
    double q3 = quaternion[3];
    double phi = atan2(2*(q0*q1+q2*q3), 1-2*(q1*q1+q2*q2));
    double theta = asin(2*(q0*q2-q3*q1));
    double psi = atan2(2*(q0*q3+q1*q2), 1-2*(q2*q2+q3*q3));

    eulerAngles[0] = phi; // x axis
    eulerAngles[1] = theta; // y axis
    eulerAngles[2] = psi; // z axis
}


////////////////////////////////////////////////////////////
/// \brief rotationMatrixToEulerAngle_viaQuaternion
/// extract Euler rotation angles from rotation matrix, via quaternion convertion
/// \param rotmat: input, rotation matrix
/// \param eulerAngles: output, [angle_x, angle_y, angle_z];
///
void rotationMatrixToEulerAngle_viaQuaternion(const double *rotmat, double *eulerAngles)
{
    double quaternion[4] = {0.0};
    rotationMatrixToQuaternion(rotmat, quaternion);
    quaternionToEulerAngle(quaternion, eulerAngles);
}


////////////////////////////////////////////////////////////
/// \brief eulerAnglesToRotationMatrix_indirect
/// compute rotation matrix indiectly, by means of quaternion
/// \param eulerAngles: input, [angle_x, angle_y, angle_z]
/// \param rotmat: output, rotation matrix
///
void eulerAnglesToRotationMatrix_indirect(const double *eulerAngles, double *rotmat)
{
    double quaternion[4] = {0.0};
    eulerAngleToQuaternion_ZYX(eulerAngles, quaternion);
    quaternionToRotationMatrix(quaternion, rotmat);
}

