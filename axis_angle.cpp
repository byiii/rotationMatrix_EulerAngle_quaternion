#include "axis_angle.h"
#include "helpers.h"

#include <cmath>
#include <cfloat>
#include <Eigen/Dense>

////////////////////////////////////////////////////////////
/// \brief axisAngleRotationMatrix
/// generate rotation matrix given rotation axis and angle
/// make sure axis is of length 3, rotmat is of length 9;
/// \param axis: input, length 3
/// \param angle: input
/// \param rotmat: output, length 9
///
void axisAngleToRotationMatrix(const double *axis, double angle, double *rotmat)
{
    double ux = axis[0];
    double uy = axis[1];
    double uz = axis[2];

    // check if the axis is an unit vector, if not, convert to unit vector
    double axis_length = std::sqrt(ux*ux+uy*uy+uz*uz);
    if(DBL_EPSILON<fabs(axis_length-1.0))
    {
        ux = ux/axis_length;
        uy = uy/axis_length;
        uz = uz/axis_length;
        axis_length = 1.0;
    }

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
/// \param axis: output, length 3;
/// \param angle: output
///
void rotationMatrixToAxisAngle(double *mat, double *axis, double &angle)
{
    double rotmat[9] = {0.0};
    for(size_t i=0; i<9; ++i)
        rotmat[i] = mat[i];

    double trace = rotmat[0]+rotmat[4]+rotmat[8];
    double cos_theta = 0.5*(trace-1);
    angle = std::acos(cos_theta);

    Eigen::Matrix3d rot = Eigen::Map<Eigen::MatrixXd>(rotmat, 3, 3).transpose();
    Eigen::Matrix3d rot_i = rot - Eigen::Matrix3d::Identity();

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(rot_i, Eigen::ComputeFullV);
    Eigen::Vector3d rot_axis = svd.matrixV().col(2);

    axis[0] = rot_axis(0);
    axis[1] = rot_axis(1);
    axis[2] = rot_axis(2);

    double uz_sin = rot(1,0)-rot(0,1);
    double uy_sin = rot(0,2)-rot(2,0);
    double ux_sin = rot(2,1)-rot(1,2);
    if(ux_sin>DBL_EPSILON)
        axis[0] = fabs(axis[0]);
    else
        axis[0] = -fabs(axis[0]);
    if(uy_sin>DBL_EPSILON)
        axis[1] = fabs(axis[1]);
    else
        axis[1] = -fabs(axis[1]);
    if(uz_sin>DBL_EPSILON)
        axis[2] = fabs(axis[2]);
    else
        axis[2] = -fabs(axis[2]);
}
