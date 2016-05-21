#include <iostream>
#include "axis_angle.h"
#include "helpers.h"

#include <Eigen/Geometry>

#include <cmath>

int main(int argc, char** argv)
{
    double axis[] = {-1, 1, 0};
    double angle = -30.0/180.0*M_PI;
    double rotmat[9] = {0.0};

    axisAngleToRotationMatrix(axis, angle, rotmat);

    std::cout << "result: " << std::endl;
    printMatrix(rotmat, 3, 3);

    //    Eigen::Vector3d rot_axis_o(axis[0], axis[1], axis[2]);
    //    Eigen::AngleAxisd rotation_true(angle, rot_axis_o.normalized());
    //    std::cout << rotation_true.matrix() << std::endl;

    double axis_2[3] = {0.0};
    double angle_2 = 0.0;
    rotationMatrixToAxisAngle(rotmat, axis_2, angle_2);
    std::cout << "result: " << std::endl;
    std::cout << "axis: ";
    printVector(axis_2, 3);
    std::cout << "angle: " << angle_2/M_PI*180 << std::endl;

    return 0;
}
