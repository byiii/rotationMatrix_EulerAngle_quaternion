#include <iostream>
#include "axis_angle.h"
#include "helpers.h"

#include <Eigen/Geometry>

#include <cmath>

int main(int argc, char** argv)
{
    double axis[] = {0, -1, 1};
    double angle = 80.0/180.0*M_PI;
    double rotmat[9] = {0.0};

    axisAngleToRotationMatrix(axis, angle, rotmat, false);

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

    double axis_3[3] = {0.0};
    double angle_3 = 0.0;
    rotationMatrixToAxisAngle_simple(rotmat, axis_3, angle_3);
    std::cout << "result: " << std::endl;
    std::cout << "axis: ";
    printVector(axis_3, 3);
    std::cout << "angle: " << angle_3/M_PI*180 << std::endl;

    double quaternion_1[4] = {0.0};
    rotationMatrixToQuaternion(rotmat, quaternion_1);
    std::cout << "quaternion result: ";
    printVector(quaternion_1, 4);
    std::cout << std::endl;

    double rotmat_2[9] = {0.0};
    quaternionToRotationMatrix(quaternion_1, rotmat_2);
    std::cout << "rotmat_2: ";
    printMatrix(rotmat_2, 3, 3);

    double quaternion_2[4] = {0.0};
    rotationMatrixToQuaternion(rotmat_2, quaternion_2);
    std::cout << "quaternion 2 result: ";
    printVector(quaternion_2, 4);
    std::cout << std::endl;

    double eulerAngles_1[3] = {0.0};
    rotationMatrixToEulerAngle_XYZ(rotmat, eulerAngles_1);
    for(size_t i=0; i<3; ++i)
    {
        eulerAngles_1[i] = eulerAngles_1[i]/M_PI*180;
    }
    std::cout << "euler angle result: ";
    printVector(eulerAngles_1, 3);
    std::cout << std::endl;

    return 0;
}
