// Convert this to c++
// def rotateMatrix(dth,rmatrix):
//     theta = np.linalg.norm(dth)
//     dth /= theta  # Ensure the axis is a unit vector
//     u_x, u_y, u_z = dth
//     cos_t = np.cos(theta)
//     sin_t = np.sin(theta)
//     one_minus_cos_t = 1 - cos_t

//     return np.array([
//         [
//             cos_t + u_x**2 * one_minus_cos_t,
//             u_x * u_y * one_minus_cos_t - u_z * sin_t,
//             u_x * u_z * one_minus_cos_t + u_y * sin_t
//         ],
//         [
//             u_y * u_x * one_minus_cos_t + u_z * sin_t,
//             cos_t + u_y**2 * one_minus_cos_t,
//             u_y * u_z * one_minus_cos_t - u_x * sin_t
//         ],
//         [
//             u_z * u_x * one_minus_cos_t - u_y * sin_t,
//             u_z * u_y * one_minus_cos_t + u_x * sin_t,
//             cos_t + u_z**2 * one_minus_cos_t
//         ]
//     ]) @ rmatrix


#include <Eigen/Dense>
#include <vector> 
#include <iostream>



Eigen::Matrix3d rotateMatrix(Eigen::Vector3d dth, Eigen::Matrix3d rmatrix) {
    double theta = dth.norm();
    dth /= theta;  // Ensure the axis is a unit vector
    double u_x = dth(0);
    double u_y = dth(1);
    double u_z = dth(2);
    double cos_t = std::cos(theta);
    double sin_t = std::sin(theta);
    double one_minus_cos_t = 1 - cos_t;

    Eigen::Matrix3d rotationMatrix;
    rotationMatrix << 
        cos_t + u_x * u_x * one_minus_cos_t,
        u_x * u_y * one_minus_cos_t - u_z * sin_t,
        u_x * u_z * one_minus_cos_t + u_y * sin_t,
        
        u_y * u_x * one_minus_cos_t + u_z * sin_t,
        cos_t + u_y * u_y * one_minus_cos_t,
        u_y * u_z * one_minus_cos_t - u_x * sin_t,
        
        u_z * u_x * one_minus_cos_t - u_y * sin_t,
        u_z * u_y * one_minus_cos_t + u_x * sin_t,
        cos_t + u_z * u_z * one_minus_cos_t;

    return rotationMatrix * rmatrix;
}
// test code
int main() {
    Eigen::Vector3d dth(0.001, 0.002, 0.003);
    Eigen::Matrix3d rmatrix = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d result = rotateMatrix(dth, rmatrix);
    std::cout << "Resulting Matrix:\n" << result << std::endl;
    return 0;
}