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

#include <array>
#include <cmath>

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;

static inline double norm3(const Vec3& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline Mat3 matMul(const Mat3& A, const Mat3& B) {
    Mat3 C{};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += A[i][k] * B[k][j];
            C[i][j] = s;
        }
    }
    return C;
}


Mat3 rotateMatrix(const Vec3& dth, const Mat3& rmatrix) {
    const double theta = norm3(dth);
    constexpr double eps = 1e-12;


    if (theta < eps) return rmatrix;

    const double u_x = dth[0] / theta;
    const double u_y = dth[1] / theta;
    const double u_z = dth[2] / theta;

    const double cos_t = std::cos(theta);
    const double sin_t = std::sin(theta);
    const double one_minus_cos_t = 1.0 - cos_t;

    Mat3 R {{
        {{
            cos_t + u_x*u_x*one_minus_cos_t,
            u_x*u_y*one_minus_cos_t - u_z*sin_t,
            u_x*u_z*one_minus_cos_t + u_y*sin_t
        }},
        {{
            u_y*u_x*one_minus_cos_t + u_z*sin_t,
            cos_t + u_y*u_y*one_minus_cos_t,
            u_y*u_z*one_minus_cos_t - u_x*sin_t
        }},
        {{
            u_z*u_x*one_minus_cos_t - u_y*sin_t,
            u_z*u_y*one_minus_cos_t + u_x*sin_t,
            cos_t + u_z*u_z*one_minus_cos_t
        }}
    }};

    return matMul(R, rmatrix);
}
