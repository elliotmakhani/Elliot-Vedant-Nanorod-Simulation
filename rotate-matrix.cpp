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
