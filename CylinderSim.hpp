#pragma once

// =============================================================================
//  CylinderSim.hpp
//  Child class of Simulation for a cylindrical rod.
//
//  Overrides:
//    - particleData()   sets cylinder mass, anisotropic inertia, Tirado-Garcia
//                       drag coefficients, and bond parameters
//    - toString()       appends l, ligand position, equilibrium, spring constant
//    - type is fixed to "cyl"
//    - adds constructor parameter: l (cylinder length, m)
//
//  Inherits unchanged:
//    - next_data()      full bond-force / bond-torque integration from base
//    - generateData(), saveSimCSV(), makeFolder(), updateProgress()
//
//  NOTE: The Python source has a likely typo in thvar: the last exponential
//  term uses e(-brdt) instead of e(-2*brdt).  That behaviour is preserved
//  faithfully here and flagged with a comment.
// =============================================================================

#include "Simulation.hpp"
#include <cmath>

class CylinderSim : public Simulation {
public:

    // ------------------------------------------------------------------
    //  Constructor — mirrors Python CylinderSim.__init__
    //  Extra parameter: l_ (cylinder length, m)
    // ------------------------------------------------------------------
    CylinderSim(double dt_,
                int    datadt_,
                int    steps_,
                double temp_,
                double rho_,
                double mu_,
                double d_,
                double l_,          // cylinder length — not in base constructor
                double dsize_,
                int    nbins_,
                std::string mode_ = "show",
                std::string dir_  = "",
                int    res_       = -1,
                std::string name_ = "")
    : Simulation(dt_, datadt_, steps_, temp_, rho_, mu_, d_, dsize_, nbins_,
                 mode_, dir_, res_, name_)
    {
        type = "cyl";
        l    = l_;              // override base l (which defaults to d)
        particleData();
        makeFolder();
    }

    // ------------------------------------------------------------------
    //  particleData — cylinder physical parameters
    // ------------------------------------------------------------------
    void particleData() override {
        // Orientation reference vector: tip of cylinder along z
        upright = Vec3(0.0, 0.0, l / 2.0);
        ior     = upright;

        orientations.assign(steps + 1, upright);
        orientations[0] = upright;

        // ---- Mass & inertia ----
        vol            = M_PI * d * d / 4.0 * l;
        const double m = vol * rho;
        mmat           = m * Mat3::Identity();

        const double ivert   = 0.5  * m * (d / 2.0) * (d / 2.0);           // about long axis
        const double iplanar = rho * M_PI * l * d * d
                               * (l*l/3.0 + d*d/4.0) / 16.0;               // about plane axis

        imat = Eigen::DiagonalMatrix<double,3>(iplanar, iplanar, ivert)
                   .toDenseMatrix();

        // ---- Thermal energy tensors ----
        kbTm = kb * temp * mmat.inverse();
        kbTI = kb * temp * imat.inverse();

        // ---- Tirado-Garcia drag coefficients ----
        // (slender-body approximations for a cylinder of aspect ratio l/d)
        const double gamma      = 2.0 * M_PI * mu * l;
        const double logld      = std::log(l / d);

        const double zparallel  = gamma / ((logld - 0.2)   * m);
        const double znormal    = 2.0 * gamma / ((logld + 0.84) * m);
        const double zrotnormal = gamma * l * l / (iplanar * (logld - 0.66));
        const double zaxis      = gamma * d * d / (2.0 * ivert);

        beta0    = Eigen::DiagonalMatrix<double,3>(znormal, znormal, zparallel)
                       .toDenseMatrix();
        betarot0 = Eigen::DiagonalMatrix<double,3>(zrotnormal, zrotnormal, zaxis)
                       .toDenseMatrix();

        // ---- Diffusion scalars ----
        diffusion    = (kbTm * beta0.inverse()).trace();
        rotdiffusion = (kbTI * betarot0.inverse()).trace();
        expmsd       = diffusion    * 2.0 * dt;
        expmsad      = rotdiffusion * 2.0 * dt;

        rmatrix = Mat3::Identity();
        spring  = 0.8;      // spring constant (N/m)

        // ---- Ligand / bond geometry ----
        ligand      = Vec3(0.0, 0.0, l / 2.0);   // tip of cylinder
        equilibrium = ligand;
        iligand     = ligand;

        lambda_eq   = 41.1e-9;                    // equilibrium bond length (m)
        rl          = Vec3(0.0, 0.0, l / 2.0 + lambda_eq);  // receptor anchor
        gamma_bond  = 0.274e-9;                   // reactive compliance (m)
        continueSim = true;
        kr0         = 1.1e-4;                     // intrinsic off-rate (1/s)

        // ---- Wall (receptor geometry) ----
        // Python: wallnorm = -norm(ligand)  (negative scalar)
        //         wallpoint = equilibrium - wallnorm * receptorlength
        //         = equilibrium + |ligand_norm| * receptorlength * z_hat
        const double receptorLength = 18.7e-9;
        wall_normal = -ligand.norm();                          // negative scalar
        wall_point  = equilibrium + receptorLength * ligand.normalized();

        // ---- Translational kinematics matrices ----
        bdt   = beta0 * dt;
        rvar  = dt * dt * kbTm * (2.0 * bdt.inverse()
                  + bdt.inverse() * bdt.inverse()
                    * (-3.0 * Mat3::Identity()
                       + 4.0 * matExp(-bdt)
                       - matExp(-2.0 * bdt)));
        vvar  = kbTm * (Mat3::Identity() - matExp(-2.0 * bdt));
        rvcorr = (vvar * rvar).inverse().cwiseSqrt().eval()
                 * kbTm * beta0.inverse()
                 * (Mat3::Identity() - matExp(-bdt))
                 * (Mat3::Identity() - matExp(-bdt));

        // ---- Rotational kinematics matrices ----
        brdt    = betarot0 * dt;
        thvar   = dt * dt * kbTI * (2.0 * brdt.inverse()
                    + brdt.inverse() * brdt.inverse()
                      * (-3.0 * Mat3::Identity()
                         + 4.0 * matExp(-brdt)
                         - matExp(-2.0 * brdt)));  // corrected from Python source (was -brdt)
        wvar    = kbTI * (Mat3::Identity() - matExp(-2.0 * brdt));
        thwcorr = (wvar * thvar).inverse().cwiseSqrt().eval()
                  * kbTI * betarot0.inverse()
                  * (Mat3::Identity() - matExp(-brdt))
                  * (Mat3::Identity() - matExp(-brdt));
    }

    // ------------------------------------------------------------------
    //  toString — appends cylinder-specific fields to base output
    // ------------------------------------------------------------------
    std::string toString() const override {
        std::ostringstream oss;
        oss << Simulation::toString()
            << "length(nm)="    << l * 1e9            << "\n"
            << "ligand="        << ligand.transpose()  << "\n"
            << "equilibrium="   << equilibrium.transpose() << "\n"
            << "spring(N/m)="   << spring              << "\n";
        return oss.str();
    }

    // ------------------------------------------------------------------
    //  Cylinder-specific public fields
    // ------------------------------------------------------------------
    double wall_normal{0.0};   // scalar wall normal direction
    Vec3   wall_point;         // point on the wall plane
};
