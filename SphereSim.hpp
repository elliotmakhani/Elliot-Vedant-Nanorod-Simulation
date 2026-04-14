#pragma once

// =============================================================================
//  SphereSim.hpp
//  Child class of Simulation for a sphere.
//
//  Overrides:
//    - particleData()  sets physical parameters analytically for a sphere
//                      (mass, inertia, Stokes drag, kinematics matrices)
//    - type is fixed to "sph"
//
//  Inherits unchanged:
//    - next_data()     full bond-force / bond-torque integration from base
//    - generateData(), saveSimCSV(), makeFolder(), updateProgress()
// =============================================================================

#include "Simulation.hpp"

class SphereSim : public Simulation {
public:

    // ------------------------------------------------------------------
    //  Constructor — mirrors Python SphereSim.__init__
    // ------------------------------------------------------------------
    SphereSim(double dt_,
              int    datadt_,
              int    steps_,
              double temp_,
              double rho_,
              double mu_,
              double d_,
              double dsize_,
              int    nbins_,
              std::string mode_ = "show",
              std::string dir_  = "",
              int    res_       = -1,
              std::string name_ = "")
    : Simulation(dt_, datadt_, steps_, temp_, rho_, mu_, d_, dsize_, nbins_,
                 mode_, dir_, res_, name_)
    {
        type = "sph";
        particleData();       // sphere-specific physical setup
        makeFolder();         // base-class folder creation (uses type = "sph")
    }

    // ------------------------------------------------------------------
    //  particleData — sphere physical parameters (overrides base stub)
    // ------------------------------------------------------------------
    void particleData() override {
        // Orientation reference vector: tip of sphere along z
        upright = Vec3(0.0, 0.0, d / 2.0);
        ior     = upright;

        orientations.assign(steps + 1, upright);
        orientations[0] = upright;

        // Mass / inertia
        vol           = M_PI * d * d * d / 6.0;
        const double m = vol * rho;
        mmat          = m * Mat3::Identity();
        imat          = (m * d * d / 10.0) * Mat3::Identity();

        // Thermal energy tensors
        kbTm = kb * temp * mmat.inverse();
        kbTI = kb * temp * imat.inverse();

        // Drag (Stokes)  β = 3πμd / m,   β_rot = π μ d³ / I
        beta0    = mmat.inverse()  * (3.0 * d * mu * M_PI);
        betarot0 = imat.inverse()  * (d * d * d * M_PI * mu);

        // Diffusion coefficients
        diffusion    = (kbTm * beta0.inverse()).trace();
        rotdiffusion = (kbTI * betarot0.inverse()).trace();
        expmsd       = diffusion    * 2.0 * dt;
        expmsad      = rotdiffusion * 2.0 * dt;

        rmatrix = Mat3::Identity();
        spring  = 1.0;

        // Translational kinematics matrices
        bdt    = beta0    * dt;
        rvar   = dt * dt * kbTm * (2.0 * bdt.inverse()
                   + bdt.inverse() * bdt.inverse()
                     * (-3.0 * Mat3::Identity()
                        + 4.0 * matExp(-bdt)
                        - matExp(-2.0 * bdt)));
        vvar   = kbTm * (Mat3::Identity() - matExp(-2.0 * bdt));
        rvcorr = (vvar * rvar).inverse().cwiseSqrt()   // element-wise on diagonal
                   .eval()
                 * kbTm * beta0.inverse()
                 * (Mat3::Identity() - matExp(-bdt))
                 * (Mat3::Identity() - matExp(-bdt));

        // Rotational kinematics matrices
        brdt    = betarot0 * dt;
        thvar   = dt * dt * kbTI * (2.0 * brdt.inverse()
                    + brdt.inverse() * brdt.inverse()
                      * (-3.0 * Mat3::Identity()
                         + 4.0 * matExp(-brdt)
                         - matExp(-2.0 * brdt)));
        wvar    = kbTI * (Mat3::Identity() - matExp(-2.0 * brdt));
        thwcorr = (wvar * thvar).inverse().cwiseSqrt()
                    .eval()
                  * kbTI * betarot0.inverse()
                  * (Mat3::Identity() - matExp(-brdt))
                  * (Mat3::Identity() - matExp(-brdt));
    }

    // next_data() is inherited from Simulation — bond force and torque
    // are fully handled there. SphereSim only needs sphere-specific
    // physical parameters, which particleData() provides above.
};
