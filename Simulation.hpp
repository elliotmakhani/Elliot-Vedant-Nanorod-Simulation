#pragma once

// =============================================================================
//  Simulation.hpp
//  Merged from two branches:
//    - Logic/completeness from the fully-implemented std::array branch
//    - Eigen types replacing all manual Vec3/Vec2/Mat3 math helpers
// =============================================================================

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#ifndef M_PI
static constexpr double M_PI = 3.14159265358979323846;
#endif
#include <Eigen/Dense>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <map>

namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
//  Convenience aliases
// ---------------------------------------------------------------------------
using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;

// ---------------------------------------------------------------------------
//  Free helper: Rodrigues rotation  (kept as free function to match original)
// ---------------------------------------------------------------------------
inline Mat3 rotateMatrix(const Vec3& dth, const Mat3& R) {
    const double th = dth.norm();
    if (th == 0.0) return R;
    const Vec3 u  = dth / th;
    const double c   = std::cos(th);
    const double s   = std::sin(th);
    const double omc = 1.0 - c;

    Mat3 Rinc;
    Rinc <<
        c + u.x()*u.x()*omc,          u.x()*u.y()*omc - u.z()*s,    u.x()*u.z()*omc + u.y()*s,
        u.y()*u.x()*omc + u.z()*s,    c + u.y()*u.y()*omc,           u.y()*u.z()*omc - u.x()*s,
        u.z()*u.x()*omc - u.y()*s,    u.z()*u.y()*omc + u.x()*s,    c + u.z()*u.z()*omc;

    return Rinc * R;
}

// theta(x,y) = atan2(x,y)  [matches Python convention]
inline double theta(double x, double y) { return std::atan2(x, y); }

// ---------------------------------------------------------------------------
//  Simulation class
// ---------------------------------------------------------------------------
class Simulation {
public:

    // ------------------------------------------------------------------
    //  Constructor  (mirrors Python __init__ + array-branch constructor)
    // ------------------------------------------------------------------
    Simulation(double dt_,
               int    datadt_,
               int    steps_,        // raw steps (Python: steps before /datadt)
               double temp_,
               double rho_,
               double mu_,
               double d_,
               double dsize_,
               int    nbins_,
               std::string mode_ = "show",
               std::string dir_  = "",
               int    res_       = -1,   // -1 => res = steps_
               std::string name_ = "")
    : type(""),
      dt(dt_),
      datadt(datadt_),
      total_raw_steps(steps_),
      steps(std::max(1, steps_ / std::max(1, datadt_))),
      last_data_idx(steps),
      totalT(dt_ * steps_),
      temp(temp_),
      rho(rho_),
      mu(mu_),
      d(d_),
      l(d_),
      dsize(dsize_),
      nbins(nbins_),
      mode(std::move(mode_)),
      dir(std::move(dir_)),
      name(std::move(name_)),
      res((res_ < 0) ? steps_ : res_)
    {
        if (!dir.empty() && dir.back() != '/') dir.push_back('/');

        // Position / kinematics storage  (size = steps+1)
        const int N = steps + 1;
        positions.resize(N, Vec3::Zero());
        velocities.resize(N, Vec3::Zero());
        thetaz.resize(N, Vec2::Zero());
        orientations.resize(N);
        for (auto& o : orientations) o = Vec3(1.0, 0.0, 0.0);
        omegas.resize(N, Vec3::Zero());
        omegas_objframe.resize(N, Vec3::Zero());
        sqdis.resize(N, 0.0);
        sqangdis.resize(N, 0.0);

        times.resize(N);
        for (int i = 0; i < N; ++i)
            times[i] = static_cast<double>(i) * (dt * datadt);

        // Particle / bond defaults
        upright     = Vec3(1.0, 0.0, 0.0);
        ior         = upright;
        ligand      = upright;
        iligand     = ligand;
        equilibrium = ligand;
        rmatrix     = Mat3::Identity();

        continueSim = true;
        rl          = Vec3::Zero();
        lambda_eq   = 1.0;
        kr0         = 0.0;
        gamma_bond  = 0.0;
        spring      = 1.0;
        sqdr             = 0.0;
        sqdth            = 0.0;
        iter_progress_pos_ = -1;

        // RNG
        std::random_device rd;
        rng.seed(rd());
    }

    virtual ~Simulation() = default;

    // ------------------------------------------------------------------
    //  makeFolder  (complete – array branch logic, fs::path used properly)
    // ------------------------------------------------------------------
    void makeFolder() {
        if (mode == "save") {
            int i = 1;
            std::string attempt = name + type;
            while (true) {
                fs::path p = fs::path(dir) / attempt;
                std::error_code ec;
                bool created = fs::create_directory(p, ec);

                if (!created || ec) {
                    ++i;
                    if (i == 2)
                        attempt = name + type + "v2";
                    else {
                        // strip last digit / suffix and append new index
                        attempt = name + type + "v" + std::to_string(i);
                    }
                } else {
                    std::ofstream file((p / "parameters.txt").string());
                    file << toString();
                    break;
                }
            }
        }
        name += type;
    }

    // ------------------------------------------------------------------
    //  particleData  – initialises dynamics matrices using Eigen
    // ------------------------------------------------------------------
    virtual void particleData() {
        upright = Vec3(1.0, 0.0, 0.0);
        ior     = upright;

        orientations.assign(steps + 1, upright);

        vol          = 0.0;
        mmat         = Mat3::Identity();
        imat         = Mat3::Identity();
        kbTm         = Mat3::Identity();
        kbTI         = Mat3::Identity();
        beta0        = Mat3::Identity();
        betarot0     = Mat3::Identity();
        diffusion    = 0.0;
        rotdiffusion = 0.0;
        expmsd       = diffusion    * 2.0 * dt;
        expmsad      = rotdiffusion * 2.0 * dt;
        rmatrix      = Mat3::Identity();

        // translational kinematics pre-computation
        bdt   = beta0 * dt;
        rvar  = computeRVar(bdt, kbTm) * (dt * dt);
        vvar  = computeVVar(bdt, kbTm);
        rvcorr = computeCorr(vvar, rvar, kbTm, beta0, bdt);

        // rotational kinematics pre-computation
        brdt    = betarot0 * dt;
        thvar   = computeRVar(brdt, kbTI) * (dt * dt);
        wvar    = computeVVar(brdt, kbTI);
        thwcorr = computeCorr(wvar, thvar, kbTI, betarot0, brdt);

        spring      = 1.0;
        ligand      = upright;
        iligand     = ligand;
        equilibrium = ligand;
    }

    // ------------------------------------------------------------------
    //  next_data  – one integration step (fully implemented, Eigen math)
    // ------------------------------------------------------------------
    virtual void next_data() {
        Vec3 bondforce = Vec3::Zero();

        if (continueSim) {
            Vec3 rr      = ir + iligand;                           // ligand position
            double ilambda = (rl - rr).norm();
            if (ilambda > 0.0)
                bondforce = spring * (1.0 - lambda_eq / ilambda) * (rl - rr);

            double delta = std::abs(lambda_eq - ilambda);
            double kr    = kr0 * std::exp(gamma_bond * spring * delta / (kb * temp));
            double p_r   = 1.0 - std::exp(-kr * dt);

            if (uniform01(rng) < p_r)
                continueSim = false;
        }

        // ---- Translational kinematics ----
        Mat3 c0 = matExp(-(bdt));
        Mat3 c1 = (Mat3::Identity() - c0) * beta0.inverse();
        Mat3 c2 = (dt * Mat3::Identity() - c1) * beta0.inverse();

        Vec3 randomr = sampleNormal3();
        randomr = diagSqrtMul(rvar, randomr);

        Vec3 randomv = sampleCorrelatedVelocity(randomr, vvar, rvcorr, rvar);

        // Rotate noise into lab frame
        randomr = rmatrix * randomr;
        randomv = rmatrix * randomv;

        Vec3 k  = mmat.inverse() * bondforce;
        Vec3 dr = randomr
                + rmatrix * c1 * rmatrix.transpose() * iv
                + rmatrix * c2 * rmatrix.transpose() * k;

        sqdr  = dr.squaredNorm();
        ir   += dr;
        iv    = rmatrix * c0 * rmatrix.transpose() * iv
              + randomv
              + rmatrix * c1 * rmatrix.transpose() * k;

        // ---- Rotational kinematics ----
        c0 = matExp(-(brdt));
        c1 = (Mat3::Identity() - c0) * betarot0.inverse();
        c2 = (dt * Mat3::Identity() - c1) * betarot0.inverse();

        Vec3 torque = iligand.cross(bondforce);

        Vec3 randomth = sampleNormal3();
        randomth = diagSqrtMul(thvar, randomth);

        Vec3 randomw = sampleCorrelatedVelocity(randomth, wvar, thwcorr, thvar);

        randomth = rmatrix * randomth;
        randomw  = rmatrix * randomw;

        Vec3 dth = randomth
                 + rmatrix * c1 * rmatrix.transpose() * iw
                 + rmatrix * c2 * imat.inverse() * rmatrix.transpose() * torque;

        sqdth = dth.squaredNorm();

        iw = rmatrix * c0 * rmatrix.transpose() * iw
           + randomw
           + rmatrix * c1 * imat.inverse() * rmatrix.transpose() * torque;

        rmatrix = rotateMatrix(dth, rmatrix);
        iligand = rmatrix * ligand;
        ior     = rmatrix * upright;
        ith     = Vec2(theta(ior.x(), ior.y()),
                       ior.z() / std::sqrt(upright.squaredNorm()));
    }

    // ------------------------------------------------------------------
    //  generateData  – main simulation loop (fully implemented)
    // ------------------------------------------------------------------
    void generateData() {
        const int rawN         = steps * datadt;
        const int progressEvery = std::max(1, rawN / 1000);

        for (int i = 1; i <= rawN; ++i) {
            next_data();

            if (!continueSim) {
                last_data_idx = std::min(steps, i / std::max(1, datadt));
                updateProgress("Stopped at " + std::to_string(i * dt) + " s");
                break;
            }

            if (i % datadt == 0) {
                int idx = i / datadt;
                if (idx >= 0 && idx <= steps) {
                    positions[idx]      = ir;
                    velocities[idx]     = iv;
                    thetaz[idx]         = ith;
                    omegas[idx]         = iw;
                    omegas_objframe[idx] = rmatrix.transpose() * iw;
                    orientations[idx]   = rmatrix * upright;
                    sqdis[idx]          = sqdr;
                    sqangdis[idx]       = sqdth;
                }
            }

            if (mode == "save" && (i % progressEvery == 0)) {
                writeIterProgress(i);
            }
        }

        last_data_idx = std::min(last_data_idx, steps);
    }

    // ------------------------------------------------------------------
    //  saveSimCSV  – writes all output to CSV files
    // ------------------------------------------------------------------
    void saveSimCSV() {
        if (mode != "save") return;

        fs::path base = fs::path(dir) / name;
        fs::create_directories(base);

        writeVec3CSV((base / "positions.csv").string(),    positions,       last_data_idx, {"x","y","z"});
        writeVec3CSV((base / "velocities.csv").string(),   velocities,      last_data_idx, {"x","y","z"});
        writeVec2CSV((base / "thetaz.csv").string(),       thetaz,          last_data_idx, {"theta","z"});
        writeVec3CSV((base / "omegas.csv").string(),       omegas,          last_data_idx, {"roll","pitch","yaw"});
        writeVec3CSV((base / "omegas_objframe.csv").string(), omegas_objframe, last_data_idx,
                     {"object_x","object_y","object_z"});
        writeVec3CSV((base / "orientations.csv").string(), orientations,    last_data_idx, {"x","y","z"});
        writeScalarCSV((base / "squared_displacement.csv").string(), sqdis, sqangdis, last_data_idx);

        updateProgress("CSV data saved");
    }

    // Append a milestone message to progress.txt.
    // Override in subclasses to redirect to stdout or a GUI.
    virtual void updateProgress(const std::string& msg) {
        if (mode != "save") return;
        fs::path p = fs::path(dir) / name / "progress.txt";
        std::ofstream f(p.string(), std::ios::app);
        f << msg << "\n";
    }

    // ------------------------------------------------------------------
    //  parseParams — read key=value pairs from a parameters.txt file
    //  Returns a map of all keys found; throws if file cannot be opened.
    // ------------------------------------------------------------------
    static std::map<std::string, std::string>
    parseParams(const std::string& path) {
        std::ifstream f(path);
        if (!f.is_open())
            throw std::runtime_error("parseParams: cannot open \"" + path + "\"");
        std::map<std::string, std::string> params;
        std::string line;
        while (std::getline(f, line)) {
            auto eq = line.find('=');
            if (eq == std::string::npos) continue;
            std::string key = line.substr(0, eq);
            std::string val = line.substr(eq + 1);
            // strip trailing whitespace/carriage-returns
            while (!val.empty() && (val.back() == ' ' || val.back() == '\r'))
                val.pop_back();
            params[key] = val;
        }
        return params;
    }

    // Overwrite only the last iteration-progress line so the file does not
    // accumulate 1000 "iteration N" entries.  A fixed-width sentinel line is
    // written on the first call; subsequent calls seek back and overwrite it.
    void writeIterProgress(int iteration) {
        if (mode != "save") return;
        fs::path p = fs::path(dir) / name / "progress.txt";

        const int    rawN = steps * datadt;
        const double pct  = 100.0 * iteration / rawN;

        std::ostringstream oss;
        oss << "iteration " << iteration
            << "  (" << std::fixed << std::setprecision(1) << pct << "%)";
        std::string line = oss.str();

        // Pad to fixed width so every overwrite covers the full previous line.
        const std::size_t LINE_WIDTH = 48;
        if (line.size() < LINE_WIDTH)
            line.append(LINE_WIDTH - line.size(), ' ');

        if (iter_progress_pos_ < 0) {
            // First call: append sentinel and remember its file position.
            std::fstream f(p.string(), std::ios::in | std::ios::out | std::ios::app);
            iter_progress_pos_ = static_cast<long>(f.tellp());
            f << line << "\n";
        } else {
            // Subsequent calls: seek back to sentinel and overwrite in-place.
            std::fstream f(p.string(), std::ios::in | std::ios::out);
            f.seekp(iter_progress_pos_);
            f << line << "\n";
        }
    }


    // ------------------------------------------------------------------
    //  Public member variables
    // ------------------------------------------------------------------
    std::string type;

    double dt;
    int    datadt;
    int    total_raw_steps;
    int    steps;
    int    last_data_idx;
    double totalT;
    double temp, rho, mu, d, l;

    Vec3 ir, iv, iw;
    Vec2 ith;

    std::vector<Vec3>   positions;
    std::vector<Vec3>   velocities;
    std::vector<Vec2>   thetaz;
    std::vector<Vec3>   orientations;
    std::vector<Vec3>   omegas;
    std::vector<Vec3>   omegas_objframe;
    std::vector<double> sqdis, sqangdis, times;

    double      dsize;
    int         nbins;
    std::string mode, dir, name;
    int         res;

    // Particle / dynamics fields
    Vec3 upright, ior, ligand, iligand, equilibrium;
    double vol{0.0};
    Mat3 mmat, imat, kbTm, kbTI, beta0, betarot0;
    double diffusion{0.0}, rotdiffusion{0.0}, expmsd{0.0}, expmsad{0.0};

    Mat3 rmatrix;
    Mat3 bdt, rvar, vvar, rvcorr;
    Mat3 brdt, thvar, wvar, thwcorr;

    double spring{1.0};

    // Bond fields
    bool   continueSim{true};
    Vec3   rl;
    double lambda_eq{1.0}, kr0{0.0}, gamma_bond{0.0};

    // Step scratch
    double sqdr{0.0}, sqdth{0.0};

    // File position of the overwriteable iteration-progress line (-1 = not yet written)
    long iter_progress_pos_{-1};

    static constexpr double kb = 1.380649e-23;

protected:
    // ------------------------------------------------------------------
    //  RNG  (protected so child classes can sample noise directly)
    // ------------------------------------------------------------------
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> uniform01{0.0, 1.0};
    std::normal_distribution<double>       normal01{0.0, 1.0};

    Vec3 sampleNormal3() {
        return Vec3(normal01(rng), normal01(rng), normal01(rng));
    }

    // ------------------------------------------------------------------
    //  Eigen-based math helpers  (protected so child classes can reuse)
    // ------------------------------------------------------------------

    // Matrix exponential for diagonal matrices (elementwise exp on diagonal)
    static Mat3 matExp(const Mat3& A) {
        Mat3 R = Mat3::Zero();
        R(0,0) = std::exp(A(0,0));
        R(1,1) = std::exp(A(1,1));
        R(2,2) = std::exp(A(2,2));
        return R;
    }

    // Multiply vector elementwise by sqrt of diagonal of a diagonal matrix
    static Vec3 diagSqrtMul(const Mat3& diagA, const Vec3& v) {
        return Vec3(
            v.x() * std::sqrt(std::max(0.0, diagA(0,0))),
            v.y() * std::sqrt(std::max(0.0, diagA(1,1))),
            v.z() * std::sqrt(std::max(0.0, diagA(2,2)))
        );
    }

    // Correlated velocity sample (diagonal assumption)
    Vec3 sampleCorrelatedVelocity(const Vec3& randomr_scaled,
                                  const Mat3& vvarD,
                                  const Mat3& corrD,
                                  const Mat3& rvarD)
    {
        Vec3 z = sampleNormal3();
        Vec3 out;
        for (int i = 0; i < 3; ++i) {
            const double vv  = vvarD(i,i);
            const double rr  = rvarD(i,i);
            const double c   = corrD(i,i);
            const double t1  = z(i) * std::sqrt(std::max(0.0, vv * (1.0 - c*c)));
            const double t2  = randomr_scaled(i) * c * std::sqrt(std::max(0.0, vv / (rr > 0 ? rr : 1)));
            out(i) = t1 + t2;
        }
        return out;
    }

    // rvar diagonal:  dt^2 * kbT * (2/b + (1/b^2)(-3 + 4e^-b - e^-2b))
    // NOTE: caller must multiply by dt^2 after this returns
    static Mat3 computeRVar(const Mat3& bdtD, const Mat3& kbTD) {
        Mat3 R = Mat3::Zero();
        for (int i = 0; i < 3; ++i) {
            const double b    = bdtD(i,i);
            const double kbT  = kbTD(i,i);
            if (b == 0.0) throw std::runtime_error("computeRVar: zero bdt diagonal");
            const double ib   = 1.0 / b;
            R(i,i) = kbT * (2.0*ib + ib*ib*(-3.0 + 4.0*std::exp(-b) - std::exp(-2.0*b)));
        }
        return R;
    }

    // vvar diagonal:  kbT * (1 - e^-2b)
    static Mat3 computeVVar(const Mat3& bdtD, const Mat3& kbTD) {
        Mat3 R = Mat3::Zero();
        for (int i = 0; i < 3; ++i)
            R(i,i) = kbTD(i,i) * (1.0 - std::exp(-2.0 * bdtD(i,i)));
        return R;
    }

    // rvcorr diagonal:  (1/sqrt(vv*rr)) * kbT * (1/beta) * (1 - e^-b)^2
    static Mat3 computeCorr(const Mat3& vvarD, const Mat3& rvarD,
                             const Mat3& kbTD,  const Mat3& betaD,
                             const Mat3& bdtD)
    {
        Mat3 R = Mat3::Zero();
        for (int i = 0; i < 3; ++i) {
            const double vv   = vvarD(i,i);
            const double rr   = rvarD(i,i);
            const double kbT  = kbTD(i,i);
            const double beta = betaD(i,i);
            const double b    = bdtD(i,i);
            if (vv <= 0.0 || rr <= 0.0 || beta == 0.0) { R(i,i) = 0.0; continue; }
            const double f = 1.0 - std::exp(-b);
            R(i,i) = (1.0 / std::sqrt(vv * rr)) * kbT * (1.0 / beta) * f * f;
        }
        return R;
    }

private:
    // ------------------------------------------------------------------
    //  CSV writers
    // ------------------------------------------------------------------
    static void writeVec3CSV(const std::string& path,
                             const std::vector<Vec3>& data,
                             int lastIdx,
                             const std::array<const char*,3>& hdr)
    {
        std::ofstream f(path);
        f << hdr[0] << "," << hdr[1] << "," << hdr[2] << "\n";
        for (int i = 0; i <= lastIdx; ++i)
            f << data[i].x() << "," << data[i].y() << "," << data[i].z() << "\n";
    }

    static void writeVec2CSV(const std::string& path,
                             const std::vector<Vec2>& data,
                             int lastIdx,
                             const std::array<const char*,2>& hdr)
    {
        std::ofstream f(path);
        f << hdr[0] << "," << hdr[1] << "\n";
        for (int i = 0; i <= lastIdx; ++i)
            f << data[i].x() << "," << data[i].y() << "\n";
    }

    static void writeScalarCSV(const std::string& path,
                               const std::vector<double>& trans,
                               const std::vector<double>& ang,
                               int lastIdx)
    {
        std::ofstream f(path);
        f << "Translational,Angular\n";
        for (int i = 0; i <= lastIdx; ++i)
            f << trans[i] << "," << ang[i] << "\n";
    }

protected:
    // ------------------------------------------------------------------
    //  toString for parameters file  (virtual so child classes can extend)
    // ------------------------------------------------------------------
    virtual std::string toString() const {
        std::ostringstream oss;
        oss << std::setprecision(17);
        oss << "type="             << type             << "\n"
            << "dt="               << dt               << "\n"
            << "datadt="           << datadt           << "\n"
            << "total_raw_steps="  << total_raw_steps  << "\n"
            << "steps(saved)="     << steps            << "\n"
            << "temp="             << temp             << "\n"
            << "rho="              << rho              << "\n"
            << "mu="               << mu               << "\n"
            << "d="                << d                << "\n"
            << "l="                << l                << "\n"
            << "mode="             << mode             << "\n"
            << "dir="              << dir              << "\n"
            << "name="             << name             << "\n"
            << "res="              << res              << "\n";
        return oss.str();
    }
};
