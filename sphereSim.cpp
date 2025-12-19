// class SphereSim(Simulation):
//     def __init__(self,dt,datadt,steps,temp,rho,mu,d,dsize,nbins,mode,dir="",res = None, name = ''):
//         super().__init__(dt,datadt,steps,temp,rho,mu,d,dsize,nbins,mode,dir,res,name)
//         self.type = 'sph'
//         self.particleData()
//         super().makeFolder()
//     def next_data(self):
//         # translational kinematics
//         randomr = np.random.normal(loc=0,scale=1,size=3) @ self.rvar**.5
//         #print(f"randomr1:{randomr}")
//         randomv = np.random.normal(loc=0, scale=1, size=3) @ (self.vvar@(np.eye(3)-self.rvcorr**2))**.5 + randomr @ self.rvcorr @ (self.vvar @ inv(self.rvar)) ** .5
//         #print(f"randomv:{randomv}")
//         randomr = np.eye(3) @ randomr
//         #print(f"randomr2:{randomr}")
//         #print(f"rmatrix:{self.rmatrix}")
//         randomr = self.rmatrix @ randomr
//         #print(f"randomr3:{randomr}")
//         randomv = self.rmatrix @ randomv
//         #print(f"randomv2:{randomv}")
//         dr = randomr + (np.eye(3)-e(-self.bdt)) @ inv(self.beta0) @ self.iv
//         #print(dr)
//         self.sqdr = dr @ dr
//         self.ir += dr
//         self.iv = e(-self.bdt) @ self.iv + randomv
        
//         # rotational kinematics 
//         randomth = np.random.normal(loc=0,scale=1,size=3) @ self.thvar**.5
//         #print(f"randomth2:{randomth}")
//         randomw = np.random.normal(loc=0, scale=1, size=3) @ (self.wvar@(np.eye(3)-self.thwcorr**2))**.5 + randomth @ self.thwcorr @ (self.wvar @ inv(self.thvar)) ** .5
//         randomth = self.rmatrix @ randomth
//         #print(f"randomth2:{randomth}")
//         randomw = self.rmatrix @ randomw
//         dth = randomth + (np.eye(3)-e(-self.brdt)) @ inv(self.betarot0) @ self.iw
//         #print(f"dth-randomth:{dth-randomth}")
//         #print(f"dth:{dth}")
//         self.sqdth = dth @ dth
//         self.iw = e(-self.brdt) @ self.iw + randomw
//         self.rmatrix = rotateMatrix(dth,self.rmatrix)
//         self.ior = self.rmatrix @ self.upright
//         #self.ith = [theta(self.ior[0],self.ior[1]),self.ior[2]/np.sqrt(self.upright@self.upright)]

//     def particleData(self):
//         self.upright = np.array([0,0,self.d/2])
//         self.ior = self.upright
//         self.orientations = np.empty((self.steps+1,3))
//         self.orientations[0] = self.upright.copy()

//         self.vol = np.pi*self.d**3/6
//         m = self.vol*self.rho
//         self.mmat = m*np.eye(3)
//         self.imat = m *self.d**2/10 * np.eye(3)
//         self.kbTm = kb*self.temp*inv(self.mmat)
//         self.kbTI = kb*self.temp*inv(self.imat)
//         self.beta0 = inv(self.mmat)*3*self.d*self.mu*np.pi
//         self.betarot0 = inv(self.imat)*self.d**3*np.pi*self.mu
//         self.diffusion = np.trace(self.kbTm @ inv(self.beta0))
//         self.rotdiffusion = np.trace(self.kbTI @ inv(self.betarot0))
//         self.expmsd = self.diffusion*2*self.dt
//         self.expmsad = self.rotdiffusion*2*self.dt
//         self.rmatrix = np.eye(3)
//         self.spring = 1 # spring constant
//         #translational kinematics data
//         self.bdt = self.beta0*self.dt
//         #print(self.bdt)
//         self.rvar = self.dt**2 * self.kbTm @ (2*inv(self.bdt) + inv(self.bdt)**2 @ (-3*np.eye(3)+4*e(-self.bdt)-e(-2*self.bdt)))
//         #print(self.rvar)
//         self.vvar = self.kbTm @ (np.eye(3)-e(-2*self.bdt))
//         #print(self.vvar)
//         self.rvcorr =  (inv(self.vvar @ self.rvar)**.5) @ self.kbTm @ inv(self.beta0) @ ((np.eye(3)-e(-1*self.bdt))**2)

//         #rotational kinematics data
//         self.brdt = self.betarot0 * self.dt 
//         self.thvar = self.dt**2 * self.kbTI @ (2*inv(self.brdt) + inv(self.brdt)**2 @ (-3*np.eye(3)+4*e(-self.brdt)-e(-2*self.brdt)))
//         self.wvar = self.kbTI @ (np.eye(3)-e(-2*self.brdt)) 
//         self.thwcorr = (inv(self.wvar @ self.thvar)**.5) @ self.kbTI @ inv(self.betarot0) @ ((np.eye(3)-e(-1*self.brdt))**2)



#include <Eigen/Dense>
#include <random>
#include <vector>
#include <cmath>

constexpr double kb = 1.380649e-23; // Boltzmann constant (J/K)


namespace {
    std::mt19937 &rng() {
        static std::mt19937 gen(std::random_device{}());
        return gen;
    }

    double normal01() {
        static std::normal_distribution<double> dist(0.0, 1.0);
        return dist(rng());
    }

    Eigen::Vector3d randn3() {
        return Eigen::Vector3d(normal01(), normal01(), normal01());
    }

    // Rodrigues rotation: rotate existing rotation matrix by "rotation vector" dth
    Eigen::Matrix3d rotateMatrix(const Eigen::Vector3d &dth,
                                 const Eigen::Matrix3d &rmatrix) {
        double theta = dth.norm();
        if (theta == 0.0) {
            return rmatrix;
        }
        Eigen::Vector3d u = dth / theta;
        double ux = u.x();
        double uy = u.y();
        double uz = u.z();
        double c = std::cos(theta);
        double s = std::sin(theta);
        double one_c = 1.0 - c;

        Eigen::Matrix3d R;
        R << c + ux * ux * one_c,      ux * uy * one_c - uz * s, ux * uz * one_c + uy * s,
             uy * ux * one_c + uz * s, c + uy * uy * one_c,      uy * uz * one_c - ux * s,
             uz * ux * one_c - uy * s, uz * uy * one_c + ux * s, c + uz * uz * one_c;

        // Apply incremental rotation on the left (same idea as your Python code)
        return R * rmatrix;
    }
} 

// class: SphereSim
class SphereSim {
public:
    double dt;
    int steps;
    double temp, rho, mu, d;

    // state
    Eigen::Vector3d ir; // position
    Eigen::Vector3d iv; // translational velocity
    Eigen::Vector3d iw; // angular velocity
    double sqdr{};      // squared displacement
    double sqdth{};     // squared angular displacement

    Eigen::Vector3d upright; // reference vector in body frame
    Eigen::Vector3d ior;     // orientation of upright in lab frame
    std::vector<Eigen::Vector3d> orientations;

    double vol{};
    double diffusion{};
    double rotdiffusion{};
    double expmsd{};
    double expmsad{};
    double spring{1.0};

    Eigen::Matrix3d rmatrix; // rotation matrix body->lab

    // Diagonal matrices stored as their 3 diagonal entries
    Eigen::Vector3d mmat_diag;
    Eigen::Vector3d imat_diag;
    Eigen::Vector3d kbTm_diag;
    Eigen::Vector3d kbTI_diag;
    Eigen::Vector3d beta0_diag;
    Eigen::Vector3d betarot0_diag;
    Eigen::Vector3d bdt_diag;
    Eigen::Vector3d brdt_diag;
    Eigen::Vector3d rvar_diag;
    Eigen::Vector3d vvar_diag;
    Eigen::Vector3d rvcorr_diag;
    Eigen::Vector3d thvar_diag;
    Eigen::Vector3d wvar_diag;
    Eigen::Vector3d thwcorr_diag;

    // constructor
    SphereSim(double dt_, int steps_, double temp_,
              double rho_, double mu_, double d_)
        : dt(dt_), steps(steps_), temp(temp_), rho(rho_), mu(mu_), d(d_) {
        ir.setZero();
        iv.setZero();
        iw.setZero();
        particleData();
    }

    void next_data();
    void particleData();
};

// particleData 
void SphereSim::particleData() {
    using Eigen::Array3d;
    constexpr double PI = 3.14159265358979323846;

    // self.upright = np.array([0,0,self.d/2])
    upright = Eigen::Vector3d(0.0, 0.0, d * 0.5);
    ior = upright;

    // self.orientations = np.empty((self.steps+1,3))
    orientations.assign(steps + 1, Eigen::Vector3d::Zero());
    if (!orientations.empty()) {
        orientations[0] = upright;
    }

    // self.vol = np.pi*self.d**3/6
    // m = self.vol*self.rho
    vol = PI * d * d * d / 6.0;
    double m = vol * rho;

    // self.mmat = m*np.eye(3)
    // self.imat = m *self.d**2/10 * np.eye(3)
    mmat_diag = Eigen::Vector3d::Constant(m);
    imat_diag = Eigen::Vector3d::Constant(m * d * d / 10.0);

    // self.kbTm = kb*self.temp*inv(self.mmat)
    // self.kbTI = kb*self.temp*inv(self.imat)
    kbTm_diag = Eigen::Vector3d::Constant(kb * temp).cwiseQuotient(mmat_diag);
    kbTI_diag = Eigen::Vector3d::Constant(kb * temp).cwiseQuotient(imat_diag);

    // self.beta0 = inv(self.mmat)*3*self.d*self.mu*np.pi
    // self.betarot0 = inv(self.imat)*self.d**3*np.pi*self.mu
    beta0_diag    = mmat_diag.cwiseInverse() * (3.0 * d * mu * PI);
    betarot0_diag = imat_diag.cwiseInverse() * (d * d * d * PI * mu);

    // self.diffusion = np.trace(self.kbTm @ inv(self.beta0))
    // self.rotdiffusion = np.trace(self.kbTI @ inv(self.betarot0))
    diffusion    = (kbTm_diag.array() / beta0_diag.array()).sum();
    rotdiffusion = (kbTI_diag.array() / betarot0_diag.array()).sum();

    // self.expmsd = self.diffusion*2*self.dt
    // self.expmsad = self.rotdiffusion*2*self.dt
    expmsd  = diffusion * 2.0 * dt;
    expmsad = rotdiffusion * 2.0 * dt;

    // self.rmatrix = np.eye(3)
    rmatrix = Eigen::Matrix3d::Identity();
    spring  = 1.0; // spring constant

    // translational kinematics data
    // self.bdt = self.beta0*self.dt
    bdt_diag = beta0_diag * dt;

    Array3d bdt          = bdt_diag.array();
    Array3d inv_bdt      = 1.0 / bdt;
    Array3d exp_neg_bdt  = (-bdt).exp();
    Array3d exp_neg_2bdt = (-2.0 * bdt).exp();

    // self.rvar = dt^2 * kbTm @ (2*inv(bdt) + inv(bdt)^2 @ (-3I + 4e(-bdt) - e(-2bdt)))
    Array3d X      = -3.0 + 4.0 * exp_neg_bdt - exp_neg_2bdt;
    Array3d inside = 2.0 * inv_bdt + inv_bdt.square() * X;
    rvar_diag      = (dt * dt) * (kbTm_diag.array() * inside).matrix();

    // self.vvar = kbTm @ (I - e(-2bdt))
    vvar_diag = (kbTm_diag.array() * (1.0 - exp_neg_2bdt)).matrix();

    // self.rvcorr = (inv(vvar @ rvar)**.5) @ kbTm @ inv(beta0) @ ((I - e(-bdt))**2)
    Array3d inv_sqrt_vr        = 1.0 / (vvar_diag.array() * rvar_diag.array()).sqrt();
    Array3d kbTm_over_beta0    = kbTm_diag.array() / beta0_diag.array();
    Array3d one_minus_exp_bdt  = 1.0 - exp_neg_bdt;
    rvcorr_diag = (inv_sqrt_vr * kbTm_over_beta0 * one_minus_exp_bdt.square()).matrix();

    //  rotational kinematics data 
    // self.brdt = self.betarot0 * self.dt
    brdt_diag = betarot0_diag * dt;

    Array3d brdt          = brdt_diag.array();
    Array3d inv_brdt      = 1.0 / brdt;
    Array3d exp_neg_brdt  = (-brdt).exp();
    Array3d exp_neg_2brdt = (-2.0 * brdt).exp();

    // self.thvar = dt^2 * kbTI @ (2*inv(brdt) + inv(brdt)^2 @ (-3I + 4e(-brdt) - e(-2brdt)))
    Array3d Xr      = -3.0 + 4.0 * exp_neg_brdt - exp_neg_2brdt;
    Array3d insider = 2.0 * inv_brdt + inv_brdt.square() * Xr;
    thvar_diag      = (dt * dt) * (kbTI_diag.array() * insider).matrix();

    // self.wvar = kbTI @ (I - e(-2brdt))
    wvar_diag = (kbTI_diag.array() * (1.0 - exp_neg_2brdt)).matrix();

    // self.thwcorr = (inv(wvar @ thvar)**.5) @ kbTI @ inv(betarot0) @ ((I - e(-brdt))**2)
    Array3d inv_sqrt_wt        = 1.0 / (wvar_diag.array() * thvar_diag.array()).sqrt();
    Array3d kbTI_over_betarot0 = kbTI_diag.array() / betarot0_diag.array();
    Array3d one_minus_exp_brdt = 1.0 - exp_neg_brdt;
    thwcorr_diag = (inv_sqrt_wt * kbTI_over_betarot0 * one_minus_exp_brdt.square()).matrix();
}

// next_data
void SphereSim::next_data() {
    using Eigen::Array3d;

    //  translational kinematics 
    // randomr = N(0,1,3) @ rvar**0.5   (diag => elementwise)
    Eigen::Vector3d rand_r  = randn3();
    Eigen::Vector3d randomr = (rand_r.array() * rvar_diag.array().sqrt()).matrix();

    // randomv Gaussian term
    // first part: N * (vvar @ (I - rvcorr**2))**0.5
    Eigen::Vector3d rand_v1 = randn3();
    Array3d vterm1_sqrt     = (vvar_diag.array() * (1.0 - rvcorr_diag.array().square())).sqrt();
    Eigen::Vector3d part1   = (rand_v1.array() * vterm1_sqrt).matrix();

    // second part: randomr @ rvcorr @ (vvar @ inv(rvar))**0.5
    Array3d ratio_vr        = (vvar_diag.array() / rvar_diag.array()).sqrt();
    Eigen::Vector3d part2   = (randomr.array() * rvcorr_diag.array() * ratio_vr).matrix();

    Eigen::Vector3d randomv = part1 + part2;

    // randomr = self.rmatrix @ randomr
    // randomv = self.rmatrix @ randomv
    randomr = rmatrix * randomr;
    randomv = rmatrix * randomv;

    // dr = randomr + (I - e(-bdt)) @ inv(beta0) @ iv
    Array3d exp_neg_bdt = (-bdt_diag.array()).exp();
    Eigen::Vector3d drift =
        ((1.0 - exp_neg_bdt) * (iv.array() / beta0_diag.array())).matrix();

    Eigen::Vector3d dr = randomr + drift;

    // self.sqdr = dr @ dr
    sqdr = dr.squaredNorm();

    // self.ir += dr
    ir += dr;

    // self.iv = e(-bdt) @ iv + randomv
    iv = (exp_neg_bdt * iv.array()).matrix() + randomv;

    //  rotational kinematics 
    // randomth = N(0,1,3) @ thvar**0.5
    Eigen::Vector3d rand_th  = randn3();
    Eigen::Vector3d randomth = (rand_th.array() * thvar_diag.array().sqrt()).matrix();

    // randomw
    Eigen::Vector3d rand_w1  = randn3();
    Array3d wterm1_sqrt      = (wvar_diag.array() * (1.0 - thwcorr_diag.array().square())).sqrt();
    Eigen::Vector3d rpart1   = (rand_w1.array() * wterm1_sqrt).matrix();

    Array3d ratio_wt         = (wvar_diag.array() / thvar_diag.array()).sqrt();
    Eigen::Vector3d rpart2   = (randomth.array() * thwcorr_diag.array() * ratio_wt).matrix();

    Eigen::Vector3d randomw  = rpart1 + rpart2;

    // randomth = rmatrix @ randomth
    // randomw = rmatrix @ randomw
    randomth = rmatrix * randomth;
    randomw  = rmatrix * randomw;

    // dth = randomth + (I - e(-brdt)) @ inv(betarot0) @ iw
    Array3d exp_neg_brdt = (-brdt_diag.array()).exp();
    Eigen::Vector3d rot_drift =
        ((1.0 - exp_neg_brdt) * (iw.array() / betarot0_diag.array())).matrix();

    Eigen::Vector3d dth = randomth + rot_drift;

    // self.sqdth = dth @ dth
    sqdth = dth.squaredNorm();

    // self.iw = e(-brdt) @ iw + randomw
    iw = (exp_neg_brdt * iw.array()).matrix() + randomw;

    // self.rmatrix = rotateMatrix(dth, self.rmatrix)
    // self.ior = self.rmatrix @ self.upright
    rmatrix = rotateMatrix(dth, rmatrix);
    ior     = rmatrix * upright;
}
