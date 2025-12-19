// class CylinderSim(Simulation):
//     def __init__(self,dt,datadt,steps,temp,rho,mu,d,l,dsize,nbins,mode,dir="",res = None,name=''):
//         super().__init__(dt,datadt,steps,temp,rho,mu,d,dsize,nbins,mode,dir,res,name)
//         self.type = 'cyl'
//         self.l = l
//         self.particleData()
//         super().makeFolder()
//     def particleData(self):
//         self.upright = np.array([0,0,self.l/2])
//         self.ior = self.upright
//         self.orientations = np.empty((self.steps+1,3))
//         self.orientations[0] = self.upright.copy()
//         self.vol = np.pi*self.d**2/4*self.l
//         m = self.vol * self.rho # mass of diameter spherical water particle
//         self.mmat = np.eye(3)*m
//         ivert = 0.5*m*(self.d/2)**2 # rotational inertia around the long axis
//         iplanar = self.rho*np.pi*self.l*self.d**2*(self.l**2/3+self.d**2/4)/16 # rotational inertia around the plane axis
//         self.imat = np.diag([iplanar,iplanar,ivert])
//         self.kbTm = kb*self.temp*inv(self.mmat)
//         self.kbTI = kb*self.temp*inv(self.imat)
//         gamma = 2*np.pi*self.mu*self.l
//         zparallel = gamma/((log(self.l/self.d)-0.2)*m) # drag coefficient for movement along the axis ('z' is a stand-in for zeta)
//         znormal = 2*gamma/((log(self.l/self.d)+0.84)*m) #  drag coefficient for movement normal to the axis
//         zrotnormal = gamma*self.l**2/(iplanar*(log(self.l/self.d)-0.66)) # drag coefficient for rotation around the center of the rod lengthwise
//         zaxis = gamma*self.d**2/(2*ivert) # drag coefficient for rotation around the axis
//         self.beta0 = np.diag([znormal,znormal,zparallel])
//         self.betarot0 = np.diag([zrotnormal,zrotnormal,zaxis])
//         self.diffusion = np.trace(self.kbTm @ inv(self.beta0))
//         self.rotdiffusion = np.trace(self.kbTI @ inv(self.betarot0))
//         self.expmsd = self.diffusion*2*self.dt
//         self.expmsad = self.rotdiffusion*2*self.dt
//         self.rmatrix = np.eye(3)
//         self.spring = 0.8 # spring constant
//         self.ligand = np.array([0,0,self.l/2])
//         # self.ligands = np.array([0,0,self.l/2],[0,0,-self.l/2],[0,self.d/2,0],[0,-self.d/2,0],[self.d/2,0,0],[-self.d/2,0,0]) # ligand positions relative to CM
//         self.bonded = [True,False,False,False,False,False] # boolean for if the bond is attached

//         self.equilibrium = self.ligand.copy()
//         self.iligand = self.ligand.copy()
//         self.lambda_eq = 41.1 * 10**-9 # equilibrium bond length
//         self.lambda_t  = self.lambda_eq # bond length, initalized at the equilibrium bond length
//         self.rl = np.array([0,0,self.l/2+self.lambda_eq])
//         self.gamma_bond = 0.274 * 10**-9 # bond reactive compliance 0.274nm
//         self.continueSim = True # boolean for if the bond is attached
//         self.attached = True # boolean for if the bond is attached
//         self.kr0 = 1.1 * 10**-4 # intrinsic bond breakage rate
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
//         self.thvar = self.dt**2 * self.kbTI @ (2*inv(self.brdt) + inv(self.brdt)**2 @ (-3*np.eye(3)+4*e(-self.brdt)-e(-self.brdt)))
//         self.wvar = self.kbTI @ (np.eye(3)-e(-2*self.brdt)) 
//         self.thwcorr = (inv(self.wvar @ self.thvar)**.5) @ self.kbTI @ inv(self.betarot0) @ ((np.eye(3)-e(-self.brdt))**2)
//         receptorlength = 18.7*10**-9 # 18.7nm length of the receptor
//         wallnorm = -np.linalg.norm(self.ligand)
//         wallpoint = self.equilibrium - wallnorm*receptorlength
//         self.wall = wallnorm, wallpoint



//     def __str__(self):
//         return super().__str__()+f'\nlength: {self.l*10**9} nm' + f'\nligand position relative to CM: {self.ligand}' + f'\nequilibrium: {self.equilibrium}' + f'\nspring constant: {self.spring} N/m' 


// simName = '051425'
// simDir = f"Simulations{os.sep}"
// simMode = 'save'
// cylindersim = CylinderSim(dt=10**-9,datadt=10**3, steps=10**6,temp=298.15,rho=997,mu=8.8891*(10**-4),d=5*10**-8,l=5*10**-7,dsize=10,nbins=50,mode=simMode,res = 10**3, name = simName,dir ="Simulations") 
// #spheresim = SphereSim(dt=10**-5,steps= 10**6,temp=298.15,rho=997,mu=8.8891*(10**-4),d=6*10**-8,dsize=10,nbins=50,mode=simMode,res = 10**5, name = simName,dir ="Simulations") 
// #readcylindersim = CylinderSim(dt=10**-5,steps=5*10**5,temp=298.15,rho=997,mu=8.8891*(10**-4),d=5*10**-8,l=5*10**-7,dsize=10,nbins=50,mode='show',res = 10**5, name = simName,dir ="Simulations") 
// #readspheresim = SphereSim(dt=10**-5,steps= 5* 10**5,temp=298.15,rho=997,mu=8.8891*(10**-4),d=6*10**-8,dsize=10,nbins=50,mode='show',res = 10**5, name = simName,dir ="Simulations") 

// #jobSim.mode = 'save'
// # sim = SphereSim(10**-6,10**5,298.15,997,8.9*(10**-4),10**-7,10,30, res= 10**4)sim.generateData()


#include <array>
#include <vector>
#include <string>
#include <optional>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>

static constexpr double kb = 1.380649e-23;
static constexpr double PI = 3.141592653589793238462643383279502884;

using Vec3  = std::array<double, 3>;
using Diag3 = std::array<double, 3>;

static inline double norm3(const Vec3& v) {
  return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline Vec3 vec_sub_scalar(const Vec3& v, double s) {
  return Vec3{v[0]-s, v[1]-s, v[2]-s};
}

static inline Diag3 diag_eye(double a = 1.0) { return Diag3{a,a,a}; }
static inline Diag3 diag_scale(const Diag3& a, double s) { return Diag3{a[0]*s, a[1]*s, a[2]*s}; }
static inline Diag3 diag_add(const Diag3& a, const Diag3& b) { return Diag3{a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
static inline Diag3 diag_sub(const Diag3& a, const Diag3& b) { return Diag3{a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
static inline Diag3 diag_mul(const Diag3& a, const Diag3& b) { return Diag3{a[0]*b[0], a[1]*b[1], a[2]*b[2]}; }

static inline Diag3 diag_inv(const Diag3& a) { return Diag3{1.0/a[0], 1.0/a[1], 1.0/a[2]}; }
static inline Diag3 diag_exp(const Diag3& a) { return Diag3{std::exp(a[0]), std::exp(a[1]), std::exp(a[2])}; }
static inline Diag3 diag_square(const Diag3& a) { return Diag3{a[0]*a[0], a[1]*a[1], a[2]*a[2]}; }
static inline Diag3 diag_sqrt(const Diag3& a) { return Diag3{std::sqrt(a[0]), std::sqrt(a[1]), std::sqrt(a[2])}; }

static inline double diag_trace(const Diag3& a) { return a[0] + a[1] + a[2]; }


class Simulation {
public:
  double dt = 0.0;
  double datadt = 0.0;
  std::size_t steps = 0;
  double temp = 0.0;
  double rho = 0.0;
  double mu  = 0.0;
  double d   = 0.0;
  int dsize  = 0;
  int nbins  = 0;

  std::string mode;
  std::string dir;
  std::optional<double> res;
  std::string name;

  Simulation(double dt_,
             double datadt_,
             std::size_t steps_,
             double temp_,
             double rho_,
             double mu_,
             double d_,
             int dsize_,
             int nbins_,
             std::string mode_,
             std::string dir_ = "",
             std::optional<double> res_ = std::nullopt,
             std::string name_ = "")
  : dt(dt_), datadt(datadt_), steps(steps_), temp(temp_), rho(rho_), mu(mu_), d(d_),
    dsize(dsize_), nbins(nbins_), mode(std::move(mode_)), dir(std::move(dir_)),
    res(res_), name(std::move(name_)) {}

  virtual ~Simulation() = default;

  virtual void makeFolder() {
   // Need to add folder creation 
  }

  virtual std::string toString() const {
    std::ostringstream oss;
    oss << "Simulation"
        << "\n  name: " << name
        << "\n  mode: " << mode
        << "\n  dir:  " << dir
        << "\n  dt:   " << std::setprecision(6) << std::scientific << dt
        << "\n  steps:" << steps
        << "\n  temp: " << temp
        << "\n  rho:  " << rho
        << "\n  mu:   " << mu
        << "\n  d:    " << d;
    return oss.str();
  }
};


class CylinderSim : public Simulation {
public:
  std::string type = "cyl";
  double l = 0.0;

  // geometry/orientation
  Vec3 upright{0,0,0};
  Vec3 ior{0,0,0};
  std::vector<Vec3> orientations;

  // physical properties
  double vol = 0.0;
  Diag3 mmat{0,0,0};
  double ivert = 0.0;
  double iplanar = 0.0;
  Diag3 imat{0,0,0};
  Diag3 kbTm{0,0,0};
  Diag3 kbTI{0,0,0};

  // drag / diffusion
  Diag3 beta0{0,0,0};
  Diag3 betarot0{0,0,0};
  double diffusion = 0.0;
  double rotdiffusion = 0.0;
  double expmsd = 0.0;
  double expmsad = 0.0;

  // rotation matrix (full 3x3 identity)
  std::array<std::array<double,3>,3> rmatrix{{ {1,0,0}, {0,1,0}, {0,0,1} }};

  // bond model
  double spring = 0.8;
  Vec3 ligand{0,0,0};
  std::array<bool,6> bonded{{true,false,false,false,false,false}};

  Vec3 equilibrium{0,0,0};
  Vec3 iligand{0,0,0};
  double lambda_eq = 41.1e-9;
  double lambda_t  = 41.1e-9;
  Vec3 rl{0,0,0};
  double gamma_bond = 0.274e-9;
  bool continueSim = true;
  bool attached = true;
  double kr0 = 1.1e-4;

  // translational kinematics
  Diag3 bdt{0,0,0};
  Diag3 rvar{0,0,0};
  Diag3 vvar{0,0,0};
  Diag3 rvcorr{0,0,0};

  // rotational kinematics
  Diag3 brdt{0,0,0};
  Diag3 thvar{0,0,0};
  Diag3 wvar{0,0,0};
  Diag3 thwcorr{0,0,0};

  // wall
  double wallnorm = 0.0;
  Vec3 wallpoint{0,0,0};

  CylinderSim(double dt_,
              double datadt_,
              std::size_t steps_,
              double temp_,
              double rho_,
              double mu_,
              double d_,
              double l_,
              int dsize_,
              int nbins_,
              std::string mode_,
              std::string dir_ = "",
              std::optional<double> res_ = std::nullopt,
              std::string name_ = "")
  : Simulation(dt_, datadt_, steps_, temp_, rho_, mu_, d_, dsize_, nbins_, std::move(mode_), std::move(dir_), res_, std::move(name_)),
    l(l_)
  {
    type = "cyl";
    particleData();
    makeFolder();
  }

  void particleData() {
    // upright/orientation history
    upright = Vec3{0.0, 0.0, l/2.0};
    ior = upright;
    orientations.assign(steps + 1, Vec3{0.0,0.0,0.0});
    orientations[0] = upright;

    // volume and mass
    vol = PI * (d*d) / 4.0 * l;
    const double m = vol * rho;
    mmat = Diag3{m,m,m};

    // inertia
    ivert   = 0.5 * m * std::pow(d/2.0, 2.0);
    iplanar = rho * PI * l * (d*d) * ( (l*l)/3.0 + (d*d)/4.0 ) / 16.0;
    imat = Diag3{iplanar, iplanar, ivert};

    // kb*T * inv(mmat/imat)
    kbTm = diag_scale(diag_inv(mmat), kb * temp);
    kbTI = diag_scale(diag_inv(imat), kb * temp);

    // drag coefficients
    const double gamma = 2.0 * PI * mu * l;
    const double log_ld = std::log(l/d);

    const double zparallel   = gamma / ((log_ld - 0.2)  * m);
    const double znormal     = 2.0 * gamma / ((log_ld + 0.84) * m);
    const double zrotnormal  = gamma * (l*l) / (iplanar * (log_ld - 0.66));
    const double zaxis       = gamma * (d*d) / (2.0 * ivert);

    beta0     = Diag3{znormal, znormal, zparallel};
    betarot0  = Diag3{zrotnormal, zrotnormal, zaxis};

    diffusion     = diag_trace(diag_mul(kbTm, diag_inv(beta0)));
    rotdiffusion  = diag_trace(diag_mul(kbTI, diag_inv(betarot0)));

    expmsd  = diffusion * 2.0 * dt;
    expmsad = rotdiffusion * 2.0 * dt;

    // identity rotation matrix 
    spring = 0.8;

    // ligand & bond
    ligand = Vec3{0.0, 0.0, l/2.0};
    bonded = {{true,false,false,false,false,false}};

    equilibrium = ligand;
    iligand     = ligand;
    lambda_eq   = 41.1e-9;
    lambda_t    = lambda_eq;
    rl          = Vec3{0.0, 0.0, l/2.0 + lambda_eq};
    gamma_bond  = 0.274e-9;
    continueSim = true;
    attached    = true;
    kr0         = 1.1e-4;

    // translational kinematics 
    bdt = diag_scale(beta0, dt);

    // rvar = dt^2 * kbTm @ ( 2*inv(bdt) + inv(bdt)^2 @ (-3*I + 4*e(-bdt) - e(-2*bdt)) )
    const Diag3 inv_bdt  = diag_inv(bdt);
    const Diag3 inv_bdt2 = diag_square(inv_bdt);

    const Diag3 I = diag_eye(1.0);
    const Diag3 bracket = diag_add(
      diag_scale(I, -3.0),
      diag_sub(
        diag_scale(diag_exp(diag_scale(bdt, -1.0)), 4.0),
        diag_exp(diag_scale(bdt, -2.0))
      )
    );

    const Diag3 term_inside = diag_add(
      diag_scale(inv_bdt, 2.0),
      diag_mul(inv_bdt2, bracket)
    );

    rvar = diag_scale(diag_mul(kbTm, term_inside), dt*dt);

    // vvar = kbTm @ (I - e(-2*bdt))
    vvar = diag_mul(kbTm, diag_sub(I, diag_exp(diag_scale(bdt, -2.0))));

    // rvcorr = (inv(vvar @ rvar)**.5) @ kbTm @ inv(beta0) @ ((I - e(-bdt))**2)
    const Diag3 vvar_rvar = diag_mul(vvar, rvar);
    const Diag3 sqrt_inv_vr = diag_sqrt(diag_inv(vvar_rvar));
    const Diag3 I_minus_exp = diag_sub(I, diag_exp(diag_scale(bdt, -1.0)));
    const Diag3 I_minus_exp_sq = diag_square(I_minus_exp);

    rvcorr = diag_mul(
      sqrt_inv_vr,
      diag_mul(kbTm, diag_mul(diag_inv(beta0), I_minus_exp_sq))
    );

    // rotational kinematics
    brdt = diag_scale(betarot0, dt);

    // thvar = dt^2 * kbTI @ ( 2*inv(brdt) + inv(brdt)^2 @ (-3*I + 4*e(-brdt) - e(-brdt)) )
    const Diag3 inv_brdt  = diag_inv(brdt);
    const Diag3 inv_brdt2 = diag_square(inv_brdt);

    const Diag3 bracket_rot = diag_add(
      diag_scale(I, -3.0),
      diag_add(
        diag_scale(diag_exp(diag_scale(brdt, -1.0)), 4.0),
        diag_scale(diag_exp(diag_scale(brdt, -1.0)), -1.0)
      )
    );

    const Diag3 term_inside_rot = diag_add(
      diag_scale(inv_brdt, 2.0),
      diag_mul(inv_brdt2, bracket_rot)
    );

    thvar = diag_scale(diag_mul(kbTI, term_inside_rot), dt*dt);

    // wvar = kbTI @ (I - e(-2*brdt))
    wvar = diag_mul(kbTI, diag_sub(I, diag_exp(diag_scale(brdt, -2.0))));

    // thwcorr = (inv(wvar @ thvar)**.5) @ kbTI @ inv(betarot0) @ ((I - e(-brdt))**2)
    const Diag3 wvar_thvar = diag_mul(wvar, thvar);
    const Diag3 sqrt_inv_wt = diag_sqrt(diag_inv(wvar_thvar));
    const Diag3 I_minus_exp_rot = diag_sub(I, diag_exp(diag_scale(brdt, -1.0)));
    const Diag3 I_minus_exp_rot_sq = diag_square(I_minus_exp_rot);

    thwcorr = diag_mul(
      sqrt_inv_wt,
      diag_mul(kbTI, diag_mul(diag_inv(betarot0), I_minus_exp_rot_sq))
    );

    // wall
    const double receptorlength = 18.7e-9;
    wallnorm  = -norm3(ligand);
    wallpoint = vec_sub_scalar(equilibrium, wallnorm * receptorlength); // mirrors Python broadcasting
  }

  std::string toString() const override {
    std::ostringstream oss;
    oss << Simulation::toString()
        << "\nlength: " << (l * 1e9) << " nm"
        << "\nligand position relative to CM: [" << ligand[0] << ", " << ligand[1] << ", " << ligand[2] << "]"
        << "\nequilibrium: [" << equilibrium[0] << ", " << equilibrium[1] << ", " << equilibrium[2] << "]"
        << "\nspring constant: " << spring << " N/m";
    return oss.str();
  }
};








