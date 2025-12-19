// Convert to cpp
// class Simulation:
//     def __init__(self,dt,datadt,steps,temp,rho,mu,d,dsize,nbins, mode= 'show',dir="",res=None,name =''):
//         self.type = ''
//         self.dt = dt # s
//         self.datadt = datadt # int representing the number of steps between data points
//         self.steps = int(steps/datadt)
//         self.laststep = self.steps
//         self.totalt = dt*steps # s
//         self.temp = temp # K
//         self.rho = rho # density, kg/m^3
//         self.mu = mu # viscosity, N*s/m^2
//         self.d = d # diameter, m
//         self.l = d 
//         self.ir = np.array([0.,0.,0.]) # position vector
//         self.iv = np.array([0.,0.,0.]) # velocity vector
//         self.ith = np.array([0.,0.]) # theta, z
//         self.iw = np.array([0.,0.,0.]) # angular velocity vector
        
//         self.positions = np.empty((self.steps+1,3))
//         self.positions[0] = self.ir.copy()

//         self.velocities = np.empty((self.steps+1,3))
//         self.velocities[0] = self.iv.copy()

//         self.thetaz = np.empty((self.steps+1,2))
//         self.thetaz[0] = self.ith.copy()  # initial position assumed the axis of the cylinder is aligned with the z axis

//         self.orientations = np.empty((self.steps+1,3))
//         self.orientations[0] = np.array([1.,0.,0.])
        
//         self.omegas = np.empty((self.steps+1,3))
//         self.omegasobjframe = np.empty((self.steps+1,3))
//         self.omegas[0] = self.iw.copy()
//         self.omegasobjframe[0] = np.array([0.,0.,0.])

        
//         self.sqdis = np.empty(self.steps+1) # squared displacements
//         self.sqangdis = np.empty(self.steps+1) # squared angular displacements
//         self.sqdis[0] = 0
//         self.sqangdis[0] = 0
    
//         self.squared_displacement = [self.sqdis.copy(),self.sqangdis.copy()]
//         self.times = np.array(range(0,self.steps+1)) * (dt * datadt)
//         self.dsize = dsize
//         self.nbins = nbins
//         self.mode = mode # either 'show' or 'save'
//         self.dir = dir
//         self.name = name
//         if (not(self.dir == "")):
//             self.dir += '/'
//         if res == None:
//             self.res = steps
//         else: 
//             self.res = res


//     def makeFolder(self):
//         if (self.mode == "save"):
//             i = 1
//             while(True):
//                 try:
//                     Path(self.dir + self.name + self.type).mkdir()
//                 except FileExistsError:
//                     i += 1
//                     if (i==2):
//                         self.name += 'v2'
//                     else:
//                         self.name = self.name[0:-1:] + str(i)
//                 else:
//                     break
            

//             with open(self.dir + self.name + self.type  + "/parameters.txt",'w') as file:
//                 file.write(str(self)) 
//         self.name += self.type      
//     def particleData(self):
//         self.upright = np.array([1,0,0])
//         self.ior = self.upright
//         self.orientations = np.empty((self.steps+1,3))
//         self.orientations[0] = self.upright.copy()
//         self.vol = 0
//         self.mmat = np.eye(3)
//         self.imat = np.eye(3)
//         self.kbTm = np.eye(3)
//         self.kbTI = np.eye(3)
//         self.beta0 = np.eye(3)
//         self.betarot0 = np.eye(3)
//         self.diffusion = 0
//         self.rotdiffusion = 0
//         self.expmsd = self.diffusion*2*self.dt
//         self.expmsad = self.rotdiffusion*2*self.dt
//         self.rmatrix = np.eye(3)
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
//         self.spring = 1
//         self.ligand = self.upright 
//         self.iligand = self.ligand.copy() # instantaneous ligand position relative to the CM
//         self.equilibrium = self.ligand # equilibrium position (aka the position of the ligand when the bond is formed)
        
//     def next_data(self):
//         #bondforce = (self.equilibrium - self.ir - self.iligand) * self.spring # position vector + arm vector - vector of equilibrium position (upright)
        
//         if self.continueSim:
//             rr = self.ir + self.iligand # position of ligand (m)
//             ilambda = np.linalg.norm(self.rl - rr) 
//             bondforce = self.spring * (1 - self.lambda_eq/ilambda) * (self.rl - rr) # force on the ligand due to the spring
//             delta = np.abs(self.lambda_eq - ilambda) # distance from the equilibrium position of the ligand
//             kr = self.kr0 * np.exp(self.gamma_bond * self.spring * delta/(kb*self.temp)) # spring constant of the bond
//             p_r = 1-np.exp(-kr*self.dt) # probability of breaking the bond
//             print(p_r)
            
//             if (np.random.rand() < p_r): # break the bond with probability p_r
//                 self.continueSim = False
//         else:
//             bondforce = np.array([0,0,0])
//         #print(f"bond force: {bondforce}")
//         #print(f"tip position:{self.ir+self.ior}")
//         #n = 2*self.ior/self.l # unit normal vector of the circular tip of the cylinder
//         #zpos = np.array([0,0,1]) # z axis, pointing towards the "wall"
//         #corner0 = zpos - zpos @ n * n # the corner of the cylinder closest to the wall
//         #corner0 = self.d*corner0/(2*np.linalg.norm(corner0)) + self.ir + self.ior # position of the corner

//         # translational kinematics
//         c0 = e(-self.bdt)
//         c1 = (np.eye(3)-c0) @ inv(self.beta0)
//         c2 = (self.dt*np.eye(3)-c1) @ inv(self.beta0)
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
        
//         k =  inv(self.mmat) @ bondforce
//         #print(f"k:{k}")
//         #print(f"randomv2:{randomv}")
//         dr = randomr + self.rmatrix @ c1 @ self.rmatrix.T @ self.iv + self.rmatrix @ c2 @ self.rmatrix.T @ k
//         #print(f"dr:{dr}")
//         self.sqdr = dr @ dr
//         self.ir += dr
//         self.iv = (self.rmatrix @ c0 @ self.rmatrix.T) @ self.iv + randomv + self.rmatrix @ c1 @ self.rmatrix.T @ k 

//         # rotational kinematics 
//         c0 = e(-self.brdt)
//         c1 = (np.eye(3)-c0) @ inv(self.betarot0)
//         c2 = (self.dt*np.eye(3)-c1) @ inv(self.betarot0)

//         #print(f"c2:{c2}")

//         torque = np.cross(self.iligand,bondforce) 
//         #print(f"torque: {torque}")
//         randomth = np.random.normal(loc=0,scale=1,size=3) @ self.thvar**.5
//         #print(f"randomth2:{randomth}")
//         randomw = np.random.normal(loc=0, scale=1, size=3) @ (self.wvar@(np.eye(3)-self.thwcorr**2))**.5 + randomth @ self.thwcorr @ (self.wvar @ inv(self.thvar)) ** .5
//         randomth = self.rmatrix @ randomth
//         #print(f"randomth2:{randomth}")
//         randomw = self.rmatrix @ randomw
//         dth = randomth + self.rmatrix @ c1 @ self.rmatrix.T @ self.iw +  self.rmatrix @ c2 @ inv(self.imat) @ self.rmatrix.T @ torque
//         #print(f"dth:{dth}")
//         #print(f"dth-randomth:{dth-randomth}")
//         #print(f"dth:{dth}")
//         #self.sqdth = dth @ dth
//         self.iw = (self.rmatrix @ c0 @ self.rmatrix.T) @ self.iw + randomw + self.rmatrix @ c1 @ inv(self.imat) @ self.rmatrix.T @ torque
//         self.rmatrix = rotateMatrix(dth,self.rmatrix)
//         self.iligand = self.rmatrix @ self.ligand
//         self.ior = self.rmatrix @ self.upright
//         self.ith = [theta(self.ior[0],self.ior[1]),self.ior[2]/np.sqrt(self.upright@self.upright)]

//         #n = 2*self.ior/self.l # checks corner position after movement
//         #cornerf = zpos - zpos @ n * n #
//         #cornerf = self.d*corner0/(2*np.linalg.norm(corner0)) + self.ir + self.ior
//         #collision = (cornerf - self.wall(1)) @ self.wall(0) < 0 # checks if the corner of the cylinder has gone past the wall


//     def generateData(self):
//         for i in range(1,self.steps*self.datadt+1):
            
//             #print(f"step {i}")
//             #print("position")
//             self.next_data()
//             if (not(self.continueSim)):
//                 self.laststep = i
//                 self.updateProgress(f"Stopped at {i*self.dt} s")
//                 break
//             #self.ir = (self.ir+graphboxsize/2) % graphboxsize - graphboxsize/2
//             # iomobj = self.rmatrix.T@self.iw
//             if (i%self.datadt == 0):
//                 idx = int(i/self.datadt)
//                 self.positions[idx] = self.ir.copy()
//                 self.velocities[idx] = self.iv.copy()
//                 self.thetaz[idx] = self.ith.copy()
//                 self.omegas[idx] = self.iw.copy()
//                 #self.omegasobjframe.append(iomobj)
//                 #self.sqdis.append(self.sqdr.copy())
//                 #self.sqangdis.append(self.sqdth.copy())
//                 self.orientations[idx] = self.rmatrix@self.upright
//                 #self.absoluteorientations = [(rotateMatrix(self.ith)@np.array([0,0,1])) @ (rotateMatrix(self.ith)@np.array([0,0,1]))]
//                 #print(self.absoluteorientations)
//                 # print('test-1')
//             if ((i)%(self.steps*self.datadt*10**-3) == 0):
//                 with open(self.dir + self.name + "/progress.txt",'w') as file:
//                     file.write(f"iteration {i}")
//           # magnitudes = {(ornt @ ornt)/(self.upright@self.upright) for ornt in self.orientations}
//           # print(f"\n{min(magnitudes)},{max(magnitudes)}")
//           self.sqdis= np.array(self.sqdis[0:self.laststep+1])
//           self.sqangdis = np.array(self.sqangdis[0:self.laststep+1])
//           #print(self.positions)
//           self.positions = pd.DataFrame(self.positions[0:self.laststep+1] , columns= ['x','y','z'])
//           # print(positions)
//           self.velocities = pd.DataFrame(self.velocities[0:self.laststep+1], columns= ['x','y','z'])
//           self.thetaz = pd.DataFrame(self.thetaz[0:self.laststep+1], columns=['theta','z'])
//           # print(thetas)
//           self.omegas = pd.DataFrame(self.omegas[0:self.laststep+1], columns=['roll','pitch','yaw'])
//           self.omegasobjframe = pd.DataFrame(self.omegasobjframe[0:self.laststep+1], columns = ['object x-axis','object y-axis','object z-axis'])
//           self.orientations = pd.DataFrame(self.orientations[0:self.laststep+1], columns = ['x','y','z'])
//           self.squared_displacement =  pd.DataFrame({'Translational': self.sqdis[0:self.laststep+1],'Angular': self.sqangdis[0:self.laststep+1]})
//           self.orientations *= 10**9
//           self.positions*= 10**9
//           self.times *= 10**9
//           self.totalt *= 10**9
//           self.dt *= 10**9 

    
//           #elif self.mode == 'save':
//               #pass
//               #plotlyfig = tls.mpl_to_plotly(self.fig)
//               #pio.write_html(plotlyfig, file = self.dir + "3d_orientations.html")


//               def saveSim(self):
//           self.datadict = {
//               'orientations': self.orientations,
//               'positions': self.positions,
//               'thetaz': self.thetaz,
//               'velocities': self.velocities,
//               'omegas': self.omegas,
//               'omegasobjframe': self.omegasobjframe,
//               'squared_displacement': self.squared_displacement
//           }
  
//           for file_name, data_object in self.datadict.items():
//               # Save to CSV
//               data_object.to_hdf(self.dir + self.name + '/' +f'{self.name}.h5', key = file_name, mode = 'a', format = 'table', complib='zlib', complevel=9)
//               # Generate a progress message based on the file name
//               progress_message = f"{file_name} saved"
//               self.updateProgress(progress_message)

                  
#include <array>
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

class Simulation {
public:

  using Vec3 = std::array<double, 3>;
  using Vec2 = std::array<double, 2>;
  using Mat3 = std::array<std::array<double, 3>, 3>;


  Simulation(double dt_,
             int datadt_,
             int steps_,         // "raw" steps (like Python's steps before dividing by datadt)
             double temp_,
             double rho_,
             double mu_,
             double d_,
             double dsize_,
             int nbins_,
             std::string mode_ = "show",
             std::string dir_ = "",
             int res_ = -1,       // -1 means "None" => res = steps
             std::string name_ = "")
  : type(""),
    dt(dt_),
    datadt(datadt_),
    total_raw_steps(steps_),
    steps(std::max(1, steps_ / std::max(1, datadt_))), // number of stored datapoints minus 1
    last_data_idx(steps),
    totalT(dt_ * steps_),
    temp(temp_),
    rho(rho_),
    mu(mu_),
    d(d_),
    l(d_), // matches your Python init: self.l = d
    ir{0.0, 0.0, 0.0},
    iv{0.0, 0.0, 0.0},
    ith{0.0, 0.0},
    iw{0.0, 0.0, 0.0},
    dsize(dsize_),
    nbins(nbins_),
    mode(std::move(mode_)),
    dir(std::move(dir_)),
    name(std::move(name_)),
    res((res_ < 0) ? steps_ : res_)
  {
    if (!dir.empty() && dir.back() != '/') dir.push_back('/');

    positions.assign(steps + 1, Vec3{0.0, 0.0, 0.0});
    velocities.assign(steps + 1, Vec3{0.0, 0.0, 0.0});
    thetaz.assign(steps + 1, Vec2{0.0, 0.0});
    orientations.assign(steps + 1, Vec3{1.0, 0.0, 0.0});
    omegas.assign(steps + 1, Vec3{0.0, 0.0, 0.0});
    omegas_objframe.assign(steps + 1, Vec3{0.0, 0.0, 0.0});

    sqdis.assign(steps + 1, 0.0);
    sqangdis.assign(steps + 1, 0.0);

    times.assign(steps + 1, 0.0);
    for (int i = 0; i <= steps; ++i) times[i] = static_cast<double>(i) * (dt * datadt);

// instatiation
    upright = Vec3{1.0, 0.0, 0.0};
    ior = upright;
    ligand = upright;
    iligand = ligand;
    equilibrium = ligand;

    rmatrix = I();


    continueSim = true;
    rl = Vec3{0.0, 0.0, 0.0};
    lambda_eq = 1.0;
    kr0 = 0.0;
    gamma_bond = 0.0;
    spring = 1.0;

    // RNG
    std::random_device rd;
    rng.seed(rd());
  }

  virtual ~Simulation() = default;


  void makeFolder() {
    if (mode == "save") {
      int i = 1;
      while (true) {
        std::filesystem::path p = std::filesystem::path(dir) / (name + type);
        std::error_code ec;
        bool created = std::filesystem::create_directory(p, ec);

        if (!created || ec) {
          ++i;
          if (i == 2) {
            name += "v2";
          } else {
            if (!name.empty()) name.pop_back();
            name += std::to_string(i);
          }
        } else {
          // write parameters
          std::ofstream file((p / "parameters.txt").string());
          file << toString();
          break;
        }
      }
    }
    name += type; 
  }

  void particleData() {
    upright = Vec3{1.0, 0.0, 0.0};
    ior = upright;

    orientations.assign(steps + 1, Vec3{0.0, 0.0, 0.0});
    orientations[0] = upright;

    vol = 0.0;
    mmat = I();         // identity by default
    imat = I();
    kbTm = I();
    kbTI = I();
    beta0 = I();
    betarot0 = I();
    diffusion = 0.0;
    rotdiffusion = 0.0;
    expmsd = diffusion * 2.0 * dt;
    expmsad = rotdiffusion * 2.0 * dt;

    rmatrix = I();

    // translational kinematics precompute
    bdt = mul(beta0, dt); // beta0 * dt 

    rvar = computeRVar(bdt, kbTm);
    vvar = computeVVar(bdt, kbTm);
    rvcorr = computeCorr(vvar, rvar, kbTm, beta0, bdt);

    // rotational kinematics precompute
    brdt = mul(betarot0, dt);
    thvar = computeRVar(brdt, kbTI);
    wvar  = computeVVar(brdt, kbTI);
    thwcorr = computeCorr(wvar, thvar, kbTI, betarot0, brdt);

    spring = 1.0;
    ligand = upright;
    iligand = ligand;
    equilibrium = ligand;
  }


  void next_data() {
    // bondforce
    Vec3 bondforce{0.0, 0.0, 0.0};
    if (continueSim) {
      Vec3 rr = add(ir, iligand); // position of ligand
      double ilambda = norm(sub(rl, rr));
      if (ilambda > 0.0) {
        bondforce = mul(sub(rl, rr), spring * (1.0 - lambda_eq / ilambda));
      }
      double delta = std::abs(lambda_eq - ilambda);
      double kr = kr0 * std::exp(gamma_bond * spring * delta / (kb * temp));
      double p_r = 1.0 - std::exp(-kr * dt);

      // std::cout << p_r << "\n";

      if (uniform01(rng) < p_r) {
        continueSim = false;
      }
    }

    // translational kinematics 
    Mat3 c0 = e(neg(bdt));
    Mat3 c1 = mul(sub(I(), c0), inv(beta0));
    Mat3 c2 = mul(sub(mul(I(), dt), c1), inv(beta0));

    Vec3 randomr = sampleNormal3();
    randomr = diagSqrtMul(rvar, randomr);

    // randomv correlated (diag assumption)
    Vec3 randomv = sampleCorrelatedVelocity(randomr, vvar, rvcorr, rvar);

    // rotate noise into lab frame
    randomr = matVec(rmatrix, randomr);
    randomv = matVec(rmatrix, randomv);

    Vec3 k = matVec(inv(mmat), bondforce);

    Vec3 dr = add(
      randomr,
      add(
        matVec(matMul(matMul(rmatrix, c1), transpose(rmatrix)), iv),
        matVec(matMul(matMul(rmatrix, c2), transpose(rmatrix)), k)
      )
    );

    sqdr = dot(dr, dr);
    ir = add(ir, dr);

    iv = add(
      matVec(matMul(matMul(rmatrix, c0), transpose(rmatrix)), iv),
      add(
        randomv,
        matVec(matMul(matMul(rmatrix, c1), transpose(rmatrix)), k)
      )
    );

    c0 = e(neg(brdt));
    c1 = mul(sub(I(), c0), inv(betarot0));
    c2 = mul(sub(mul(I(), dt), c1), inv(betarot0));

    Vec3 torque = cross(iligand, bondforce);

    Vec3 randomth = sampleNormal3();
    randomth = diagSqrtMul(thvar, randomth);

    Vec3 randomw = sampleCorrelatedVelocity(randomth, wvar, thwcorr, thvar);

    randomth = matVec(rmatrix, randomth);
    randomw  = matVec(rmatrix, randomw);

    Vec3 dth = add(
      randomth,
      add(
        matVec(matMul(matMul(rmatrix, c1), transpose(rmatrix)), iw),
        matVec(
          matMul(
            matMul(rmatrix, c2),
            matMul(inv(imat), transpose(rmatrix))
          ),
          torque
        )
      )
    );

    // store angular increment magnitude^2 (Python had it commented; this is useful)
    sqdth = dot(dth, dth);

    iw = add(
      matVec(matMul(matMul(rmatrix, c0), transpose(rmatrix)), iw),
      add(
        randomw,
        matVec(
          matMul(matMul(rmatrix, c1), matMul(inv(imat), transpose(rmatrix))),
          torque
        )
      )
    );

    rmatrix = rotateMatrix(dth, rmatrix); // assumes left-multiply (lab-frame increment)
    iligand = matVec(rmatrix, ligand);
    ior = matVec(rmatrix, upright);

    ith = Vec2{ theta(ior[0], ior[1]),
                ior[2] / std::sqrt(dot(upright, upright)) };
  }

  void generateData() {
    const int rawN = steps * datadt; // matches Python loop upper bound
    const int progressEvery = std::max(1, rawN / 1000);

    for (int i = 1; i <= rawN; ++i) {
      next_data();

      if (!continueSim) {
        // stop at raw i; last stored datapoint is floor(i/datadt)
        last_data_idx = std::min(steps, i / std::max(1, datadt));
        updateProgress("Stopped at " + std::to_string(i * dt) + " s");
        break;
      }

      if (i % datadt == 0) {
        int idx = i / datadt;
        if (idx >= 0 && idx <= steps) {
          positions[idx] = ir;
          velocities[idx] = iv;
          thetaz[idx] = ith;
          omegas[idx] = iw;

          // object-frame omega 
          omegas_objframe[idx] = matVec(transpose(rmatrix), iw);

          orientations[idx] = matVec(rmatrix, upright);

          // squared displacement measures at stored times
          sqdis[idx] = sqdr;
          sqangdis[idx] = sqdth;
        }
      }

      if (mode == "save" && (i % progressEvery == 0)) {
        std::filesystem::path p = std::filesystem::path(dir) / name;
        std::ofstream file((p / "progress.txt").string());
        file << "iteration " << i;
      }
    }

    // ensure last_data_idx is full length, if simulation stops early
    last_data_idx = std::min(last_data_idx, steps);
  }


  void saveSimCSV() {
    if (mode != "save") return;

    std::filesystem::path base = std::filesystem::path(dir) / name;
    std::filesystem::create_directories(base);

    writeVec3CSV((base / "positions.csv").string(), positions, last_data_idx, {"x","y","z"});
    writeVec3CSV((base / "velocities.csv").string(), velocities, last_data_idx, {"x","y","z"});
    writeVec2CSV((base / "thetaz.csv").string(), thetaz, last_data_idx, {"theta","z"});
    writeVec3CSV((base / "omegas.csv").string(), omegas, last_data_idx, {"roll","pitch","yaw"});
    writeVec3CSV((base / "omegas_objframe.csv").string(), omegas_objframe, last_data_idx,
                 {"object_x","object_y","object_z"});
    writeVec3CSV((base / "orientations.csv").string(), orientations, last_data_idx, {"x","y","z"});
    writeScalarCSV((base / "squared_displacement.csv").string(), sqdis, sqangdis, last_data_idx);

    updateProgress("CSV data saved");
  }

  // update progress
  virtual void updateProgress(const std::string& msg) {
    std::cout << msg << std::endl;
  }


  std::string type;

  double dt;                 // s
  int datadt;                // steps between saved datapoints
  int total_raw_steps;       // raw steps passed in
  int steps;                 // number of saved intervals (arrays are steps+1)
  int last_data_idx;         // last valid saved index
  double totalT;             // total time in seconds
  double temp;               // K
  double rho;                // kg/m^3
  double mu;                 // N*s/m^2
  double d;                  // m
  double l;                  // m (set to d in your init)

  Vec3 ir;                   // position
  Vec3 iv;                   // velocity
  Vec2 ith;                  // theta,z
  Vec3 iw;                   // angular velocity

  std::vector<Vec3> positions;
  std::vector<Vec3> velocities;
  std::vector<Vec2> thetaz;
  std::vector<Vec3> orientations;
  std::vector<Vec3> omegas;
  std::vector<Vec3> omegas_objframe;
  std::vector<double> sqdis;
  std::vector<double> sqangdis;
  std::vector<double> times;

  double dsize;
  int nbins;
  std::string mode;
  std::string dir;
  std::string name;
  int res;

  // particle / dynamics fields from particleData()
  Vec3 upright{1.0, 0.0, 0.0};
  Vec3 ior{1.0, 0.0, 0.0};

  double vol{0.0};
  Mat3 mmat{I()};
  Mat3 imat{I()};
  Mat3 kbTm{I()};
  Mat3 kbTI{I()};
  Mat3 beta0{I()};
  Mat3 betarot0{I()};
  double diffusion{0.0};
  double rotdiffusion{0.0};
  double expmsd{0.0};
  double expmsad{0.0};

  Mat3 rmatrix{I()};

  Mat3 bdt{I()};
  Mat3 rvar{I()};
  Mat3 vvar{I()};
  Mat3 rvcorr{I()};

  Mat3 brdt{I()};
  Mat3 thvar{I()};
  Mat3 wvar{I()};
  Mat3 thwcorr{I()};

  double spring{1.0};
  Vec3 ligand{1.0, 0.0, 0.0};
  Vec3 iligand{1.0, 0.0, 0.0};
  Vec3 equilibrium{1.0, 0.0, 0.0};

  // bond-related (used in next_data)
  bool continueSim{true};
  Vec3 rl{0.0, 0.0, 0.0};
  double lambda_eq{1.0};
  double kr0{0.0};
  double gamma_bond{0.0};

  // last-step scratch
  double sqdr{0.0};
  double sqdth{0.0};

  // constants
  static constexpr double kb = 1.380649e-23;

private:
  // RNG
  std::mt19937_64 rng;
  std::uniform_real_distribution<double> uniform01{0.0, 1.0};
  std::normal_distribution<double> normal01{0.0, 1.0};

  std::string toString() const {
    std::ostringstream oss;
    oss << std::setprecision(17);
    oss << "type=" << type << "\n";
    oss << "dt=" << dt << "\n";
    oss << "datadt=" << datadt << "\n";
    oss << "total_raw_steps=" << total_raw_steps << "\n";
    oss << "steps(saved)=" << steps << "\n";
    oss << "temp=" << temp << "\n";
    oss << "rho=" << rho << "\n";
    oss << "mu=" << mu << "\n";
    oss << "d=" << d << "\n";
    oss << "l=" << l << "\n";
    oss << "mode=" << mode << "\n";
    oss << "dir=" << dir << "\n";
    oss << "name=" << name << "\n";
    oss << "res=" << res << "\n";
    return oss.str();
  }

  static Mat3 I() {
    return Mat3{{
      {{1.0, 0.0, 0.0}},
      {{0.0, 1.0, 0.0}},
      {{0.0, 0.0, 1.0}}
    }};
  }

  static Mat3 neg(const Mat3& A) {
    Mat3 R = A;
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) R[i][j] = -R[i][j];
    return R;
  }

  static Mat3 add(const Mat3& A, const Mat3& B) {
    Mat3 R{};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) R[i][j] = A[i][j] + B[i][j];
    return R;
  }

  static Mat3 sub(const Mat3& A, const Mat3& B) {
    Mat3 R{};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) R[i][j] = A[i][j] - B[i][j];
    return R;
  }

  static Mat3 mul(const Mat3& A, double s) {
    Mat3 R{};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) R[i][j] = A[i][j] * s;
    return R;
  }

  static Mat3 mul(double s, const Mat3& A) { return mul(A, s); }

  static Mat3 matMul(const Mat3& A, const Mat3& B) {
    Mat3 R{};
    for (int i=0;i<3;++i) {
      for (int j=0;j<3;++j) {
        double sum = 0.0;
        for (int k=0;k<3;++k) sum += A[i][k] * B[k][j];
        R[i][j] = sum;
      }
    }
    return R;
  }

  static Vec3 matVec(const Mat3& A, const Vec3& v) {
    return Vec3{
      A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2],
      A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2],
      A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2]
    };
  }

  static Mat3 transpose(const Mat3& A) {
    Mat3 R{};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j) R[i][j] = A[j][i];
    return R;
  }

  // Assumes diagonal matrices 
  static Mat3 e(const Mat3& diagA) {
    Mat3 R = Mat3{{ {{0,0,0}}, {{0,0,0}}, {{0,0,0}} }};
    for (int i=0;i<3;++i) R[i][i] = std::exp(diagA[i][i]);
    return R;
  }

  static Mat3 inv(const Mat3& diagA) {
    Mat3 R = Mat3{{ {{0,0,0}}, {{0,0,0}}, {{0,0,0}} }};
    for (int i=0;i<3;++i) {
      if (diagA[i][i] == 0.0) throw std::runtime_error("inv(): zero diagonal");
      R[i][i] = 1.0 / diagA[i][i];
    }
    return R;
  }

  static Mat3 mul(const Mat3& A, const Mat3& B) { // matrix-matrix
    return matMul(A, B);
  }

  static Mat3 mul(const Mat3& diagA, const Mat3& diagB, bool /*diagOnly*/) {
    return matMul(diagA, diagB);
  }


  static Mat3 mul(const Mat3& A, double dt) {
    return mul(A, dt);
  }

  static Vec3 add(const Vec3& a, const Vec3& b) { return {a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
  static Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
  static Vec3 mul(const Vec3& a, double s) { return {a[0]*s, a[1]*s, a[2]*s}; }

  static double dot(const Vec3& a, const Vec3& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
  static double norm(const Vec3& a) { return std::sqrt(dot(a,a)); }

  static Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3{
      a[1]*b[2] - a[2]*b[1],
      a[2]*b[0] - a[0]*b[2],
      a[0]*b[1] - a[1]*b[0]
    };
  }

  // theta(x,y)=arctan2(x,y)
  static double theta(double x, double y) {
    return std::atan2(x, y);
  }

  Vec3 sampleNormal3() {
    return Vec3{ normal01(rng), normal01(rng), normal01(rng) };
  }

  // Multiply vector by sqrt(diagonal matrix) (assumes diagonal)
  static Vec3 diagSqrtMul(const Mat3& diagA, const Vec3& v) {
    return Vec3{
      v[0] * std::sqrt(std::max(0.0, diagA[0][0])),
      v[1] * std::sqrt(std::max(0.0, diagA[1][1])),
      v[2] * std::sqrt(std::max(0.0, diagA[2][2]))
    };
  }

  // randomv formula under diagonal assumption
  Vec3 sampleCorrelatedVelocity(const Vec3& randomr_scaled,
                                const Mat3& vvarDiag,
                                const Mat3& corrDiag,
                                const Mat3& rvarDiag)
  {
    Vec3 z = sampleNormal3();

    Vec3 out{};
    for (int i=0;i<3;++i) {
      const double vv = vvarDiag[i][i];
      const double rr = rvarDiag[i][i];
      const double c  = corrDiag[i][i];

      const double term1 = z[i] * std::sqrt(std::max(0.0, vv * (1.0 - c*c)));
      const double term2 = randomr_scaled[i] * c * std::sqrt(std::max(0.0, vv / rr));
      out[i] = term1 + term2;
    }
    return out;
  }

  static Mat3 computeRVar(const Mat3& bdtDiag, const Mat3& kbTDiag) {
    Mat3 R = Mat3{{ {{0,0,0}}, {{0,0,0}}, {{0,0,0}} }};
    for (int i=0;i<3;++i) {
      const double b = bdtDiag[i][i];
      const double kbT = kbTDiag[i][i];
      if (b == 0.0) throw std::runtime_error("computeRVar(): zero bdt diagonal");
      const double invb = 1.0 / b;
      const double expr = 2.0*invb + (invb*invb) * (-3.0 + 4.0*std::exp(-b) - std::exp(-2.0*b));
      R[i][i] = (/*dt^2 handled outside by caller via bdt=beta*dt*/ 1.0) * kbT * expr;
     // applied dt^2 in caller by passing bdt computed as beta*dt and multiplying by dt^2 there.
    }
    return R;
  }

  static Mat3 computeVVar(const Mat3& bdtDiag, const Mat3& kbTDiag) {
    Mat3 R = Mat3{{ {{0,0,0}}, {{0,0,0}}, {{0,0,0}} }};
    for (int i=0;i<3;++i) {
      const double b = bdtDiag[i][i];
      const double kbT = kbTDiag[i][i];
      R[i][i] = kbT * (1.0 - std::exp(-2.0*b));
    }
    return R;
  }

  static Mat3 computeCorr(const Mat3& vvarDiag,
                          const Mat3& rvarDiag,
                          const Mat3& kbTDiag,
                          const Mat3& betaDiag,
                          const Mat3& bdtDiag)
  {
    Mat3 R = Mat3{{ {{0,0,0}}, {{0,0,0}}, {{0,0,0}} }};
    for (int i=0;i<3;++i) {
      const double vv = vvarDiag[i][i];
      const double rr = rvarDiag[i][i];
      const double kbT = kbTDiag[i][i];
      const double beta = betaDiag[i][i];
      const double b = bdtDiag[i][i];

      if (vv <= 0.0 || rr <= 0.0 || beta == 0.0) {
        R[i][i] = 0.0;
        continue;
      }
      const double inv_sqrt_vvrr = 1.0 / std::sqrt(vv * rr);
      const double factor = (1.0 - std::exp(-b));
      R[i][i] = inv_sqrt_vvrr * kbT * (1.0 / beta) * (factor * factor);
    }
    return R;
  }

  // Rodrigues incremental rotation from axis-angle vector dth; returns Rinc * R
  static Mat3 rotateMatrix(const Vec3& dth, const Mat3& R) {
    const double th = norm(dth);
    if (th == 0.0) return R;

    const Vec3 u{ dth[0]/th, dth[1]/th, dth[2]/th };
    const double ux=u[0], uy=u[1], uz=u[2];
    const double c = std::cos(th);
    const double s = std::sin(th);
    const double omc = 1.0 - c;

    Mat3 Rinc{{
      {{ c + ux*ux*omc,      ux*uy*omc - uz*s, ux*uz*omc + uy*s }},
      {{ uy*ux*omc + uz*s,   c + uy*uy*omc,    uy*uz*omc - ux*s }},
      {{ uz*ux*omc - uy*s,   uz*uy*omc + ux*s, c + uz*uz*omc    }}
    }};

    return matMul(Rinc, R);
  }

  // CSV writers
  static void writeVec3CSV(const std::string& path,
                           const std::vector<Vec3>& data,
                           int lastIdx,
                           const std::array<std::string,3>& header)
  {
    std::ofstream f(path);
    f << header[0] << "," << header[1] << "," << header[2] << "\n";
    for (int i=0; i<=lastIdx; ++i) {
      f << data[i][0] << "," << data[i][1] << "," << data[i][2] << "\n";
    }
  }

  static void writeVec2CSV(const std::string& path,
                           const std::vector<Vec2>& data,
                           int lastIdx,
                           const std::array<std::string,2>& header)
  {
    std::ofstream f(path);
    f << header[0] << "," << header[1] << "\n";
    for (int i=0; i<=lastIdx; ++i) {
      f << data[i][0] << "," << data[i][1] << "\n";
    }
  }

  static void writeScalarCSV(const std::string& path,
                             const std::vector<double>& trans,
                             const std::vector<double>& ang,
                             int lastIdx)
  {
    std::ofstream f(path);
    f << "Translational,Angular\n";
    for (int i=0; i<=lastIdx; ++i) {
      f << trans[i] << "," << ang[i] << "\n";
    }
  }


  //   rvar = computeRVar(...) ; then rvar = rvar * (dt*dt)
  // Same for thvar.
  //implemented multiply-in-place for diagonal matrices:
  static Mat3 mulDiag(const Mat3& A, double s) {
    Mat3 R = Mat3{{ {{0,0,0}}, {{0,0,0}}, {{0,0,0}} }};
    for (int i=0;i<3;++i) R[i][i] = A[i][i] * s;
    return R;
  }

  // correct dt^2 factor
public:
  void particleData_with_dt2_fix() {
    particleData();

    rvar  = mulDiag(rvar,  dt*dt);
    thvar = mulDiag(thvar, dt*dt);
  }
};
