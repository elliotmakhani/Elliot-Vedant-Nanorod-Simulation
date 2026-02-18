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

                  
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

class Simulation {
public:
    std::string type;
    double dt; // s
    int datadt; // int representing the number of steps between data points
    int steps;
    int laststep;
    double totalt; // s
    double temp; // K
    double rho; // density, kg/m^3
    double mu; // viscosity, N*s/m^2
    double d; // diameter, m
    double l;
    Eigen::Vector3d ir; // position vector
    Eigen::Vector3d iv; // velocity vector
    Eigen::Vector2d ith; // theta, z
    Eigen::Vector3d iw; // angular velocity vector
    Eigen::MatrixXd positions; // all positions
    Eigen::MatrixXd velocities; // all velocities
    Eigen::MatrixXd thetaz; // all thetas
    Eigen::MatrixXd orientations; // all orientations
    Eigen::MatrixXd omegas; // all omegas
    Eigen::MatrixXd omegasobjframe; // all omegas in object frame
    Eigen::VectorXd sqdis; // squared displacements
    Eigen::VectorXd sqangdis; // squared angular displacements
    std::vector<std::pair<double, double>> squared_displacement; // combined squared displacements
    Eigen::VectorXd times;
    double dsize;
    int nbins;
    std::string mode; // either 'show' or 'save'
    std::string dir;
    int res;
    std::string name;
    

    Simulation(double dt, int datadt, int steps, double temp, double rho, double mu, double d, double dsize, int nbins, std::string mode = "show", std::string dir = "", int res = -1, std::string name = ""): dt(dt), datadt(datadt), steps(steps), temp(temp), rho(rho), mu(mu), d(d), dsize(dsize), nbins(nbins), mode(mode), dir(dir), res(res), name(name) {
        totalt = dt * steps; // s
        laststep = steps;
        ir = Eigen::Vector3d::Zero(); // position vector
        iv = Eigen::Vector3d::Zero(); // velocity vector
        ith = Eigen::Vector2d::Zero(); // theta, z
        iw = Eigen::Vector3d::Zero(); // angular velocity vector
        positions = Eigen::MatrixXd::Zero(steps + 1, 3);
        positions.row(0) = ir;
        velocities = Eigen::MatrixXd::Zero(steps + 1, 3);
        velocities.row(0) = iv;
        thetaz = Eigen::MatrixXd::Zero(steps + 1, 2);
        thetaz.row(0) = ith;
        orientations = Eigen::MatrixXd::Zero(steps + 1, 3);
        orientations.row(0) = Eigen::Vector3d(1.0, 0.0, 0.0);
        omegas = Eigen::MatrixXd::Zero(steps + 1, 3);
        omegasobjframe = Eigen::MatrixXd::Zero(steps + 1, 3);
        omegas.row(0) = iw;
        sqdis = Eigen::VectorXd::Zero(steps + 1); // squared displacements
        sqangdis = Eigen::VectorXd::Zero(steps + 1); // squared angular displacements
        squared_displacement.push_back(std::make_pair(sqdis(0), sqangdis(0)));
        times = Eigen::VectorXd::LinSpaced(steps + 1, 0, totalt);
        if (!dir.empty() && dir.back() != '/') {
            dir += '/';
        }
        if (res == -1) {
            this->res = steps;
        }
    }
    void makeFolder(){
        if (mode == "save") {
            int i = 1;
            std::string full_path = dir;
            full_path /= name + type;
            while (true) {                
                try {
                    if(fs::create_directory(full_path)) {
                        break;
                    } else {
                        throw std::runtime_error("Directory already exists");
                    }
                } catch (std::runtime_error& e) {
                    i += 1;
                    full_path = dir;
                    full_path /= name + type + "v" + std::to_string(i);
                }
            }
            // Write parameters to file logic here
        }
    }
    void particleData();
    void next_data();
    void generateData();
    void saveSim();
};
