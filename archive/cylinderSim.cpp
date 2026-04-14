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
