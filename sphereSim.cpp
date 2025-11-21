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
