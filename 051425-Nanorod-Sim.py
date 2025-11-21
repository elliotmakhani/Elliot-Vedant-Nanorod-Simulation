import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from math import log
from matplotlib.widgets import RangeSlider
from matplotlib import cm
from pathlib import Path
import time
import os
# constants/static functions

start = time.time()
graphboxsize = 10**-5 # m
kb = 1.380649 * (10 ** -23) # boltzmann constant, J/K
def e(matrix):
    return np.diag(np.exp(np.diag(matrix)))

def inv(matrix):
    return np.diag(1 / np.diag(matrix))

def theta(x, y):
    return np.arctan2(x, y)


def rotateMatrix(dth,rmatrix):
    theta = np.linalg.norm(dth)
    dth /= theta  # Ensure the axis is a unit vector
    u_x, u_y, u_z = dth
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    one_minus_cos_t = 1 - cos_t

    return np.array([
        [
            cos_t + u_x**2 * one_minus_cos_t,
            u_x * u_y * one_minus_cos_t - u_z * sin_t,
            u_x * u_z * one_minus_cos_t + u_y * sin_t
        ],
        [
            u_y * u_x * one_minus_cos_t + u_z * sin_t,
            cos_t + u_y**2 * one_minus_cos_t,
            u_y * u_z * one_minus_cos_t - u_x * sin_t
        ],
        [
            u_z * u_x * one_minus_cos_t - u_y * sin_t,
            u_z * u_y * one_minus_cos_t + u_x * sin_t,
            cos_t + u_z**2 * one_minus_cos_t
        ]
    ]) @ rmatrix

def closeOpen():
    if not bool(plt.get_fignums()):
        plt.close()
class Simulation:
    def __init__(self,dt,datadt,steps,temp,rho,mu,d,dsize,nbins, mode= 'show',dir="",res=None,name =''):
        self.type = ''
        self.dt = dt # s
        self.datadt = datadt # int representing the number of steps between data points
        self.steps = int(steps/datadt)
        self.laststep = self.steps
        self.totalt = dt*steps # s
        self.temp = temp # K
        self.rho = rho # density, kg/m^3
        self.mu = mu # viscosity, N*s/m^2
        self.d = d # diameter, m
        self.l = d 
        self.ir = np.array([0.,0.,0.]) # position vector
        self.iv = np.array([0.,0.,0.]) # velocity vector
        self.ith = np.array([0.,0.]) # theta, z
        self.iw = np.array([0.,0.,0.]) # angular velocity vector
        
        self.positions = np.empty((self.steps+1,3))
        self.positions[0] = self.ir.copy()

        self.velocities = np.empty((self.steps+1,3))
        self.velocities[0] = self.iv.copy()

        self.thetaz = np.empty((self.steps+1,2))
        self.thetaz[0] = self.ith.copy()  # initial position assumed the axis of the cylinder is aligned with the z axis

        self.orientations = np.empty((self.steps+1,3))
        self.orientations[0] = np.array([1.,0.,0.])
        
        self.omegas = np.empty((self.steps+1,3))
        self.omegasobjframe = np.empty((self.steps+1,3))
        self.omegas[0] = self.iw.copy()
        self.omegasobjframe[0] = np.array([0.,0.,0.])

        
        self.sqdis = np.empty(self.steps+1) # squared displacements
        self.sqangdis = np.empty(self.steps+1) # squared angular displacements
        self.sqdis[0] = 0
        self.sqangdis[0] = 0
    
        self.squared_displacement = [self.sqdis.copy(),self.sqangdis.copy()]
        self.times = np.array(range(0,self.steps+1)) * (dt * datadt)
        self.dsize = dsize
        self.nbins = nbins
        self.mode = mode # either 'show' or 'save'
        self.dir = dir
        self.name = name
        if (not(self.dir == "")):
            self.dir += '/'
        if res == None:
            self.res = steps
        else: 
            self.res = res


    def makeFolder(self):
        if (self.mode == "save"):
            i = 1
            while(True):
                try:
                    Path(self.dir + self.name + self.type).mkdir()
                except FileExistsError:
                    i += 1
                    if (i==2):
                        self.name += 'v2'
                    else:
                        self.name = self.name[0:-1:] + str(i)
                else:
                    break
            

            with open(self.dir + self.name + self.type  + "/parameters.txt",'w') as file:
                file.write(str(self)) 
        self.name += self.type      
    def particleData(self):
        self.upright = np.array([1,0,0])
        self.ior = self.upright
        self.orientations = np.empty((self.steps+1,3))
        self.orientations[0] = self.upright.copy()
        self.vol = 0
        self.mmat = np.eye(3)
        self.imat = np.eye(3)
        self.kbTm = np.eye(3)
        self.kbTI = np.eye(3)
        self.beta0 = np.eye(3)
        self.betarot0 = np.eye(3)
        self.diffusion = 0
        self.rotdiffusion = 0
        self.expmsd = self.diffusion*2*self.dt
        self.expmsad = self.rotdiffusion*2*self.dt
        self.rmatrix = np.eye(3)
        #translational kinematics data
        self.bdt = self.beta0*self.dt
        #print(self.bdt)
        self.rvar = self.dt**2 * self.kbTm @ (2*inv(self.bdt) + inv(self.bdt)**2 @ (-3*np.eye(3)+4*e(-self.bdt)-e(-2*self.bdt)))
        #print(self.rvar)
        self.vvar = self.kbTm @ (np.eye(3)-e(-2*self.bdt))
        #print(self.vvar)
        self.rvcorr =  (inv(self.vvar @ self.rvar)**.5) @ self.kbTm @ inv(self.beta0) @ ((np.eye(3)-e(-1*self.bdt))**2)
        #rotational kinematics data
        self.brdt = self.betarot0 * self.dt 
        self.thvar = self.dt**2 * self.kbTI @ (2*inv(self.brdt) + inv(self.brdt)**2 @ (-3*np.eye(3)+4*e(-self.brdt)-e(-2*self.brdt)))
        self.wvar = self.kbTI @ (np.eye(3)-e(-2*self.brdt)) 
        self.thwcorr = (inv(self.wvar @ self.thvar)**.5) @ self.kbTI @ inv(self.betarot0) @ ((np.eye(3)-e(-1*self.brdt))**2)
        self.spring = 1
        self.ligand = self.upright 
        self.iligand = self.ligand.copy() # instantaneous ligand position relative to the CM
        self.equilibrium = self.ligand # equilibrium position (aka the position of the ligand when the bond is formed)
        
    def next_data(self):
        #bondforce = (self.equilibrium - self.ir - self.iligand) * self.spring # position vector + arm vector - vector of equilibrium position (upright)
        
        if self.continueSim:
            rr = self.ir + self.iligand # position of ligand (m)
            ilambda = np.linalg.norm(self.rl - rr) 
            bondforce = self.spring * (1 - self.lambda_eq/ilambda) * (self.rl - rr) # force on the ligand due to the spring
            delta = np.abs(self.lambda_eq - ilambda) # distance from the equilibrium position of the ligand
            kr = self.kr0 * np.exp(self.gamma_bond * self.spring * delta/(kb*self.temp)) # spring constant of the bond
            p_r = 1-np.exp(-kr*self.dt) # probability of breaking the bond
            print(p_r)
            
            if (np.random.rand() < p_r): # break the bond with probability p_r
                self.continueSim = False
        else:
            bondforce = np.array([0,0,0])
        #print(f"bond force: {bondforce}")
        #print(f"tip position:{self.ir+self.ior}")
        #n = 2*self.ior/self.l # unit normal vector of the circular tip of the cylinder
        #zpos = np.array([0,0,1]) # z axis, pointing towards the "wall"
        #corner0 = zpos - zpos @ n * n # the corner of the cylinder closest to the wall
        #corner0 = self.d*corner0/(2*np.linalg.norm(corner0)) + self.ir + self.ior # position of the corner

        # translational kinematics
        c0 = e(-self.bdt)
        c1 = (np.eye(3)-c0) @ inv(self.beta0)
        c2 = (self.dt*np.eye(3)-c1) @ inv(self.beta0)
        randomr = np.random.normal(loc=0,scale=1,size=3) @ self.rvar**.5
        #print(f"randomr1:{randomr}")
        randomv = np.random.normal(loc=0, scale=1, size=3) @ (self.vvar@(np.eye(3)-self.rvcorr**2))**.5 + randomr @ self.rvcorr @ (self.vvar @ inv(self.rvar)) ** .5
        #print(f"randomv:{randomv}")
        randomr = np.eye(3) @ randomr
        #print(f"randomr2:{randomr}")
        #print(f"rmatrix:{self.rmatrix}")
        randomr = self.rmatrix @ randomr
        #print(f"randomr3:{randomr}")
        randomv = self.rmatrix @ randomv
        
        k =  inv(self.mmat) @ bondforce
        #print(f"k:{k}")
        #print(f"randomv2:{randomv}")
        dr = randomr + self.rmatrix @ c1 @ self.rmatrix.T @ self.iv + self.rmatrix @ c2 @ self.rmatrix.T @ k
        #print(f"dr:{dr}")
        self.sqdr = dr @ dr
        self.ir += dr
        self.iv = (self.rmatrix @ c0 @ self.rmatrix.T) @ self.iv + randomv + self.rmatrix @ c1 @ self.rmatrix.T @ k 

        # rotational kinematics 
        c0 = e(-self.brdt)
        c1 = (np.eye(3)-c0) @ inv(self.betarot0)
        c2 = (self.dt*np.eye(3)-c1) @ inv(self.betarot0)

        #print(f"c2:{c2}")

        torque = np.cross(self.iligand,bondforce) 
        #print(f"torque: {torque}")
        randomth = np.random.normal(loc=0,scale=1,size=3) @ self.thvar**.5
        #print(f"randomth2:{randomth}")
        randomw = np.random.normal(loc=0, scale=1, size=3) @ (self.wvar@(np.eye(3)-self.thwcorr**2))**.5 + randomth @ self.thwcorr @ (self.wvar @ inv(self.thvar)) ** .5
        randomth = self.rmatrix @ randomth
        #print(f"randomth2:{randomth}")
        randomw = self.rmatrix @ randomw
        dth = randomth + self.rmatrix @ c1 @ self.rmatrix.T @ self.iw +  self.rmatrix @ c2 @ inv(self.imat) @ self.rmatrix.T @ torque
        #print(f"dth:{dth}")
        #print(f"dth-randomth:{dth-randomth}")
        #print(f"dth:{dth}")
        #self.sqdth = dth @ dth
        self.iw = (self.rmatrix @ c0 @ self.rmatrix.T) @ self.iw + randomw + self.rmatrix @ c1 @ inv(self.imat) @ self.rmatrix.T @ torque
        self.rmatrix = rotateMatrix(dth,self.rmatrix)
        self.iligand = self.rmatrix @ self.ligand
        self.ior = self.rmatrix @ self.upright
        self.ith = [theta(self.ior[0],self.ior[1]),self.ior[2]/np.sqrt(self.upright@self.upright)]

        #n = 2*self.ior/self.l # checks corner position after movement
        #cornerf = zpos - zpos @ n * n #
        #cornerf = self.d*corner0/(2*np.linalg.norm(corner0)) + self.ir + self.ior
        #collision = (cornerf - self.wall(1)) @ self.wall(0) < 0 # checks if the corner of the cylinder has gone past the wall


    def generateData(self):
        for i in range(1,self.steps*self.datadt+1):
            
            #print(f"step {i}")
            #print("position")
            self.next_data()
            if (not(self.continueSim)):
                self.laststep = i
                self.updateProgress(f"Stopped at {i*self.dt} s")
                break
            #self.ir = (self.ir+graphboxsize/2) % graphboxsize - graphboxsize/2
            # iomobj = self.rmatrix.T@self.iw
            if (i%self.datadt == 0):
                idx = int(i/self.datadt)
                self.positions[idx] = self.ir.copy()
                self.velocities[idx] = self.iv.copy()
                self.thetaz[idx] = self.ith.copy()
                self.omegas[idx] = self.iw.copy()
                #self.omegasobjframe.append(iomobj)
                #self.sqdis.append(self.sqdr.copy())
                #self.sqangdis.append(self.sqdth.copy())
                self.orientations[idx] = self.rmatrix@self.upright
                #self.absoluteorientations = [(rotateMatrix(self.ith)@np.array([0,0,1])) @ (rotateMatrix(self.ith)@np.array([0,0,1]))]
                #print(self.absoluteorientations)
                # print('test-1')
            if ((i)%(self.steps*self.datadt*10**-3) == 0):
                with open(self.dir + self.name + "/progress.txt",'w') as file:
                    file.write(f"iteration {i}")
        # magnitudes = {(ornt @ ornt)/(self.upright@self.upright) for ornt in self.orientations}
        # print(f"\n{min(magnitudes)},{max(magnitudes)}")
        self.sqdis= np.array(self.sqdis[0:self.laststep+1])
        self.sqangdis = np.array(self.sqangdis[0:self.laststep+1])
        #print(self.positions)
        self.positions = pd.DataFrame(self.positions[0:self.laststep+1] , columns= ['x','y','z'])
        # print(positions)
        self.velocities = pd.DataFrame(self.velocities[0:self.laststep+1], columns= ['x','y','z'])
        self.thetaz = pd.DataFrame(self.thetaz[0:self.laststep+1], columns=['theta','z'])
        # print(thetas)
        self.omegas = pd.DataFrame(self.omegas[0:self.laststep+1], columns=['roll','pitch','yaw'])
        self.omegasobjframe = pd.DataFrame(self.omegasobjframe[0:self.laststep+1], columns = ['object x-axis','object y-axis','object z-axis'])
        self.orientations = pd.DataFrame(self.orientations[0:self.laststep+1], columns = ['x','y','z'])
        self.squared_displacement =  pd.DataFrame({'Translational': self.sqdis[0:self.laststep+1],'Angular': self.sqangdis[0:self.laststep+1]})
        self.orientations *= 10**9
        self.positions*= 10**9
        self.times *= 10**9
        self.totalt *= 10**9
        self.dt *= 10**9 

    def __str__(self):
        return f'Simulation {self.name} parameters:\ntime step: {self.dt:.2e} ns\n# of steps: {self.steps*self.datadt:.2e}\nsimulation duration: {self.dt*self.steps*self.datadt:.2e} (ns)\ntemperature: {self.temp}K\nfluid density: {self.rho} kg/m^3\nviscosity: {self.mu:.4e}\ndot size: {self.dsize}\n# of bins: {self.nbins}\ndata resolution: {self.res:.2e}\ndiameter: {self.d*10**9} nm'

    def graphxv(self):
        closeOpen()
        self.fig, self.ax = plt.subplots(nrows=2,ncols=3,figsize=(18,7))
        plt.subplots_adjust(wspace=1, hspace=1)
        #position histogram
        counts, bins,patches= self.ax[0][0].hist(self.positions['x'],bins=self.nbins )
        
        # bin_width = np.diff(bins)[0]
        self.ax[0][0].axhline(np.mean(counts), color='red', linestyle='--')
        self.ax[0][0].set_title('Distribution of particle x position')
        self.ax[0][0].set_xlabel("Position (nm)")
        self.ax[0][0].set_ylabel("Frequency")
        #gauss = self.steps*bin_width/(np.std(self.positions['x'])*(np.pi*2)**.5) * np.exp(((self.positions['x']-np.mean(self.positions['x']))/np.std(self.positions['x']))**2*(-1/2))
        #self.ax[0][0].scatter(self.positions['x'],gauss,s=self.dsize)

        counts, bins, patches= self.ax[0][1].hist(self.positions['y'],bins=self.nbins )
        # bin_width = np.diff(bins)[0]
        self.ax[0][1].axhline(np.mean(counts), color='red', linestyle='--')
        self.ax[0][1].set_title('Distribution of particle y position')
        self.ax[0][1].set_xlabel("Position (nm)")
        self.ax[0][1].set_ylabel("Frequency")
        #gauss = self.steps*bin_width/(np.std(self.positions['y'])*(np.pi*2)**.5) * np.exp(((self.positions['y']-np.mean(self.positions['y']))/np.std(self.positions['y']))**2*(-1/2))
        #self.ax[0][1].scatter(self.positions['y'],gauss,s=self.dsize)
        
        counts, bins, patches= self.ax[0][2].hist(self.positions['z'],bins=self.nbins )
        # bin_width = np.diff(bins)[0]
        self.ax[0][2].axhline(np.mean(counts), color='red', linestyle='--')
        self.ax[0][2].set_title('Distribution of particle z position')
        self.ax[0][2].set_xlabel("Position (nm)")
        self.ax[0][2].set_ylabel("Frequency")
        #gauss = self.steps*bin_width/(np.std(self.positions['z'])*(np.pi*2)**.5) * np.exp(((self.positions['z']-np.mean(self.positions['z']))/np.std(self.positions['z']))**2*(-1/2))
        #self.ax[0][2].scatter(self.positions['z'],gauss,s=self.dsize)

        #velocity histogram
        vcounts, vbins,vpatches= self.ax[1][0].hist(self.velocities['x'],bins=self.nbins)
        vbin_width = np.diff(vbins)[0]
        self.ax[1][0].set_title('Distribution of particle x-axis Velocity')
        self.ax[1][0].set_xlabel("Velocity (nm/ns)")
        self.ax[1][0].set_ylabel("Frequency")
        gaussv = self.steps*vbin_width/(np.std(self.velocities['x'])*(np.pi*2)**.5) * np.exp(((self.velocities['x']-np.mean(self.velocities['x']))/np.std(self.velocities['x']))**2*(-1/2))
        self.ax[1][0].scatter(self.velocities['x'],gaussv,  color='red', s=self.dsize)

        vcounts, vbins,vpatches= self.ax[1][1].hist(self.velocities['y'],bins=self.nbins)
        vbin_width = np.diff(vbins)[0]
        self.ax[1][1].set_title('Distribution of particle y-axis Velocity')
        self.ax[1][1].set_xlabel("Velocity (nm/ns)")
        self.ax[1][1].set_ylabel("Frequency")
        gaussv = self.steps*vbin_width/(np.std(self.velocities['y'])*(np.pi*2)**.5) * np.exp(((self.velocities['y']-np.mean(self.velocities['y']))/np.std(self.velocities['y']))**2*(-1/2))
        self.ax[1][1].scatter(self.velocities['y'],gaussv,  color='red', s=self.dsize)
    
        vcounts, vbins,vpatches= self.ax[1][2].hist(self.velocities['z'],bins=self.nbins)
        vbin_width = np.diff(vbins)[0]
        self.ax[1][2].set_title('Distribution of particle z-axis Velocity')
        self.ax[1][2].set_xlabel("Velocity (nm/ns)")
        self.ax[1][2].set_ylabel("Frequency")
        gaussv = self.steps*vbin_width/(np.std(self.velocities['z'])*(np.pi*2)**.5) * np.exp(((self.velocities['z']-np.mean(self.velocities['z']))/np.std(self.velocities['z']))**2*(-1/2))
        self.ax[1][2].scatter(self.velocities['z'],gaussv,  color='red', s=self.dsize)
        plt.tight_layout(pad=2.0)
        # Display the plot

        plt.savefig( self.dir + self.name + '/' +  "xvgraph.png")    
    def graphthw(self):
        closeOpen()
        self.fig, self.ax = plt.subplots(nrows=3,ncols=3,figsize=(18,10))
        plt.subplots_adjust(wspace=1, hspace=1)
        #angle histogram
        counts, bins,patches= self.ax[0][0].hist(self.thetaz['theta'],bins=self.nbins )
        # bin_width = np.diff(bins)[0]
        self.ax[0][0].axhline(np.mean(counts), color='red', linestyle='--')

        self.ax[0][0].set_title('Distribution of particle yaw angle')
        self.ax[0][0].set_xlabel("Angle (rad)")
        self.ax[0][0].set_ylabel("Frequency")
        #gauss = self.steps*bin_width/(np.std(self.thetaz['roll'])*(np.pi*2)**.5) * np.exp(((self.thetaz['roll']-np.mean(self.thetaz['roll']))/np.std(self.thetaz['roll']))**2*(-1/2))
        #self.ax[0][0].scatter(self.thetaz['roll'],gauss,s=self.dsize)

        counts, bins, patches= self.ax[0][1].hist(self.thetaz['z'],bins=self.nbins )
        self.ax[0][1].axhline(np.mean(counts), color='red', linestyle='--')
        # bin_width = np.diff(bins)[0]
        self.ax[0][1].set_title('Distribution of the z-component of the unit orientation vector')
        self.ax[0][1].set_xlabel("\u1E91")
        self.ax[0][1].set_ylabel("Frequency")
        #gauss = self.steps*bin_width/(np.std(self.thetaz['pitch'])*(np.pi*2)**.5) * np.exp(((self.thetaz['pitch']-np.mean(self.thetaz['pitch']))/np.std(self.thetaz['pitch']))**2*(-1/2))
        #self.ax[0][1].scatter(self.thetaz['pitch'],gauss,s=self.dsize)
        

        #gauss = self.steps*bin_width/(np.std(self.thetaz['yaw'])*(np.pi*2)**.5) * np.exp(((self.thetaz['yaw']-np.mean(self.thetaz['yaw']))/np.std(self.thetaz['yaw']))**2*(-1/2))
        #self.ax[0][2].scatter(self.thetaz['yaw'],gauss,s=self.dsize)

        #velocity histogram
        vcounts, vbins,vpatches= self.ax[1][0].hist(self.omegas['roll'],bins=self.nbins)
        vbin_width = np.diff(vbins)[0]
        self.ax[1][0].set_title('Distribution of particle roll angular velocity (fixed frame of ref.)')
        self.ax[1][0].set_xlabel("Angular Velocity (rad/s)")
        self.ax[1][0].set_ylabel("Frequency")
        gaussv = self.steps*vbin_width/(np.std(self.omegas['roll'])*(np.pi*2)**.5) * np.exp(((self.omegas['roll']-np.mean(self.omegas['roll']))/np.std(self.omegas['roll']))**2*(-1/2))
        self.ax[1][0].scatter(self.omegas['roll'],gaussv,  color='red', s=self.dsize)

        vcounts, vbins,vpatches= self.ax[1][1].hist(self.omegas['pitch'],bins=self.nbins)
        vbin_width = np.diff(vbins)[0]
        self.ax[1][1].set_title('Distribution of particle pitch angular velocity (fixed frame of ref.)')
        self.ax[1][1].set_xlabel("Angular Velocity (rad/s)")
        self.ax[1][1].set_ylabel("Frequency")
        gaussv = self.steps*vbin_width/(np.std(self.omegas['pitch'])*(np.pi*2)**.5) * np.exp(((self.omegas['pitch']-np.mean(self.omegas['pitch']))/np.std(self.omegas['pitch']))**2*(-1/2))
        self.ax[1][1].scatter(self.omegas['pitch'],gaussv,  color='red', s=self.dsize)
    
        vcounts, vbins,vpatches= self.ax[1][2].hist(self.omegas['yaw'],bins=self.nbins)
        vbin_width = np.diff(vbins)[0]
        self.ax[1][2].set_title('Distribution of particle yaw angular velocity (fixed frame of ref.)')
        self.ax[1][2].set_xlabel("Angular Velocity (rad/s)")
        self.ax[1][2].set_ylabel("Frequency")
        gaussv = self.steps*vbin_width/(np.std(self.omegas['yaw'])*(np.pi*2)**.5) * np.exp(((self.omegas['yaw']-np.mean(self.omegas['yaw']))/np.std(self.omegas['yaw']))**2*(-1/2))
        self.ax[1][2].scatter(self.omegas['yaw'],gaussv,  color='red', s=self.dsize)
        plt.tight_layout(pad=2.0)
        '''
        vcounts, vbins,vpatches= self.ax[2][0].hist(self.omegasobjframe['object x-axis'],bins=self.nbins)
        vbin_width = np.diff(vbins)[0]
        self.ax[2][0].set_title('Distribution of particle roll angular velocity (object frame of ref.)')
        self.ax[2][0].set_xlabel("Angular Velocity (rad/s)")
        self.ax[2][0].set_ylabel("Frequency")
        gaussv = self.steps*vbin_width/(np.std(self.omegasobjframe['object x-axis'])*(np.pi*2)**.5) * np.exp(((self.omegasobjframe['object x-axis']-np.mean(self.omegasobjframe['object x-axis']))/np.std(self.omegasobjframe['object x-axis']))**2*(-1/2))
        self.ax[2][0].scatter(self.omegasobjframe['object x-axis'],gaussv,  color='red', s=self.dsize)

        vcounts, vbins,vpatches= self.ax[2][1].hist(self.omegasobjframe['object y-axis'],bins=self.nbins)
        vbin_width = np.diff(vbins)[0]
        self.ax[2][1].set_title('Distribution of particle pitch angular velocity (object frame of ref.)')
        self.ax[2][1].set_xlabel("Angular Velocity (rad/s)")
        self.ax[2][1].set_ylabel("Frequency")
        gaussv = self.steps*vbin_width/(np.std(self.omegasobjframe['object y-axis'])*(np.pi*2)**.5) * np.exp(((self.omegasobjframe['object y-axis']-np.mean(self.omegasobjframe['object y-axis']))/np.std(self.omegasobjframe['object y-axis']))**2*(-1/2))
        self.ax[2][1].scatter(self.omegasobjframe['object y-axis'],gaussv,  color='red', s=self.dsize)

        vcounts, vbins,vpatches= self.ax[2][2].hist(self.omegasobjframe['object z-axis'],bins=self.nbins)
        vbin_width = np.diff(vbins)[0]
        self.ax[2][2].set_title('Distribution of particle yaw angular velocity (object frame of ref.)')
        self.ax[2][2].set_xlabel("Angular Velocity (rad/s)")
        self.ax[2][2].set_ylabel("Frequency")
        gaussv = self.steps*vbin_width/(np.std(self.omegasobjframe['object z-axis'])*(np.pi*2)**.5) * np.exp(((self.omegasobjframe['object z-axis']-np.mean(self.omegasobjframe['object z-axis']))/np.std(self.omegasobjframe['object z-axis']))**2*(-1/2))
        self.ax[2][2].scatter(self.omegasobjframe['object z-axis'],gaussv,  color='red', s=self.dsize)
        '''
        # Display the plot
        plt.savefig( self.dir + self.name + '/' +  "thwgraph.png")

    def graph_rt_msd(self):
        closeOpen()
        mean = np.sum(self.squared_displacement['Translational'])/self.steps
        self.fig, self.ax = plt.subplots(nrows=2,figsize=(10,12))
        counts, bins, patches = self.ax[0].hist(self.squared_displacement['Translational'],bins=self.nbins)
        threshold = self.steps * 5* 10**-4
        counts = np.where(counts<threshold,0,counts)
        lastbin = 0
        for i in range(0,len(counts)):
            if counts[i] == 0:
                lastbin = i
                break
        self.ax[0].cla()
        self.ax[0].bar(bins[:-1],counts,width = np.diff(bins)) 
        self.ax[0].set_xlim(0, bins[lastbin-1])
        self.ax[0].set_title('Distribution of squared displacment')
        self.ax[0].set_xlabel("Displacement squared in a single timestep (m^2)")
        self.ax[0].set_ylabel("Frequency")
        self.ax[0].plot([mean]*2,[*self.ax[0].get_ylim()],label = f"mean squared displacement ={mean:4e}")
        self.ax[0].plot([self.expmsd]*2,[*self.ax[0].get_ylim()],label = f'2*D*dt = {self.expmsd:4e}') # \n z-score = {np.sqrt(self.steps/2)*(mean/self.expmsd-1):3f}
        self.ax[0].legend(loc='upper right', bbox_to_anchor=(1, 1))

        self.ax[1].scatter(self.times, self.positions['x'].values, s=self.dsize,label = 'x position (nm)')
        self.ax[1].scatter(self.times, self.positions['y'].values, color='red',s=self.dsize,label = 'y position (nm)')
        self.ax[1].scatter(self.times, self.positions['z'].values, color='green',s=self.dsize,label = 'z position (nm)') 
        self.ax[1].set_xlim(0,self.totalt)
        self.ax[1].set_ylim(self.positions.min().min(),self.positions.max().max())
        self.ax[1].set_title('Position vs time')
        self.ax[1].set_xlabel("Time (ns)")
        self.ax[1].set_ylabel("Position (nm)")
        self.ax[1].legend(loc='upper right', bbox_to_anchor=(1, 1))

        plt.tight_layout(pad=2.0)
        # Display the plot

        plt.savefig( self.dir + self.name + '/' +  "rt_msd.png")

    def graph_tht_msad(self):
        closeOpen()
        mean = np.sum(self.squared_displacement['Angular'])/self.steps
        self.fig, self.ax = plt.subplots(nrows=2,figsize=(10,12))
        counts, bins, patches = self.ax[0].hist(self.squared_displacement['Angular'],bins=self.nbins)
        threshold = self.steps * 5* 10**-4
        counts = np.where(counts<threshold,0,counts)
        lastbin = 0
        for i in range(0,len(counts)):
            if counts[i] == 0:
                lastbin = i
                break
            
        self.ax[0].cla()
        self.ax[0].bar(bins[:-1],counts,width = np.diff(bins)) 
        self.ax[0].set_xlim(0, bins[lastbin-1])
        self.ax[0].set_title('Distribution of squared angular displacment')
        self.ax[0].set_xlabel("Displacement squared in a single timestep (rad^2)")
        self.ax[0].set_ylabel("Frequency")
        self.ax[0].plot([mean]*2,[*self.ax[0].get_ylim()],label = f"mean squared angular displacement ={mean:.4e}")
        self.ax[0].plot([self.expmsad]*2,[*self.ax[0].get_ylim()],label = f'2*Dr*dt = {self.expmsad:4e}') # \n z-score = {np.sqrt(self.steps/2)*(mean/self.expmsad-1):3f}
        self.ax[0].legend(loc='upper right', bbox_to_anchor=(1, 1))


        self.ax[1].scatter(self.times, self.thetaz['theta'].values, s=self.dsize,label = 'yaw angle (rad)')
        self.ax[1].scatter(self.times, self.thetaz['z'].values, color='red',s=self.dsize,label = '\u1E91')
        self.ax[1].set_xlim(0,self.totalt)
        self.ax[1].set_ylim(self.thetaz.min().min(),self.thetaz.max().max())
        self.ax[1].set_title('Angular position vs time')
        self.ax[1].set_xlabel("Time (ns)")
        self.ax[1].set_ylabel("Position (nm)")
        self.ax[1].legend(loc='upper right', bbox_to_anchor=(1, 1))

        plt.tight_layout(pad=2.0)
    
        # Display the plot
        plt.savefig( self.dir + self.name + '/' +  "tht_msad_graph.png")
    def graphPositions(self):
        n  = int(self.steps/self.res)
        self.fig = plt.figure(figsize=(8, 6))
        self.ax = self.fig.add_subplot(111,projection='3d')
        sc = self.ax.scatter(self.positions['x'].values[::n],self.positions['y'].values[::n],self.positions['z'].values[::n],c=self.times[::n]/10**9,cmap='jet')
        # self.ax.quiver(positions['x'].values[::n],positions['y'].values[::n],positions['z'].values[::n],velocities['x'].values[::n],velocities['y'].values[::n],velocities['y'].values[::n],length=10,normalize=True)
        self.ax.set_title('Particle Brownian Motion Positions in 3D')
        self.ax.set_xlabel('x (nm)')
        self.ax.set_ylabel('y (nm)')
        self.ax.set_zlabel('z (nm)')
        boxsize = self.positions.abs().max().max()*2
        self.ax.set_xlim(-boxsize/2,boxsize/2)
        self.ax.set_ylim(-boxsize/2,boxsize/2)
        self.ax.set_zlim(-boxsize/2,boxsize/2)
        self.ax.set_box_aspect([1, 1, 1])
        # Adjust the layout to make room for the slider
        self.fig.subplots_adjust(bottom=0.16)

        # Create the slider axes (place it below the plot)
        range_slider_ax = self.fig.add_axes([0.2, 0.05, 0.65, 0.03])

        # Create the range slider
        time_slider = RangeSlider(range_slider_ax, 'Time', 0, self.totalt, valinit=(0, self.totalt))
        self.cbar = plt.colorbar(sc, ax=self.ax, shrink=0.6)
        self.cbar.set_label('color map (s)')
        # Update function for the slider
        def update(val):
            # Get the current slider values
            min_time, max_time = time_slider.val
            
            # Convert time to indices
            min_index = int(min_time / (self.dt*self.datadt))
            max_index = int(max_time / (self.dt*self.datadt))
            
            # Update the scatter plot with the new range
            if self.cbar:
                self.cbar.remove()
            self.ax.clear()
            sc = self.ax.scatter(self.positions.iloc[min_index:max_index]['x'].values[::n], 
                    self.positions.iloc[min_index:max_index]['y'].values[::n], 
                    self.positions.iloc[min_index:max_index]['z'].values[::n], 
                    c=self.times[min_index:max_index:n]/10**9, cmap='jet')
            
            # Redraw the plot
            self.ax.set_title('Particle Brownian Motion in 3D')
            self.ax.set_xlabel('x (nm)')
            self.ax.set_ylabel('y (nm)')
            self.ax.set_zlabel('z (nm)')
            self.ax.set_xlim(-boxsize / 2, boxsize / 2)
            self.ax.set_ylim(-boxsize / 2, boxsize / 2)
            self.ax.set_zlim(-boxsize / 2, boxsize / 2)
            self.fig.canvas.draw_idle()
            self.cbar = plt.colorbar(sc, ax=self.ax, shrink=0.6)
            self.cbar.set_label('color map (s)')

        # Connect the slider to the update function
        time_slider.on_changed(update)

        # Display the plot
        if self.mode == 'show':
            plt.show()
            return 
        #elif self.mode == 'save':
            #print(type(self.fig))
            #for ax in self.fig.axes:
                #print(type(ax), ax.get_frame_on)
            #plotlyfig = tls.mpl_to_plotly(self.fig)
            #pio.write_html(plotlyfig, file = self.dir + "3d_positions.html")
    def graphRod(self):
        n  = int(self.steps/self.res)
        self.fig = plt.figure(figsize=(8, 6))
        self.ax = self.fig.add_subplot(111,projection='3d')
        lowertip = self.positions - self.orientations
        #print(lowertip)
        uppertip = 2*self.orientations
        print(lowertip+uppertip)
        sc = self.ax.quiver(lowertip['x'].values[::n],lowertip['y'].values[::n],lowertip['z'].values[::n],uppertip['x'].values[::n],uppertip['y'].values[::n],uppertip['z'].values[::n],color =cm.jet(self.times[::n]/self.times[-1]),arrow_length_ratio=0.1)
        # plt.quiver(positions['x'].values[::n],positions['y'].values[::n],positions['z'].values[::n],velocities['x'].values[::n],velocities['y'].values[::n],velocities['y'].values[::n],length=10,normalize=True)
        self.ax.set_title('Particle Orientation in 3D')
        self.ax.set_xlabel('x (nm)')
        self.ax.set_ylabel('y (nm)')
        self.ax.set_zlabel('z (nm)')
        minpoint = abs(lowertip.min().min())
        mintip = abs((lowertip+uppertip).max().max())

        maxpoint = abs(lowertip.max().max())
        maxtip = abs((lowertip+uppertip).max().max())
        l = max(minpoint,mintip,maxpoint,maxtip)*1.1

        self.ax.set_xlim(-l,l)
        self.ax.set_ylim(-l,l)
        self.ax.set_zlim(-l,l)

        self.ax.set_box_aspect([1, 1, 1])
        self.fig.subplots_adjust(bottom=0.16)

        # Create the slider axes (place it below the plot)
        range_slider_ax = self.fig.add_axes([0.2, 0.05, 0.65, 0.03])

        # Create the range slider
        time_slider = RangeSlider(range_slider_ax, 'Time', 0, self.totalt, valinit=(0, self.totalt))
        self.cbar = plt.colorbar(sc, ax=self.ax, shrink=0.6)
        self.cbar.set_label('color map (s)')
        # Update function for the slider
        def update(val):
            # Get the current slider values
            min_time, max_time = time_slider.val
            
            # Convert time to indices
            min_index = int(min_time /  (self.dt*self.datadt))
            max_index = int(max_time / (self.dt*self.datadt))
            
            # Update the scatter plot with the new range
            if self.cbar:
                self.cbar.remove()
            self.ax.clear()
            sc = self.ax.quiver(lowertip.iloc[min_index:max_index]['x'].values[::n], 
                    lowertip.iloc[min_index:max_index]['y'].values[::n], 
                    lowertip.iloc[min_index:max_index]['z'].values[::n],uppertip.iloc[min_index:max_index]['x'].values[::n], 
                    uppertip.iloc[min_index:max_index]['y'].values[::n], 
                    uppertip.iloc[min_index:max_index]['z'].values[::n], 
                    color=cm.jet(self.times[min_index:max_index:n]/self.times[-1]),arrow_length_ratio=0.1)
            
            # Redraw the plot
            self.ax.set_title('Particle Orientation in 3D')
            self.ax.set_xlabel('x (nm)')
            self.ax.set_ylabel('y (nm)')
            self.ax.set_zlabel('z (nm)')
            self.ax.set_xlim(-l,l)
            self.ax.set_ylim(-l,l)
            self.ax.set_zlim(-l,l)

            self.cbar = plt.colorbar(sc, ax=self.ax, shrink=0.6)
            self.cbar.set_label('color map (s)')

            self.fig.canvas.draw_idle()

        # Connect the slider to the update function
        time_slider.on_changed(update)
        # Display the plot
        if self.mode == 'show':
            plt.show()


        #elif self.mode == 'save':
            #pass
            #plotlyfig = tls.mpl_to_plotly(self.fig)
            #pio.write_html(plotlyfig, file = self.dir + "3d_orientations.html")
    def graphOrientations(self):
        n  = int(self.steps/self.res)
        self.fig = plt.figure(figsize=(8, 6))
        self.ax = self.fig.add_subplot(111,projection='3d')
        sc = self.ax.scatter(self.orientations['x'].values[::n],self.orientations['y'].values[::n],self.orientations['z'].values[::n],c=self.times[::n]/10**9,cmap='jet')
        # plt.quiver(positions['x'].values[::n],positions['y'].values[::n],positions['z'].values[::n],velocities['x'].values[::n],velocities['y'].values[::n],velocities['y'].values[::n],length=10,normalize=True)
        self.ax.set_title('Particle Orientation in 3D')
        self.ax.set_xlabel('x (nm)')
        self.ax.set_ylabel('y (nm)')
        self.ax.set_zlabel('z (nm)')
        l = self.l * 10**9
        self.ax.set_xlim(-l/2,l/2)
        self.ax.set_ylim(-l/2,l/2)
        self.ax.set_zlim(-l/2,l/2)
        self.ax.set_box_aspect([1, 1, 1])
        self.fig.subplots_adjust(bottom=0.16)

        # Create the slider axes (place it below the plot)
        range_slider_ax = self.fig.add_axes([0.2, 0.05, 0.65, 0.03])

        # Create the range slider
        time_slider = RangeSlider(range_slider_ax, 'Time', 0, self.totalt, valinit=(0, self.totalt))
        self.cbar = plt.colorbar(sc, ax=self.ax, shrink=0.6)
        self.cbar.set_label('color map (s)')
        # Update function for the slider
        def update(val):
            # Get the current slider values
            min_time, max_time = time_slider.val
            
            # Convert time to indices
            min_index = int(min_time /  (self.dt*self.datadt))
            max_index = int(max_time / (self.dt*self.datadt))
            
            # Update the scatter plot with the new range
            if self.cbar:
                self.cbar.remove()
            self.ax.clear()
            sc = self.ax.scatter(self.orientations.iloc[min_index:max_index]['x'].values[::n], 
                    self.orientations.iloc[min_index:max_index]['y'].values[::n], 
                    self.orientations.iloc[min_index:max_index]['z'].values[::n], 
                    c=self.times[min_index:max_index:n]/10**9, cmap='jet')
            
            # Redraw the plot
            self.ax.set_title('Particle Orientation in 3D')
            self.ax.set_xlabel('x (nm)')
            self.ax.set_ylabel('y (nm)')
            self.ax.set_zlabel('z (nm)')
            self.ax.set_xlim(-l / 2, l / 2)
            self.ax.set_ylim(-l / 2, l / 2)
            self.ax.set_zlim(-l / 2, l / 2)

            self.cbar = plt.colorbar(sc, ax=self.ax, shrink=0.6)
            self.cbar.set_label('color map (s)')

            self.fig.canvas.draw_idle()

        # Connect the slider to the update function
        time_slider.on_changed(update)
        # Display the plot
        if self.mode == 'show':
            plt.show()


        #elif self.mode == 'save':
            #pass
            #plotlyfig = tls.mpl_to_plotly(self.fig)
            #pio.write_html(plotlyfig, file = self.dir + "3d_orientations.html")
    def saveSim(self):
        self.datadict = {
            'orientations': self.orientations,
            'positions': self.positions,
            'thetaz': self.thetaz,
            'velocities': self.velocities,
            'omegas': self.omegas,
            'omegasobjframe': self.omegasobjframe,
            'squared_displacement': self.squared_displacement
        }

        for file_name, data_object in self.datadict.items():
            # Save to CSV
            data_object.to_hdf(self.dir + self.name + '/' +f'{self.name}.h5', key = file_name, mode = 'a', format = 'table', complib='zlib', complevel=9)
            # Generate a progress message based on the file name
            progress_message = f"{file_name} saved"
            self.updateProgress(progress_message)
    def readData(self):
        # Define the mapping between file names and class attributes
        file_to_attribute = {
            'orientations': 'orientations',
            'positions': 'positions',
            'thetaz': 'thetaz',
            'velocities': 'velocities',
            'omegas': 'omegas',
           'omegasobjframe': 'omegasobjframe',
            'squared_displacement': 'squared_displacement'
        }

        self.datadict = {}  # Initialize an empty dictionary to store the updated DataFrames
        #print('test')
        for file_name, attribute_name in file_to_attribute.items():
            # Read the data from the CSV
 
            data = pd.read_hdf(self.dir + self.name + '/' +f'{self.name}.h5', key = file_name)
            self.updateProgress(f'{file_name} read')
            #print(data.head())

            # Update both the dictionary and the corresponding attribute
            setattr(self, attribute_name, data)
            self.datadict[file_name] = data
        self.times *= 10**9
        self.totalt *= 10**9
        self.dt *= 10**9 

            
    def updateProgress(self,msg):
        with open(self.dir + self.name + "/progress.txt",'a') as file:
            file.write('\n' + msg) 
class SphereSim(Simulation):
    def __init__(self,dt,datadt,steps,temp,rho,mu,d,dsize,nbins,mode,dir="",res = None, name = ''):
        super().__init__(dt,datadt,steps,temp,rho,mu,d,dsize,nbins,mode,dir,res,name)
        self.type = 'sph'
        self.particleData()
        super().makeFolder()
    def next_data(self):
        # translational kinematics
        randomr = np.random.normal(loc=0,scale=1,size=3) @ self.rvar**.5
        #print(f"randomr1:{randomr}")
        randomv = np.random.normal(loc=0, scale=1, size=3) @ (self.vvar@(np.eye(3)-self.rvcorr**2))**.5 + randomr @ self.rvcorr @ (self.vvar @ inv(self.rvar)) ** .5
        #print(f"randomv:{randomv}")
        randomr = np.eye(3) @ randomr
        #print(f"randomr2:{randomr}")
        #print(f"rmatrix:{self.rmatrix}")
        randomr = self.rmatrix @ randomr
        #print(f"randomr3:{randomr}")
        randomv = self.rmatrix @ randomv
        #print(f"randomv2:{randomv}")
        dr = randomr + (np.eye(3)-e(-self.bdt)) @ inv(self.beta0) @ self.iv
        #print(dr)
        self.sqdr = dr @ dr
        self.ir += dr
        self.iv = e(-self.bdt) @ self.iv + randomv
        
        # rotational kinematics 
        randomth = np.random.normal(loc=0,scale=1,size=3) @ self.thvar**.5
        #print(f"randomth2:{randomth}")
        randomw = np.random.normal(loc=0, scale=1, size=3) @ (self.wvar@(np.eye(3)-self.thwcorr**2))**.5 + randomth @ self.thwcorr @ (self.wvar @ inv(self.thvar)) ** .5
        randomth = self.rmatrix @ randomth
        #print(f"randomth2:{randomth}")
        randomw = self.rmatrix @ randomw
        dth = randomth + (np.eye(3)-e(-self.brdt)) @ inv(self.betarot0) @ self.iw
        #print(f"dth-randomth:{dth-randomth}")
        #print(f"dth:{dth}")
        self.sqdth = dth @ dth
        self.iw = e(-self.brdt) @ self.iw + randomw
        self.rmatrix = rotateMatrix(dth,self.rmatrix)
        self.ior = self.rmatrix @ self.upright
        #self.ith = [theta(self.ior[0],self.ior[1]),self.ior[2]/np.sqrt(self.upright@self.upright)]

    def particleData(self):
        self.upright = np.array([0,0,self.d/2])
        self.ior = self.upright
        self.orientations = np.empty((self.steps+1,3))
        self.orientations[0] = self.upright.copy()

        self.vol = np.pi*self.d**3/6
        m = self.vol*self.rho
        self.mmat = m*np.eye(3)
        self.imat = m *self.d**2/10 * np.eye(3)
        self.kbTm = kb*self.temp*inv(self.mmat)
        self.kbTI = kb*self.temp*inv(self.imat)
        self.beta0 = inv(self.mmat)*3*self.d*self.mu*np.pi
        self.betarot0 = inv(self.imat)*self.d**3*np.pi*self.mu
        self.diffusion = np.trace(self.kbTm @ inv(self.beta0))
        self.rotdiffusion = np.trace(self.kbTI @ inv(self.betarot0))
        self.expmsd = self.diffusion*2*self.dt
        self.expmsad = self.rotdiffusion*2*self.dt
        self.rmatrix = np.eye(3)
        self.spring = 1 # spring constant
        #translational kinematics data
        self.bdt = self.beta0*self.dt
        #print(self.bdt)
        self.rvar = self.dt**2 * self.kbTm @ (2*inv(self.bdt) + inv(self.bdt)**2 @ (-3*np.eye(3)+4*e(-self.bdt)-e(-2*self.bdt)))
        #print(self.rvar)
        self.vvar = self.kbTm @ (np.eye(3)-e(-2*self.bdt))
        #print(self.vvar)
        self.rvcorr =  (inv(self.vvar @ self.rvar)**.5) @ self.kbTm @ inv(self.beta0) @ ((np.eye(3)-e(-1*self.bdt))**2)

        #rotational kinematics data
        self.brdt = self.betarot0 * self.dt 
        self.thvar = self.dt**2 * self.kbTI @ (2*inv(self.brdt) + inv(self.brdt)**2 @ (-3*np.eye(3)+4*e(-self.brdt)-e(-2*self.brdt)))
        self.wvar = self.kbTI @ (np.eye(3)-e(-2*self.brdt)) 
        self.thwcorr = (inv(self.wvar @ self.thvar)**.5) @ self.kbTI @ inv(self.betarot0) @ ((np.eye(3)-e(-1*self.brdt))**2)

class CylinderSim(Simulation):
    def __init__(self,dt,datadt,steps,temp,rho,mu,d,l,dsize,nbins,mode,dir="",res = None,name=''):
        super().__init__(dt,datadt,steps,temp,rho,mu,d,dsize,nbins,mode,dir,res,name)
        self.type = 'cyl'
        self.l = l
        self.particleData()
        super().makeFolder()
    def particleData(self):
        self.upright = np.array([0,0,self.l/2])
        self.ior = self.upright
        self.orientations = np.empty((self.steps+1,3))
        self.orientations[0] = self.upright.copy()
        self.vol = np.pi*self.d**2/4*self.l
        m = self.vol * self.rho # mass of diameter spherical water particle
        self.mmat = np.eye(3)*m
        ivert = 0.5*m*(self.d/2)**2 # rotational inertia around the long axis
        iplanar = self.rho*np.pi*self.l*self.d**2*(self.l**2/3+self.d**2/4)/16 # rotational inertia around the plane axis
        self.imat = np.diag([iplanar,iplanar,ivert])
        self.kbTm = kb*self.temp*inv(self.mmat)
        self.kbTI = kb*self.temp*inv(self.imat)
        gamma = 2*np.pi*self.mu*self.l
        zparallel = gamma/((log(self.l/self.d)-0.2)*m) # drag coefficient for movement along the axis ('z' is a stand-in for zeta)
        znormal = 2*gamma/((log(self.l/self.d)+0.84)*m) #  drag coefficient for movement normal to the axis
        zrotnormal = gamma*self.l**2/(iplanar*(log(self.l/self.d)-0.66)) # drag coefficient for rotation around the center of the rod lengthwise
        zaxis = gamma*self.d**2/(2*ivert) # drag coefficient for rotation around the axis
        self.beta0 = np.diag([znormal,znormal,zparallel])
        self.betarot0 = np.diag([zrotnormal,zrotnormal,zaxis])
        self.diffusion = np.trace(self.kbTm @ inv(self.beta0))
        self.rotdiffusion = np.trace(self.kbTI @ inv(self.betarot0))
        self.expmsd = self.diffusion*2*self.dt
        self.expmsad = self.rotdiffusion*2*self.dt
        self.rmatrix = np.eye(3)
        self.spring = 0.8 # spring constant
        self.ligand = np.array([0,0,self.l/2])
        # self.ligands = np.array([0,0,self.l/2],[0,0,-self.l/2],[0,self.d/2,0],[0,-self.d/2,0],[self.d/2,0,0],[-self.d/2,0,0]) # ligand positions relative to CM
        self.bonded = [True,False,False,False,False,False] # boolean for if the bond is attached

        self.equilibrium = self.ligand.copy()
        self.iligand = self.ligand.copy()
        self.lambda_eq = 41.1 * 10**-9 # equilibrium bond length
        self.lambda_t  = self.lambda_eq # bond length, initalized at the equilibrium bond length
        self.rl = np.array([0,0,self.l/2+self.lambda_eq])
        self.gamma_bond = 0.274 * 10**-9 # bond reactive compliance 0.274nm
        self.continueSim = True # boolean for if the bond is attached
        self.attached = True # boolean for if the bond is attached
        self.kr0 = 1.1 * 10**-4 # intrinsic bond breakage rate
        #translational kinematics data
        self.bdt = self.beta0*self.dt
        #print(self.bdt)
        self.rvar = self.dt**2 * self.kbTm @ (2*inv(self.bdt) + inv(self.bdt)**2 @ (-3*np.eye(3)+4*e(-self.bdt)-e(-2*self.bdt)))
        #print(self.rvar)
        self.vvar = self.kbTm @ (np.eye(3)-e(-2*self.bdt))
        #print(self.vvar)
        self.rvcorr =  (inv(self.vvar @ self.rvar)**.5) @ self.kbTm @ inv(self.beta0) @ ((np.eye(3)-e(-1*self.bdt))**2)
        #rotational kinematics data
        self.brdt = self.betarot0 * self.dt 
        self.thvar = self.dt**2 * self.kbTI @ (2*inv(self.brdt) + inv(self.brdt)**2 @ (-3*np.eye(3)+4*e(-self.brdt)-e(-self.brdt)))
        self.wvar = self.kbTI @ (np.eye(3)-e(-2*self.brdt)) 
        self.thwcorr = (inv(self.wvar @ self.thvar)**.5) @ self.kbTI @ inv(self.betarot0) @ ((np.eye(3)-e(-self.brdt))**2)
        receptorlength = 18.7*10**-9 # 18.7nm length of the receptor
        wallnorm = -np.linalg.norm(self.ligand)
        wallpoint = self.equilibrium - wallnorm*receptorlength
        self.wall = wallnorm, wallpoint



    def __str__(self):
        return super().__str__()+f'\nlength: {self.l*10**9} nm' + f'\nligand position relative to CM: {self.ligand}' + f'\nequilibrium: {self.equilibrium}' + f'\nspring constant: {self.spring} N/m' 


simName = '051425'
simDir = f"Simulations{os.sep}"
simMode = 'save'
cylindersim = CylinderSim(dt=10**-9,datadt=10**3, steps=10**6,temp=298.15,rho=997,mu=8.8891*(10**-4),d=5*10**-8,l=5*10**-7,dsize=10,nbins=50,mode=simMode,res = 10**3, name = simName,dir ="Simulations") 
#spheresim = SphereSim(dt=10**-5,steps= 10**6,temp=298.15,rho=997,mu=8.8891*(10**-4),d=6*10**-8,dsize=10,nbins=50,mode=simMode,res = 10**5, name = simName,dir ="Simulations") 
#readcylindersim = CylinderSim(dt=10**-5,steps=5*10**5,temp=298.15,rho=997,mu=8.8891*(10**-4),d=5*10**-8,l=5*10**-7,dsize=10,nbins=50,mode='show',res = 10**5, name = simName,dir ="Simulations") 
#readspheresim = SphereSim(dt=10**-5,steps= 5* 10**5,temp=298.15,rho=997,mu=8.8891*(10**-4),d=6*10**-8,dsize=10,nbins=50,mode='show',res = 10**5, name = simName,dir ="Simulations") 

#jobSim.mode = 'save'
# sim = SphereSim(10**-6,10**5,298.15,997,8.9*(10**-4),10**-7,10,30, res= 10**4)sim.generateData()
def runSim(sim):
    sim.generateData()
    sim.updateProgress("data generated")
    sim.saveSim()
    ''''''
#runSim(spheresim)
#runSim(cylindersim)
def graphSim(sim):
    sim.readData()
    sim.updateProgress('data read')
    #print(sim.positions)
    sim.graphPositions()
    #print(sim.ori)
    sim.graphRod()
    sim.graphOrientations()
    #sim.graph_rt_msd()
    #sim.updateProgress('rtmsd graph saved')
    sim.graphthw()
    sim.updateProgress('thw graph saved')
    sim.graphxv()
    sim.updateProgress('xv graph saved')
    #sim.graph_tht_msad()
    #sim.updateProgress('thtmsad graph saved')
def main(sim):
    if sim.mode == 'save':
        runSim(sim)
    elif sim.mode == 'show':
        graphSim(sim)


'''spheresim.readData()
spheresim.saveSim()
cylindersim.readData()
cylindersim.saveSim()'''

#main(spheresim)
main(cylindersim)
# main(jobSimcyl)
#main(cylindersim)
elapsed = time.time()- start
cylindersim.updateProgress(f"Elapsed time: {elapsed:.4f} seconds")
