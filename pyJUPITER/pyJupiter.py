from numpy import *
import matplotlib
#matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

class output: #name of the class

    def __init__(self, N):
        '''
        init method with N as the argument
        '''
        self.N = N
        self.folder = 'output%.5d'%N #name of the output folder outputXXXXX

        #read the descriptor file
        with open(self.folder+'/Descriptor%.1d.dat'%N,'r') as f:
             self.lines = f.readlines()

        #highest refinement level
        self.maxreflvl=int(self.lines[-3][-2:-1]) #number of grid refinement levels
        #create empty containers for grid info
        self.xdim = empty(self.maxreflvl+1,dtype=int) #numper of cells along first dimension (typically azimuth)
        self.ydim = empty(self.maxreflvl+1,dtype=int) #numper of cells along second dimension (typically radius)
        self.zdim = empty(self.maxreflvl+1,dtype=int) #numper of cells along second dimension (typically polar)
        self.xlim = empty((self.maxreflvl+1,2)) #lower and upper boundaries of the first dimension [min,max]
        self.ylim = empty((self.maxreflvl+1,2)) #lower and upper boundaries of the second dimension [min,max]
        self.zlim = empty((self.maxreflvl+1,2)) #lower and upper boundaries of the third dimension [min,max]

        for lvl in range(self.maxreflvl+1): #here we loop over all the grid refinement levels.
            dim_line = lvl*11+6 #current refinement level
            #fill in the grid onvo of the current refinement level
            dim = self.lines[dim_line]
            self.xdim[lvl] = int(dim.split(" ")[0])
            self.xlim[lvl,0] = float(self.lines[dim_line+2].split(" ")[2])
            self.xlim[lvl,1] = float(self.lines[dim_line+2].split(" ")[-4])
            self.ydim[lvl] = int(dim.split(" ")[1])
            self.ylim[lvl,0] = float(self.lines[dim_line+3].split(" ")[2])
            self.ylim[lvl,1] = float(self.lines[dim_line+3].split(" ")[-4])
            self.zdim[lvl] = int(dim.split(" ")[2])
            self.zlim[lvl,0] = float(self.lines[dim_line+4].split(" ")[2])
            self.zlim[lvl,1] = float(self.lines[dim_line+4].split(" ")[-4])

        # check whether the simulation is 2D or 3D
        if self.zdim[lvl] == 1:
            self.dim=2
        else:
            self.dim=3

        #create some epmty instance variables
        self.field = "" #name of the field
        self.azimuth = 0.0
        self.polar = False
        self.velocityfield = False
        self.verticalplot = False
        self.midplaneplot = False

        self.cmap=plt.get_cmap('viridis') #standard colormap


    def readdata(self,field,reflvl):
        '''
        reads in .dat file and returns z*y*x dimensional array with simulation data
        '''
        #read .dat file
        data = fromfile(self.folder+'/%s%.1d_%.1d_%.1d.dat'%(field,self.N,reflvl,reflvl))
        #rearange
        xdim = self.xdim[reflvl]
        ydim = self.ydim[reflvl]
        zdim = self.zdim[reflvl]
        if (field[-8:]=='velocity'):
            self.velocityfield=True #the field is a vector field
        if self.velocityfield==True:
            fielddata = reshape(data,(3,zdim,ydim,xdim))
        else:
            fielddata = reshape(data,(zdim,ydim,xdim))

        if (field=='dustdensity'): # I chose a different colormap for dust
            self.cmap=plt.get_cmap('inferno')
        else:
            self.cmap=plt.get_cmap('viridis')

        return fielddata

    def upddate_maxmin(self,field,reflvls,vel=0,z=0,azigrid=0,avg=True,vert=False):
        '''
        computes the maximum and minimum values over all the refinement levels for the colormap
        '''
        self.midmax = 0.0
        self.vertmax = 0.0
        self.midmin = 1.0e10
        self.vertmin = 1.0e10
        for lvl in reflvls:
            if self.velocityfield==True:
                    middata = self.readdata(field,lvl)[vel,-1-z,:,:]
                    vertdata = self.readdata(field,lvl)[vel,:,:,int(self.xdim[lvl]/2)]
                    data = self.readdata(field,lvl)[vel,:,:]
            else:

                    middata = self.readdata(field,lvl)[-1-z,:,:]
                    vertdata = self.readdata(field,lvl)[:,:,int(self.xdim[lvl]/2)]
                    data = self.readdata(field,lvl)
            midmaximum = amax(middata)
            midminimum = amin(middata)
            vertmaximum = amax(vertdata)
            vertminimum = amin(vertdata)

            if midmaximum > self.midmax:
                self.midmax = midmaximum
            if midminimum < self.midmin:
                self.midmin = midminimum
            if vertmaximum > self.vertmax:
                self.vertmax = vertmaximum
            if vertminimum < self.vertmin:
                self.vertmin = vertminimum

            if (self.midmin<1e-50):
                self.midmin = 1e-50

            if vert==True:
                if avg==True:
                    self.vertmax = amax(data.mean(-1))
                    self.vertmin = amin(data.mean(-1))
                else:
                    self.vertmax = amax(data[:,:,azigrid])
                    self.vertmin = amin(data[:,:,azigrid])
                if (self.vertmin<1e-50):
                    self.vertmin = 1e-50
        return None


    def plot_setlimits(self,polar=False,vertical=False):
        if self.polar:
            if vertical:
                self.ax.set_thetamin(80)
                self.ax.set_thetamax(100)
                self.ax.set_rorigin(-0.2)
                self.ax.set_rmin(0.2)
                self.ax.set_rmax(2.6)
                self.ax.set_theta_zero_location('S')
            else:
                self.ax.set_theta_zero_location('E')

            self.ax.set_xlabel('r [au]')
            self.ax.set_ylabel('colatidute [deg]')
        else:
            if self.verticalplot:
                self.ax.set_ylabel('polar [rad]')
                self.ax.set_xlabel(r'r [$r_p$]')
            if self.midplaneplot:
                self.ax.set_ylabel(r'r [$r_p$]')
                self.ax.set_xlabel('azimuth [rad]')
        if vertical:
            None
            #self.ax.set_ylim(pi/2-0.01,pi/2+0.01)
            #self.ax.set_xlim(0.99,1.01)
            #self.ax.set_xlim(self.zlim[0][0],(self.zlim[0][1]-self.zlim[0][0])+self.zlim[0][1]+0.01)
            #self.ax.set_ylim(self.zlim[0][0],(self.zlim[0][1]-self.zlim[0][0])+self.zlim[0][1]+0.01)

        if self.azimuth != 0.0:
            self.ax.set_title(self.folder+' '+self.field+', azimuth=%.2f'%self.azimuth)
        else:
            self.ax.set_title(self.folder+' '+self.field+' ')

        return None


    def plotmidlvl(self,lvl):
        slcdata = self.slice
        xcsize = (self.xlim[lvl,1]-self.xlim[lvl,0])/self.xdim[lvl]
        ycsize = (self.ylim[lvl,1]-self.ylim[lvl,0])/self.ydim[lvl]
        xplot = linspace(self.xlim[lvl,0],self.xlim[lvl,1],self.xdim[lvl]+1)
        yplot = linspace(self.ylim[lvl,0],self.ylim[lvl,1],self.ydim[lvl]+1)
        xplot_2d, yplot_2d = meshgrid(xplot, yplot)
        if self.velocityfield==True:
            cx=self.ax.pcolormesh(xplot_2d,yplot_2d,slcdata,shading='flat',cmap=self.cmap)
        else:
            cx=self.ax.pcolormesh(xplot_2d,yplot_2d,slcdata,norm=colors.LogNorm(vmin=self.midmin, vmax=self.midmax),shading='flat',cmap=self.cmap)

        if lvl==0:
            self.fig.colorbar(cx)
        return None

    def plotmidslice(self,field='gasdensity',reflvls=-1,polar=False,vel=0,z=0):
        self.field=field
        self.polar = polar
        self.midplaneplot = True
        if (field[-8:]=='velocity'):
            self.velocityfield=True
        if reflvls==-1: #plot all the refinement levels
            reflvls = arange(0,self.maxreflvl+1)
        self.upddate_maxmin(field,reflvls,vel,z,avg=False,vert=False)

        if polar:
            figsz = (9,9)
        else:
            figsz=(9,5)
        self.fig=plt.figure(figsize=figsz)
        self.ax = self.fig.add_subplot(111,polar=polar)
        for lvl in reflvls:
            if self.velocityfield==True:
                self.slice = self.readdata(field,lvl)[vel,-1-z,:,:]
            else:
                self.slice = self.readdata(field,lvl)[-1-z,:,:] #midplane data of reflevel 'lvl'
            self.plotmidlvl(lvl) #plot midplane data of reflevel 'lvl'
        self.plot_setlimits(polar,vertical=False)
        plt.show(block=False)
        return None


    def plotvertlvl(self,lvl):
        slcdata = self.slice
        ycsize = (self.ylim[lvl,1]-self.ylim[lvl,0])/self.ydim[lvl]
        zcsize = (self.zlim[lvl,1]-self.zlim[lvl,0])/self.zdim[lvl]
        yplot = linspace(self.ylim[lvl,0],self.ylim[lvl,1],self.ydim[lvl]+1)
        zplot = linspace(self.zlim[lvl,0],self.zlim[lvl,1]+(self.zlim[lvl,1]-self.zlim[lvl,0]),self.zdim[lvl]*2+1)
        yplot_2d, zplot_2d = meshgrid(yplot, zplot)

        if self.polar:
            cont = transpose(zplot_2d)
            zplot_2d = transpose(yplot_2d)
            yplot_2d = cont
            slcdata = transpose(slcdata)
        if self.velocityfield==True:
            cx = self.ax.pcolormesh(yplot_2d,zplot_2d,slcdata,
                          shading='flat',cmap=self.cmap)
        else:
            cx = self.ax.pcolormesh(yplot_2d,zplot_2d,slcdata,norm=colors.LogNorm(vmin=self.vertmin, vmax=self.vertmax),
                          shading='flat',cmap=self.cmap,alpha=1.0)

        return cx

    def plotvertslice(self,field='gasdensity',reflvls=-1,polar=False,azi=0.0,vel=0):
        if self.dim!=3:
            print('no vertical slice for 2 dimensional array')
            return None
        self.field=field
        self.polar = polar
        self.azimuth = azi
        self.verticalplot = True
        if (field[-8:]=='velocity'):
            self.velocityfield=True
        if reflvls==-1:
            reflvls = arange(0,self.maxreflvl+1)
        self.upddate_maxmin(field,reflvls,vel,z=0,avg=False,vert=True)
        if polar:
            figsz = (9,9)
        else:
            figsz=(12,4)
        fig=plt.figure(figsize=figsz)
        self.ax = fig.add_subplot(111,polar=polar)
        if azi != 0.0:
            reflvls = arange(0,1);
        for lvl in reflvls:
            azigrid = int(self.xdim[lvl]/2)
            if (lvl==0):
                self.upddate_maxmin(field,[lvl],vel,z=0,azigrid=azigrid,avg=False,vert=True)
            if azi != 0.0:
                azigrid = int(self.xdim[lvl]*(azi + pi)/(2.*pi))
            if self.velocityfield==True:
                self.slice = self.readdata(field,lvl)[vel,:,:,azigrid]
            else:
                self.slice = self.readdata(field,lvl)[:,:,azigrid]
            self.slice = concatenate((self.slice,self.slice[::-1,:]),axis=0)
            cx = self.plotvertlvl(lvl) #plot midplane data of reflevel 'lvl'
        fig.colorbar(cx)
        self.plot_setlimits(polar,vertical=True)
        plt.show(block=False)
        return None

    def plotvertavg(self,field='gasdensity',polar=False):
        self.verticalplot = True
        if self.dim!=3:
            print('no vertical average for 2 dimensional array')
            return None
        self.field=field
        self.polar = polar
        lvl = ([0])
        self.upddate_maxmin(field,lvl,avg=True,vert=True)

        if polar:
            figsz = (9,9)
        else:
            figsz=(12,4)
        fig=plt.figure(figsize=figsz)
        self.ax = fig.add_subplot(111,polar=polar)

        self.slice = average(self.readdata(field,0),axis=2)
        self.slice = concatenate((self.slice,self.slice[::-1,:]),axis=0)
        cx = self.plotvertlvl(0) #plot midplane data of reflevel 'lvl'
        fig.colorbar(cx)
        self.plot_setlimits(polar,vertical=True)
        plt.show(block=False)
        return None

    def surfdens(self,field):
        lvl = 0
        zcsize = (self.zlim[lvl,1]-self.zlim[lvl,0])/self.zdim[lvl]
        y = zcsize*linspace(self.ylim[lvl,0],self.ylim[lvl,1],self.ydim[lvl])
        dz = reshape(tile(y,self.xdim[0]),(self.xdim[0],self.ydim[0]))
        dz = transpose(dz)
        fielddata = 2.*sum(self.readdata(field,0),axis=0)*dz
        self.surfdens = fielddata
        return fielddata

    def plotgapprofile(self,field='gasdensity',ylog=True,avg=False,azi=0.0):
        lvl=0
        self.azimuth=azi
        self.field=field
        if avg:
            if self.dim == 2:
                data = average(self.readdata(field,0)[0,:,:],axis=1)
            if self.dim == 3:
                data = average(self.surfdens(field),axis=1)
        else:
            azigrid = int(self.xdim[lvl]*(pi+ self.azimuth)/(2.*pi))



            if self.dim == 2:
                data = self.readdata(field,0)[0,:,azigrid]
            if self.dim == 3:
                data = self.surfdens(field)[:,azigrid]

        yplot = linspace(self.ylim[lvl,0],self.ylim[lvl,1],self.ydim[lvl])
        plt.figure(figsize=(8,4))
        plt.plot(yplot,data,color = 'k')
        if ylog:
            plt.yscale('log')
        plt.ylabel(r'$\Sigma$ [$M_{*}/au^2$]')
        plt.xlabel('r [au]')
        if avg:
            plt.title('surafce density profile (azim avg), '+self.folder+' '+self.field+'')
        else:
            plt.title('surafce density profile, azimuth=%.2f'%self.azimuth)
        plt.show(block=False)
        return None
