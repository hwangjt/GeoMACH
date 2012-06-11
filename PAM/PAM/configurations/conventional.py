from __future__ import division
import numpy
import sys
from PAM.components import fullplate, halfbody, fulljunction, halfjunction, fuse_sections
from PAM.configurations import configuration


class conventional(configuration):

    def __init__(self):
        fuse = halfbody([70,10,10,20,10,10,10,50,10,25,10,10],[25,25,25,25],[15])
        fuse.translatePoints(0,0,0)
        wing = fullplate([10,10,10,50],[10,20,10,10])
        wing.translatePoints(4,0,0)
        tail = fullplate([30],[25])
        tail.translatePoints(8,0,0)
        nacelle = halfbody([50,10,20,30],[30],[20,10,10,20],full=True)
        nacelle.translatePoints(12,0,0)
        pylon = fullplate([30],[20],opentip=True)
        pylon.translatePoints(16,0,0)
        fin = fullplate([30],[25],half=True)
        fin.translatePoints(20,0,0)
        wingfuse = fulljunction(wing, 0, fuse, 2, [2,1], [3,6])
        tailfuse = fulljunction(tail, 0, fuse, 2, [1,8], [2,10])
        pylonnacelle = fulljunction(pylon, 0, nacelle, 1, [1,1], [2,3])
        pylonwing = fulljunction(pylon, 1, wing, 1, [2,1], [0,2])
        finfuse = halfjunction(fin, 0, fuse, 1, [0,8], [0,10])

        self.components = []
        self.components.append(fuse)
        self.components.append(wing)
        self.components.append(tail)
        self.components.append(nacelle)
        self.components.append(pylon)
        self.components.append(fin)
        self.components.append(wingfuse)
        self.components.append(tailfuse)
        self.components.append(pylonnacelle)
        self.components.append(pylonwing)
        self.components.append(finfuse)

        self.assembleComponents()

        self.components[0].setSections(2,fuse_sections.rounded2)
        self.components[0].setSections(3,fuse_sections.rounded2)
        self.components[0].setSections(4,fuse_sections.rounded2)
        self.components[0].setSections(5,fuse_sections.rounded2)
        self.components[0].props['posx'].set([0,10],[0,1])
        self.components[0].props['posy'].set([0.3,0.5,0.5],[0,0.15,1],w=[1.0,0,0],d=[1,0,0])
        self.components[0].props['ry'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])
        self.components[0].props['rz'].set([0.1,0.5,0.5,0.1],[0,0.15,0.75,1.0],w=[0.9985,0,0,0],d=[1,0,0,0])

        self.components[1].offset[:] = [3.75, 0.3, 0.5]
        self.components[1].setAirfoil("rae2822.dat")
        self.components[1].props['posx'].set([0,3.2,4],[0,0.8,1],w=[0.4,1,0])
        self.components[1].props['posy'].set([0,0.9,2.1],[0,0.8,1],w=[0.5,1,0])
        self.components[1].props['posz'].set([0,4.5,5],[0,0.8,1],w=[0,1,0])
        self.components[1].props['prpx'].set([1,1],[0,1])
        self.components[1].props['prpy'].set([0,0],[0,1])
        self.components[1].props['chord'].set([2,0.25],[0,1])

        self.components[2].offset[:] = [8.5, 0.5, 0.35]
        self.components[2].props['posx'].set([0,1.6],[0,1],w=[0.2,0])
        self.components[2].props['posy'].set([0,0.3],[0,1],w=[0,0])
        self.components[2].props['posz'].set([0,1.7],[0,1])
        self.components[2].props['roty'].set([10,0],[0,1])
        self.components[2].props['prpx'].set([0,0],[0,1])
        self.components[2].props['prpy'].set([0,0],[0,1])
        self.components[2].props['chord'].set([0.85,0.15],[0,1])

        e = numpy.zeros((11,4))
        e[0,:] = [0.20, 0.03, 1, 1]
        e[1,:] = [0.30, 0.10, 0, 0]
        e[2,:] = [0.30, 0.29, 0.5, 0]
        e[3,:] = [0.00, 0.29, 1, 1]
        e[4,:] = [0.00, 0.31, 1, 1]
        e[5,:] = [0.50, 0.35, 1, 0]
        e[6,:] = [1.00, 0.30, 0, 0]
        e[7,:] = [1.00, 0.27, 0, 0]
        e[8,:] = [0.65, 0.27, 0, 0]
        e[9,:] = [0.65, 0.20, 1, 0]
        e[10,:] = [1.2, 0.05, 0, 0]
        e[:,:2] *= 0.8
        l = numpy.linspace(0,1,e.shape[0])
        self.components[3].offset[:] = [3.3, -0.05, 1.5]
        self.components[3].props['posx'].set(e[:,0],l)
        self.components[3].props['posy'].set([0,0],[0,1])
        self.components[3].props['ry'].set(e[:,1],l,e[:,2],e[:,3])
        self.components[3].props['rz'].set(e[:,1],l,e[:,2],e[:,3])

        self.components[4].offset[:] = [3.8, 0.2, 1.5]
        self.components[4].setAirfoil("naca0015")
        self.components[4].props['posx'].set([0,0.2],[0,1])
        self.components[4].props['posy'].set([0,0.08],[0,1])
        self.components[4].props['posz'].set([0,0],[0,1])
        self.components[4].props['prpx'].set([1,1],[0,1])
        self.components[4].props['prpy'].set([0,0],[0,1])
        self.components[4].props['chord'].set([0.75,0.75],[0,1])

        self.components[5].offset[:] = [8.5, 0.85, 0]
        self.components[5].props['posx'].set([0,1.5],[0,1])
        self.components[5].props['posy'].set([0,1.4],[0,1])
        self.components[5].props['posz'].set([0,0],[0,1])
        self.components[5].props['rotz'].set([10,0],[0,1])
        self.components[5].props['prpx'].set([1,1],[0,1])
        self.components[5].props['prpy'].set([0,0],[0,1])
        self.components[5].props['chord'].set([1,0.2],[0,1])

        self.computePoints()

if __name__ == '__main__':

    aircraft = conventional()
    aircraft.oml0.write2Tec('conventional')
    aircraft.oml0.write2TecC('conventionalC')
    aircraft.plot()
