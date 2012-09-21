from __future__ import division
import numpy, time
import sys
from PAM.components import Wing, Body, FullInterface, HalfInterface
from PAM.configurations import Configuration


class Supersonic(Configuration):

    def __init__(self):
        super(Supersonic,self).__init__() 

        self.addComp('fuse', Body([90,10,160,10,10,35,10,10],[40,40,10,10],[15]))
        self.addComp('wing', Wing([90],[160]))
        self.addComp('tail', Wing([30],[35]))
        self.addComp('nacelle', Body([50,10,20,30],[31],[20,10,10,20],full=True))
        #self.addComp('pylon', Wing([30],[20],opentip=True))
        #self.addComp('fin', Wing([30],[25],half=True))

        self.separateComps()

        self.addComp('wingfuse', FullInterface(self.comps, 'wing', 0, 'fuse', 2, [0,1], [1,3]))
        self.addComp('tailfuse', FullInterface(self.comps, 'tail', 0, 'fuse', 2, [0,4], [1,6]))
        #self.addComp('pylonnacelle', FullInterface(self.comps, 'pylon', 0, 'nacelle', 1, [1,1], [2,3]))
        #self.addComp('pylonwing', FullInterface(self.comps, 'pylon', 1, 'wing', 1, [2,1], [0,2]))
        #self.addComp('finfuse', HalfInterface(self.comps, 'fin', 0, 'fuse', 1, [0,8], [0,10]))

        self.assembleComponents()

        c = self.comps
        #c['fuse'].setSections(sections=[2,3,4,5], t1L=0.35, t2L=0.65)
        c['fuse'].props['posx'].set([0,10],[0,1])
        c['fuse'].props['posy'].set([0.2,0.5,0.5],[0,0.5,1],w=[1.0,0,0],d=[1,0,0])
        c['fuse'].props['ry'].set([0.01,0.2,0.2,0.20,0.05],[0,0.4,0.6,0.75,1.0],w=[0.9985,0,0,1,0],d=[1,0,0,0,0])
        c['fuse'].props['rz'].set([0.01,0.3,0.3,0.20,0.05],[0,0.4,0.6,0.75,1.0],w=[0.9985,0,0,1,0],d=[1,0,0,0,0])

        c['wing'].offset[:] = [3.35, 0.55, 0.25]
        c['wing'].setAirfoil("naca0003.dat")
        c['wing'].props['posx'].set([0,2.5,5],[0,0.4,1],w=[0,0,0])
        c['wing'].props['posy'].set([0,0],[0,1],w=[0,0])
        c['wing'].props['posz'].set([0,0.7,1.8],[0,0.5,1],w=[0,0,0])
        c['wing'].props['prpx'].set([1,1],[0,1])
        c['wing'].props['prpy'].set([0,0],[0,1])
        c['wing'].props['chord'].set([4.75,2.2,0.1],[0,0.4,1])

        c['tail'].offset[:] = [8.3, 0.55, 0.20]
        c['tail'].setAirfoil("naca0003.dat")
        c['tail'].props['posx'].set([0,1.9],[0,1],w=[0,0])
        c['tail'].props['posy'].set([0,1],[0,1],w=[0,0])
        c['tail'].props['posz'].set([0,1],[0,1])
        c['tail'].props['roty'].set([5,0],[0,1])
        c['tail'].props['prpx'].set([0,0],[0,1])
        c['tail'].props['prpy'].set([0,0],[0,1])
        c['tail'].props['chord'].set([1.3,0.08],[0,1])

        e = numpy.zeros((11,4))
        e[0,:] = [0.20, 0.08, 1, 1]
        e[1,:] = [0.30, 0.10, 0, 0]
        e[2,:] = [0.30, 0.29, 0.5, 0]
        e[3,:] = [0.00, 0.29, 1, 1]
        e[4,:] = [0.00, 0.31, 1, 1]
        e[5,:] = [0.50, 0.35, 1, 0]
        e[6,:] = [1.00, 0.30, 0, 0]
        e[7,:] = [1.00, 0.27, 0, 0]
        e[8,:] = [0.65, 0.27, 0, 0]
        e[9,:] = [0.65, 0.20, 1, 0]
        e[10,:] = [1.2, 0.03, 0, 0]
        e[:,:2] *= 0.8
        l = numpy.linspace(0,1,e.shape[0])
        c['nacelle'].offset[:] = [7.6, 0.4, 0.7]
        c['nacelle'].props['posx'].set(1.3*e[:,0],l)
        c['nacelle'].props['posy'].set([0,0],[0,1])
        c['nacelle'].props['ry'].set(0.4*e[:,1],l,0.4*e[:,2],0.4*e[:,3])
        c['nacelle'].props['rz'].set(0.4*e[:,1],l,0.4*e[:,2],0.4*e[:,3])

        #c['pylon'].offset[:] = [3.8, 0.2, 1.5]
        #c['pylon'].setAirfoil("naca0010")
        #c['pylon'].props['posx'].set([0,0.2],[0,1])
        #c['pylon'].props['posy'].set([0,0.08],[0,1])
        #c['pylon'].props['posz'].set([0,0],[0,1])
        #c['pylon'].props['prpx'].set([1,1],[0,1])
        #c['pylon'].props['prpy'].set([0,0],[0,1])
        #c['pylon'].props['chord'].set([0.75,0.75],[0,1])

        #c['fin'].offset[:] = [8.25, 0.85, 0]
        #c['fin'].props['posx'].set([0,1.5],[0,1])
        #c['fin'].props['posy'].set([0,1.4],[0,1])
        #c['fin'].props['posz'].set([0,0],[0,1])
        #c['fin'].props['rotz'].set([10,0],[0,1])
        #c['fin'].props['prpx'].set([1,1],[0,1])
        #c['fin'].props['prpy'].set([0,0],[0,1])
        #c['fin'].props['chord'].set([1,0.2],[0,1])

        self.computePoints()
  
        c['fuse'].addMembers('Longerons', 2, 1, 12, 15, A1=[0,0,0.95], C1=[1,0,1], A2=[0,1,0.95], C2=[1,1,1])
        c['fuse'].addMembers('Frames', 2, 2, 16, 11, A1=[0,0,0.85], C1=[0,1,1], A2=[1,0,0.85], C2=[1,1,1])

        c['wing'].addMembers('RibsLE', 1, 1, 13, 1, A1=[0,0,0], C1=[0,0.125,1], A2=[1,0,0], C2=[1,0.125,1])
        c['wing'].addMembers('Ribs', 1, 2, 13, 5, A1=[0,0.125,0], C1=[0,0.75,1], A2=[1,0.125,0], C2=[1,0.75,1])
        c['wing'].addMembers('RibsTE', 1, 1, 13, 1, A1=[0,0.75,0], C1=[0,1,1], A2=[1,0.75,0], C2=[1,1,1])
        c['wing'].addMembers('Spars', 1, 2, 2, 12, A1=[0,0.125,0], C1=[1,0.125,1], A2=[0,0.75,0], C2=[1,0.75,1])
        c['wing'].addMembers('Ustiff', 1, 1, 4, 12, A1=[0,0.25,0.9], C1=[1,0.25,1], A2=[0,0.625,0.9], C2=[1,0.625,1])
        c['wing'].addMembers('UstiffL', 1, 1, 4, 12, A1=[0,0.25,0.9], B1=[0,0.255,0.9], C1=[1,0.255,0.9], D1=[1,0.25,0.9], A2=[0,0.625,0.9], B2=[0,0.63,0.9], C2=[1,0.63,0.9], D2=[1,0.625,0.9])
        c['wing'].addMembers('Lstiff', 1, 1, 4, 12, A1=[0,0.25,0], C1=[1,0.25,0.1], A2=[0,0.625,0], C2=[1,0.625,0.1])
        c['wing'].addMembers('LstiffL', 1, 1, 4, 12, A1=[0,0.25,0.1], B1=[0,0.255,0.1], C1=[1,0.255,0.1], D1=[1,0.25,0.1], A2=[0,0.625,0.1], B2=[0,0.63,0.1], C2=[1,0.63,0.1], D2=[1,0.625,0.1])

        c['tail'].addMembers('RibsLE', 1, 1, 7, 1, A1=[0,0,0], C1=[0,0.125,1], A2=[1,0,0], C2=[1,0.125,1])
        c['tail'].addMembers('Ribs', 1, 2, 7, 5, A1=[0,0.125,0], C1=[0,0.75,1], A2=[1,0.125,0], C2=[1,0.75,1])
        c['tail'].addMembers('RibsTE', 1, 1, 7, 1, A1=[0,0.75,0], C1=[0,1,1], A2=[1,0.75,0], C2=[1,1,1])
        c['tail'].addMembers('Spars', 1, 2, 2, 6, A1=[0,0.125,0], C1=[1,0.125,1], A2=[0,0.75,0], C2=[1,0.75,1])
        c['tail'].addMembers('Ustiff', 1, 1, 4, 6, A1=[0,0.25,0.9], C1=[1,0.25,1], A2=[0,0.625,0.9], C2=[1,0.625,1])
        c['tail'].addMembers('UstiffL', 1, 1, 4, 6, A1=[0,0.25,0.9], B1=[0,0.255,0.9], C1=[1,0.255,0.9], D1=[1,0.25,0.9], A2=[0,0.625,0.9], B2=[0,0.63,0.9], C2=[1,0.63,0.9], D2=[1,0.625,0.9])
        c['tail'].addMembers('Lstiff', 1, 1, 4, 6, A1=[0,0.25,0], C1=[1,0.25,0.1], A2=[0,0.625,0], C2=[1,0.625,0.1])
        c['tail'].addMembers('LstiffL', 1, 1, 4, 6, A1=[0,0.25,0.1], B1=[0,0.255,0.1], C1=[1,0.255,0.1], D1=[1,0.25,0.1], A2=[0,0.625,0.1], B2=[0,0.63,0.1], C2=[1,0.63,0.1], D2=[1,0.625,0.1])

        t0 = time.time()
        self.computePoints()        
        print time.time() - t0


if __name__ == '__main__':

    name = 'Supersonic'
    aircraft = Supersonic()
    aircraft.buildStructure()
    aircraft.writeStructure(name)
    aircraft.export.write2Tec(name)
    aircraft.export.write2TecC(name+'_C')
    #aircraft.export.write2IGES(name)
    #aircraft.export.write2EGADS(name+'_EGADS')
    aircraft.plot()
