from __future__ import division
import numpy, time
import sys
from PAM.components import Wing, Body, FullInterface, HalfInterface
from PAM.configurations import Configuration

class HWB(Configuration):
    
    def __init__(self):
        super(HWB,self).__init__() 

        self.addComp('wing',Wing([10,10,10,10,150],[50,10,30,10,10]))
        #self.addComp('pylon1',Wing([30],[30],half=True,opentip=True))
        self.addComp('pylon2',Wing([30],[30],opentip=True))
        #self.addComp('nac1',Body([50,10,30,10,50],[20],[10,10],full=False))
        self.addComp('nac2',Body([50,10,30,10,50],[20],[10,10,10,10],full=True))

        self.separateComps()

        #self.addComp('wing_pylon1',HalfInterface(self.comps,'pylon1',0,'wing',0,[3,0],[1,0]))
        self.addComp('wing_pylon2',FullInterface(self.comps,'pylon2',0,'wing',0,[3,2],[1,3]))
        #self.addComp('pylon_nac1' ,HalfInterface(self.comps,'pylon1',1,'nac1',3,[1,3],[1,1]))
        self.addComp('pylon_nac2' ,FullInterface(self.comps,'pylon2',1,'nac2',3,[2,3],[1,1]))

        self.assembleComponents()

        #return

        c = self.comps

        c['wing'].setAirfoil("rae2822.dat")
        c['wing'].props['posx'].set([0,0.4,1.2,2.31,3],[0,0.08,0.3,0.7,1.0],w=[1,0.4,0,0,0],d=[0,1,0,0,0])
        c['wing'].props['posy'].set([0,0,0.7],[0,0.9,1])
        c['wing'].props['posz'].set([0,3.3,3.4],[0,0.9,1],w=[0,0,0])
        c['wing'].props['prpx'].set([1,1],[0,1])
        c['wing'].props['prpy'].set([0,0],[0,1])
        c['wing'].props['chord'].set([2.95,2.5,1.3,0.85,0.55,0.2],[0,0.08,0.3,0.5,0.7,1.0],w=[1,0.5,0.5,0,0,0],d=[0,1,1,0,0,0])

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
        e[:,:2] *= 0.5
        l = numpy.linspace(0,1,e.shape[0])
        #c['nac1'].offset[:] = [1.3, 0.4, 0.0]
        #c['nac1'].props['posx'].set(e[:,0],l)
        #c['nac1'].props['posy'].set([0,0],[0,1])
        #c['nac1'].props['ry'].set(e[:,1],l,e[:,2],e[:,3])
        #c['nac1'].props['rz'].set(e[:,1],l,e[:,2],e[:,3])

        c['nac2'].offset[:] = [2.15, 0.405, 0.4]
        c['nac2'].props['posx'].set(e[:,0],l)
        c['nac2'].props['posy'].set([0,0],[0,1])
        c['nac2'].props['ry'].set(e[:,1],l,e[:,2],e[:,3])
        c['nac2'].props['rz'].set(e[:,1],l,e[:,2],e[:,3])

        #c['pylon1'].offset[:] = [1.3, 0.1, 0.0]
        #c['pylon1'].setAirfoil("rae2822.dat")
        #c['pylon1'].props['posx'].set([0,0],[0,1.0])
        #c['pylon1'].props['posy'].set([0,0.15],[0,1])
        #c['pylon1'].props['posz'].set([0,0],[0,1])
        #c['pylon1'].props['prpx'].set([1,1],[0,1])
#        c['pylon1'].props['prpy'].set([0,0],[0,1])
        #c['pylon1'].props['chord'].set([0.3,0.3],[0,1.0])

        c['pylon2'].offset[:] = [2.15, 0.105, 0.4]
        c['pylon2'].setAirfoil("rae2822.dat")
        c['pylon2'].props['posx'].set([0,0],[0,1.0])
        c['pylon2'].props['posy'].set([0,0.15],[0,1])
        c['pylon2'].props['posz'].set([0,0],[0,1])
        c['pylon2'].props['prpx'].set([1,1],[0,1])
        c['pylon2'].props['prpy'].set([0,0],[0,1])
        c['pylon2'].props['rotz'].set([10,0],[0,1])
        c['pylon2'].props['chord'].set([0.3,0.3],[0,1.0])

        self.computePoints()

        c['wing'].addMembers('RibsLE', 1, 1, 13, 1, A1=[0,0,0], C1=[0,0.125,1], A2=[1,0,0], C2=[1,0.125,1])
        c['wing'].addMembers('Ribs', 1, 2, 13, 5, A1=[0,0.125,0], C1=[0,0.75,1], A2=[1,0.125,0], C2=[1,0.75,1])
        c['wing'].addMembers('RibsTE', 1, 1, 13, 1, A1=[0,0.75,0], C1=[0,1,1], A2=[1,0.75,0], C2=[1,1,1])
        c['wing'].addMembers('Spars', 1, 2, 2, 12, A1=[0,0.125,0], C1=[1,0.125,1], A2=[0,0.75,0], C2=[1,0.75,1])
        c['wing'].addMembers('Ustiff', 1, 1, 4, 12, A1=[0,0.25,0.9], C1=[1,0.25,1], A2=[0,0.625,0.9], C2=[1,0.625,1])
        c['wing'].addMembers('UstiffL', 1, 1, 4, 12, A1=[0,0.25,0.9], B1=[0,0.255,0.9], C1=[1,0.255,0.9], D1=[1,0.25,0.9], A2=[0,0.625,0.9], B2=[0,0.63,0.9], C2=[1,0.63,0.9], D2=[1,0.625,0.9])
        c['wing'].addMembers('Lstiff', 1, 1, 4, 12, A1=[0,0.25,0], C1=[1,0.25,0.1], A2=[0,0.625,0], C2=[1,0.625,0.1])
        c['wing'].addMembers('LstiffL', 1, 1, 4, 12, A1=[0,0.25,0.1], B1=[0,0.255,0.1], C1=[1,0.255,0.1], D1=[1,0.25,0.1], A2=[0,0.625,0.1], B2=[0,0.63,0.1], C2=[1,0.63,0.1], D2=[1,0.625,0.1])

if __name__ == '__main__':

    name = 'HWB'
    aircraft = HWB()
    aircraft.export.write2Tec(name)
    aircraft.export.write2TecC(name+'_C')
#    aircraft.oml0.write2Tec(name)
    aircraft.plot()
