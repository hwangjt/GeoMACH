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
        c['wing'].props['posx'].set([0,1.4,2.2],[0,0.5,1.0],w=[1,0,0])
        c['wing'].props['posy'].set([0,0],[0,1])
        c['wing'].props['posz'].set([0,3.4],[0,1],w=[0,0])
        c['wing'].props['prpx'].set([0,0],[0,1])
        c['wing'].props['prpy'].set([0,0],[0,1])
        c['wing'].props['chord'].set([2.0,0.5,0.2],[0,0.5,1.0],w=[1,0,0])

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

        c['nac2'].offset[:] = [1.3, 0.4, 0.4]
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

        c['pylon2'].offset[:] = [1.3, 0.1, 0.4]
        c['pylon2'].setAirfoil("rae2822.dat")
        c['pylon2'].props['posx'].set([0,0],[0,1.0])
        c['pylon2'].props['posy'].set([0,0.15],[0,1])
        c['pylon2'].props['posz'].set([0,0],[0,1])
        c['pylon2'].props['prpx'].set([1,1],[0,1])
#        c['pylon2'].props['prpy'].set([0,0],[0,1])
        c['pylon2'].props['chord'].set([0.3,0.3],[0,1.0])

        self.computePoints()

if __name__ == '__main__':

    name = 'HWB'
    aircraft = HWB()
    aircraft.plot()
#    aircraft.oml0.write2Tec(name)
