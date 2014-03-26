from __future__ import division
from tecplot import Tecplot
import numpy, time
import sys

from GeoMACH.PGM.components import Wing, Body, Shell, Junction, Cone
from GeoMACH.PGM.configurations import Configuration
from GeoMACH.PSM import Airframe


class Conventional(Configuration):

    def __init__(self):
        super(Conventional,self).__init__()
        c = self.comps

        self.addComp('HWB', Wing(nx=5, nz=8))
        self.addComp('Prop', Wing(nz=2, left=0, right=0))

        self.separateComps()

        self.addComp('JL', Junction(c['HWB'], 0, 1, [1,1], c['Prop'], mSide=1))
        self.addComp('JR', Junction(c['HWB'], 0, 3, [1,1], c['Prop'], mSide=0))

        self.assembleComponents()

        c['HWB'].setm(0,0,[4,4,7,4,50])
        c['HWB'].setm(0,1,[30,4,4,10,10,4,4,30])
        c['Prop'].setm(0,1,[15,15])

        self.update()

#        c['HWB'].params['pos'].setP([[0,0,-50],[0,5,50]])
        c['HWB'].params['pos'].setP([[0,0,0],[0,0,0]])
        c['HWB'].params['ogn'].setP([0,0,0])
        c['HWB'].params['scl'].setP([3.5])
        c['HWB'].params['nor'].setP([0])
        #c['HWB'].addParam('pos1','pos',[7,3],P=[[50,0,0],[20,0,0],[5,0,0],[0,0,0],[5,0,0],[20,0,0],[50,0,0]],T=[0,0.35,0.46,0.5,0.54,0.65,1],B=[[True,False,False] for j in range(7)], D=[[-3/3.5,0,0],[-3/3.5,0,0],[-5/2.5,0,0],[0,0,0],[5/2.5,0,0],[3/3.5,0,0],[3/3.5,0,0]])
        #c['HWB'].addParam('pos1','pos',[7,3],P=[[50,0,-50],[20,0,-15],[5,0,-4],[0,0,0],[5,0,4],[20,0,15],[50,0,50]],T=[0,0.35,0.46,0.5,0.54,0.65,1],B=[[True,False,True] for j in range(7)], D=[[-3,0,3.5],[-3,0,3.5],[-5,0,2.5],[0,0,0],[5,0,2.5],[3,0,3.5],[3,0,3.5]])
        c['HWB'].addParam('pos1','pos',[7,3],P=[[50,0,-50],[20,0,-15],[-5,0,-4],[-10,0,0],[-5,0,4],[20,0,15],[50,0,50]],T=[0,0.35,0.46,0.5,0.54,0.65,1],B=[[True,False,True] for j in range(7)], D=[[-3,0,3.5],[-3,0,3.5],[-5,0,2.5],[0,0,0],[5,0,2.5],[3,0,3.5],[3,0,3.5]])
        c['HWB'].addParam('wingletL','pos',[2,3],P=[[0,0,0],[3,10,0]],B=[[True for i in range(3)] for j in range(2)],D=[[0,0,5],[0,5,0]],T=[0.97,1])
        c['HWB'].addParam('wingletR','pos',[2,3],P=[[3,10,0],[0,0,0]],B=[[True for i in range(3)] for j in range(2)],D=[[0,0,5],[0,5,0]],T=[0,0.03])
        c['HWB'].addParam('wingletL2','scl',[2,1],P=[0,-1.5],T=[0.97,1])
        c['HWB'].addParam('wingletR2','scl',[2,1],P=[-1.5,0],T=[0,0.03])
        c['HWB'].addParam('scl1','scl',[7,1],P=[0,30,55,60,55,30,0],T=[0,0.35,0.46,0.5,0.54,0.65,1],B=[True for i in range(7)], D=[3,3,5,0,-5,-3,-3])
        c['HWB'].addParam('scl2','scl',[5,1],P=[0,-17,-5,-17,0],T=[0,0.26,0.5,0.74,1],B=[True for i in range(5)],D=[-0.65,0,0,0,0.65])
        #c['HWB'].addParam('scl2','scl',[7,1],P=[0,3.5,4.6,5,4.6,3.5,0],T=[0,0.35,0.46,0.5,0.54,0.65,1])

        c['Prop'].params['pos'].setP([[0,0,0],[0,0,0]])
        c['Prop'].params['ogn'].setP([0,0,0])
        c['Prop'].params['scl'].setP([0])
        c['Prop'].addParam('scl1','scl',[1,3],[5,2.5,0])
        c['Prop'].addParam('pos1','pos',[4,3],[[33,1.1,-14],[33,2.3,-14],[33,2.3,14],[33,1.1,14]],B=[[True,True,True] for i in range(4)],D=[[0,0.3,0],[0,0,0.3],[0,0,0.3],[0,-0.3,0]])
        c['Prop'].addParam('rot1','rot',[4,3],[[0,7,0],[0,0,0],[0,0,0],[0,-7,0]])


        #c['rw'].params['pos'].setP([[0,0,0],[0,3,30]])
        #c['rw'].params['ogn'].setP([0,0,0])
        #c['rw'].addParam('offset','pos',[1,3],P=[18,-1,3])
        #c['rw'].addParam('pos1','pos',[3,3],P=[[0,0,0],[18,0,25],[22,0,29]],T=[0,0.9,1.0])
        #c['rw'].params['scl'].setP([1])
        #c['rw'].addParam('scl1','scl',[3,1],P=[10,5,0.8],T=[0,0.35,1.0])

        self.computePoints()

    def meshStructure(self):
        afm = Airframe(self, 1)

        idims = numpy.linspace(0.25,0.65,7)
        jdims = numpy.linspace(0.02,0.36,16)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('Mlw_1:'+str(i)+':'+str(j),'HWB',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('Mrw_1:'+str(i)+':'+str(j),'HWB',[idims[i],1-jdims[j]],[idims[i+1],1-jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                if i is 0 or i is idims.shape[0]-1:
                    afm.addVertFlip('Mlw_2:'+str(i)+':'+str(j),'HWB',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                    afm.addVertFlip('Mrw_2:'+str(i)+':'+str(j),'HWB',[idims[i],1-jdims[j]],[idims[i],1-jdims[j+1]])
                else:
                    afm.addVertFlip('Mlw_2a:'+str(i)+':'+str(j),'HWB',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1,0.85])
                    afm.addVertFlip('Mlw_2b:'+str(i)+':'+str(j),'HWB',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[0.15,0])
                    afm.addVertFlip('Mrw_2a:'+str(i)+':'+str(j),'HWB',[idims[i],1-jdims[j]],[idims[i],1-jdims[j+1]],w=[1,0.85])
                    afm.addVertFlip('Mrw_2b:'+str(i)+':'+str(j),'HWB',[idims[i],1-jdims[j]],[idims[i],1-jdims[j+1]],w=[0.15,0])

        ind = self.inds
        def run(v1,v2,v3):
            a = 0.2
            b = 0.75
            jdims = numpy.linspace(a,b,7)
            for j in range(jdims.shape[0]-1):
                afm.addMember(name,[
                        [[ind['HWB'], 0,0], [1,jdims[j],v1], [1,jdims[j+1],v1], [0.8,jdims[j],v2], [0.8,jdims[j+1],v2]],
                        [[ind['HWB'], 1,0], [0,1-jdims[j],v1], [0,1-jdims[j+1],v1], [0.2,1-jdims[j],v2], [0.2,1-jdims[j+1],v2]],
                        ])
                afm.addMember(name,[
                        [[ind['HWB'], 0,0], [1,jdims[j],v3+0], [1,jdims[j+1],v3+0], [0.8,jdims[j],v2], [0.8,jdims[j+1],v2]],
                        [[ind['HWB'], 1,0], [0,1-jdims[j],v3+0], [0,1-jdims[j+1],v3+0], [0.2,1-jdims[j],v2], [0.2,1-jdims[j+1],v2]],
                        ])
                afm.addMember(name,[
                        [[ind['HWB'], 0,0], [0,jdims[j],v1], [0,jdims[j+1],v1], [0.2,jdims[j],v2], [0.2,jdims[j+1],v2]],
                        [[ind['HWB'], 1,0], [1,1-jdims[j],v1], [1,1-jdims[j+1],v1], [0.8,1-jdims[j],v2], [0.8,1-jdims[j+1],v2]],
                        ])
                afm.addMember(name,[
                        [[ind['HWB'], 0,0], [0,jdims[j],v3+0], [0,jdims[j+1],v3+0], [0.2,jdims[j],v2], [0.2,jdims[j+1],v2]],
                        [[ind['HWB'], 1,0], [1,1-jdims[j],v3+0], [1,1-jdims[j+1],v3+0], [0.8,1-jdims[j],v2], [0.8,1-jdims[j+1],v2]],
                        ])
                afm.addMember(name,[
                        [[ind['HWB'], 0,0], [0.8,jdims[j],v2], [0.8,jdims[j+1],v2], [0.2,jdims[j],v2], [0.2,jdims[j+1],v2]],
                        [[ind['HWB'], 1,0], [0.2,1-jdims[j],v2], [0.2,1-jdims[j+1],v2], [0.8,1-jdims[j],v2], [0.8,1-jdims[j+1],v2]],
                        ])

        #run(0.36,0.37,0.387)
        #run(0.43,0.45,0.477)
        #run(0.53,0.55,0.577)
        #run(0.62,0.63,0.647)
        run(0.36,0.37,0.38)
        run(0.43,0.45,0.47)
        run(0.53,0.55,0.57)
        run(0.62,0.63,0.64)

        afm.preview('HWB_pvw.dat')
        afm.mesh()
        afm.computeMesh('HWB_str.dat')

if __name__ == '__main__':

#    import cProfile
#    cProfile.run('Conventional()')

    name = 'HWB'
    aircraft = Conventional()

    der = aircraft.getDerivatives('HWB', 'scl', (0,0), FD=False)
    aircraft.oml0.addVars(['der'])
    #aircraft.oml0.addVars(['dx','dy','dz'])
    #aircraft.oml0.P0[:,6] = aircraft.oml0.exportPjtn(der)
    aircraft.oml0.Q[:,6] = numpy.sum(der*der,1)**0.5
    aircraft.oml0.computePoints()
    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name)
    aircraft.meshStructure()

    #aircraft.oml0.plot()
