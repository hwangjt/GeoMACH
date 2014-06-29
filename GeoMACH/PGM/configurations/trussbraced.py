from __future__ import division
from tecplot import Tecplot
import numpy, time
import sys

from GeoMACH.PGM.components import Wing, Body, Shell, Junction, Cone, Tip
from GeoMACH.PGM.configurations import Configuration
from GeoMACH.PSM import Airframe


class Trussbraced(Configuration):

    def define_primitive_comps(self):
        comps = {
            'fu': Body(nx=17, ny=6, nz=4),
            'lw': Wing(nx=7, nz=7, right=0),
            'rw': Wing(nx=7, nz=7, left=0),
            'ls': Wing(nx=4, nz=4, left=0, right=0),
            'lv': Wing(left=0, right=0),
            'vt': Wing(nx=4,nz=4,right=0),
            'lt': Wing(right=0),
        }
        return comps

    def define_interpolant_comps(self):
        c = self.comps
        comps = {
            'fu_n': Cone(c['fu'], 0),
            'fu_t': Cone(c['fu'], 1),
            'lw_T': Tip(c['lw'], 1),
            'rw_T': Tip(c['rw'], 0),
            'lt_T': Tip(c['lt'], 1),
            'vt_T': Tip(c['vt'], 1),
            'lw_fu': Junction(c['fu'], 'lft', 0, [0,1], c['lw'], mSide=0),
            'rw_fu': Junction(c['fu'], 'rgt', 2, [0,7], c['rw'], mSide=1),
            'ls_fu': Junction(c['fu'], 'lft', 0, [4,2], c['ls'], mSide=0),
            'ls_lw': Junction(c['lw'], 'low', 3, [4,1], c['ls'], mSide=1),
            'lv_lw': Junction(c['lw'], 'low', 3, [1,3], c['lv'], mSide=1),
            'lv_ls': Junction(c['ls'], 'upp', 3, [1,0], c['lv'], mSide=0),
            'vt_fu': Junction(c['fu'], 'top', 0, [1,11], c['vt'], mSide=0),
            'lt_vt': Junction(c['vt'], 'low', 1, [0,0], c['lt'], mSide=0),
            #'lp_ln': Junction(c['ln'], 'tp0', 2, [1,0], c['lp'], mSide=1),
            #'rp_rw': Junction(c['rw'], 'low', 1, [1,0], c['rp'], mSide=0),
            #'rp_rn': Junction(c['rn'], 'tp0', 2, [1,0], c['rp'], mSide=1),
            #'lw_fu': Junction(c['fu'], 'lft', 0, [2,1], c['lw'], mSide=0),
            #'rw_fu': Junction(c['fu'], 'rgt', 2, [2,5], c['rw'], mSide=1),
            #'lt_fu': Junction(c['fu'], 'lft', 0, [1,9], c['lt'], mSide=0),
            #'rt_fu': Junction(c['fu'], 'rgt', 2, [1,0], c['rt'], mSide=1),
                 }
        return comps

    def define_dvs(self):
        self.add_dv('jtn', [1], val=0.0, lower=-2, upper=3)
        #self.add_dv('fuL', [25], val=0.0, lower=0, upper=1)
        self.add_dv('lwU', [100], val=0.0, lower=0, upper=1)
        self.add_dv('lwL', [100], val=0.0, lower=0, upper=1)
        self.add_dv('lsU', [64], val=0.0, lower=0, upper=1)
        self.add_dv('lsL', [64], val=0.0, lower=0, upper=1)
        self.add_dv('lvU', [25], val=0.0, lower=0, upper=1)
        self.add_dv('lvL', [25], val=0.0, lower=0, upper=1)
        self.add_dv('ltU', [64], val=0.0, lower=0, upper=1)
        self.add_dv('ltL', [64], val=0.0, lower=0, upper=1)

    def apply_dvs(self):
        Das, Dis, Djs = [], [], []

        def addShape(Das, Dis, Djs, 
                     dv_name, comp_name, prop_name, param_name, 
                     sgn, ni, nj):
            Da = numpy.ones(ni*nj) * sgn
            Di = numpy.zeros(ni*nj, int)
            Dj = numpy.zeros(ni*nj, int)

            param = self.comps[comp_name].props[prop_name].params[param_name]
            dv = self.dvs[dv_name]

            ind = 0
            for j in xrange(nj):
                for i in xrange(ni):
                    param.param_vec[i,j] = sgn*dv.vec[ind]
                    Di[ind] = param.param_ind[i,j]
                    Dj[ind] = dv.ind[ind]
                    ind += 1
                    
            Das.append(Da)
            Dis.append(Di)
            Djs.append(Dj)

        def addNonZero(a, i, j):
            Das.append(a)
            Dis.append(i)
            Djs.append(j)

        #addShape(Das, Dis, Djs, 'fuL', 'fu', ('shX', 'lft'), 's', -1, 5, 5)
        addShape(Das, Dis, Djs, 'lwU', 'lw', ('shY', 'upp'), 's',  1, 10, 10)
        addShape(Das, Dis, Djs, 'lwL', 'lw', ('shY', 'low'), 's', -1, 10, 10)
        addShape(Das, Dis, Djs, 'lsU', 'ls', ('shY', 'upp'), 's',  1,  8,  8)
        addShape(Das, Dis, Djs, 'lsL', 'ls', ('shY', 'low'), 's', -1,  8,  8)
        addShape(Das, Dis, Djs, 'lvU', 'lv', ('shY', 'upp'), 's',  1,  5,  5)
        addShape(Das, Dis, Djs, 'lvL', 'lv', ('shY', 'low'), 's', -1,  5,  5)
        addShape(Das, Dis, Djs, 'ltU', 'lt', ('shY', 'upp'), 's',  1,  8,  8)
        addShape(Das, Dis, Djs, 'ltL', 'lt', ('shY', 'low'), 's', -1,  8,  8)

        par_lw = self.comps['lw'].props['pos'].params['pos1']
        par_ls = self.comps['ls'].props['pos'].params['pos1']
        par_lv = self.comps['lv'].props['pos'].params['offset']
        dv = self.dvs['jtn']

        par_lw.param_vec[1,2,0] = 15 + dv.vec[0]
        par_ls.param_vec[1,2,0] = (15 + dv.vec[0]) * 0.9
        par_lv.param_vec[0,2,0] = 2 + (15 + dv.vec[0]) * 6.8/15

        addNonZero(1, par_lw.param_ind[1,2,0], dv.ind[0])
        addNonZero(0.9, par_ls.param_ind[1,2,0], dv.ind[0])
        addNonZero(6.8/15, par_lv.param_ind[0,2,0], dv.ind[0])

        return Das, Dis, Djs

    def set_oml_resolution(self):
        comps = self.comps
        comps['lv'].faces['upp'].setm(0,[4])
        comps['lv'].faces['upp'].setm(1,[4])
        comps['lv'].faces['low'].setm(0,[4])
        comps['lv'].faces['low'].setm(1,[4])
        comps['lt'].faces['upp'].setm(0,[4])
        comps['lt'].faces['upp'].setm(1,[4])
        comps['lt'].faces['low'].setm(0,[4])
        comps['lt'].faces['low'].setm(1,[4])

        comps['fu'].faces['rgt'].setm(1,[42,4,4,4,7,4,4,4,4,4,33,4,4,12,4,4,4])
        comps['fu'].faces['lft'].setm(0,[4,4,8,8,4,4])
        comps['lw'].faces['upp'].setm(1,[8,4,4,5,4,4,22])
        comps['rw'].faces['upp'].setm(1,[8,4,4,5,4,4,22][::-1])
#        comps['rw'].faces['low'].setm(0,[4,4,4,7,4,4,4])
        comps['ls'].faces['upp'].setm(1,[10,4,4,10])

        comps['lw_T'].faces['def'].setm(0,[10])
        comps['lt_T'].faces['def'].setm(0,[10])
        comps['vt_T'].faces['def'].setm(0,[10])

        #comps['fu'].faces['rgt'].setm(1,[18,4,4,4,4,8,4,15,4,4,10,4])
        #comps['fu'].faces['rgt'].setm(0,[4,4,4,4])
        #comps['fu'].faces['top'].setm(0,[8,8])
        #comps['lw'].faces['upp'].setm(1,[6,4,4,20])

    def define_oml_parameters(self):
        c = self.comps

        c['lw'].setAirfoil('rae2822.dat')
        c['rw'].setAirfoil('rae2822.dat')

        c['fu'].props['nor'].addParam('nor1',[1,3],P=[1.0,0.0,1.0])
        c['fu'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[38,0,0]])
        c['fu'].props['pos'].addParam('nose',[2,3],P=[[0,-1.1,0],[0,0,0]],T=[0,0.13],B=[False,True])
        c['fu'].props['pos'].addParam('tail',[2,3],P=[[0,0,0],[0,1.6,0]],T=[0.75,1.0],B=[False,False])
        c['fu'].props['scl'].addParam('rad',[4,1],P=[0.65,2,2,0.4],T=[0,0.14,0.75,1.0],B=[False,True,False,False])
        c['fu'].props['scl'].addParam('tail',[2,3],P=[[0,0,0],[-0.3,0,0]],T=[0.75,1.0])
        c['fu'].props['flt'].addParam('flt1',[2,4],P=[[0,0,0.9,0.9],[0,0,0.9,0.9]],T=[0.33,0.50])
        c['fu'].props['flt'].addParam('flt2',[2,4],P=[[0.9,0.9,0,0],[0.9,0.9,0,0]],T=[0.32,0.55])

        c['lw'].props['pos'].addParam('offset',[1,3],P=[14,1.5,2.1])
        c['lw'].props['pos'].addParam('pos1',[3,3],P=[[0,0,0],[1.9,0.75,15],[3.8,1.5,30]])
        c['lw'].props['scl'].addParam('scl1',[2,1],P=[4.9,0.9])
        c['lw'].props['shY','upp'].addParam2('s',[10,10],Tu=[0.05,0.95])
        c['lw'].props['shY','low'].addParam2('s',[10,10],Tu=[0.05,0.95])

        c['rw'].props['pos'].addParam('offset',[1,3],P=[14,1.5,-2.1])
        c['rw'].props['pos'].addParam('pos1',[2,3],P=[[3.8,1.5,-30],[0,0,0]])
        c['rw'].props['scl'].addParam('scl1',[2,1],P=[0.9,4.9])

        c['ls'].props['pos'].addParam('offset',[1,3],P=[15,-1.6,2.3])
        c['ls'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[1,3.0,15*0.9]])
        c['ls'].props['pos'].addParam('curv',[2,3],P=[[0,0,0],[0,0.6,0]],B=[True,False],T=[0.8,1.0])
        c['ls'].props['scl'].addParam('scl1',[2,1],P=[2.5,1.9])
        c['ls'].props['nor'].addParam('nor',[1,1],P=[1.0])
        c['ls'].props['rot'].addParam('rot',[2,3],P=[[0,0,0],[-60,0,0]],T=[0.7,1])
        c['ls'].props['shY','upp'].addParam2('s',[8,8],Tu=[0.05,0.95])
        c['ls'].props['shY','low'].addParam2('s',[8,8],Tu=[0.05,0.95])

        c['lv'].props['pos'].addParam('offset',[1,3],P=[15.5,0,8.8])
        c['lv'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[0,1.58,0]])
        c['lv'].props['scl'].addParam('scl1',[2,1],P=[1,1.4])
        c['lv'].props['nor'].addParam('nor',[1,1],P=[1.0])
        c['lv'].props['shY','upp'].addParam2('s',[5,5],Tu=[0.05,0.95])
        c['lv'].props['shY','low'].addParam2('s',[5,5],Tu=[0.05,0.95])

        c['vt'].props['pos'].addParam('offset',[1,3],P=[31.3,2.1,0])
        c['vt'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[3,6,0]])
        c['vt'].props['scl'].addParam('scl1',[2,1],P=[5.8,3.3])
        c['vt'].props['nor'].addParam('nor',[1,3],P=[1.0,0.0,0.0])

        c['lt'].props['pos'].addParam('offset',[1,3],P=[33.9,6.7,0.25])
        c['lt'].props['pos'].addParam('pos1',[2,3],P=[[0,0,0],[3.3,0,5]])
        c['lt'].props['scl'].addParam('scl1',[2,1],P=[2.9,1])
        c['lt'].props['shY','upp'].addParam2('s',[8,8],Tu=[0.05,0.95])
        c['lt'].props['shY','low'].addParam2('s',[8,8],Tu=[0.05,0.95])

        c['fu_n'].props['fC1'].params['fC1'].setP([4])
        c['fu_t'].props['fC1'].params['fC1'].setP([4])
        c['lw_T'].props['fC1'].params['fC1'].setP([0.5])
        c['vt_T'].props['fC1'].params['fC1'].setP([1.0])
        c['lt_T'].props['fC1'].params['fC1'].setP([0.5])
        c['lw_fu'].props['mC1'].params['mC1'].setP([0.1])
        c['lw_fu'].props['fC1'].params['fC1'].setP([1.0])
        c['ls_fu'].props['mC1'].params['mC1'].setP([0.02])
        c['ls_fu'].props['fC1'].params['fC1'].setP([0.02])
        c['ls_lw'].props['mC1'].params['mC1'].setP([0.1])
        c['ls_lw'].props['fC1'].params['fC1'].setP([0.1])
        c['lv_ls'].props['mC1'].params['mC1'].setP([0.1])
        c['lv_ls'].props['fC1'].params['fC1'].setP([0.1])
        c['lv_lw'].props['mC1'].params['mC1'].setP([0.1])
        c['lv_lw'].props['fC1'].params['fC1'].setP([0.1])
        c['vt_fu'].props['mC1'].params['mC1'].setP([0.1])
        c['vt_fu'].props['fC1'].params['fC1'].setP([0.1])
        c['lt_vt'].props['mC1'].params['mC1'].setP([0.1])
        c['lt_vt'].props['fC1'].params['fC1'].setP([0.1])

    def meshStructure(self):
        afm = Airframe(self, 0.4)

        idims = numpy.linspace(0.3,0.85,7)
        jdims = numpy.linspace(0,1,17)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('Mlw1:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                if i is 0 or i is idims.shape[0]-1:
                    afm.addVertFlip('Mlw2:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                else:
                    afm.addVertFlip('Mlw2a:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1,0.85])
                    afm.addVertFlip('Mlw2b:'+str(i)+':'+str(j),'lw',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[0.15,0])

        idims = numpy.linspace(0.3,0.85,7)
        jdims = numpy.linspace(0,0.9,17)
        for i in range(idims.shape[0]):
            if i is 0 or i is idims.shape[0]-1:
                afm.addCtrVert('Mcw2:'+str(i)+':'+str(j),'lw','rw',idims[i])
            else:
                afm.addCtrVert('Mcw2a:'+str(i)+':'+str(j),'lw','rw',idims[i],w=[1,0.85])
                afm.addCtrVert('Mcw2b:'+str(i)+':'+str(j),'lw','rw',idims[i],w=[0.15,0])
        for i in range(idims.shape[0]-1):
            afm.addCtr('Mcwu:','lw','rw',0,[idims[i],idims[i+1]])
        for i in range(idims.shape[0]-1):
            afm.addCtr('Mcwl:','lw','rw',1,[1-idims[i],1-idims[i+1]])
        afm.addCtrVert('Mcwsec:'+str(i)+':'+str(j),'lw','rw',0.18)

        idims = numpy.linspace(0.25,0.65,2)
        jdims = numpy.linspace(0,1,11)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVertFlip('Mlt1:'+str(i)+':'+str(j),'lt',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('Mvt1:'+str(i)+':'+str(j),'vt',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('Mls1:'+str(i)+':'+str(j),'ls',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
                afm.addVertFlip('Mlv1:'+str(i)+':'+str(j),'lv',[idims[i],jdims[j]],[idims[i+1],jdims[j]])
        for i in range(idims.shape[0]):
            for j in range(jdims.shape[0]-1):
                afm.addVertFlip('Mlt2:'+str(i)+':'+str(j),'lt',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                afm.addVertFlip('Mvt2:'+str(i)+':'+str(j),'vt',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                afm.addVertFlip('Mls2:'+str(i)+':'+str(j),'ls',[idims[i],jdims[j]],[idims[i],jdims[j+1]])
                afm.addVertFlip('Mlv2:'+str(i)+':'+str(j),'lv',[idims[i],jdims[j]],[idims[i],jdims[j+1]])

        idims = numpy.linspace(0,1,5)
        jdims = numpy.linspace(0,1,21)
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]):
                afm.addVert('MfuF1:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[0,2])
                afm.addVert('MfuF2:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[1,3])
                afm.addVert('MfuF3:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[2,0])
                afm.addVert('MfuF4:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i+1],jdims[j]],w=[1.0,0.94],i=[3,1])
        for i in range(idims.shape[0]-1):
            for j in range(jdims.shape[0]-1):
                afm.addVert('MfuL1:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[0,2])
                afm.addVert('MfuL2:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[1,3])
                afm.addVert('MfuL3:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[2,0])
                afm.addVert('MfuL4:'+str(i)+':'+str(j),'fu',[idims[i],jdims[j]],[idims[i],jdims[j+1]],w=[1.0,0.97],i=[3,1])
        for j in range(jdims.shape[0]-1):
            afm.addVertFlip('Mfu0a:'+str(j),'fu',[0.4,jdims[j]],[0.4,jdims[j+1]],w=[1.0,0.5],i=[0,2])
            afm.addVertFlip('Mfu0b:'+str(j),'fu',[0.4,jdims[j]],[0.4,jdims[j+1]],w=[0.5,0.0],i=[0,2])

        afm.preview('trussbraced_pvw.dat')
        afm.mesh()
        afm.computeMesh('trussbraced_str.dat')


if __name__ == '__main__':

    import cProfile

    name = 'trussbraced'
    aircraft = Trussbraced()

    #der = aircraft.get_derivatives('lw', 'scl1', (0,0), useFD=False)
    #aircraft.oml0.addVars(['der'])
    ##aircraft.oml0.addVars(['dx','dy','dz'])
    ##aircraft.oml0.P0[:,6] = aircraft.oml0.exportPjtn(der)
    #aircraft.oml0.Q[:,6] = numpy.sum(der*der,1)**0.5
    aircraft.oml0.computePoints()
    aircraft.oml0.write2Tec(name)
    aircraft.oml0.write2TecC(name)
    aircraft.oml0.write2IGES(name)

    aircraft.comps['lw'].add_thk_con('name', 
                                     numpy.linspace(0.1,0.9,10),
                                     numpy.linspace(0.1,0.9,10),
                                     0.25)

    #aircraft.test_derivatives_dv()
    aircraft.meshStructure()

    #aircraft.test_derivatives()
    #aircraft.oml0.plot()
    #aircraft.meshStructure()

#    cProfile.run('aircraft.meshStructure()')
