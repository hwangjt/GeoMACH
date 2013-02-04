from __future__ import division
from tecplot import Tecplot
from conventional import Conventional
import numpy, time
import sys
from PAM.components import Wing, Body, Shell, Junction, Cone
from PAM.configurations import Configuration



name = 'conventional'
aircraft = Conventional()
aircraft.oml0.addVars(['der'])

t = Tecplot()
counter = 100   
for c in aircraft.comps.keys():
    comp = aircraft.comps[c]
    for p in comp.params.keys():
        par = comp.params[p]
        for j in range(par.P.shape[1]):
            for i in range(par.P.shape[0]):
                print c, p, i, j
                der = aircraft.getDerivatives(c, p, (i,j), clean=False, FD=False)
                aircraft.oml0.P0[:,6] = aircraft.oml0.exportPjtn(der)
                aircraft.oml0.write2Tec(name+str(counter))
                t.importDataSet(name+str(counter),aircraft.oml0.var,True)
                t.setTransparency(False)
                t.setCamera(66.0214, 133.832, -140.403, -251, 268, 170.3)
                t.plotContours(7)
                t.writeImage('temp%03d'%(counter))
                counter += 1

t.runTecplot()
t.makeVideo()

