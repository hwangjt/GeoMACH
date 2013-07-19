from tecplot import Tecplot
import os

if 1:
    counter = 200
    t = Tecplot()
    for i in range(counter):
        t.importDataSet('temp'+str(i),True)
        t.createMirror(1,165,3)
        t.importDataSet('temp'+str(i)+'_str')
        t.setTranslucency(1,2617,90)
        t.setTranslucency(166,330,0)
#        t.setTransparency(False)
        t.setRelCameraPosition(-20, 10, 20, -140)
        t.writeImage('temp%03d'%(i))
    t.runTecplot()
    for i in range(counter,counter+30):
        os.system('cp temp%03d.jpg temp%03d.jpg'%(counter-1,i))
    t.makeVideo()
