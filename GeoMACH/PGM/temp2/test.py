from tecplot import Tecplot

t = Tecplot()
t.loadFile('macro1.mcr')
t.writeImage('test.png')
t.runTecplot()
