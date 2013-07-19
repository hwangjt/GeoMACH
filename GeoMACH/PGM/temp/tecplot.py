from __future__ import division
import numpy, os


class Tecplot(object):

    def __init__(self):
        self.lines = []
        self.lines.append('#!MC 1300')
        self.lines.append('# Created by Tecplot 360 build 13.1.0.15185')
        self.lines.append('$!VarSet |pwd| = \'' + os.getcwd() + '\'')

    def importDataSet(self, filename, new=False):
        if not filename[-4:]=='.dat':
            filename = filename + '.dat'
        self.lines.append('$!READDATASET  \'"|pwd|/' + filename + '" \'')
        if new:
            self.lines.append('  READDATAOPTION = NEW')
            self.lines.append('  RESETSTYLE = YES')
        else:
            self.lines.append('  READDATAOPTION = APPEND')
            self.lines.append('  RESETSTYLE = NO')
        self.lines.append('  INCLUDETEXT = NO')
        self.lines.append('  INCLUDEGEOM = NO')
        self.lines.append('  INCLUDECUSTOMLABELS = NO')
        self.lines.append('  VARLOADMODE = BYNAME')
        self.lines.append('  ASSIGNSTRANDIDS = YES')
        self.lines.append('  INITIALPLOTTYPE = CARTESIAN3D')
        self.lines.append('  VARNAMELIST = \'"x" "y" "z"\'')

    def setTransparency(self, transparent):
        if transparent:
            self.lines.append('$!FIELDLAYERS USETRANSLUCENCY = YES')
        else:
            self.lines.append('$!FIELDLAYERS USETRANSLUCENCY = NO')

    def setTranslucency(self, s1, s2, val):
        self.lines.append('$!FIELDMAP ['+str(s1)+'-'+str(s2)+']  EFFECTS{SURFACETRANSLUCENCY = '+str(val)+'}')

    def createMirror(self, s1, s2, axis):
        self.lines.append('$!CREATEMIRRORZONES')
        self.lines.append('  SOURCEZONES =  ['+str(s1)+'-'+str(s2)+']')
        self.lines.append('  MIRRORVARS =  ['+str(axis)+']')

    def setRelCameraPosition(self, x, y, z, twist):
        def arctan(x,y):
            if x==0:
                if y > 0:
                    t = numpy.pi/2.0
                elif y < 0:
                    t = 3*numpy.pi/2.0
            elif y==0:
                if x > 0:
                    t = 0
                elif x < 0:
                    t = numpy.pi
            elif x<0:
                t = numpy.arctan(y/x) + numpy.pi
            elif y<0:
                t = numpy.arctan(y/x) + 2*numpy.pi
            elif y>0:
                t = numpy.arctan(y/x)
            else:
                t = 0
            if t>numpy.pi:
                t -= 2*numpy.pi
            return t*180.0/numpy.pi

        psi = arctan(z,(x**2+y**2)**0.5)
        theta = arctan(x,y)

        n = numpy.array([x,y,z])
        e2 = numpy.array([1,0,0])
        e3 = numpy.array([0,0,1])
        v = e2 - n*numpy.dot(e2,n)/numpy.dot(n,n)
        alpha = arctan(v[2],(v[0]**2+v[1]**2)**0.5)
        if numpy.dot(n,numpy.cross(e3,v)) > 0:
            alpha *= -1           

        self.lines.append('$!THREEDVIEW')
        self.lines.append('  PSIANGLE = '+str(psi))
        self.lines.append('  THETAANGLE = '+str(theta))
        self.lines.append('  ALPHAANGLE = '+str(twist))
        self.lines.append('  VIEWERPOSITION')
        self.lines.append('    {')
        self.lines.append('    X = '+str(x))
        self.lines.append('    Y = '+str(y))
        self.lines.append('    Z = '+str(z))
        self.lines.append('    }')
        

    def writeMcr(self, filename):
        if not filename[-4:]=='.mcr':
            filename = filename + '.mcr'
        f = open(filename,'w')
        for line in self.lines:
            f.write(line + '\n')
        f.close()

    def writeImage(self, filename, imgwidth=2000, aa=2):
        if not filename[-4:]=='.png':
            filename = filename + '.png'
        self.lines.append('$!VIEW CENTER')
        self.lines.append('$!VIEW DATAFIT')
        self.lines.append('$!REDRAWALL')
        self.lines.append('$!EXPORTSETUP EXPORTFORMAT = PNG')
        self.lines.append('$!EXPORTSETUP IMAGEWIDTH = '+str(imgwidth))
        self.lines.append('$!EXPORTSETUP USESUPERSAMPLEANTIALIASING = YES')
        self.lines.append('$!EXPORTSETUP SUPERSAMPLEFACTOR = '+str(aa))
        self.lines.append('$!EXPORTSETUP EXPORTFNAME = \'|pwd|/'+filename+'\'')
        self.lines.append('$!EXPORT')
        self.lines.append('  EXPORTREGION = CURRENTFRAME')

    def runTecplot(self):
        self.writeMcr('macro.mcr')
        os.system('tec360 -mesa -b macro.mcr')
        #os.system('rm macro.mcr')
        #os.system('rm batch.log')

    def makeVideo(self, fps=10):
#	os.system('jpegoptim *.jpg --max=50')
        os.system('mencoder mf://*.png -mf fps='+str(fps)+':type=png -ovc x264 -x264encopts bitrate=32000 -o output.avi')
        #os.system('mencoder mf://*.jpg -mf w=800:h=600:fps=15:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi')
if __name__ == '__main__':

    t = Tecplot()
    t.importDataSet('conventional')
    t.setTransparency(False)
    t.createMirror(1,165,3)
    t.setRelCameraPosition(-20, 10, 20, -140)
    t.writeImage('test0')
