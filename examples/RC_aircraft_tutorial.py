# RC aircraft tutorial
# Ney Secco
# March, 2015

#IMPORTS
from __future__ import division

from GeoMACH.PGM.core import PGMconfiguration, PGMparameter, PGMdv
from GeoMACH.PGM.components import PGMwing, PGMbody, PGMshell
from GeoMACH.PGM.components import PGMjunction, PGMtip, PGMcone
from GeoMACH.PSM import Airframe
import numpy

##CREATE AIRCRAFT OBJECT
class test(PGMconfiguration):

	#DEFINE COMPONENTS
	def _define_comps(self):
		self.comps['fuse'] = PGMbody(num_x=12, num_y=4, num_z=4) #Creates the fuselage component with 12 surfaces on the longitudinal direction, 4 surfaces on y direction and 4 surfaces on z direction
	        self.comps['cone_f'] = PGMcone(self, 'fuse', 'front',2) #Adds a cone to close the front of the fuselage. The last parameter controls how the cone uses the tangential component of the primary component (fuselage) to generate the junction
	        self.comps['cone_r'] = PGMcone(self, 'fuse', 'rear',0.1) #Adds a cone to close the rear of the fuselage. The last parameter controls how the cone uses the tangential component of the primary component (fuselage) to generate the junction

		self.comps['wing_l'] = PGMwing(num_x=3,num_z=4,left_closed=True) #Creates the left wing component. We only need to close the left tip
		self.comps['wingtip_l'] = PGMtip(self, 'wing_l', 'left', 0.1) #Creates the left wing tip
		self.comps['junc_fus_wing_l'] = PGMjunction(self, 'fuse', 'lft', 'E', [0,3], 'wing_l', 'right') #We will connect the 'right' end of the left wing ('wing_l') to the left face ('lft') of the fuselage ('fuse'). We look at the face on the female component in the orientation that has the upper face of the male component facing up. The 'v' direction of the female face is pointing to the right, so we specify 'E' for east. The northwest-most surface of the junction has coordinates (i=0,j=3), where the i axis points down, and the j axis points right; i.e., (i=0,j=0) is the northwest-most surface on the female face.
		self.comps['wing_r'] = PGMwing(num_x=3,num_z=4,right_closed=True) #Creates the right wing component. We only need to close the right tip
		self.comps['wingtip_r'] = PGMtip(self, 'wing_r', 'right', 0.1) #Creates the right wing tip
		self.comps['junc_fus_wing_r'] = PGMjunction(self, 'fuse', 'rgt', 'W', [0,4], 'wing_r', 'left') #We will connect the 'left' end to the right wing ('wing_r') to the right face ('rgt') of the fuselage ('fuse').

		self.comps['htail_l'] = PGMwing(num_x=1, num_z=3, left_closed=True) #This creates the left horizontal tail component
		self.comps['htailtip_l'] = PGMtip(self, 'htail_l', 'left', 0.1) #Creates the left tail tip
		self.comps['junc_fus_htail_l'] = PGMjunction(self, 'fuse', 'lft', 'E', [1,9], 'htail_l', 'right') #We will connect the 'right' end to the left tail ('htail_l') to the left face ('lft') of the fuselage ('fuse').
		self.comps['htail_r'] = PGMwing(num_x=1, num_z=3, right_closed=True) #This creates the right horizontal tail component
		self.comps['htailtip_r'] = PGMtip(self, 'htail_r', 'right', 0.1) #Creates the right tail tip
		self.comps['junc_fus_htail_r'] = PGMjunction(self, 'fuse', 'rgt', 'W', [1,0], 'htail_r', 'left') #We will connect the 'left' end to the right tail ('htail_r') to the right face ('rgt') of the fuselage ('fuse').
		self.comps['vtail'] = PGMwing(num_x=1, num_z=3, left_closed=True) #This creates the vertical tail component. We will define this component with the upper face pointing to the right.
		self.comps['vtailtip'] = PGMtip(self, 'vtail', 'left', 0.1) #Creates the vertical tail tip
		self.comps['junc_fus_vtail'] = PGMjunction(self, 'fuse', 'top', 'E', [1,9], 'vtail', 'right') #We will connect the 'right' end to the vertical tail ('vtail') to the top face ('top') of the fuselage ('fuse').

	#DEFINE PARAMETERS OF EACH PROPERTY
	def _define_params(self):
		fuse_props = self.comps['fuse'].props #Get the reference to the body properties
		fuse_props['pos'].params[''] = PGMparameter(4,3, order_u=3) #This command states that we will have 4 points (linearly spaced between u=0.0 and u=1.0) to specify the fuselage center-line control points. In addition, the interpolation between the sections will be parabolic due to the order_u=3 option.
		fuse_props['nor'].params[''] = PGMparameter(1,1) #This controls the normal of the faces with respect to the center-line
		fuse_props['scl'].params[''] = PGMparameter(4,3) #We will specify the radius along 4 linearly spaced sections and in 3 directions (x,y,z, although the z one is ignored) to make an ellipse. Linear interpolation will be used for the radius of intermediate sections.
		fuse_props['flt'].params[''] = PGMparameter(1,1) #This parameter steadily turns elliptical cross-sections into rectangular cross-sections.

		wing_l_props = self.comps['wing_l'].props #Get the reference to the left wing parameters
		wing_l_props['pos'].params[''] = PGMparameter(2,3) #We need to set 3 coordinates (x,y,z) for the root section and the tip section (root + tip, hence the 2)
		wing_l_props['scl'].params[''] = PGMparameter(2,1) #We will specify the chord of the root and the tip with the scaling parameter (1 because the same scaling is applied to x,y,z, though z is irrelevant)
		wing_r_props = self.comps['wing_r'].props #Get the reference to the right wing parameters
		wing_r_props['pos'].params[''] = PGMparameter(2,3)
		wing_r_props['scl'].params[''] = PGMparameter(2,1)

		htail_l_props = self.comps['htail_l'].props #Get the reference to the left tail parameters
		htail_l_props['pos'].params[''] = PGMparameter(2,3)
		htail_l_props['scl'].params[''] = PGMparameter(2,3)
		htail_l_props['rot'].params[''] = PGMparameter(2,3) #Now we use 3 parameters to set rotations around 3 axes (x,y,z) for the root and tip sections
		htail_r_props = self.comps['htail_r'].props #Get the reference to the right tail parameters
		htail_r_props['pos'].params[''] = PGMparameter(2,3)
		htail_r_props['scl'].params[''] = PGMparameter(2,3)
		htail_r_props['rot'].params[''] = PGMparameter(2,3)
		vtail_props = self.comps['vtail'].props #Get the reference to the vertical tail parameters
		vtail_props['pos'].params[''] = PGMparameter(2,3)
		vtail_props['scl'].params[''] = PGMparameter(2,3)
		vtail_props['nor'].params[''] = PGMparameter(2,3) #We will turn on the normal condition for the root and tip sections
		vtail_props['rot'].params[''] = PGMparameter(2,3)

	#DEFINE THE PARAMETERS VALUES
	def _compute_params(self):
		fuse_props = self.comps['fuse'].props #Get the reference to the body properties
		fuse_props['pos'].params[''].data[:,:] = [[0,0,0],[0.210,0,0],[0.565,0,0],[1.340,0,0]] #Sets the position of the control points
		fuse_props['nor'].params[''].data[0] = 1.0 #Set the faces to be perpendicular to the center line axis
		fuse_props['scl'].params[''].data[:,:] = [[0.015,0.015,1],[0.105,0.063,1],[0.105,0.063,1],[0.001,0.001,1]] #Sets the radius for each direction (x,y,z). We use scaling in z as 1 to avoid any type of scaling in this axis.
		fuse_props['flt'].params[''].data[:,:] = 0.98 #1 turns the section in a complete square, and 0 remains elliptical. We used 0.98 to get almost rectangular corners.

		wing_l_props = self.comps['wing_l'].props #Get the reference to the left wing parameters
		wing_l_props['pos'].params[''].data[:,:] = [[0.228,0.025,0.105],[0.228,0.025,0.725]] #We need to set 3 coordinates (x,y,z) for the root and the tip. The root is exactly on the surface of the fuselage.
		wing_l_props['scl'].params[''].data[:,:] = [[0.256],[0.256]] #We will specify the chord of each section
		wing_r_props = self.comps['wing_r'].props #Get the reference to the right wing parameters
		wing_r_props['pos'].params[''].data[:,:] = [[0.228,0.025,-0.725],[0.228,0.025,-0.105]]
		wing_r_props['scl'].params[''].data[:,:] = [[0.256],[0.256]]

		htail_l_props = self.comps['htail_l'].props #Get the reference to the left tail parameters
		htail_l_props['pos'].params[''].data[:,:] = [[1.059,0,0.040],[1.118,0,0.219]]
		htail_l_props['scl'].params[''].data[:,:] = [[0.238,0.03,1],[0.154,0.03,1]] #We will specify the chord and thickness of each section
		htail_l_props['rot'].params[''].data[:,:] = [[0,7.943,0],[0,0,0]] #We will only use rotation around the y axis for the root section. The rotation is about the leading edge.
		htail_r_props = self.comps['htail_r'].props #Get the reference to the right tail parameters
		htail_r_props['pos'].params[''].data[:,:] = [[1.118,0,-0.219],[1.059,0,-0.040]]
		htail_r_props['scl'].params[''].data[:,:] = [[0.154,0.03,1],[0.238,0.03,1]]
		htail_r_props['rot'].params[''].data[:,:] = [[0,0,0],[0,-7.943,0]]
		vtail_props = self.comps['vtail'].props #Get the reference to the vertical tail parameters
		vtail_props['pos'].params[''].data[:,:] = [[1.106,0.021,0],[1.1466,0.271,0]]
		vtail_props['scl'].params[''].data[:,:] = [[0.164,0.03,1],[0.106,0.03,1]]
		vtail_props['nor'].params[''].data[:,:] = [[1,0,0],[1,0,0]] #Make the sections perpendicular to the vertical tail center-line
		vtail_props['rot'].params[''].data[:,:] = [[0,4.574,0],[0,0,0]]

		return [], [], [] #No gradient needs to be returned because we have no design variables

##MAIN ROUTINE
if __name__ == '__main__':

	#INITIALIZE GEOMETRY
	pgm = test() #Create the aircraft class
	bse = pgm.initialize() #Initialize the aircraft and return a pointer to the BSE object 

	#ADD AIRFOILS
	pgm.comps['wing_l'].set_airfoil() #This sets a NACA0012 to the left wing
	pgm.comps['wing_r'].set_airfoil() #This sets a NACA0012 to the right wing
	pgm.comps['htail_l'].set_airfoil() #This sets a NACA0012 to the left tail
	pgm.comps['htail_r'].set_airfoil() #This sets a NACA0012 to the right tail
	pgm.comps['vtail'].set_airfoil() #This sets a NACA0012 to the vertical tail

	#COMPUTE GEOMETRY
	pgm.compute_all()

	#EXPORT GEOMETRY
	bse.vec['pt_str']._hidden[:] = False #Show both sides of the aircraft
	bse.vec['pt_str'].export_tec_str()
