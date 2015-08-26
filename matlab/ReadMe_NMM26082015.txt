First Run directions:

1) set up paths
-open the add_paths.m file
-change pathnames to those appropriate to your operating system (‘/‘ for mac, ‘\’ for windows) and directory where you have placed the gitHub directory
-run add_paths

2) parameter definitions
- The parameters for the simulation are defined in the CCMParams.m file which is a class definition. There are 3 types of parameters in this file
	a) there are single value “properties,” with units, which can be changed 
	b) properties which are protected and cannot be changed (Constant) tagged
	c) dependent parameters, which are calculated from other parameters defined in the file
	d)abstract parameters which are calculated in other files (this allows you to define them differently depending on whether you have a compartment or not), many of these are non-dimensional — which was useful for calculating analytic solutions.

- The non-dimensional abstract parameters are calculated using either CCMParams_Csome.m for when there is a compartment, and CCMParams_NoCsome.m when there is no compartment. These two class structures are super classes which include CCMParams, and are what are actually used by the rest of the code. Go here fore information about classdefs in matlab: 
http://www.mathworks.com/help/matlab/object-oriented-programming-in-matlab.html

3) numerical code
- 

