This is mhd2d, a simple tutorial 2-dim MHD finite volume code.

To compile type 'make'.
The code should compile with only one or two warnings about unused functions.


To run, type './mhd2d'.
It should produce a number of vtk data files in the subdirectory 'data'.


To inspect the data, start 'visit'.
Select 'File -> Open file...', browse into the 'data' directory and open 
the 'MHD_*.vtk' database.
Then, select 'Plots -> Pseudocolor -> rho' to display the mass density rho.
You can select other variables for display under the 'Variables' menu.
Cycle through the time steps using the time slider and/or the forward/backward buttons.


To change simulation parameters, modify the file 'global.cc', *save the file*,
'make' and run the code again.
You can also change the name stem for the output files to create multiple file sets.

To modify the code, change the corresponding *.cc and/or *.h files and redo the
'make'-run cycle.
