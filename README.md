# R3B_Neutron_Analyzer

This simple macro is used to derive output from the simulation that can be compared directly to the data in a fast way. It supports simulation involving Land or Neuland.

Before to run, source the config.sh of your R3BRoot install.

You might also want to check some paths.

To run :

>make

>./invariant_mass --file=sim.root --NWall=Land //with LAND

OR

>./invariant_mass --file=sim.root --NWall=Neuland //with Neuland

--file= should be followed by the name of the file to analyze in the input/ directory
