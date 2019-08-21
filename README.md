# rdf-com-xyz
This program calculates COM radial distributions for molecules in XYZ file.
# cp2k-rdf
This fortran program reads xyz file and calculates center of mass radial distributions.
The center of mass coordinates are calculated automatically. The rdf.par file is required 
to include a number of molecular types, a number of atoms for each molecule
in order of xyz file, a distance withing the radial distributions should be estimated,
a number of bins, and box dimensions.
# compile
  In order to compile the program: gfortran rdf_com.f
