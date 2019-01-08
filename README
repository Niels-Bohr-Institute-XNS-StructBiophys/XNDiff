/*
   folder structure:
   XNDiff/ : 		XNDiff.cpp mtrand.cpp mtrand.h compile_XNDiff.sh Makefile
   XNDiff/in/ :		syminfo.lib (optional) f0_InterTables.dat DeBe_NeutronNews.dat atomweights.dat cif_core.dic 
                	further user-provided input files like cif-files (less common pdb-files), symmetry operation files, parameter files (less common, parameters can be provided upon call of XNDiff on cmd)
   XNDiff/out/ : 	XNDiff generated files for SAXS and SANS patterns, log-files, further output files for testing (distributions, (modified)cif-files, ...)

*/


/* 
   mandatory and optional external libraries / installations to run XNDiff on GNU Linux and (with minor code modifcations) on other POSIX systems (FreeBSD, ...):

   -MTRand PRNG (MANDATORY, might be replaced with a PRNG from GSL later),
    http://www.bedaux.net/mtrand/
    just place mtrand.cpp and mtrand.h at the same location of XNDiff.cpp before compilation

   -GSL (MANDATORY, tested with v1.15)
    http://www.gnu.org/software/gsl/ 
    comes along with most GNU Linux distros and is installed on most HPC systems by default

   -CCP4 (OPTIONAL, tested with v6.1.13, when using -sym <ccp4 ...> option)
    http://www.ccp4.ac.uk/
    Note: CCP4 can be binary installation, since CCP4 is called via opening a shell
 
   -Armadillo (OPTIONAL, tested with v3.0.3, mostly for test purposes, when using -bt <bt_twin, bt_double> option)
*/


/* compilation using g++, assuming everything is installed / put in the right place:

   -without using optional Armadillo library, Optimization: DEBUG / O3-Level, no further architecture-specific optimization
    e.g. ./compile_XNDiff.sh g++ DEBUG 0 NONE
    e.g. ./compile_XNDiff.sh g++ 3 0 NONE
   -using optional Armadillo library, Optimization: DEBUG / O3-Level, no further architecture-specific optimization
    e.g. ./compile_XNDiff.sh g++ DEBUG 1 NONE
    e.g. ./compile_XNDiff.sh g++ 3 1 NONE

*/


/*
   further mandatory / optional parameter-files with physical constants and symmetry operations etc to run XNDiff:

   -download from link below and save as syminfo.lib to in/
    http://www.ccp4.ac.uk/cvs/viewvc.cgi/libccp4/data/syminfo.lib?view=co
    (OPTIONAL, symmetry operations, required when using -sym <symopfile> or -sym <ccp4 ...> option without providing a symopfile, see source code for further details)
    (further reading http://www.ccp4.ac.uk/html/symlib.html)

   -download from link below and save as f0_InterTables.dat to in/
    http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/f0_InterTables.dat
    (MANDATORY, Cromer-Mann coefficients for X-ray scattering)

   -download from link below and save as DeBe_NeutronNews.dat to in/
    http://ftp.esrf.eu/pub/scisoft/DabaxFiles/DeBe_NeutronNews.dat
    (MANDATORY, neutron coherent scattering lengths etc)

   -save as atomweights.dat the ASCII table from link below into in/ and make some modifications (has been already done and included in your download of XNDiff!!!)
    http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii&isotype=some
    (MANDATORY, atom/isotope masses)

   -download from link below the current version and save as cif_core.dic to in/
    ftp://ftp.iucr.org/pub/cif_core.dic
    (MANDATORY, cif-dictionary required to read cif-files, tested with v2.4.3 from 2012-05-16)
*/
