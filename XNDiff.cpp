/***************************************************************************
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


/*** History: ***/
/*
  The development of a program for the evaluation of very complicated functions
  was started in Oktober 2000 by T. Unruh. This program was called "func_calc".
  The main goal was to write a program for the calculation of diffraction patterns
  of very complicated structures numerically.
  After first qualitatively successful tests in 2000 first quantitatively convincing
  results could be achieved in 2007. At this stage the program was renamed to "XDiff"
  and made available by Sourceforge.net ( http://sourceforge.net/projects/xdiff/).

  15.10.2000: func_calc Version 0.02
              av_corr applied to S not to E !
  28.10.2000: Version 0.03
              karthesian coordinates of a vectors corrected
	      hkl averaging (mode 1) bugs corrected
	      dispersion medium correction added
  25.11.2000: func_calc Version 0.04
              some speed optimizations
	      pid included in outputs
  27.11.2000: func_calc Version 0.05
              some changes of output data
              log file name changed to <mfname>.log
   3.11.2000: func_calc Version 0.06
              V = V * n1 * n2 * n3!
              some minor bugs
   4.11.2000: func_calc Version 0.07
              minor bug fix in output file
  21.12.2000: func_calc Version 0.08
              minor speed optimizations 
  11.01.2005: func_calc Version 0.09
              background bug fixes
  21.02.2006: func_calc Added stabilizer scattering
  17.03.2006: func_calc Added new crystal (plus stabilizer) overlap test
  26.06.2007: XDiff version 0.01 extracted from func_calc
              (cf. J. Appl. Cryst. submitted)
              The code is very preliminary and only intended to simulate very special
	      systems so far
  26.06.2009: Begin of further developement of XNDiff
		- add neutron diffraction support
		- optimization for speed
		- add crystallographic stuff such as space groups and symmetry operations
		- add a tool for the generation of structure parameter files with complete molecules
		  from e.g. CSD (cif) or PDB file format

  20140406 minor bug fix, for New crystal: the printed w3 and wt3 vectors were wrong w(t)3[0][k][:] was printed (i.e. always n3=1) but not w(t)3[n3[k]-1][k][:] for the correct n3[k] thickness
           only printed results in the log-files were wrong but not the computations
*/



/**********************************************/
/* short description of the stackcpp function */
/**********************************************/

/* (to be added later) */





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




#include <iostream>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <vector>
#include <ctype.h>
#include <time.h>

/* directory exists http://bytes.com/topic/c/answers/584434-check-directory-exists-c */
#include <dirent.h>

/* mkdir http://www.imb-jena.de/~gmueller/kurse/c_c++/c_file.html */
#include <sys/stat.h>

/* Unix-Default-Lib-Header (Windows C++ stdlib.h ) */
#include <unistd.h>

/* for typeid */
#include <typeinfo>

/**************************************/
/* OpenMP support for parallelization */
/**************************************/
#include <omp.h> 

/* print variable name http://ubuntuforums.org/showthread.php?t=506402 */
#define var_name(var) #var

/**************************************************************************/
/* COMPLEX NUMBER HEADER OR CLASSES http://en.wikipedia.org/wiki/Complex.h*/
/**************************************************************************/
#include <complex>

/***************************************************************************/
/* PRNG Mersenne-Twister MTRand for MC mode (much better than C++'s rand())*/
/*                    http://www.bedaux.net/mtrand/                        */
/***************************************************************************/
#include "mtrand.h"

/*****************************************************************/
/* GNU Scientific Library (GSL), PRNG for discrete distributions */
/*             http://www.gnu.org/software/gsl/                  */
/*****************************************************************/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>




using namespace std;

/****************************************/
/* ARMADILLO LIBRARY FOR LINEAR ALGEBRA */
/*      http://arma.sourceforge.net/    */
/****************************************/
#if (ARMADILLO)
#include <armadillo>
#endif
#if (ARMADILLO)
using namespace arma;
#endif



typedef complex<double> dcmplx;




class XNDiff
{

  public:
	/* program information */  
	const char* XNDIFF_VERSION ; /* see constructor */
	const char* XNDIFF_DATE ; /* see constructor */
	const char* XNDIFF_AUTHOR ; /* see constructor */
	double DR ; /* see constructor */
	double RD ; /* see constructor */

	/* List of available functions for XPPSA Method */
	/* struct function, with func[i]->name and func[i].nop */
	/* List of available functions for basis transformations */
	static const unsigned int NOF = 1 ;
	static const unsigned int NOBT = 2 ;
	bool bt_flag ;
	char bt[1024] ; /* bt function name */
	struct function
	{
		const char *name ; /* Function name */
		unsigned int nop ; /* Number of parameters */
	} func_XPPSA[NOF], func_bt[NOBT] ;

	struct par_entries
	{
		const char *THICKNESS_INLAY ; 
		const char *THICKNESS_OUTLAY ;
		const char *XRAY_RHO_INLAY ;
		const char *XRAY_RHO_OUTLAY ;
		const char *XRAY_RHO_DM ;
		const char *NEUT_SLD_INLAY ;
		const char *NEUT_SLD_OUTLAY ;
		const char *NEUT_SLD_DM ;
		const char *NEUT_ICSD_INLAY ;
		const char *NEUT_ICSD_OUTLAY ;
		const char *NEUT_ICSD_DM ;
		const char *CONC_CRY ;
		const char *CONC_DM ;
	} par_keyword ;

	/* process ID http://www.tutorials.de/forum/linux-unix/158373-datentypen-pid_t-dev_t-ssize_t.html */
	pid_t pid ;
	unsigned int time_flag ; /* no protocolling(0), only at beginning/FFF/terminating(1) or detailed(2) */
	unsigned int mem_flag ; /* no protocolling(0), only at beginning(1) or detailed(2) */
	char mem_filename[64] ;

	/* log file and output files */
	bool singlefilesoutput_flag ;

	bool extension_flag ;

	bool intermediatesteps_flag ;
	unsigned int imsteps ;

	int log_flag;
	FILE *logfile ;

	char mfname[1024] ; /* Master file name of output files */	

	char **evarg ; /* Export function call parameters */
	int *ecarg ; /* Export function call parameters */


	/* XNDiff PRNG for particles, ... */
	#define R1 16807
	#define R2 2147483647
	#define R3 (1.0/R2)
	#define R4 127773
	#define R5 2836
	#define R6 32
	#define R7 (1+(R2-1)/R6)
	#define R8 1.1e-7
	#define R9 (1.0-R8)


	/* Parameters are contained in the structure par of type struct stack_parameters */
	struct stack_parameters
	{       /* shared parameters */
		double n11, n12 ;
		double n21, n22 ;
		double n31, n32 ;
		double d11, d12 ;
		double d21, d22 ;
		double D1, D2 ;
		int td1 ;
		int td2 ;
		int stackmode ;
		double rs1, rs2 ;
		unsigned int np ;
		unsigned int nr ;
		double h, k, l ;
		double hs, ks, ls ;
		unsigned int nsp ; /* number of layers in single particles e.g. 10 */
		unsigned int nms ; /* maximum number of particles in stacks e.g. 5 */

		int av_mode ;
		double av_pol[3] ;
		double av_azi[3] ;
		double r31, r32 ;
		double av_r3[2] ; /* range of r3-rotation in [ av_r3[0] , av_r3[1] ] */
		unsigned int n_r3 ; /* number of grid points for r3-rotation between -90° and +90°, e.g. 180 */
		unsigned int n_MC ; /* number of points for MC mode e.g. 180*24 = 4320 */
	} *par ;

	static const int MAX_AV_MODE = 3 ; /* see +av flag */

	unsigned int nav_pc, nav_ac, nav_r3, nav_MC ;
	double av_r3d ;
	double ds, s1 ;
	
	vector <double> av_pol_ang ;
	vector <double> av_azi_ang ;

	char par_file[1024] ;		/* filename of the .par-file which provides the necessary information about the stabilizer layer and dispersion media */
	char ccp4_config_file[1024] ;	/* configuration file with path information for an installed CCP4-library */

	char ccp4_bash_file[1024] ; /* filename of the bash-script used for running the CCP4 (pdbset) */
	const char *def_ccp4_bash_file ;

	char cif_dic[1024] ;		/* cif-dictionary filename */
	char cif_file[1024] ;		/* cif-file filename */

	char cif_pdb_file[1024] ;
	const char *def_cif_pdb_file ;

	bool jmol_cif_flag ; /* write cif-file for molecular viewers like Jmol,PyMol,Rasmol,... */
		
	bool sym_pdb_mode ;
	char pdb_file[1024] ;		/* pdb-filename either by default (def_pdb_file) generated by CCP4 from cif_pdb_file or user-supplied */
	const char *def_pdb_file ;

	bool pcr_flag ;			/* flag for writing a file with atomic coordinates in the pcr-file format for FullProf */
	char pcr_file[1024] ;

	struct cif_dic_entries		/* structure of entries in the cif-dictionary */
	{
		bool isnumeric;
		char* units;
		char* category;
		char** names;
		int num_names;
	} *cde;
	unsigned int num_cif_dic_entries;	/* number of entries in the cif-dictionary */

	struct cif_file_entries		/* structure of entries in the cif-file */
	{
		bool isnumeric ;
		char* units ;
		char* name ;
		unsigned int data_length ;
		char** text_data ;
		double* numeric_data ;
		double* numeric_stderr_data ; /* default 0.0 */
	} *cfe ;
	unsigned int num_cif_file_entries ;	/* number of entries in the cif-file */

	bool sym_ccp4_mode ;
	bool sym_cif_mode ;

	bool symop_userdef ; /* flag whether to use a user-supplied symmetry operations file */
	char symop_file[1024] ; /* filename of the user-supplied symmetry operations file */

	char** symopxyz; /* string array of symmetry operations */
	unsigned int num_symop; /* amount of symmetry operations */
	char** cenopxyz; /* string array of center operations */
	unsigned int num_cenop; /* amount of center operations */

	vector<double> cis ; /* c_i's fitted from diluted systems (\sum_i c_i will be normalized to 1 after read in) */
	vector<double> nis ; /* n_i = (c_i/i)/[\sum_i(c_i/i)] in case cis are provided, else if nis provided normaliyed to 1 */
	bool cis_isdef ;
	bool nis_isdef ;

	vector< vector<unsigned int> > bond_site12 ; /* connectivity table derived from cif-file */
	vector< vector<string> > bond_symmetry12 ;
	vector<double> bond_distance ;

	/* Standard Atomic Weight for elements and Relative Atomic Mass for isotopes */
	struct StAtWtandRelAtMa
	{
		char* nuclei_name;
		double mass; /* Standard Atomic Weight for elements and Relative Atomic Mass for isotopes */
	} *atom_mass;
	int num_atom_mass_entries;

	/* NEutron SCattering LEngth and CRoss SECtions */
	struct nesclecrsec
	{
		char* nuclei_name;
		double nat_abd; /* natural abundance [0-100 %] (for isotopes the half-life is given instead) */
		dcmplx coh_sc_length; /* complex bound coherent scattering length [fm] */
		dcmplx incoh_sc_length; /* complex bound incoherent scattering length [fm] */ 
		double coh_sc_cross; /* bound coherent scattering cross section [barn] */
		double incoh_sc_cross; /* bound incoherent scattering cross section [barn] */ 
		double total_sc_cross; /* total bound scattering cross section [barn] */ 
		double abs_cross; /* absorption cross section for 2200 m/s neutrons [barn] */
		double mass; /* Standard Atomic Weight for elements and Relative Atomic Mass for isotopes */
	} *n_atom_sc;
	int num_n_atom_sc_entries;

	double cell_par[6] ; 	       /* Unit cell parameters, edge lengths and angles of the triclinic unit cell in [nm] */
	struct atom_entries            /* Atom parameters */
	{
		char label[10] ;       /* Atom label */
		char type[5] ;         /* Atom type (cf. get_atomic_scattering_factor */
		double coord[3] ;      /* Atom fractional coordinates */
		double mult ;          /* Multiplicity */
		double temp[7] ;       /* 1 isotropic and 6 anisotropic temperature factors (optional) */
		dcmplx coh_sc_length ; /* coherent scattering length */
		double incoh_sc_cross ;/* incoherent scattering cross section */
		int f0ind ; /* index to a table with the pre-computed scattering amplitudes as a function of s, f_0(s) */
	} *atom ;
	unsigned int noa ;                      /* Number of atoms in the unit cell */

	double** f0table ; /* numf0 x par->np table containing f_0^j(s), j=0,...,numf0 */
	int numf0 ;

	bool par_num_def ;
	double *rho_dm ; /* electron density of dispersion medium [electrons/nm^3] */
	double *sld_dm ; /* neutron SLD of dispersion medium [1e-6 A^(-2)] */
	double *icsd_dm ; /* neutron ICSD of the dispersion medium [1/cm] */

	double thickness_isl ; /* Thickness of inner stabilizer layer [nm] */
	double *rho_isl ; /* electron density of inner stabilizer layer [electrons/nm^3] */
	double *sld_isl ; /* neutron SLD of inner stabilizer layer [1e-6 A^(-2)] */
	double *icsd_isl ; /* neutron ICSD of inner stabilizer layer [1/cm] */

	double thickness_osl ; /* Thickness of outer stabilizer layer [nm] */
	double *rho_osl ; /* electron density of outer stabilizer layer [electrons/nm^3] */
	double *sld_osl ; /* neutron SLD of outer stabilizer layer [1e-6 A^(-2)] */
	double *icsd_osl ; /* neutron ICSD of outer stabilizer layer [1/cm] */

	double rho_cry ; /* X-ray electron density for the crystal [electrons/nm^3] */
	double icsd_cry ; /* neutron incoherent cross section density for the crystal [1/cm] */
	double sld_cry ; /* neutron scattering length density for the crystal [10^(-6) 1/A^2] */

	bool conc_userdef ; /* flag whether to use default concentrations for conc_dm and conc_cry ( false, 1.0 and 0.01) or those specified in the par-file (true) */
	double conc_dm ; /* vol% concentration of dispersion medium in the nanosuspension */
	double conc_cry ; /* vol% concentration of crystals in the nanosuspension */

	struct rho_multi_struct
	{
		bool multiq ;
		unsigned int n ;
		char* XRAY_RHO_INLAY_file ;
		unsigned int XRAY_RHO_INLAY_col ;
		char* XRAY_RHO_OUTLAY_file ;
		unsigned int XRAY_RHO_OUTLAY_col ;
		char* XRAY_RHO_DM_file ;
		unsigned int XRAY_RHO_DM_col ;
	} rho_multi ;

	struct sld_multi_struct
	{
		bool multiq ;
		unsigned int n ;
		char* NEUT_SLD_INLAY_file ;
		unsigned int NEUT_SLD_INLAY_col ;
		char* NEUT_SLD_OUTLAY_file ;
		unsigned int NEUT_SLD_OUTLAY_col ;
		char* NEUT_SLD_DM_file ;
		unsigned int NEUT_SLD_DM_col ;
		char* NEUT_ICSD_INLAY_file ;
		unsigned int NEUT_ICSD_INLAY_col ;
		char* NEUT_ICSD_OUTLAY_file ;
		unsigned int NEUT_ICSD_OUTLAY_col ;
		char* NEUT_ICSD_DM_file ;
		unsigned int NEUT_ICSD_DM_col ;
	} sld_multi ;

	double r_e ; /* see constructor */
	double sc_scl ; /* see constructor */
	double sc_rho ; /* see constructor */
	double sc_sld ; /* see constructor */
	double sc_icsd ; /* see constructor */

	double sc_irr ; /* see constructor */

	const char *unit_S_X, *unit_S_n ; /* see constructor */
	const char *unit_rho, *unit_sld, *unit_icsd ; /* see constructor */
	const char *unit_F_X, *unit_F_n ; /* see constructor */
	const char *unit_Yc_X[10], *unit_Yc_n[10], *unit_Yi_n[4] ; /* see constructor */
	const char *unit_V_CRY, *unit_V_ISL, *unit_V_OSL, *unit_V_DM, *unit_V_IRR ; /* see constructor */

	int Xn ;  /* flag, that decides whether X-ray (Xn=-1), neutron (Xn=+1) or both patterns (average of both, Xn=0) should be computed */
	dcmplx ***FFF_X ; /* structure amplitude for X-rays [nm] */
	dcmplx ***FFF_n ; /* structure amplitude for neutrons [nm] */

	bool orav_flag ;

	bool openmp_flag ;

	bool Fr_flag ;
	int num_Fr_flag ;
	char Fr_file[2][1024] ;
	bool Fw_flag ;
	bool FNoMem_flag ;

	bool distr_flag ;

	bool weight_flag ;

	int init_n_rand ;
	bool init_n_rand_userdef ;
	long unsigned int init_d_rand ;
	bool init_d_rand_userdef ;
	unsigned int init_MTRand ;
	bool init_MTRand_userdef ;

	double ***eQvj ; /* scalar products for all three a_j's with all the unit vectors of Q depending on theta,phi/r3 */

	double a1[3], a2[3], a3[3] ; /* Lattice parameters in a Cartesian coordinate system */
	double b1[3], b2[3], b3[3] ; /* Reciprocal unit cell vectors in kartesian coordinate system */
	double G0[3] ; /* Reciprocal lattice unit vector parallel to (s-s0) */
	double kg1[3], kg2[3] ; /* Two vectors in (hkl)-plane */	
	double Gs[3] ; /* Reciprocal lattice unit vector of stacking direction */
	double k1[3], k2[3] ; /* Two vectors in (hs ks ls)-plane */
	double V_uc ; /* Volume of unit cell = det|a_1 a_2 a_3| [nm^3] */

	char startingtime[64];

	const char *def_ifo, *def_ofo ;
	char ifo[1024], ofo[1024] ;

	bool overwrite_flag ;
	bool silent_flag ;

	bool test_XnEquiv_flag ;
	bool test_I_at_s_equal_zero_flag ;

	public:


	/* constructor */
	XNDiff()
	{
		/* program informations */
		this->func_XPPSA[0].name = "stackcpp" ;
		this->func_XPPSA[0].nop = 27 ;
		this->func_bt[0].name = "bt_twin" ;
		this->func_bt[0].nop = 0 ;
		this->func_bt[1].name = "bt_double" ;
		this->func_bt[1].nop = 0 ;
		this->XNDIFF_VERSION = "0.4.1" ;
		this->XNDIFF_DATE = "8.2000-04.2014" ;
		this->XNDIFF_AUTHOR = "Tobias Unruh and Martin Schmiele" ;

		/* class constants */
		this->DR = M_PI/180.0 ;
		this->RD = (180.0/M_PI) ;

		this->r_e = 2.8179402894 ;  /* classical electron radius in [fm] */
		this->sc_scl = 1.0e-6 ; /* scattering lengths / r_e in [fm] -> [nm] for F_X/n and sc_rho */
		this->sc_rho = r_e * sc_scl ; /* [ 1/nm^3 * ( r_e = 2.82 fm ) ] -> [1/nm^2] */
		this->sc_sld = 1.0e-4 ; /* [10^(-6) 1/A^2 = 10^10 1/cm^2] -> [1/nm^2] */
 		this->sc_icsd = 1.0e-7 ; /* incoherent cross section density (ICSD) in [1/cm], volume in [nm^3], aim is [nm^2] for B_n */
		this->sc_irr = 1.0e-7 ; /* [nm^3] for V -> [nm^2 * cm] for Vtot_irr */

		/* units */
		this->unit_S_X = "1/cm" ; /* absolute units */
		this->unit_S_n = "1/cm" ; /* absolute units */
		this->unit_rho = "1/nm^3" ; /* user-friendly unit for rho, note that r_e * rho can be expressed as for neutrons in [10^10 1/cm^2] */
		this->unit_sld = "10^(-6) 1/A^2" ; /* user-friendly unit for sld, equivalent to [10^10 1/cm^2] */
		this->unit_icsd = "1/cm" ; /* absolute units */
		this->unit_F_X = "nm" ; /* atom form factor F */
		this->unit_F_n = "nm" ;

		this->unit_V_CRY = "nm^3" ;
		this->unit_V_ISL = "nm^3" ;
		this->unit_V_OSL = "nm^3" ;
		this->unit_V_DM = "nm^3" ;
		this->unit_V_IRR = "nm^2 cm" ;

		this->unit_Yc_X[0] = "1/cm" ; /* absolute units */
		for ( int q=1; q<7; ++q) { this->unit_Yc_X[q] = "nm^6/cm" ; }  /* absolute units */
		for ( int q=7; q<10; ++q) { this->unit_Yc_X[q] = "nm^3/cm" ; }  /* absolute units */

		this->unit_Yc_n[0] = "1/cm" ; /* absolute units */
		for ( int q=1; q<7; ++q) { this->unit_Yc_n[q] = "nm^4/cm" ; }  /* absolute units */
		for ( int q=7; q<10; ++q) { this->unit_Yc_n[q] = "nm^2/cm" ; }  /* absolute units */

		this->unit_Yi_n[0] = "1/cm" ; /* absolute units */
		for ( int q=1; q<4; ++q) { this->unit_Yi_n[q] = "1" ; }  /* absolute units */

		this->par_keyword.THICKNESS_INLAY = "THICKNESS_INLAY" ;
		this->par_keyword.THICKNESS_OUTLAY = "THICKNESS_OUTLAY" ;
		this->par_keyword.XRAY_RHO_INLAY = "XRAY_RHO_INLAY" ;
		this->par_keyword.XRAY_RHO_OUTLAY = "XRAY_RHO_OUTLAY" ;
		this->par_keyword.XRAY_RHO_DM = "XRAY_RHO_DM" ;
		this->par_keyword.NEUT_SLD_INLAY = "NEUT_SLD_INLAY" ;
		this->par_keyword.NEUT_SLD_OUTLAY = "NEUT_SLD_OUTLAY" ;
		this->par_keyword.NEUT_SLD_DM = "NEUT_SLD_DM" ;
		this->par_keyword.NEUT_ICSD_INLAY = "NEUT_ICSD_INLAY" ;
		this->par_keyword.NEUT_ICSD_OUTLAY = "NEUT_ICSD_OUTLAY" ;
		this->par_keyword.NEUT_ICSD_DM = "NEUT_ICSD_DM" ;
		this->par_keyword.CONC_CRY = "CONC_CRY" ;
		this->par_keyword.CONC_DM = "CONC_DM" ;

		this->num_symop = 0 ;
		this->num_cenop = 0 ;
		this->noa = 0 ;
	}


	/* check if a directory path exists */
	bool directory_exists( const char* path)
	{
		if ( path == NULL) return false ;
		
		DIR *dir ;
		bool ret = false ;
		
		dir = opendir(path) ;
		
		if ( dir != NULL )
		{
			ret = true ;
			(void) closedir(dir) ;
		}
		return ret ;
	}


	/* Test computer code if it returns the same results for neutrons as for X-rays when replacing neutron parameters with X-ray ones
	   replaced will be neutron slds by rhos and scaling factors
	   FFF_n will be defacto copied from FFF_X elsewhere in the code, see:

		compute_structure_amplitude()
		read_structure_amplitude()
		stackcpp()
	*/
	void test_XnEquiv()
	{
		/* test works only for X-ray AND neutron scattering */
		if ( Xn != 0 ) { XNDIFF_ERROR(82) ; }

		/* test should better not be applied when using multiple rho/sld/icsd's '*/
		if ( rho_multi.multiq || sld_multi.multiq ) { XNDIFF_ERROR(83) ; }

		/* coherent neutron scattering -> X */
		this->sc_sld = this->sc_rho ;
		sld_isl[0] = rho_isl[0] ;
		sld_osl[0] = rho_osl[0] ;
		sld_dm[0] = rho_dm[0] ;

		/* no incoherent scattering */
		this->sc_icsd = 0.0 ; 
		icsd_isl[0] = 0.0 ;
		icsd_osl[0] = 0.0 ;
		icsd_dm[0] = 0.0 ;

		/* units cannot be replaced because the of const char* */
	}


	/* Test computer code if it returns the correct intensity at I(s=0)

	   Use :
	   conc_cry = 1.0 (and optionally conc_dm = 0.0)
	   to normalize only by volume of crystals (without stabilizer layer)
	   Note that in XNDiff Vtot_irr ~ Vtot_cry / conc_cry, hence only conc_cry counts, 
	   but not conc_dm, which is only relevant for Vtot_dm (via VN_dm) ~ Vtot_cry * conc_dm / conc_cry 
	   and by this for B_n and Yi_n
	   Hence with conc_cry = 1.0 it holds Vtot_irr = Vtot_cry

	   Generally normalization is done in XNDiff via :
	   S_X/n -> S_X/n / MULT(coh) / SUM_WEIGHT_CORR / Vtot_irr
	   Yc_X/n -> Yc_X/n / MULT(coh) / SUM_WEIGHT_CORR / Vtot_irr
	   and 
	   B_n -> B_n / MULT(inc) / Vtot_irr
	   Yi_n -> Yi_n / MULT(inc) / Vtot_irr
	  
	   In exported files for S_X/n and Yc_X/n etc 
	   V_CRY ~ Vtot_cry ~ volume of crystal
	   V_ISL ~ Vtot_isl ~ volume of crystal AND inner stabilizer layer
	   V_OSL ~ Vtot_osl ~ volume of crystal AND inner AND outer stabilizer layers
	   V_DM and V_IRR correspond to Vtot_dm and Vtot_irr as given from above via Vtot_cry
	   
	   V_CRY corresponds physically to the volume of crystal and therefore for XNDiff or later the correct conc_cry should be applied !!!
	   ( V_ISL - V_CRY )
	   and 
	   ( V_OSL - V_ISL )
	   correspond to volumes of stabilizer layer covering the parallelepipedal particles. Ideally they should be also related to the concentration of the 
	   stabilizer layer (e.g. phospholipid + sodium glycocholate etc ). However they must not, cause the stabilizer layer must not be densely packed and dispersion
	   medium and crystal substance may contribute. Their contributions might be later estimated by the rho/sld and from the difference between the expected volume
	   concentration of the stabilizer (sample composition) and the one from the simulation/fit.
	
	*/
	void test_I_at_s_equal_zero()
	{
		/* NOT IMPLEMENTED YET */

	}


	#if (ARMADILLO)
	/* basis transformation for twinned triglyceride crystals */
	void bt_twin()
	{
		char sdummy[1024] ;
		int idummy ;

		vec a1_vec = zeros<vec>(3) ;
		vec a2_vec = zeros<vec>(3) ;
		vec a3_vec = zeros<vec>(3) ;

		vec delta_vec = zeros<vec>(3) ;

		vec ap1_vec = zeros<vec>(3) ;
		vec ap2_vec = zeros<vec>(3) ;
		vec ap3_vec = zeros<vec>(3) ;

		mat A = zeros<mat>(3,3) ; /* A=(a1,a2,a3) */
		mat Ap = zeros<mat>(3,3); /* Ap=(ap1,ap2,ap3) */
		mat T = zeros<mat>(3,3); /* inv(A) */

		vec x_vec = zeros<vec>(3) ; /* x_vec=(x,y,z) */
		vec xi_vec = zeros<vec>(3) ; /* xi_vec=(xi,eta,zeta) */

		/* setup a1_vec,...,a3_vec, delta_vec and ap1_vec,...,ap3_vec */
		for ( unsigned int i = 0; i<3; ++i)
		{
			a1_vec(i) = a1[i] ;
			a2_vec(i) = a2[i] ;
			a3_vec(i) = a3[i] ;
		}

		delta_vec(0) = 0.0 ;
		delta_vec(1) = 0.0 ;
		delta_vec(2) = 2.0 * a3_vec(2) ;

		ap1_vec = a1_vec ;
		ap2_vec = a2_vec ;
		ap3_vec = delta_vec ;

		if ( log_flag ) { fprintf( logfile, "\t a1_vec = ( %15.10G, %15.10G, %15.10G)\n", a1_vec[0], a1_vec[1], a1_vec[2]) ; fflush (logfile) ; }
		if ( log_flag ) { fprintf( logfile, "\t a2_vec = ( %15.10G, %15.10G, %15.10G)\n", a2_vec[0], a2_vec[1], a2_vec[2]) ; fflush (logfile) ; }
		if ( log_flag ) { fprintf( logfile, "\t a3_vec = ( %15.10G, %15.10G, %15.10G)\n", a3_vec[0], a3_vec[1], a3_vec[2]) ; fflush (logfile) ; }
		if ( log_flag ) { fprintf( logfile, "\n") ; }
		if ( log_flag ) { fprintf( logfile, "\tap1_vec = ( %15.10G, %15.10G, %15.10G)\n", ap1_vec[0], ap1_vec[1], ap1_vec[2]) ; fflush (logfile) ; }
		if ( log_flag ) { fprintf( logfile, "\tap2_vec = ( %15.10G, %15.10G, %15.10G)\n", ap2_vec[0], ap2_vec[1], ap2_vec[2]) ; fflush (logfile) ; }
		if ( log_flag ) { fprintf( logfile, "\tap3_vec = ( %15.10G, %15.10G, %15.10G)\n", ap3_vec[0], ap3_vec[1], ap3_vec[2]) ; fflush (logfile) ; }

		/* setup matrix A=(a1,a2,a3) and Ap=(ap1,ap2,ap3) */
		A.col(0) = a1_vec ;
		A.col(1) = a2_vec ;
		A.col(2) = a3_vec ;

		Ap.col(0) = ap1_vec ;
		Ap.col(1) = ap2_vec ;
		Ap.col(2) = ap3_vec ;

		/* set up transformation matrix from old coordinates/basis to new ones, T = inv(Ap) */
		T = Ap.i(true) ;

		/**********************/
		/* update cif-entries */
		/**********************/

		/* update symmetry */
		strcpy( sdummy, "_symmetry_space_group_name_H-M") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 )
		{
			strcpy( sdummy, "Unknown") ;
			cfe[idummy].text_data[0] = (char *) realloc ( cfe[idummy].text_data[0], ( strlen(sdummy) + 1 ) * sizeof(char)) ;
			strcpy( cfe[idummy].text_data[0], sdummy);
		}

		strcpy( sdummy, "_symmetry_space_group_name_Hall") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 )
		{
			strcpy( sdummy, "Unknown") ;
			cfe[idummy].text_data[0] = (char *) realloc ( cfe[idummy].text_data[0], ( strlen(sdummy) + 1 ) * sizeof(char)) ;
			strcpy( cfe[idummy].text_data[0], sdummy);
		}

		strcpy( sdummy, "_symmetry_equiv_pos_as_xyz") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) 
		{
			for ( unsigned int i=0; i<cfe[idummy].data_length; ++i)
			{
				free( cfe[idummy].text_data[i] ) ;
			}
			free( cfe[idummy].text_data ) ;

			cfe[idummy].data_length = 1 ;
			cfe[idummy].text_data = (char **) calloc( 1, sizeof(char *)) ;
			strcpy( sdummy, "+x,+y,+z") ;
			cfe[idummy].text_data[0] = (char *) calloc( strlen(sdummy)+1, sizeof(char)) ;
			strcpy( cfe[idummy].text_data[0], sdummy);
		}

		/* cell, update c and alpha, beta */
		strcpy( sdummy, "_cell_length_c") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { cfe[idummy].numeric_data[0] = 10.0 * norm( ap3_vec, 2) ; }

		strcpy( sdummy, "_cell_angle_alpha") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { cfe[idummy].numeric_data[0] = 90.0 ; }

		strcpy( sdummy, "_cell_angle_beta") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { cfe[idummy].numeric_data[0] = 90.0 ; }

		/* update atom with coordinates */
		atom = (atom_entries *) realloc (atom, 2*noa * sizeof(atom_entries)) ;
		for ( unsigned int i = 0; i<noa; ++i)
		{
			/* first unit cell remains the same */
			x_vec(0) = atom[i].coord[0] ;
			x_vec(1) = atom[i].coord[1] ;
			x_vec(2) = atom[i].coord[2] ;

			/* define position vector with old basis and coordinates x, y, z */
			x_vec = A * x_vec ;
			/* transform to new basis coordinates xi, eta, zeta */
			xi_vec = T * x_vec ;

			atom[i].coord[0] = xi_vec(0) ;
			atom[i].coord[1] = xi_vec(1) ;
			atom[i].coord[2] = xi_vec(2) ;

			/* (001)-mirrored coordinates for second unit cell, use x_vec from above,
			   define position vector with old basis and coordinates x, y, z */

			/* x_vec hasn't changed, use it again */
			x_vec(2) *= -1 ;
			x_vec += delta_vec ;

			/* transform to new basis coordinates xi, eta, zeta */
			xi_vec = T * x_vec ;

			atom[i+noa].coord[0] = xi_vec(0) ;
			atom[i+noa].coord[1] = xi_vec(1) ;
			atom[i+noa].coord[2] = xi_vec(2) ;
		}
		/* update noa and atom labels and types */
		noa *= 2 ;
		derive_atom_labels_and_types_from_cif() ;
	}
	#endif

	#if (ARMADILLO)
	/* basis transformation, test scenario with a doubled unit cell along a3-direction for triglyceride crystals */
	void bt_double()
	{
		char sdummy[1024] ;
		int idummy ;

		vec a1_vec = zeros<vec>(3) ;
		vec a2_vec = zeros<vec>(3) ;
		vec a3_vec = zeros<vec>(3) ;

		vec ap1_vec = zeros<vec>(3) ;
		vec ap2_vec = zeros<vec>(3) ;
		vec ap3_vec = zeros<vec>(3) ;

		mat A = zeros<mat>(3,3) ; /* A=(a1,a2,a3) */
		mat Ap = zeros<mat>(3,3); /* Ap=(ap1,ap2,ap3) */
		mat T = zeros<mat>(3,3); /* inv(A) */

		vec x_vec = zeros<vec>(3) ; /* x_vec=(x,y,z) */
		vec xi_vec = zeros<vec>(3) ; /* xi_vec=(xi,eta,zeta) */

		/* setup a1_vec,...,a3_vec, delta_vec and ap1_vec,...,ap3_vec */
		for ( unsigned int i = 0; i<3; ++i)
		{
			a1_vec(i) = a1[i] ;
			a2_vec(i) = a2[i] ;
			a3_vec(i) = a3[i] ;
		}

		ap1_vec = a1_vec ;
		ap2_vec = a2_vec ;
		ap3_vec = 2.0 * a3_vec ;

		if ( log_flag ) { fprintf( logfile, "\t a1_vec = ( %15.10G, %15.10G, %15.10G)\n", a1_vec[0], a1_vec[1], a1_vec[2]) ; fflush (logfile) ; }
		if ( log_flag ) { fprintf( logfile, "\t a2_vec = ( %15.10G, %15.10G, %15.10G)\n", a2_vec[0], a2_vec[1], a2_vec[2]) ; fflush (logfile) ; }
		if ( log_flag ) { fprintf( logfile, "\t a3_vec = ( %15.10G, %15.10G, %15.10G)\n", a3_vec[0], a3_vec[1], a3_vec[2]) ; fflush (logfile) ; }
		if ( log_flag ) { fprintf( logfile, "\n") ; }
		if ( log_flag ) { fprintf( logfile, "\tap1_vec = ( %15.10G, %15.10G, %15.10G)\n", ap1_vec[0], ap1_vec[1], ap1_vec[2]) ; fflush (logfile) ; }
		if ( log_flag ) { fprintf( logfile, "\tap2_vec = ( %15.10G, %15.10G, %15.10G)\n", ap2_vec[0], ap2_vec[1], ap2_vec[2]) ; fflush (logfile) ; }
		if ( log_flag ) { fprintf( logfile, "\tap3_vec = ( %15.10G, %15.10G, %15.10G)\n", ap3_vec[0], ap3_vec[1], ap3_vec[2]) ; fflush (logfile) ; }

		/* setup matrix A=(a1,a2,a3) and Ap=(ap1,ap2,ap3) */
		A.col(0) = a1_vec ;
		A.col(1) = a2_vec ;
		A.col(2) = a3_vec ;

		Ap.col(0) = ap1_vec ;
		Ap.col(1) = ap2_vec ;
		Ap.col(2) = ap3_vec ;

		/* set up transformation matrix from old coordinates/basis to new ones, T = inv(Ap) */
		T = Ap.i(true) ;


		/**********************/
		/* update cif-entries */
		/**********************/

		/* update symmetry */
		strcpy( sdummy, "_symmetry_space_group_name_H-M") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 )
		{
			strcpy( sdummy, "Unknown") ;
			cfe[idummy].text_data[0] = (char *) realloc ( cfe[idummy].text_data[0], ( strlen(sdummy) + 1 ) * sizeof(char)) ;
			strcpy( cfe[idummy].text_data[0], sdummy);
		}

		strcpy( sdummy, "_symmetry_space_group_name_Hall") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 )
		{
			strcpy( sdummy, "Unknown") ;
			cfe[idummy].text_data[0] = (char *) realloc ( cfe[idummy].text_data[0], ( strlen(sdummy) + 1 ) * sizeof(char)) ;
			strcpy( cfe[idummy].text_data[0], sdummy);
		}

		strcpy( sdummy, "_symmetry_equiv_pos_as_xyz") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) 
		{
			for ( unsigned int i=0; i<cfe[idummy].data_length; ++i)
			{
				free( cfe[idummy].text_data[i] ) ;
			}
			free( cfe[idummy].text_data ) ;

			cfe[idummy].data_length = 1 ;
			cfe[idummy].text_data = (char **) calloc( 1, sizeof(char *)) ;
			strcpy( sdummy, "+x,+y,+z") ;
			cfe[idummy].text_data[0] = (char *) calloc( strlen(sdummy)+1, sizeof(char)) ;
			strcpy( cfe[idummy].text_data[0], sdummy);
		}

		/* cell, update c */
		strcpy( sdummy, "_cell_length_c") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { cfe[idummy].numeric_data[0] = 10.0 * norm( ap3_vec, 2) ; }

		/* update atom with coordinates */
		atom = (atom_entries *) realloc (atom, 2*noa * sizeof(atom_entries)) ;
		for ( unsigned int i = 0; i<noa; ++i)
		{
			/* first unit cell remains the same */
			x_vec(0) = atom[i].coord[0] ;
			x_vec(1) = atom[i].coord[1] ;
			x_vec(2) = atom[i].coord[2] ;

			/* define position vector with old basis and coordinates x, y, z */
			x_vec = A * x_vec ;
			/* transform to new basis coordinates xi, eta, zeta */
			xi_vec = T * x_vec ;

			atom[i].coord[0] = xi_vec(0) ;
			atom[i].coord[1] = xi_vec(1) ;
			atom[i].coord[2] = xi_vec(2) ;

			/* second unit cell is the same as the first one, but shifted by a3, use x_vec from above */
			/* define position vector with old basis and coordinates x, y, z */

			/* x_vec' = A * x_vec hasn't changed, use it again */
			/* x_vec = x_vec' + a3_vec */
			x_vec += a3_vec ;

			/* transform to new basis coordinates xi, eta, zeta */
			xi_vec = T * x_vec ;

			atom[i+noa].coord[0] = xi_vec(0) ;
			atom[i+noa].coord[1] = xi_vec(1) ;
			atom[i+noa].coord[2] = xi_vec(2) ;
		}
		/* update noa and atom labels and types */
		noa *= 2 ;
		derive_atom_labels_and_types_from_cif() ;
	}
	#endif


	/* Not yet implemented

	   The Multiplicity for each atom in a cif-file seems to 1.0 -> set 1.0
	
	   The 7 parameters for thermal displacement ( -> Debye Waller factor) Uiso (1x) and Uani (6x) are set to 0.0
	   However, they might be read from the cif file, if they are provided there
	   (_atom_site_U_iso_or_equiv for isotropic and data_atom_site_aniso_U_ij for anisotropic, or
	    _atom_site_B_iso_or_equiv for isotropic and data_atom_site_aniso_B_ij for anisotropic, ... )
	*/
	void derive_Uiso_Uani_and_Mult_from_cif()
	{
		for ( unsigned int i=0; i<noa; ++i)
		{
			atom[i].mult = 1.0 ;
			for ( unsigned int j=0; j<7; ++j) { atom[i].temp[j] = 0.0 ; }
		}
	}


	/* updates the atom struct labels and types with the entries given in the cif-file.
	   intended for using when noa is a multiple of noa(cif) what is checked.
	   if so then for( i=0:noa/noa(cif-1) for( j=0:noa(cif-1) atom(i*noa(cif)+j).label='cif-label(j)'+'ext_i' ) )
	   types will be derived from initial letters of label or from cif-file (if provided)

	   D-isotopes are assigned by type 2H
	*/
	void derive_atom_labels_and_types_from_cif()
	{
		if ( log_flag ) { fprintf( logfile, "\tderive_atom_labels_and_types_from_cif()\n") ; }

		string sdummy, sdummy2, sdummy3 ;
		ostringstream ssdummy ;
		
		int i_label = find_index_in_cfe("_atom_site_label") ;
		/* look also for type in cif-file, if it does not exist only warn and derive type from label */
		int i_type = find_index_in_cfe("_atom_site_type_symbol", 0) ;
		if ( i_type < 0 )
		{ 
			if ( log_flag )
			{ 
				fprintf( logfile, "\t\tCould not find the entry '_atom_site_type_symbol' in cif-file.\n") ;
				fprintf( logfile, "\t\tAtom types will be derived from atom labels.\n") ;
			}
		}
		unsigned int noa_cif = cfe[i_label].data_length ;

		if ( noa % noa_cif != 0 ) { XNDIFF_ERROR(65) ; }

		unsigned int noa_ratio = noa / noa_cif ;

		for ( unsigned int i=0; i<noa_cif; ++i)
		{
			/* read cif-label */
			sdummy = cfe[i_label].text_data[i] ;
			/* type, use _atom_site_type_symbol if included in cif-file,
			   otherwise derive type from label 
			*/
			if ( i_type > 0 )
			{
				sdummy2 = cfe[i_type].text_data[i] ;
			}
			else
			{
				sdummy2 = sdummy.at(0) ;
				/* symbol has at least one letter, if more read them, ignore numbers */
				for ( unsigned int j=1; j<sdummy.length(); ++j)
				{
					if ( !isdigit(sdummy[j]) ) { sdummy2 += sdummy.at(j) ; }
					else { break ; }
				}
			}
			/* furthermore replace deuterium symbol D with 2H */
			if ( sdummy2.compare("D") == 0 ) { sdummy2 = "2H" ; }

			/* update all labels and types for the bunch of i+j*noa_cif where j runs over all noa_ratio */
			for ( unsigned int j=0; j<noa_ratio; ++j)
			{
				/* extend label */
				ssdummy << j ;
				sdummy3 = sdummy + "_" + ssdummy.str() ;
				ssdummy.str("");
				ssdummy.clear();

				/* update label and type */
				strcpy( atom[i+j*noa_cif].label, sdummy3.c_str() ) ;
				strcpy( atom[i+j*noa_cif].type, sdummy2.c_str() ) ;
			}
		}
		if ( log_flag ) { fprintf( logfile, "\tdone\n") ; }
	}


	/* similar as for atom labels and types, update bond_site12, bond_symmetry12, bond_distance */
	void derive_connectivity_table_from_cif()
	{
		/* number of atoms (noa) and bonds (nob) in cif-file */
		int i_label = find_index_in_cfe("_atom_site_label") ;
		unsigned int noa_cif = cfe[i_label].data_length ;
		unsigned int nob_cif = bond_site12.size() ;

		if ( noa % noa_cif != 0 ) { XNDIFF_ERROR(65) ; }

		unsigned int noa_ratio = noa / noa_cif ;
		vector<unsigned int> uidummy ;

		/* update only everything above noa_cif */
		for ( unsigned int j=1; j<noa_ratio; ++j)
		{
			for ( unsigned int i=0; i<nob_cif; ++i)
			{
				uidummy.assign( 2, j*noa_cif) ;
				uidummy[0] += bond_site12[i][0] ;
				uidummy[1] += bond_site12[i][1] ;
				bond_site12.push_back( uidummy ) ;
				bond_symmetry12.push_back( bond_symmetry12[i] ) ;

				bond_distance.push_back( bond_distance[i] ) ;
			}
		}
	}


	/* check if ifo and ofo folders exist, if not create them */
	void check_ifo_ofo()
	{
		/* mkdir returns 0 on success, -1 if an error occurred.
		   http://home.fhtw-berlin.de/~junghans/cref/MAN/mkdir.htm
		   http://stackoverflow.com/questions/675039/how-can-i-create-directory-tree-in-c-linux 
		*/
		if ( !directory_exists(ifo) )
		{
			if ( mkdir( ifo, 0777) != 0 ) { XNDIFF_ERROR(43) ; } ;
		}
		if ( !directory_exists(ofo) )
		{
			if ( mkdir( ofo, 0777) != 0 ) { XNDIFF_ERROR(44) ; } ;
		}
	}


	/* returns string about remaining time as computed from the time 
	   elapsed up to now, et [s], and the progress achieved within this time, pr [%]. */
	void remtime (long int et, double pr, char* str)
	{
		long int rt ;

		if ( pr < 1.0e-3 )  /* at least 0.1 % for a forecast */
			strcpy( str , "?" ) ;
		else
		{
			rt = (long int)( ( 100.0 / pr - 1.0 ) * ( (double) et ) ) ;
			if ( rt < 60 )
				sprintf ( str, "%02ld (ss)", rt) ;
			else
			{
				if ( rt < 3600 )
					sprintf ( str, "%02ld:%02ld (mm:ss)", ( rt - (rt % 60) ) / 60 , rt % 60) ;
				else
				{
					if ( rt < 86400 )
						sprintf ( str, "%02ld:%02ld:%02ld (hh:mm:ss)", ( rt - (rt % 3600) ) / 3600 , ( (rt % 3600) - ( (rt % 3600) % 60) ) / 60 , rt % 60) ;
					else 
						sprintf ( str, "%02ld:%02ld:%02ld:%02ld (dd:hh:mm:ss)", ( rt - (rt % 86400) ) / 86400, ( (rt % 86400) - ( (rt % 86400) % 3600) ) / 3600 , ( ( (rt % 86400) % 3600) - ( ( (rt % 86400) % 3600) % 60) ) / 60 , rt % 60) ;
					/* not more than 99 days */
				}
			}
		}
	}


	/* Neutron Scattering Lengths and cross sections into variable n_atom_sc */
	/* Neutron News, Vol. 3, No. 3, 1992, pp. 29-37 */
	/* http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/DeBe_NeutronNews.dat */
	void read_neutron_scattering_lengths_and_cross_sections()
	{
		FILE *inpf ; 
		const int numsigns = 1024 ;
		char sdummy[numsigns] ;
		char data[numsigns] ;
		char *strptr1, *strptr2, *strptr3;
		int i;
		int max_entries=9; /*#N 8 or 9 */

		/* Open the neutron scattering data file */ 
		strcpy( sdummy, ifo) ;
		strcat( sdummy, "DeBe_NeutronNews.dat") ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(31) ; }

		num_n_atom_sc_entries=0;
		while (fgets (sdummy, numsigns, inpf) != NULL)
		{
			strptr1 = sdummy ;
			while (*strptr1 == ' ') 
				++strptr1 ;

			if (*strptr1 == '\n')
				continue ;

			if (*strptr1 == '#')
			{
				++strptr1;
				if (*strptr1 == 'S')
				{
					/* new entry of the form '#S  1    H\t\n' */
					if (num_n_atom_sc_entries==0)
						n_atom_sc=(nesclecrsec *) calloc( 1, sizeof(nesclecrsec) ) ;
					else
						n_atom_sc=(nesclecrsec *) realloc ( n_atom_sc, (num_n_atom_sc_entries+1) * sizeof(nesclecrsec) ) ;

					++strptr1;
					/* skip the number */
					while (*strptr1 == ' ') 
						++strptr1 ;

					strptr2=strchr (strptr1, ' ');
					strptr1=strptr2;
					while (*strptr1 == ' ') 
						++strptr1 ;

					/* read nuclei name */
					if ((strptr2=strchr (strptr1, '\r'))==NULL)
						strptr2=strchr (strptr1, '\n');

					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;
					delallspc(data);

					n_atom_sc [num_n_atom_sc_entries].nuclei_name=(char *) calloc( strlen(data)+1, sizeof(char));
					strcpy(n_atom_sc[num_n_atom_sc_entries].nuclei_name, data);
				}
				else if (*strptr1 == 'N')
				{
					/* read number of entries, either 8 or 9 */
					/* #N 9\r\n */
					++strptr1;
					while (*strptr1 == ' ') 
						++strptr1 ;

					/* read nuclei name */
					if ((strptr2=strchr (strptr1, '\r'))==NULL)
						strptr2=strchr (strptr1, '\n');

					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;
					delallspc(data);

					max_entries=strtol(data,NULL,10);

				}
				else if (*strptr1 == 'L')
				{
					/* read data in the next line after #L */
					if ( fgets( sdummy, numsigns, inpf) == NULL)
						break;

					strptr1 = sdummy ;
					
					int num_entries=0;
					char delimiter='\t';

					if (max_entries==9)
						delimiter=' ';


					while ( (strchr (strptr1, '\r') != NULL) || (strchr (strptr1, delimiter) != NULL) )
					{
						if ((strptr2=strchr (strptr1, delimiter))==NULL)
							strptr2=strchr (strptr1, '\r');

						++num_entries;

						i = (int)(strptr2-strptr1) ;
						memmove (data, strptr1, i) ;
						data[i] = 0 ;

						/* read natural abundance or half life */
						if (num_entries==1)
						{
							if (!strcmp(data,"1"))
							{
								n_atom_sc[num_n_atom_sc_entries].nat_abd = 100.0 ;
							}
							else if (data[0]=='(')
							{
								/* strptr1 still points to position of first character of the current entry i.e. '(' */
								++strptr1;
								strptr3=strptr1;
								
								while ( (*strptr3!='a') && (*strptr3!=')') )
									++strptr3;

								i = (int)( strptr3-strptr1 ) ;
								memmove (data, strptr1, i) ;
								data[i] = 0 ;

								n_atom_sc[num_n_atom_sc_entries].nat_abd = strtod(data,NULL) ;
							}
							else
								n_atom_sc[num_n_atom_sc_entries].nat_abd = strtod(data,NULL) ;
						}
						else
						{
							/* skip (+/-) and '>' or '<' in front and e.s.d. values (...) at the end of an entry */
							if ( !strncmp(strptr1, "(+/-)", 5) )
								strptr1 += 5;

							if ( (*strptr1=='<') || (*strptr1=='>') )
								++strptr1;

							if ((strptr3=strchr (strptr1, '(')) == NULL)
								strptr3=strptr2;


							i = (int)( strptr3-strptr1 ) ;
							memmove (data, strptr1, i) ;
							data[i] = 0 ;

							/* there's no Im(b_inc) for N=8, i.e. atoms */
							/* for elements (i.e. not isotopes) there are only N=8 entries 
							   Re(b_coh), Im(b_coh) typically 0, Re(b_inc)=9999 (not defined for elements)
							   and Im(b_inc) is missing, then coh_cs, incoh_cs, scat_cs, abs_cs  */
							if ( max_entries==8 && num_entries>4 )
								++num_entries;

							switch (num_entries)
							{
								case 2:
									n_atom_sc[num_n_atom_sc_entries].coh_sc_length = dcmplx( strtod(data,NULL), 0.0);
									break;

								case 3:
									n_atom_sc[num_n_atom_sc_entries].coh_sc_length += dcmplx( 0.0, strtod(data,NULL));
									break;
								case 4:
									n_atom_sc[num_n_atom_sc_entries].incoh_sc_length = dcmplx( strtod(data,NULL), 0.0);
									break;
								case 5:
									n_atom_sc[num_n_atom_sc_entries].incoh_sc_length += dcmplx( 0.0, strtod(data,NULL));
									break;
								case 6:
									n_atom_sc[num_n_atom_sc_entries].coh_sc_cross = strtod(data,NULL);
									break;
								case 7:
									n_atom_sc[num_n_atom_sc_entries].incoh_sc_cross = strtod(data,NULL);
									break;
								case 8:
									n_atom_sc[num_n_atom_sc_entries].total_sc_cross = strtod(data,NULL);
									break;
								case 9:
									n_atom_sc[num_n_atom_sc_entries].abs_cross = strtod(data,NULL);
									break;
								default:
									break;
							}

							/* bring num_entries back to the right value */
							if ( max_entries==8 && num_entries>(4+1) )
								--num_entries;
						
						}

						if (num_entries==max_entries)	
							break;

						/* update the pointer strptr1 on the last position */
						strptr1=strptr2+1;
						while (*strptr1 == ' ')
							++strptr1 ;
					}

					if (num_entries!=max_entries)
					{
						/* print a warning message and continue*/
						if ( log_flag ) { fprintf( logfile, "\tWarning. Found only %d data entries for nuclei %s. Expected %d. Continuing.\n", num_entries, n_atom_sc[num_n_atom_sc_entries].nuclei_name, max_entries) ; }
					}

					/* now increase number of entries */
					++num_n_atom_sc_entries;
				}
				else
					continue ;
			}
		}
		/* close the neutron scattering data file */
		fclose(inpf);
	}


	/* clear n_atom_sc */
	void free_neutron_scattering_lengths_and_cross_sections()
	{
		for (int i=0; i<num_n_atom_sc_entries; ++i)
		{
			free(n_atom_sc[i].nuclei_name);
		}
		free(n_atom_sc);
	}


	/* returns the index i (0<=i<=num_n_atom_sc_entries-1) of an entry in the  */
	/* Neutron Scattering Length and Cross Sections array n_atom_sc by comparing with its element/isotope names */
	/* if an entry is not found it returns -1 */
	int find_index_in_neutron_scattering_lengths_and_cross_sections(const char* entry)
	{
		int i=-1;
		
		for (int j=0; j<num_n_atom_sc_entries; ++j)
		{
			if (!strcmp(entry,n_atom_sc[j].nuclei_name))
			{
				i=j;
				break;
			}
		}

		if ( i < 0 )
		{
			if ( log_flag )
				fprintf( logfile, "\tWarning: Could not find the element/isotope %s in the neutron scattering database.\n", entry) ;
		}

		return i;
	}


	dcmplx get_coherent_scattering_length(const char* aatom)
	{
		/* check on charge labels and remove them, only the isotopic composition is relevant for neutron scattering */
		/* e.g. Ca2+ -> Ca or 2H1- -> 2H */
		int i;
		int strlenaatom = strlen(aatom);
		char* entry = (char *) calloc( strlenaatom+1, sizeof(char) ) ;

		strcpy(entry, aatom);
		if ( entry[strlenaatom-1]=='+' || entry[strlenaatom-1]=='-' )
		{
			entry[strlenaatom-1] = 0 ;
			i = strlenaatom-2 ;
			while ( isdigit(entry[i]) )
			{
				entry[i] = 0 ;
				--i;
			}
		}
	
		i = find_index_in_neutron_scattering_lengths_and_cross_sections(entry) ;
		free(entry);

		if ( i < 0 )
		{
			fprintf( logfile, "\tWarning: The coherent scattering length for the missing element/isotope %s has been set to 0.0 + 0.0 i.\n", entry) ;
			return dcmplx(0.0,0.0);
		}
		else
		{
			return n_atom_sc[i].coh_sc_length ;
		}
	}


	double get_incoherent_scattering_crosssec(const char* aatom)
	{
		/* check on charge labels and remove them, only the isotopic composition is relevant for neutron scattering */
		/* e.g. Ca2+ -> Ca or 2H1- -> 2H */
		int i;
		int strlenaatom = strlen(aatom);
		char* entry = (char *) calloc( strlenaatom+1, sizeof(char) ) ;

		strcpy(entry, aatom);
		if ( entry[strlenaatom-1]=='+' || entry[strlenaatom-1]=='-' )
		{
			entry[strlenaatom-1] = 0 ;
			i = strlenaatom-2 ;
			while ( isdigit(entry[i]) ) 
			{
				entry[i] = 0 ;
				--i;
			}
		}
	
		i = find_index_in_neutron_scattering_lengths_and_cross_sections(entry) ;
		free(entry);

		if ( i < 0 )
		{
			fprintf( logfile, "\tWarning: The incoherent scattering cross section for the missing element/isotope %s has been set to 0.0.\n", entry) ;
			return 0.0;
		}
		else
		{
			return n_atom_sc[i].incoh_sc_cross;
		}
	}


	/* computes the X-ray electron density in the crystal's unit cell in [electrons/nm^3] */
	void compute_electron_density()
	{
		rho_cry = 0.0 ;
		for (unsigned int l=0; l<noa; ++l)
			rho_cry += get_atomic_scattering_factor ( atom[l].type, 0.0 ) ; /* [electrons] */

		rho_cry /= V_uc ; rho_cry *= 1.0 ; /* V_uc [nm^3], to get [electrons/nm^3] */
	}


	/* computes the incoherent scattering cross section density (ICSD) in the crystal's unit cell in [1/cm] */
	void compute_incoherent_cross_section_density()
	{
		double ddummy ;

		icsd_cry = 0.0 ;
		/* spin-incoherent contribution */
		ddummy = 0.0 ;
		for ( unsigned int l=0; l<noa; ++l)
			ddummy += atom[l].incoh_sc_cross / (4.0 * M_PI ) ; /* [barn=10^(-24) cm^2] */
		
		icsd_cry += ddummy ;

		/* composition-incoherent contribution */
		ddummy = 0.0 ;
		for ( unsigned int l=0; l<noa; ++l)
			ddummy += pow( real( atom[l].coh_sc_length ), 2.0 ) ; /* [fm^2] */

		ddummy *= 1.0e-2 ; /* (b_coh)^2 [fm^2], to get [barn=10^(-24) cm^2] */ 
		icsd_cry += ddummy ;

		ddummy = 0.0 ;
		for ( unsigned int l=0; l<noa; ++l)
			ddummy += real( atom[l].coh_sc_length ) ; /* [fm] */

		ddummy = pow( ddummy, 2.0 ) / ( (double) noa ) ;
		ddummy *= 1.0e-2 ; /* (b_coh)^2 [fm^2], to get [barn=10^(-24) cm^2] */
		icsd_cry -= ddummy ; /* subtract ! */

		icsd_cry /= V_uc ; icsd_cry *= 1.0e-3 ; /* V_uc [nm^3], to get [cm^(-1)] for ICSD */
	}


	/* computes the real part of the scattering legth density (SLD) in the crystal's unit cell in [10^(-6) 1/A^2 = 10^10 1/cm^2] */
	void compute_scattering_length_density()
	{
		sld_cry = 0.0 ;
		for ( unsigned  int l=0; l<noa; ++l)
			sld_cry += ( real(atom[l].coh_sc_length) ) ; /* [fm] */

		sld_cry /= V_uc ; sld_cry *= 1.0e-2 ; /* V_uc [nm^3], to get [10^(-6) 1/A^2] for SLD */
	}


	/* Atomic Weights and Isotopic Compositions for All Elements 		    */
	/* http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii&isotype=some */
	/* http://en.wikipedia.org/wiki/Relative_atomic_mass 			    */
	/* for normal atoms in a compound take the Standard Atomic Weight value     */
	/* for specific isotopes take the Relative Atomic Mass directly 	    */ 
	/*									    */
	/*              Relative           Isotopic      Standard                   */
	/* Isotope      Atomic Mass        Composition   Atomic Weight    Notes     */
	/* _________________________________________________________________________*/
	/* 1   H   1    1.0078250321(4 99.9885(70)   1.00794(7)       g,m,r,c,w     */
	/*     D   2    2.0141017780(4)    0.0115(70)  				    */
	/*     T   3    3.0160492675(11)               				    */
	/*     H   4    4.02783(12)                    				    */
	/*         5    5.03954(102)                   				    */
	/*         6    6.04494(28)    						    */
	/*          1111111111222222222233333333334444444444555555555566666666667777*/ 
	/* 1234567890123456789012345678901234567890123456789012345678901234567890123*/
	/*									    */
	/* for elements: 							    */
	/* - read element symbol from pos 5-7 after a new section begins    	    */
	/* - read Standard Atomic Weight from pos. 47 or 48 on up to (...)	    */
	/*   (for elements in []-brackets e.g. [209] for Po read the value	    */
	/*    within the []-brackets and give awarning to the logfile)		    */
	/*									    */
	/* for isotopes: 							    */
	/* - use element symbol from pos 5-7 from the element entry above 	    */
	/* - read isotopic number from pos. 9-11 				    */
	/* - read Relative Atomic Mass from pos. 14 on up to (...)		    */
	void read_standard_atomic_weight_and_relative_atomic_mass()
	{
		FILE *inpf ; 
		const int numsigns=1024;
		char sdummy[numsigns] ;
		char data[numsigns] ;
		char *strptr1, *strptr2 ;
		int i, strlen_current_sym ;
		char current_sym[4]; /* "Uuu\0" */
		bool new_entry=false;

		/* Open the neutron scattering data file */
		strcpy( sdummy, ifo) ;
		strcat( sdummy, "atomweights.dat") ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(32) ; }

		num_atom_mass_entries=0;
		while (fgets (sdummy, numsigns, inpf) != NULL)
		{
			strptr1 = sdummy ;
			if ( strlen(sdummy)<72 )
			{ 
				for (int j=strlen(sdummy)+1; j<73; ++j)
					sdummy[j]=' ';
			}	

			/* check for new entries after lines with 73 underscore signs */
			if (*strptr1 == '_')
			{
				new_entry = true ;
				for (int j=0; j<72; ++j)
				{
					if ( *( ++strptr1 ) != '_')
					{
						new_entry = false ;
						break;
					}
				}
				continue ;
			}

			if ( new_entry )
			{
				/* read data for elements from pos. 47 or 48 on */
				/* the first case is indicated by a Standard Atomic Weight not written in []-brackets */

				if (num_atom_mass_entries==0)
				{
					if ( isdigit(sdummy[46]) || sdummy[46]=='[' )
						atom_mass = (StAtWtandRelAtMa *) calloc( 2, sizeof(StAtWtandRelAtMa) );
					else
						atom_mass = (StAtWtandRelAtMa *) calloc( 1, sizeof(StAtWtandRelAtMa) );			
				}
				else
				{
					if ( isdigit(sdummy[46]) || sdummy[46]=='[' )
						atom_mass = (StAtWtandRelAtMa *) realloc ( atom_mass, (num_atom_mass_entries+2) * sizeof(StAtWtandRelAtMa) );
					else
						atom_mass = (StAtWtandRelAtMa *) realloc ( atom_mass, (num_atom_mass_entries+1) * sizeof(StAtWtandRelAtMa) );
				}

				if ( isdigit(sdummy[46]) || sdummy[46]=='[' )
				{
					/* read current element symbol from pos 5 on */
					strptr1 = strptr1 + 4 ;
					strptr2 = strchr ( strptr1, ' ') ;
					strlen_current_sym = i = (int)(strptr2-strptr1) ;
					memmove (current_sym, strptr1, i) ;
					current_sym[i] = 0 ;

					strptr1 = sdummy ;

					/* elements */

					atom_mass[num_atom_mass_entries].nuclei_name = (char *) calloc( strlen_current_sym+1, sizeof(char));
					strcpy(atom_mass[num_atom_mass_entries].nuclei_name, current_sym);

					/* read current element Standard Atomic Weight from pos 47 ( or 48 for []-brackets ) on */
					strptr2 = strptr1 = strptr1 + 46 ;
					if ( sdummy[46]=='[' )
					{
						++strptr1;
						++strptr2;
					}
					while ( isdigit(*strptr2) || *strptr2 == '.' ) 
						++strptr2;

					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;

					atom_mass[num_atom_mass_entries].mass = strtod( data, NULL) ;

					if ( sdummy[46]=='[' )
					{
						if ( log_flag )
							fprintf( logfile, "\tWarning: No precise Standard Atomic Weight for element %s is defined. Using the given value %d.\n", atom_mass[num_atom_mass_entries].nuclei_name, (int)(atom_mass[num_atom_mass_entries].mass) ) ;	
					}
		
					++num_atom_mass_entries;

					strptr1 = sdummy ;
				}

				/* isotopes */

				/* read current isotope number from pos 9 on */
				strptr1 = strptr1 + 8 ;

				if ( isdigit(*strptr1) )
				{
					strptr2 = strptr1 ;
					while ( isdigit(*strptr2) ) 
						++strptr2;
	
					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;
	
					atom_mass[num_atom_mass_entries].nuclei_name = (char *) calloc( strlen_current_sym + i + 1, sizeof(char) ) ;
					/* copy isotope number */
					for (int j=0; j<i; ++j)
						atom_mass[num_atom_mass_entries].nuclei_name[j] = data[j] ;
	
					/* copy element symbol + terminating 0 */
					for (int j=0; j<strlen_current_sym + 1; ++j)
						atom_mass[num_atom_mass_entries].nuclei_name[i+j] = current_sym[j] ;
	
					/* read isotope Relative Atomic Mass from pos 14 on */	
					strptr2 = strptr1 = strptr1 + 5 ;
					while ( isdigit(*strptr2) || *strptr2 == '.' ) 
						++strptr2;

					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;

					atom_mass[num_atom_mass_entries].mass = strtod( data, NULL) ;

					++num_atom_mass_entries;
					/* continue; */
				}
				else
				{
					/* if there are no new isotopic numbers finish this entry */
					new_entry = false ;
					/* continue; */
				}
			}
		}
		/* close the neutron scattering data file */
		fclose(inpf);
	}


	/* clears atom_mass */
	void free_standard_atomic_weight_and_relative_atomic_mass()
	{
		for (int i=0; i<num_atom_mass_entries; ++i)
		{
			free(atom_mass[i].nuclei_name);
		}
		free(atom_mass);
	}


	/* returns the index i (0<=i<=num_atom_mass_entries-1) of an entry in the  */
	/* Standard Atomic Weight and Relative Atom Mass array atom_mass by comparing with its element/isotope names */
	/* if an entry is not found it returns -1 */
	int find_index_in_standard_atomic_weight_and_relative_atomic_mass(const char* entry)
	{
		int i=-1;

		for (int j=0; j<num_atom_mass_entries; ++j)
		{
			if ( !strcmp(entry,atom_mass[j].nuclei_name) )
			{
				i = j ;
				break ;
			}
		}

		return i;
	}

	/* tries for each element/isotope entry in the n_atom_sc database 
	   to find the corresponding element/isotope mass in the atom_mass database and to include it */
	void add_atom_mass_to_neutron_scattering_database()
	{
		int num_fail=0;
		for (int i=0; i<num_n_atom_sc_entries; ++i)
		{
			int j;
			for (j=0; j<num_atom_mass_entries; ++j)
			{
				if ( !strcmp(atom_mass[j].nuclei_name,n_atom_sc[i].nuclei_name) )
					break;
			}
			/* print message to logfile if an entry was not found and set the mass of this entry to -1.0 */
			if (j==num_atom_mass_entries)
			{
				if ( log_flag )
					fprintf( logfile, "\tWarning: Could not find an element/isotope mass for %s in element/isotope mass database\n", n_atom_sc[i].nuclei_name) ;

				n_atom_sc[i].mass = -1.0 ;
				++num_fail;
			}
			else
			{
				n_atom_sc[i].mass = atom_mass[j].mass ;
			}
	
		}

		if (num_fail>0) { fprintf( logfile, "\tThe masses of the missing element(s)/isotope(s) have been set to -1.0\n") ; }

	}

	/* reads cif-dictionary (option +cif cif-file cif-dictionary) information into variable cde */
	void read_cif_dictionary()
	{

		FILE *inpf ; /* FILE pointer, for the cif-dictionary */
		const int numsigns = 1024 ; /* in principle 80/81 should be sufficient, because cif files are restricted to 80 characters in text-width */
		int i,j;
		char sdummy[numsigns] ;
		char data[numsigns] ;
		bool commentblock;
		bool categoryblock, dicversionblock;
		bool loopnameblock;
		char *strptr1, *strptr2 ;
		long int lineindex = 0;


		/* Open dictionary */
		strcpy( sdummy, ifo) ;
		strcat( sdummy, cif_dic) ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(21) ; }

		/* read line by line */
		commentblock=false;
		categoryblock=false;
		loopnameblock=false;
		dicversionblock=false;
		
		num_cif_dic_entries = 0 ;
		while ( fgets( sdummy, numsigns, inpf) != NULL)
		{
			++lineindex;

			/* only '#'' , ';', empty line, 'data_' allowed at the beginning, additional white space in front will be always ignored  */
			strptr1 = sdummy ;
			while (*strptr1 == ' ') 
				++strptr1 ;

			if (commentblock)
			{ 
				if(*strptr1 == ';')
					commentblock=false;
				
				continue;
			}

			if (*strptr1 == ';')
			{
				commentblock=true;
				continue ;
			}

			if (*strptr1 == '#')
				continue ;
			if (*strptr1 == '\n')
				continue ;

			/* CHECK for 'data' */
			/* http://www.cplusplus.com/reference/clibrary/cstring/strncmp/ */
			if ( !strncmp(strptr1, "data_", 5) )
			{
				/* new data_(category,dictionary-version) block */
				dicversionblock=false;
				categoryblock=false;

				/* ignore the category entries, that are ending with data_...._[] */
				/* find substring http://www.cplusplus.com/reference/clibrary/cstring/strstr/ */
				if ( strstr(strptr1, "_[]") != NULL)
				{
					categoryblock=true;
					continue;
				}

				/* also skip the information on this dictionary version */
				if ( strstr(strptr1, "on_this_dictionary") != NULL)
				{
					dicversionblock=true;
					continue;
				}

				if ( num_cif_dic_entries == 0 )
				{
					cde = (cif_dic_entries *) calloc( 1, sizeof(cif_dic_entries)) ;
				}
				else
				{
					cde = (cif_dic_entries *) realloc (cde, (num_cif_dic_entries + 1) * sizeof(cif_dic_entries)) ;
				}
				++num_cif_dic_entries;

				/* make default unit '1' */
				cde[num_cif_dic_entries-1].units = (char *) calloc( 1, sizeof(char)) ;
				cde[num_cif_dic_entries-1].units[0]='1';

				continue;
			}

			/* ignore category and dicversion data blocks */
			/* this will take as long as a new data_ entry is found above */
			if (categoryblock || dicversionblock )
				continue;
		
			/* if a data block is found read its entries */
			/* check for _name or loop_ _name, _category, _type, _units  */

			/* read the entries until the first white space after them */
			/* strptr2 will be overwritten from above, strptr1 same position as above */

			if ((strptr2=strchr (strptr1, ' ')) == NULL)
			{
				strptr2=strchr (strptr1, '\n');

				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;

				if (data[0]!='\'')
				{
					loopnameblock=false;
					continue;
				}
			}
			else
			{
				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;
	
				/* strptr1 points now to the first white space after entryname */
				strptr1=strptr2;

			}

			/* single _name field */
			if (!strcmp(data,"_name"))
			{
				loopnameblock=false;

				/* read value of the '_name' entry until the EOL into data */	
				strptr2=strchr (strptr1, '\n');
				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;

				/* expecting the form ' 'entryname'', delete the white spaces in data */
				/* (REMIND: entrynames itself are not allowed to contain whitespaces) */
				/* count the amount of characters between the apostrophes */
				/* (for security restrict j to numsigns) */
				/* allocate the names-field and write the entryname inside it */
				delallspc(data);
				j=1;
				while ( (data[j] != '\'') && (j < numsigns) ) 
					++j;
				
				cde[num_cif_dic_entries-1].num_names=1;
				cde[num_cif_dic_entries-1].names=(char **) calloc( 1, sizeof(char *)) ;
				cde[num_cif_dic_entries-1].names[0]=(char *) calloc( (j), sizeof(char)) ;

				for (i=1; i<j; ++i)
				{
					cde[num_cif_dic_entries-1].names[0][i-1]=data[i];
				}
				cde[num_cif_dic_entries-1].names[0][j-1]=0;
				continue;
			}

			/* _category field */
			if (!strcmp(data,"_category"))
			{
				loopnameblock = false ;

				/* read value of the '_category' entry until the EOL into data */
				strptr2=strchr (strptr1, '\n');
				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;
				/* REMIND: category names do not contain any white spaces itself by format specification */
				delallspc(data);

				j=0;
				while ( (data[j] != '\0') && (j < numsigns) ) 
					++j;
				
				cde[num_cif_dic_entries-1].category=(char *) calloc( (j), sizeof(char)) ;

				for (i=0; i<j; ++i)
				{
					cde[num_cif_dic_entries-1].category[i]=data[i];
				}
				continue;
			}

			/* _type field */
			if (!strcmp(data,"_type"))
			{
				loopnameblock = false ;

				cde[num_cif_dic_entries-1].isnumeric=false;

				/* read value of the '_name' entry until the EOL into data */
				strptr2=strchr (strptr1, '\n');
				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;
				delallspc(data);

				/* numb or char ? */
				if (!strcmp(data,"numb"))
					cde[num_cif_dic_entries-1].isnumeric=true;

				continue;
			}

			/* _units field */
			if (!strcmp(data,"_units"))
			{
				loopnameblock = false ;

				/* read value of the '_units' entry until the EOL into data */	
				strptr2=strchr (strptr1, '\n');
				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;
				delallspc(data);

				j=0;
				while ( (data[j] != '\0') && (j < numsigns) ) { ++j; }

				cde[num_cif_dic_entries-1].units=(char *) realloc ( cde[num_cif_dic_entries-1].units, (j) * sizeof(char)) ;

				for (i=0; i<j; ++i)
				{
					cde[num_cif_dic_entries-1].units[i]=data[i];
				}
				continue;
			}

			/* multiple 'loop_ _name' field */
			if (loopnameblock)
			{
				/* read first value of the '_name' entries until the EOL into data */	

				/* expecting the form ' 'entryname'', delete the white spaces in data */
				/* (REMIND: entrynames itself are not allowed to contain whitespaces) */
				/* count the amount of characters between the apostrophes */
				/* (for security restrict j to numsigns) */
				/* allocate the names-field and write the entryname inside it */
				delallspc(data);

				j=1;
				while ( (data[j] != '\'') && (j < numsigns) ) { ++j; }

				cde[num_cif_dic_entries-1].num_names += 1;
				cde[num_cif_dic_entries-1].names=(char **) realloc ( cde[num_cif_dic_entries-1].names , cde[num_cif_dic_entries-1].num_names * sizeof(char *)) ;
				cde[num_cif_dic_entries-1].names[cde[num_cif_dic_entries-1].num_names-1]=(char *) calloc( (j), sizeof(char)) ;

				for (i=1; i<j; ++i)
				{
					cde[num_cif_dic_entries-1].names[cde[num_cif_dic_entries-1].num_names-1][i-1]=data[i];
				}
				cde[num_cif_dic_entries-1].names[cde[num_cif_dic_entries-1].num_names-1][j-1]=0;
				continue;
			}

			/* multiple 'loop_ _name' field, but also 'loop_ _enumerates', ... */
			/* http://www.cplusplus.com/reference/clibrary/cstring/strncmp/ */
			if ( !strncmp(data, "loop_", 5) )
			{

				/* skip possible white spaces */
				while (*strptr1 == ' ')
					++strptr1 ;

				if ((strptr2=strchr (strptr1, ' ')) == NULL)
				{
					// fprintf(stdout, "Expecting control sequence after 'loop_ ...' in cif-dictionary in line %ld .\n", lineindex);
					continue;
				}

				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;

				/* strptr1 points now to the first white space after 'loop_ _entryname' */
				strptr1=strptr2;

				if (!strcmp(data,"_name"))
				{
					
					/* loopnameblock will true as long as _units, _type, _category, _name will be found */
					/* or the string starts not with an apostrophe ' (after ignoring the white space) */
					loopnameblock=true;

					/* read first value of the '_name' entries until the EOL into data */	
					strptr2=strchr (strptr1, '\n');
					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;

					/* expecting the form ' 'entryname'', delete the white spaces in data */
					/* (REMIND: entrynames itself are not allowed to contain whitespaces) */
					/* count the amount of characters between the apostrophes */
					/* (for security restrict j to numsigns) */
					/* allocate the names-field and write the entryname inside it */
					delallspc(data);
					j=1;
					while ( (data[j] != '\'') && (j < numsigns) ) { ++j; }

					cde[num_cif_dic_entries-1].num_names=1;
					cde[num_cif_dic_entries-1].names=(char **) calloc( 1, sizeof(char *)) ;
					cde[num_cif_dic_entries-1].names[0]=(char *) calloc( (j), sizeof(char)) ;

					for (i=1; i<j; ++i)
					{
						cde[num_cif_dic_entries-1].names[0][i-1]=data[i];
					}
					cde[num_cif_dic_entries-1].names[0][j-1]=0;

				}
				continue;
			}
		}
		/* close the file */
		fclose(inpf) ;
	}


	/* For the given molecular structure in the cif-file the connection table is generated by using :
		_geom_bond_atom_site_label_1
		_geom_bond_atom_site_label_2
	   Furthermore also read :
		_geom_bond_site_symmetry_1
		_geom_bond_site_symmetry_2
		_geom_bond_distance
	   Connection table has structure { {1,5},{2,3},... } where the numbers indicate the atoms.
	   Hence dimensions are #-of-bonds x 2.
	 */
	int read_connectivity_table_from_cif()
	{
		/* if an entry can't be found return -1 as an error flag */
		int i_label_1 = find_index_in_cfe( "_geom_bond_atom_site_label_1", 0) ;
		if ( i_label_1 < 0 ) { return -1 ; }
		int i_label_2 = find_index_in_cfe( "_geom_bond_atom_site_label_2", 0) ;
		if ( i_label_2 < 0 ) { return -1 ; }
		int i_symmetry_1 = find_index_in_cfe( "_geom_bond_site_symmetry_1", 0) ;
		if ( i_symmetry_1 < 0 ) { return -1 ; }
		int i_symmetry_2 = find_index_in_cfe( "_geom_bond_site_symmetry_2", 0) ;
		if ( i_symmetry_2 < 0 ) { return -1 ; }
		int i_distance = find_index_in_cfe( "_geom_bond_distance", 0) ;
		if ( i_distance < 0 ) { return -1 ; }

		int i_label = find_index_in_cfe( "_atom_site_label", 0) ;
		if ( i_label < 0 ) { return -1 ; }

		unsigned int noa_cif = cfe[i_label].data_length ;

		unsigned int n_bonds = cfe[i_label_1].data_length ;

		for ( unsigned int i=0; i<n_bonds; ++i)
		{
			bond_site12.push_back( vector<unsigned int>() ) ;
			for ( unsigned int j=0; j<noa_cif; ++j)
			{
				if ( !strcmp( cfe[i_label_1].text_data[i], cfe[i_label].text_data[j]) )
				{
					bond_site12[i].push_back(j) ;
					break ;
				}
			}
			for ( unsigned int j=0; j<noa_cif; ++j)
			{
				if ( !strcmp( cfe[i_label_2].text_data[i], cfe[i_label].text_data[j]) )
				{
					bond_site12[i].push_back(j) ;
					break ;
				}
			}
			bond_symmetry12.push_back( vector<string>() ) ;
			bond_symmetry12[i].push_back( cfe[i_symmetry_1].text_data[i] ) ;
			bond_symmetry12[i].push_back( cfe[i_symmetry_2].text_data[i] ) ;

			bond_distance.push_back( cfe[i_distance].numeric_data[i] ) ;
		}
		return 0 ;
	}


	/* write a cif-file dedicated for Jmol, PyMol, Rasmol, ... 
	   useful after a bt of the atom-coordinates,
	   simply for symmetry operations its not necessary because this can be easily done also in Jmol
	   
	   use entries the (more or less essential entries ):
		_symmetry_space_group_name_H-M
		_symmetry_space_group_name_Hall
		loop_
			_symmetry_equiv_pos_as_xyz
		_cell_length_a
		_cell_length_b
		_cell_length_c
		_cell_angle_alpha
		_cell_angle_beta
		_cell_angle_gamma
		_cell_formula_units_Z
		loop_
			_atom_site_label
			_atom_site_fract_x
			_atom_site_fract_y
			_atom_site_fract_z
			_atom_site_U_iso_or_equiv (optionally, if provided in cif-file)
			_atom_site_thermal_displace_type (optionally, if provided in cif-file)
			_atom_site_calc_flag (optionally, if provided in cif-file)
			_atom_site_calc_attached_atom (optionally, if provided in cif-file)
		loop_
			_geom_bond_atom_site_label_1 (optionally, if provided in cif-file)
			_geom_bond_atom_site_label_2 (optionally, if provided in cif-file)
			_geom_bond_site_symmetry_1 (optionally, if provided in cif-file)
			_geom_bond_site_symmetry_2 (optionally, if provided in cif-file)
			_geom_bond_distance (optionally, if provided in cif-file)

	   _symmetry and _cell entries are read from cif-file if they should differ, these entries must be modified
	   in one of the dedicated bt-functions or manually.

	   atom coordinates are used from atom struct.
	   label and bond information are derived from cif-entries if they exist there, it is assumed that by symmetry operations and bt
	   the noa is a multiple of noa in cif-file, in order to use the same labels and bonds for the equivalent molecules
	*/
	void write_cif_file_for_Jmol()
	{
		FILE *outf ;
		char sdummy[1024] ;
		int idummy ;

		unsigned int noa_cif = cfe[find_index_in_cfe("_atom_site_label")].data_length ;

		/* Open file */
		strcpy( sdummy, ofo) ;
		strcat( sdummy, mfname) ;
		strcat( sdummy, "_jmol_") ;
		strcat( sdummy, cif_file) ;
		if ( ( outf = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(64) ; }

		/* write obligatory data entry */
		fprintf( outf, "data_XNDiff\n" ) ;

		/* symmetry */
		strcpy( sdummy, "_symmetry_space_group_name_H-M") ;
		idummy = find_index_in_cfe(sdummy, 0) ;
		if ( idummy >= 0 ) { fprintf( outf, "%s     '%s'\n", sdummy, cfe[idummy].text_data[0] ) ; }

		strcpy( sdummy, "_symmetry_space_group_name_Hall") ;
		idummy = find_index_in_cfe(sdummy, 0) ;
		if ( idummy >= 0 ) { fprintf( outf, "%s     '%s'\n", sdummy, cfe[idummy].text_data[0] ) ; }

		strcpy( sdummy, "_symmetry_equiv_pos_as_xyz") ;
		idummy = find_index_in_cfe(sdummy, 0) ;
		if ( idummy >= 0 ) 
		{ 
			fprintf( outf, "loop_\n    %s\n", sdummy ) ;
			for ( unsigned int i=0; i<cfe[idummy].data_length; ++i)
			{
				fprintf( outf, "    %s\n", cfe[idummy].text_data[i] ) ;
			}
		}

		/* cell */
		strcpy( sdummy, "_cell_length_a") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { fprintf( outf, "%s     %.4lf\n", sdummy, cfe[idummy].numeric_data[0] ) ; }

		strcpy( sdummy, "_cell_length_b") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { fprintf( outf, "%s     %.4lf\n", sdummy, cfe[idummy].numeric_data[0] ) ; }

		strcpy( sdummy, "_cell_length_c") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { fprintf( outf, "%s     %.4lf\n", sdummy, cfe[idummy].numeric_data[0] ) ; }

		strcpy( sdummy, "_cell_angle_alpha") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { fprintf( outf, "%s     %.4lf\n", sdummy, cfe[idummy].numeric_data[0] ) ; }

		strcpy( sdummy, "_cell_angle_beta") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { fprintf( outf, "%s     %.4lf\n", sdummy, cfe[idummy].numeric_data[0] ) ; }

		strcpy( sdummy, "_cell_angle_gamma") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { fprintf( outf, "%s     %.4lf\n", sdummy, cfe[idummy].numeric_data[0] ) ; }

		/* atoms */
		fprintf( outf, "loop_\n") ;
		fprintf( outf, "    _atom_site_label\n") ;
		fprintf( outf, "    _atom_site_fract_x\n") ;
		fprintf( outf, "    _atom_site_fract_y\n") ;
		fprintf( outf, "    _atom_site_fract_z\n") ;

		/* optionally add thermal displacements etc */
		strcpy( sdummy, "_atom_site_U_iso_or_equiv") ;
		if ( find_index_in_cfe(sdummy, 0) >= 0 ) { fprintf( outf, "    %s\n", sdummy) ; }
		strcpy( sdummy, "_atom_site_thermal_displace_type") ;
		if ( find_index_in_cfe(sdummy, 0) >= 0 ) { fprintf( outf, "    %s\n", sdummy) ; }
		strcpy( sdummy, "_atom_site_calc_attached_atom") ;
		if ( find_index_in_cfe(sdummy, 0) >= 0 ) { fprintf( outf, "    %s\n", sdummy) ; }
		strcpy( sdummy, "_atom_site_calc_flag") ;
		if ( find_index_in_cfe(sdummy, 0) >= 0 ) { fprintf( outf, "    %s\n", sdummy) ; }

		for ( unsigned int i=0; i<noa; ++i)
		{
			fprintf( outf, "    %s %.5lf %.5lf %.5lf", atom[i].label, atom[i].coord[0], atom[i].coord[1], atom[i].coord[2] ) ;

			strcpy( sdummy, "_atom_site_U_iso_or_equiv") ;
			idummy = find_index_in_cfe(sdummy, 0) ;
			if ( idummy >= 0 ) { fprintf( outf, " %.3lf", cfe[idummy].numeric_data[i%noa_cif]) ; }
			strcpy( sdummy, "_atom_site_thermal_displace_type") ;
			idummy = find_index_in_cfe(sdummy, 0) ;
			if ( idummy >= 0 ) { fprintf( outf, " %s", cfe[idummy].text_data[i%noa_cif]) ; }
			strcpy( sdummy, "_atom_site_calc_attached_atom") ;
			idummy = find_index_in_cfe(sdummy, 0) ;
			if ( idummy >= 0 ) { fprintf( outf, " %s", cfe[idummy].text_data[i%noa_cif]) ; }
			strcpy( sdummy, "_atom_site_calc_flag") ;
			idummy = find_index_in_cfe(sdummy, 0) ;
			if ( idummy >= 0 ) { fprintf( outf, " %s", cfe[idummy].text_data[i%noa_cif]) ; }

			fprintf( outf, "\n" ) ;
		}

		/* bonds */
		if ( log_flag ) { fprintf( logfile, "\tread_connectivity_table_from_cif()\n") ; } 
		if ( read_connectivity_table_from_cif() < 0 )
		{
			if ( log_flag ) { fprintf( logfile, "\tread_connectivity_table_from_cif() failed due to at least one missing cif-entry. Bonds will be omitted.\n") ; } 
		}
		else
		{
			if ( log_flag ) { fprintf( logfile, "\tdone\n\tderive_connectivity_table_from_cif()\n") ; } 
			derive_connectivity_table_from_cif() ;
			if ( log_flag ) { fprintf( logfile, "\tdone\n") ; } 

			fprintf( outf, "loop_\n") ;
			fprintf( outf, "    _geom_bond_atom_site_label_1\n") ;
			fprintf( outf, "    _geom_bond_atom_site_label_2\n") ;
			fprintf( outf, "    _geom_bond_site_symmetry_1\n") ;
			fprintf( outf, "    _geom_bond_site_symmetry_2\n") ;
			fprintf( outf, "    _geom_bond_distance\n") ;

			for ( unsigned int i=0; i<bond_site12.size(); ++i)
			{
				fprintf( outf, "    %s %s %s %s %.3lf\n", atom[bond_site12[i][0]].label, atom[bond_site12[i][1]].label, bond_symmetry12[i][0].c_str(), bond_symmetry12[i][1].c_str(), bond_distance[i] ) ;
			}
		}

		fclose(outf) ;
	}


	/* Add the coherent scattering lengths and incoherent scattering cross sections 
	   from n_atom_sc to atom structure if neutron scattering is requested. */
	void add_neutron_coh_sc_length_and_inc_sc_crossec_to_atom()
	{
		/* add coherent scattering lengths and incoherent scattering cross sections from n_atom_sc to atom */
		for ( unsigned int i=0; i<noa; ++i)
		{
			atom[i].coh_sc_length = get_coherent_scattering_length (atom[i].type) ;
			atom[i].incoh_sc_cross = get_incoherent_scattering_crosssec (atom[i].type) ;
		}
	}


	int memmove_ptr(char *data, char *strptr1, char *strptr2 )
	{
		int i = (int)(strptr2-strptr1) ;
		if ( i > 0 )
		{
			memmove( data, strptr1, i) ;
			data[i] = 0 ;
		}
		return i ;
	}

	
	void vec_mult( vector<double> &x, vector<double> &a)
	{ 
		for ( unsigned int i=0; i<x.size(); ++i)
		{
			x.at(i) *= a.at(i) ;
		}
	}

	void vec_scale( vector<double> &x, double &a)
	{ 
		for ( unsigned int i=0; i<x.size(); ++i)
		{
			x.at(i) *= a ;
		}
	}


	/* Update atom[noa_ratio*noa_cif...(noa_ratio+1)*noa_cif].coord[coo] with the formula
	   The formula may comprise of several expressions separated by addition/ subtraction 
	   e.g. x'= -3/2 + 2.1*x or z'=-x or y'=-x+y+z*3/2  
	*/
	void formula_interpreter( char* formula, vector<double> &x, vector<double> &y, vector<double> &z, unsigned int coo, unsigned int noa_cif, unsigned int noa_ratio)
	{
		const int numsigns = 1024 ;
		char data[numsigns] ;
		char sdummy[numsigns] ;
		char *strptr1, *strptr2, *strptr3 ;
		int i ;
		double number ;
		vector<double> tmp ;

		bool next_expression ;
		bool number_def = false ;
		bool x_def = false ;
		bool x_def_all = false ;

		/* set atom[:+noa_ratio*noa_cif][coo] to 0.0 and add later each expression in the formula */
		for ( unsigned int j=0; j<noa_cif; ++j)
		{
			atom[j+noa_ratio*noa_cif].coord[coo] = 0.0 ;
		}

		strptr1 = formula ;

		do 
		{
			next_expression = false ;
			/* for each expression in a formula there might be a number and coordinate x,y,z */
			number_def = false ;
			x_def = false ;

			/* delete + at the beginning of each expression e.g. +x+2/3, +1.2*..., y+1/2 */
			if ( *strptr1 == '+' ) { ++strptr1 ; }

			/* if a + is found -> more than one expression e.g. x+1/2
			   similarly for - but ignore the first sign in formula e.g. -x+1/2
			*/
			if ( ( strptr2 = strchr( strptr1, '+') ) != NULL ) { next_expression = true ; }
			if ( ( strptr3 = strchr( strptr1+1, '-') ) != NULL ) { next_expression = true ; }
			/* if multiple expressions occur e.g. 1/2-x+y -> choose min for strptr2 */
			if ( strptr2 != NULL && strptr3 != NULL )
			{
				if ( (int)(strptr3-strptr2) < 0 ) { strptr2 = strptr3 ; }
			}
			/* if only - were found */
			if ( strptr2 == NULL && strptr3 != NULL ) { strptr2 = strptr3 ; strptr3 = NULL ; }

			/* only one expression (neither + nor - was found) */
			if ( next_expression == false ) { strptr2 = strchr( strptr1, '\0') ; }
			memmove_ptr( data, strptr1, strptr2) ;

			strptr3 = strptr2 ;

			/* check for numeric expression at beginning (constant addition or multiplication) */
			strptr1 = data ;
			strptr2 = strptr1 ;
			while ( isdigit(*strptr2) || *strptr2 == '.' || *strptr2 == '-' ) { ++strptr2; }

			i = memmove_ptr( sdummy, strptr1, strptr2) ;
			if ( i > 0 ) 
			{
				/* -, 2, 1.2394, -54.34 */
				if ( ( i == 1 ) && *strptr1 == '-' ) { number = -1.0 ; }
				else { number = strtod( sdummy, NULL) ; }

				strptr1 = strptr2 ;

				if ( *strptr1 == '/' )
				{
					/* 1/2, 2/4, -5/3
					   accept only integers > 0 */
					++strptr1 ;
					strptr2 = strptr1 ;
					while ( isdigit(*strptr2) ) { ++strptr2; }

					i = memmove_ptr( sdummy, strptr1, strptr2) ;
					if ( ( i > 0 ) && ( *strptr1 != '0' ) ) 
					{
						/* > 0 integers */
						number /= strtod( sdummy, NULL) ;
						strptr1 = strptr2 ;
					}
					else { XNDIFF_ERROR(76) ; }
				}

				tmp.assign( noa_cif, number) ;
				number_def = true ;
			}
			else
			{
				/* x, y, z */
				if ( *strptr1 == 'x') {  tmp = x ; }
				else if ( *strptr1 == 'y') {  tmp = y ; }
				else if ( *strptr1 == 'z') {  tmp = z ; }
				else { XNDIFF_ERROR(76) ; }
				++strptr1 ;
				x_def = true ;
			}

			/* evaluate second part of expression or update-continue */
			if ( ( *strptr1 != '+' ) && ( *strptr1 != '-' ) && ( *strptr1 != '\0' ) )
			{
				if ( ( *strptr1 == '*' ) && ( number_def ) )
				{
					/* ...*x, ...*y, ...*z */
					++strptr1 ;
					if ( *strptr1 == 'x') {  vec_mult( tmp, x) ; }
					else if ( *strptr1 == 'y') {  vec_mult( tmp, y) ; }
					else if ( *strptr1 == 'z') {  vec_mult( tmp, z) ; }
					else { XNDIFF_ERROR(76) ; }
					x_def = true ;
					++strptr1 ;
				}
				else if ( *(strptr1-1) == '-' )
				{
					/* -x, -y, -z, number=-1.0 */
					if ( *strptr1 == 'x') { vec_mult( tmp, x) ; }
					else if ( *strptr1 == 'y') { vec_mult( tmp, y) ; }
					else if ( *strptr1 == 'z') { vec_mult( tmp, z) ; }
					else { XNDIFF_ERROR(76) ; }
					x_def = true ;
					++strptr1 ;
				}
				else if ( ( *strptr1 == '*' ) && ( x_def ) )
				{
					/* x*..., y*..., z*..., read only positive numbers */

					++strptr1 ;
					strptr2 = strptr1 ;
					while ( isdigit(*strptr2) || *strptr2 == '.' ) { ++strptr2; }

					i = memmove_ptr( sdummy, strptr1, strptr2) ;
					if ( i > 0 ) 
					{
						/* 2, 1.2394 */
						number = strtod( sdummy, NULL) ;

						strptr1 = strptr2 ;

						if ( *strptr1 == '/' )
						{
							/* 1/2, 2/4
							   accept only integers > 0 */
							++strptr1 ;
							strptr2 = strptr1 ;
							while ( isdigit(*strptr2) ) { ++strptr2; }

							i = memmove_ptr( sdummy, strptr1, strptr2) ;
							if ( ( i > 0 ) && ( *strptr1 != '0' ) ) 
							{
								/* > 0 integers */
								number /= strtod( sdummy, NULL) ;
								strptr1 = strptr2 ;
							}
							else { XNDIFF_ERROR(76) ; }
						}
					}
					else { XNDIFF_ERROR(76) ; }

					vec_scale( tmp, number) ;
					number_def = true ;

// 					++strptr1 ;
				}
				else { XNDIFF_ERROR(76) ; }

				/* final check if end / next expression is now reached */
				if ( ( *strptr1 != '+' ) && ( *strptr1 != '-' ) && ( *strptr1 != '\0' ) ) { XNDIFF_ERROR(76) ; }
			}
			/* now update atom with first expression and continue */
			for ( unsigned int j=0; j<noa_cif; ++j)
			{
				atom[j+noa_ratio*noa_cif].coord[coo] += tmp.at(j) ;
			}
			/* if next_expression == false <==> *strptr1 == \0 
			   or for the next expression in the formula 
			   if next_expression == true <==> *strptr1 == + || *strptr1 == -
			   prepare strptr1 for the next expression
			*/
			strptr1 = strptr3 ;

			/* if coordinate x,y,z was included */
			if ( x_def ) { x_def_all = true ; } 
		}
		while ( next_expression ) ;

		/* each formula must contain at least one expression with x, y, z */
		if ( x_def_all == false ) { XNDIFF_ERROR(76) ; }
		return ;
	}


	/* use symmetry operations (e.g. "-x,-y,-z" or "x+1/2,-y+3/2,1-z") from:
	   -userdefined symmetry operation file, e.g. symop.txt (symop and cenop must be taken together)
	   -cif-file (symop and cenop are taken together by default in cif-files)
	   -syminfo.lib (symop and cenop will be read and saved in a ;-separated list 
	    e.g. x,-y,-z;x+1/2,y+1/2,1-z what is equivalent to x+1/2,-y+3/2,1-z )
	  
	   What formulas can be read:
	   x,y+1/2,z-1/2
	   -x+1,-y+1/2,-z+3/2
	   1.0*x,1.0*y+0.5,z*1.0-0.5
	   -1.0*x+1.0,-1.0*y+0.5,-1.0*z+1.5
	   x,y+1/4+1/4,-0.25+z-1/4
	   +x,+y,z
	   x*2.3333+1/10,y*2+0.5,z*1.0-0.5
	   -5.0*x+5.5,-1.0*y+0.5,-1.0*z+1.5
	   -5.0*x+2*y-5.5,-1.0*z+0.5,-1.0*x+1.5
	   -x-2.5*y,z+1/2,-x
	   
	   Of course also lists of them work:
	   -x,-y,-z;-x-2.5*y,z+1/2,-x
	   -x,-y,-z;2*x+0.5,y-1/2,z*2.0+5/3

	   Limitations:
	   an arbitrary number of linear terms in coordinates x,y,z AND offset constants are possible !
	   e.g. -1.2*x+y-2/3*z+1/4,-1/6-y+1/3,y+z+z  
	   
	   -x+y*2.5,z+1/2,-x works well
	   -x-y*2.5,z+1/2,-x does not work !!! -> -x-2.5*y,z+1/2,-x works
	   hence factors can only applied on the rhs of a coordinate x,y,z if there is no minus sign
	   
	   
	   Make a copy of current atom entries, what is necessary if several formulas are applied (e.g. symop AND cenop)
	   Applied on the coordinates given in the cif-file and appended to the atom structure
	*/
	void compute_coordinates_by_formula_interpreter( char* formula)
	{
		if ( log_flag ) { fprintf( logfile, "\tcompute_coordinates_by_formula_interpreter(\"%s\")\n", formula) ; fflush (logfile) ; }

		const int numsigns = 1024 ;
		char data[numsigns] ;
		char sep ;
		char *strptr1, *strptr2 ;
		int i ;
		bool next_formula_exists = true ;

		int i_label = find_index_in_cfe("_atom_site_label") ;
		int i_x = find_index_in_cfe("_atom_site_fract_x") ;
		int i_y = find_index_in_cfe("_atom_site_fract_y") ;
		int i_z = find_index_in_cfe("_atom_site_fract_z") ;

		vector<double> x, y, z ;

		unsigned int noa_cif = cfe[i_label].data_length ;
		unsigned int noa_ratio ;

		/* for each block of formula create a new bunch of coordinates in atom */
		if ( noa == 0 ) { noa += noa_cif ; atom = (atom_entries *) calloc( noa, sizeof(atom_entries)) ; }
		else { noa += noa_cif ; atom = (atom_entries *) realloc( atom, noa * sizeof(atom_entries)) ; }

		noa_ratio = noa / noa_cif - 1 ;
			
		/* initialize atom with coordinates from cif-file
		   they will be updated subsequently by all formulas given in formula 
		*/
		for ( unsigned int j=0; j<noa_cif; ++j)
		{
			atom[j+noa_ratio*noa_cif].coord[0] = cfe[i_x].numeric_data[j]  ;
			atom[j+noa_ratio*noa_cif].coord[1] = cfe[i_y].numeric_data[j] ;
			atom[j+noa_ratio*noa_cif].coord[2] = cfe[i_z].numeric_data[j] ; 
		}

		delallspc(formula) ;
		strptr1 = formula ;
		do 
		{
			/* for each formula in formula we need a copy of atom from the previous update
			   before we overwrite atom again
			   furthermore x,y,z need to be cleared if several formulas are used (;-separated list)
			*/
			if ( !(x.empty()) ) { x.clear() ; }
			if ( !(y.empty()) ) { y.clear() ; }
			if ( !(z.empty()) ) { z.clear() ; }
			for ( unsigned int j=0; j<noa_cif; ++j)
			{
				x.push_back( atom[j+noa_ratio*noa_cif].coord[0] ) ;
				y.push_back( atom[j+noa_ratio*noa_cif].coord[1] ) ;
				z.push_back( atom[j+noa_ratio*noa_cif].coord[2] ) ;
			}

			for ( unsigned int j=0; j<3; ++j)
			{
				if (j<2) { sep = ',' ; }
				else { sep = ';' ; }

				if ( ( strptr2 = strchr( strptr1, sep) ) == NULL )
				{
					if ( j<2 ) { XNDIFF_ERROR(76) ; }
					else
					{
						/* j==2 */
						next_formula_exists = false ;
						if ( ( strptr2 = strchr( strptr1, '\0') ) == NULL ) { XNDIFF_ERROR(76) ; }
					}
				}

				i = (int)(strptr2-strptr1) ;
				memmove( data, strptr1, i) ;
				data[i] = 0 ;
				/* data contains now strings for x'_j e.g. -y+1/2 */

				/* update atom[...].coord[j] with the current formula for x'_j and x,y,z */
				formula_interpreter( data, x, y, z, j, noa_cif, noa_ratio) ;

				if ( j < 2 || next_formula_exists ) { strptr1 = strptr2 + 1 ; }
			}
		}
		while ( next_formula_exists ) ; 

		if ( log_flag ) { fprintf( logfile, "\tdone\n") ; fflush(logfile) ; }
	}


	/* -sym option */
	void compute_atoms_by_symmetry()
	{
		char sdummy[1024] ;
		int idummy ;

		/* all modes are fully complementary !!! */ 
		if ( sym_cif_mode)
		{
			strcpy( sdummy, "_symmetry_equiv_pos_as_xyz") ;
			idummy = find_index_in_cfe(sdummy) ;
			if ( idummy >= 0 )
			{
				if ( log_flag ) { fprintf( logfile, "\tUse symmetry operations (_symmetry_equiv_pos_as_xyz) from cif-file %s:\n", cif_file) ; fflush (logfile) ; }
				for ( unsigned int i=0; i<cfe[idummy].data_length; ++i)
				{
					compute_coordinates_by_formula_interpreter( cfe[idummy].text_data[i] );
				}
			}
			else { XNDIFF_ERROR(68) ; }
		}

		if ( sym_ccp4_mode )
		{
			if ( log_flag ) { fprintf( logfile, "\tUse CCP4 (pdbset):\n") ; fflush (logfile) ; }
			if ( log_flag ) { fprintf( logfile, "\twrite_pdb_file(...) : %s\n", cif_pdb_file) ; fflush (logfile) ; }
			/*
				int index_atom_site_label=find_index_in_cfe("_atom_site_fract_x") ;
				int index_atom_site_fract_x=find_index_in_cfe("_atom_site_fract_x") ;
				int index_atom_site_fract_y=find_index_in_cfe("_atom_site_fract_y") ;
				int index_atom_site_fract_z=find_index_in_cfe("_atom_site_fract_y") ;
			*/
			write_pdb_file( cfe[find_index_in_cfe("_atom_site_fract_x")].data_length, cfe[find_index_in_cfe("_atom_site_fract_x")].numeric_data, cfe[find_index_in_cfe("_atom_site_fract_y")].numeric_data, cfe[find_index_in_cfe("_atom_site_fract_z")].numeric_data, cfe[find_index_in_cfe("_atom_site_label")].text_data, NULL) ;

			if ( log_flag ) { fprintf( logfile, "\tdone\n\n") ; fprintf( logfile, "\twrite_ccp4_bash_script() : %s\n", ccp4_bash_file) ; fflush (logfile) ; }
			write_ccp4_bash_script() ;
			if ( log_flag ) { fprintf( logfile, "\tdone\n\n") ; fprintf( logfile, "\trun_ccp4_bash_script()\n") ; fflush (logfile) ; }
			run_ccp4_bash_script() ;
			if ( log_flag ) { fprintf( logfile, "\tdone\n\n") ; fflush (logfile) ; }
			/* Read atomic informations from CCP4-generated pdb-file into structure atom */
			if ( log_flag ) { fprintf( logfile, "\tRead CCP4-generated pdb-file %s in read_pdb_file()\n", pdb_file) ; fflush (logfile) ; }
			read_pdb_file() ;
			if ( log_flag ) { fprintf( logfile, "\tdone\n") ; }
		}

		if ( sym_pdb_mode )
		{
			/* Read atomic informations from user-supplied pdb-file into structure atom */
			if ( log_flag ) { fprintf( logfile, "\tRead user-supplied pdb-file %s in read_pdb_file()\n", pdb_file) ; fflush (logfile) ; }
			read_pdb_file() ;
			if ( log_flag ) { fprintf( logfile, "\tdone\n") ; }
		}

		/* read symmetry operations from syminfo.lib or a user-defined symop-file
		   for the former a space_group name (or number, not yet supported) from the cif-file must be provided
		   Use symop and cenop entries in syminfo.lib, but nothing more like basisop and no limit control !
		*/
		if ( !sym_ccp4_mode && !sym_cif_mode && !sym_pdb_mode )
		{
			if ( symop_userdef )
			{
				read_symopxyz_from_file() ;
				for ( unsigned int i=0; i<num_symop; ++i)
				{
					compute_coordinates_by_formula_interpreter( symopxyz[i] );
				}
			}
			else
			{
				vector<string> cif_symmetry_keywords ;
				cif_symmetry_keywords.push_back("_symmetry_space_group_name_H-M") ;
				cif_symmetry_keywords.push_back("_symmetry_space_group_name_Hall") ;
				/* cif_symmetry_keywords.push_back("_symmetry_Int_Tables_number") ; */

				vector<string> syminfolib_symmetry_keywords ;
				syminfolib_symmetry_keywords.push_back("symbol xHM") ;
				syminfolib_symmetry_keywords.push_back("symbol Hall") ;
				/* syminfolib_symmetry_keywords.push_back("number") ; */

				bool SYMMETRY_SUCCESS = false ;
				bool SYMINOLIB_ENTRY_FOUND = false ;
				char pattern[1024] ;
				for ( unsigned int i=0; i<cif_symmetry_keywords.size(); ++i)
				{
					idummy = find_index_in_cfe( cif_symmetry_keywords.at(i).c_str() ) ;
					sprintf( pattern, "'%s'", cfe[idummy].text_data[0]) ;
					if ( idummy >= 0 )
					{
						if ( log_flag ) { fprintf( logfile, "\tread_symopxyz_from_syminfolib( %s, %s)\n", syminfolib_symmetry_keywords.at(i).c_str(), pattern) ; fflush (logfile) ; }
						SYMINOLIB_ENTRY_FOUND = read_symopxyz_and_cenopxyz_from_syminfolib( syminfolib_symmetry_keywords.at(i).c_str(), pattern ) ;

						if ( SYMINOLIB_ENTRY_FOUND == false ) 
						{ 
							if ( log_flag ) { fprintf( logfile, "\tfailed\n") ; fflush (logfile) ; }
							continue ;
						}

						if ( log_flag ) { fprintf( logfile, "\tdone\n") ; fflush (logfile) ; }

						for ( unsigned int j=0; j<num_cenop; ++j)
						{
							for ( unsigned int i=0; i<num_symop; ++i)
							{
								/* 0\ will be automatically appended by sprintf */
								sprintf( sdummy, "%s;%s", symopxyz[i], cenopxyz[j]) ;
								compute_coordinates_by_formula_interpreter( sdummy ) ;
							}
						}

						SYMMETRY_SUCCESS = true ;
						break ;
					}
				}
				if ( SYMMETRY_SUCCESS == false ) { XNDIFF_ERROR(69); }
			}
		}

		/* derive thermal displacements and multiplicity with the cif-file
		   presently not yet implemented,they are set to:
		   atom[:].mult = 1.0 (reasonable)
		   atom[:].temp[:] = 0.0

		   Note that temp and mult values from a read_pdb_file() will be overwritten -> doesn't matter 
		*/
		derive_Uiso_Uani_and_Mult_from_cif() ;

		/* derive atom labels and types with the cif-file */
		derive_atom_labels_and_types_from_cif() ;

		free_symopxyz_and_cenopxyz() ;
	}


	/* reads cif-file information with the help of the cif-dictionary cde into variable cfe.
	   also reads the cell parameters a,b,c,alpha,beta,gamma and HM-symbol from the cif-file entries cfe and scales a,b,c from [A] to [nm].

		not supported yet:
		for both _single and loop_ entries
		;
			multi-line			
			enviroments with 
			starting and terminating ';'
			are ignored and not saved as text entries !
		;

		for a mix of invalid + valid entries in a loop_, ALL are ignored !
		as long ALL entries in a loop_ are EITHER valid OR invalid there's no problem.

		in some cif-files entries from different dictionaries (powder diffraction, _pd; macro molecules, ...) can be present, therefore it would be good to allow the rading of several dictionaries.
		dictionaries:
		core 			http://www.iucr.org/resources/cif/dictionaries/cif_core
		powder diffraction 	http://www.iucr.org/resources/cif/dictionaries/cif_pd
		macromolecules 		http://www.iucr.org/resources/cif/dictionaries/cif_mm
		...


	*/
	void read_cif_file()
	{
		FILE *inpf ; /* FILE pointer, for the cif-dictionary */
		const int numsigns = 1024 ; /* in principle 80/81 should be sufficient, because cif files are restricted to 80 characters in text-width */
		int i ;
		char sdummy[numsigns] ;
		char data[numsigns] ;

		bool entry_in_dic ;
		char *strptr1, *strptr2 ;
		long int lineindex = 0 ;

		/* Open cif-file */
		strcpy( sdummy, ifo) ;
		strcat( sdummy, cif_file) ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(22) ; }

		/* read line by line */
		bool loop_flag = false ;
		bool loop_data_flag = false ;
		bool loop_search_next_line = false ;
		bool single_flag = false ;
		bool search_next_line = false ;


		int num_loop_entries = 0 ;
		int entries_before_loop = 0 ;
		int current_entry = 0 ;
		
		num_cif_file_entries = 0 ;
		while ( fgets( sdummy, numsigns, inpf) != NULL)
		{
			++lineindex;

			/* ignore comments starting with '#'' ,empty lines, additional white space in front of them will be always ignored  */
			strptr1 = sdummy ;
			while (*strptr1 == ' ')
				++strptr1 ;

			if (*strptr1 == '#')
				continue ;
			if (*strptr1 == '\n')
				continue ;



			/* CHECK for single entries starting with '_' */
			if ( *strptr1 == '_' )
			{
				/* if a '_' is found if loop_data_flag=true, a new single data block starts */
				if (loop_data_flag)
				{
					loop_flag=false;
					loop_data_flag=false;
				}

				if (!loop_flag)
				{
					if ((strptr2=strchr (strptr1, ' ')) == NULL )
					{
						strptr2=strchr (strptr1, '\n');
						search_next_line=true;
					}

					entry_in_dic=false;
					/* read entry label */
					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;

					if (!search_next_line)
					{
						strptr1=strptr2;
		
						while (*strptr1 == ' ')
							++strptr1 ;
					}
	
					delallspc(data);
	
					/* compare entry label with dictionary */
					/* if none is found ignore it and further proceed */
					for ( unsigned int k=0; k<num_cif_dic_entries; ++k)
					{
						for ( int l=0; l<cde[k].num_names; ++l)
						{
							if (!strcmp(data,cde[k].names[l]))
							{
								entry_in_dic=true;
								single_flag=true;
	
								if ( num_cif_file_entries == 0 )
								{
									cfe = (cif_file_entries *) calloc( 1, sizeof(cif_file_entries)) ;
								}
								else
								{
									cfe = (cif_file_entries *) realloc (cfe, (num_cif_file_entries + 1) * sizeof(cif_file_entries)) ;
								}
								++num_cif_file_entries;
	
								/* copy informations from cde to cfe */
								cfe[num_cif_file_entries-1].units = (char *) calloc( strlen(cde[k].units)+1, sizeof(char)) ;
								strcpy(cfe[num_cif_file_entries-1].units,cde[k].units) ;
				
								cfe[num_cif_file_entries-1].name = (char *) calloc( strlen(cde[k].names[l])+1, sizeof(char)) ;
								strcpy(cfe[num_cif_file_entries-1].name,cde[k].names[l]) ;
								
								cfe[num_cif_file_entries-1].isnumeric = cde[k].isnumeric ;
								cfe[num_cif_file_entries-1].data_length = 1 ;
	
								break;
							}
						}
						if (entry_in_dic)
							break;
					}
					if (!entry_in_dic)
					{
						/* error the entry label was not found in the cif_dic */
						if ( log_flag ) { fprintf( logfile, "\tThe label '%s' in line %ld of the cif-file could not be found in the cif-dictionary. Ignoring it.\n", data, lineindex) ; }
						search_next_line=false;
						continue;
					}
				}
				else
				{
					++num_loop_entries;

					strptr2=strchr (strptr1, '\n');

					/* read entry label */
					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;

					delallspc(data) ;

					entry_in_dic=false ;
					/* compare entry label with dictionary */
					/* if none is found ignore it and further proceed */
					for (unsigned int k=0; k<num_cif_dic_entries; ++k)
					{
						for (int l=0; l<cde[k].num_names; ++l)
						{
							if (!strcmp(data,cde[k].names[l]))
							{
								entry_in_dic=true;
								if ( num_cif_file_entries == 0 )
								{
									cfe = (cif_file_entries *) calloc( 1, sizeof(cif_file_entries)) ;
								}
								else
								{
									cfe = (cif_file_entries *) realloc (cfe, (num_cif_file_entries + 1) * sizeof(cif_file_entries)) ;
								}
								++num_cif_file_entries;


								/* copy informations from cde to cfe */
								cfe[num_cif_file_entries-1].units = (char *) calloc( strlen(cde[k].units)+1, sizeof(char)) ;
								strcpy(cfe[num_cif_file_entries-1].units,cde[k].units);
				
								cfe[num_cif_file_entries-1].name = (char *) calloc( strlen(cde[k].names[l])+1, sizeof(char)) ;
								strcpy(cfe[num_cif_file_entries-1].name,cde[k].names[l]);
								
								cfe[num_cif_file_entries-1].isnumeric = cde[k].isnumeric;
								cfe[num_cif_file_entries-1].data_length=0;

								break;
							}
						}
						if (entry_in_dic)
							break;
					}
					if (!entry_in_dic)
					{
						/* error the entry label was not found in the cif_dic */
						if ( log_flag ) { fprintf( logfile, "\tThe label '%s' (inside a loop) in line %ld of the cif-file could not be found in the cif-dictionary. Ignoring it. All other labels in this loop are ignored too!\n", data, lineindex) ; }

						/* reset/correct values and free already allocated cfe[i] text fields and cfe[i] from valid entries */
						for (unsigned int k=0; k< ( num_cif_file_entries - entries_before_loop ); ++k)
						{
							free(cfe[k].name);
							free(cfe[k].units);
							/* free(cfe[k]); (or reallocate) not possible / necessary */
						}
						num_cif_file_entries = entries_before_loop ;
						loop_flag = false ;
					}

					/* get the next line */
					continue;
				}
			}

			/* CHECK for 'loop_' */
			if ( strstr(strptr1, "loop_") != NULL)
			{
				loop_flag = true ;
				loop_data_flag = false ;
				num_loop_entries = 0 ;
				/* reset current_entry to 0 when next loop_ starts */
				current_entry = 0 ;
				entries_before_loop = num_cif_file_entries ;
				continue ;
			}

			/* read data from a loop if no further entry starting with '_' below loop_ is found */
			if (loop_flag)
			{
				loop_data_flag=true;
				loop_search_next_line=false;

				for (int k=current_entry; k<num_loop_entries; ++k)
				{
					++cfe[entries_before_loop+k].data_length;

					/* Initialize strptr2 with the NULL Pointer */
					strptr2=NULL;

					if ( k < num_loop_entries-1 )
					{
						if ( (!cfe[entries_before_loop+k].isnumeric) && (*(strptr1) == '\'') )
						{
							/* text/char starting with ' */
							strptr2=strptr1+1;
							while ( !((*strptr2 == '\'') && (*(strptr2+1) == ' ')) )
							{
								++strptr2;
								if (*strptr2=='\n')
								{
									/* if "'\n" occurs instead of "' " */
									strptr2=NULL;
									break;
								}
							}

							if (strptr2!=NULL)
								++strptr2;
						}
						else
						{
							/* numeric and text/char not starting with ' */
							strptr2=strchr (strptr1, ' ');
						}

						/* if no expected entry was found */
						if (strptr2==NULL)
						{
							/* data entries missing on this line check the next one and so on */
							strptr2=strchr (strptr1, '\n');
							loop_search_next_line=true;
						}
					}
					else
					{
						strptr2=strchr (strptr1, '\n');
					}

					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;

					if ( cfe[entries_before_loop+k].isnumeric )
					{
						if ( cfe[entries_before_loop+k].data_length == 1 )
						{
							cfe[entries_before_loop+k].numeric_data=(double* ) calloc( 1, sizeof(double));
							cfe[entries_before_loop+k].numeric_stderr_data=(double* ) calloc( 1, sizeof(double));
						}
						else
						{
							cfe[entries_before_loop+k].numeric_data=(double* ) realloc ( cfe[entries_before_loop+k].numeric_data, cfe[entries_before_loop+k].data_length * sizeof(double));
							cfe[entries_before_loop+k].numeric_stderr_data=(double* ) realloc ( cfe[entries_before_loop+k].numeric_stderr_data, cfe[entries_before_loop+k].data_length * sizeof(double));
						}

						str2numeric( data, cfe[entries_before_loop+k].numeric_data[cfe[entries_before_loop+k].data_length-1], cfe[entries_before_loop+k].numeric_stderr_data[cfe[entries_before_loop+k].data_length-1] ); 	
					}
					else
					{
						/* char/text data */
						if ( cfe[entries_before_loop+k].data_length == 1 )
						{
							cfe[entries_before_loop+k].text_data=(char** ) calloc( 1, sizeof(char *));
						}
						else
						{
							cfe[entries_before_loop+k].text_data=(char** ) realloc ( cfe[entries_before_loop+k].text_data, cfe[entries_before_loop+k].data_length * sizeof(char *)) ;
						}

						readtext(inpf, lineindex, data, cfe[entries_before_loop+k].text_data[cfe[entries_before_loop+k].data_length-1]);
					}

					if (loop_search_next_line)
					{
						current_entry=k+1;
						break;
						/* continue with next line and reading from current_entry on */
					}

					/* before next entry set pointer strptr1 to the position of strptr2 and delete possibly white spaces */
					strptr1=strptr2;
					while (*strptr1 == ' ') 
						++strptr1 ;

				}
				/* end of for-loop, continue */

				/* reset current_entry to 0 when next data line starts */
				if (!loop_search_next_line)
					current_entry=0;
			}

			/* read data from single lines */
			if (single_flag)
			{
				if (search_next_line)
				{
					search_next_line=false;
					continue;
				}
				else

				/* read the entry until EOL */
				strptr2=strchr (strptr1, '\n');

				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;

				if ( cfe[num_cif_file_entries-1].isnumeric )
				{
					/* numeric input possibly with standard error in (..) */
					cfe[num_cif_file_entries-1].numeric_data=(double* ) calloc( 1, sizeof(double));
					cfe[num_cif_file_entries-1].numeric_stderr_data=(double* ) calloc( 1, sizeof(double));

					/* extend later with error flag coming from str2numeric */
					str2numeric( data, cfe[num_cif_file_entries-1].numeric_data[0], cfe[num_cif_file_entries-1].numeric_stderr_data[0]);
				}
				else
				{
					/* char/text data */
					cfe[num_cif_file_entries-1].text_data=(char **) calloc( 1, sizeof(char *)) ;

					readtext(inpf, lineindex, data, cfe[num_cif_file_entries-1].text_data[0]);
				}

				single_flag=false;
			}

		/* end of outer while loop */
		}
		/* close the file */
		fclose(inpf) ;

		/* print number of found cfe */
		if ( log_flag ) { fprintf( logfile, "\tFound %d entries in cif-file.\n", num_cif_file_entries) ; fflush (logfile) ; }	

		/* read cell parameters a,b,c,alpha,beta,gamma from cfe and scale from [A] to [nm] */
		cell_par[0]=read_single_numeric_data_from_cfe("_cell_length_a") / 10.0;
		cell_par[1]=read_single_numeric_data_from_cfe("_cell_length_b") / 10.0;
		cell_par[2]=read_single_numeric_data_from_cfe("_cell_length_c") / 10.0;

		cell_par[3]=read_single_numeric_data_from_cfe("_cell_angle_alpha");
		cell_par[4]=read_single_numeric_data_from_cfe("_cell_angle_beta");
		cell_par[5]=read_single_numeric_data_from_cfe("_cell_angle_gamma");

	/* end of read_cif_file() */
	}


	/* returns true if data contains a number representation (0,1,...,9,+,-,.,e,E)
	   e.g. -1.056, 2, 2.56E+07 
	   but also for e.g. E-0.4+4.098E98....-+e
	*/
	bool is_numeric(char* data, bool with_exp=true)
	{
		bool bdummy = true ;
		unsigned int i = 0 ;
		while ( data[i] != 0 )
		{
			if ( isdigit(data[i]) ) { ++i ; continue ; }
			if ( data[i] == '.' ) { ++i ; continue ; }
			if ( ( data[i] == '-' ) || ( data[i] == '+' ) ) { ++i ; continue ; }
			if ( with_exp ) {  if ( ( data[i] == 'e' ) || ( data[i] == 'E' ) ) { ++i ; continue ; } }
			bdummy = false ;
			break ;
		}
		return bdummy ;
	}


	/* returns true if data contains a positive integer number
	   e.g. 90, +87, 0, +0009876 
	*/
	bool is_unsigned_int(char* data)
	{
		bool bdummy = true ;
		unsigned int i = 0 ;
		if ( data[0] == '+' ) { ++i ; }
		while ( data[i] != 0 )
		{
			if ( isdigit(data[i]) ) { ++i ; continue ; }
			bdummy = false ;
			break ;
		}
		return bdummy ;
	}


	/* convert string from cif-file to number and its e.s.d. */
	void str2numeric(char* str, double &number, double &stderr) 
	{

		/* delete possible white space after the numeric input */
		delallspc(str);

		const int numsigns=82;
		char *strptr1, *strptr3;
		bool is_point=false;
		bool is_exp=false;
		bool is_bracket=false;
		char data[numsigns];

		int pos_point=0, pos_exp=0, pos_bracket=0, pos_end=0 ;
		long int exponent;

		strptr1=str;
		strptr3=str;

		while (*strptr3 != '\0')
		{
			if (*strptr3 == '.')
			{
				is_point = true;
				pos_point = (int)(strptr3-strptr1) ;
			}
			else if ((*strptr3 == 'e') || (*strptr3 == 'E'))
			{
				is_exp = true;
				pos_exp = (int)(strptr3-strptr1) ;
			}
			else if (*strptr3 == '(')
			{
				is_bracket = true;
				pos_bracket = (int)(strptr3-strptr1) ;
			}
			else if (*strptr3 == ')')
			{
				pos_end = (int)(strptr3-strptr1) ;
			}
			++strptr3;
		}

		if (is_bracket)
		{
			/* expecting 1.9876(4), -.0098765(32), 37(2), 1.567E3(3), 1E-4(3), ... */
			/* rewrite data up to the position of '(' and read the number */
			memmove (data, strptr1, pos_bracket) ;
			data[pos_bracket] = 0 ;
			number=strtod(data,NULL);

			/* read e.s.d. between (...) */
			memmove (data, strptr1+pos_bracket+1, pos_end-1-pos_bracket) ;
			data[pos_end-1-pos_bracket] = 0 ;
			stderr=strtod(data,NULL);

			if (is_exp)
			{
				if (is_point)
					stderr *= pow(10.0,(double)(-(pos_exp-1-pos_point)));

				memmove (data, strptr1+pos_exp+1, pos_bracket-1-pos_exp) ;
				data[pos_bracket-1-pos_exp] = 0 ;
				exponent=strtol(data,NULL,10) ;

				stderr *= pow(10.0,(double)(exponent)) ;
			}
			else
			{
				if (is_point)
					stderr *= pow(10.0,(double)(-(pos_bracket-1-pos_point)));
			}

		}
		else
		{
			/* expecting 1.9876, -.0098765, 1.567E4, ... */
			/* set e.s.d. to 0.0 and read number from data from above */ 
			stderr=0.0;
			number=strtod(str,NULL);
		}
	}


	/* read text from cif-file, call with parameters by reference */
	void readtext(FILE* inpf, long int& lineindex, char* str, char*& text)
	{
		const int numsigns=82;

		/* str does not contain any white spaces in front, but possibly in the back */
		/* last sign of str should be '\0' and not '\n' */
		if (str[0]=='\'')
		{
		
			int j=1;
			while ( (str[j] != '\'') && ( (str[j+1] != ' ') || (str[j+1] != '\0') ) && (j < numsigns) )
				++j;

			text=(char *) calloc( (j-1), sizeof(char)) ;
			for (int i=1; i<j; ++i)
			{
				text[i-1]=str[i];
			}
		}
		else if (str[0]==';')
		{
			char sdummy[numsigns];
			char data[numsigns];
			char *strptr1, *strptr2;
			int line=0;
			int length=0;
			int totallength=0;

			while (fgets (sdummy, numsigns, inpf) != NULL)
			{
				++lineindex;
	
				/* ignore comments starting with '#'' ,empty lines, additional white space in front of them will be always ignored  */
				strptr1 = sdummy ;
				while (*strptr1 == ' ')
					++strptr1 ;
	
				if (*strptr1 == '#')
					continue ;
				if (*strptr1 == '\n')
					continue ;

				if (*strptr1==';')
				{
					/* write the terminating NULL signs before leaving */
					if (totallength>0)
						text[totallength-1]='\0';

					break;
				}

				/* read the text until EOL */
				strptr2=strchr (strptr1, '\n');

				int i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;

				length = strlen(data)+1 ;

				if (line==0)
					text=(char *) calloc( length , sizeof(char) ) ;
				else
					text=(char *) realloc ( text, (totallength + length) * sizeof(char) ) ;

				for (int i=0; i<length; ++i)
				{
					text[totallength+i]=data[i];
				}
				/* write EOL sign */
				totallength += length;
				text[totallength-1]='\n';

				++line;
			}
		}
		else
		{
			delallspc(str);
			text=(char *) calloc( strlen(str)+1, sizeof(char)) ;
			strcpy(text,str);
		}
	}


	/* clears cde */
	void free_cif_dictionary()
	{
		for ( unsigned int i=0; i<num_cif_dic_entries; ++i)
		{
			free(cde[i].category);
			free(cde[i].units);
			for (int j=0; j<cde[i].num_names; ++j)
				free(cde[i].names[j]);
	
		}
		free(cde);
	}


	/* clears cfe */
	void free_cif_file()
	{
		for ( unsigned int i=0; i<num_cif_file_entries; ++i)
		{
			free(cfe[i].name);
			free(cfe[i].units);

			if (cfe[i].isnumeric)
			{
				free(cfe[i].numeric_data);
				free(cfe[i].numeric_stderr_data);
			}
			else
			{
				for ( unsigned int j=0; j<cfe[i].data_length; ++j)
					free(cfe[i].text_data[j]);
			}
		}
		free(cfe);
	}


	/* returns the index i (0<=i<=num_cif_file_entries) of entry in the cif-file database cfe 
	   if entry is not found i=-1. If warn_or_exit >= 0 only a warning message will be written, 
	   for warn_or_exit < 0 (default value) program will terminate */
	int find_index_in_cfe( const char* entry, int warn_or_exit = -1 )
	{
		int i=-1;

		for ( unsigned int j=0; j<num_cif_file_entries; ++j)
		{
			if ( !strcmp( entry, cfe[j].name) )
			{
				i = (int)j ;
				break ;
			}
		}

		if ( i<0 )
		{
			if ( warn_or_exit < 0 ) 
			{ 
				if ( log_flag ) { fprintf( logfile, "Error: Could not find the entry '%s' in cif-file. Exit.\n", entry) ; fflush(logfile) ; }
				XNDIFF_ERROR(30); 
			}
			else
			{ 
				if ( log_flag ) { fprintf( logfile, "Warning: Could not find the entry '%s' in cif-file.\n", entry) ; fflush(logfile) ; }
			}
		}

		return i;
	}


	/* returns numerical data for a specific entry in the cfe */
	double read_single_numeric_data_from_cfe(const char* entry, int warn_or_exit = -1 )
	{
		int index=find_index_in_cfe( entry, warn_or_exit);

		if ( !cfe[index].isnumeric )
		{
			if ( log_flag ) { fprintf( logfile, "\tError in read_single_numeric_data_from_cfe(). The data type for the entry %s is not numeric.\n", entry) ; }
			XNDIFF_ERROR(79) ;
		}
		return cfe[index].numeric_data[0];
	}


	/* returns text data for a specific entry in the cfe */
	void read_single_text_data_from_cfe(const char* entry, char** text_data, int warn_or_exit = -1 )
	{
		int index=find_index_in_cfe( entry, warn_or_exit);

		if ( cfe[index].isnumeric )
		{
			if ( log_flag ) { fprintf( logfile, "\tError in read_single_text_data_from_cfe(). The data type for the entry %s is not text.\n", entry) ; }
			XNDIFF_ERROR(80) ;
		}
		*text_data=(char *) calloc( strlen(cfe[index].text_data[0])+1, sizeof(char) );
		strcpy( *text_data, cfe[index].text_data[0]);
	}


	/* Writes the necessary information from the cif file (atomic positions,...) to a pdb file cif_pdb_file.
	   Later on, this file will be processed by the CCP4 routine pdbset specified in a bash script.

	   A Multiplicity of 1 is assumed for every atom in the cif-file !!!
	*/
	void write_pdb_file( unsigned int N, double* x, double* y, double* z, char** label, char** type, unsigned int mode=0)
	{
		FILE *inpf ;
		char sdummy[1024], symmetry_space_group_name_HM[1024] ;
		int idummy ;

		/* write the bash-script */
		strcpy( sdummy, ofo) ;
		strcat( sdummy, cif_pdb_file) ;

		if ( ( inpf = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(28) ; }

		fprintf(inpf,"REMARK   1 PDB file created from XNDIFF                               \n");

		/* PDB Format used here for entry CRYST1, there are 70 columns in a line.
		Here a virtual cubic crystal with a=b=c=1.000 is used.
		Therefore the fractional coordinates coming from the cif-file are identical
		to the orthogonal coordinates (Transformation Matrix is unit matrix).
		PDBSET usually expects orthogonal coordinates, transforms them with their
		information given in CRYST1 (a,b,c,alpha,beta,gamma) to fractional coordinates,
		and applies then SYMGEN on the fractional coordinates. After this the inverse 
		transformation is done to obtain the orthogonal coordinates. In short:
	
		x_frac'=(M)^(-1).M.SYMGEN.(M)^(-1).M.x_frac;
		( where  x_orth=M.x_frac;  PDBSET=M.SYMGEN.(M)^(-1); )
	
		To avoid the calculation of M (and because it's senseless) set up virtually
		a unit sized cubic crystal such that x_orth=x_frac;
	
				1111111111222222222233333333334444444444555555555566666666667\n
			1234567890123456789012345678901234567890123456789012345678901234567890\n
			CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 SPACE_GROUP    \n   
			
			where SPACE_GROUP must be 'left justified' !
		*/

		strcpy( sdummy, "_symmetry_space_group_name_H-M") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { strcpy( symmetry_space_group_name_HM, cfe[idummy].text_data[0]) ; }
		else { strcpy( symmetry_space_group_name_HM, "Unknown") ; }

		/* don't use cell parameters from cif-file, see remarks above */
		fprintf(inpf,"CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 %s", symmetry_space_group_name_HM );
		/* fill right hand side with additional white spaces if necessary */
		for ( unsigned int j=0; j < ( 11 - (unsigned int)strlen(symmetry_space_group_name_HM) ); ++j)
			fprintf(inpf," ");

		fprintf(inpf,"    \n");


		/* PDB Format used here for entry ATOM, there are 70 columns in a line.
		The fractional coordinates from the cif-file have the form +X.XXYYY(...).
		For best approxination in the PDB format scale them with x100 such that 
		they suit to the PDB format +XXX.YYY. Later rescale them by 0.01 !
	
				1111111111222222222233333333334444444444555555555566666666667\n
			1234567890123456789012345678901234567890123456789012345678901234567890\n
			ATOM  ##### TYPE   - -   0    +XXX.YYY+XXX.YYY+XXX.YYY  1.00  0.00    \n   
	
		note for the current version 3.2 of the PDB-format there are 80 columns.
		see e.g. http://www.wwpdb.org/documentation/format32/sect9.html

				1	   2         3         4	 5         6         7	       8\n
			12345678901234567890123456789012345678901234567890123456789012345678901234567890\n
			ATOM  ##### TYPE   - -   0    +XXX.YYY+XXX.YYY+XXX.YYY  1.00  0.00    	    EL+2\n   

		where there is in addition the Element symbol and the charge */

		for ( unsigned int j=0; j<N; ++j)
		{
			/* print ATOM and atom index */
			fprintf( inpf,"ATOM  %5u ", (j+1));

			/* print atom type */
			/* use the last 3 out of the 4 places, but adjust them on the left side as ' C2a' */
			fprintf( inpf," %-3.3s ",label[j]);     
			fprintf( inpf,"  - -   0    "); 

			/* print fractional coordinates*100 */
			/* |+XXX.YYY| -> minimum 8 signs and 3 signs after '.' */
			if (fabs(x[j])>=10.0)
			{
				if ( log_flag ) { fprintf( logfile, "\tError: _atom_site_fract_x[%d]|>10 will lead to PDB-format errors.\n",j) ; }
				XNDIFF_ERROR(81) ;
			}
			fprintf( inpf,"%8.3lf",100.0*x[j]);

			if (fabs(y[j])>=10.0)
			{
				if ( log_flag ) { fprintf( logfile, "\tError: _atom_site_fract_y[%d]|>10 will lead to PDB-format errors.\n",j) ; }
				XNDIFF_ERROR(81) ;
			}
			fprintf( inpf,"%8.3lf",100.0*y[j]);

			if (fabs(z[j])>=10.0)
			{
				if ( log_flag ) { fprintf( logfile, "\tError: _atom_site_fract_z[%d]|>10 will lead to PDB-format errors.\n",j) ; }
				XNDIFF_ERROR(81) ;
			}
			fprintf( inpf,"%8.3lf",100.0*z[j]);

			/* print rest up to EOL */
			if ( mode == 0 ) { fprintf( inpf,"  1.00  0.00    \n"); }
			else if ( mode == 1 ) { fprintf( inpf,"  1.00  0.00      -    %s\n", type[j]); }
		}
		/*terminate ATOM list by TER only if mode == 0 */
		if ( mode == 0 ) { fprintf(inpf,"TER   %5u        - -   0                                            \n", N+1); }

		/* write END */
		fprintf( inpf,"END\n");

		/* close the PDB-file */
		fclose(inpf) ;
	}


	/* read a specific column from an ascii file and return via reference a vector
	   -T can be e.g. double and int
	   -sep e.g. ' ' (includes also multiple empty spaces), ',' or ';'
	   -col, column number, col >= 1
	   -max_data == 0 -> no limits; otherwise it describes the maximum number of data rows to be read
	   -append == false, vector Tvec will be filled from the beginning on and increased in size only if necessary
	    append == true, elements will be appended to the vector Tvec
	*/
	template <class T>
	void read_column_from_ascii_file( string filename, char sep, unsigned int col, vector<T>& Tvec, bool append=false, unsigned long int max_data = 0)
	{
		FILE *inpf ;
		char sdummy[1024] ;
		char data[1024] ;
		int i ;
		char *strptr1, *strptr2 ;

		int data_type = -1 ;
		unsigned int line_index = 0 ;
		unsigned int data_index = 0 ;
		T Tval ;

		if ( typeid(Tvec) == typeid(vector<double>) ) { data_type = 0 ; }
		else if ( typeid(Tvec) == typeid(vector<float>) ) { data_type = 1 ; }
		else if ( typeid(Tvec) == typeid(vector<int>) ) { data_type = 2 ; }
		else if ( typeid(Tvec) == typeid(vector<long int>) ) { data_type = 3 ; }
		else if ( typeid(Tvec) == typeid(vector<unsigned int>) ) { data_type = 4 ; }
		else if ( typeid(Tvec) == typeid(vector<unsigned long int>) ) { data_type = 5 ; }

		/* Open file */
		if ( log_flag ) { fprintf( logfile, "\tread_column_from_ascii_file(%s,%u,%c,...)\n", filename.c_str(), col, sep) ; }

		/* get type of array buf and compare with type of the tif-file */
		if ( data_type < 0 )
		{
			fprintf( stdout, "\t\tError: Unsupported data type. Exit.\n") ;
			exit(1) ;
		}

		strcpy( sdummy, ifo) ;
		strcat( sdummy, filename.c_str()) ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) 
		{
			fprintf( stdout, "\t\tError: Could not open %s. Exit.\n", filename.c_str()) ;
			exit(1) ;
		}

		while ( fgets( sdummy, 1024, inpf) != NULL )
		{
			++line_index ;

			/* skip comments ('#') and empty lines, additional white space in front will be always ignored */
			strptr1 = sdummy ;
			while (*strptr1 == ' ') { ++strptr1 ; }

			if (*strptr1 == '#') { continue ; }
			if (*strptr1 == '\n') { continue ; }

			for ( unsigned int j=0; j<col; ++j)
			{
				/* skip all previous colons */
				if ( j < (col-1) )
				{
					if ( ( strptr2=strchr( strptr1, sep)) == NULL )
					{
						fprintf( stdout, "\t\tError: Found only %u colons in line %u. Exit.\n", j, line_index) ;
						exit(1) ;
					}
					strptr1 = strptr2 + 1 ;
					while (*strptr1 == ' ') { ++strptr1 ; }
				}
				else
				{
					/* read data from colon col */
					if ( ( strptr2=strchr( strptr1, sep)) == NULL )
					{
						if ( ( strptr2=strchr( strptr1, '\n')) == NULL )
						{
							fprintf( stdout, "\t\tError: Found only %u colons in line %u. Exit.\n", j, line_index) ;
							exit(1) ;
						}
					}
					i = (int)(strptr2-strptr1) ;
					memmove( data, strptr1, i) ;
					data[i] = 0 ;

					switch (data_type)
					{
						case 0: Tval = strtod(data,NULL) ; break ; /* double */
						case 1: Tval = strtof(data,NULL) ; break ; /* float */
						case 2: Tval = (int)strtol(data,NULL,10) ; break ; /* int */
						case 3: Tval = strtol(data,NULL,10) ; break ; /* long int */
						case 4: Tval = (unsigned int)strtoul(data,NULL,10) ; break ; /* unsigned int */
						case 5: Tval = strtoul(data,NULL,10) ; break ; /* unsigned long int */
					}

					if ( append ) { Tvec.push_back(Tval) ; }
					else
					{
						if ( Tvec.size() > data_index ) { Tvec[data_index] = Tval ; }
						else { Tvec.push_back(Tval) ; }
					}
				}
			}
			++data_index ;
			if ( max_data != 0 && data_index >= max_data ) { break ; }
			/* continue with next line */
		}

		/* close the file */
		fclose(inpf) ;

		if ( log_flag ) { fprintf( logfile, "\tdone\n") ; }
	}


	/* read the user-supplied symmetry operations from a file, symmetry operations must include symop and cenop !!
	
	   if -sym <symopfile> 
	   use traditional xyz representation
	   symop x,y+1/2,z-1/2
	   symop -x+1,-y+1/2,-z+3/2
	
	   if -sym ccp4 <ccp4config> <symopfile> 
	   use the symmetry operations for the CCP4 routine pdbset, translations must be scaled by 100 !!
	   symop x,y+100/2,z-100/2
	   symop -x+100,-y+100/2,-z+300/2
	*/
	void read_symopxyz_from_file()
	{
		FILE *inpf ;
		char sdummy[1024] ;
		char data[1024] ;
		int i ;

		char *strptr1, *strptr2 ;

		/* Open file */ 
		if ( log_flag ) { fprintf( logfile, "\tReading symmetry operations file: %s\n", symop_file) ; }

		strcpy( sdummy, ifo) ;
		strcat( sdummy, symop_file) ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(24) ; }

		num_symop = 0 ;
		while (fgets (sdummy, 1024, inpf) != NULL)
		{
			/* skip comments ('#') and empty lines, additional white space in front will be always ignored  */
			strptr1 = sdummy ;
			while (*strptr1 == ' ') 
				++strptr1 ;

			if (*strptr1 == '#')
				continue ;
			if (*strptr1 == '\n')
				continue ;


			/* CHECK for 'symop' tag */
			if ( (strptr2=strchr (strptr1, ' '))!=NULL )
			{
				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;

				if (!strcmp(data,"symop"))
				{
					++num_symop;

					strptr1=strptr2;
					while (*strptr1 == ' ') 
						++strptr1 ;

					strptr2=strchr (strptr1, '\n');
					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;

					if ( num_symop == 1 )
					{
						symopxyz=(char**) calloc( num_symop, sizeof(char*) );
						symopxyz[0]=(char*) calloc( (i+1), sizeof(char) );
					}
					else 
					{
						symopxyz=(char**) realloc ( symopxyz,  num_symop * sizeof(char*) );
						symopxyz[num_symop-1]=(char*) calloc( (i+1), sizeof(char) );
					}

					strcpy(symopxyz[num_symop-1], data);
				}
				else 
				{
					XNDIFF_ERROR(84) ;
				}
			}
		}

		/* close the file */
		fclose(inpf) ;

		if ( log_flag ) { fprintf( logfile, "\tdone\n") ; }
	}


	/* read the symop and cenop entrties from syminfo.lib */
	bool read_symopxyz_and_cenopxyz_from_syminfolib( const char* syminfolib_symmetry_keyword, char* pattern)
	{
		/* Typical entry in a cif-file :
			...
			_symmetry_cell_setting orthorhombic
			_symmetry_space_group_name_H-M 'A 21 a m'
			_symmetry_Int_Tables_number 36
			loop_
			_symmetry_equiv_pos_site_id
			_symmetry_equiv_pos_as_xyz
			1 x,y,z
			2 x,y,-z
			3 1/2+x,-y,-z
			4 1/2+x,-y,z
			5 x,1/2+y,1/2+z
			6 x,1/2+y,1/2-z
			7 1/2+x,1/2-y,1/2-z
			8 1/2+x,1/2-y,1/2+z
			...

		   Corresponding entry in syminfo.lib:

			begin_spacegroup
			number  36
			basisop z,y,-x
			symbol ccp4 0
			symbol Hall ' C 2c -2 (z,y,-x)'
			symbol xHM  'A 21 a m'
			symbol old  ''
			symbol laue '-P 2 2' 'mmm'
			symbol patt '-A 2 2' 'mmm'
			symbol pgrp ' P -2 2' 'mm2'
			hklasu ccp4 'h>=0 and k>=0 and l>=0'
			mapasu ccp4 0<=x<-1; 0<=y<-1; 0<=z<-1
			mapasu zero 0<=x<1; 0<=y<=1/4; 0<=z<=1/2
			mapasu nonz 0<=x<1; 0<=y<=1/4; 0<=z<=1/2
			cheshire 0<=x<=1/2; 0<=y<=1/2; 0<=z<=0
			symop x,y,z
			symop x+1/2,-y,-z
			symop x,y,-z
			symop x+1/2,-y,z
			cenop x,y,z
			cenop x,y+1/2,z+1/2
			end_spacegroup

		   Should be extented later for a more general treatment or better use CCP4 !
		*/

		FILE *inpf ;
		const int numsigns = 1024 ;
		char sdummy[numsigns] ;
		char data[numsigns] ;
		char *strptr1, *strptr2 ;
		int i ;
		long int lineindex = 0 ;

		unsigned int syminfolib_symmetry_keyword_length = strlen(syminfolib_symmetry_keyword) ;
		unsigned int pattern_length = strlen(pattern) ;
		bool FOUND = false ;

		/* open syminfo.lib */
		strcpy( sdummy, ifo) ;
		strcat( sdummy, "syminfo.lib") ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(75) ; }

		while ( fgets( sdummy, numsigns, inpf) != NULL )
		{
			++lineindex;

			/* ignore comments starting with '#'', empty lines, additional white space in front of them will be always ignored  */
			strptr1 = sdummy ;

			while (*strptr1 == ' ') { ++strptr1 ; }
			if (*strptr1 == '#') { continue ; }
			if (*strptr1 == '\n') { continue ; }

			if ( !strncmp( strptr1, syminfolib_symmetry_keyword, syminfolib_symmetry_keyword_length) && !FOUND )
			{
				/* found an syminfo.lib keyword like "symbol Hall" or "symbol xHM"
				   read entry and check if the entry fits the pattern 
				*/
				for ( unsigned int j=0; j<syminfolib_symmetry_keyword_length; ++j) { ++strptr1 ; }
				while (*strptr1 == ' ') { ++strptr1 ; }

				strptr2 = strchr( strptr1, '\n') ;

				i = (int)(strptr2-strptr1) ;
				memmove( data, strptr1, i) ;
				data[i] = 0 ;

				if ( !strncmp( data, pattern, pattern_length) )
				{
					/* the correct spacegroup has been found, now read symop and cenop */
					FOUND = true ;
					if ( log_flag ) { fprintf( logfile, "\t\tFound %s %s in syminfo.lib\n", syminfolib_symmetry_keyword, pattern) ; }
					continue ;
				}
				/* else { continue ; } */
			}
			/* else { continue ; } */

			if ( FOUND )
			{
				/* read all symop and cenop until end of selected space group */
				while (*strptr1 == ' ') { ++strptr1 ; }

				if ( !strncmp( strptr1, "end_spacegroup", 14) ) { break ; }

				if ( ( strptr2 = strchr( strptr1, ' ') ) != NULL )
				{
					i = (int)(strptr2-strptr1) ;
					memmove( data, strptr1, i) ;
					data[i] = 0 ;

					if ( !strcmp( data, "symop") ) 
					{
						strptr1 = strptr2 ;
						while (*strptr1 == ' ') { ++strptr1 ; }

						strptr2 = strchr( strptr1, '\n') ;

						i = (int)(strptr2-strptr1) ;
						memmove( data, strptr1, i) ;
						data[i] = 0 ;

						++num_symop;

						if ( num_symop == 1 )
						{
							symopxyz=(char**) calloc( num_symop, sizeof(char*) );
							symopxyz[0]=(char*) calloc( (i+1), sizeof(char) );
						}
						else 
						{
							symopxyz=(char**) realloc ( symopxyz,  num_symop * sizeof(char*) );
							symopxyz[num_symop-1]=(char*) calloc( (i+1), sizeof(char) );
						}

						strcpy( symopxyz[num_symop-1], data) ;
						if ( log_flag ) { fprintf( logfile, "\t\tRead symop %s from syminfo.lib\n", symopxyz[num_symop-1]) ; }
					}
					else if ( !strcmp( data, "cenop") ) 
					{
						strptr1 = strptr2 ;
						while (*strptr1 == ' ') { ++strptr1 ; }

						strptr2 = strchr( strptr1, '\n') ;

						i = (int)(strptr2-strptr1) ;
						memmove( data, strptr1, i) ;
						data[i] = 0 ;

						++num_cenop;

						if ( num_cenop == 1 )
						{
							cenopxyz=(char**) calloc( num_cenop, sizeof(char*) );
							cenopxyz[0]=(char*) calloc( (i+1), sizeof(char) );
						}
						else 
						{
							cenopxyz=(char**) realloc ( cenopxyz,  num_cenop * sizeof(char*) );
							cenopxyz[num_cenop-1]=(char*) calloc( (i+1), sizeof(char) );
						}

						strcpy( cenopxyz[num_cenop-1], data) ;
						if ( log_flag ) { fprintf( logfile, "\t\tRead cenop %s from syminfo.lib\n", cenopxyz[num_cenop-1]) ; }
					}
				}
			}
		}
		fclose(inpf) ;

		return FOUND ;
	}


	/* clears arrays of symopxyz and cenopxyz */
	void free_symopxyz_and_cenopxyz()
	{
		if ( num_symop > 0 )
		{
			for ( unsigned int i=0; i<num_symop; ++i) { free(symopxyz[i]) ; }
			free(symopxyz) ;
		}
		if ( num_cenop > 0 )
		{
			for ( unsigned int i=0; i<num_cenop; ++i) { free(cenopxyz[i]) ; }
			free(cenopxyz) ;
		}
	}


	/* Create a bash script for CCP4-routine pdbset that will be executed later on (see run_ccp4_bash_script() ).
	   First, read the configuration of the CCP4 (user-supplied paths of the CCP4 installation).
           Second, create the ccp4_bash_file script, use either the HM space group symbol from the cif-file or
	   read directly symmetry operations from a user supplied file
	*/
	void write_ccp4_bash_script()
	{
		/* example bash-file


		#!/bin/bash
		#
		source /home/mschmiele/CCP4-bin/setup-scripts/sh/ccp4.setup
		source /home/mschmiele/CCP4-bin/setup-scripts/sh/ccp4-others.setup
		#
		/home/mschmiele/CCP4-bin/ccp4-6.1.13/bin/pdbset xyzin cif.pdb xyzout pdbset.pdb <<eof-1
		symgen 'P -1'
		eof-1 


		Note that core-cif cannot be read by CCP4, primarly only PDB and mmCIF, same for output.
		mmCIF can be generated by additional commands at the end, however mmCIF isn't better than PDB

		output cif
		end

		http://www.ccp4.ac.uk/html/pdbset.html
		*/

		FILE *inpf ;
		char sdummy[1024], symmetry_space_group_name_HM[1024] ;
		char data[1024] ;
		int i, idummy ;

		char** paths_setup;
		char* path_pdbset;

		int num_paths_setup=0;

		char *strptr1, *strptr2 ;


		/* Open CCP4 configuration file */ 
		if ( log_flag ) { fprintf( logfile, "\tReading CCP4-configuration file: %s\n", ccp4_config_file) ; }

		strcpy( sdummy, ifo) ;
		strcat( sdummy, ccp4_config_file) ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(26) ; }

		while ( fgets( sdummy, 1024, inpf) != NULL)
		{
			/* skip comments ('#') and empty lines, additional white space in front will be always ignored  */
			strptr1 = sdummy ;
			while (*strptr1 == ' ') 
				++strptr1 ;

			if (*strptr1 == '#')
				continue ;
			if (*strptr1 == '\n')
				continue ;

			/* CHECK for 'setup' or 'pdbset' tag */
			if ( (strptr2=strchr (strptr1, ' '))!=NULL )
			{		
				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;

				if (!strcmp(data,"setup"))
				{
					++num_paths_setup;

					strptr1=strptr2;
					while (*strptr1 == ' ') 
						++strptr1 ;

					strptr2=strchr (strptr1, '\n');
					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;

					if (num_paths_setup==1)
					{
						paths_setup=(char**) calloc( num_paths_setup, sizeof(char*) );
						paths_setup[0]=(char*) calloc( (i+1), sizeof(char) );
					}
					else 
					{
						paths_setup=(char**) realloc ( paths_setup,  num_paths_setup * sizeof(char*) );
						paths_setup[num_paths_setup-1]=(char*) calloc( (i+1), sizeof(char) );
					}

					strncpy(paths_setup[num_paths_setup-1],data, i);
				}
				else if ( !strcmp( data, "pdbset") )
				{
					strptr1 = strptr2 ;
					while (*strptr1 == ' ') 
						++strptr1 ;

					strptr2 = strchr (strptr1, '\n') ;
					i = (int)(strptr2-strptr1) ;
					memmove( data, strptr1, i) ;
					data[i] = 0 ;

					path_pdbset=(char*) calloc( (i+1), sizeof(char) );

					strncpy( path_pdbset, data, i);
				}
			}
		}

		/* close the CCP4-configuration file */
		fclose(inpf) ;

		if ( log_flag ) { fprintf( logfile, "\tdone\n\n") ; }

		if ( symop_userdef ) { read_symopxyz_from_file() ; }

		/* write the bash-script */
		strcpy( sdummy, ofo) ;
		strcat( sdummy, ccp4_bash_file) ;
		if ( ( inpf = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(27) ; }

		fprintf( inpf, "#!/bin/bash\n") ;
		fprintf( inpf, "#\n") ;

		for (int j=0; j<num_paths_setup; ++j) { fprintf( inpf, "source %s\n",paths_setup[j]) ; }

		fprintf( inpf,"#\n") ;
		fprintf( inpf,"%s xyzin %s%s xyzout %s%s <<eof-1\n",path_pdbset, ofo, cif_pdb_file, ofo, pdb_file) ;

		/* write symmetry operations either as user-defined x,y,z form */
		/* or use the spacegroup name from cif-file in HM notation */
		if (symop_userdef)
		{
			for (unsigned int j=0; j<num_symop; ++j) { fprintf( inpf, "symgen %s\n", symopxyz[j]) ; }
		}
		else
		{
			strcpy( sdummy, "_symmetry_space_group_name_H-M") ;
			idummy = find_index_in_cfe(sdummy) ;
			if ( idummy >= 0 ) { strcpy( symmetry_space_group_name_HM, cfe[idummy].text_data[0]) ; }
			else { strcpy( symmetry_space_group_name_HM, "Unknown") ; }

			fprintf( inpf, "symgen \'%s\'\n", symmetry_space_group_name_HM) ; 
		}

		fprintf( inpf, "eof-1\n") ;

		/* close the bash-script */
		fclose(inpf) ;

		/* free path string fields */
		for (int j=0; j<num_paths_setup; ++j) { free(paths_setup[j]) ; }

		free(paths_setup) ;
		free(path_pdbset) ;
	}


	/* runs the bash script ccp4_bash_file, i.e. 
	   applies the symmetry operations according to the CCP4 routine pdbset.
	   Output file will be as specified in pdb_file */
	void run_ccp4_bash_script()
	{
		/* popen  http://www.unix.com/high-level-programming/77863-running-shell-commands-c-c.html */

		FILE *inpf ;
		char buff[1024] ;
		char sdummy[1024] ;
		
		/* popen creates a pipe so the output can be read of the program that is invoked */
		/* change the executable rights of the CCP4-bash file */

		/* print the output of popen to logfile or commandline */
		if ( log_flag ) { fprintf( logfile, "\tChanging the executable right of bash script\n") ; }

		sprintf( sdummy, "chmod u+rwx %s%s", ofo, ccp4_bash_file) ;
		if ( ( inpf = popen( sdummy, "r")) == NULL) { XNDIFF_ERROR(23); }

		while ( fgets( buff, sizeof(buff), inpf) != NULL )
		{
			if ( log_flag ) { fprintf( logfile, "%s", buff) ; }
		}

		/* close the pipe */
		pclose(inpf);

		if ( log_flag ) { fprintf( logfile, "\tdone\n\n") ; }

		/* execute bash-file and write output of pdbset to logfile or the commandline */
		if ( log_flag ) { fprintf( logfile, "\tOutput from CCP4 (pdbset):\n") ; }

		sprintf( sdummy, "./%s%s", ofo, ccp4_bash_file) ;
		if ( ( inpf = popen( sdummy, "r")) == NULL) { XNDIFF_ERROR(23) ; }

		while ( fgets( buff, sizeof(buff), inpf) != NULL ) 
		{
			if ( log_flag ) { fprintf( logfile, "\t%s", buff) ; }
		}

		/* close the file */
		fclose(inpf) ;

		if ( log_flag )
		{
			fprintf( logfile, "\t\n") ;
			fprintf( logfile, "\tdone\n") ;
		}
	}


	/* write atomic coordinates to a pcr-file with format that can be copy and pasted to a real FullProf pcr file
	   The atomic coordinates and types will be used from atom. The labels are derived from cif-file assuming that the
	   atoms read from the pdb-file are repeated blocks of those from the cif-file with symmetry operations applied on it.
	   The beta_ij values from the anisotropic atoms (Uani vs Uiso) are similarly derived from the cif-file.
	   To transform the U_ij to beta_ij see J. Appl. Cryst. (2002). 35, 477-480 :

	   http://dx.doi.org/10.1107/S0021889802008580

	   a*,b*,c* lengths of reciprocal lattice vectors in [nm^-1]
           N = diag(a*,b*,c*)
	   U* = N U_{cif} N^T 
	   beta = 2\pi^2 * U* / 100 [A^{-2}]
           ( factor 100 due to A<->nm )

	   or see 

	   http://journals.chester.iucr.org/iucr-top/cif/ddlm/DDLm_Aug08/TEST_DIC/core_struc.dic

	   ...

	   The symmetric anisotropic atomic displacement matrix beta(IJ)
	   appears in a structure factor expression as:
	
	   t = exp -       ( beta11 h h + ............ 2 beta23 k l )
	
	   It is related to the adp matrices U(IJ) and B(IJ) as follows:
	
	   t = exp -2pi**2 ( U11    h h a* a* + ...... 2 U23    k l b* c* ) 
	   t = exp - 0.25  ( B11    h h a* a* + ...... 2 B23    k l b* c* ) 

           ...

	   Note that for isotropic atoms (H,D) (Uiso flag) there will be no beta's generated.
	*/
	void write_pcr_file()
	{
		FILE *outf ;
		char sdummy[1024] ;

		double a_st, b_st, c_st ; 
		int noa_1, Uaniso ;

		a_st = betrag(b1) ;
		b_st = betrag(b2) ;
		c_st = betrag(b3) ;

		double beta[6] ;
		double Biso ;

		char atomtype[10] ;

		int index_atom_site_label = find_index_in_cfe("_atom_site_label") ;

		int index_atom_site_thermal_displace_type = find_index_in_cfe("_atom_site_thermal_displace_type", 1) ;
		if ( index_atom_site_thermal_displace_type < 0 )
		{
			/* if no Uiso or Uani entries and values are provided, use isotropic Biso = 0.0 for all atoms */
			Uaniso = 0 ; Biso = 0.0 ;

			if ( log_flag ) { fprintf( logfile, "	Using Biso = 0.0 for all atoms.\n") ; }
		}
		
		int index_atom_site_U_iso_or_equiv = find_index_in_cfe("_atom_site_U_iso_or_equiv", 1) ;

		int index_atom_site_aniso_label = find_index_in_cfe("_atom_site_aniso_label", 1) ;
		int index_atom_site_aniso_U_11 = find_index_in_cfe("_atom_site_aniso_U_11", 1) ;
		int index_atom_site_aniso_U_22 = find_index_in_cfe("_atom_site_aniso_U_22", 1) ;
		int index_atom_site_aniso_U_33 = find_index_in_cfe("_atom_site_aniso_U_33", 1) ;
		int index_atom_site_aniso_U_12 = find_index_in_cfe("_atom_site_aniso_U_12", 1) ;
		int index_atom_site_aniso_U_13 = find_index_in_cfe("_atom_site_aniso_U_13", 1) ;
		int index_atom_site_aniso_U_23 = find_index_in_cfe("_atom_site_aniso_U_23", 1) ;

		noa_1 = cfe[index_atom_site_label].data_length ;
		
		/* Open file */
		strcpy( sdummy, ofo) ;
		strcat( sdummy, pcr_file) ;
		if ( ( outf = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(34) ; }

		/* export data for the first symmetry operation */

		/*
		!Atom   Typ       X        Y        Z     Biso       Occ     In Fin N_t Spc /Codes
		!    beta11   beta22   beta33   beta12   beta13   beta23  /Codes
		C1     C       0.30610  0.40960  0.98145  0.00000   1.00000   0   0   2    0                                                                                  
				  0.00     0.00     0.00     0.00      0.00
		      0.02266  0.00532  0.00045 -0.00094  0.00073  -0.00004
			 0.00     0.00     0.00     0.00     0.00      0.00
		*/
		for ( int i=0; i<noa_1; ++i)
		{
			/* distinguish between Uiso and Uani */
			if ( index_atom_site_thermal_displace_type >= 0 )
			{
				if ( !strcmp(cfe[index_atom_site_thermal_displace_type].text_data[i],"Uani") )
				{
					Uaniso = 2 ; Biso = 0.0;
				}
				else if ( !strcmp(cfe[index_atom_site_thermal_displace_type].text_data[i],"Uiso") )
				{
					Uaniso = 0 ; Biso = 8.0 * pow(M_PI,2.0) * cfe[index_atom_site_U_iso_or_equiv].numeric_data[i];
				}
				else
				{
					/*error*/
				}
			}

			/* allow for element name only 4 signs */
			if ( !strcmp(atom[i].type,"2H") )
				strcpy( atomtype, "D") ;
			else
				strcpy( atomtype, atom[i].type) ;

			fprintf( outf, "%-7.4s%-6s %8.5lf %8.5lf %8.5lf %8.5lf  %8.5lf%4d%4d%4d%5d                                                                                  \n", cfe[index_atom_site_label].text_data[i], atomtype, atom[i].coord[0], atom[i].coord[1], atom[i].coord[2], Biso, 1.0, 0, 0, Uaniso, 0 ) ;
			fprintf( outf, "                  0.00     0.00     0.00     0.00      0.00\n" ) ;

			if ( Uaniso > 0 )
			{
				/* anisotropic beta_ij factors */

				/* search for the right element */
				for ( unsigned int j=0; j<cfe[index_atom_site_aniso_label].data_length; ++j)
				{
					if ( !strcmp(cfe[index_atom_site_aniso_label].text_data[j],cfe[index_atom_site_label].text_data[i]) )
					{
						beta[0] = cfe[index_atom_site_aniso_U_11].numeric_data[j] * a_st * a_st * 2.0 * pow(M_PI,2.0) / 100.0 ;
						beta[1] = cfe[index_atom_site_aniso_U_22].numeric_data[j] * b_st * b_st * 2.0 * pow(M_PI,2.0) / 100.0 ;
						beta[2] = cfe[index_atom_site_aniso_U_33].numeric_data[j] * c_st * c_st * 2.0 * pow(M_PI,2.0) / 100.0 ;
						beta[3] = cfe[index_atom_site_aniso_U_12].numeric_data[j] * a_st * b_st * 2.0 * pow(M_PI,2.0) / 100.0 ;
						beta[4] = cfe[index_atom_site_aniso_U_13].numeric_data[j] * a_st * c_st * 2.0 * pow(M_PI,2.0) / 100.0 ;
						beta[5] = cfe[index_atom_site_aniso_U_23].numeric_data[j] * b_st * c_st * 2.0 * pow(M_PI,2.0) / 100.0 ;
						break;
					}
				}

				fprintf( outf, "     %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf  %8.5lf\n", beta[0], beta[1], beta[2], beta[3], beta[4], beta[5] ) ;
				fprintf( outf, "         0.00     0.00     0.00     0.00     0.00      0.00\n" ) ;
			}
		}


		/* close pdb-file */
		fclose(outf) ;
	}


	/* Reads the atomic coordinates from the PDB file as specified in pdb_file and generated 
	   from CCP4 routine pdbset into the format used in the program later.
	   
	   Multiplicity, Temperature coefficients could be read from pdb however the former are 1.0 by write_pdb_file(),
	   the latter are set to 0.0 though they might be derived from the cif-file !!!
	*/
	void read_pdb_file()
	{
		/* make use of the fixed column numbers in the PDB-format specification */

		FILE *inpf ; /* FILE pointer, for the structure file */
		char sdummy[1024] ;
		char data[1024] ;

		char *strptr1, *strptr2 ;
		int i ;

		/* Open file sym_pdb_mode or sym_ccp4_mode */
		if ( sym_pdb_mode ) { strcpy( sdummy, ifo) ; }
		else { strcpy( sdummy, ofo) ; }
		strcat( sdummy, pdb_file) ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(29) ; }

		/* read atom information */
		while ( fgets( sdummy, 1024, inpf) != NULL )
		{
			/* read line by line, delete white space and read next line, if a # is detected */
			strptr1 = sdummy ;

			while (*strptr1 == ' ') 
				++strptr1 ;

			if (*strptr1 == '#')
				continue ;
	
			if ((strptr2=strchr (strptr1, ' ')) != NULL)
			{
				i = (int)(strptr2-strptr1) ;
				memmove (data, strptr1, i) ;
				data[i] = 0 ;

				if (!strcmp(data,"ATOM"))
				{
					if ( noa == 0 )
						atom = (atom_entries *) calloc( 1, sizeof(atom_entries)) ;
					else
						atom = (atom_entries *) realloc (atom, (noa+1) * sizeof(atom_entries)) ;

					/* read atom label */
					strptr1 = sdummy ;

					memmove (atom[noa].label, strptr1+12, 4) ;
					atom[noa].label[4] = 0 ;
					delallspc(atom[noa].label) ;


					/* read atomic coordinates */
					/* !!! now scale them with 0.01 to their right values !!! */
					for (int j=0; j<3; ++j)
					{
						memmove (data, strptr1+30+j*8, 8) ;
						data[8] = 0 ;
						delallspc (data) ;

						atom[noa].coord[j] = 0.01 * strtod (data,NULL) ;
					}

					/* read multiplicity */
					memmove (data, strptr1+54, 6) ;
					data[6] = 0 ;
					delallspc (data) ;

					atom[noa].mult = strtod( data, NULL) ;

					/* set temparature coefficients to 0.0 */
					/* however, they might be read from the cif file, if provided there 
					   (_atom_site_U_iso_or_equiv for isotropic and data_atom_site_aniso_U_ij for anisotropic, or
					    _atom_site_B_iso_or_equiv for isotropic and data_atom_site_aniso_B_ij for anisotropic, ... ) */
					/* the ATOM statement in the pdb file can only contain the isotropic B value or the equivalent U-isotropic U(eq)
					   for further information see PDB-format specification */
					for (int j=0; j<7; ++j)
						atom[noa].temp[j] = 0.0 ;

					/* read atom symbol/type from position 77+78 */
					/* http://www.wwpdb.org/documentation/format32/sect9.html */
					memmove (atom[noa].type, strptr1+76, 2) ;
					atom[noa].type[2] = 0 ;
					delallspc(atom[noa].type) ;

					/* replace deuterium symbol D with 2H */
					if ( !strcmp(atom[noa].type,"D") )
						strcpy(atom[noa].type,"2H") ;

					++noa;
				}
				/* else (REMARK,CRYST1, SCALE, ...) continue */
			}
		}

		/* close pdb-file */
		fclose(inpf) ;
	}


	/* write some information about unit cell, stabilizer layer and dispersion medium, 
	   used for logfile and data output in write_header_in_logfile(), write_data_in_file(), write_Y() */
	void write_common_header(FILE *outf, const char *presp, int Xn_flag)
	{
		char sdummy[1024], symmetry_space_group_name_HM[1024] ;
		int idummy ;

		strcpy( sdummy, "_symmetry_space_group_name_H-M") ;
		idummy = find_index_in_cfe(sdummy) ;
		if ( idummy >= 0 ) { strcpy( symmetry_space_group_name_HM, cfe[idummy].text_data[0]) ; }
		else { strcpy( symmetry_space_group_name_HM, "Unknown") ; }
		
		fprintf( outf, "%s********************\n", presp) ;
		fprintf( outf, "%sCrystal information:\n", presp) ;
		fprintf( outf, "%s********************\n", presp) ;
		fprintf( outf, "%sConcentration: %lf [x100 vol%%]\n", presp, conc_cry) ;
		fprintf( outf, "%sSpace group symbol [HM]: %s\n", presp, symmetry_space_group_name_HM) ;
		fprintf( outf, "%sLattice constants:     a=%10.7G [nm],    b=%10.7G [nm],     c=%10.7G [nm]\n", presp, cell_par[0], cell_par[1], cell_par[2]) ;
		fprintf( outf, "%s                   alpha=%10.7G [ °], beta=%10.7G [ °], gamma=%10.7G [ °]\n", presp, cell_par[3], cell_par[4], cell_par[5]) ;
		fprintf( outf, "%sVolume unit cell: %lf [nm^3]\n", presp, V_uc) ;
		if ( Xn_flag >= 0 )
		{
			fprintf( outf, "%sElectron density: %lf [%s]\n", presp, rho_cry, unit_rho) ;
		}
		if ( Xn_flag <= 0 )
		{
			fprintf( outf, "%sScattering length density: %lf [%s]\n", presp, sld_cry, unit_sld) ;
			fprintf( outf, "%sIncoherent cross section density: %lf [%s]\n", presp, icsd_cry, unit_icsd) ;
		}
		fprintf( outf, "%s\n", presp) ;

		fprintf( outf, "%s***********************************\n", presp) ;
		fprintf( outf, "%sInner stabilizer layer information:\n", presp) ;     
		fprintf( outf, "%s***********************************\n", presp) ;
		fprintf( outf, "%sThickness: %lf [nm]\n", presp, thickness_isl) ;
		if ( Xn_flag >= 0 )
		{
			fprintf( outf, "%sElectron density: see Table\n", presp) ;
		}
		if ( Xn_flag <= 0 )
		{
			fprintf( outf, "%sScattering length density: see Table\n", presp) ;
			fprintf( outf, "%sIncoherent cross section density: see Table\n", presp) ;
		}
		fprintf( outf, "%s\n", presp) ;

		fprintf( outf, "%s***********************************\n", presp) ;
		fprintf( outf, "%sOuter stabilizer layer information:\n", presp) ;
		fprintf( outf, "%s***********************************\n", presp) ;
		fprintf( outf, "%sThickness: %lf [nm]\n", presp, thickness_osl) ;
		if ( Xn_flag >= 0 )
		{
			fprintf( outf, "%sElectron density: see Table\n", presp) ;
		}
		if ( Xn_flag <= 0 )
		{
			fprintf( outf, "%sScattering length density: see Table\n", presp) ;
			fprintf( outf, "%sIncoherent cross section density: see Table\n", presp) ;
		}
		fprintf( outf, "%s\n", presp) ;

		fprintf( outf, "%s******************************\n", presp) ;
		fprintf( outf, "%sDispersion medium information:\n", presp) ;
		fprintf( outf, "%s******************************\n", presp) ;
		fprintf( outf, "%sConcentration: %lf [x100 vol%%]\n", presp, conc_dm) ;
		if ( Xn_flag >= 0 )
		{
			fprintf( outf, "%sElectron density: see Table\n", presp) ;
		}
		if ( Xn_flag <= 0 )
		{
			fprintf( outf, "%sScattering length density: see Table\n", presp) ;
			fprintf( outf, "%sIncoherent cross section density: see Table\n", presp) ;
		}
		fprintf( outf, "%s\n", presp) ;
	}


	/* tabulate rho, SLD and ICSD's in a user-specified file.
	   presp denotes a prefix attached before each line. Xn_flag will be 
	   typically Xn, if Xn_flag>0 only rho will be written, if Xn_flag<0
	   only SLD and ICSD will be written, for Xn_flag=0 both will be written */
	void write_rho_sld( FILE *file, const char* presp, int Xn_flag)
	{
		char sdummy[1024], sdummy2[1024] ;
		char sdummy3[128], sdummy4[128] ;
		char pattern[1024], rho_pattern[1024], sldicsd_pattern[1024] ;
		int min_num = 15 ;
		char sep[64] = "  " ;
		int rho_len[3], sldicsd_len[6] ;

		unsigned int q_max=0, q_min=0 ;

		fprintf( file, "%s\n", presp) ;
		sprintf( pattern, "%s%15s", presp, "index") ;
		if ( Xn_flag >= 0 )
		{
			for ( unsigned int i=0; i<3; ++i)
			{
				switch (i)
				{
					case 0: strcpy( sdummy3, var_name(rho_isl)) ; break;
					case 1: strcpy( sdummy3, var_name(rho_osl)) ; break;
					case 2: strcpy( sdummy3, var_name(rho_dm)) ; break;
				}
				get_var_unit_length( sdummy3, unit_rho, rho_len[i], min_num) ;
				sprintf( rho_pattern, "%%s%%%ds", rho_len[i]) ; sprintf( sdummy2, "%s [%s]", sdummy3, unit_rho) ; 
				sprintf( sdummy, rho_pattern, sep, sdummy2) ; strcat( pattern, sdummy) ;
			}
		}

		if ( Xn_flag <= 0 )
		{
			for ( unsigned int i=0; i<6; ++i)
			{
				switch (i)
				{
					case 0: strcpy( sdummy3, var_name(sld_isl)) ; break; 
					case 1: strcpy( sdummy3, var_name(sld_osl)) ; break;
					case 2: strcpy( sdummy3, var_name(sld_dm)) ; break;
					case 3: strcpy( sdummy3, var_name(icsd_isl)) ; break;
					case 4: strcpy( sdummy3, var_name(icsd_osl)) ; break;
					case 5: strcpy( sdummy3, var_name(icsd_dm)) ; break;
				}
				if ( i<3 ) { strcpy( sdummy4, unit_sld) ; }
				else { strcpy( sdummy4, unit_icsd) ; }

				get_var_unit_length( sdummy3, sdummy4, sldicsd_len[i], min_num) ;
				sprintf( sldicsd_pattern, "%%s%%%ds", sldicsd_len[i]) ; sprintf( sdummy2, "%s [%s]", sdummy3, sdummy4) ; 
				sprintf( sdummy, sldicsd_pattern, sep, sdummy2) ; strcat( pattern, sdummy) ;
			}
		}
 		fprintf( file, "%s\n", pattern) ;

		/* by default rho/sld_multi.n = 0, if otherwise set. always > 0 (single = 1, multi >= 1).
		   Xn > 0 -> rho_multi.n > 0 and sld_multi.n = 0
		   Xn < 0 -> rho_multi.n = 0 and sld_multi.n > 0
		   Xn = 0 -> rho_multi.n > 0 and sld_multi.n > 0 */

		if ( Xn_flag > 0 ) { q_min = rho_multi.n ; q_max = q_min ; }
		if ( Xn_flag < 0 ) { q_min = sld_multi.n ; q_max = q_min ; }
		if ( Xn_flag == 0 ) { q_min = min(rho_multi.n,sld_multi.n) ; q_max = max(rho_multi.n,sld_multi.n) ; }

		/* if min_q==0 the following code will have no influence ! */
		if ( Xn_flag >= 0 )
		{
			strcpy( rho_pattern, "") ;
			for ( unsigned int i=0; i<3; ++i)
			{
				sprintf( sdummy, "%s%%%d.10G", sep, rho_len[i]) ;
				strcat( rho_pattern, sdummy) ;
			}
		}
		if ( Xn_flag <= 0 )
		{
			strcpy( sldicsd_pattern, "") ;
			for ( unsigned int i=0; i<6; ++i)
			{ 
				sprintf( sdummy, "%s%%%d.10G", sep, sldicsd_len[i]) ;
				strcat( sldicsd_pattern, sdummy) ;
			}
		}

		for ( unsigned int q=0; q<q_min; ++q)
		{
			sprintf( sdummy, "%s%15d", presp, q) ;
			if ( Xn_flag >= 0 ) 
			{
				sprintf( pattern, rho_pattern, rho_isl[q], rho_osl[q], rho_dm[q]) ; 
				strcat( sdummy, pattern) ;
			}
			if ( Xn_flag <= 0 ) 
			{ 
				sprintf( pattern, sldicsd_pattern, sld_isl[q], sld_osl[q], sld_dm[q], icsd_isl[q], icsd_osl[q], icsd_dm[q]) ; 
				strcat( sdummy, pattern) ;
			}
			fprintf( file, "%s\n", sdummy) ;
		}

		/* if rho_multi.n==sld_multi.n the following code will have no influence ! */
		if ( Xn_flag >= 0 )
		{
			strcpy( rho_pattern, "") ;
			for ( unsigned int i=0; i<3; ++i)
			{
				if ( rho_multi.n >= sld_multi.n ) 
				{
					sprintf( sdummy, "%s%%%d.10G", sep, rho_len[i]) ;
				}
				else
				{
					sprintf( sdummy, "%s%%%ds", sep, rho_len[i]) ;
				}
				strcat( rho_pattern, sdummy) ;
			}
			
		}
		if ( Xn_flag <= 0 )
		{
			strcpy( sldicsd_pattern, "") ;
			for ( unsigned int i=0; i<6; ++i)
			{
				if ( rho_multi.n <= sld_multi.n ) 
				{
					sprintf( sdummy, "%s%%%d.10G", sep, sldicsd_len[i]) ;
				}
				else
				{
					sprintf( sdummy, "%s%%%ds", sep, sldicsd_len[i]) ;
				}
				strcat( sldicsd_pattern, sdummy) ;
			}
		}

		for ( unsigned int q=q_min; q<q_max; ++q)
		{
			sprintf( sdummy, "%s%15d", presp, q) ;
			if ( Xn_flag >= 0 )
			{
				if ( rho_multi.n >= sld_multi.n ) 
				{
					sprintf( pattern, rho_pattern, rho_isl[q], rho_osl[q], rho_dm[q]) ;
				}
				else
				{
					sprintf( pattern, rho_pattern, " ", " ", " ") ;
				}
				strcat( sdummy, pattern) ;
			}
			if ( Xn_flag <= 0 )
			{
				if ( rho_multi.n <= sld_multi.n ) 
				{
					sprintf( pattern, sldicsd_pattern, sld_isl[q], sld_osl[q], sld_dm[q], icsd_isl[q], icsd_osl[q], icsd_dm[q]) ;
				}
				else
				{
					sprintf( pattern, sldicsd_pattern, " ", " ", " ", " ", " ", " ") ;
				}
				strcat( sdummy, pattern) ;
			}
			fprintf( file, "%s\n", sdummy) ;
		}

		fprintf( file, "%s\n%s\n", presp, presp) ;
	}


	/* write an abstract of important structural data from the par+cif/pdb-file to the logfile */
	void write_header_in_logfile( string presp )
	{
		fprintf( logfile, "%sFractional atomic coordinates:\n", presp.c_str()) ;

		fprintf( logfile, "%s   %10s %5s %5s %10s %10s %10s", presp.c_str(), "Label", "Type", "Mult", "x", "y", "z") ;
		if ( Xn >= 0 ) { fprintf( logfile, " %10s", "Z=f0(s=0)") ; }
		if ( Xn <= 0 ) { fprintf( logfile, " %21s %10s", "b_coh [fm] (Re, Im)", "ics [barn]") ; }
		fprintf( logfile, " %7s %7s %7s %7s %7s %7s %7s\n", "Uiso", "Uani_11", "Uani_22", "Uani_33", "Uani_12", "Uani_13", "Uani_23") ;

		for ( unsigned int j=0; j<noa; ++j)
		{
			fprintf( logfile, "%s   %10s %5s %5.2G %10.7G %10.7G %10.7G", presp.c_str(), atom[j].label, atom[j].type, atom[j].mult, atom[j].coord[0], atom[j].coord[1], atom[j].coord[2]) ; 
			if ( Xn >= 0 ) { fprintf( logfile, " %10.7G", get_atomic_scattering_factor( atom[j].type, 0.0 ) ) ; }
			if ( Xn <= 0 ) { fprintf( logfile, " %10.7G %10.7G %10.7G", real(atom[j].coh_sc_length), imag(atom[j].coh_sc_length), atom[j].incoh_sc_cross ) ; }
			fprintf( logfile, " %7.4G %7.4G %7.4G %7.4G %7.4G %7.4G %7.4G\n", atom[j].temp[0], atom[j].temp[1], atom[j].temp[2], atom[j].temp[3], atom[j].temp[4], atom[j].temp[5], atom[j].temp[6]) ;
		}

		fprintf( logfile, "%s\n", presp.c_str()) ;

		write_common_header( logfile, presp.c_str(), Xn) ;

		write_rho_sld( logfile, presp.c_str(), Xn) ;
	}


	/* returns by reference the string length of variable name + its unit */
	void get_var_unit_length( const char *var, const char* unit, int& len, int min_len=0)
	{
		char sdummy[1024] ;
		sprintf ( sdummy, "%s [%s]", var, unit) ; 
		len = strlen(sdummy) ;

		if ( len < min_len ) { len = min_len ; }
	}


	/* lists call command and time stamps in a file */
	void write_call_arguments( FILE *outf, const char *presp)
	{
		fprintf( outf, "%sXNDiff --- A calculation program for X-ray and neutron powder diffraction patterns\n", presp) ;
		fprintf( outf, "%sVersion %s    %s by %s \n", presp, XNDIFF_VERSION, XNDIFF_DATE, XNDIFF_AUTHOR) ;
		fprintf( outf, "%s----------------------------------------------------------------------------------\n", presp) ;
		fprintf( outf, "%sProgram called with:\n", presp) ;
		fprintf( outf, "%s", presp) ;
		/* Export function call parameters ecarg */
		for ( int i=0; i<*ecarg; ++i) { fprintf( outf, "%s ", evarg[i]) ; }
		fprintf( outf, "\n") ;
		fprintf( outf, "%s\n", presp) ;
		fprintf( outf, "%s\n", presp) ;
		
		fprintf( outf, "%sProgram started at %s", presp, startingtime ) ;
		if ( time_flag > 0 )
		{
			time_t now ;
			struct tm * timeinfo ;

			time(&now) ;
			timeinfo = localtime(&now) ;

			fprintf( outf, "%sStart writing this file at %s", presp, asctime(timeinfo) ) ;
		}
		fprintf( outf, "%s\n", presp) ;
	}


	/* generate output files with (d)S_X/n and (d)B_n */
	/* for (d)B_n use MULT_INC = MULT/(nav_pc*nav_ac) instead of MULT as for (d)S_n
	   because for the incoherent scattering no orientational average was performed */
	void write_data_in_file( char *oname, const char *presp, char X_or_n, int C_REP, int MULT, char *MULT_str, double SUM_WEIGHT, double V_CRY, double V_ISL, double V_OSL, double V_DM, double V_IRR, const char *S_name, double **S, const char *B_name=NULL, double *B=NULL )
	{
		int min_num = 15 ;
		char sep[64] = "  " ;
		int indent = 0 ;
		int nq ;
		double cs ;

		int S_len = 0 ;
		int B_len = 0 ;

		int idummy = 0 ;

		unsigned int MULT_INC ;
		double SUM_WEIGHT_CORR ;

		char sdummy[1024] ;
		char pattern[1024] ;
		char S_str[256] ;
		char B_str[256] ;

		FILE *outf;

		strcpy( sdummy, ofo) ;
		strcat( sdummy, oname) ;
		if ( ( outf = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(19) ; }

		write_call_arguments( outf, presp) ;

		/* write some information about unit cell, stabilizer layer and dispersion medium */
		if ( X_or_n == 'X') { write_common_header( outf, presp, 1) ; }
		if ( X_or_n == 'n') { write_common_header( outf, presp, -1) ; }

		/* write (total) volumes of V_CRY, V_OSL, V_ISL, V_DM and V_IRR */
		fprintf( outf, "%s%s = %-15.10G [%s]\n", presp, var_name(V_CRY), V_CRY, unit_V_CRY) ;
		fprintf( outf, "%s%s = %-15.10G [%s]\n", presp, var_name(V_ISL), V_ISL, unit_V_ISL) ;
		fprintf( outf, "%s%s = %-15.10G [%s]\n", presp, var_name(V_OSL), V_OSL, unit_V_OSL) ;
		fprintf( outf, "%s%s = %-15.10G [%s]\n", presp, var_name(V_DM), V_DM, unit_V_DM) ;
		fprintf( outf, "%s%s = %-15.10G [%s]\n", presp, var_name(V_IRR), V_IRR, unit_V_IRR) ;

		/* write MULTIPLICITY and (normalized) SUM_WEIGHT_FACTOR information */
		fprintf( outf, "%s", MULT_str) ;
		fprintf( outf, "%sC_REP + 1 = %d\n", presp, C_REP + 1 ) ;
		fprintf( outf, "%sSUM_WEIGHT_FACTOR = %15.10G\n", presp, SUM_WEIGHT) ;
		SUM_WEIGHT_CORR = SUM_WEIGHT / (double)MULT / (double)( C_REP + 1 ) ;
		fprintf( outf, "%sSUM_WEIGHT_CORR = %15.10G\n", presp, SUM_WEIGHT_CORR ) ;
		if ( X_or_n == 'n' ) 
		{
			MULT_INC = 1 ;
			fprintf( outf, "%sINCOHERENT MULTIPLICITY = %d\n", presp, MULT_INC) ;
			fprintf( outf, "%s\n", presp) ;
		}

		indent = min_num ;
		if ( X_or_n == 'X' )
		{
			nq = rho_multi.n ;

			get_var_unit_length( var_name(rho_isl), unit_rho, idummy) ;
			if ( idummy > indent ) { indent = idummy ; }
			get_var_unit_length( var_name(rho_osl), unit_rho, idummy) ;
			if ( idummy > indent ) { indent = idummy ; }
			get_var_unit_length( var_name(rho_dm), unit_rho, idummy) ;
			if ( idummy > indent ) { indent = idummy ; }

			sprintf ( S_str, "%s [%s]", S_name, unit_S_X) ;
			S_len = strlen(S_str) ;
			if ( S_len < min_num )
				S_len = min_num ;

		}
		if ( X_or_n == 'n' )
		{
			if ( ( B_name == NULL ) || ( B == NULL ) )
				XNDIFF_ERROR( 40 ) ;

			nq = sld_multi.n ;

			get_var_unit_length( var_name(sld_isl), unit_sld, idummy) ;
			if ( idummy > indent ) { indent = idummy ; }
			get_var_unit_length( var_name(sld_osl), unit_sld, idummy) ;
			if ( idummy > indent ) { indent = idummy ; }
			get_var_unit_length( var_name(sld_dm), unit_sld, idummy) ;
			if ( idummy > indent ) { indent = idummy ; }

			get_var_unit_length( var_name(icsd_isl), unit_icsd, idummy) ;
			if ( idummy > indent ) { indent = idummy ; }
			get_var_unit_length( var_name(icsd_osl), unit_icsd, idummy) ;
			if ( idummy > indent ) { indent = idummy ; }
			get_var_unit_length( var_name(icsd_dm), unit_icsd, idummy) ;
			if ( idummy > indent ) { indent = idummy ; }

			sprintf ( S_str, "%s [%s]", S_name, unit_S_n) ;
			S_len = strlen(S_str) ;
			if ( S_len < min_num )
				S_len = min_num ;

			sprintf ( B_str, "%s [%s]", B_name, unit_S_n) ;
			B_len = strlen(B_str) ;
			if ( B_len < min_num )
				B_len = min_num ;

		}
		sprintf( pattern, "%s%%%ds", presp, indent) ;
		fprintf( outf, pattern, "index") ;

		if ( X_or_n == 'X' ) 
		{
			/* write indices q */
			sprintf( pattern, "%s%%%dd", sep, S_len) ; 
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, q) ;

			fprintf( outf, "\n#\n") ;

			/* write rho's */
			sprintf( pattern, "%s%%%ds", presp, indent) ; 
			sprintf( sdummy, "%s [%s]", var_name(rho_isl), unit_rho) ;
			fprintf( outf, pattern, sdummy ) ;
			sprintf( pattern, "%s%%%d.10G", sep, S_len) ; 
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, rho_isl[q]) ;

			fprintf( outf, "\n") ;

			sprintf( pattern, "%s%%%ds", presp, indent) ; 
			sprintf( sdummy, "%s [%s]", var_name(rho_osl), unit_rho) ;
			fprintf( outf, pattern, sdummy ) ;
			sprintf( pattern, "%s%%%d.10G", sep, S_len) ; 
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, rho_osl[q]) ;

			fprintf( outf, "\n") ;

			sprintf( pattern, "%s%%%ds", presp, indent) ; 
			sprintf( sdummy, "%s [%s]", var_name(rho_dm), unit_rho) ;
			fprintf( outf, pattern, sdummy ) ;
			sprintf( pattern, "%s%%%d.10G", sep, S_len) ; 
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, rho_dm[q]) ;

			fprintf( outf, "\n") ;

			/* write column variables names */
			sprintf( pattern, "%s%%%ds", presp, indent) ;
			fprintf( outf, pattern, "s [1/nm]") ;
			sprintf( pattern, "%s%%%ds", sep, S_len) ; 
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, S_str) ;

			fprintf( outf, "\n") ;

			/* write data */
			indent += strlen(presp) ;

			cs = s1 ;
			for ( unsigned int i=0; i<par->np; ++i)
			{
				cs += ds ;
				sprintf( pattern, "%%%d.10G", indent) ; 
				fprintf( outf, pattern, cs) ;
				sprintf( pattern, "%s%%%d.10G", sep, S_len) ; 
				for ( int q=0; q<nq; ++q)
					fprintf( outf, pattern, S[q][i] / ( (double)MULT * SUM_WEIGHT_CORR * V_IRR ) ) ; /* [1/cm] */
	
				fprintf( outf, "\n") ;
			}
		}
		if ( X_or_n == 'n' ) 
		{
			/* write indices q */
		 	sprintf( pattern, "%s%%%ds%s%%%dd", sep, S_len, sep, B_len) ; 
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, " ", q) ;

			fprintf( outf, "\n#\n") ;

			/* write SLD's and ICSD's */
			sprintf( pattern, "%s%%%ds", presp, indent) ; 
			sprintf( sdummy, "%s [%s]", var_name(sld_isl), unit_sld) ;
			fprintf( outf, pattern, sdummy ) ;
			sprintf( pattern, "%s%%%ds%s%%%d.10G", sep, S_len, sep, B_len) ;
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, " ",  sld_isl[q]) ;

			fprintf( outf, "\n") ;

			sprintf( pattern, "%s%%%ds", presp, indent) ; 
			sprintf( sdummy, "%s [%s]", var_name(sld_osl), unit_sld) ;
			fprintf( outf, pattern, sdummy ) ;
			sprintf( pattern, "%s%%%ds%s%%%d.10G", sep, S_len, sep, B_len) ;
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, " ", sld_osl[q]) ;

			fprintf( outf, "\n") ;

			sprintf( pattern, "%s%%%ds", presp, indent) ; 
			sprintf( sdummy, "%s [%s]", var_name(sld_dm), unit_sld) ;
			fprintf( outf, pattern, sdummy ) ;
			sprintf( pattern, "%s%%%ds%s%%%d.10G", sep, S_len, sep, B_len) ;
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, " ", sld_dm[q]) ;

			fprintf( outf, "\n") ;

			sprintf( pattern, "%s%%%ds", presp, indent) ; 
			sprintf( sdummy, "%s [%s]", var_name(icsd_isl), unit_icsd) ;
			fprintf( outf, pattern, sdummy ) ;
			sprintf( pattern, "%s%%%ds%s%%%d.10G", sep, S_len, sep, B_len) ;
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, " ", icsd_isl[q]) ;

			fprintf( outf, "\n") ;

			sprintf( pattern, "%s%%%ds", presp, indent) ; 
			sprintf( sdummy, "%s [%s]", var_name(icsd_osl), unit_icsd) ;
			fprintf( outf, pattern, sdummy ) ;
			sprintf( pattern, "%s%%%ds%s%%%d.10G", sep, S_len, sep, B_len) ;
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, " ", icsd_osl[q]) ;

			fprintf( outf, "\n") ;

			sprintf( pattern, "%s%%%ds", presp, indent) ; 
			sprintf( sdummy, "%s [%s]", var_name(icsd_dm), unit_icsd) ;
			fprintf( outf, pattern, sdummy ) ;
			sprintf( pattern, "%s%%%ds%s%%%d.10G", sep, S_len, sep, B_len) ;
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, " ", icsd_dm[q]) ;

			fprintf( outf, "\n") ;

			/* write column variables names */
			sprintf( pattern, "%s%%%ds", presp, indent) ; 
			fprintf( outf, pattern, "s [1/nm]") ;
			sprintf( pattern, "%s%%%ds%s%%%ds", sep, S_len, sep, B_len) ;
			for ( int q=0; q<nq; ++q)
				fprintf( outf, pattern, S_str, B_str) ;

			fprintf( outf, "\n") ;

			/* write data */
			indent += strlen(presp) ;

			cs = s1 ;
			for ( unsigned int i=0; i<par->np; ++i)
			{
				cs += ds ;
				sprintf( pattern, "%%%d.10G", indent) ;
				fprintf( outf, pattern, cs) ;
				sprintf( pattern, "%s%%%d.10G%s%%%d.10G", sep, S_len, sep, B_len) ;
				for ( int q=0; q<nq; ++q)
					fprintf( outf, pattern, S[q][i] / ( (double)MULT * SUM_WEIGHT_CORR * V_IRR ), B[q] / ( (double)MULT_INC * V_IRR ) ) ; /* [1/cm] */
	
				fprintf( outf, "\n") ;
			}
		}
		fclose(outf) ;
	}


	/* generate output files with Yc(i)_X/n */
	/* for (d)Yi_n use MULT_INC = MULT/(nav_pc*nav_ac) instead of MULT as for (d)Yc_n
	   because for the incoherent scattering no orientational average was performed */
	void write_Y( char *oname, const char *presp, char X_or_n, int C_REP, unsigned int MULT, char *MULT_str, double SUM_WEIGHT, double V_CRY, double V_ISL, double V_OSL, double V_DM, double V_IRR, const char *Yc_name, double **Yc, const char *Yi_name=NULL, double *Yi=NULL )
	{
		int min_num = 15 ;
		char sep[64] = "  " ;
		double cs ;

		int Y_len ;

		int indent ;

		int idummy ;
		
		unsigned int MULT_INC ;
		double SUM_WEIGHT_CORR ;

		char sdummy[1024] ;
		char pattern[1024] ;
		char Y_str[256] ;
		const char s_str[256] = "s [1/nm]" ;

		FILE *outf;

		strcpy( sdummy, ofo) ;
		strcat( sdummy, oname) ;
		if ( ( outf = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(19) ; }

		write_call_arguments( outf, presp) ;

		/* write some information about unit cell, stabilizer layer and dispersion medium */
		if ( X_or_n == 'X') { write_common_header( outf, presp, 1) ; }
		if ( X_or_n == 'n') { write_common_header( outf, presp, -1) ; }

		/* write (total) volumes of V_CRY, V_OSL, V_ISL, V_DM and V_IRR */
		fprintf( outf, "%s%s = %-15.10G [%s]\n", presp, var_name(V_CRY), V_CRY, unit_V_CRY) ;
		fprintf( outf, "%s%s = %-15.10G [%s]\n", presp, var_name(V_ISL), V_ISL, unit_V_ISL) ;
		fprintf( outf, "%s%s = %-15.10G [%s]\n", presp, var_name(V_OSL), V_OSL, unit_V_OSL) ;
		fprintf( outf, "%s%s = %-15.10G [%s]\n", presp, var_name(V_DM), V_DM, unit_V_DM) ;
		fprintf( outf, "%s%s = %-15.10G [%s]\n", presp, var_name(V_IRR), V_IRR, unit_V_IRR) ;

		/* write MULTIPLICITY and (normalized) SUM_WEIGHT_FACTOR information */
		fprintf( outf, "%s", MULT_str) ;
		fprintf( outf, "%sC_REP + 1 = %d\n", presp, C_REP + 1 ) ;		
		fprintf( outf, "%sSUM_WEIGHT_FACTOR = %15.10G\n", presp, SUM_WEIGHT) ;
		SUM_WEIGHT_CORR = SUM_WEIGHT / (double)MULT / (double)( C_REP + 1 ) ;
		fprintf( outf, "%sSUM_WEIGHT_CORR = %15.10G\n", presp, SUM_WEIGHT_CORR ) ;
		if ( X_or_n == 'n' ) 
		{
			MULT_INC = 1 ;
			fprintf( outf, "%sINCOHERENT MULTIPLICITY = %d\n", presp, MULT_INC) ;
			fprintf( outf, "%s\n", presp) ;
		}

		Y_len = min_num ;
		if ( X_or_n == 'X' )
		{
			for ( int q=0; q<10; ++q)
			{
				sprintf( Y_str, "%s[%d] [%s]", Yc_name, q, unit_Yc_X[q]) ;
				idummy = strlen(Y_str) ;
				if ( idummy > Y_len ) { Y_len = idummy ; }
			}
		}
		if ( X_or_n == 'n' )
		{
			if ( ( Yi_name == NULL ) || ( Yi == NULL ) ) { XNDIFF_ERROR( 42 ) ; }

			for ( int q=0; q<10; ++q)
			{
				sprintf( Y_str, "%s[%d] [%s]", Yc_name, q, unit_Yc_n[q]) ;
				idummy = strlen(Y_str) ;
				if ( idummy > Y_len ) { Y_len = idummy ; }
			}

			for ( int q=0; q<4; ++q)
			{
				sprintf( Y_str, "%s[%d] [%s]", Yi_name, q, unit_Yi_n[q]) ;
				idummy = strlen(Y_str) ;
				if ( idummy > Y_len ) { Y_len = idummy ; }
			}

		}

		if ( X_or_n == 'X' ) 
		{
			/* write column variables names */
			indent = strlen(s_str) ;
			if ( indent < min_num ) { indent = min_num ; }
			sprintf( pattern, "%s%%%ds", presp, indent) ;
			fprintf( outf, pattern, s_str) ;

			sprintf( pattern, "%s%%%ds", sep, Y_len) ;
			for ( int q=0; q<10; ++q)
			{
				sprintf( Y_str, "%s[%d] [%s]", Yc_name, q, unit_Yc_X[q]) ;
				fprintf( outf, pattern, Y_str) ;
			}

			fprintf( outf, "\n") ;

			/* write data */
			indent += strlen(presp) ;

			cs = s1 ;
			for ( unsigned int i=0; i<par->np; ++i)
			{
				cs += ds ;
				sprintf( pattern, "%%%d.10G", indent) ;
				fprintf( outf, pattern, cs) ;
				sprintf( pattern, "%s%%%d.10G", sep, Y_len) ;
				for ( int q=0; q<10; ++q)
					fprintf( outf, pattern, Yc[q][i] / ( (double)MULT * SUM_WEIGHT_CORR * V_IRR ) ) ; /* [ unit_Yc_X[q] ] */

				fprintf( outf, "\n") ;
			}
		}
		if ( X_or_n == 'n' ) 
		{
			/* write column variables names */
			indent = strlen(s_str) ;
			if ( indent < min_num ) { indent = min_num ; }
			sprintf( pattern, "%s%%%ds", presp, indent) ;
			fprintf( outf, pattern, s_str) ;

			sprintf( pattern, "%s%%%ds", sep, Y_len) ;
			for ( int q=0; q<10; ++q)
			{
				sprintf( Y_str, "%s[%d] [%s]", Yc_name, q, unit_Yc_n[q]) ;
				fprintf( outf, pattern, Y_str) ;
			}
			for ( int q=0; q<4; ++q)
			{
				sprintf( Y_str, "%s[%d] [%s]", Yi_name, q, unit_Yi_n[q]) ;
				fprintf( outf, pattern, Y_str) ;
			}

			fprintf( outf, "\n") ;

			/* write data */
			indent += strlen(presp) ;

			cs = s1 ;
			for ( unsigned int i=0; i<par->np; ++i)
			{
				cs += ds ;
				sprintf( pattern, "%%%d.10G", indent) ;
				fprintf( outf, pattern, cs) ;
				sprintf( pattern, "%s%%%d.10G", sep, Y_len) ;
				for ( int q=0; q<10; ++q)
					fprintf( outf, pattern, Yc[q][i] / ( (double)MULT * SUM_WEIGHT_CORR * V_IRR ) ) ; /* [ unit_Yc_n[q] ] */

				for ( int q=0; q<4; ++q)
					fprintf( outf, pattern, Yi[q] / ( (double)MULT_INC * V_IRR ) ) ; /* [ unit_Yi_n[q] ] */

				fprintf( outf, "\n") ;
			}
		}
		fclose(outf) ;
	}


	void read_column_from_file( char *file, unsigned int col, unsigned int& n, double **arr)
	{
		FILE *inpf ; 
		char sdummy[1024] ;
		char data[1024] ;
		unsigned int counter = 0 ;
		unsigned int curr_col ;
		bool bdummy = false ;

		char *strptr1, *strptr2 ;
		int i, j ;

		/* Open multi-file */
		strcpy( sdummy, ifo) ;
		strcat( sdummy, file) ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(38) ; }

		while ( fgets( sdummy, 1024, inpf) != NULL )
		{
			strptr1 = sdummy ;
			/* read line by line, delete white space and read next line, if a # is detected */
			while ( *strptr1 == ' ' ) 
				++strptr1 ;

			if ( *strptr1 == '#' || *strptr1 == '\n' )
				continue ;

			strptr1 = sdummy ;
			curr_col = 1 ;
			while ( curr_col < col && ( strptr2 = strchr( strptr1, ',') ) != NULL )
			{
				strptr1 = strptr2 ;
				/* skip comma in order to check for the next one */
				++strptr1 ;

				++curr_col ;
			}
			if ( curr_col != col )
				XNDIFF_ERROR(39) ;

			while (*strptr1 == ' ') 
				++strptr1 ;

			if ( ( strptr2 = strchr( strptr1, ',') ) == NULL )
				strptr2 = strchr( strptr1, '\n') ;

			i = (int)(strptr2-strptr1) ;
			/* test if there's still data in the column (i>0) or not i==0.
			   In the latter case return to read_rho_and_sld. */
			if ( i == 0 )
				break ;

			memmove (data, strptr1, i) ;
			data[i] = 0 ;

			delallspc(data) ;

			/* Furthermore check for valid numerical data, otherwise return */
			j = 0 ;
			while ( data[j] != 0 )
			{
				if ( isdigit(data[j]) ) { ++j ; continue ; }
				if ( data[j] == '.' ) { ++j ; continue ; }
				if ( ( data[j] == '-' ) || ( data[j] == '+' ) ) { ++j ; continue ; }
				if ( ( data[j] == 'e' ) || ( data[j] == 'E' ) ) { ++j ; continue ; }

				bdummy = true ;
				break ;
			}
			if ( bdummy ) { break ; }

			if ( counter == 0 )
				*arr = (double *) calloc( 1, sizeof(double)) ;
			else
				*arr = (double *) realloc ( *arr, ( counter + 1 ) * sizeof(double)) ;

			(*arr)[counter] = strtod( data, NULL) ;
			++counter;
		}
		/* close multi-file */
		fclose(inpf) ;

		n = counter ;
	}


	void read_rho_and_sld( char *str, char **file, unsigned int& col, double **arr, char X_or_n)
	{
		char *strptr1, *strptr2 ;
		char data[1024] ;
		int i ;
		unsigned int n = 0 ;

		strptr1 = str ;
		strptr2 = strchr( strptr1, ' ') ;
		if ( strptr2 != NULL )
		{
			i = (int)(strptr2-strptr1) ;
			memmove (data, strptr1, i) ;
			data[i] = 0 ;
 			*file = (char *) calloc( i+1, sizeof (char)) ;
			strcpy( *file, data) ;

			strptr1 = strptr2 ;
			while (*strptr1 == ' ') 
				++strptr1 ;

			strptr2 = strchr( strptr1, '\n') ;
			i = (int)(strptr2-strptr1) ;
			memmove (data, strptr1, i) ;
			data[i] = 0 ;
			col = (unsigned int)strtol( data, NULL, 10) ;

			/* set rho/sld_multi = true */
			if ( X_or_n == 'X' )
			{
				rho_multi.multiq = true ;
				read_column_from_file( *file, col, n, arr) ;

				if ( rho_multi.n == 0 )
				{
					rho_multi.n = n ;
				}
				else
				{
					if ( n != rho_multi.n ) { XNDIFF_ERROR(37) ; }
				}
			}
			else if ( X_or_n == 'n' )
			{
				sld_multi.multiq = true ;
				read_column_from_file( *file, col, n, arr) ;

				if ( sld_multi.n == 0 )
				{
					sld_multi.n = n ;
				}
				else
				{
					if ( n != sld_multi.n ) { XNDIFF_ERROR(37) ; }
				}
			}
		}
		else
		{
			/* check by rho/sld_multi.n > 1 that no mixing occurs between multi and single.
			   allow rho/sld_multi.n = 1 with rho/sld_multi.multiq = true */
			if ( X_or_n == 'X' )
			{
				if ( rho_multi.n == 0 )
				{
					rho_multi.n = 1 ;
				}
				else
				{
					if ( rho_multi.n != 1 )
						XNDIFF_ERROR(37) ;
				}
			}
			else if ( X_or_n == 'n' )
			{
				if ( sld_multi.n == 0 )
				{
					sld_multi.n = 1 ;
				}
				else
				{
					if ( sld_multi.n != 1 )
						XNDIFF_ERROR(37) ;
				}
			}

			i = strlen(par_file) ;
 			*file = (char *) calloc( i+1, sizeof (char)) ;
			strcpy( *file, par_file) ;

			strptr2 = strchr( strptr1, '\n') ;
			i = (int)(strptr2-strptr1) ;
			memmove (data, strptr1, i) ;
			data[i] = 0 ;
			*arr = (double *) calloc( 1, sizeof(double)) ;
			(*arr)[0] = strtod( data, NULL) ;
		}
	}


	/* +par */
	void read_par_file()
	{
		FILE *inpf ; /* FILE pointer, for the structure file */	
		char sdummy[1024] ;
		char item[1024] ;
		char data[1024] ;
		int lineindex = 0 ;

		char *strptr0, *strptr1, *strptr2 ;
		int i ;

		/* Open par-file */
		if ( par_num_def ) { strcpy( sdummy, ofo) ; }
		else { strcpy( sdummy, ifo) ; }
		strcat( sdummy, par_file) ;
		if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(8) ; }

		/* initialize some necessary rho/sld_multi values */
		rho_multi.multiq = false ;
		rho_multi.n = 0 ;
		sld_multi.multiq = false ;
		sld_multi.n = 0 ;


		bool rho_isdef[3] = { false, false , false} ;
		bool sld_isdef[3] = { false, false , false} ;
		bool icsd_isdef[3] = { false, false , false} ;
		bool thickness_isdef[2] = { false, false} ;
		bool conc_isdef[2] = { false, false} ;

		while (fgets (sdummy, 1024, inpf) != NULL)
		{
			++lineindex;

			/* read line by line, delete white space and read next line, if a # is detected */
			strptr0 = sdummy ;
			strptr1 = sdummy ;

			while (*strptr1 == ' ') { ++strptr1 ; }

			if (*strptr1 == '#') { continue ; }
	
			if ((strptr2=strchr(strptr1, ' ')) != NULL)
			{
				i = (int)(strptr2-strptr1) ;
				memmove (item, strptr1, i) ;
				item[i] = 0 ;

				strptr1 = strptr2 ;
				while (*strptr1 == ' ') { ++strptr1 ; }


				if ( !strcmp( item, par_keyword.THICKNESS_INLAY) )
				{
					/* thickness [nm] of the inner stabilizer layer */
					strptr2 = strchr (strptr1, '\n') ;
					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;

					thickness_isl = strtod (data,NULL) ;
					thickness_isdef[0] = true ;
				}
				else if ( !strcmp( item, par_keyword.THICKNESS_OUTLAY) )
				{
					/* thickness [nm] of the outer stabilizer layer */
					strptr2 = strchr (strptr1, '\n') ;
					i = (int)(strptr2-strptr1) ;
					memmove (data, strptr1, i) ;
					data[i] = 0 ;

					thickness_osl = strtod (data,NULL) ;
					thickness_isdef[1] = true ;
				}
				else if ( ( !strcmp( item, par_keyword.XRAY_RHO_INLAY) ) && ( Xn >= 0 ) )
				{
					/* electron density [e-/nm^3] of the inner stabilizer layer */
 					read_rho_and_sld( &sdummy[(int)(strptr1-strptr0)], &rho_multi.XRAY_RHO_INLAY_file, rho_multi.XRAY_RHO_INLAY_col, &rho_isl, 'X') ;
					rho_isdef[0] = true ;
				}
				else if ( ( !strcmp( item, par_keyword.XRAY_RHO_OUTLAY) ) && ( Xn >= 0 ) )
				{
					/* electron density [e-/nm^3] of the outer stabilizer layer */
 					read_rho_and_sld( &sdummy[(int)(strptr1-strptr0)], &rho_multi.XRAY_RHO_OUTLAY_file, rho_multi.XRAY_RHO_OUTLAY_col, &rho_osl, 'X') ;
					rho_isdef[1] = true ;
				}
				else if ( ( !strcmp( item, par_keyword.XRAY_RHO_DM) ) && ( Xn >= 0 ) )
				{
					/* electron density of dispersion medium [e-/nm^3] */
					read_rho_and_sld( &sdummy[(int)(strptr1-strptr0)], &rho_multi.XRAY_RHO_DM_file, rho_multi.XRAY_RHO_DM_col, &rho_dm, 'X') ;
					rho_isdef[2] = true ;
				}
				else if ( ( !strcmp( item, par_keyword.NEUT_SLD_INLAY) ) && ( Xn <= 0 ) )
				{
					/* neutron SLD of the inner stabilizer layer [1e-6 A^(-2)] */
 					read_rho_and_sld( &sdummy[(int)(strptr1-strptr0)], &sld_multi.NEUT_SLD_INLAY_file, sld_multi.NEUT_SLD_INLAY_col, &sld_isl, 'n') ;
					sld_isdef[0] = true ;
				}
				else if ( ( !strcmp( item, par_keyword.NEUT_SLD_OUTLAY) ) && ( Xn <= 0 ) )
				{
					/* neutron SLD of the outer stabilizer layer [1e-6 A^(-2)] */
 					read_rho_and_sld( &sdummy[(int)(strptr1-strptr0)], &sld_multi.NEUT_SLD_OUTLAY_file, sld_multi.NEUT_SLD_OUTLAY_col, &sld_osl, 'n') ;
					sld_isdef[1] = true ;
				}
				else if ( ( !strcmp( item, par_keyword.NEUT_SLD_DM) ) && ( Xn <= 0 ) )
				{
					/* neutron SLD of dispersion medium [1e-6 A^(-2)] */
 					read_rho_and_sld( &sdummy[(int)(strptr1-strptr0)], &sld_multi.NEUT_SLD_DM_file, sld_multi.NEUT_SLD_DM_col, &sld_dm, 'n') ;
					sld_isdef[2] = true ;
				}
				else if ( ( !strcmp( item, par_keyword.NEUT_ICSD_INLAY) ) && ( Xn <= 0 ) )
				{
					/* neutron neutron ICSD of the inner stabilizer layer [cm^(-1)] */
 					read_rho_and_sld( &sdummy[(int)(strptr1-strptr0)], &sld_multi.NEUT_ICSD_INLAY_file, sld_multi.NEUT_ICSD_INLAY_col, &icsd_isl, 'n') ;
					icsd_isdef[0] = true ;
				}
				else if ( ( !strcmp( item, par_keyword.NEUT_ICSD_OUTLAY) ) && ( Xn <= 0 ) )
				{
					/* neutron neutron ICSD of the outer stabilizer layer [cm^(-1)] */
 					read_rho_and_sld( &sdummy[(int)(strptr1-strptr0)], &sld_multi.NEUT_ICSD_OUTLAY_file, sld_multi.NEUT_ICSD_OUTLAY_col, &icsd_osl, 'n') ;
					icsd_isdef[1] = true ;
				}
				else if ( ( !strcmp( item, par_keyword.NEUT_ICSD_DM) ) && ( Xn <= 0 ) )
				{
					/* neutron neutron ICSD of the dispersion medium [cm^(-1)] */
 					read_rho_and_sld( &sdummy[(int)(strptr1-strptr0)], &sld_multi.NEUT_ICSD_DM_file, sld_multi.NEUT_ICSD_DM_col, &icsd_dm, 'n') ;
					icsd_isdef[2] = true ;
				}
				else if ( conc_userdef && !strcmp( item, par_keyword.CONC_CRY) ) 
				{
					/* overwrite default value for conc_cry, if conc_userdef == true */

					/* vol% concentration of crystals in the nanosuspension */
					strptr2 = strchr (strptr1, '\n') ;
					i = (int)(strptr2-strptr1) ;
					memmove( data, strptr1, i) ;
					data[i] = 0 ;

					conc_cry = strtod (data,NULL) ;
					conc_isdef[0] = true ;
				}
				else if ( conc_userdef && !strcmp( item, par_keyword.CONC_DM) )
				{
					/* overwrite default value for conc_dm, if conc_userdef == true */

					/* vol% concentration of dispersion medium in the nanosuspension */
					strptr2 = strchr (strptr1, '\n') ;
					i = (int)(strptr2-strptr1) ;
					memmove( data, strptr1, i) ;
					data[i] = 0 ;

					conc_dm = strtod (data,NULL) ;
					conc_isdef[1] = true ;
				}
				/* else continue */
			}
		}

		/* close par-file */  
		fclose(inpf) ;

		/* finally apply checks if all necessary parameters were found in the par-file */
		if ( thickness_isdef[0] == false ) { XNDIFF_ERROR(49) ; }
		if ( thickness_isdef[1] == false ) { XNDIFF_ERROR(50) ; }

		if ( ( Xn >= 0 ) && ( rho_isdef[0] == false ) ) { XNDIFF_ERROR(51) ; }
		if ( ( Xn >= 0 ) && ( rho_isdef[1] == false ) ) { XNDIFF_ERROR(52) ; }
		if ( ( Xn >= 0 ) && ( rho_isdef[2] == false ) ) { XNDIFF_ERROR(53) ; }

		if ( ( Xn <= 0 ) && ( sld_isdef[0] == false ) ) { XNDIFF_ERROR(54) ; }
		if ( ( Xn <= 0 ) && ( sld_isdef[1] == false ) ) { XNDIFF_ERROR(55) ; }
		if ( ( Xn <= 0 ) && ( sld_isdef[2] == false ) ) { XNDIFF_ERROR(56) ; }

		if ( ( Xn <= 0 ) && ( icsd_isdef[0] == false ) ) { XNDIFF_ERROR(57) ; }
		if ( ( Xn <= 0 ) && ( icsd_isdef[1] == false ) ) { XNDIFF_ERROR(58) ; }
		if ( ( Xn <= 0 ) && ( icsd_isdef[2] == false ) ) { XNDIFF_ERROR(59) ; }

		if ( conc_userdef && ( conc_isdef[0] == false ) ) { XNDIFF_ERROR(60) ; }
		if ( conc_userdef && ( conc_isdef[1] == false ) ) { XNDIFF_ERROR(61) ; }
	}


	/* calculates the lattice vectors a1,a2,a3 and reciprocal lattice vectors b1,b2,b3 as well as 
	   the vectors G0,kg1,kg2 and Gs,k1,k2 in reciprocal space */
	void compute_lattice_vectors()
	{
		double ddummy1;

		/* Transform a1, a2, a3 to kartesian coordinates and calculate reciprocal unit cell vectors */
		/* x direction along vector a1 */
		a1[0] = cell_par[0] ;
		a1[1] = 0.0 ;
		a1[2] = 0.0 ;
  
		a2[0] = cell_par[1] * cos (DR*cell_par[5]) ;
		a2[1] = cell_par[1] * sin (DR*cell_par[5]) ;
		a2[2] = 0.0 ;
  
		a3[0] = cell_par[2] * cos (DR*cell_par[4]) ;
		a3[1] = cell_par[2] * (cos (DR*cell_par[3]) - cos (DR*cell_par[4]) * cos (DR*cell_par[5])) / sin (DR*cell_par[5]) ;
		a3[2] = sqrt (SQUARE(cell_par[2]) - SQUARE(a3[0]) - SQUARE(a3[1])) ;

		/* Volume of the unit cell by parallelepipedial product (determinant of |a1 a2 a3|) ) */
		V_uc = a1[0] * (a2[1]*a3[2] - a2[2]*a3[1]) + a1[1] * (a2[2]*a3[0] - a2[0]*a3[2]) + a1[2] * (a2[0]*a3[1] - a2[1]*a3[0]) ;

		/* reciprocal lattice vectors bi=1/V_uc*(bj x bk) */
		b1[0] = (a2[1]*a3[2] - a2[2]*a3[1]) / V_uc ;
		b1[1] = (a2[2]*a3[0] - a2[0]*a3[2]) / V_uc ;
		b1[2] = (a2[0]*a3[1] - a2[1]*a3[0]) / V_uc ;
  
		b2[0] = (a3[1]*a1[2] - a3[2]*a1[1]) / V_uc ;
		b2[1] = (a3[2]*a1[0] - a3[0]*a1[2]) / V_uc ;
		b2[2] = (a3[0]*a1[1] - a3[1]*a1[0]) / V_uc ;
  
		b3[0] = (a1[1]*a2[2] - a1[2]*a2[1]) / V_uc ;
		b3[1] = (a1[2]*a2[0] - a1[0]*a2[2]) / V_uc ;
		b3[2] = (a1[0]*a2[1] - a1[1]*a2[0]) / V_uc ;

		/* Calculate G0 and Gs */
		/* where G0 is the reciprocal lattice unit vector parallel to (s-s0) ??? and */
		/* Gs is the reciprocal lattice unit vector of stacking direction */
		/* both G0 and Gs should not be a zero vector */
		if ( par->h==0.0 && par->k==0.0 && par->l==0.0 )
		{
			XNDIFF_ERROR(10) ;
		}
		if ( par->hs==0.0 && par->ks==0.0 && par->ls==0.0 )
		{
			XNDIFF_ERROR(11) ;
		}
		/* G0=h*b1+k*b2+l*b3 */  
		G0[0] = par->h*b1[0] + par->k*b2[0] + par->l*b3[0] ;
		G0[1] = par->h*b1[1] + par->k*b2[1] + par->l*b3[1] ;
		G0[2] = par->h*b1[2] + par->k*b2[2] + par->l*b3[2] ;
		ddummy1 = betrag(G0) ;
		G0[0] /= ddummy1 ;
		G0[1] /= ddummy1 ;
		G0[2] /= ddummy1 ;
  
		Gs[0] = par->hs*b1[0] + par->ks*b2[0] + par->ls*b3[0] ;
		Gs[1] = par->hs*b1[1] + par->ks*b2[1] + par->ls*b3[1] ;
		Gs[2] = par->hs*b1[2] + par->ks*b2[2] + par->ls*b3[2] ;
		ddummy1 = betrag (Gs) ;
		Gs[0] /= ddummy1 ;
		Gs[1] /= ddummy1 ;
		Gs[2] /= ddummy1 ;

		/* write unit cell, reciprocal unit cell and G0, Gs into logfile */
		if ( log_flag )
		{
			fprintf( logfile, "\tUnit cell vectors in kartesian coordinates:\n") ;
			fprintf( logfile, "\ta1 = (%G %G %G)\n", a1[0], a1[1], a1[2]) ;
			fprintf( logfile, "\ta2 = (%G %G %G)\n", a2[0], a2[1], a2[2]) ;
			fprintf( logfile, "\ta3 = (%G %G %G)\n", a3[0], a3[1], a3[2]) ;
			fprintf( logfile, "\tReciprocal unit cell vectors in kartesian coordinates:\n") ;
			fprintf( logfile, "\tb1 = (%G %G %G)\n", b1[0], b1[1], b1[2]) ;
			fprintf( logfile, "\tb2 = (%G %G %G)\n", b2[0], b2[1], b2[2]) ;
			fprintf( logfile, "\tb3 = (%G %G %G)\n", b3[0], b3[1], b3[2]) ;
			fprintf( logfile, "\tUnit vector in direction of scattering vector:\n") ;
			fprintf( logfile, "\tG0 = (%G %G %G)\n", G0[0], G0[1], G0[2]) ;
			fprintf( logfile, "\tUnit vector in stacking direction:\n") ;
			fprintf( logfile, "\tGs = (%G %G %G)\n", Gs[0], Gs[1], Gs[2]) ;
			fflush (logfile) ;
		}

		/* Calculate Debye-Waller-factors */

		/* Calculate scattering function */

		/* Some definitions */

  

		/* kg1 and kg2: two vectors perpendicular to G */
  
		/* logic table, if h,k,l occur explicitely, they are non-zero */
		/* hkl		kg1	kg2	*/
		/* -----------------------------*/  
		/* 00l		done	done	*/
		/* 0kl,0k0	done	-	*/
		/* h0l,h00	done	-	*/
		/* hkl,hk0	done	-	*/ 
		if (par->h == 0 || par->k == 0)
		{
			/* 0kl,0k0 or h0l,h00 or 00l (000 is forbidden) */

			if (par->h == 0)
			{
				/* 0kl,0k0 or 00l */
				/* therefore kg1 along a1 */
				kg1[0] = a1[0] ;
				kg1[1] = a1[1] ;
				kg1[2] = a1[2] ;
				if (par->k == 0)
				{
					/* hkl=00l */ 
					/* therefore kg2 along a2 */
					/* but kg2=a2 is in general not perpendicular to a1 yet */
					kg2[0] = a2[0] ;
					kg2[1] = a2[1] ;
					kg2[2] = a2[2] ;
				}
			}
			else 
			{
				/* h0l,h00 */
				/* therefore kg1 along a1 ??? should it not be a2 ? -> yes, changed */
				/* if condition below seems to be redundant by the outer two if conditions */
// 				if (par->k == 0)
// 				{ 
// 					kg1[0] = a1[0] ;
// 					kg1[1] = a1[1] ;
// 					kg1[2] = a1[2] ;
// 				}
				if (par->k == 0)
				{ 
					kg1[0] = a2[0] ;
					kg1[1] = a2[1] ;
					kg1[2] = a2[2] ;
				}
			}
		}
		else
		{
			/* hkl,hk0 */
			kg1[0] = (a1[0]/par->h - a2[0]/par->k) ;
			kg1[1] = (a1[1]/par->h - a2[1]/par->k) ;
			kg1[2] = (a1[2]/par->h - a2[2]/par->k) ;
		}
		ddummy1 = betrag (kg1) ;
		kg1[0] /= ddummy1 ;
		kg1[1] /= ddummy1 ;
		kg1[2] /= ddummy1 ;
  

		/* logic table, if h,k,l occur explicitely, they are non-zero */
		/* hkl		kg1	kg2	*/
		/* -----------------------------*/  
		/* 00l		-	-	*/
		/* 0kl,0k0	-	done	*/
		/* h0l,h00	-	done	*/
		/* hkl,hk0	-	done	*/ 
		if (par->h != 0 || par->k != 0)
		{
			/* hkl,hk0 or 0kl,0k0 or h0l,h00 */

			if (par->l == 0)
			{
				/* hk0 or 0k0 or h00 */
				/* therefore kg2 along a3 */
				kg2[0] = a3[0] ;
				kg2[1] = a3[1] ;
				kg2[2] = a3[2] ;
			}
			else 
			{
				/* hkl or 0kl or h0l */
				if (par->k == 0)
				{
					/* h0l */
					kg2[0] = (a3[0]/par->l - a1[0]/par->h) ;
					kg2[1] = (a3[1]/par->l - a1[1]/par->h) ;
					kg2[2] = (a3[2]/par->l - a1[2]/par->h) ;
				}
				else 
				{
					/* hkl or 0kl */
					kg2[0] = (a2[0]/par->k - a3[0]/par->l) ;
					kg2[1] = (a2[1]/par->k - a3[1]/par->l) ;
					kg2[2] = (a2[2]/par->k - a3[2]/par->l) ;
				}
			}
		}
		ddummy1 = betrag (kg2) ;
		kg2[0] /= ddummy1 ;
		kg2[1] /= ddummy1 ;
		kg2[2] /= ddummy1 ;

		/* make kg2 perpendicular to kg1 */
		ddummy1 = sp (kg1, kg2) ;
		kg2[0] = (kg1[0] - 1/ddummy1 * kg2[0]) ;
		kg2[1] = (kg1[1] - 1/ddummy1 * kg2[1]) ;
		kg2[2] = (kg1[2] - 1/ddummy1 * kg2[2]) ;

		ddummy1 = betrag (kg2) ;
		kg2[0] /= ddummy1 ;
		kg2[1] /= ddummy1 ;
		kg2[2] /= ddummy1 ;
		/* from now on kg1,kg2,G0 perpendicular on each other and they are unit vectors */
  

		/* k1 and k2: two vectors perpendicular to Gs */
		/* exactly the same procedure for Gs as for G0 */
		/* use only par->hs,ks,ls instead of hkl for G0 */
		/* same error ??? as above */
		if (par->hs == 0 || par->ks == 0)
		{
			if (par->hs == 0)
			{
				k1[0] = a1[0] ;
				k1[1] = a1[1] ;
				k1[2] = a1[2] ;
				if (par->ks == 0)
				{
					k2[0] = a2[0] ;
					k2[1] = a2[1] ;
					k2[2] = a2[2] ;
				}
			}
			else
			{
				if (par->ks == 0)
				{
					k1[0] = a1[0] ;
					k1[1] = a1[1] ;
					k1[2] = a1[2] ;
				}
			}
		}
		else
		{
			k1[0] = (a1[0]/par->hs - a2[0]/par->ks) ;
			k1[1] = (a1[1]/par->hs - a2[1]/par->ks) ;
			k1[2] = (a1[2]/par->hs - a2[2]/par->ks) ;
		}
		ddummy1 = betrag (k1) ;
		k1[0] /= ddummy1 ;
		k1[1] /= ddummy1 ;
		k1[2] /= ddummy1 ;

		if (par->hs != 0 && par->ks != 0)
		{
			if (par->ls == 0)
			{
				k2[0] = a3[0] ;
				k2[1] = a3[1] ;
				k2[2] = a3[2] ;
			}
			else
			{
				if (par->ks == 0)
				{
					k2[0] = (a3[0]/par->ls - a1[0]/par->hs) ;
					k2[1] = (a3[1]/par->ls - a1[1]/par->hs) ;
					k2[2] = (a3[2]/par->ls - a1[2]/par->hs) ;
				}
				else
				{
					k2[0] = (a2[0]/par->ks - a3[0]/par->ls) ;
					k2[1] = (a2[1]/par->ks - a3[1]/par->ls) ;
					k2[2] = (a2[2]/par->ks - a3[2]/par->ls) ;
				}
			}
		}
		ddummy1 = betrag (k2) ;
		k2[0] /= ddummy1 ;
		k2[1] /= ddummy1 ;
		k2[2] /= ddummy1 ;

		/* k2 should be perpendicular to k1 */
		ddummy1 = sp (k1, k2) ;
		k2[0] = (k1[0] - 1/ddummy1 * k2[0]) ;
		k2[1] = (k1[1] - 1/ddummy1 * k2[1]) ;
		k2[2] = (k1[2] - 1/ddummy1 * k2[2]) ;
		ddummy1 = betrag (k2) ;
		k2[0] /= ddummy1 ;
		k2[1] /= ddummy1 ;
		k2[2] /= ddummy1 ;

		if ( log_flag )
		{
			fprintf( logfile, "\tVectors kg1 and kg2 perpendicular to G:\n") ;
			fprintf( logfile, "\tkg1 = (%G %G %G)\n", kg1[0], kg1[1], kg1[2]) ;
			fprintf( logfile, "\tkg2 = (%G %G %G)\n", kg2[0], kg2[1], kg2[2]) ;
			fprintf( logfile, "\tVectors k1 and k2 perpendicular to Gs:\n") ;
			fprintf( logfile, "\tk1  = (%G %G %G)\n", k1[0], k1[1], k1[2]) ;
			fprintf( logfile, "\tk2  = (%G %G %G)\n", k2[0], k2[1], k2[2]) ;
			fflush (logfile) ;
		}
		/* from now on k1,k2,Gs perpendicular on each other and they are unit vectors */
	}


	/* determine used memory during runtime */
	void get_current_memory()
	{
		FILE *file ;
		char buff[1024] ;
		char sdummy[1024] ;

		#ifdef __FreeBSD__
		sprintf( mem_filename, "%smem_%s_PID=%d.sh", ofo, "BSD", pid) ;
		#elif __linux__
		sprintf( mem_filename, "%smem_%s_PID=%d.sh", ofo, "GLX", pid) ;
		#endif

		double mem_size = -1.0 ;

		/* write bash script if it does not exist yet in order to get the current memory */

		/* (over)write file */
		if ( log_flag ) { fprintf( logfile, "%sWriting bash script %s\n", "\t", mem_filename) ; fflush (logfile) ; }
		if ( ( file = fopen( mem_filename, "w")) == NULL ) { XNDIFF_ERROR(72) ; }

		fprintf( file, "#!/bin/bash\n#\n") ;

// 		#ifdef __FreeBSD__
// 		fprintf( file, "string=`ps -o vsz -p %d`\n", pid) ;
// 		fprintf( file, "echo ${string#*VSZ*}\n") ;

// 		#elif __linux__
		fprintf( file, "string=`cat /proc/%d/status | grep VmSize`\n", pid) ;
		fprintf( file, "echo $string | awk '{print $2}'\n") ;
// 		fprintf( file, "echo ${string/VmSize:/}\n") ;

// 		#endif

		fclose(file) ;
		if ( log_flag ) { fprintf( logfile, "%sdone\n\n", "\t") ; fflush (logfile) ; }

		/* change executable rights */
		if ( log_flag ) { fprintf( logfile, "%sChange executable rights for %s\n", "\t", mem_filename) ; fflush (logfile) ; }
		sprintf( sdummy,"chmod +x %s", mem_filename) ;
		if ( ( file = popen( sdummy, "r") ) == NULL ) { XNDIFF_ERROR(73) ; }
		pclose(file) ;
		if ( log_flag ) { fprintf( logfile, "%sdone\n\n", "\t") ; fflush (logfile) ; }

		/* run the bash script and read out the current memory, the
		   number comes in KByte for FreeBSD and Linux, convert to MByte */
		if ( log_flag ) { fprintf( logfile, "%sRunning bash script %s\n", "\t", mem_filename) ; fflush (logfile) ; }
		sprintf( sdummy,"./%s", mem_filename) ;
		if ( ( file = popen( sdummy, "r")) == NULL ) { XNDIFF_ERROR(74) ; }

		while ( fgets( buff, sizeof(buff), file) != NULL )
		{
			if ( log_flag ) { fprintf( logfile, "%sbash output = \"%s\"\n", "\t", buff) ; fflush (logfile) ; }
			delallspc(buff) ;
			/* remove \n at the end of buff what may cause errors in strtod */
			for ( unsigned int j=0; j<sizeof(buff); ++j)
			{
				if ( buff[j] == '\n' ) { buff[j] = 0 ; break ; }
			}
			if ( is_numeric(buff) )
			{
				mem_size = strtod( buff, NULL) ;
				mem_size /= 1024.0 ;
			}
			else
			{
				if ( log_flag ) { fprintf( logfile, "%sCannot convert output to number.\n", "\t") ; fflush (logfile) ; }
			}
		}

		/* close the pipe and output result */
		pclose(file) ;
		if ( log_flag ) { fprintf( logfile, "%sdone\n\n", "\t") ; fflush (logfile) ; }

		if ( log_flag ) { fprintf( logfile, "%sMemory Usage: %.3f [MByte]\n", "\t", mem_size) ; }
	}


	/* tool, that allows to compute and export the orientational averaged structure factor |F(s)|^2 for neutron and X-ray scattering.
	   mode is an optional parameter and is set by default to par->av_mode. For X-ray |F(s)|^2 [nm^2] and same for 
	   neutrons |F(s)|^2 [nm^2] */
	/* DO NOT call this function before FFF_X/n are computed, this is fulfilled ...
	   -if FNoMem_flag is not set, after calling compute_structure_amplitude()
	   -if FNoMem_flag is set, after c_rep=0 has been passed
	*/
	void compute_orientational_averaged_structure_factor(dcmplx ***FFF_Xn, char X_or_n, int mode = -1 )
	{
 		FILE* outf ; /* FILE pointer for output file */
		char sdummy[1024] ;
		double ABS_F, cs ;
		unsigned int MULTIPLICITY ;
		int idummy ;

		double weight_factor, sum_weight_factor, SUM_WEIGHT_CORR ;

		/* if mode = -1 isn't overwritten by a user-supplied mode parameter, use by default par->av_mode */
		if ( mode == -1 ) { mode = par->av_mode ; }

		if ( X_or_n == 'n') { sprintf( sdummy, "%s%s%s.dat", ofo, mfname, "_F2_n") ; }
		else if ( X_or_n == 'X') { sprintf( sdummy, "%s%s%s.dat", ofo, mfname, "_F2_X") ; }
		if ( ( outf = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(19) ; }

		fprintf( outf, "# Orientationally averaged structure amplitude F(s)\n") ;
		if ( X_or_n == 'n') { fprintf( outf, "# for neutron scattering\n") ; }
		else if ( X_or_n == 'X') { fprintf( outf, "# for X-ray scattering\n") ; }

		if ( ( mode == 0 ) || ( mode == 2 ) || ( mode == 3 ) )
		{
			MULTIPLICITY = nav_pc * nav_ac ;
			sum_weight_factor = 0.0 ;
			for ( unsigned int k=0; k<nav_ac; ++k)
			{
				for ( unsigned int j=0; j<nav_pc; ++j)
				{
					if ( weight_flag ) { weight_factor = sin( DR * av_pol_ang.at(j) ) ; }
					else { weight_factor = 1.0 ; }
					sum_weight_factor += weight_factor ;
				}
			}
			SUM_WEIGHT_CORR = sum_weight_factor / MULTIPLICITY ;

			fprintf( outf, "# MULTIPLICITY = %d ( nav_ac = %d, nav_pc = %d)\n", MULTIPLICITY, nav_ac, nav_pc) ;
			fprintf( outf, "# SUM_WEIGHT_FACTOR = %15.10G\n", sum_weight_factor) ;
			fprintf( outf, "# SUM_WEIGHT_CORR = %15.10G\n", SUM_WEIGHT_CORR ) ;
			fprintf( outf, "# %15s", "s [1/nm]" ) ;
			if ( X_or_n == 'n' )
			{
				sprintf( sdummy, "|F[s]|^2 [%s^2]", unit_F_n) ;
				fprintf( outf, " %15s\n", sdummy) ;
			}
			else if ( X_or_n == 'X' )
			{
				sprintf( sdummy, "|F[s]|^2 [%s^2]", unit_F_X) ;
				fprintf( outf, " %15s\n", sdummy) ;
			}

			cs = s1 ;
			for ( unsigned int i=0; i<par->np; ++i)
			{
				ABS_F = 0.0;
				cs += ds;

				/* orientational average over azimuthal and polar angles */
				for ( unsigned int j=0; j<nav_pc; ++j)
				{
					/* optionally apply a weighting factor ~ sin / (\sum sin ) */
					if ( weight_flag ) { weight_factor = sin( DR * av_pol_ang.at(j) ) ; }
					else { weight_factor = 1.0 ; }

					for ( unsigned int k=0; k<nav_ac; ++k)
					{
						ABS_F += fabs( real( FFF_Xn[i][j][k] * conj(FFF_Xn[i][j][k]) ) ) * weight_factor ;
					}
				}

				fprintf( outf, "  %15.10G %15.10G\n", cs, ABS_F / ( (double)MULTIPLICITY * SUM_WEIGHT_CORR )  ) ;
			}
		}
		else if ( mode == 1 )
		{
			/* better avoid this mode when using FNoMem_flag, errors will be likely */

			/* When FNoMem_flag is set only nav_pc*par->nms out of the total nav_pc*nav_r3 entries of FFF are non-zero if c_rep=0 was limited !!
			   This due to the fact that for each particle in the stack (c_rep=0) one random orientation is chosen, which are hopefully all distinct.
			   Hence use par->nms as normalization instead of nav_r3.
			   However when randomly multiple occurences of a r3-rotation occured in c_rep=0 the result will not be correct in this mode !!
			*/

			/* If FNoMem_flag is not set, the computation should be fine, cause typically FFF_X/n is finely sampled over the relevant orientations in this mode, 
			   allowing later "all" possible random orientations 
			*/
			if ( !FNoMem_flag ) { idummy = nav_r3 ; }
			else { idummy = par->nms ; }

			MULTIPLICITY = nav_pc * idummy ;
			sum_weight_factor = 0.0 ;
			for ( unsigned int j=0; j<nav_pc; ++j)
			{
				if ( weight_flag ) { weight_factor = sin( DR * av_pol_ang.at(j) ) ; }
				else { weight_factor = 1.0 ; }

				for ( int k=0; k<idummy; ++k)
				{
					sum_weight_factor += weight_factor ;
				}
			}
			SUM_WEIGHT_CORR = sum_weight_factor / MULTIPLICITY ;

			fprintf( outf, "# MULTIPLICITY = %d ( nav_r3 = %d, nav_pc = %d)\n", MULTIPLICITY, idummy, nav_pc) ;
			fprintf( outf, "# SUM_WEIGHT_FACTOR = %15.10G\n", sum_weight_factor) ;
			fprintf( outf, "# SUM_WEIGHT_CORR = %15.10G\n", SUM_WEIGHT_CORR ) ;
			fprintf( outf, "# %15s", "s [1/nm]" ) ;
			if ( X_or_n == 'n' )
			{
				sprintf( sdummy, "|F[s]|^2 [%s^2]", unit_F_n) ;
				fprintf( outf, " %15s\n", sdummy) ;
			}
			else if ( X_or_n == 'X' )
			{
				sprintf( sdummy, "|F[s]|^2 [%s^2]", unit_F_X) ;
				fprintf( outf, " %15s\n", sdummy) ;
			}

			cs = s1 ;
			for ( unsigned int i=0; i<par->np; ++i)
			{
				ABS_F = 0.0;
				cs += ds;

				/* orientational average over r3-rotation and polar angles */
				for ( unsigned int j=0; j<nav_pc; ++j)
				{
					/* optionally apply a weighting factor ~ sin / (\sum sin ) */
					if ( weight_flag ) { weight_factor = sin( DR * av_pol_ang.at(j) ) ; }
					else { weight_factor = 1.0 ; }

					for ( unsigned int k=0; k<nav_r3; ++k)
					{
						ABS_F += fabs( real( FFF_Xn[i][j][k] * conj(FFF_Xn[i][j][k]) ) ) * weight_factor ;
					}
				}
				fprintf( outf, "%15.10G %15.10G\n", cs, ABS_F / ( (double) MULTIPLICITY * SUM_WEIGHT_CORR ) ) ;
			}
		}
		fclose(outf) ;
	}


	/* tool, that allows to compute and export the orientational averaged complex structure amplitude F(s) for neutron and X-ray scattering.
	   mode is an optional parameter and is set by default to par->av_mode. For X-ray and neutrons F(s) has units [nm]
	*/
	/* DO NOT call this function before FFF_X/n are computed, this is fulfilled ...
	   -if FNoMem_flag is not set, after calling compute_structure_amplitude()
	   -if FNoMem_flag is set, after c_rep=0 has been passed
	*/
	void compute_orientational_averaged_structure_amplitude(dcmplx ***FFF_Xn, char X_or_n, int mode = -1 )
	{
 		FILE* outf ; /* FILE pointer for output file */
		char sdummy[1024] ;
		double cs ;
		dcmplx F ;
		unsigned int MULTIPLICITY ;
		int idummy ;

		double weight_factor, sum_weight_factor, SUM_WEIGHT_CORR ;

		/* if mode = -1 isn't overwritten by a user-supplied mode parameter, use by default par->av_mode */
		if ( mode == -1 ) { mode = par->av_mode ; }

		if ( X_or_n == 'n') { sprintf( sdummy, "%s%s%s.dat", ofo, mfname, "_F_n") ; }
		else if ( X_or_n == 'X') { sprintf( sdummy, "%s%s%s.dat", ofo, mfname, "_F_X") ; }
		if ( ( outf = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(19) ; }

		fprintf( outf, "# Orientationally averaged structure amplitude F(s)\n") ;
		if ( X_or_n == 'n') { fprintf( outf, "# for neutron scattering\n") ; }
		else if ( X_or_n == 'X') { fprintf( outf, "# for X-ray scattering\n") ; }

		if ( ( mode == 0 ) || ( mode == 2 ) || ( mode == 3 ) )
		{
			MULTIPLICITY = nav_pc * nav_ac ;
			sum_weight_factor = 0.0 ;
			for ( unsigned int k=0; k<nav_ac; ++k)
			{
				for ( unsigned int j=0; j<nav_pc; ++j)
				{
					if ( weight_flag ) { weight_factor = sin( DR * av_pol_ang.at(j) ) ; }
					else { weight_factor = 1.0 ; }
					sum_weight_factor += weight_factor ;
				}
			}
			SUM_WEIGHT_CORR = sum_weight_factor / MULTIPLICITY ;

			fprintf( outf, "# MULTIPLICITY = %d ( nav_ac = %d, nav_pc = %d)\n", MULTIPLICITY, nav_ac, nav_pc) ;
			fprintf( outf, "# SUM_WEIGHT_FACTOR = %15.10G\n", sum_weight_factor) ;
			fprintf( outf, "# SUM_WEIGHT_CORR = %15.10G\n", SUM_WEIGHT_CORR ) ;
			fprintf( outf, "# %15s", "s [1/nm]" ) ;
			if ( X_or_n == 'n' )
			{
				sprintf( sdummy, "Re(F[s]) [%s]", unit_F_n) ;
				fprintf( outf, " %15s", sdummy) ;
				sprintf( sdummy, "Im(F[s]) [%s]", unit_F_n) ;
				fprintf( outf, " %15s\n", sdummy) ;
			}
			else if ( X_or_n == 'X' )
			{
				sprintf( sdummy, "Re(F[s]) [%s]", unit_F_X) ;
				fprintf( outf, " %15s", sdummy) ;
				sprintf( sdummy, "Im(F[s]) [%s]", unit_F_X) ;
				fprintf( outf, " %15s\n", sdummy) ;
			}

			cs = s1 ;
			for ( unsigned int i=0; i<par->np; ++i)
			{
				F = dcmplx(0.0,0.0) ;
				cs += ds;

				/* orientational average over azimuthal and polar angles */
				for ( unsigned int j=0; j<nav_pc; ++j)
				{
					/* optionally apply a weighting factor ~ sin / (\sum sin ) */
					if ( weight_flag ) { weight_factor = sin( DR * av_pol_ang.at(j) ) ; }
					else { weight_factor = 1.0 ; }

					for ( unsigned int k=0; k<nav_ac; ++k)
					{
						F += FFF_Xn[i][j][k] * (dcmplx)weight_factor ;
					}
				}
				fprintf( outf, "  %15.10G %15.10G %15.10G\n", cs, real(F) / ( (double)MULTIPLICITY * SUM_WEIGHT_CORR ), imag(F) / ( (double)MULTIPLICITY * SUM_WEIGHT_CORR ) ) ;
			}
		}
		else if ( mode == 1 )
		{
			/* better avoid this mode when using FNoMem_flag, errors will be likely */

			/* When FNoMem_flag is set only nav_pc*par->nms out of the total nav_pc*nav_r3 entries of FFF are non-zero if c_rep=0 was limited !!
			   This due to the fact that for each particle in the stack (c_rep=0) one random orientation is chosen, which are hopefully all distinct.
			   Hence use par->nms as normalization instead of nav_r3.
			   However when randomly multiple occurences of a r3-rotation occured in c_rep=0 the result will not be correct in this mode !!
			*/

			/* If FNoMem_flag is not set, the computation should be fine, cause typically FFF_X/n is finely sampled over the relevant orientations in this mode, 
			   allowing later "all" possible random orientations 
			*/
			if ( !FNoMem_flag ) { idummy = nav_r3 ; }
			else { idummy = par->nms ; }

			MULTIPLICITY = nav_pc * idummy ;
			sum_weight_factor = 0.0 ;
			for ( int k=0; k<idummy; ++k)
			{
				for ( unsigned int j=0; j<nav_pc; ++j)
				{
					if ( weight_flag ) { weight_factor = sin( DR * av_pol_ang.at(j) ) ; }
					else { weight_factor = 1.0 ; }
					sum_weight_factor += weight_factor ;
				}
			}
			SUM_WEIGHT_CORR = sum_weight_factor / MULTIPLICITY ;

			fprintf( outf, "# MULTIPLICITY = %d ( nav_r3 = %d, nav_pc = %d)\n", MULTIPLICITY, idummy, nav_pc) ;
			fprintf( outf, "# SUM_WEIGHT_FACTOR = %15.10G\n", sum_weight_factor) ;
			fprintf( outf, "# SUM_WEIGHT_CORR = %15.10G\n", SUM_WEIGHT_CORR ) ;
			fprintf( outf, "# %15s", "s [1/nm]" ) ;
			if ( X_or_n == 'n' )
			{
				sprintf( sdummy, "Re(F[s]) [%s]", unit_F_n) ;
				fprintf( outf, " %15s", sdummy) ;
				sprintf( sdummy, "Im(F[s]) [%s]", unit_F_n) ;
				fprintf( outf, " %15s\n", sdummy) ;
			}
			else if ( X_or_n == 'X' )
			{
				sprintf( sdummy, "Re(F[s]) [%s]", unit_F_X) ;
				fprintf( outf, " %15s", sdummy) ;
				sprintf( sdummy, "Im(F[s]) [%s]", unit_F_X) ;
				fprintf( outf, " %15s\n", sdummy) ;
			}

			cs = s1 ;
			for ( unsigned int i=0; i<par->np; ++i)
			{
				F = dcmplx(0.0,0.0) ;
				cs += ds;
	
				/* orientational average over r3-rotation and polar angles */
				for ( unsigned int j=0; j<nav_pc; ++j)
				{
					/* optionally apply a weighting factor ~ sin / (\sum sin ) */
					if ( weight_flag ) { weight_factor = sin( DR * av_pol_ang.at(j) ) ; }
					else { weight_factor = 1.0 ; }

					for ( unsigned int k=0; k<nav_r3; ++k)
					{
						F += FFF_Xn[i][j][k] * (dcmplx)weight_factor ;
					}
				}
				fprintf( outf, "  %15.10G %15.10G %15.10G\n", cs, real(F) / ( (double)MULTIPLICITY * SUM_WEIGHT_CORR ), imag(F) / ( (double)MULTIPLICITY * SUM_WEIGHT_CORR ) ) ;
			}
		}
		fclose(outf) ;
	}


	/* Tool, that allows to calculate and export the orientational averaged structure factor |F(s)|^2
	   for neutron and X-ray scattering using the analytical solution (for b_coh with zero imaginary part).
	   For X-ray and for neutrons |F(s)|^2 has units [nm^2].

	   Note that before calling this procedure,
	   add_neutron_coh_sc_length_and_inc_sc_crossec_to_atom()
	   compute_f0table_and_add_f0ind_to_atom()
	   should have been executed.
	*/
	void calculate_orientational_averaged_structure_factor()
	{
 		FILE *outf_X, *outf_n ; /* FILE pointers for output files */
		char sdummy[1024] ;
		double F2_X, F2_n, cs ;
		double dummy1[3], dummy2[3], dummy3[3] ;
		double sinc ;

		if ( Xn >= 0 )
		{
			sprintf( sdummy, "%s%s%s.dat", ofo, mfname, "_F2_X_exact") ;
			if ( ( outf_X = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(19) ; }

			fprintf( outf_X, "# Averaged structure factor |F(s)|^2 versus s \n") ;
			fprintf( outf_X, "# for X-ray scattering \n") ;
			sprintf( sdummy, "|F[s]|^2 [%s^2]", unit_F_X) ; fprintf( outf_X, "# %15s %15s\n", "s [1/nm]", sdummy) ;
		}
		if ( Xn <= 0 )
		{
			sprintf( sdummy, "%s%s%s.dat", ofo, mfname, "_F2_n_exact") ;
			if ( ( outf_n = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(19) ; }

			fprintf( outf_n, "# Averaged structure factor |F(s)|^2 versus s \n") ;
			fprintf( outf_n, "# for neutron scattering \n") ;
			sprintf( sdummy, "|F[s]|^2 [%s^2]", unit_F_n) ; fprintf( outf_n, "# %15s %15s\n", "s [1/nm]", sdummy) ;
		}

		cs = s1 ;
		for ( unsigned int i=0; i<par->np; ++i)
		{
			F2_X = 0.0 ;
			F2_n = 0.0 ;

			/* l == j == 0 */
			if ( Xn >= 0 ) { F2_X += f0table[atom[0].f0ind][i] * f0table[atom[0].f0ind][i] ; }
			if ( Xn <= 0 ) { F2_n += real( atom[0].coh_sc_length ) * real( atom[0].coh_sc_length ) ; }

			cs += ds;
			for ( unsigned int l=1; l<noa; ++l)
			{
				/* terms with l == j and l>0 */
				if ( Xn >= 0 ) { F2_X += f0table[atom[l].f0ind][i] * f0table[atom[l].f0ind][i] ; }
				if ( Xn <= 0 ) { F2_n += real( atom[l].coh_sc_length ) * real( atom[l].coh_sc_length ) ; }

				dummy2[0] = dummy2[1] = dummy2[2] = 0.0 ;
				csp( dummy1, atom[l].coord[0], a1) ;
				vadd( dummy2, dummy1) ;
				csp( dummy1, atom[l].coord[1], a2) ;
				vadd( dummy2, dummy1) ;
				csp( dummy1, atom[l].coord[2], a3) ;
				vadd( dummy2, dummy1) ;

				/* cross terms j<l */
				for ( unsigned int j=0; j<l; ++j)
				{
					dummy3[0] = dummy3[1] = dummy3[2] = 0.0 ;
					csp( dummy1, atom[j].coord[0], a1) ;
					vadd( dummy3, dummy1) ;
					csp( dummy1, atom[j].coord[1], a2) ;
					vadd( dummy3, dummy1) ;
					csp( dummy1, atom[j].coord[2], a3) ;
					vadd( dummy3, dummy1) ;

					vsub ( dummy3, dummy2) ;
		
					sinc = sin (2.0 * M_PI * cs * betrag (dummy3) ) / (2.0 * M_PI * cs * betrag (dummy3) ) ;
					sinc *= atom[l].mult ;
					if ( Xn >= 0 ) { F2_X += 2.0 * sinc * f0table[atom[j].f0ind][i] * f0table[atom[l].f0ind][i] ; /* [electrons] */ }

					if ( Xn <= 0 ) { F2_n += 2.0 * sinc * ( real( atom[j].coh_sc_length ) * real( atom[l].coh_sc_length ) ) ; /* [fm] */ }
				}
			}
			if ( Xn >= 0 ) { fprintf( outf_X, "  %15.10G %15.10G\n", cs, F2_X * pow( sc_rho, 2) ) ; } /* [nm^2] */
			if ( Xn <= 0 ) { fprintf( outf_n, "  %15.10G %15.10G\n", cs, F2_n * pow( sc_scl, 2) ) ; } /* [nm^2] */

		}
		if ( Xn >= 0 ) { fclose(outf_X) ; }
		if ( Xn <= 0 ) { fclose(outf_n) ; }
	}



	/* Tool, that allows to calculate and export the orientational averaged structure factor |F(s)|^2
	   for X-ray scattering using the analytical solution. |F(s)|^2 has units [nm^2].

	   Note that before calling this procedure,
	   add_neutron_coh_sc_length_and_inc_sc_crossec_to_atom()
	   compute_f0table_and_add_f0ind_to_atom()
	   should have been executed.
	*/
	double ff_sphere( double cs, double r_sph)
	{
		return 3.0 * ( sin( 2.0 * M_PI * cs * r_sph) - (2.0 * M_PI * cs * r_sph) * cos( 2.0 * M_PI * cs * r_sph) ) / pow( 2.0 * M_PI * cs * r_sph, 3.0 ) ; 
	}

	void calculate_orientational_averaged_structure_factor_in_solvent()
	{
 		FILE *outf_X ; /* FILE pointers for output files */
		char sdummy[1024] ;
		double F2_X, cs ;
		double dummy1[3], dummy2[3], dummy3[3] ;
		double sinc ;
		double rho_0, rho_l, rho_j, Z_solv, r_sph ; 

		/* v3 
		Z_solv = ( 2.0 * get_atomic_scattering_factor( "H", 0.0) + get_atomic_scattering_factor( "O", 0.0) ) ;
		r_sph = 0.193 ;
		*/

		// v4
		Z_solv = ( 2.0 * get_atomic_scattering_factor( "H", 0.0) + get_atomic_scattering_factor( "O", 0.0) ) ;

		if ( Xn >= 0 )
		{
			sprintf( sdummy, "%s%s%s.dat", ofo, mfname, "_F2_X_H2O_exact") ;
			if ( ( outf_X = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(19) ; }

			fprintf( outf_X, "# Averaged structure factor |F(s)|^2 versus s \n") ;
			fprintf( outf_X, "# for X-ray scattering \n") ;
			sprintf( sdummy, "|F[s]|^2 [%s^2]", unit_F_X) ; fprintf( outf_X, "# %15s %15s\n", "s [1/nm]", sdummy) ;

			cs = s1 ;
			for ( unsigned int i=0; i<par->np; ++i)
			{
				cs += ds;
				F2_X = 0.0 ;

				/* l == j == 0 */
				// F2_X += f0table[atom[0].f0ind][i] * f0table[atom[0].f0ind][i] ;
				// v1 rho_0 = f0table[atom[0].f0ind][i] - 2.0 * get_atomic_scattering_factor( "H", cs) - get_atomic_scattering_factor( "O", cs) ;
				
				/* v2 
				if ( !strcmp( atom[0].type, "H" ) ) { rho_0 = 0.0 ; } 
				else { rho_0 = f0table[atom[0].f0ind][i] - get_atomic_scattering_factor( "O", cs) ; }
				*/ 
				// v3 rho_0 = f0table[atom[0].f0ind][i] - Z_solv * ff_sphere( cs, r_sph) ;
				
				/* v4 */
				if ( !strcmp( atom[0].type, "H" ) ) { r_sph = 0.032 ; }
 				else if ( !strcmp( atom[0].type, "C" ) ) { r_sph = 0.077 ; }
 				else { r_sph = 0.066 ; } /* "O" */
				rho_0 = f0table[atom[0].f0ind][i] - Z_solv * ff_sphere( cs, r_sph) ;

				F2_X += rho_0 * rho_0 ;

				for ( unsigned int l=1; l<noa; ++l)
				{
					/* terms with l == j and l>0 */
					// F2_X += f0table[atom[l].f0ind][i] * f0table[atom[l].f0ind][i] ;
					// v1 rho_l = f0table[atom[l].f0ind][i] - 2.0 * get_atomic_scattering_factor( "H", cs) - get_atomic_scattering_factor( "O", cs) ;
					/* v2
					if ( !strcmp( atom[l].type, "H" ) ) { rho_l = 0.0 ; } 
					else { rho_l = f0table[atom[l].f0ind][i] - get_atomic_scattering_factor( "O", cs) ; }
					*/
					// v3 rho_l = f0table[atom[l].f0ind][i] - Z_solv * ff_sphere( cs, r_sph) ;

					/* v4 */
					if ( !strcmp( atom[l].type, "H" ) ) { r_sph = 0.032 ; }
					else if ( !strcmp( atom[l].type, "C" ) ) { r_sph = 0.077 ; }
					else { r_sph = 0.066 ; } /* "O" */
					rho_l = f0table[atom[l].f0ind][i] - Z_solv * ff_sphere( cs, r_sph) ;

					F2_X += rho_l * rho_l ;

					dummy2[0] = dummy2[1] = dummy2[2] = 0.0 ;
					csp( dummy1, atom[l].coord[0], a1) ;
					vadd( dummy2, dummy1) ;
					csp( dummy1, atom[l].coord[1], a2) ;
					vadd( dummy2, dummy1) ;
					csp( dummy1, atom[l].coord[2], a3) ;
					vadd( dummy2, dummy1) ;

					/* cross terms j<l */
					for ( unsigned int j=0; j<l; ++j)
					{
						dummy3[0] = dummy3[1] = dummy3[2] = 0.0 ;
						csp( dummy1, atom[j].coord[0], a1) ;
						vadd( dummy3, dummy1) ;
						csp( dummy1, atom[j].coord[1], a2) ;
						vadd( dummy3, dummy1) ;
						csp( dummy1, atom[j].coord[2], a3) ;
						vadd( dummy3, dummy1) ;

						vsub ( dummy3, dummy2) ;
		
						sinc = sin (2.0 * M_PI * cs * betrag (dummy3) ) / (2.0 * M_PI * cs * betrag (dummy3) ) ;
						// F2_X += 2.0 * sinc * f0table[atom[j].f0ind][i] * f0table[atom[l].f0ind][i] ; /* [electrons] */
						// v1 rho_j = f0table[atom[j].f0ind][i] - 2.0 * get_atomic_scattering_factor( "H", cs) - get_atomic_scattering_factor( "O", cs) ;
						/* v2
						if ( !strcmp( atom[j].type, "H" ) ) { rho_j = 0.0 ; } 
						else { rho_j = f0table[atom[j].f0ind][i] - get_atomic_scattering_factor( "O", cs) ; }
						*/
						
						// v3 rho_j = f0table[atom[j].f0ind][i] - Z_solv * ff_sphere( cs, r_sph) ;

						/* v4 */
						if ( !strcmp( atom[j].type, "H" ) ) { r_sph = 0.032 ; }
		 				else if ( !strcmp( atom[j].type, "C" ) ) { r_sph = 0.077 ; }
		 				else { r_sph = 0.066 ; } /* "O" */
						rho_j = f0table[atom[j].f0ind][i] - Z_solv * ff_sphere( cs, r_sph) ;

						F2_X += 2.0 * sinc * rho_j * rho_l ; /* [electrons] */
					}
				}
				fprintf( outf_X, "  %15.10G %15.10G\n", cs, F2_X * pow( sc_rho, 2) ) ; /* [nm^2] */
			}
			fclose(outf_X) ;
		}
	}




	/* Tool, that allows to calculate and export the exact orientational averaged structure amplitude F(s),
	   for neutron (using only real part of coherent scattering length) and X-ray scattering, using the analytical solution (for b_coh with zero imaginary part).
	   For X-ray and for neutrons F(s) has units [nm].

	   Note that before calling this procedure,
	   add_neutron_coh_sc_length_and_inc_sc_crossec_to_atom()
	   compute_f0table_and_add_f0ind_to_atom()
	   should have been executed.

	   Note that -b_coh in Mathematica -> should not have any inluence
	*/
	void calculate_orientational_averaged_structure_amplitude()
	{
 		FILE *outf_X, *outf_n ; /* FILE pointers for output files */
		double F_X, F_n, cs ;
		char sdummy[1024] ;
		double dummy1[3], dummy2[3] ;
		double sinc ;

		if ( Xn >= 0 )
		{
			sprintf( sdummy, "%s%s%s.dat", ofo, mfname, "_F_X_exact") ;
			if ( ( outf_X = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(19) ; }

			fprintf( outf_X, "# Averaged structure amplitude F(s) versus s \n") ;
			fprintf( outf_X, "# for X-ray scattering \n") ;
			sprintf( sdummy, "F[s] [%s]", unit_F_X) ; fprintf( outf_X, "# %15s %15s\n", "s [1/nm]", sdummy) ;
		}
		if ( Xn <= 0 )
		{
			sprintf( sdummy, "%s%s%s.dat", ofo, mfname, "_F_n_exact") ;
			if ( ( outf_n = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(19) ; }

			fprintf( outf_n, "# Averaged structure amplitude F(s) versus s \n") ;
			fprintf( outf_n, "# for neutron scattering \n") ;
			sprintf( sdummy, "F[s] [%s]", unit_F_n) ; fprintf( outf_n, "# %15s %15s\n", "s [1/nm]", sdummy) ;
		}

		cs = s1 ;
		for ( unsigned int i=0; i<par->np; ++i)
		{
			F_X = 0.0;
			F_n = 0.0;

			cs += ds;

			for ( unsigned int l=0; l<noa; ++l)
			{
				dummy2[0] = dummy2[1] = dummy2[2] = 0.0 ;
				csp( dummy1, atom[l].coord[0], a1) ;
				vadd( dummy2, dummy1) ;
				csp( dummy1, atom[l].coord[1], a2) ;
				vadd( dummy2, dummy1) ;
				csp( dummy1, atom[l].coord[2], a3) ;
				vadd( dummy2, dummy1) ;
	
				sinc = sin (2.0 * M_PI * cs * betrag (dummy2) ) / (2.0 * M_PI * cs * betrag (dummy2) ) ;
				sinc *= atom[l].mult ;
				/* compute orientational averaged structure amplitude F for X-ray with the atomic scattering factors */
				if ( Xn >= 0 ) { F_X += sinc * f0table[atom[l].f0ind][i] ; } /* [electrons] */ 

				/* the same for neutrons with the real part of the coherent neutron scattering lengths
				   don't use minus sign convention for b_coh
				*/
				if ( Xn <= 0 ) { F_n += sinc * real( atom[l].coh_sc_length ) ; } /* [fm] */
			}
			if ( Xn >= 0 ) { fprintf( outf_X, "  %15.10G %15.10G\n", cs, F_X * sc_rho ) ; } /* [nm] */
			if ( Xn <= 0 ) { fprintf( outf_n, "  %15.10G %15.10G\n", cs, F_n * sc_scl) ; } /* [nm] */

		}
		if ( Xn >= 0 ) { fclose(outf_X) ; }
		if ( Xn <= 0 ) { fclose(outf_n) ; }
	}


	/* read the structure amplitudes FFF_X and FFF_n. First check if the setups were the 
           same as it is specified in the "Identifiers" section. Automatically check which file is for 
	   neutrons and which for X-ray */
	bool read_structure_amplitude()
	{
		int mode = par->av_mode ;

		/* The implementation of reading / writing of FFF_X/n in the case of random orientations makes little sense */
		if ( ( mode == 2 ) || ( mode == 3 ) )
		{
			if ( log_flag ) { fprintf( logfile, "\tWarning: Reading FFF_X/n for mode=%d is not supported !\n", mode) ; fflush (logfile) ; }
			return false ;
		}

		/* The implementation of reading of FFF_X/n in the case when FNoMem_flag is set makes little sense */
		if ( FNoMem_flag )
		{
			if ( log_flag ) { fprintf( logfile, "\tWarning: Reading FFF_X/n when the FNoMem_flag is set, is not supported !\n") ; fflush (logfile) ; }
			return false ;
		}

		if ( Xn != 0 && num_Fr_flag != 1 ) { if (log_flag ) { fprintf( logfile, "\tExpected 1 structure amplitude file\n") ; fflush (logfile) ; } ; return false ; }
		if ( Xn == 0  && num_Fr_flag != 2 ) { if (log_flag ) { fprintf( logfile, "\tExpected 2 structure amplitude files\n") ; fflush (logfile) ; } ; return false ; }


		unsigned int nav_rot = 0 ;

		if ( mode == 0 ) { nav_rot = nav_ac ; }
		else if ( mode == 1 ) { nav_rot = nav_r3 ; }

		FILE *inpf ;
		const int numsigns = 1024 ;
		char sdummy[numsigns] ;
		char *sdummy2 ; 
		char sdummy3[100] ; 
		double ddummy = 0.0 ;
		char data1[numsigns], data2[numsigns] ;
		char *strptr1, *strptr2 ;
		int X_file = -1 ;
		int n_file = -1 ;
		int num_id ;

		/* check files and refer the right file to either X or n */
		try
		{
			for (int i=0; i<num_Fr_flag; ++i)
			{
				strcpy( sdummy, ifo) ;
				strcat( sdummy, Fr_file[i]) ;
				if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(33) ; }

				num_id = 0 ;
				while (fgets (sdummy, numsigns, inpf) != NULL)
				{
					/* read line by line, delete white space and read next line, if a # is detected */
					strptr1 = sdummy ;
		
					while (*strptr1 == ' ') 
						++strptr1 ;
		
					if (*strptr1 == '#')
						continue ;
			
					if (*strptr1 == '%')
					{
						if ( !strncmp( sdummy, "% Identifiers", strlen(sdummy) - 1 ) )
						{
							if (fgets (sdummy, numsigns, inpf) != NULL)
							{
								num_id += 1 ;
								delallspc(sdummy) ;
								if ( !strcmp( sdummy, "X\n") ) { X_file = i ; }
								else if ( !strcmp( sdummy, "n\n") ) { n_file = i ; }
								else  { throw (i) ; }
							}
							if (fgets (sdummy, numsigns, inpf) != NULL)
							{
								num_id += 1 ;
								delallspc(sdummy) ;
								read_single_text_data_from_cfe( "_chemical_name_common", &sdummy2) ;
								if ( strncmp( sdummy, sdummy2, strlen(sdummy) - 1 ) ) { throw (i) ; }
							}
							if (fgets (sdummy, numsigns, inpf) != NULL)
							{
								num_id += 1 ;
								delallspc(sdummy) ;
								if ( mode != (int) strtol ( sdummy, NULL, 10) ) { throw (i) ; }
							}
							if (fgets (sdummy, numsigns, inpf) != NULL)
							{
								num_id += 1 ;
								delallspc(sdummy) ;
								/* in write_structure_amplitude use %.7G therefore 1e-7, actually due to
								   rounding effects 0.5e-7 would be a sharper criterion
								*/
								if ( abs( par->rs1 - strtod ( sdummy, NULL) ) > 1e-7 ) { throw (i) ; }
							}
							if (fgets (sdummy, numsigns, inpf) != NULL)
							{
								num_id += 1 ;
								delallspc(sdummy) ;
								if ( abs( par->rs2 - strtod ( sdummy, NULL) ) > 1e-7 ) { throw (i) ; }
							}
							if (fgets (sdummy, numsigns, inpf) != NULL)
							{
								num_id += 1 ;
								delallspc(sdummy) ;
								if ( par->np != (unsigned int) strtol( sdummy, NULL, 10) ) { throw (i) ; }
							}
							for (int j=0; j<3; ++j)
							{
								if (fgets (sdummy, numsigns, inpf) != NULL)
								{
									num_id += 1 ;
									delallspc(sdummy) ;
									if ( abs( par->av_pol[j] - strtod ( sdummy, NULL) ) > 1e-7 ) { throw (i) ; }
								}
							}

							for (int j=0; j<2; ++j)
							{
								if (fgets (sdummy, numsigns, inpf) != NULL)
								{
									num_id += 1 ;
									delallspc(sdummy);
									if ( mode == 0 ) { ddummy = par->av_azi[j] ; }
									else if  ( mode == 1 ) { ddummy = par->av_r3[j] ; }

									if ( abs( ddummy - strtod ( sdummy, NULL) ) > 1e-7 ) { throw (i) ; }
								}
							}
							if (fgets (sdummy, numsigns, inpf) != NULL)
							{
								num_id += 1 ;
								delallspc(sdummy) ;
								if ( nav_rot != (unsigned int) strtol ( sdummy, NULL, 10) ) { throw (i) ; }
							}

							break ;
						}
					}
				}

				/* close Fr-file */  
				fclose(inpf) ;

				if ( num_id != 12 ) { throw (-i) ; }
			}
		}
		catch (int i)
		{
			if ( log_flag ) 
			{
				if ( i < 0 )
					fprintf( logfile, "\tIn Fr_file %s file reading terminated before reading all Identifiers (num_id = %d)\n", Fr_file[-i], num_id) ;
				else
					fprintf( logfile, "\tIn Fr_file %s the Identifier no. %d does not match\n", Fr_file[i], num_id) ;

				fflush (logfile) ;
			}
			return false ;
		}

		if ( Xn >= 0 && X_file < 0 ) { if (log_flag ) { fprintf( logfile, "\tNo structure amplitude file found for X-rays\n") ; fflush (logfile) ; } ; return false ; }
		if ( Xn <= 0 && n_file < 0 ) { if (log_flag ) { fprintf( logfile, "\tNo structure amplitude file found for neutrons\n") ; fflush (logfile) ; } ; return false ; }

		/* Allocate memory */
		if ( mode == 0 ) /* s,theta, phi */
		{	
			if ( Xn >= 0 ) { FFF_X = (dcmplx ***) calloc( par->np, sizeof(dcmplx **)) ; }
			if ( Xn <= 0 ) { FFF_n = (dcmplx ***) calloc( par->np, sizeof(dcmplx **)) ; }
			
			for ( unsigned int i=0; i<par->np; ++i)
			{
				if ( Xn >= 0 )
				{
					FFF_X[i] = (dcmplx **) calloc( nav_pc, sizeof(dcmplx *)) ;
					for ( unsigned int j=0; j<nav_pc; ++j) { FFF_X[i][j] = (dcmplx *) calloc( nav_ac, sizeof(dcmplx)) ; }
				}
				if ( Xn <= 0 )
				{
					FFF_n[i] = (dcmplx **) calloc( nav_pc, sizeof(dcmplx *)) ;
					for ( unsigned int j=0; j<nav_pc; ++j) { FFF_n[i][j] = (dcmplx *) calloc( nav_ac, sizeof(dcmplx)) ; }
				}
			}
		}
		else if ( mode == 1 ) /* s,theta, only one phi, rotation angle (0,r32) */
		{
			if ( Xn >= 0 ) { FFF_X = (dcmplx ***) calloc( par->np, sizeof(dcmplx **)) ; }
			if ( Xn <= 0 ) { FFF_n = (dcmplx ***) calloc( par->np, sizeof(dcmplx **)) ; }
			
			for ( unsigned int i=0; i<par->np; ++i)
			{
				if ( Xn >= 0 )
				{
					FFF_X[i] = (dcmplx **) calloc( nav_pc, sizeof(dcmplx *)) ;
					for ( unsigned int j=0; j<nav_pc; ++j) { FFF_X[i][j] = (dcmplx *) calloc( nav_r3, sizeof(dcmplx)) ; }
				}
				if ( Xn <= 0 )
				{
					FFF_n[i] = (dcmplx **) calloc( nav_pc, sizeof(dcmplx *)) ;
					for ( unsigned int j=0; j<nav_pc; ++j) { FFF_n[i][j] = (dcmplx *) calloc( nav_r3, sizeof(dcmplx)) ; }
				}
			}
		}


		try
		{
			bool data_flag = false ;
			unsigned int ii, jj, kk ;
			int i ;

			if ( Xn >= 0 )
			{
				sprintf ( sdummy3, "%% Re(FFF) [%s] Im(FFF) [%s]\n", unit_F_X, unit_F_X) ;

				strcpy( sdummy, ifo) ;
				strcat( sdummy, Fr_file[X_file]) ;
				if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(33) ; }
	
				ii = jj = kk = 0 ;
				while (fgets (sdummy, numsigns, inpf) != NULL)
				{
					/* read line by line, delete white space and read next line, if a # is detected */
					strptr1 = sdummy ;
		
					while (*strptr1 == ' ') 
						++strptr1 ;
		
					if (*strptr1 == '#')
					{
						if ( data_flag )
						{
							if ( kk % nav_rot != 0 ) { throw(X_file) ; }
						}
						continue ;
					}
					if (*strptr1 == '%')
					{
						if ( !strncmp( sdummy, sdummy3, strlen(sdummy) - 1 ) )
							data_flag = true ;

						continue ;
					}
					if ( data_flag )
					{
						if ( kk > nav_rot - 1 )
						{
							kk = 0 ; ++jj ; 
							if ( jj > nav_pc - 1 ) { jj = 0 ; ++ii ; }
						}

						/* read FFF_X[i][j][k] */
						if ( ( strptr2 = strchr( strptr1, ' ') ) != NULL )
						{
							i = (int)(strptr2-strptr1) ;
							memmove( data1, strptr1, i) ;
							data1[i] = 0 ;

							strptr1 = strptr2 ;
							while (*strptr1 == ' ') 
								++strptr1 ;
			
							strptr2 = strchr( strptr1, '\n');
							i = (int)(strptr2-strptr1) ;
							memmove( data2, strptr1, i) ;
							data2[i] = 0 ;

							FFF_X[ii][jj][kk] = dcmplx( strtod( data1, NULL), strtod( data2, NULL) ) ;
						}
						else { throw(-X_file) ; }

						++ kk ;
					}
				}
				++jj ; ++ii ;
				/* perform completeness check */
				if ( ii != par->np || jj != nav_pc || kk != nav_rot ) { throw(X_file) ; }

				/* close Fr-file */
				fclose(inpf) ;
			}

			data_flag = false ;
			if ( Xn <= 0 )
			{
				sprintf ( sdummy3, "%% Re(FFF) [%s] Im(FFF) [%s]\n", unit_F_n, unit_F_n) ;	
	
				strcpy( sdummy, ifo) ;
				strcat( sdummy, Fr_file[n_file]) ;
				if ( ( inpf = fopen( sdummy, "r")) == NULL) { XNDIFF_ERROR(33) ; }

				ii = jj = kk = 0 ;
				while (fgets (sdummy, numsigns, inpf) != NULL)
				{
					/* read line by line, delete white space and read next line, if a # is detected */
					strptr1 = sdummy ;
		
					while (*strptr1 == ' ') 
						++strptr1 ;
		
					if (*strptr1 == '#')
					{
						if ( data_flag )
						{
							if ( kk % nav_rot != 0 ) { throw(n_file) ; }
						}
						continue ;
					}
					if (*strptr1 == '%')
					{
						if ( !strncmp( sdummy, sdummy3, strlen(sdummy) - 1 ) )
							data_flag = true ;

						continue ;
					}
					if ( data_flag )
					{
						if ( kk > nav_rot - 1 )
						{
							kk = 0 ; ++jj ; 
							if ( jj > nav_pc - 1 ) { jj = 0 ; ++ii ; }
						}

						/* read FFF_n[i][j][k] */
						if ( ( strptr2 = strchr( strptr1, ' ') ) != NULL )
						{
							i = (int)(strptr2-strptr1) ;
							memmove( data1, strptr1, i) ;
							data1[i] = 0 ;

							strptr1 = strptr2 ;
							while (*strptr1 == ' ') 
								++strptr1 ;
			
							strptr2 = strchr( strptr1, '\n');
							i = (int)(strptr2-strptr1) ;
							memmove( data2, strptr1, i) ;
							data2[i] = 0 ;

							FFF_n[ii][jj][kk] = dcmplx( strtod( data1, NULL), strtod( data2, NULL) ) ;

							/* Xn=0 mus tbe fulfilled, FFF_X has been already read */
							if ( test_XnEquiv_flag ) { FFF_n[ii][jj][kk] = FFF_X[ii][jj][kk] ; }
						}
						else { throw(-n_file) ; }

						++ kk ;
					}
				}
				++jj ; ++ii ;
				/* perform completeness check */
				if ( ii != par->np || jj != nav_pc || kk != nav_rot ) { throw(n_file) ; }

				/* close Fr-file */  
				fclose(inpf) ;
			}
		}
		catch (int i)
		{
			if (log_flag ) 
			{
				if ( i < 0 )
					fprintf( logfile, "\tMissing second column in %s\n", Fr_file[-i]) ; 
				else
					fprintf( logfile, "\tDimension mismatch while reading data in %s\n", Fr_file[i]) ;
 
 				fflush (logfile) ;
			}
			free_structure_amplitude() ;
			return false ;
		}
		/* exit normally read_structure_amplitude () */
		return true ;
	}


	/* export the structure amplitudes FFF_X and FFF_n */
	bool write_structure_amplitude(int mode = -1, char X_or_n = ' ', double rot_1 = 0.0, double rot_2 = 0.0, unsigned int nav_rot = 0, dcmplx ***FFF_Xn = NULL, char *oname = NULL )
	{
		if ( mode == -1 )
		{
			mode = par->av_mode ;
			/* The implementation of reading / writing of FFF_X/n in the case of random orientations makes little sense */
			if ( ( mode == 2 ) || ( mode == 3 ) )
			{
				if ( log_flag ) { fprintf( logfile, "\tWarning: Writing FFF_X/n for mode=%d is not supported !\n", mode) ; fflush (logfile) ; }
				return false ;
			}

			char ext[50];
			if ( Xn >= 0 )
			{
				sprintf (ext, "_FFF_X.dat") ;
				oname = (char *) realloc ( oname, ( strlen(mfname) + strlen(ext) + 1 ) * sizeof(char) ) ;
				sprintf (oname, "%s%s", mfname, ext) ;
				if ( mode == 0) { write_structure_amplitude( mode, 'X', par->av_azi[0], par->av_azi[1], nav_ac, &FFF_X[0], oname) ; }
				else if ( mode == 1) { write_structure_amplitude( mode, 'X', par->av_r3[0], par->av_r3[1], nav_r3, &FFF_X[0], oname) ; }
			}
			if ( Xn <= 0 ) 
			{
				sprintf (ext, "_FFF_n.dat") ;
				oname = (char *) realloc ( oname, ( strlen(mfname) + strlen(ext) + 1 ) * sizeof(char) ) ;
				sprintf (oname, "%s%s", mfname, ext) ;
				if ( mode == 0) { write_structure_amplitude( mode, 'n', par->av_azi[0], par->av_azi[1], nav_ac, &FFF_n[0], oname) ; }
				else if ( mode == 1) { write_structure_amplitude( mode, 'n', par->av_r3[0], par->av_r3[1], nav_r3, &FFF_n[0], oname) ; }
			}
			free(oname) ;
	
		}
		else
		{
			double cs, av_pc ;
			char sdummy[1024] ;
			char *sdummy2 ;
			FILE *outf ;

			strcpy( sdummy, ofo) ;
			strcat( sdummy, oname) ;
			if ( ( outf = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(19) ; }

			write_call_arguments( outf, "# ") ;

			fprintf( outf, "%% Identifiers\n") ;
			fprintf( outf, "%10c\n", X_or_n ) ;
			read_single_text_data_from_cfe( "_chemical_name_common", &sdummy2) ;
			fprintf( outf, "%10s\n", sdummy2 );
			fprintf( outf, "%10d\n", mode ) ;
			fprintf( outf, "%10.7G\n", par->rs1 ) ;
			fprintf( outf, "%10.7G\n", par->rs2 ) ;
			fprintf( outf, "%10d\n", par->np ) ;
			fprintf( outf, "%10.7G\n", par->av_pol[0] ) ;
			fprintf( outf, "%10.7G\n", par->av_pol[1] ) ;
			fprintf( outf, "%10.7G\n", par->av_pol[2] ) ;
			fprintf( outf, "%10.7G\n", rot_1 ) ;
			fprintf( outf, "%10.7G\n", rot_2 ) ;
			fprintf( outf, "%10d\n", nav_rot ) ;

			const char *unit ;
			if ( X_or_n == 'X' ) { unit = &unit_F_X[0] ; }
			else if ( X_or_n == 'n' ) { unit = &unit_F_n[0] ; }

			fprintf( outf, "%% Re(FFF) [%s] Im(FFF) [%s]\n", unit, unit) ;	
			cs = s1 ;
			for ( unsigned int i=0; i<par->np; ++i)
			{
				cs += ds ;

				for ( unsigned int j=0; j<nav_pc; ++j)
				{
					av_pc = av_pol_ang.at(j) ; 
					fprintf( outf, "%% s=%.10G theta=%.10G\n", cs, av_pc) ;

					for (unsigned int k=0; k<nav_rot; ++k)
					{
						fprintf( outf, "%15.10G %15.10G\n", real( FFF_Xn[i][j][k] ), imag( FFF_Xn[i][j][k] ) ) ;
					}	
				}
			}
			free(sdummy2) ;
			fclose(outf) ;
		}
		return true ;
	}


	/* compute for all types of atoms in the crystal the scattering amplitudes as a function of s, f_0(s) 
	   and save an index in structure atom */
	void compute_f0table_and_add_f0ind_to_atom()
	{
		int f0ind = 0 ;
		double cs ;

		int* atomtypeindlist = (int *) calloc( f0ind + 1 , sizeof(int)) ; /* stores index of first occurence of an atom type in atom */
		f0table = (double **) calloc( f0ind + 1 , sizeof(double *)) ; /* f0ind x par->np table containing f_0^j(s), j=0,...,f0ind */
		f0table[f0ind] = (double *) calloc( par->np , sizeof(double)) ;

		/* extra treatment of the first atom's type */
		atomtypeindlist[f0ind] = 0;
		atom[0].f0ind = f0ind ;
		cs = s1 ;
		for ( unsigned int i=0; i<par->np; ++i)
		{
			cs += ds ;
			f0table[f0ind][i] = get_atomic_scattering_factor( atom[0].type , cs ) ;
		}

		for ( unsigned int i=1; i<noa; ++i)
		{
			for (int j=0; j<=f0ind; ++j)
			{
				if (!strcmp(atom[atomtypeindlist[j]].type,atom[i].type))
				{
					atom[i].f0ind = atom[atomtypeindlist[j]].f0ind ;
					break ;
				}
				else
				{
					/* create a new entry if none of the previous types did match */
					if (j == f0ind)
					{
						++f0ind ;
						atomtypeindlist = (int *) realloc ( atomtypeindlist , (f0ind + 1) * sizeof(int)) ;
						atomtypeindlist[f0ind] = i;

						f0table = (double **) realloc ( f0table , (f0ind + 1) * sizeof(double *)) ;
						f0table[f0ind] = (double *) calloc( par->np , sizeof(double)) ;
						cs = s1 ;
						for ( unsigned int k=0; k<par->np; ++k)
						{
							cs += ds ;
							f0table[f0ind][k] = get_atomic_scattering_factor( atom[i].type , cs ) ;
						}
						atom[i].f0ind = f0ind ;

						/* break */
					}
				}
			}
		}
		numf0 = f0ind ;
		free(atomtypeindlist) ;
	}


	/* deallocate FFF_X and FFF_n, does not depend on par->av_mode */
	void free_structure_amplitude()
	{
		if ( Xn >= 0 )
		{		
			for ( unsigned int i=0; i<par->np; ++i)
			{
				for ( unsigned int j=0; j<nav_pc; ++j) { free(FFF_X[i][j]) ; }
				free(FFF_X[i]);
			}
			free(FFF_X);
		}
		if ( Xn <= 0 )
		{		
			for ( unsigned int i=0; i<par->np; ++i)
			{
				for ( unsigned int j=0; j<nav_pc; ++j) { free(FFF_n[i][j]); }
				free(FFF_n[i]);
			}
			free(FFF_n);
		}
	}


	/* Computes structure amplitude for X-ray and neutron scattering and keeps them in memory.
	   X-ray: F [nm] (r_e is now included in the calculation), Neutrons: F [nm] */
	/* Note, that it must be called after compute_eQvj(int mode)  */
	void compute_structure_amplitude(int mode = -1)
	{
		double cs;
		dcmplx FFF_cdummy;

		/* if mode = -1 isn't overwritten by a user-supplied mode parameter, use by default par->av_mode */
		if ( mode == -1 ) { mode = par->av_mode ; }

		/* use eQvj */

		double pr ; /* progress indicator */
		char rtstr[30] ; /* time string filled by remtime */
		time_t start, now ; /* timer events */

		if ( time_flag > 1 ) { time(&start) ; }

		#pragma omp parallel if ( openmp_flag )
		{
			if ( omp_in_parallel() )
			{
				#pragma omp single
				if ( log_flag ) { fprintf( logfile ,"%sParallelized Run: #-CPU=%d, #-Threads=%d\n", "\t", omp_get_num_procs(), omp_get_num_threads() ) ; fflush( logfile ) ; }
			}
			else
			{
				if ( log_flag ) { fprintf( logfile ,"%sSerialized Run\n", "\t") ; fflush( logfile ) ; }
			}
		}

		/* loop variables must be ALL declared outside the parallel block in order to get OpenMP working */
		unsigned int jj, kk, ll ;
		unsigned int ii, ii_count ;
		dcmplx FFF_X_sum, FFF_n_sum ;

		/* parallelize ii-loop or inner jj-loop has no influence on the performance !!! */

		/* s,theta,phi, all rotations (r11,r12),...,(r31,r32) are zero */
		/* Polar-Azimuthal Grid and MC mode */
		if ( ( mode == 0 ) || ( mode == 2 ) || ( mode == 3 ) )
		{
			if ( Xn >= 0 ) { FFF_X = (dcmplx ***) calloc( par->np, sizeof(dcmplx **)) ; }
			if ( Xn <= 0 ) { FFF_n = (dcmplx ***) calloc( par->np, sizeof(dcmplx **)) ; }
			
			for ( unsigned int i=0; i<par->np; ++i)
			{
				if ( Xn >= 0 )
				{
					FFF_X[i] = (dcmplx **) calloc( nav_pc, sizeof(dcmplx *)) ;
					for ( unsigned int j=0; j<nav_pc; ++j)
						FFF_X[i][j] = (dcmplx *) calloc( nav_ac, sizeof(dcmplx)) ;

				}
				if ( Xn <= 0 )
				{
					FFF_n[i] = (dcmplx **) calloc( nav_pc, sizeof(dcmplx *)) ;
					for ( unsigned int j=0; j<nav_pc; ++j)
						FFF_n[i][j] = (dcmplx *) calloc( nav_ac, sizeof(dcmplx)) ;
				}
			}

			if ( !FNoMem_flag )
			{
				ii_count = 0 ;
				#pragma omp parallel private( jj, kk, ll, cs, FFF_cdummy, FFF_X_sum, FFF_n_sum) shared( ii_count, pr, now, start, rtstr, stdout ) default(none) if ( openmp_flag )
				{
					#pragma omp for schedule(static)
					for ( ii=0; ii<par->np; ++ii)
					{
						/* reduction(+:ii_count) clause does not work, must be shared */
						ii_count += 1 ;

						cs = s1 + (double)( ii + 1 ) * ds ;

						/* set FFF to zero, usually not necessary cause of calloc instead of malloc */
						for ( jj=0; jj<nav_pc; ++jj)
						{
							for ( kk=0; kk<nav_ac; ++kk)
							{
								if ( Xn >= 0 ) { FFF_X[ii][jj][kk] = dcmplx(0.0,0.0) ; }
								if ( Xn <= 0 ) { FFF_n[ii][jj][kk] = dcmplx(0.0,0.0) ; }
							}
						}

						if ( time_flag > 1 )
						{
							// /* all threads update pr while computing FFF */
							// #pragma omp critical
							/* as for single thread allow update in parallel mode only for master thread */
							if ( omp_get_thread_num() == 0 )
							{
								pr = 100.0*((double) ii_count)/((double) par->np ) ;
								time(&now) ;
								remtime( (long int)(difftime( now, start)), pr, &rtstr[0]) ;
								fprintf( stdout ,"PID %d | Computing F | %5.1f %% | %-30s\r", pid, pr, rtstr) ; fflush( stdout ) ;
							}
						}

						for ( jj=0; jj<nav_pc; ++jj)
						{
							for ( kk=0; kk<nav_ac; ++kk)
							{
								if ( Xn >= 0 ) { FFF_X_sum = dcmplx( 0.0, 0.0) ; }
								if ( Xn <= 0 ) { FFF_n_sum = dcmplx( 0.0, 0.0) ; }

								for ( ll=0; ll<noa; ++ll)
								{
									FFF_cdummy = exp (dcmplx( 0.0, 2.0 * M_PI * cs * ( atom[ll].coord[0] * eQvj[jj][kk][0] + atom[ll].coord[1] * eQvj[jj][kk][1] + atom[ll].coord[2] * eQvj[jj][kk][2] ) ) ) ;
									FFF_cdummy *= atom[ll].mult ;

									/* compute structure amplitude F_k for X-ray with the atomic scattering factors */
									if ( Xn >= 0 ) { FFF_X_sum += FFF_cdummy * f0table[atom[ll].f0ind][ii] * ( (dcmplx) sc_rho ) ; } /* [nm] */
									/* the same for neutrons with the coherent neutron scattering lengths
									   don't use any minus sign convention for b_coh, use just b_coh,
									   the sld will also not include later minus signs
									*/
									if ( Xn <= 0 ) { FFF_n_sum += FFF_cdummy * atom[ll].coh_sc_length * ( (dcmplx) sc_scl ) ; } /* [nm] */

									if ( test_XnEquiv_flag ) { FFF_n_sum = FFF_X_sum ; }
								}
								#pragma omp critical
								{
									if ( Xn >= 0 ) { FFF_X[ii][jj][kk] += FFF_X_sum ; } /* [nm] */
									if ( Xn <= 0 ) { FFF_n[ii][jj][kk] += FFF_n_sum ; } /* [nm] */
								}
							}
						}
					}
				}
				if ( time_flag > 1 ) 
				{
					fprintf( stdout, "PID %d | Computing F | %-40s\r", pid, "Done") ; fflush (stdout) ;
				}
			}
		}
		else if ( mode == 1 ) /* s, theta, only one phi, rotation angle (0,r32) */
		{
			/* use eQvj */

			if ( Xn >= 0 ) { FFF_X = (dcmplx ***) calloc( par->np, sizeof(dcmplx **)) ; }
			if ( Xn <= 0 ) { FFF_n = (dcmplx ***) calloc( par->np, sizeof(dcmplx **)) ; }

			for ( unsigned int i=0; i<par->np; ++i)
			{
				if ( Xn >= 0 )
				{
					FFF_X[i] = (dcmplx **) calloc( nav_pc, sizeof(dcmplx *)) ;
					for ( unsigned int j=0; j<nav_pc; ++j) { FFF_X[i][j] = (dcmplx *) calloc( nav_r3, sizeof(dcmplx)) ; }
				}
				if ( Xn <= 0 )
				{
					FFF_n[i] = (dcmplx **) calloc( nav_pc, sizeof(dcmplx *)) ;
					for ( unsigned int j=0; j<nav_pc; ++j) { FFF_n[i][j] = (dcmplx *) calloc( nav_r3, sizeof(dcmplx)) ; }
				}
			}

			if ( !FNoMem_flag )
			{
				ii_count = 0 ;
				#pragma omp parallel private( jj, kk, ll, cs, FFF_cdummy, FFF_X_sum, FFF_n_sum) shared( ii_count, pr, now, start, rtstr, stdout ) default(none) if ( openmp_flag )
				{
					#pragma omp for schedule(static)
					for ( ii=0; ii<par->np; ++ii)
					{
						cs = s1 + (double)( ii + 1 ) * ds ;

						/* set FFF to zero, usually not necessary cause of calloc instead of malloc */
						for ( jj=0; jj<nav_pc; ++jj)
						{
							for ( kk=0; kk<nav_ac; ++kk)
							{
								if ( Xn >= 0 ) { FFF_X[ii][jj][kk] = dcmplx(0.0,0.0) ; }
								if ( Xn <= 0 ) { FFF_n[ii][jj][kk] = dcmplx(0.0,0.0) ; }
							}
						}

						if ( time_flag > 1 )
						{
							/* as for single thread allow update in parallel mode only for master thread */
							if ( omp_get_thread_num() == 0 )
							{
								pr = 100.0*((double) ii_count)/((double) par->np ) ;
								time(&now) ;
								remtime( (long int)(difftime( now, start)), pr, &rtstr[0]) ;
								fprintf( stdout ,"PID %d | Computing F | %5.1f %% | %-30s\r", pid, pr, rtstr) ; fflush( stdout ) ;
							}
						}

						for ( jj=0; jj<nav_pc; ++jj)
						{
							/* fixed av_ac value */
							for ( kk=0; kk<nav_r3; ++kk)
							{
								if ( Xn >= 0 ) { FFF_X_sum = dcmplx( 0.0, 0.0) ; }
								if ( Xn <= 0 ) { FFF_n_sum = dcmplx( 0.0, 0.0) ; }

								for ( ll=0; ll<noa; ++ll)
								{
									/* use eQvj */
									FFF_cdummy = exp (dcmplx( 0.0, 2.0 * M_PI * cs * ( atom[ll].coord[0] * eQvj[jj][kk][0] + atom[ll].coord[1] * eQvj[jj][kk][1] + atom[ll].coord[2] * eQvj[jj][kk][2] ) ) ) ;
									FFF_cdummy *= atom[ll].mult ;

									/* compute structure amplitude F_k for X-ray with the atomic scattering factors */
									if ( Xn >= 0 ) { FFF_X_sum += FFF_cdummy * f0table[atom[ll].f0ind][ii] * ( (dcmplx) sc_rho ) ; } /* [nm] */

									/* the same for neutrons with the coherent neutron scattering lengths
									   don't use any minus sign convention for b_coh, use just b_coh,
									   the sld will also not include later minus signs
									*/
									if ( Xn <= 0 ) { FFF_n_sum += FFF_cdummy * atom[ll].coh_sc_length * ( (dcmplx) sc_scl ) ; } /* [nm] */

									if ( test_XnEquiv_flag ) { FFF_n_sum = FFF_X_sum ; }
								}
								#pragma omp critical
								{
									if ( Xn >= 0 ) { FFF_X[ii][jj][kk] += FFF_X_sum ; } /* [nm] */
									if ( Xn <= 0 ) { FFF_n[ii][jj][kk] += FFF_n_sum ; } /* [nm] */
								}
							}
						}
						/* reduction(+:ii_count) clause does not work, must be shared */
						ii_count += 1 ;
					}
				}
				if ( time_flag > 1 ) 
				{
					fprintf( stdout, "PID %d | Computing F | %-40s\r", pid, "Done") ; fflush (stdout) ;
				}
			}
		}
	}

	/* compute array of to-be-used polar and azimuthal angles for the Powder Average */
	/* also computes nav_pc, av_p1, av_p2, av_pd and nav_ac, av_a1, av_a2, av_ad as well as nav_r3, av_r3d and nav_MC and s1, ds */
	void compute_pol_azi_ang()
	{
		double av_p1, av_p2, av_pd, av_a1, av_a2, av_ad ;
		double av_pc, av_ac ;

		if ( par->av_mode == 0 )
		{
			/* theta, phi grid */
			av_p1 = par->av_pol[0] ;
			av_p2 = par->av_pol[1] ;
			av_pd = par->av_pol[2] ;

			av_a1 = par->av_azi[0] ;
			av_a2 = par->av_azi[1] ;
			av_ad = par->av_azi[2] ;

			/* av_pc starts effectively with av_p1 + av_pd/2.0 and similar for av_ac */
			if ( ( av_p2 > av_p1 ) && ( av_pd > 0.0 ) )
			{
				nav_pc = 0 ;
				av_pc = av_p1 + av_pd/2.0 ;
				do { ++nav_pc; av_pol_ang.push_back(av_pc) ; av_pc += av_pd ; } while ( !(av_pc > av_p2) ) ;
			}
			else if ( ( av_p2 == av_p1 ) && ( av_pd == 0.0 ) )
			{
				nav_pc = 1 ;
				av_pol_ang.push_back(av_p1) ;
			}

			if ( ( av_a2 > av_a1 ) && ( av_ad > 0.0 ) )
			{
				nav_ac = 0 ;
				av_ac = av_a1 + av_ad/2.0 ;
				do { ++nav_ac; av_azi_ang.push_back(av_ac) ; av_ac += av_ad ; } while ( !(av_ac > av_a2) ) ;
			}
			else if ( ( av_a2 == av_a1 ) && ( av_ad == 0.0 ) )
			{
				nav_ac = 1 ;
				av_azi_ang.push_back(av_ac) ;
			}
		}
		else if ( par->av_mode == 1 )
		{
			/* theta grid, only one phi, random r3 rotation */
			nav_r3 = par->n_r3 ;
			av_r3d = (par->av_r3[1]-par->av_r3[0])/((double)(nav_r3-1)) ;
			
			av_p1 = par->av_pol[0] ;
			av_p2 = par->av_pol[1] ;
			av_pd = par->av_pol[2] ;

			av_a1 = par->av_azi[0] ;
			av_a2 = par->av_azi[1] ;
			av_ad = par->av_azi[2] ;

			/* av_pc starts effectively with av_p1 + av_pd/2.0 and similar for av_ac */
			if ( ( av_p2 > av_p1 ) && ( av_pd > 0.0 ) )
			{
				nav_pc = 0 ;
				av_pc = av_p1 + av_pd/2.0 ;
				do { ++nav_pc; av_pol_ang.push_back(av_pc) ; av_pc += av_pd ; } while ( !(av_pc > av_p2) ) ;
			}
			else if ( ( av_p2 == av_p1 ) && ( av_pd == 0.0 ) )
			{
				nav_pc = 1 ;
				av_pol_ang.push_back(av_p1) ;
			}
			
			nav_ac = 1 ;
			av_azi_ang.push_back( av_a1 + av_ad/2.0 ) ;
		}
		else if ( ( par->av_mode == 2 ) || ( par->av_mode == 3 ) )
		{
			double x1, x2, x, y, z ;

			/* MC mode */
		
			if ( !init_MTRand_userdef )
			{
				/* initialize the PRNG with CPU time, later d_r() calls will provide pseudo random numbers within the range [0,1) */
				time_t seconds = time (NULL) ;
				init_MTRand = (unsigned int)seconds ;
				if ( log_flag ) { fprintf( logfile, "\tInitializing MTRand PRNG for MC Powder Average with seed %u\n", init_MTRand) ; fflush (logfile) ; }
			}
			else
			{
				/* initialize the PRNG with CPU time, later d_r() calls will provide pseudo random numbers within the range [0,1) */
				if ( log_flag ) { fprintf( logfile, "\tInitializing MTRand PRNG for MC Powder Average with user-defined seed %u\n", init_MTRand) ; fflush (logfile) ; }
			}
			MTRand d_r( init_MTRand ) ;
			if ( log_flag ) { fprintf( logfile, "\tdone\n") ; fflush (logfile) ; }

			nav_MC = par->n_MC ;
			av_pol_ang.resize( nav_MC ) ;
			av_azi_ang.resize( nav_MC ) ;

			for ( unsigned int j=0; j<av_pol_ang.size(); ++j)
			{
				/* apply only simple restrictions for range of polar and azimuthal angle */
				do
				{
					if ( par->av_mode == 2 )
					{
						/* MC ArcCos */
						/* In principle instead of using the simple restriction method by just dicing until the restrictions are fulfilled, 
						   for the ArcCos method there is a more efficient scheme:

						   For theta1=par->av_pol[0] < theta2=par->av_pol[1] and phi1=par->av_azi[0] < phi2=par->av_azi[1] use:
						   av_pc = RD *  acos( cos(theta2) + ( cos(theta1) - cos(theta2) ) * d_r() )
						   and
						   av_ac = phi1 + ( phi2 - phi1 ) * d_r()

						   This reduces for theta1=0 and theta2=pi/2 to the known case:
						   av_pc = RD *  acos( d_r() )
						   Similarly for theta1=0 and theta2=pi to:
						   av_pc = RD *  acos( 2 * d_r() - 1 )
						   And for phi1=0 and phi2=2*pi:
						   av_ac = 360.0 * d_r() ;
						*/
						for ( unsigned int j=0; j<av_pol_ang.size(); ++j)
						{
							av_pc = RD *  acos( d_r() ) ;
							av_ac = 360.0 * d_r() ;
						}
					}
					else /* par->av_mode == 3 */
					{
						/* MC Marsaglia */
						/* WARNING: The criteria ( x1 * x1 + x2 * x2 ) >= 0.5 to assure that z>0 instead of ( x1 * x1 + x2 * x2 ) >= 1.0 is applied here !!!
						            If the general case (full instead of only the upper half sphere) shall be considered replace 0.5 by 1.0 !!!
						*/
						/* To design x1 and x2 of the kind x1 = a_1 + b_1 * x, x2 = a_2 + b_2 * y such that they fulfill:
						   phi_1 <= ArcTan[x2/x1] ( + pi ) <= phi2
						   and
						   theta_1 <= ArcCos[1-2*( x1^2 + x2^2 )] < theta_2
						   with at least
						   x1^2 + x2^2 <= 1
						   holding, and of course x,y ~ U_[0,1)
						   seems not to be such an easy task to solve !!
						*/
						do 
						{
							x1 = 2.0 * d_r() - 1.0 ;
							x2 = 2.0 * d_r() - 1.0 ;
						} while ( ( x1 * x1 + x2 * x2 ) >= 0.5 ) ;
						
						x = 2.0 * x1 * sqrt( 1.0 - x1 * x1 - x2 * x2 ) ;
						y = 2.0 * x2 * sqrt( 1.0 - x1 * x1 - x2 * x2 ) ;
						z = 1.0 - 2.0 *  ( x1 * x1 + x2 * x2 ) ;

						// r = sqrt( x * x + y * y + z * z ) ; r == 1 by construction -> skip r

						av_pc = RD * acos(z) ;

						av_ac = atan2 ( y, x) ;
						if ( av_ac < 0 ) { av_ac += 2 * M_PI ; } 
						av_ac *= RD ;
					}
				} while ( ( av_pc < par->av_pol[0] ) || ( av_pc > par->av_pol[1] ) || ( av_ac < par->av_azi[0] ) || ( av_ac > par->av_azi[1] ) ) ;
				
				av_pol_ang.at(j) = av_pc ;
				av_azi_ang.at(j) = av_ac ;
			}
			nav_pc = nav_MC ;
			nav_ac = 1 ;
		}

		/* for the s-range */
		ds = (par->rs2 - par->rs1) / ((double)(par->np-1)) ;
		s1 = par->rs1-ds ;

		/* if -distr flag is set export the angles */
		if ( distr_flag )
		{
			FILE *outf ;
			char oname[1024]; /* name for output file */
			/* clear .distr file */
			sprintf( oname, "%s%s.pol_azi_ang", ofo, mfname) ;
			if ( ( outf = fopen( oname, "w")) == NULL) { XNDIFF_ERROR(18) ; }

			fprintf( outf, "# MTRand PRN\n" ) ;
			fprintf( outf, "# av_mode: %d\n", par->av_mode ) ;
			fprintf( outf, "# nav_pc: %d\n", nav_pc ) ;
			fprintf( outf, "# nav_ac: %d\n", nav_ac ) ;
			if ( par->av_mode == 1 ) { fprintf( outf, "# nav_r3: %d\n", nav_r3 ) ; }
			if ( ( par->av_mode == 2 ) || ( par->av_mode == 3 ) ) { fprintf( outf, "# nav_MC: %d\n", nav_MC ) ; }
			fprintf( outf, "# index %15s %15s\n", var_name(av_pol_ang), var_name(av_azi_ang)) ;

			for ( unsigned int j=0; j<max(av_pol_ang.size(),av_azi_ang.size()); ++j)
			{
				if ( j<min(av_pol_ang.size(),av_azi_ang.size()) )
				{
					fprintf( outf, "%7u %15.10G %15.10G\n", j, av_pol_ang.at(j), av_azi_ang.at(j)) ;
				}
				else if ( j<av_pol_ang.size() && j>av_azi_ang.size() )
				{
					fprintf( outf, "%7u %15.10G\n", j, av_pol_ang.at(j)) ;
				} 
				else if ( j>av_pol_ang.size() && j<av_azi_ang.size() )
				{
					fprintf( outf, "%7u %15s %15.10G\n", j, "", av_azi_ang.at(j)) ;
				} 
			}
			fclose(outf) ;
		}
	}




	/* compute all scalar products between the three a_j's and all unit vectors eQ */
	void compute_eQvj(int mode = -1 )
	{
		/* if mode = -1 isn't overwritten by a user-supplied mode parameter, use by default par->av_mode */
		if ( mode == -1 ) { mode = par->av_mode ; }

		double av_pc = 0.0, av_ac = 0.0, av_r3c = 0.0 ;
		double ddummy1,ddummy2,ddummy3,ddummy4,ddummy5;
		double eQ[3],v1[3],v2[3],v3[3];
		double sinav_ac = 0.0,cosav_ac = 0.0 ;
		double sinav_pc = 0.0,cosav_pc = 0.0 ;

		if ( ( mode == 0 ) || ( mode == 2 ) || ( mode == 3 ) )
		{
			/* theta, phi, all rotations (r11,r12),...,(r31,r32) are zero (r1=r2=r3=0) !!! */
			eQvj = (double ***) calloc( nav_pc, sizeof(double **)) ;	
			for ( unsigned int j=0; j<nav_pc; ++j)
			{
				av_pc = av_pol_ang.at(j) ;
				
				eQvj[j] = (double **) calloc( nav_ac, sizeof(double *)) ;

				sinav_pc=sin(DR * av_pc) ;
				cosav_pc=cos(DR * av_pc) ;

				for ( unsigned int k=0; k<nav_ac; ++k)
				{
					if ( mode == 0 ) { av_ac = av_azi_ang.at(k) ; }
					else if  ( ( mode == 2 ) || ( mode == 3 ) ) { av_ac = av_azi_ang.at(j) ; }
					
					eQvj[j][k] = (double *) calloc( 3, sizeof(double)) ;

					ddummy1 = sinav_pc * cos(DR * av_ac) ;
					ddummy2 = sinav_pc * sin(DR * av_ac) ;
					ddummy3 = cosav_pc ;

					eQ[0] = ddummy1 * kg1[0] + ddummy2 * kg2[0] + ddummy3 * G0[0] ;
					eQ[1] = ddummy1 * kg1[1] + ddummy2 * kg2[1] + ddummy3 * G0[1] ;
					eQ[2] = ddummy1 * kg1[2] + ddummy2 * kg2[2] + ddummy3 * G0[2] ;
					ddummy1 = betrag (eQ) ;
					if (ddummy1 > 1e-10)
					{
						eQ[0] /= ddummy1 ;
						eQ[1] /= ddummy1 ;
						eQ[2] /= ddummy1 ;
					}

					eQvj[j][k][0]=sp(eQ,a1);
					eQvj[j][k][1]=sp(eQ,a2);
					eQvj[j][k][2]=sp(eQ,a3);
				}
			}
		}
		else if ( mode == 1 )
		{
			/* theta, only one phi, rotation angle r3!=0 */
			/* fixed av_ac value */
			av_ac = av_azi_ang.at(0) ;
			sinav_ac = sin (DR * av_ac) ;
			cosav_ac = cos (DR * av_ac) ;

			eQvj = (double ***) calloc( nav_pc, sizeof(double **)) ;
			for ( unsigned int j=0; j<nav_pc; ++j)
			{
				av_pc = av_pol_ang.at(j) ;
				sinav_pc = sin (DR * av_pc) ;
				cosav_pc = cos (DR * av_pc) ;

				av_r3c = par->av_r3[0] - av_r3d ;

				eQvj[j] = (double **) calloc( nav_r3, sizeof(double *)) ;
				for ( unsigned int k=0; k<nav_r3; ++k)
				{
					av_r3c += av_r3d ;
					eQvj[j][k] = (double *) calloc( 3, sizeof(double)) ;

					ddummy1 = sinav_pc * cosav_ac ;
					ddummy2 = sinav_pc * sinav_ac ;
					ddummy3 = cosav_pc ;
					eQ[0] = ddummy1 * kg1[0] + ddummy2 * kg2[0] + ddummy3 * G0[0] ;
					eQ[1] = ddummy1 * kg1[1] + ddummy2 * kg2[1] + ddummy3 * G0[1] ;
					eQ[2] = ddummy1 * kg1[2] + ddummy2 * kg2[2] + ddummy3 * G0[2] ;
					ddummy1 = betrag (eQ) ;
					if (ddummy1 > 1e-10)
					{
						eQ[0] /= ddummy1 ;
						eQ[1] /= ddummy1 ;
						eQ[2] /= ddummy1 ;
					}

					/* Rotate crystal k on Gs -> e1,e2,e3 */
					/* e_l=<a_l,G_s>G_s+(a_l-<a_l,G_s>G_s)*cos(r3[k])+(a_l x G_s)sin(r3[k]) for l=1,2,3 */
					/* a_l are the lattice vectors of length a_l */
					/* G_s unit vector */
					/* e_l are rotated vectors of length a_l */
					/* due to r1,r2=0 it holds v_l=e_l */
					ddummy1 = Gs[0] ;
					ddummy2 = Gs[1] ;
					ddummy3 = Gs[2] ;
					ddummy4 = cos(DR*av_r3c) ;
					ddummy5 = sin(DR*av_r3c) ;
					v1[0] = a1[0] * (ddummy4+SQUARE(ddummy1)*(1.0-ddummy4)) +
							a1[1] * (ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a1[2] * (-ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) ;
					v1[1] = a1[0] * (-ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a1[1] * (ddummy4+SQUARE(ddummy2)*(1.0-ddummy4)) +
							a1[2] * (ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) ;
					v1[2] = a1[0] * (ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) +
							a1[1] * (-ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) +
							a1[2] * (ddummy4+SQUARE(ddummy3)*(1.0-ddummy4)) ;
					v2[0] = a2[0] * (ddummy4+SQUARE(ddummy1)*(1.0-ddummy4)) +
							a2[1] * (ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a2[2] * (-ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) ;
					v2[1] = a2[0] * (-ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a2[1] * (ddummy4+SQUARE(ddummy2)*(1.0-ddummy4)) +
							a2[2] * (ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) ;
					v2[2] = a2[0] * (ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) +
							a2[1] * (-ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) +
							a2[2] * (ddummy4+SQUARE(ddummy3)*(1.0-ddummy4)) ;
					v3[0] = a3[0] * (ddummy4+SQUARE(ddummy1)*(1.0-ddummy4)) +
							a3[1] * (ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a3[2] * (-ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) ;
					v3[1] = a3[0] * (-ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a3[1] * (ddummy4+SQUARE(ddummy2)*(1.0-ddummy4)) +
							a3[2] * (ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) ;
					v3[2] = a3[0] * (ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) +
							a3[1] * (-ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) +
							a3[2] * (ddummy4+SQUARE(ddummy3)*(1.0-ddummy4)) ;

					eQvj[j][k][0]=sp(eQ,v1);
					eQvj[j][k][1]=sp(eQ,v2);
					eQvj[j][k][2]=sp(eQ,v3);
				}
			}
		}
	}




	int stackcpp ()
	{
		FILE* outf ; /* FILE pointer for output files */
		char oname[1024]; /* names for output files */

		double ddummy1, ddummy2, ddummy3, ddummy4, ddummy5 ;
		double sinav_pc = 0.0, cosav_pc = 0.0, sinav_ac = 0.0, cosav_ac = 0.0 ;

		const gsl_rng_type* GSL_RNG_T ;
		gsl_rng* GSL_RNG ;
		gsl_ran_discrete_t* d_rand ;
		double* nis_arr ;

		int idummy;
		char sdummy[1024], sdummy2[1024];

		dcmplx exp_N1Qv1, exp_N2Qv2, exp_Qv1, exp_Qv2, exp_Qv3 ;
		dcmplx exp_N1Qw1, exp_N2Qw2, exp_N1Qwt1, exp_N2Qwt2 ;
		dcmplx *exp_N3Qv3, *exp_N3Qw3, *exp_N3Qwt3 ;
		dcmplx *exp_RQ ;
		dcmplx exp_DRQ ;
		unsigned int k_phi ;
		
		double piQv1, piQv2, piQv3 ;
		double piQwt1, piQwt2, piQwt3 ;
		double piQw1, piQw2, piQw3 ;

		double G[3] ;					/* Reciprocal lattice unit vector parallel to (s-s0) */

		dcmplx **GA ;					/* lattice amplitude GA */
		dcmplx GA_N1_N2 ;

		dcmplx Bg_X ;					/* electric field amplitude for the dispersion medium */
		dcmplx Bg_n ; 

		double ***S_X ;					/* nspsp x rho/sld_multi.n x par->np, pointer to scattering patterns S */
		double ***S_n ;

		dcmplx ***dE_X ;				/* par->nsp x par->nms x rho/sld_multi.n, pointer to incremental scattered electric field dE */
		dcmplx ***dE_n ;

		dcmplx dEstack_X ;				 /* dummy variables for adding up scattering in a stack */
		dcmplx dEstack_n ;
		double dBstack_n ;

		double ****dS_X ;				/* nspsp x nms/1 x rho/sld_multi.n x par->np pointer to incremental scattering intensity dS */
		double ****dS_n ;				/* single particles 0...par->nsp-1, stacks par->nsp...nspsp-1 */ 
								/* single particles 0...par->nms, stacks 0 */   
								/* rho/sld_multi.n applies for dS_X/dS_n respectively*/
		double ****dS_X_parallel_thread ;
		double ****dS_n_parallel_thread ;

		double ***dB_n ;				/* nspsp x nms/1 x sld_multi.n, no par->np dependence necessary */
		double **B_n ;					/* nspsp x sld_multi.n */

		double ***Yc_X ;                                /* par->nsp x 10 x par->np, Y[m][0...9][i] pointer to the 10 coherent structure parts */
		double ***Yc_X_parallel_thread ;
		double dYc_X_stack[10] ;
		double ***Yc_n ;				/* par->nsp x 10 x par->np, Y[m][0...9][i] pointer to the 10 coherent structure parts */
		double ***Yc_n_parallel_thread ;
		double dYc_n_stack[10] ;
		double **Yi_n ;					/* par->nsp x 4, Y[m][q] pointer to the 4 incoherent structure parts */

		dcmplx **P ;					/* par->nsp x par->nms, pointer for crystal term V * sc_rho/sld * "cdummy" */
		dcmplx **P_isl ;				/* par->nsp x par->nms, pointer for crystal+isl term V_isl * sc_rho/sld * "cdummy" */
		dcmplx **P_osl ;				/* par->nsp x par->nms, pointer for crystal+isl+osl term V_osl * sc_rho/sld * "cdummy" */

		dcmplx cdummy11, cdummy12, cdummy13 ;		/* Complex dummy variables */
		dcmplx cdummy21, cdummy22, cdummy23 ;		/* Complex dummy variables */
		dcmplx cdummy31, cdummy32, cdummy33 ;		/* Complex dummy variables */
		dcmplx *cdummy41, *cdummy42, *cdummy43;		/* Complex dummy variables */
		dcmplx FFF_cdummy ;				/* Complex dummy variables */

		double l1, l2, l3 ;				/* variables for calculations of dispersion medium and stabilizer scattering */
		double ol1, ol2, ol3 ;				/* variables for calculations of dispersion medium and stabilizer scattering */
		double *kf1, *kf2, *kf3;
		double *okf1, *okf2, *okf3;

		unsigned int ii ;
		unsigned int jj, kk, ll, mm, pp, qq ;
		unsigned int kk_jj_count ;

		double cs, s[3] ;				/* current s value cs, vector Q ~ (s-s0) */

		unsigned int* n1 ;           /* Number of unit cells in a1 direction of ith crystal */
		unsigned int* n2 ;           /* Number of unit cells in a2 direction of ith crystal */
		unsigned int* n3 ;           /* Number of unit cells in a3 direction of ith crystal */
		double* d1 ;             /* Displacement along k1 in nm */
		double* d2 ;             /* Displacement along k2 in nm */
		double* r3 ;             /* Rotation in degree */
		double* D ;              /* Distance to next crystal along Gs in nm */
		double** R ;           /* Vectors to the crystals origins */
		double DMIN ; /* variable for minimum spacing D between two particles */

		double** v1 ;          /* Unit cell vector a1 of ith crystal */
		double** v2 ;          /* Unit cell vector a2 of ith crystal */
		double** v3 ;          /* Unit cell vector a3 of ith crystal */

		double** w1 ;          /* Unit cell vector a1 of ith crystal including inner stabilizer layer */
		double** w2 ;          /* Unit cell vector a2 of ith crystal including inner stabilizer layer */
		double*** w3 ;          /* Unit cell vector a3 of ith crystal including inner stabilizer layer */

		double** wt1 ;         /* Unit cell vector a1 of ith crystal including total stabilizer layer */
		double** wt2 ;         /* Unit cell vector a2 of ith crystal including total stabilizer layer */
		double*** wt3 ;         /* Unit cell vector a3 of ith crystal including total stabilizer layer */

		double** V ; 		/* V(crystal) = V_uc, unit cell of crystal [nm^3] */
		double** V_isl ;        /* V_isl = V(crystal) + V(isl), unit cell of crystal and inner stabilizer layer [nm^3] */
		double** V_osl ;        /* V_osl = V(crystal)+V(isl+osl), unit cell of crystal and total stabilizer layer [nm^3] */

		double** VN ; 		/* par->nsp x par->nms, V(crystal) = (n1*n2*n3) * V_uc, crystal volume [nm^3] */
		double** VN_isl ;        /* par->nsp x par->nms, V_isl = V(crystal) + V(isl), volume of crystal and inner stabilizer layer [nm^3] */
		double** VN_osl ;        /* par->nsp x par->nms, V_osl = V(crystal)+V(isl+osl), volume of crystal and total stabilizer layer [nm^3] */

		double** VN_irr ; /* par->nsp x par->nms + par->nms x 1, irradiated sample volume for single particle and stacks */
		double** VN_dm ; /* par->nsp x par->nms + par->nms x 1, proportional volume of the dispersion medium for single particle and stacks */

		double *Vtot_cry ; /* total crystals volume, nspsp pointer */
		double *Vtot_isl ; /* total isl volume, nspsp pointer */
		double *Vtot_osl ; /* total osl volume, nspsp pointer */

		double *Vtot_dm ; /* total dm volume, nspsp pointer, derived from V_tot with conc_cry and conc_dm */
		double *Vtot_irr ; /* total irradiated sample volume, nspsp pointer, derived from V_tot with conc_cry */

		double *drerho_dm_osl ; /* stores r_e * ( rho_dm - rho_osl ) in [1/nm^2] */
		double *drerho_osl_isl ; /* stores r_e * ( rho_osl - rho_isl ) in [1/nm^2] */
		double *drerho_isl ; /* stores r_e * rho_isl in [1/nm^2] */

		double *dsld_dm_osl ; /* stores ( sld_dm - sld_osl ) in [1/nm^2] */
		double *dsld_osl_isl ; /* stores ( sld_osl - sld_isl ) in [1/nm^2] */
		double *dsld_isl ; /* stores sld_isl in [1/nm^2] */

		/* polar averaging angles */
		double av_pc = 0.0 ;
		/* azimuthal averaging angles */
		double av_ac = 0.0  ;
		/* weighting factor */
		double weight_factor, sum_weight_factor, single_sum_weight_factor ;

		unsigned int* ind_r3; /* index of r3[] on the r3-grid */

		unsigned int nspsp ;
		unsigned int* MULTIPLICITY ;
		const char presp[64] = "# " ;

		time_t start, now, last; /* timer events */

		char extensiondecision[1024];
		char extensionrepeats[1024];
		unsigned int repeatsincrease;


		/* get process ID and start timer */
		pid = getpid () ;
		if ( log_flag ) { fprintf( logfile, "PID = %d\n\n", pid) ; fflush (logfile) ; }

		if ( time_flag > 0 ) { time(&start) ; }

		if ( log_flag ) { fprintf( logfile, "Use distribution types: td1=%d, td2=%d\n", par->td1, par->td2) ; fflush (logfile) ; }

		/* initialize discrete PRNG for par->td2==2 */
		if ( par->td2 == 2 )
		{
			/* setup uniform PRNG */
			gsl_rng_env_setup() ;
			GSL_RNG_T = gsl_rng_default ;
			GSL_RNG = gsl_rng_alloc(GSL_RNG_T) ;

			if ( log_flag ) { fprintf( logfile, "Initializing d_rand PRNG (n3).\n") ; fflush (logfile) ; }

			if ( !init_d_rand_userdef )
			{
				/* default random initialization */
				if ( log_flag ) { fprintf( logfile, "\tRandom Number Generator Initialization with seed: %lu\n", gsl_rng_default_seed) ; fflush (logfile) ; }
			}
			else
			{
				/* use user supplied initialization number */
				if ( log_flag ) { fprintf( logfile, "\tRandom Number Generator Initialization with user-defined seed: %lu\n", init_d_rand) ; fflush (logfile) ; }
				gsl_rng_set( GSL_RNG, init_d_rand ) ;
			}

			// if cis provided, normalize cis to 1 and compute nis
			// here we use the approach n_i -> c_i / i, which is strictly only true if disl=dosl=0
			// it might be worth to implement n_i -> c_i / ( i * dhkl + 2 * (disl+dosl) ) to account for platelet volumes includinmg the stabilzer shells, since the cis are related to VOSL in case of VOSL normalization
			if( cis_isdef )
			{ 
				if ( log_flag )
				{
					fprintf( logfile, "\tcis = ( ") ;
					for( unsigned int i=0; i<par->nsp-1; ++i){ fprintf( logfile, " %.3lf,", cis[i]) ; }
					fprintf( logfile, " %.3lf)\n", cis[par->nsp-1]) ;
				}
			
				/* normalize cis to 1 */
				ddummy1 = 0.0 ;
				for( unsigned int i=0; i<par->nsp; ++i) { ddummy1 += cis[i] ; }
				for( unsigned int i=0; i<par->nsp; ++i) { cis[i] /= ddummy1 ; }

				if ( log_flag )
				{
					fprintf( logfile, "\tcis = ( ") ;
					for( unsigned int i=0; i<par->nsp-1; ++i){ fprintf( logfile, " %.3lf,", cis[i]) ; }
					fprintf( logfile, " %.3lf) (normalized to 1)\n", cis[par->nsp-1]) ;
				}

				/* calculate unnormalized nis from cis */
				nis = cis ;
				for( unsigned int i=0; i<par->nsp; ++i) { nis[i] /= (double)(i+1) ; }
			}

			/* normalized nis to 1 */ 
			ddummy1 = 0.0 ;
			for( unsigned int i=0; i<par->nsp; ++i) { ddummy1 += nis[i] ; }
			for( unsigned int i=0; i<par->nsp; ++i) { nis[i] /= ddummy1 ; }

			// assign nis to nis_arr, print nis_arr
			nis_arr = (double*) calloc( par->nsp, sizeof(double)) ;
			for( unsigned int i=0; i<par->nsp; ++i) { nis_arr[i] = nis[i] ; }

			if ( log_flag )
			{
				fprintf( logfile, "\tnis = ( ") ;
				for( unsigned int i=0; i<par->nsp-1; ++i){ fprintf( logfile, " %.3lf,", nis_arr[i]) ; }
				fprintf( logfile, " %.3lf)\n", nis_arr[par->nsp-1]) ;
				fprintf( logfile, "done\n\n") ; 
				fflush( logfile ) ;
			}


			/* now setup discrete PRNG for numbers 0...par->nsp-1 with weight according to nis */
			d_rand = gsl_ran_discrete_preproc( par->nsp, nis_arr) ;
		}

		/* initialize always n_rand PRNG, seed will be written from within n_rand */
		if ( log_flag ) { fprintf( logfile, "Initializing n_rand PRNG (n1,n2,n3,d1,d2,D,r3)\n") ; fflush (logfile) ; }
		if ( !init_n_rand_userdef )
		{
			/* random initialization using time */
			n_rand( -1, 0.0, 0.0) ;
		}
		else
		{
			/* use user supplied initialization number */
			n_rand( -1, 0.0, 0.0, init_n_rand) ;
		}

		/* also computes nav_pc, nav_ac, nav_r3, nav_MC and s1, ds */
		if ( log_flag ) { fprintf( logfile, "done\n\n") ; fprintf( logfile, "Compute polar and azimuthal angles for Powder Average\n") ; fflush (logfile) ; }
		compute_pol_azi_ang() ;

		if ( log_flag ) { fprintf( logfile, "done\n\n") ; fprintf( logfile, "reading element and isotope masses from file: %s\n", "atomweights.dat") ; fflush (logfile) ;}
		read_standard_atomic_weight_and_relative_atomic_mass() ;
		if ( Xn <= 0 )
		{
			if ( log_flag ) { fprintf( logfile, "done\n\n") ; fprintf( logfile, "read_neutron_scattering_lengths_and_cross_sections() : %s\n", "DeBe_NeutronNews.dat") ; fflush (logfile) ;}	
			read_neutron_scattering_lengths_and_cross_sections() ;
			if ( log_flag ) { fprintf( logfile, "done\n\n") ; fprintf( logfile, "Adding element and isotope masses to the neutron scattering database\n") ; fflush (logfile) ;}
			add_atom_mass_to_neutron_scattering_database() ;
		}
		if ( log_flag ) { fprintf( logfile, "done\n\n") ; fprintf( logfile, "read_cif_dictionary() : %s\n",cif_dic) ; fflush (logfile) ;}
		read_cif_dictionary() ;
		if ( log_flag ) { fprintf( logfile, "done\n\n") ; fprintf( logfile, "read_cif_file() : %s\n",cif_file) ; fflush (logfile) ;} 
		read_cif_file() ;
		if ( log_flag ) { fprintf( logfile, "done\n\n") ; fprintf( logfile, "Compute vectors {a_i}, {b_i}, {k1, k2, Gs}, {kg1, kg2, G0} and V_uc\n") ; fflush (logfile) ;}	
		compute_lattice_vectors() ;
		/* compute all scalar products between the three a_j's and all unit vectors eQ */
		if (log_flag ) { fprintf( logfile, "done\n\n") ; fprintf( logfile, "compute_eQvj()\n") ; fflush (logfile) ; }
		if ( time_flag > 0 ) { time(&last) ; }
		compute_eQvj() ; 
		if ( time_flag > 0 ) { time(&now) ; }
		if (log_flag )
		{
			if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
			else { fprintf( logfile, "done\n\n") ; }
			fflush (logfile) ;
		}

		/* atom will be used explicitely later (maybe even if there is a -Fr option) in :
		   - write_pcr_file()
		   - compute_electron_density()
		   - compute_scattering_length_density()
		   - compute_incoherent_cross_section_density()
		   - write_header_in_logfile()
		   - calculate_orientational_averaged_structure_factor()
		   - calculate_orientational_averaged_structure_amplitude()
		   - compute_f0table_and_add_f0ind_to_atom()
		   - add_neutron_coh_sc_length_and_inc_sc_crossec_to_atom()
		   - compute_coordinates_by_formula_interpreter()
		   - derive_atom_labels_and_types_from_cif()
		   - compute_structure_amplitude(...)
		*/
		if ( log_flag ) { fprintf( logfile, "compute_atoms_by_symmetry()\n") ; fflush (logfile) ; }
		compute_atoms_by_symmetry() ;
		if ( log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }

		/* basis transformation mode */
		/* read_pdb_file() must have been passed before */
		if ( bt_flag )
		{
			#if (ARMADILLO)
			bool func_found = false ;
			int func_num ;
			for ( unsigned int l=0; l<NOBT; ++l)
			{
				/* test bt on func_bt[l].name */
				if (strcmp (func_bt[l].name, bt) == 0)
				{
					func_found = true ;
					func_num = l ;
					if ( log_flag ) { fprintf( logfile, "basis transformation with function %s()\n", bt) ; fflush (logfile) ; }
					switch (func_num)
					{
						case 0: /* bt_twin function */
							bt_twin() ;
							break ;
						case 1: /* bt_double function */
							bt_double() ;
							break ;
					}
					if ( log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }
					break ;
				}
			}
			if ( func_found == false ) { XNDIFF_ERROR(4) ; }
			#else
			if ( log_flag ) { fprintf( logfile, "basis transformation with function %s() failed.\nArmadillo library must be installed and linked via -larmadillo.\nFurthermore g++ option -DARMADILLO must be used.\n", bt) ; fflush (logfile) ; exit(0) ; }
			#endif 
		}

		/* write cif-file for Jmol with current atom structure */
		if ( jmol_cif_flag ) 
		{
			if ( log_flag ) { fprintf( logfile, "write_cif_file_for_Jmol()\n") ; fflush (logfile) ; }
			write_cif_file_for_Jmol() ;
			if ( log_flag ) { fprintf( logfile, "done\n\n") ; }
		}

		if ( bt_flag )
		{
			/* always terminate */
			exit(0) ;
		}

		/* read the par-file (contrasts & thicknesses for stabilizer layer and dispersion medium, ...) */
		if ( log_flag ) { fprintf( logfile, "read_par_file() : %s\n",par_file) ; fflush (logfile) ; }
		read_par_file() ;
		if ( log_flag ) { fprintf( logfile, "done\n\n") ; }

		/* Similarly always compute first f0table and link the entries in atom to it via f0ind f0table will be used in :
		   - compute_structure_amplitude(...)
		   - calculate_orientational_averaged_structure_factor()
		   - calculate_orientational_averaged_structure_amplitude()
		   Same has to done with the neutron coherent scattering lengths and incoherent cross sections. 
		*/
		if ( Xn >= 0 )
		{
			if ( log_flag ) { fprintf( logfile, "compute_f0table_and_add_f0ind_to_atom()\n") ; fflush (logfile) ; }
			if ( time_flag > 0 ) { time(&last) ; }
			compute_f0table_and_add_f0ind_to_atom() ;
			if ( time_flag > 0 ) { time(&now) ; }
			if ( log_flag )
			{
				if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
				else { fprintf( logfile, "done\n\n") ; }
				fflush (logfile) ;
			}
		}
		if ( Xn <= 0 )
		{
			if ( log_flag ) { fprintf( logfile, "add_neutron_coh_sc_length_and_inc_sc_crossec_to_atom()\n") ; fflush (logfile) ; }
			add_neutron_coh_sc_length_and_inc_sc_crossec_to_atom();
			if ( log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }
		}

		/* test_XnEquiv() */
		if ( test_XnEquiv_flag ) { test_XnEquiv() ; }


		/* read or compute complex structure amplitudes FFF_X and FFF_n and keep them in memory.
		   If an error occurs compute FFF by setting Fr_flag = false */
		if ( Fr_flag == true )
		{
			if ( log_flag ) { fprintf( logfile, "read_structure_amplitude()\n") ; fflush (logfile) ; }
			if ( time_flag > 0 ) { time(&last) ; }
			if ( !read_structure_amplitude() )
			{ 
				Fr_flag = false ;
				if ( log_flag ) { fprintf( logfile, "failed\n\n") ; fflush (logfile) ; }
				if ( log_flag ) { fprintf( logfile, "Try to compute structure amplitude.\n\n") ; fflush (logfile) ; }
			}
			else
			{
				if ( time_flag > 0 ) { time(&now) ; }
				if ( log_flag )
				{
					if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
					else { fprintf( logfile, "done\n\n") ; }
					fflush (logfile) ;
				}
			}
		}

		/* if -Fr option was not succesful compute FFF */
		if ( Fr_flag == false ) /* don't use else ! */
		{
			if ( log_flag ) { fprintf( logfile, "compute_structure_amplitude()\n") ; fflush (logfile) ; }
			if ( time_flag > 0 ) { time(&last) ; }
			compute_structure_amplitude() ;
			if ( time_flag > 0 ) { time(&now) ; }
			if ( log_flag )
			{
				if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
				else { fprintf( logfile, "done\n\n") ; }
				fflush (logfile) ;
			}
		}

		/* create pcr file for FullProf Suite with crystallographic information */
		if ( pcr_flag )
		{
			if ( log_flag ) { fprintf( logfile, "write_pcr_file() : %s\n", pcr_file) ; fflush (logfile) ; }	
			write_pcr_file() ; 
			if ( log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }
		}

		if ( Xn >= 0 )
		{
			if ( log_flag ) { fprintf( logfile, "Compute electron density for crystal\n") ; fflush (logfile) ;}
			compute_electron_density() ;
			if ( log_flag ) { fprintf( logfile, "done\n\n") ; }
		}
		if ( Xn <= 0 )
		{
			if ( log_flag ) { fprintf( logfile, "Compute incoherent cross section density for crystal\n") ; fflush (logfile) ;}
			compute_incoherent_cross_section_density() ;
			if ( log_flag ) { fprintf( logfile, "done\n\n") ; fprintf( logfile, "Compute scattering length density for crystal\n") ; fflush (logfile) ;}
			compute_scattering_length_density() ;
			if ( log_flag ) { fprintf( logfile, "done\n\n") ; }
		}
		if ( log_flag ) { fprintf( logfile, "Summarize important parameters\n") ; fflush (logfile) ; }	
		if ( log_flag ) { write_header_in_logfile("\t") ; }


		/* define sum of stacks and single particles -1 */
		nspsp = par->nms + par->nsp - 1 ;


		/* write structure amplitude */
		if ( Fw_flag && !FNoMem_flag )
		{
			if ( log_flag ) { fprintf( logfile, "write_structure_amplitude()\n") ; fflush (logfile) ; }
			if ( time_flag > 0 ) { time(&last) ; }
			if ( write_structure_amplitude() )
			{
				if ( time_flag > 0 ) { time(&now) ; }
				if ( log_flag )
				{
					if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
					else { fprintf( logfile, "done\n\n") ; }
					fflush (logfile) ;
				}
			}
			else
			{
				if ( time_flag > 0 ) { time(&now) ; }
				if ( log_flag )
				{
					if ( time_flag > 0 ) { fprintf( logfile, "failed (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
					else { fprintf( logfile, "failed\n\n") ; }
					fflush (logfile) ;
				}
			}
		}
		/* option to perform orientational averages <F_X>, <F_n>, <|F_X|^2> and <|F_n|^2> */
		if ( orav_flag )
		{
			if ( !FNoMem_flag )
			{
				if ( Xn >= 0 )
				{
					if ( log_flag ) { fprintf( logfile, "compute_orientational_averaged_structure_factor( &FFF_X[0], 'X') for X-ray\n") ; fflush (logfile) ; }
					if ( time_flag > 0 ) { time(&last) ; }
					compute_orientational_averaged_structure_factor( &FFF_X[0], 'X') ;
					if ( time_flag > 0 ) { time(&now) ; }
					if ( log_flag )
					{
						if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
						else { fprintf( logfile, "done\n\n") ; }
						fflush (logfile) ;
					}
					if ( log_flag ) { fprintf( logfile, "compute_orientational_averaged_structure_amplitude( &FFF_X[0], 'X') for X-ray\n") ; fflush (logfile) ; }
					if ( time_flag > 0 ) { time(&last) ; }
					compute_orientational_averaged_structure_amplitude( &FFF_X[0], 'X') ;
					if ( time_flag > 0 ) { time(&now) ; }
					if ( log_flag )
					{
						if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
						else { fprintf( logfile, "done\n\n") ; }
						fflush (logfile) ;
					}
				}

				if ( Xn <= 0 )
				{
					if ( log_flag ) { fprintf( logfile, "compute_orientational_averaged_structure_factor( &FFF_n[0], 'n') for neutrons\n") ; fflush (logfile) ;} 
					if ( time_flag > 0 ) { time(&last) ; }
					compute_orientational_averaged_structure_factor( &FFF_n[0], 'n') ;
					if ( time_flag > 0 ) { time(&now) ; }
					if ( log_flag )
					{
						if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
						else { fprintf( logfile, "done\n\n") ; }
						fflush (logfile) ;
					}
					if ( log_flag ) { fprintf( logfile, "compute_orientational_averaged_structure_amplitude( &FFF_n[0], 'n') for neutrons\n") ; fflush (logfile) ;} 
					if ( time_flag > 0 ) { time(&last) ; }
					compute_orientational_averaged_structure_amplitude( &FFF_n[0], 'n') ;
					if ( time_flag > 0 ) { time(&now) ; }
					if ( log_flag )
					{
						if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
						else { fprintf( logfile, "done\n\n") ; }
						fflush (logfile) ;
					}
				}
			}

			/* analytical calculations */ 
			if ( log_flag ) { fprintf( logfile, "calculate_orientational_averaged_structure_amplitude()\n") ; fflush (logfile) ; }
			if ( time_flag > 0 ) { time(&last) ; }
			calculate_orientational_averaged_structure_amplitude() ;
			if ( time_flag > 0 ) { time(&now) ; }
			if ( log_flag )
			{
				if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
				else { fprintf( logfile, "done\n\n") ; }
				fflush (logfile) ;
			}

			if ( log_flag ) { fprintf( logfile, "calculate_orientational_averaged_structure_factor()\n") ; fflush (logfile) ; }
			// if ( log_flag ) { fprintf( logfile, "calculate_orientational_averaged_structure_factor_in_solvent()\n") ; fflush (logfile) ; }
			if ( time_flag > 0 ) { time(&last) ; }
			calculate_orientational_averaged_structure_factor() ;
			// for test purposes calculate_orientational_averaged_structure_factor_in_solvent() ;
			if ( time_flag > 0 ) { time(&now) ; }
			if ( log_flag )
			{
				if ( time_flag > 0 ) { fprintf( logfile, "done (%ld seconds)\n\n",(long int)(difftime( now, last)) ) ; }
				else { fprintf( logfile, "done\n\n") ; }
				fflush (logfile) ;
			}
		}


		/* allocate memory for variables used only (locally) in stackcpp */
		if ( log_flag ) { fprintf( logfile, "Allocating memory\n") ; fflush (logfile) ; }

		n1=(unsigned int *) calloc( par->nms, sizeof(unsigned int));
		n2=(unsigned int *) calloc( par->nms, sizeof(unsigned int));
		n3=(unsigned int *) calloc( par->nms, sizeof(unsigned int));
		r3=(double *) calloc( par->nms, sizeof(double));
		ind_r3=(unsigned int *) calloc( par->nms, sizeof(unsigned int));
		v1=(double **) calloc( par->nms, sizeof(double *));
		v2=(double **) calloc( par->nms, sizeof(double *));
		v3=(double **) calloc( par->nms, sizeof(double *));
		w1=(double **) calloc( par->nms, sizeof(double *));
		w2=(double **) calloc( par->nms, sizeof(double *));
		wt1=(double **) calloc( par->nms, sizeof(double *));
		wt2=(double **) calloc( par->nms, sizeof(double *));

		for ( unsigned int i=0; i<par->nms; ++i)
		{
			v1[i]=(double *) calloc( 3, sizeof(double));
			v2[i]=(double *) calloc( 3, sizeof(double));
			v3[i]=(double *) calloc( 3, sizeof(double));
			w1[i]=(double *) calloc( 3, sizeof(double));
			w2[i]=(double *) calloc( 3, sizeof(double));
			wt1[i]=(double *) calloc( 3, sizeof(double));
			wt2[i]=(double *) calloc( 3, sizeof(double));
		}
		kf1=(double *) calloc( par->nms, sizeof(double));
		kf2=(double *) calloc( par->nms, sizeof(double));
		okf1=(double *) calloc( par->nms, sizeof(double));
		okf2=(double *) calloc( par->nms, sizeof(double));

		kf3=(double *) calloc( par->nsp, sizeof(double));
		okf3=(double *) calloc( par->nsp, sizeof(double));

		w3=(double ***) calloc( par->nsp, sizeof(double **));
		wt3=(double ***) calloc( par->nsp, sizeof(double **));

		if ( Xn >= 0 )
		{
			drerho_dm_osl = (double *) calloc( rho_multi.n, sizeof(double)) ;
			drerho_osl_isl = (double *) calloc( rho_multi.n, sizeof(double)) ;
			drerho_isl = (double *) calloc( rho_multi.n, sizeof(double)) ;
		}
		if ( Xn <= 0 )
		{
			dsld_dm_osl = (double *) calloc( sld_multi.n, sizeof(double)) ;
			dsld_osl_isl = (double *) calloc( sld_multi.n, sizeof(double)) ;
			dsld_isl = (double *) calloc( sld_multi.n, sizeof(double)) ;
		}

		Vtot_cry = (double *) calloc( nspsp, sizeof(double)) ;
		for ( unsigned int i=0; i<nspsp; ++i) { Vtot_cry[i] = 0.0 ; }
		Vtot_isl = (double *) calloc( nspsp, sizeof(double)) ;
		for ( unsigned int i=0; i<nspsp; ++i) { Vtot_isl[i] = 0.0 ; }
		Vtot_osl = (double *) calloc( nspsp, sizeof(double)) ;
		for ( unsigned int i=0; i<nspsp; ++i) { Vtot_osl[i] = 0.0 ; }
		Vtot_irr = (double *) calloc( nspsp, sizeof(double)) ;
		for ( unsigned int i=0; i<nspsp; ++i) { Vtot_irr[i] = 0.0 ; }
		Vtot_dm = (double *) calloc( nspsp, sizeof(double)) ;
		for ( unsigned int i=0; i<nspsp; ++i) { Vtot_dm[i] = 0.0 ; }

		V = (double **) calloc( par->nsp, sizeof(double *)) ;
		V_isl = (double **) calloc( par->nsp, sizeof(double *)) ;
		V_osl = (double **) calloc( par->nsp, sizeof(double *)) ;

		for ( unsigned int i=0; i<par->nsp; ++i)
		{
			w3[i]=(double **) calloc( par->nms, sizeof(double *));
			wt3[i]=(double **) calloc( par->nms, sizeof(double *));
			for ( unsigned int j=0; j<par->nms; ++j)
			{
				w3[i][j]=(double *) calloc( 3, sizeof(double)) ;
				wt3[i][j]=(double *) calloc( 3, sizeof(double)) ;
			}

			V[i] = (double *) calloc( par->nms, sizeof(double)) ;
			V_isl[i] = (double *) calloc( par->nms, sizeof(double)) ;
			V_osl[i] = (double *) calloc( par->nms, sizeof(double)) ;
		}

		VN = (double **) calloc( nspsp, sizeof(double *)) ;
		VN_isl = (double **) calloc( nspsp, sizeof(double *)) ;
		VN_osl = (double **) calloc( nspsp, sizeof(double *)) ;
		VN_dm = (double **) calloc( nspsp, sizeof(double *)) ;
		VN_irr = (double **) calloc( nspsp, sizeof(double *)) ;
		for ( unsigned int i=0; i<par->nsp; ++i)
		{
			VN[i] = (double *) calloc( par->nms, sizeof(double)) ;
			VN_isl[i] = (double *) calloc( par->nms, sizeof(double)) ;
			VN_osl[i] = (double *) calloc( par->nms, sizeof(double)) ;
			VN_dm[i] = (double *) calloc( par->nms, sizeof(double)) ;
			VN_irr[i] = (double *) calloc( par->nms, sizeof(double)) ;
		}
		for ( unsigned int i=par->nsp; i<nspsp; ++i)
		{
			VN[i] = (double *) calloc( 1, sizeof(double)) ;
			VN_isl[i] = (double *) calloc( 1, sizeof(double)) ;
			VN_osl[i] = (double *) calloc( 1, sizeof(double)) ;
			VN_dm[i] = (double *) calloc( 1, sizeof(double)) ;
			VN_irr[i] = (double *) calloc( 1, sizeof(double)) ;
		}

		R=(double **) calloc( par->nms, sizeof(double *)) ;
		for ( unsigned int i=0; i<par->nms; ++i) { R[i]=(double *) calloc(3, sizeof(double)) ; }

		d1=(double *) calloc( par->nms, sizeof(double)) ;
		d2=(double *) calloc( par->nms, sizeof(double)) ;
		D=(double *) calloc( par->nms, sizeof(double)) ;

		MULTIPLICITY = (unsigned int *) calloc( nspsp, sizeof(unsigned int)) ;

		/* S is a double matrix of size ( nspsp x rho/sld_multi.n x par->np ) */
		if ( Xn >= 0 )
		{
			S_X = (double ***) calloc( nspsp, sizeof(double **)) ;
			dS_X = (double ****) calloc( nspsp, sizeof(double ***)) ;
		}
		if ( Xn <= 0 )
		{
			S_n = (double ***) calloc( nspsp, sizeof(double **)) ;
			dS_n = (double ****) calloc( nspsp, sizeof(double ***)) ;
		}
		for ( unsigned int i=0; i<nspsp; ++i)
		{
			/* S_X/n, initialize matrices to zero */
			if ( Xn >= 0 )
			{
				S_X[i] = (double **) calloc( rho_multi.n, sizeof(double *)) ;
				for ( unsigned int q=0; q<rho_multi.n; ++q)
				{
					S_X[i][q] = (double *) calloc( par->np, sizeof(double)) ;
					for ( unsigned int j=0; j<par->np; ++j) { S_X[i][q][j] = 0.0 ; }
				}
			}
			if ( Xn <= 0 )
			{
				S_n[i] = (double **) calloc( sld_multi.n, sizeof(double *)) ;
				for ( unsigned int q=0; q<sld_multi.n; ++q)
				{
					S_n[i][q] = (double *) calloc( par->np, sizeof(double)) ;
					for ( unsigned int j=0; j<par->np; ++j) { S_n[i][q][j] = 0.0 ; }
				}
			}

			/*  dS_X/n, is set to zero in the beginning of each repeat step */
			if ( i < par->nsp )
			{
				/* dS single particle */
				if ( Xn >= 0 )
				{
					dS_X[i] = (double ***) calloc( par->nms, sizeof(double **)) ;
					for ( unsigned int j=0; j<par->nms; ++j)
					{
						dS_X[i][j] = (double **) calloc( rho_multi.n, sizeof(double *)) ;
						for ( unsigned int q=0; q<rho_multi.n; ++q)
							dS_X[i][j][q] = (double *) calloc( par->np, sizeof(double)) ;
					}
				}
				if ( Xn <= 0 )
				{
					dS_n[i] = (double ***) calloc( par->nms, sizeof(double **)) ;
					for ( unsigned int j=0; j<par->nms; ++j)
					{
						dS_n[i][j] = (double **) calloc( sld_multi.n, sizeof(double *)) ;
						for ( unsigned int q=0; q<sld_multi.n; ++q)
							dS_n[i][j][q] = (double *) calloc( par->np, sizeof(double)) ;

					}
				}
			}
			else
			{
				/* dS stacks */
				if ( Xn >= 0 )
				{
					dS_X[i] = (double ***) calloc( 1, sizeof(double **)) ;
					dS_X[i][0] = (double **) calloc( rho_multi.n, sizeof(double *)) ;
					for ( unsigned int q=0; q<rho_multi.n; ++q)
						dS_X[i][0][q] = (double *) calloc( par->np, sizeof(double)) ;

				}
				if ( Xn <= 0 )
				{
					dS_n[i] = (double ***) calloc( 1, sizeof(double **)) ;
					dS_n[i][0] = (double **) calloc( sld_multi.n, sizeof(double *)) ;
					for ( unsigned int q=0; q<sld_multi.n; ++q)
						dS_n[i][0][q] = (double *) calloc( par->np, sizeof(double)) ;

				}
			}
		}

		/* Allocate matrices B_n ( nspsp x sld_multi.n ) */
		if ( Xn <= 0 )
		{
			dB_n = (double ***) calloc( nspsp, sizeof(double **)) ;
			B_n = (double **) calloc( nspsp, sizeof(double *)) ;

			/* Initialize vector B_n to zero */
			for ( unsigned int j=0; j<nspsp; ++j)
			{
				B_n[j] = (double *) calloc( sld_multi.n, sizeof(double)) ;
				for ( unsigned int q=0; q<sld_multi.n; ++q) { B_n[j][q] = 0.0 ; }

				if (j<par->nsp)
				{
					dB_n[j] = (double **) calloc( par->nms, sizeof(double *)) ;
					for ( unsigned int k=0; k<par->nms; ++k)
						dB_n[j][k] = (double *) calloc( sld_multi.n, sizeof(double)) ;
				}
				else
				{
					dB_n[j] = (double **) calloc( 1, sizeof(double *)) ;
					dB_n[j][0] = (double *) calloc( sld_multi.n, sizeof(double)) ;
				}
			}
		}

		/* Yc(i)_X/n are double matrices of size ( nspsp x 10(4) x par->np ) */
		if ( Xn >= 0 ) { Yc_X = (double ***) calloc( nspsp, sizeof(double **)) ; }
		if ( Xn <= 0 ) 
		{
			Yc_n = (double ***) calloc( nspsp, sizeof(double **)) ;
			Yi_n = (double **) calloc( nspsp, sizeof(double **)) ;
		}

		for ( unsigned int i=0; i<nspsp; ++i)
		{
			if ( Xn >= 0 )
			{
				Yc_X[i] = (double **) calloc( 10, sizeof(double *)) ;
				for (int j=0; j<10; ++j) 
				{
					Yc_X[i][j] = (double *) calloc( par->np, sizeof(double)) ;
					for ( unsigned int k=0; k<par->np; ++k) { Yc_X[i][j][k] = 0.0 ; }
				}
			}
			if ( Xn <= 0 )
			{
				Yc_n[i] = (double **) calloc( 10, sizeof(double *)) ;
				for (int j=0; j<10; ++j)
				{
					Yc_n[i][j] = (double *) calloc( par->np, sizeof(double)) ;
					for ( unsigned int k=0; k<par->np; ++k) { Yc_n[i][j][k] = 0.0 ; }
				}

				Yi_n[i] = (double *) calloc( 4, sizeof(double)) ;
			}
		}

		if (log_flag) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }

		if ( mem_flag > 0 )
		{
			if (log_flag) { fprintf( logfile, "Determine used memory: get_current_memory()\n") ; fflush (logfile) ; }
			get_current_memory() ;
			if (log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }
		}

		/* prepare drerho_... and dsld_... , apply scale factor such that each of them is in [1/nm^2] */
		if ( Xn >= 0 )
		{
			for ( unsigned int q=0; q<rho_multi.n; ++q)
			{
				drerho_dm_osl[q] = sc_rho * ( rho_dm[q] - rho_osl[q] ) ;
				drerho_osl_isl[q] = sc_rho * ( rho_osl[q] - rho_isl[q] ) ;
				drerho_isl[q] = sc_rho * rho_isl[q] ;
			}
		}
		if ( Xn <= 0 )
		{
			for ( unsigned int q=0; q<sld_multi.n; ++q)
			{
				dsld_dm_osl[q] = sc_sld * ( sld_dm[q] - sld_osl[q] ) ;
				dsld_osl_isl[q] = sc_sld * ( sld_osl[q] - sld_isl[q] ) ;
				dsld_isl[q] = sc_sld * sld_isl[q] ;
			}
		}


		/*******************************/
		/* BEGIN OF ENSEMBLE AVERAGING */
		/*******************************/

		if ( log_flag )
		{
			fprintf( logfile, "Ensemble and Powder Average\n\n") ;
			if ( par->av_mode == 0 )
			{
				fprintf( logfile, "%sPowder Average using Polar-Azimuthal grid with nav_pc x nav_ac = %d x %d:\n", "\t", nav_pc, nav_ac) ;
				fprintf( logfile, "%sPolar: av_p1=%10.5G, av_p2=%10.5G, av_pd=%10.5G\n", "\t", par->av_pol[0],  par->av_pol[1],  par->av_pol[2]) ;
				fprintf( logfile, "%sAzimuthal: av_a1=%10.5G, av_a2=%10.5G, av_ad=%10.5G\n\n", "\t", par->av_azi[0],  par->av_azi[1],  par->av_azi[2]) ;
			}
			else if ( par->av_mode == 1 )
			{
				fprintf( logfile, "%sPowder Average using Polar grid and random r3-rotation with nav_pc x nav_ac = %d x 1:\n", "\t", nav_pc) ;
				fprintf( logfile, "%sPolar: av_p1=%10.5G, av_p2=%10.5G, av_pd=%10.5G\n", "\t", par->av_pol[0],  par->av_pol[1],  par->av_pol[2]) ;
				fprintf( logfile, "%sr3-rotation: av_r3[0]=%10.5G, av_r3[1]=%10.5G, av_r3d=%10.5G\n\n", "\t", par->av_r3[0], par->av_r3[1], av_r3d) ;
			}
			else if ( ( par->av_mode == 2 ) || ( par->av_mode == 3 ) )
			{
				fprintf( logfile, "%sPowder Average using nav_MC = %d MC points:\n", "\t", nav_MC) ;
				fprintf( logfile, "%sPolar: av_p1=%10.5G, av_p2=%10.5G\n", "\t", par->av_pol[0],  par->av_pol[1]) ;
				fprintf( logfile, "%sAzimuthal: av_a1=%10.5G, av_a2=%10.5G\n\n", "\t", par->av_azi[0],  par->av_azi[1]) ;
			}
			#pragma omp parallel if ( openmp_flag )
			{
				if ( omp_in_parallel() )
				{
					#pragma omp single
					if ( log_flag ) { fprintf( logfile ,"%sParallelized Run: #-CPU=%d, #-Threads=%d\n", "\t", omp_get_num_procs(), omp_get_num_threads() ) ; fflush( logfile ) ; }
				}
				else
				{
					if ( log_flag ) { fprintf( logfile ,"%sSerialized Run\n", "\t") ; fflush( logfile ) ; }
				}
			}
			fprintf( logfile, "\n") ;
			fflush (logfile) ;
		}

		/* Progress-bar for averaging */
		if ( time_flag > 1 ) { fprintf( stdout ,"\n" ) ; }
		double pr ; /* progress indicator */
		char rtstr[30] ; /* time string filled by remtime */

		time_t rep_start ; /* starting point of averaging */
		if ( time_flag > 1 )
		{
			time(&rep_start) ; /* for overall progress of repeats */
			time(&last) ; /* now and last for a single repeat step c_rep */
		}

		/* if the XNDiff PRN shall be saved, first write a header to the file */
		if ( distr_flag )
		{
			sprintf( oname, "%s%s.distr", ofo, mfname) ;
			if ( ( outf = fopen( oname, "w")) == NULL) { XNDIFF_ERROR(18) ; }

			fprintf( outf, "# Random numbers in XNDiff\n" ) ;
			fprintf( outf, "# Use distribution types: td1=%d, td2=%d\n", par->td1, par->td2) ;

			if ( ( par->av_mode == 0 ) || ( par->av_mode == 2 ) || ( par->av_mode == 3 ) )
			{
				fprintf( outf, "# %6s %6s %6s %15s %15s %15s\n", var_name(n1), var_name(n2), var_name(n3), var_name(d1), var_name(d2), var_name(D) ) ;
			}
			else if ( par->av_mode == 1 ) { fprintf( outf, "# %6s %6s %6s %15s %15s %15s %15s\n", var_name(n1), var_name(n2), var_name(n3), var_name(d1), var_name(d2), var_name(D), var_name(r3) ) ; }

			fclose(outf) ;
		}

		/* FIRST loop ensemble averaging */
		sum_weight_factor = 0.0 ;
		for ( unsigned int c_rep=0; c_rep<par->nr; ++c_rep)
		{
			if ( time_flag > 1 )
			{
				pr = 100.0*((double) c_rep)/((double) par->nr) ;
				time(&now) ;
				remtime( (long int)(difftime( now, rep_start)), pr, &rtstr[0]) ;
			}

			/* write event of a new stack to logfile */
			if ( log_flag ) { fprintf( logfile, "%s*** -> New Stack: c_rep: %d\n", "\t", c_rep) ; fflush(logfile) ; }


			/* protocol current memory */
			if ( mem_flag > 1 ) { get_current_memory() ; }

			/* reset sum over weight factors for a single particle to zero */
			single_sum_weight_factor = 0.0 ;

			/* initialize dS_X/n, dB_n always to zero for every c_rep round */
			if ( Xn >= 0 )
			{
// 				for ( long int i=0; i<par->np; ++i)
// 				{
// 					for ( int j=0; j<nspsp; ++j)
// 					{
// 						if ( j < par->nsp )
// 						{
// 							for ( int k=0; k<par->nms; ++k)
// 								for ( int q=0; q<rho_multi.n; ++q)
// 									dS_X[i][j][k][q] = 0.0 ;
// 						}
// 						else
// 						{
// 							for (int q=0; q<rho_multi.n; ++q)
// 								dS_X[i][j][0][q] = 0.0 ;
// 						}
// 					}
// 				}

				for ( unsigned int j=0; j<nspsp; ++j)
				{
					if ( j < par->nsp )
					{
						for ( unsigned int k=0; k<par->nms; ++k)
						{
							for ( unsigned int q=0; q<rho_multi.n; ++q)
							{
								for ( unsigned int i=0; i<par->np; ++i) { dS_X[j][k][q][i] = 0.0 ; }
							}
						}
					}
					else
					{
						for ( unsigned int q=0; q<rho_multi.n; ++q)
						{
							for ( unsigned int i=0; i<par->np; ++i) { dS_X[j][0][q][i] = 0.0 ; }
						}
					}
				}
			}
			if ( Xn <= 0 )
			{
// 				for ( long int i=0; i<par->np; ++i)
// 				{
// 					for ( int j=0; j<nspsp; ++j)
// 					{
// 						if ( j < par->nsp )
// 						{
// 							for ( int k=0; k<par->nms; ++k)
// 								for ( int q=0; q<sld_multi.n; ++q)
// 									dS_n[i][j][k][q] = 0.0 ;
// 						}
// 						else
// 						{
// 							for ( int q=0; q<sld_multi.n; ++q)
// 								dS_n[i][j][0][q] = 0.0 ;
// 						}
// 					}
// 				}

				for ( unsigned int j=0; j<nspsp; ++j)
				{
					if ( j < par->nsp )
					{
						for ( unsigned int k=0; k<par->nms; ++k)
						{
							for ( unsigned int q=0; q<sld_multi.n; ++q)
							{
								for ( unsigned int i=0; i<par->np; ++i) { dS_n[j][k][q][i] = 0.0 ; }
							}
						}
					}
					else
					{
						for ( unsigned int q=0; q<sld_multi.n; ++q)
						{
							for ( unsigned int i=0; i<par->np; ++i) { dS_n[j][0][q][i] = 0.0 ; }
						}
					}
				}

// 				for ( int j=0; j<nspsp; ++j)
// 				{
// 					if ( j < par->nsp )
// 					{
// 						for ( int k=0; k<par->nms; ++k)
// 							for ( int q=0; q<sld_multi.n; ++q) { dB_n[j][k][q] = 0.0 ; }
// 					}
// 					else
// 					{
// 						for ( int q=0; q<sld_multi.n; ++q) { dB_n[j][0][q] = 0.0 ; }
// 					}
// 				}
			}

			/* SECOND.1 loop for the nms crystalline particles in one stack */ 
			for ( unsigned int k=0; k<par->nms; ++k)
			{
				/****************************************/
				/* RANDOM NUMBERS GENERATION FOR STACKS */
				/****************************************/

				/* create random distributed crystalline particle of size n_i*a_i in the i-th unit cell direction i=1,2,3 */
				/* i.e. n_i is the number of unit cells along a_i direction i=1,2,3 */
				/* ni1,ni2 are the mean values and variances of the number of unit cells along a_i direction */

				/* setup n1 and n2 with par->td1 */
				do { ddummy1 = floor( n_rand( par->td1, par->n11, par->n12) / cell_par[0] + 0.5) ; }
				while ( ddummy1 < 1 ) ;
				n1[k] = (unsigned int) ddummy1 ;

				do { ddummy1 = floor( n_rand( par->td1, par->n21, par->n22) / cell_par[1] + 0.5) ; }
				while ( ddummy1 < 1 ) ;
				n2[k] = (unsigned int) ddummy1 ;

                                /* 
				n1[k] = (unsigned int) floor (n_rand (par->td1, par->n11, par->n12) / cell_par[0] + 0.5) ;
				if (n1[k] < 1) { n1[k] = 1 ; }

				n2[k] = (unsigned int) floor (n_rand (par->td1, par->n21, par->n22) / cell_par[1] + 0.5) ;
				if (n2[k] < 1) { n2[k] = 1 ; }
                                */

				/* setup n3 with par->td2 */
				if ( par->td2 == 2 )
				{
					/* discrete distribution according to nis */
					n3[k] = (unsigned int) gsl_ran_discrete( GSL_RNG, d_rand) + 1 ;
				}
				else
				{
					/* (log)normal distribution according to par->td2 and par->n31, par->n32 */
					do
					{
						ddummy1 = floor( n_rand( par->td2, par->n31, par->n32) / cell_par[2] + 0.5) ;
					}
					while ( ( ddummy1 < 1 ) || ( ddummy1 > par->nsp ) ) ;
					n3[k] = (unsigned int) ddummy1 ;
				}

				if (k == 0)
				{
					/* set the first crystalline particle in a stack into the origin */
					D[k] = d1[k] = d2[k] = R[0][0] = R[0][1] = R[0][2] = 0.0 ;
				}
				else
				{
					/* D and d1,d2 are the displacements between the crystalline particles in nm */
					/* along Gs and k1,k2, where all vectors are perpendicular on each other */
					/* d1i,d2i are the mean displacement and variance of crystals in stack along k1,k2 */

					/* setup DMIN without stabilizer layers */
					DMIN  = n3[k-1] * abs( sp( Gs, a3) ) / betrag (Gs) ;

					if ( par->stackmode == 0 )
					{
						/* 
						   for normal stacks, use distribution function, apply no collision check, only that D is positive
						*/
						while ((D[k] = n_rand (par->td1, par->D1, par->D2)) < 0.0 ) ;
					}
					else if ( par->stackmode == 1 )
					{
						/* 
						   for normal stacks, use distribution function,
						   check that D[k] > thickness of particle k-1 (i.e. no collisions),
						   user should provide reasonable parameters for nsp=1,2,3,... , distribution for n3
						   thickness_isl/osl and the distribution for D, 
						   i.e. <D> ~ <n3> * <Gs,a3> + 2 *(thickness_isl+thickness_osl)
						   and var(D) must be broad enough to include also D ~ max(nsp) * <Gs,a3> + 2 *(thickness_isl+thickness_osl)
						*/
						DMIN += 2.0 * ( thickness_isl + thickness_osl ) ;
						while ((D[k] = n_rand (par->td1, par->D1, par->D2)) < DMIN) ;
					}
					else if ( par->stackmode == 2 )
					{
						/* 
						   for DNA stacks add n3 times the projection of a3 on Gs to D, that is the n3-thickness,
						   leave out the soft stabilizer layer, this is however questionable 
						*/
						/*
						   old
						   while ( (D[k] = n_rand (par->td, par->D1, par->D2)) < 0) ;
						   D[k] += n3[k-1] * abs( sp( Gs, a3) ) / betrag (Gs) ;
						   new
						   D[k] = n3[k-1] * abs( sp( Gs, a3) ) / betrag (Gs) ; ~ DMIN
						*/
						D[k] = DMIN ;

					}
					else if ( par->stackmode == 3 )
					{
						/* for DNA/DMPC stacks add n3 times the projection of a3 on Gs to D, that is the n3-thickness */
						/* include also both soft stabilizer layer */
						D[k] = DMIN + 2.0 * ( thickness_isl + thickness_osl ) ;
					}
					else if ( par->stackmode == 4 )
					{
						/* if DMIN + both soft stabilizer layers >  D[k] -> D[k] = DMIN */
						/* reasonable if difference btw max(DMIN) ~ <D[k]> but not max(DMIN) >> <D[k]> */
						DMIN += 2.0 * ( thickness_isl + thickness_osl ) ;
						while ((D[k] = n_rand (par->td1, par->D1, par->D2)) < 0.0 ) ;
						if ( DMIN > D[k] ) { D[k] = DMIN ; } 
					}

					d1[k] = n_rand (par->td1, par->d11, par->d12) ;
					d2[k] = n_rand (par->td1, par->d21, par->d22) ;

					/* R is the vector to origin of the crystalline particles, see figure in paper */
					/* R[k]=Gs*D[k]+k1*d1[k]+k2*d2[k] */
					/* note k1,k2 unit vectors, Gs in general not, therefore distance along Gs is |Gs|*D[k] */
		
					R[k][0] = Gs[0] * D[k] + d1[k] * k1[0] + d2[k] * k2[0] ;
					R[k][1] = Gs[1] * D[k] + d1[k] * k1[1] + d2[k] * k2[1] ;
					R[k][2] = Gs[2] * D[k] + d1[k] * k1[2] + d2[k] * k2[2] ;

					/* skip other overlap checks and add immediately the previous R[k-1] */
					vadd (R[k], R[k-1]) ;
				}
	
				/********************************/
				/* COLLISION CONTROL FOR STACKS */ 
				/********************************/

				/* Old collision control algorithm has been been removed since:
				- only r3-rotation around Gs (av_mode=1) applies for single particles
				- R=D*Gs+d1*k1+d2*k2 shifts the single particles sufficiently large away from each 
				- Later on it would be good to use a more simpler collision control which controls if 
				D > Thickness of the current particle + epsilon, nevertheless the current 
				distribution of R is quite narrow (as was also observed experimentally by EM)
				*/


				/***************************************************************/
				/* LATTICE VECTORS AND UNIT CELL VOLUMES FOR CURRENT PARTICLES */
				/***************************************************************/

				/* Calculate lattice vectors v_l from a_l with r3=0 (av_mode=0.2.3) and randomly chosen angles r3 (av_mode=1)
				r1,r2=0 holds for all modes i.e. all particles remain parallel to each other */
				if ( ( par->av_mode == 0 ) || ( par->av_mode == 2 ) || ( par->av_mode == 3 ) )
				{
					/* if r1=r2=r3=0 then v_l=a_l */
					v1[k][0]=a1[0];
					v1[k][1]=a1[1];
					v1[k][2]=a1[2];
	
					v2[k][0]=a2[0];
					v2[k][1]=a2[1];
					v2[k][2]=a2[2];
	
					v3[k][0]=a3[0];
					v3[k][1]=a3[1];
					v3[k][2]=a3[2];
				}
				else if ( par->av_mode == 1 )
				{
					/* compute random (azimuthal) rotation around the Gs-axis in ° */
// 					r3[k] = n_rand (par->td, par->r31, par->r32) ; 
// 
// 					if (r3[k]>par->av_r3[1] || r3[k]<par->av_r3[0])
// 					{	/* restrict to range between av_r3[0] and av_r3[1] */
// 						if (r3[k]>par->av_r3[1])
// 						{
// 							ind_r3[k]=nav_r3-1;
// 						}
// 						else
// 						{
// 							ind_r3[k]=0;
// 						}
// 					}
// 					else
// 					{
// 						/* round to r3-grid and compute index ind_r3[k] */
// 						ind_r3[k]=(int)((r3[k]-par->av_r3[0])/av_r3d);
// 					}
					do 
					{
						r3[k] = n_rand( par->td1, par->r31, par->r32) ; 
					}
					while ( ( r3[k]<par->av_r3[0] ) || ( r3[k]>par->av_r3[1] ) ) ;

					/* round to r3-grid and compute index ind_r3[k] */
					ind_r3[k]=(unsigned int)((r3[k]-par->av_r3[0])/av_r3d);
					r3[k]=par->av_r3[0]+(double)(ind_r3[k])*av_r3d;

					/* Calculate new unit vectors after rotation */

					/* Rotate crystal k on Gs -> e1,e2,e3 */
					/* e_l=<a_l,G_s>G_s+(a_l-<a_l,G_s>G_s)*cos(r3[k])+(a_l x G_s)sin(r3[k]) for l=1,2,3 */
					/* a_l are the lattice vectors of length a_l */
					/* G_s unit vector */
					/* e_l are rotated vectors of length a_l */
					/* due to r1,r2=0 v_l=e_l */
					ddummy1 = Gs[0] ;
					ddummy2 = Gs[1] ;
					ddummy3 = Gs[2] ;
					ddummy4 = cos (DR*r3[k]) ;
					ddummy5 = sin (DR*r3[k]) ;
					v1[k][0] = a1[0] * (ddummy4+SQUARE(ddummy1)*(1.0-ddummy4)) +
							a1[1] * (ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a1[2] * (-ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) ;
					v1[k][1] = a1[0] * (-ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a1[1] * (ddummy4+SQUARE(ddummy2)*(1.0-ddummy4)) +
							a1[2] * (ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) ;
					v1[k][2] = a1[0] * (ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) +
							a1[1] * (-ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) +
							a1[2] * (ddummy4+SQUARE(ddummy3)*(1.0-ddummy4)) ;
					v2[k][0] = a2[0] * (ddummy4+SQUARE(ddummy1)*(1.0-ddummy4)) +
							a2[1] * (ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a2[2] * (-ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) ;
					v2[k][1] = a2[0] * (-ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a2[1] * (ddummy4+SQUARE(ddummy2)*(1.0-ddummy4)) +
							a2[2] * (ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) ;
					v2[k][2] = a2[0] * (ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) +
							a2[1] * (-ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) +
							a2[2] * (ddummy4+SQUARE(ddummy3)*(1.0-ddummy4)) ;
					v3[k][0] = a3[0] * (ddummy4+SQUARE(ddummy1)*(1.0-ddummy4)) +
							a3[1] * (ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a3[2] * (-ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) ;
					v3[k][1] = a3[0] * (-ddummy3*ddummy5+ddummy1*ddummy2*(1.0-ddummy4)) +
							a3[1] * (ddummy4+SQUARE(ddummy2)*(1.0-ddummy4)) +
							a3[2] * (ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) ;
					v3[k][2] = a3[0] * (ddummy2*ddummy5+ddummy1*ddummy3*(1.0-ddummy4)) +
							a3[1] * (-ddummy1*ddummy5+ddummy2*ddummy3*(1.0-ddummy4)) +
							a3[2] * (ddummy4+SQUARE(ddummy3)*(1.0-ddummy4)) ;
				}

				/* in both cases the v_j vectors have the length of a_j */


				/* Calculate volumes of (virtual) unit cells in the crystal, crystal+inner stabilizer layer,
				crystal+total stabilizer layer. Use them for scaling with (differences of) rho's, sld's.
				Calculate the volumes of the crystal, crystal+inner stabilizer layer, crystal+total stabilizer layer.
				Finally calculate the inc. cross sections by scaling with ICSD's  */

				/* the coherent contributions Vrho, Vrho_isl, Vrho_osl, Vsld, Vsld_isl, Vsld_osl
				are obtained by multiplying the respective unit cell volume with the (difference)
				of the electron/scattering length densities. This is due the fact that only
				V/(N1*N2*N3)*rho according to eq. (2) and (3) in the paper are necessary. 
				However the incoherent contribution is obtained by computing the volume of the 	
				complete crystals V and multiplying with the ICSD */

				/* 1. crystal unit cell */
				for ( unsigned int m=0; m<par->nsp; ++m)
				{
					V[m][k] = V_uc ; /* [nm^3] */

					/* make VN the volume of the crystal */
					VN[m][k] = V[m][k] * (double)( n1[k]*n2[k]*(m+1) ) ; /* [nm^3] */
				}


				/* 2. virtual unit cell of crystal k + inner stabilizer shell */
				/* thickness_isl, thickness_osl are the thicknesses of the outer and inner stabilizer layer in nm, 
				inverse projection into direction of the lattice vectors, gives the thickness of
				the stabilizer layer in direction of the a_l for l=1,2,3 */
				/* assume that betrag(v3[k])=betrag(a3[k])=c therefore kf3,okf3 are independent on k */
				/* THE FOLLOWING THREE COMPUTATIONS MIGTH BE PULLED OUT BEFORE THE AVERAGING */
				l1 = thickness_isl / sin (DR * cell_par[3]) ;
				l2 = thickness_isl / sin (DR * cell_par[4]) ;
				l3 = thickness_isl / sin (DR * cell_par[5]) ;
				/* kf_l is a scaling factor with regard to a_l=betrag(v_l[k]) in nm and the number 
				of unit cells in a_l direction, n_l, for l=1,2,3, therefore kf_l*n_l[k]*a_l=n_l[k]*a_l+2*l_l
				describes the length of particle k and its inner stabilizer layer */
				kf1[k] = 1.0 + 2.0 * l1 / (double)( n1[k] ) / betrag(v1[k]) ;
				kf2[k] = 1.0 + 2.0 * l2 / (double)( n2[k] ) / betrag(v2[k]) ;

				/* scale rotated lattice vectors v_l of length a_l with kf_l for l=1,2,3 and
				save them as w_l, where w_l is now the rotated lattice vector of the
				crystalline particle and the inner stabilzer layer, i.e. w_l=a_l+2*l_l/nl[k]
				such that 2*l_l/nl[k] is added onto the lattice constant a_l */
				csp (w1[k], kf1[k], v1[k]) ;
				csp (w2[k], kf2[k], v2[k]) ;
				for ( unsigned int m=0; m<par->nsp; ++m)
				{
					kf3[m] = 1.0 + 2.0 * l3 / ( (double)(m+1)) / betrag(v3[k]) ;
					csp (w3[m][k], kf3[m], v3[k]) ;
					V_isl[m][k] = w1[k][0] * (w2[k][1]*w3[m][k][2] - w2[k][2]*w3[m][k][1]) + 
							w1[k][1] * (w2[k][2]*w3[m][k][0] - w2[k][0]*w3[m][k][2]) + 
							w1[k][2] * (w2[k][0]*w3[m][k][1] - w2[k][1]*w3[m][k][0]) ;

					/* make VN_isl the volume of crystal + inner stabilizer layer */
					VN_isl[m][k] = V_isl[m][k] * (double)( n1[k]*n2[k]*(m+1) ) ; /* [nm^3] */
				}

				/* 3. virtual unit cell of crystal k and the total stabilizer shell */
				/* THE FOLLOWING THREE COMPUTATIONS MIGTH BE PULLED OUT BEFORE THE AVERAGING */
				ol1 = ( thickness_osl + thickness_isl ) / sin (DR * cell_par[3]) ;
				ol2 = ( thickness_osl + thickness_isl ) / sin (DR * cell_par[4]) ;
				ol3 = ( thickness_osl + thickness_isl ) / sin (DR * cell_par[5]) ;
	
				okf1[k] = 1.0 + 2.0 * ol1 / (double)( n1[k] ) / betrag(v1[k]) ;
				okf2[k] = 1.0 + 2.0 * ol2 / (double)( n2[k] ) / betrag(v2[k]) ;

				csp (wt1[k], okf1[k], v1[k]) ;
				csp (wt2[k], okf2[k], v2[k]) ;
				for ( unsigned int m=0; m<par->nsp; ++m)
				{
					okf3[m] = 1.0 + 2.0 * ol3 / ( (double)(m+1)) / betrag(v3[k]) ;
					csp (wt3[m][k], okf3[m], v3[k]) ;
					V_osl[m][k] = wt1[k][0] * (wt2[k][1]*wt3[m][k][2] - wt2[k][2]*wt3[m][k][1]) + 
						wt1[k][1] * (wt2[k][2]*wt3[m][k][0] - wt2[k][0]*wt3[m][k][2]) + 
						wt1[k][2] * (wt2[k][0]*wt3[m][k][1] - wt2[k][1]*wt3[m][k][0]) ;

					/* make VN_osl the volume of crystal + stabilizer layer */
					VN_osl[m][k] = V_osl[m][k] * (double)( n1[k]*n2[k]*(m+1) ) ; /* [nm^3] */
				}


				/* write characteristic data of the k-th cristalline particle to the log file ... */	
				/* for example data of the position and size and rotated lattice vectors (with and without the stabilizer layer) */
				if ( log_flag )
				{
					fprintf( logfile, "%sCrystal %3d: n1= %10d, n2=%10d, n3=%10d\n", "\t", k, n1[k], n2[k], n3[k]) ;
					fprintf( logfile, "%s              D= %15.10G, d1= %15.10G, d2= %15.10G\n", "\t", D[k], d1[k], d2[k]) ;
					fprintf( logfile, "%s              R=(%15.10G, %15.10G, %15.10G)\n", "\t", R[k][0], R[k][1], R[k][2]) ;					
					fprintf( logfile, "%s             r3=%15.10G\n", "\t", r3[k]) ;
					fprintf( logfile, "%s             v1=(%15.10G, %15.10G, %15.10G)\n", "\t", v1[k][0], v1[k][1], v1[k][2]) ;
					fprintf( logfile, "%s             v2=(%15.10G, %15.10G, %15.10G)\n", "\t", v2[k][0], v2[k][1], v2[k][2]) ;
					fprintf( logfile, "%s             v3=(%15.10G, %15.10G, %15.10G)\n", "\t", v3[k][0], v3[k][1], v3[k][2]) ;
					fprintf( logfile, "%s             w1=(%15.10G, %15.10G, %15.10G)\n", "\t", w1[k][0], w1[k][1], w1[k][2]) ;
					fprintf( logfile, "%s             w2=(%15.10G, %15.10G, %15.10G)\n", "\t", w2[k][0], w2[k][1], w2[k][2]) ;
					fprintf( logfile, "%s             w3=(%15.10G, %15.10G, %15.10G)\n", "\t", w3[n3[k]-1][k][0], w3[n3[k]-1][k][1], w3[n3[k]-1][k][2]) ;
					fprintf( logfile, "%s            wt1=(%15.10G, %15.10G, %15.10G)\n", "\t", wt1[k][0], wt1[k][1], wt1[k][2]) ;
					fprintf( logfile, "%s            wt2=(%15.10G, %15.10G, %15.10G)\n", "\t", wt2[k][0], wt2[k][1], wt2[k][2]) ;
					fprintf( logfile, "%s            wt3=(%15.10G, %15.10G, %15.10G)\n", "\t", wt3[n3[k]-1][k][0], wt3[n3[k]-1][k][1], wt3[n3[k]-1][k][2]) ;
					fflush (logfile) ;
				}
	
			/* END SECOND.1 loop of the par->nms particles in the stack */
			}

			/* option to export the generated random numbers, volumes for the current particles by appending them to a .distr file */
			if ( distr_flag )
			{
				sprintf( oname, "%s%s.distr", ofo, mfname) ;
				if ( ( outf = fopen( oname, "a")) == NULL) { XNDIFF_ERROR(18) ; }
	
				for ( unsigned int k=0; k<par->nms; ++k)
				{
					if ( ( par->av_mode == 0 ) || ( par->av_mode == 2 ) || ( par->av_mode == 3 ) ) { fprintf( outf, "  %6d %6d %6d %15.10G %15.10G %15.10G\n", n1[k], n2[k], n3[k], d1[k], d2[k], D[k]) ; }
					else if ( par->av_mode == 1 ) { fprintf( outf, "  %6d %6d %6d %15.10G %15.10G %15.10G %15.10G\n", n1[k], n2[k], n3[k], d1[k], d2[k], D[k], r3[k]) ; }
				}
				fclose(outf) ;
			}



			/*********************************/
			/* VOLUMES FOR CURRENT PARTICLES */
			/*********************************/

			/* After all random numbers and volumes for the par->nms particles in the stack have been computed in the SECOND.1 loop, 
			update V..._irr and V..._dm-variables that are used later for the incoherent scattering contribution and
			for an absolute unit scale, for both single particles and stacks */

			/* VN_irr and VN_dm are used solely for singlefilesoutput, therefore
			for each of the par->nms single particles with thickness (m+1) and
			for each stack the volumes must be computed 
			*/

			/* irradiated sample volume and volume of the dm for single particles ... */
			for ( unsigned int m=0; m<par->nsp; ++m)
			{
				for ( unsigned int k=0; k<par->nms; ++k)
				{
					VN_irr[m][k] = VN[m][k] / conc_cry ; /* [nm^3] */
					VN_dm[m][k] = VN[m][k] * ( conc_dm / conc_cry ) ; /* [nm^3] */
				}
			}

			/* ... and for stacks */
			ddummy1 = VN[n3[0]-1][0] ;
			ddummy2 = VN_isl[n3[0]-1][0] ;
			ddummy3 = VN_osl[n3[0]-1][0] ;
			for ( unsigned int m=par->nsp; m<nspsp; ++m)
			{
				/* idummy =1,2,3,...,(par->nms - 1) for the 2th,3th,4th.,...par->nms th particle for a stack,
				idummy-th particle with random n3[idummy] and idummy-th set of random sizes n1,n2 */	
				idummy = m - par->nsp + 1 ;

				ddummy1 += VN[n3[idummy]-1][idummy] ;
				ddummy2 += VN_isl[n3[idummy]-1][idummy] ;
				ddummy3 += VN_osl[n3[idummy]-1][idummy] ;

				VN[m][0] = ddummy1 ; /* [nm^3] */
				VN_isl[m][0] = ddummy2 ; /* [nm^3] */
				VN_osl[m][0] = ddummy3 ; /* [nm^3] */

				VN_dm[m][0] = ddummy1 * ( conc_dm / conc_cry ) ; /* [nm^3] */

				VN_irr[m][0] = ddummy1 / conc_cry * sc_irr ; /* !!! [nm^2 * cm] !!! */
			}

			/*****************************************/
			/* TOTAL VOLUMES (FOR THE FULL ENSEMBLE) */
			/*****************************************/

			/* single particles, sum over all par->nms */
			/* update Vtot_cry, Vtot_isl, Vtot_osl */
			for ( unsigned int m=0; m<par->nsp; ++m)
			{
				for ( unsigned int k=0; k<par->nms; ++k)
				{
					Vtot_cry[m] += VN[m][k] ; /* [nm^3] */
					Vtot_isl[m] += VN_isl[m][k] ; /* [nm^3] */
					Vtot_osl[m] += VN_osl[m][k] ; /* [nm^3] */
				}
			}
			/* stacks */
			/* update Vtot_cry, Vtot_isl, Vtot_osl */
			/* VN[m=par->nsp...nspsp-1][k=0] are already the correct volumes of the stacks, as computed above !!! */
			for ( unsigned int m=par->nsp; m<nspsp; ++m)
			{
				Vtot_cry[m] += VN[m][0] ; /* [nm^3] */
				Vtot_isl[m] += VN_isl[m][0] ; /* [nm^3] */
				Vtot_osl[m] += VN_osl[m][0] ; /* [nm^3] */
			}

			/* finally update also Vtot_dm and Vtot_irr with conc_dm and conc_cry (both in vol%) */
			/* Vtot_irr important for both, X/n for absolute units, Vtot_dm only for neutron incoherent scattering */
			for ( unsigned int m=0; m<nspsp; ++m)
			{
				Vtot_irr[m] = Vtot_cry[m] / conc_cry * sc_irr ; /* !!! [nm^2 * cm] !!! */
				Vtot_dm[m] = Vtot_cry[m] * ( conc_dm / conc_cry ) ; /* [nm^3] */
			}



			/*********************************************************************/
			/* INCOHERENT SCATTERING FOR CURRENT PARTICLES AND THE FULL ENSEMBLE */
			/*********************************************************************/

			/* update neutron incoherent scattering contribution in B_n as well as in Yi_n */
			if ( Xn <= 0 )
			{
				/* Placing the computation inside the SECOND.2 and THIRD loop
				over polar and azimuthal averaging has the advantage that the 
				MULTIPLICITY applied for coherent scattering is the same as for 
				the incoherent scattering contribution. the computed values of 
				dB and B apply for each single value of s.
				The computation should not be placed inside the FOURTH loop over 
				the s-range next to dE and dS, because incoherent scattering is 
				independent of s (index [i]) and a further normalization by 
				par->np would be necessary.
				Placing here between the FIRST and SECOND.2 loop has the advantage that 
				the computation is minimal faster (saving the Powder Average) but
				different MULTIPLICITY are be necessary due to the missing orientational 
				average.
				*/

				for ( unsigned int q=0; q<sld_multi.n; ++q)
				{
					/* single particles */
					for ( unsigned int m=0; m<par->nsp; ++m)
					{
						for ( unsigned int k=0; k<par->nms; ++k)
						{
							/* sum over incoherent scattering of crystal, outer and inner stabilizer layer 
							and the dispersion medium. */

							/* crystal incoherent scattering */
							dB_n[m][k][q] = VN[m][k] * icsd_cry * sc_icsd ; /* [nm^2] */
							/* add isl's incoherent scattering */
							dB_n[m][k][q] += ( VN_isl[m][k] - VN[m][k] ) * icsd_isl[q] * sc_icsd ; /* [nm^2] */
							/* add osl's incoherent scattering */
							dB_n[m][k][q] += ( VN_osl[m][k] - VN_isl[m][k] ) * icsd_osl[q] * sc_icsd ; /* [nm^2] */
							/* add dm's incoherent scattering */
							dB_n[m][k][q] += VN_dm[m][k] * icsd_dm[q] * sc_icsd ; /* [nm^2] */
						}
					}

					/* stacks, add consecutively the incoherent scattering contributions starting 
					from particle [0] in the stack up to particle par->nms-1 to dBstack_n and
					add dBstack_n to dB */
					/* so the incoherent scattering of the stack is just the sum over the 
					single particles and their proportional dm in the stack*/
					dBstack_n = dB_n[n3[0]-1][0][q] ;
					for ( unsigned int m=par->nsp; m<nspsp; ++m)
					{
						idummy = m - par->nsp + 1 ;
				
						dBstack_n += dB_n[n3[idummy]-1][idummy][q] ;
						// FALSCHES += in dB_n[m][0][q] += dBstack_n ; /* [nm^2] */
						dB_n[m][0][q] = dBstack_n ; /* [nm^2] */
					}
				}

				/* for B_n computations apply for single particles as well as for stacks */
				for ( unsigned int m=0; m<nspsp; ++m)
				{
					for ( unsigned int q=0; q<sld_multi.n; ++q)
					{
						/* first set crystals incoherent scattering */
						B_n[m][q] = icsd_cry * Vtot_cry[m] * sc_icsd ; /* [nm^2] */
						/* add isl's incoherent scattering */
						B_n[m][q] += icsd_isl[q] * ( Vtot_isl[m] - Vtot_cry[m] ) * sc_icsd ; /* [nm^2] */
						/* add osl's incoherent scattering */
						B_n[m][q] += icsd_osl[q] * ( Vtot_osl[m] - Vtot_isl[m] ) * sc_icsd ; /* [nm^2] */
						/* add dm's incoherent scattering */
						B_n[m][q] += icsd_dm[q] * Vtot_dm[m] * sc_icsd ; /* [nm^2] */
					}
				}

				/* Yi_n computations applies for single particles and stacks, use total volumes from above */
				for ( unsigned int m=0; m<nspsp; ++m)
				{
					/* crystals incoherent scattering */
					Yi_n[m][0] = icsd_cry * Vtot_cry[m] * sc_icsd ; /* [nm^2] */
					/* isl's incoherent scattering */
					Yi_n[m][1] = ( Vtot_isl[m] - Vtot_cry[m] ) * sc_icsd ; /* [nm^2 * cm] */
					/* osl's incoherent scattering */
					Yi_n[m][2] = ( Vtot_osl[m] - Vtot_isl[m] ) * sc_icsd ; /* [nm^2 * cm] */
					/* dm's incoherent scattering */
					Yi_n[m][3] = Vtot_dm[m] * sc_icsd ; /* [nm^2 * cm] */
				}
			}

			/*******************************************************************************************/
			/* ORIENTATIONAL AVERAGE / COHERENT SCATTERING FOR CURRENT PARTICLES AND THE FULL ENSEMBLE */
			/*******************************************************************************************/
			/* For security default(none) clause is used -> each variable in the parallel block must be declared explicitely 
			   either as private or shared (without clause default would be shared) -> prevent errors .
			   Outer parallel loop variable must not be declared private, but nested (inner) loop counters have to be declared private in C(++) !!!
			   Similar settings are applied in compute_structure_amplitude().
			*/
			#pragma omp parallel private( ii, jj, ll, mm, pp, qq, av_ac, sinav_ac, cosav_ac, av_pc, sinav_pc, cosav_pc, weight_factor, ddummy1, ddummy2, ddummy3, G, cs, s, piQv1, piQv2, piQv3, piQw1, piQw2, piQw3, piQwt1, piQwt2, piQwt3, exp_Qv1, exp_Qv2, exp_Qv3, exp_N1Qv1, exp_N2Qv2, exp_N3Qv3, exp_N1Qw1, exp_N2Qw2, exp_N3Qw3, exp_N1Qwt1, exp_N2Qwt2, exp_N3Qwt3, GA_N1_N2, cdummy11, cdummy12, cdummy13, cdummy21, cdummy22, cdummy23, cdummy31, cdummy32, cdummy33, cdummy41, cdummy42, cdummy43, GA, P, P_isl, P_osl, FFF_cdummy, Bg_X, dE_X, Bg_n, dE_n, idummy, exp_RQ, k_phi, dEstack_X, dEstack_n, dYc_X_stack, dYc_n_stack, exp_DRQ, dS_X_parallel_thread, dS_n_parallel_thread, Yc_X_parallel_thread, Yc_n_parallel_thread) shared( dS_X, dS_n, Yc_X, Yc_n, Yi_n, kk_jj_count, nspsp, rtstr, pr, c_rep, stdout, ind_r3, v1, v2, v3, kf1, kf2, kf3, okf1, okf2, okf3, n1, n2, n3, R, V, V_isl, V_osl, drerho_dm_osl, drerho_osl_isl, drerho_isl, dsld_dm_osl, dsld_osl_isl, dsld_isl) reduction(+:single_sum_weight_factor) reduction(+:sum_weight_factor) default(none) if ( openmp_flag )
			{
				/* prepare parallel region, allocate dynamic arrays that should be used as private variables
				   GA, P, P_isl, P_osl, dE_X, dE_n, exp_N3Qv3, exp_N3Qw3, exp_N3Qwt3, cdummy41, cdummy42, cdummy43, exp_RQ,
				   dS_X_parallel_thread, dS_n_parallel_thread, Yc_X_parallel_thread, Yc_n_parallel_thread
				   Note that variables allocated/defined within a parallel section are by default private.
				   dE_X/n and related variables are private for each thread, their data will be written to newly introduced private dummy arrays like
				   dS_X/n_parallel_thread, dYc_X/n_parallel_thread and related ones. They will be at the end written to dS_X/n, dYc_X/n etc.
				*/
				GA = (dcmplx **) calloc( par->nsp, sizeof(dcmplx *)) ;
				P = (dcmplx **) calloc( par->nsp, sizeof(dcmplx *)) ;
				P_isl = (dcmplx **) calloc( par->nsp, sizeof(dcmplx *)) ;
				P_osl = (dcmplx **) calloc( par->nsp, sizeof(dcmplx *)) ;
				for ( mm=0; mm<par->nsp; ++mm) 
				{
					GA[mm] = (dcmplx *) calloc( par->nms, sizeof(dcmplx)) ;
					P[mm] = (dcmplx *) calloc( par->nms, sizeof(dcmplx)) ;
					P_isl[mm] = (dcmplx *) calloc( par->nms, sizeof(dcmplx)) ;
					P_osl[mm] = (dcmplx *) calloc( par->nms, sizeof(dcmplx)) ;
				}

				/* dE_X/n is a dcmplx matrix of size ( par->nsp x par->nms x rho/sld_multi.n ) */
				if ( Xn >= 0 ) { dE_X = (dcmplx ***) calloc( par->nsp, sizeof(dcmplx **)) ; }
				if ( Xn <= 0 ) { dE_n = (dcmplx ***) calloc( par->nsp, sizeof(dcmplx **)) ; }

				for ( mm=0; mm<par->nsp; ++mm) 
				{
					if ( Xn >= 0 )
					{
						dE_X[mm] = (dcmplx **) calloc( par->nms, sizeof(dcmplx *)) ;
						for ( pp=0; pp<par->nms; ++pp) { dE_X[mm][pp] = (dcmplx *) calloc( rho_multi.n, sizeof(dcmplx)) ; }
					}
					if ( Xn <= 0 )
					{
						dE_n[mm] = (dcmplx **) calloc( par->nms, sizeof(dcmplx *)) ;
						for ( pp=0; pp<par->nms; ++pp) { dE_n[mm][pp] = (dcmplx *) calloc( sld_multi.n, sizeof(dcmplx)) ; }
					}
				}

				exp_N3Qv3 = (dcmplx *) calloc( par->nsp, sizeof(dcmplx)) ;
				exp_N3Qw3 = (dcmplx *) calloc( par->nsp, sizeof(dcmplx)) ;
				exp_N3Qwt3 = (dcmplx *) calloc( par->nsp, sizeof(dcmplx)) ;
				cdummy41 = (dcmplx *) calloc( par->nsp, sizeof(dcmplx)) ;
				cdummy42 = (dcmplx *) calloc( par->nsp, sizeof(dcmplx)) ;
				cdummy43 = (dcmplx *) calloc( par->nsp, sizeof(dcmplx)) ;

				exp_RQ = (dcmplx *) calloc( par->nms, sizeof(dcmplx)) ;

				/* dS_X/n_parallel_thread*/
				if ( Xn >= 0 ) { dS_X_parallel_thread = (double ****) calloc( nspsp, sizeof(double ***)) ; }
				if ( Xn <= 0 ) { dS_n_parallel_thread = (double ****) calloc( nspsp, sizeof(double ***)) ; }
				for ( mm=0; mm<nspsp; ++mm)
				{
					if ( mm < par->nsp )
					{
						/* single particle */
						if ( Xn >= 0 )
						{
							dS_X_parallel_thread[mm] = (double ***) calloc( par->nms, sizeof(double **)) ;
							for ( pp=0; pp<par->nms; ++pp)
							{
								dS_X_parallel_thread[mm][pp] = (double **) calloc( rho_multi.n, sizeof(double *)) ;
								for ( qq=0; qq<rho_multi.n; ++qq) { dS_X_parallel_thread[mm][pp][qq] = (double *) calloc( par->np, sizeof(double)) ; }
							}
						}
						if ( Xn <= 0 )
						{
							dS_n_parallel_thread[mm] = (double ***) calloc( par->nms, sizeof(double **)) ;
							for ( pp=0; pp<par->nms; ++pp)
							{
								dS_n_parallel_thread[mm][pp] = (double **) calloc( sld_multi.n, sizeof(double *)) ;
								for ( qq=0; qq<sld_multi.n; ++qq) { dS_n_parallel_thread[mm][pp][qq] = (double *) calloc( par->np, sizeof(double)) ; }
							}
						}
					}
					else
					{	/* dS stacks */
						if ( Xn >= 0 )
						{
							dS_X_parallel_thread[mm] = (double ***) calloc( 1, sizeof(double **)) ;
							dS_X_parallel_thread[mm][0] = (double **) calloc( rho_multi.n, sizeof(double *)) ;
							for ( qq=0; qq<rho_multi.n; ++qq) { dS_X_parallel_thread[mm][0][qq] = (double *) calloc( par->np, sizeof(double)) ; }

						}
						if ( Xn <= 0 )
						{
							dS_n_parallel_thread[mm] = (double ***) calloc( 1, sizeof(double **)) ;
							dS_n_parallel_thread[mm][0] = (double **) calloc( sld_multi.n, sizeof(double *)) ;
							for ( qq=0; qq<sld_multi.n; ++qq) { dS_n_parallel_thread[mm][0][qq] = (double *) calloc( par->np, sizeof(double)) ; }
						}
					}
				}

				/* Yc_X/n_parallel_thread */
				if ( Xn >= 0 ) { Yc_X_parallel_thread = (double ***) calloc( nspsp, sizeof(double **)) ; }
				if ( Xn <= 0 ) { Yc_n_parallel_thread = (double ***) calloc( nspsp, sizeof(double **)) ; }
				for ( mm=0; mm<nspsp; ++mm)
				{
					if ( Xn >= 0 )
					{
						Yc_X_parallel_thread[mm] = (double **) calloc( 10, sizeof(double *)) ;
						for ( pp=0; pp<10; ++pp) { Yc_X_parallel_thread[mm][pp] = (double *) calloc( par->np, sizeof(double)) ; }
					}
					if ( Xn <= 0 )
					{
						Yc_n_parallel_thread[mm] = (double **) calloc( 10, sizeof(double *)) ;
						for ( pp=0; pp<10; ++pp) { Yc_n_parallel_thread[mm][pp] = (double *) calloc( par->np, sizeof(double)) ; }
					}
				}

				/* for time display */
				kk_jj_count = 0 ;

				/* SECOND.2 loop over azimuthal range (cycled only once if par->av_mode==1, i.e. r3-average )*/
				#pragma omp for schedule(static)
				for ( kk=0; kk<nav_ac; ++kk)
				{
					/* GA, P_osl, P_isl, P, cdummy41, cdummy42, cdummy43, exp_N3Qv3, exp_N3Qw3, exp_N3Qwt3, dE_X, dE_n will be overwritten */
					/* dS_X_parallel_thread, dS_n_parallel_thread, Yc_X_parallel_thread, Yc_n_parallel_thread will be updated (+=) within each thread and finally written to the corresponding variables */

					/* reset variables to 0.0 */
// 					/* par->nsp x par->nms */
// 					for ( mm=0; mm<par->nsp; ++mm)
// 					{
// 						for ( pp=0; pp<par->nms; ++pp)
// 						{
// 							GA[mm][pp] = dcmplx( 0.0, 0.0) ;
// 							P[mm][pp] = dcmplx( 0.0, 0.0) ;
// 							P_isl[mm][pp] = dcmplx( 0.0, 0.0) ;
// 							P_osl[mm][pp] = dcmplx( 0.0, 0.0) ;
// 						}
// 					}
// 
// 					/* dE_X/n which are a double complex matrices of size ( par->nsp x par->nms x rho/sld_multi.n ) */
// 					if ( Xn >= 0 )
// 					{
// 						for ( mm=0; mm<par->nsp; ++mm) 
// 						{
// 							for ( pp=0; pp<par->nms; ++pp) 
// 							{ 
// 								for ( qq=0; qq<rho_multi.n; ++qq) 
// 								{ 
// 									dE_X[mm][pp][qq] = dcmplx( 0.0, 0.0) ;
// 								}
// 							}
// 						}
// 					}
// 					if ( Xn <= 0 )
// 					{
// 						for ( mm=0; mm<par->nsp; ++mm) 
// 						{
// 							for ( pp=0; pp<par->nms; ++pp) 
// 							{ 
// 								for ( qq=0; qq<sld_multi.n; ++qq) 
// 								{ 
// 									dE_n[mm][pp][qq] = dcmplx( 0.0, 0.0) ;
// 								}
// 							}
// 						}
// 					}
// 
// 					for ( mm=0; mm<par->nsp; ++mm) 
// 					{
// 						exp_N3Qv3[mm] = dcmplx( 0.0, 0.0) ;
// 						exp_N3Qw3[mm] = dcmplx( 0.0, 0.0) ;
// 						exp_N3Qwt3[mm] = dcmplx( 0.0, 0.0) ;
// 						cdummy41[mm] = dcmplx( 0.0, 0.0) ;
// 						cdummy42[mm] = dcmplx( 0.0, 0.0) ;
// 						cdummy43[mm] = dcmplx( 0.0, 0.0) ;
// 					}
// 
// 					for ( pp=0; pp<par->nms; ++pp) 
// 					{ 
// 						exp_RQ[pp] = dcmplx( 0.0, 0.0) ;
// 					}
// 
// 					/* free dS_X/n_thread_parallel */
// 					if ( Xn >= 0 )
// 					{
// 						for ( mm=0; mm<nspsp; ++mm)
// 						{
// 							if ( mm<par->nsp )
// 							{
// 								for ( pp=0; pp<par->nms; ++pp)
// 								{
// 									for ( qq=0; qq<rho_multi.n; ++qq) { free(dS_X_parallel_thread[mm][pp][qq]) ; }
// 									free(dS_X_parallel_thread[mm][pp]) ;
// 								}
// 							}
// 							else
// 							{
// 								for ( qq=0; qq<rho_multi.n; ++qq) { free(dS_X_parallel_thread[mm][0][qq]) ; }
// 								free(dS_X_parallel_thread[mm][0]) ;
// 							}
// 							free(dS_X_parallel_thread[mm]) ;
// 						}
// 						free(dS_X_parallel_thread) ;
// 					}
// 					if ( Xn <= 0 )
// 					{
// 						for ( mm=0; mm<nspsp; ++mm)
// 						{
// 							if ( mm<par->nsp )
// 							{
// 								for ( pp=0; pp<par->nms; ++pp)
// 								{
// 									for ( qq=0; qq<sld_multi.n; ++qq) { free(dS_n_parallel_thread[mm][pp][qq]) ; }
// 									free(dS_n_parallel_thread[mm][pp]) ;
// 								}
// 							}
// 							else
// 							{
// 								for ( qq=0; qq<sld_multi.n; ++qq) { free(dS_n_parallel_thread[mm][0][qq]) ; }
// 								free(dS_n_parallel_thread[mm][0]) ;
// 							}
// 							free(dS_n_parallel_thread[mm]) ;
// 						}
// 						free(dS_n_parallel_thread) ;
// 					}
// 
// 					/* free Yc_X/n_parallel_thread */
// 					if ( Xn >= 0 )
// 					{
// 						for ( mm=0; mm<nspsp; ++mm) 
// 						{
// 							for ( pp=0; pp<10; ++pp) { free(Yc_X_parallel_thread[mm][pp]) ; }
// 							free(Yc_X_parallel_thread[mm]) ;
// 						}
// 						free(Yc_X_parallel_thread);
// 					}
// 					if ( Xn <= 0 )
// 					{
// 						for ( mm=0; mm<nspsp; ++mm) 
// 						{
// 							for ( pp=0; pp<10; ++pp) { free(Yc_n_parallel_thread[mm][pp]) ; }
// 							free(Yc_n_parallel_thread[mm]) ;
// 						}
// 						free(Yc_n_parallel_thread);
// 					}




					/* Grid based Powder Averages */
					if ( ( par->av_mode == 0 ) || ( par->av_mode == 1 ) )
					{
						av_ac = av_azi_ang.at(kk) ;
						sinav_ac = sin (DR * av_ac) ;
						cosav_ac = cos (DR * av_ac) ;
					}

					// av_pc = av_p1 - av_pd/2.0 ;
					/* THIRD loop over polar angle range */
					for ( jj=0; jj<nav_pc; ++jj)
					{
						kk_jj_count += 1 ;
						if ( time_flag > 1 ) 
						{
							/* update only for master thread */
// 							if ( omp_get_thread_num() == 0 )
// 							{
								fprintf( stdout ,"PID %d | R %4d | %5.1f %% | T %5.1f %% | %-30s\r", pid, c_rep, 100.0 * ((double) ( kk_jj_count ) )/((double) (( nav_pc )*( nav_ac )) ), pr, rtstr ) ; fflush( stdout ) ;
// 							}
						}

						/* MC mode */
						if ( ( par->av_mode == 2 ) || ( par->av_mode == 3 ) )
						{
							av_ac = av_azi_ang.at(jj) ;
							sinav_ac = sin (DR * av_ac) ;
							cosav_ac = cos (DR * av_ac) ;
						}

						av_pc = av_pol_ang.at(jj) ;
						sinav_pc = sin (DR * av_pc) ;
						cosav_pc = cos (DR * av_pc) ;


						/* optionally apply a weighting factor ~ sin / (\sum sin ) */
						if ( weight_flag ) { weight_factor = sinav_pc ; }
						else { weight_factor = 1.0 ; }
						single_sum_weight_factor += weight_factor ;
						sum_weight_factor += weight_factor ;


						/* set up G in the kg1,kg2,G0 system, if G!=0, normalize G to a (dimensionless) unit vector */
						ddummy1 = sinav_pc * cosav_ac ;
						ddummy2 = sinav_pc * sinav_ac ;
						ddummy3 = cosav_pc ;
						G[0] = ddummy1 * kg1[0] + ddummy2 * kg2[0] + ddummy3 * G0[0] ;
						G[1] = ddummy1 * kg1[1] + ddummy2 * kg2[1] + ddummy3 * G0[1] ;
						G[2] = ddummy1 * kg1[2] + ddummy2 * kg2[2] + ddummy3 * G0[2] ;
						ddummy1 = betrag (G) ;
						if (ddummy1 > 1e-10)
						{
							G[0] /= ddummy1 ;
							G[1] /= ddummy1 ;
							G[2] /= ddummy1 ;
						}

						/* s1 = par->rs1 - par->ds ! */
						cs = s1 ;
						/* FOURTH loop over the s-range */
						for ( ii=0; ii<par->np; ++ii)
						{
							cs += ds ;

							/* vector s will be used for stacks <s,R> */
							s[0] = G[0] * cs ;
							s[1] = G[1] * cs ;
							s[2] = G[2] * cs ;

							if ( ( par->av_mode == 0 ) || ( par->av_mode == 2 ) || ( par->av_mode == 3 ) )
							{
								/* use index kk in eQvj[jj][kk] due to azimuthal averaging */
								/* over all n1[p], n2[p] from p=1...nms */

								piQv1 = M_PI * cs * eQvj[jj][kk][0] ;
								piQv2 = M_PI * cs * eQvj[jj][kk][1] ;
								piQv3 = M_PI * cs * eQvj[jj][kk][2] ;

								/* can be pulled out from p-for, but only for par->av_mode=0 */
								exp_Qv1 = exp ( dcmplx( 0.0, 2.0 * piQv1 ) ) - 1.0 ;
								exp_Qv2 = exp ( dcmplx( 0.0, 2.0 * piQv2 ) ) - 1.0 ;
								exp_Qv3 = exp ( dcmplx( 0.0, 2.0 * piQv3 ) ) - 1.0 ;


								/* can be pulled out from pp-for, but only for par->av_mode=0 */
								for ( mm=0; mm<par->nsp; ++mm)
								{
									piQwt3 = okf3[mm] * piQv3 ;
									piQw3 = kf3[mm] * piQv3 ;

									exp_N3Qv3[mm] = exp ( dcmplx( 0.0, 2.0 * piQv3 * (double)(mm+1) ) ) - 1.0 ;
									/* w_j = kf_j * v_j */
									exp_N3Qw3[mm] = exp ( dcmplx( 0.0, 2.0 * piQw3 * (double)(mm+1) ) ) - 1.0 ;
									/* wt_j = okf_j * v_j */
									exp_N3Qwt3[mm] = exp ( dcmplx( 0.0, 2.0 * piQwt3 * (double)(mm+1)  ) ) - 1.0 ;

									if (fabs(piQwt3) > 1e-10) { cdummy41[mm] = exp_N3Qwt3[mm] / ( dcmplx( 0.0, 2.0 * piQwt3) ) ; }
									else { cdummy41[mm] = dcmplx( (double)(mm+1), 0.0) ; }

									if (fabs(piQw3) > 1e-10) { cdummy42[mm] = exp_N3Qw3[mm] / ( dcmplx( 0.0, 2.0 * piQw3) ); }
									else { cdummy42[mm] = dcmplx( (double)(mm+1), 0.0) ; }

									if (fabs(piQv3) > 1e-10) { cdummy43[mm] = exp_N3Qv3[mm] / ( dcmplx( 0.0, 2.0 * piQv3) ); }
									else { cdummy43[mm] = dcmplx( (double)(mm+1), 0.0) ; }
								}

								/* over all par->nms N1,N2 configurations */
								for ( pp=0; pp<par->nms; pp++)
								{
									/* compute FFF each time but still save it in FFF[ii][jj][kk], due to Grid / MC Grid all entries of FFF are filled  */
									if ( FNoMem_flag )
									{
										if ( Xn >= 0 ) { FFF_X[ii][jj][kk] = dcmplx(0.0,0.0) ; }
										if ( Xn <= 0 ) { FFF_n[ii][jj][kk] = dcmplx(0.0,0.0) ; }

										for ( ll=0; ll<noa; ++ll)
										{
											FFF_cdummy = exp( dcmplx( 0.0, 2.0 * M_PI * ( s[0]*(atom[ll].coord[0]*v1[pp][0]+atom[ll].coord[1]*v2[pp][0]+atom[ll].coord[2]*v3[pp][0]) +
											s[1]*(atom[ll].coord[0]*v1[pp][1]+atom[ll].coord[1]*v2[pp][1]+atom[ll].coord[2]*v3[pp][1]) +
											s[2]*(atom[ll].coord[0]*v1[pp][2]+atom[ll].coord[1]*v2[pp][2]+atom[ll].coord[2]*v3[pp][2]) ) ) ) ;

											FFF_cdummy *= atom[ll].mult ;

											/* compute structure amplitude F_k for X-ray with the atomic scattering factors */
											if ( Xn >= 0 ) { FFF_X[ii][jj][kk] += FFF_cdummy * f0table[atom[ll].f0ind][ii] * ( (dcmplx) sc_rho ) ; } /* [nm] */

											/* the same for neutrons with the coherent neutron scattering lengths
											   don't use any minus sign convention for b_coh, use just b_coh,
											   the sld will also not include later minus signs
											*/
											if ( Xn <= 0 ) { FFF_n[ii][jj][kk] += FFF_cdummy * atom[ll].coh_sc_length * ( (dcmplx) sc_scl ) ; } /* [nm] */

											if ( test_XnEquiv_flag ) { FFF_n[ii][jj][kk] = FFF_X[ii][jj][kk] ; }
										}
									}

									piQwt1 = okf1[pp] * piQv1 ;
									piQwt2 = okf2[pp] * piQv2 ;
		
									piQw1 = kf1[pp] * piQv1 ;
									piQw2 = kf2[pp] * piQv2 ;

									exp_N1Qv1 = exp ( dcmplx( 0.0, 2.0 * piQv1 * (double)n1[pp]) ) - 1.0 ;
									exp_N2Qv2 = exp ( dcmplx( 0.0, 2.0 * piQv2 * (double)n2[pp]) ) - 1.0 ;

									/* w_j = kf_j * v_j */
									exp_N1Qw1 = exp ( dcmplx( 0.0, 2.0 * piQw1 * (double)n1[pp] ) ) - 1.0 ;
									exp_N2Qw2 = exp ( dcmplx( 0.0, 2.0 * piQw2 * (double)n2[pp] ) ) - 1.0 ;
									/* wt_j = okf_j * v_j */
									exp_N1Qwt1 = exp ( dcmplx( 0.0, 2.0 * piQwt1 * (double)n1[pp] ) ) - 1.0 ;
									exp_N2Qwt2 = exp ( dcmplx( 0.0, 2.0 * piQwt2 * (double)n2[pp] ) ) - 1.0 ;

									/* for GA */
									GA_N1_N2 = dcmplx( 1.0, 0.0) ;

									if ( fabs( real( exp_Qv1 * conj(exp_Qv1) ) ) > 1e-10) { GA_N1_N2 *= ( exp_N1Qv1 / exp_Qv1 ) ; }
									else { GA_N1_N2 *= dcmplx( ( (double)n1[pp]), 0.0) ; }

									if ( fabs( real( exp_Qv2 * conj(exp_Qv2) ) ) > 1e-10) { GA_N1_N2 *= ( exp_N2Qv2 / exp_Qv2 ) ; }
									else { GA_N1_N2 *= dcmplx( ( (double)n2[pp]), 0.0) ; }

									/* for background and dispersion medium S_D */
									if (fabs(piQwt1) > 1e-10) { cdummy21 = exp_N1Qwt1 / ( dcmplx( 0.0, 2.0 * piQwt1 ) ) ; }
									else { cdummy21 = dcmplx( ( (double)n1[pp]), 0.0) ; }

									if (fabs(piQwt2) > 1e-10) { cdummy31 = exp_N2Qwt2 / ( dcmplx( 0.0, 2.0 * piQwt2) ) ; }
									else { cdummy31 = dcmplx( ( (double)n2[pp]), 0.0) ; }


									if (fabs(piQw1) > 1e-10) { cdummy22 = exp_N1Qw1 / ( dcmplx( 0.0, 2.0 * piQw1) ); }
									else { cdummy22 = dcmplx( ( (double)n1[pp]), 0.0) ; }

									if (fabs(piQw2) > 1e-10) { cdummy32 = exp_N2Qw2 / ( dcmplx( 0.0, 2.0 * piQw2) ); }
									else { cdummy32 = dcmplx( ( (double)n2[pp]), 0.0) ; }


									if (fabs(piQv1) > 1e-10) { cdummy23 = exp_N1Qv1 / ( dcmplx( 0.0, 2.0 * piQv1) ); }
									else { cdummy23 = dcmplx( ( (double)n1[pp]), 0.0) ; }

									if (fabs(piQv2) > 1e-10) { cdummy33 = exp_N2Qv2 / ( dcmplx( 0.0, 2.0 * piQv2) ); }
									else { cdummy33 = dcmplx( ( (double)n2[pp]), 0.0) ; }


									/* compute scattering contribution for all single nanoparticles with n3 ~ m + 1 = 1...nsp */
									/* save in dE[mm][pp] and add |dE[mm][pp]|^2 to dS[ii][mm][pp] and S[ii][mm] */
									for ( mm=0; mm<par->nsp; ++mm)
									{
										/* lattice amplitude GA */
										/* GA=(exp_N1Qv1/exp_Qv1)*(exp_N2Qv2/exp_Qv2)*(exp_N3Qv3/exp_Qv3) */

										if ( fabs( real( exp_Qv3 * conj(exp_Qv3) ) ) > 1e-10) { GA[mm][pp] = GA_N1_N2 * ( exp_N3Qv3[mm] / exp_Qv3 ) ; }
										else { GA[mm][pp] = GA_N1_N2 * dcmplx( (double)(mm+1), 0.0) ; }

										/* Scattering of excluded dispersion medium and stabilizer layer */
										/* according to eq.(3)  S_D^k=S^{k''}(rho_dm)-S^{k''}(rho_osl) + S^{k'}(rho_osl)-S^{k'}(rho_isl) + S^{k}(rho_isl) */
										/* where S^{k}(rho) is given as in eq.(2)  S^{k}(rho)=V_u^{k}*rho*[ \prod_{j=1}^3 (exp(I*N_j^k*<a_j^k,Q>)-1)/(I*<a_j^k,Q>)] */
										/* at this V_u^{k} is the volume of the unit cell and a_j^{k} the lattice vectors see table below and */
										/* k ~ crystal, k' ~ crystal with inner stabilizer layer, k'' ~ crystal with total stabilizer layer */
										/*  ___________________________  */
										/* |   |  V_u^{k}  |  a_j^{k}  | */
										/* |___|___________|___________| */
										/* |k  |  V        |  v_j[k]   | */
										/* |k' |  V_sl[k]  |  w_j[k]   | */
										/* |k''|  V_osl[k] |  wt_j[k]  | */
										/* |___|___________|___________| */
										/* the 5 terms in eq.(3) can be conflated into 3 ones, due to two pairs of terms with k'' and k', i.e. */
										/* S_D^k=S^{k''}(rho_dm-rho_osl) + S^{k'}(rho_osl-rho_isl) + S^{k}(rho_isl) */
										/* compute in 3 steps Bg=S_D^k */

										cdummy11 = cdummy21 * cdummy31 * cdummy41[mm] ;
										cdummy12 = cdummy22 * cdummy32 * cdummy42[mm] ;
										cdummy13 = cdummy23 * cdummy33 * cdummy43[mm] ;

										P_osl[mm][pp] = cdummy11 * ((dcmplx) V_osl[mm][pp]) ; /* [nm^3] */
										P_isl[mm][pp] = cdummy12 * ((dcmplx) V_isl[mm][pp]) ;
										P[mm][pp] = cdummy13 * ((dcmplx) V[mm][pp]) ;

										if ( Xn >= 0 )
										{
											for ( qq=0; qq<rho_multi.n; ++qq)
											{
												Bg_X = dcmplx(0.0,0.0) ;  /* [nm] */
												Bg_X += ((dcmplx) drerho_dm_osl[qq] ) * P_osl[mm][pp] ;
												Bg_X += ((dcmplx) drerho_osl_isl[qq] ) * P_isl[mm][pp] ;
												Bg_X += ((dcmplx) drerho_isl[qq] ) * P[mm][pp] ;

												dE_X[mm][pp][qq] = FFF_X[ii][jj][kk] * GA[mm][pp] - Bg_X ;  /* [nm] */
			
	// 											dS_X[m][p][q][i] += fabs( real( dE_X[m][p][q] * conj(dE_X[m][p][q]) ) ) ; /* [nm^2] */
												dS_X_parallel_thread[mm][pp][qq][ii] += cabs2( dE_X[mm][pp][qq] ) * weight_factor ; /* [nm^2] */
											}
										}

										if ( Xn <= 0 )
										{
											for ( qq=0; qq<sld_multi.n; ++qq)
											{
												Bg_n = dcmplx(0.0,0.0) ; /* [nm] */
												Bg_n += ((dcmplx) dsld_dm_osl[qq] ) * P_osl[mm][pp] ;
												Bg_n += ((dcmplx) dsld_osl_isl[qq] ) * P_isl[mm][pp] ;
												Bg_n += ((dcmplx) dsld_isl[qq] ) * P[mm][pp] ;

												dE_n[mm][pp][qq] = FFF_n[ii][jj][kk] * GA[mm][pp] - Bg_n ; /* [nm] */

												dS_n_parallel_thread[mm][pp][qq][ii] += cabs2( dE_n[mm][pp][qq] ) * weight_factor ; /* [nm^2] */
											}
										}

										/* add-up the terms for Yc_X/n[mm][0...9][i] and for Yi_n[mm][0...3] */
										/* !!! it is summed automatically over all pp !!! */
										/* note that Yc_X[mm][qq][ii] ~ Yc_n[mm][qq][ii] for qq=1...6, because F isn't involved into it */
										if ( Xn >= 0 )
										{
											Yc_X_parallel_thread[mm][0][ii] += cabs2( FFF_X[ii][jj][kk] * GA[mm][pp] ) * weight_factor ; /* [nm^2] */

											Yc_X_parallel_thread[mm][1][ii] += cabs2( P_osl[mm][pp] ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
											Yc_X_parallel_thread[mm][2][ii] += cabs2( P_isl[mm][pp] ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
											Yc_X_parallel_thread[mm][3][ii] += cabs2( P[mm][pp] ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */

											Yc_X_parallel_thread[mm][4][ii] += 2.0 * real ( P_osl[mm][pp] * conj( P_isl[mm][pp] ) ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
											Yc_X_parallel_thread[mm][5][ii] += 2.0 * real ( P_osl[mm][pp] * conj( P[mm][pp] ) ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
											Yc_X_parallel_thread[mm][6][ii] += 2.0 * real ( P_isl[mm][pp] * conj( P[mm][pp] ) ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */

											Yc_X_parallel_thread[mm][7][ii] += -2.0 * real ( FFF_X[ii][jj][kk] * GA[mm][pp] * conj( P_osl[mm][pp] ) ) * sc_rho * weight_factor ; /* [nm^5] */
											Yc_X_parallel_thread[mm][8][ii] += -2.0 * real ( FFF_X[ii][jj][kk] * GA[mm][pp] * conj( P_isl[mm][pp] ) ) * sc_rho * weight_factor ; /* [nm^5] */ 
											Yc_X_parallel_thread[mm][9][ii] += -2.0 * real ( FFF_X[ii][jj][kk] * GA[mm][pp] * conj( P[mm][pp] ) ) * sc_rho * weight_factor ; /* [nm^5] */
										}
										if ( Xn <= 0 )
										{
											Yc_n_parallel_thread[mm][0][ii] += cabs2( FFF_n[ii][jj][kk] * GA[mm][pp] ) * weight_factor ; /* [nm^2] */

											Yc_n_parallel_thread[mm][1][ii] += cabs2( P_osl[mm][pp] ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
											Yc_n_parallel_thread[mm][2][ii] += cabs2( P_isl[mm][pp] ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
											Yc_n_parallel_thread[mm][3][ii] += cabs2( P[mm][pp] ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */

											Yc_n_parallel_thread[mm][4][ii] += 2.0 * real ( P_osl[mm][pp] * conj( P_isl[mm][pp] ) ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
											Yc_n_parallel_thread[mm][5][ii] += 2.0 * real ( P_osl[mm][pp] * conj( P[mm][pp] ) ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
											Yc_n_parallel_thread[mm][6][ii] += 2.0 * real ( P_isl[mm][pp] * conj( P[mm][pp] ) ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */

											Yc_n_parallel_thread[mm][7][ii] += -2.0 * real ( FFF_n[ii][jj][kk] * GA[mm][pp] * conj( P_osl[mm][pp] ) ) * sc_sld * weight_factor ; /* [nm^4] */
											Yc_n_parallel_thread[mm][8][ii] += -2.0 * real ( FFF_n[ii][jj][kk] * GA[mm][pp] * conj( P_isl[mm][pp] ) ) * sc_sld * weight_factor ; /* [nm^4] */
											Yc_n_parallel_thread[mm][9][ii] += -2.0 * real ( FFF_n[ii][jj][kk] * GA[mm][pp] * conj( P[mm][pp] ) ) * sc_sld * weight_factor ; /* [nm^4] */
										}
									} /* mm = 0 ... par->nsp-1 */
								} /* pp = 0 ... par->nms-1 */
							}
							else if ( par->av_mode == 1 )
							{
								/* av_mode==1 random r3-rotation */
								/* use index ind_r3[p] in eQvj[j][ind_r3[p]] due to random r3-rotation averaging */
								/* now the particles have distinct v_j and so on, due to random r3-rotation ! */
								/* not so much computations can be pulled out to simplify in this mode */
								for ( pp=0; pp<par->nms; pp++)
								{
									/* compute FFF each time but still save it in FFF[i][j][k], only par->nms r3-angles are used for one c_rep step, later FFF_X/n[i][j][ind_r3[p]] */
									/* if some of the par->nms r3 angles have a multiple occurence FFF will just be overwritten with the same value ! */
									if ( FNoMem_flag )
									{
										if ( Xn >= 0 ) { FFF_X[ii][jj][ind_r3[pp]] = dcmplx(0.0,0.0) ; }
										if ( Xn <= 0 ) { FFF_n[ii][jj][ind_r3[pp]] = dcmplx(0.0,0.0) ; }

										for ( ll=0; ll<noa; ++ll)
										{
											FFF_cdummy = exp( dcmplx( 0.0, 2.0 * M_PI * ( s[0]*(atom[ll].coord[0]*v1[pp][0]+atom[ll].coord[1]*v2[pp][0]+atom[ll].coord[2]*v3[pp][0]) +
											s[1]*(atom[ll].coord[0]*v1[pp][1]+atom[ll].coord[1]*v2[pp][1]+atom[ll].coord[2]*v3[pp][1]) +
											s[2]*(atom[ll].coord[0]*v1[pp][2]+atom[ll].coord[1]*v2[pp][2]+atom[ll].coord[2]*v3[pp][2]) ) ) ) ;

											FFF_cdummy *= atom[ll].mult ;

											/* compute structure amplitude F_k for X-ray with the atomic scattering factors */
											if ( Xn >= 0 ) { FFF_X[ii][jj][ind_r3[pp]] += FFF_cdummy * f0table[atom[ll].f0ind][ii] * ( (dcmplx) sc_rho ) ; } /* [nm] */

											/* the same for neutrons with the coherent neutron scattering lengths
											   don't use any minus sign convention for b_coh, use just b_coh,
											   the sld will also not include later minus signs
											*/
											if ( Xn <= 0 ) { FFF_n[ii][jj][ind_r3[pp]] += FFF_cdummy * atom[ll].coh_sc_length * ( (dcmplx) sc_scl ) ; } /* [nm] */

											if ( test_XnEquiv_flag ) { FFF_n[ii][jj][ind_r3[pp]] = FFF_X[ii][jj][ind_r3[pp]] ; }
										}
									}

									piQv1 = M_PI * cs * eQvj[jj][ind_r3[pp]][0] ;
									piQv2 = M_PI * cs * eQvj[jj][ind_r3[pp]][1] ;
									piQv3 = M_PI * cs * eQvj[jj][ind_r3[pp]][2] ;

									piQwt1 = okf1[pp] * piQv1 ;
									piQwt2 = okf2[pp] * piQv2 ;

									piQw1 = kf1[pp] * piQv1 ;
									piQw2 = kf2[pp] * piQv2 ;

									exp_N1Qv1 = exp ( dcmplx( 0.0, 2.0 * piQv1 * (double)n1[pp]) ) - 1.0 ;
									exp_N2Qv2 = exp ( dcmplx( 0.0, 2.0 * piQv2 * (double)n2[pp]) ) - 1.0 ;

									exp_Qv1 = exp ( dcmplx( 0.0, 2.0 * piQv1 ) ) - 1.0 ;
									exp_Qv2 = exp ( dcmplx( 0.0, 2.0 * piQv2 ) ) - 1.0 ;
									exp_Qv3 = exp ( dcmplx( 0.0, 2.0 * piQv3 ) ) - 1.0 ;


									/* for GA */
									GA_N1_N2 = dcmplx( 1.0, 0.0) ;

									if ( fabs( real( exp_Qv1 * conj(exp_Qv1) ) ) > 1e-10) { GA_N1_N2 *= ( exp_N1Qv1 / exp_Qv1 ) ; }
									else { GA_N1_N2 *= dcmplx( ( (double)n1[pp]), 0.0) ; }

									if ( fabs( real( exp_Qv2 * conj(exp_Qv2) ) ) > 1e-10) { GA_N1_N2 *= ( exp_N2Qv2 / exp_Qv2 ) ; }
									else { GA_N1_N2 *= dcmplx( ( (double)n2[pp]), 0.0) ; }


									/* for background and dispersion medium S_D */

									/* w_j = kf_j * v_j */
									exp_N1Qw1 = exp ( dcmplx( 0.0, 2.0 * piQw1 * (double)n1[pp] ) ) - 1.0 ;
									exp_N2Qw2 = exp ( dcmplx( 0.0, 2.0 * piQw2 * (double)n2[pp] ) ) - 1.0 ;
									/* wt_j = okf_j * v_j */
									exp_N1Qwt1 = exp ( dcmplx( 0.0, 2.0 * piQwt1 * (double)n1[pp] ) ) - 1.0 ;
									exp_N2Qwt2 = exp ( dcmplx( 0.0, 2.0 * piQwt2 * (double)n2[pp] ) ) - 1.0 ;


									if (fabs(piQwt1) > 1e-10) { cdummy21 = exp_N1Qwt1 / ( dcmplx( 0.0, 2.0 * piQwt1 ) ) ; }
									else { cdummy21 = dcmplx( ( (double)n1[pp]), 0.0) ; }

									if (fabs(piQwt2) > 1e-10) { cdummy31 = exp_N2Qwt2 / ( dcmplx( 0.0, 2.0 * piQwt2) ) ; }
									else { cdummy31 = dcmplx( ( (double)n2[pp]), 0.0) ; }


									if (fabs(piQw1) > 1e-10) { cdummy22 = exp_N1Qw1 / ( dcmplx( 0.0, 2.0 * piQw1) ) ; }
									else { cdummy22 = dcmplx( ( (double)n1[pp]), 0.0) ; }

									if (fabs(piQw2) > 1e-10) { cdummy32 = exp_N2Qw2 / ( dcmplx( 0.0, 2.0 * piQw2) ) ; }
									else { cdummy32 = dcmplx( ( (double)n2[pp]), 0.0) ; }


									if (fabs(piQv1) > 1e-10) { cdummy23 = exp_N1Qv1 / ( dcmplx( 0.0, 2.0 * piQv1) ) ; }
									else { cdummy23 = dcmplx( ( (double)n1[pp]), 0.0) ; }

									if (fabs(piQv2) > 1e-10) { cdummy33 = exp_N2Qv2 / ( dcmplx( 0.0, 2.0 * piQv2) ) ; }
									else { cdummy33 = dcmplx( ( (double)n2[pp]), 0.0) ; }

									/* compute scattering contribution for all single nanoparticles with n3 ~ m = 1 ... par->nsp */
									/* save in dE[mm] and add |dE[mm]|^2 to S[ii][mm] */
									for ( mm=0; mm<par->nsp; ++mm)
									{
										piQwt3 = okf3[mm] * piQv3 ;
										piQw3 = kf3[mm] * piQv3 ;
		
										/* lattice amplitude GA */
										exp_N3Qv3[0] = exp ( dcmplx( 0.0, 2.0 * piQv3 * (double)(mm+1) ) ) - 1.0 ;
										/* w_j = kf_j * v_j */
										exp_N3Qw3[0] = exp ( dcmplx( 0.0, 2.0 * piQw3 * (double)(mm+1) ) ) - 1.0 ;
										/* wt_j = okf_j * v_j */
										exp_N3Qwt3[0] = exp ( dcmplx( 0.0, 2.0 * piQwt3 * (double)(mm+1) ) ) - 1.0 ;
			
										/* lattice amplitude GA */
										/* GA=(exp_N1Qv1/exp_Qv1)*(exp_N2Qv2/exp_Qv2)*(exp_N3Qv3/exp_Qv3) */
										if ( fabs( real( exp_Qv3 * conj(exp_Qv3) ) ) > 1e-10) { GA[mm][pp] = GA_N1_N2 * ( exp_N3Qv3[0] / exp_Qv3 ) ; }
										else { GA[mm][pp] = GA_N1_N2 * dcmplx( (double)(mm+1), 0.0) ; }
			
										/* background Bg */
										if (fabs(piQwt3) > 1e-10) { cdummy41[0] = exp_N3Qwt3[0] / ( dcmplx( 0.0, 2.0 * piQwt3) ) ; }
										else { cdummy41[0] = dcmplx( (double)(mm+1), 0.0) ; }

										if (fabs(piQw3) > 1e-10) { cdummy42[0] = exp_N3Qw3[0] / ( dcmplx( 0.0, 2.0 * piQw3) ) ; }
										else { cdummy42[0] = dcmplx( (double)(mm+1), 0.0) ; }

										if (fabs(piQv3) > 1e-10) { cdummy43[0] = exp_N3Qv3[0] / ( dcmplx( 0.0, 2.0 * piQv3) ) ; }
										else { cdummy43[0] = dcmplx( (double)(mm+1), 0.0) ; }

										cdummy11 = cdummy21 * cdummy31 * cdummy41[0] ;
										cdummy12 = cdummy22 * cdummy32 * cdummy42[0] ;
										cdummy13 = cdummy23 * cdummy33 * cdummy43[0] ;

										P_osl[mm][pp] = cdummy11 * ((dcmplx) V_osl[mm][pp]) ; /* [nm^3] */
										P_isl[mm][pp] = cdummy12 * ((dcmplx) V_isl[mm][pp]) ;
										P[mm][pp] = cdummy13 * ((dcmplx) V[mm][pp]) ;

										if ( Xn >= 0 )
										{
											for ( qq=0; qq<rho_multi.n; ++qq)
											{
												Bg_X = dcmplx(0.0,0.0) ;  /* [nm] */
												Bg_X += ((dcmplx) drerho_dm_osl[qq] ) * P_osl[mm][pp] ;
												Bg_X += ((dcmplx) drerho_osl_isl[qq] ) * P_isl[mm][pp] ;
												Bg_X += ((dcmplx) drerho_isl[qq] ) * P[mm][pp] ;

												dE_X[mm][pp][qq] = FFF_X[ii][jj][ind_r3[pp]] * GA[mm][pp] - Bg_X ;  /* [nm] */

												dS_X_parallel_thread[mm][pp][qq][ii] += cabs2( dE_X[mm][pp][qq] ) * weight_factor ; /* [nm^2] */
											}
										}

										if ( Xn <= 0 )
										{
											for ( qq=0; qq<sld_multi.n; ++qq)
											{
												Bg_n = dcmplx(0.0,0.0) ; /* [nm] */
												Bg_n += ((dcmplx) dsld_dm_osl[qq] ) * P_osl[mm][pp] ;
												Bg_n += ((dcmplx) dsld_osl_isl[qq] ) * P_isl[mm][pp] ;
												Bg_n += ((dcmplx) dsld_isl[qq] ) * P[mm][pp] ;

												dE_n[mm][pp][qq] = FFF_n[ii][jj][ind_r3[pp]] * GA[mm][pp] - Bg_n ; /* [nm] */

												dS_n_parallel_thread[mm][pp][qq][ii] += cabs2( dE_n[mm][pp][qq] ) * weight_factor ; /* [nm^2] */
											}
										}


										/* add-up the terms for Yc_X/n[mm][0...9][ii] and for Yi_n[mm][0...3] */
										/* !!! it is summed automatically over all pp !!! */
										/* note that Yc_X[mm][qq][ii] ~ Yc_n[mm][qq][ii] for qq=1...6, because F isn't involved into it */
										if ( Xn >= 0 )
										{
											Yc_X_parallel_thread[mm][0][ii] += cabs2( FFF_X[ii][jj][ind_r3[pp]] * GA[mm][pp] ) * weight_factor ; /* [nm^2] */

											Yc_X_parallel_thread[mm][1][ii] += cabs2( P_osl[mm][pp] ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
											Yc_X_parallel_thread[mm][2][ii] += cabs2( P_isl[mm][pp] ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
											Yc_X_parallel_thread[mm][3][ii] += cabs2( P[mm][pp] ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */

											Yc_X_parallel_thread[mm][4][ii] += 2.0 * real ( P_osl[mm][pp] * conj( P_isl[mm][pp] ) ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
											Yc_X_parallel_thread[mm][5][ii] += 2.0 * real ( P_osl[mm][pp] * conj( P[mm][pp] ) ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
											Yc_X_parallel_thread[mm][6][ii] += 2.0 * real ( P_isl[mm][pp] * conj( P[mm][pp] ) ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */

											Yc_X_parallel_thread[mm][7][ii] += -2.0 * real ( FFF_X[ii][jj][ind_r3[pp]] * GA[mm][pp] * conj( P_osl[mm][pp] ) ) * sc_rho * weight_factor ; /* [nm^5] */
											Yc_X_parallel_thread[mm][8][ii] += -2.0 * real ( FFF_X[ii][jj][ind_r3[pp]] * GA[mm][pp] * conj( P_isl[mm][pp] ) ) * sc_rho * weight_factor ; /* [nm^5] */
											Yc_X_parallel_thread[mm][9][ii] += -2.0 * real ( FFF_X[ii][jj][ind_r3[pp]] * GA[mm][pp] * conj( P[mm][pp] ) ) * sc_rho * weight_factor ; /* [nm^5] */
										}
										if ( Xn <= 0 )
										{
											Yc_n_parallel_thread[mm][0][ii] += cabs2( FFF_n[ii][jj][ind_r3[pp]] * GA[mm][pp] ) * weight_factor ; /* [nm^2] */

											Yc_n_parallel_thread[mm][1][ii] += cabs2( P_osl[mm][pp] ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
											Yc_n_parallel_thread[mm][2][ii] += cabs2( P_isl[mm][pp] ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
											Yc_n_parallel_thread[mm][3][ii] += cabs2( P[mm][pp] ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
												
											Yc_n_parallel_thread[mm][4][ii] += 2.0 * real ( P_osl[mm][pp] * conj( P_isl[mm][pp] ) ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
											Yc_n_parallel_thread[mm][5][ii] += 2.0 * real ( P_osl[mm][pp] * conj( P[mm][pp] ) ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
											Yc_n_parallel_thread[mm][6][ii] += 2.0 * real ( P_isl[mm][pp] * conj( P[mm][pp] ) ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */

											Yc_n_parallel_thread[mm][7][ii] += -2.0 * real ( FFF_n[ii][jj][ind_r3[pp]] * GA[mm][pp] * conj( P_osl[mm][pp] ) ) * sc_sld * weight_factor ; /* [nm^4] */
											Yc_n_parallel_thread[mm][8][ii] += -2.0 * real ( FFF_n[ii][jj][ind_r3[pp]] * GA[mm][pp] * conj( P_isl[mm][pp] ) ) * sc_sld * weight_factor ; /* [nm^4] */
											Yc_n_parallel_thread[mm][9][ii] += -2.0 * real ( FFF_n[ii][jj][ind_r3[pp]] * GA[mm][pp] * conj( P[mm][pp] ) ) * sc_sld * weight_factor ; /* [nm^4] */
										}
									} /* mm = 0 ... par->nsp-1 */
								} /* pp = 0 ... par->nms-1 */
							} /* if av_mode == 0,2,3 or 1 */


							/* compute scattering contribution for the randomly chosen stacks */


							/* Actually there's no need for exp_RQ[0], computation might be omitted */
							/* exp_RQ[0] = exp (dcmplx( 0.0, 2.0 * M_PI * sp (R[0], s))) ; */
							exp_RQ[0] = dcmplx( 1.0, 0.0)  ;
							/* now compute the necessary exponentials for R[1] ... R[par->nms-1] */
							for ( mm=par->nsp; mm<nspsp; ++mm)
							{
								/* idummy = 1,2,3,...,(par->nms - 1) for the 2nd,3rd,4th.,...par->nms th particle in a stack,
								Note that for idummy = 0, i.e. the 1st particle, exp_RQ = 1. */
								idummy = mm - par->nsp + 1 ;
								exp_RQ[idummy] = exp( dcmplx( 0.0, 2.0 * M_PI * sp( R[idummy], s) ) ) ;
							}

							if ( Xn >= 0 )
							{
								for ( qq=0; qq<rho_multi.n; ++qq)
								{
									/* first particle in each the stack, where R[0]={0,0,0}, skip the first exponential exp_RQ[0] = 1.0 */
									/* dEstack_X/n = exp( dcmplx( 0.0, 2.0 * M_PI * sp (R[0], s) ) ) * dE_X/n[n3[0]-1][0] =  dE_X/n[n3[0]-1][0] */
									dEstack_X = dE_X[n3[0]-1][0][qq] ;
									/* Add up the contributions of all other particles in the stack,
									the whole procedure is done iteratively, i.e. the data for the 2-stack
									is used by all thicker >2-stacks and so on.
									Just to remind dE[n3[idummy]-1][idummy] is the complex scattered amplitude for the
									idummy-th particle in the stack that has random sizes of n1[idummy],n2[idummy],n3[idummy]
									and if av_mode=1 also a random r3[idummy] azimuthal rotation around Gs */
									for ( mm=par->nsp; mm<nspsp; ++mm)
									{
										/* idummy = 1,2,3,...,(par->nms - 1) for the 2nd,3rd,4th.,...par->nms th particle in a stack,
										Note that idummy = 0, i.e. the 1st particle has been already done above. */
										idummy = mm - par->nsp + 1 ;
										dEstack_X += exp_RQ[idummy] * dE_X[n3[idummy]-1][idummy][qq] ;

										dS_X_parallel_thread[mm][0][qq][ii] += fabs( real( dEstack_X * conj(dEstack_X) ) ) * weight_factor ;
									}
								}
							}
							if ( Xn <= 0 )
							{
								for ( qq=0; qq<sld_multi.n; ++qq)
								{
									dEstack_n = dE_n[n3[0]-1][0][qq] ;

									for ( mm=par->nsp; mm<nspsp; ++mm)
									{
										idummy = mm - par->nsp + 1 ;
		
										dEstack_n += exp_RQ[idummy] * dE_n[n3[idummy]-1][idummy][qq] ;

										dS_n_parallel_thread[mm][0][qq][ii] += fabs( real( dEstack_n * conj(dEstack_n) ) ) * weight_factor ;
									}
								}
							}


							/* add-up the terms for Yc_X/n[mm][0...9][ii] and for Yi_n[mm][0...3] */
							/* note that Yc_X[mm][qq][ii] ~ Yc_n[mm][qq][ii] for qq=1...6, because F isn't involved into it */
							if ( Xn >= 0 )
							{
								/* clear dYc_X_stack */
								for ( mm=0; mm<10; ++mm)
								{
									dYc_X_stack[mm] = 0.0 ;
								}

								for ( pp=0; pp<par->nms; ++pp)
								{
									/* pp = 0,1,2,3,...,(par->nms - 1) for the 1st,2nd,3rd,4th.,...par->nms th particle in a stack */

									/* respect correct index for phi valid for all av_mode */
									if ( ( par->av_mode == 0 ) || ( par->av_mode == 2 ) || ( par->av_mode == 3 ) ) { k_phi = kk ; }
									else if ( par->av_mode == 1 ) { k_phi = ind_r3[pp] ; }

									/* scattering contributions without interparticle interferences */
									dYc_X_stack[0] += cabs2( FFF_X[ii][jj][k_phi] * GA[n3[pp]-1][pp] ) * weight_factor ; /* [nm^2] */
									dYc_X_stack[1] += cabs2( P_osl[n3[pp]-1][pp] ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
									dYc_X_stack[2] += cabs2( P_isl[n3[pp]-1][pp] ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
									dYc_X_stack[3] += cabs2( P[n3[pp]-1][pp] ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */

									dYc_X_stack[4] += 2.0 * real ( P_osl[n3[pp]-1][pp] * conj( P_isl[n3[pp]-1][pp] ) ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
									dYc_X_stack[5] += 2.0 * real ( P_osl[n3[pp]-1][pp] * conj( P[n3[pp]-1][pp] ) ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
									dYc_X_stack[6] += 2.0 * real ( P_isl[n3[pp]-1][pp] * conj( P[n3[pp]-1][pp] ) ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */

									dYc_X_stack[7] += -2.0 * real ( FFF_X[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P_osl[n3[pp]-1][pp] ) ) * sc_rho * weight_factor ; /* [nm^5] */
									dYc_X_stack[8] += -2.0 * real ( FFF_X[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P_isl[n3[pp]-1][pp] ) ) * sc_rho * weight_factor ; /* [nm^5] */ 
									dYc_X_stack[9] += -2.0 * real ( FFF_X[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P[n3[pp]-1][pp] ) ) * sc_rho * weight_factor ; /* [nm^5] */

									/* scattering contributions from interparticle interferences, applies only from the 2nd particle on (pp>0) */
									for ( mm=0; mm<pp; ++mm)
									{
										exp_DRQ = exp_RQ[pp] * conj( exp_RQ[mm] ) ;
										dYc_X_stack[0] += 2.0 * cabs2( FFF_X[ii][jj][k_phi] ) * real ( GA[n3[pp]-1][pp] * conj( GA[n3[mm]-1][mm] ) * exp_DRQ ) * weight_factor ; /* [nm^2] */

										dYc_X_stack[1] += 2.0 * real ( P_osl[n3[pp]-1][pp] * conj( P_osl[n3[mm]-1][mm] ) * exp_DRQ ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
										dYc_X_stack[2] += 2.0 * real ( P_isl[n3[pp]-1][pp] * conj( P_isl[n3[mm]-1][mm] ) * exp_DRQ ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
										dYc_X_stack[3] += 2.0 * real ( P[n3[pp]-1][pp] * conj( P[n3[mm]-1][mm] ) * exp_DRQ ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */

										dYc_X_stack[4] += 2.0 * real ( ( P_osl[n3[pp]-1][pp] * conj( P_isl[n3[mm]-1][mm] ) + P_isl[n3[pp]-1][pp] * conj( P_osl[n3[mm]-1][mm] ) ) * exp_DRQ ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
										dYc_X_stack[5] += 2.0 * real ( ( P_osl[n3[pp]-1][pp] * conj( P[n3[mm]-1][mm] ) + P[n3[pp]-1][pp] * conj( P_osl[n3[mm]-1][mm] ) ) * exp_DRQ ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */
										dYc_X_stack[6] += 2.0 * real ( ( P_isl[n3[pp]-1][pp] * conj( P[n3[mm]-1][mm] ) + P[n3[pp]-1][pp] * conj( P_isl[n3[mm]-1][mm] ) ) * exp_DRQ ) * pow( sc_rho, 2) * weight_factor ; /* [nm^8] */

										dYc_X_stack[7] += -2.0 * real ( ( FFF_X[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P_osl[n3[mm]-1][mm] ) + P_osl[n3[pp]-1][pp] * conj( FFF_X[ii][jj][k_phi] * GA[n3[mm]-1][mm] ) ) * exp_DRQ ) * sc_rho * weight_factor ; /* [nm^5] */
										dYc_X_stack[8] += -2.0 * real ( ( FFF_X[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P_isl[n3[mm]-1][mm] ) + P_isl[n3[pp]-1][pp] * conj( FFF_X[ii][jj][k_phi] * GA[n3[mm]-1][mm] ) ) * exp_DRQ ) * sc_rho * weight_factor ; /* [nm^5] */
										dYc_X_stack[9] += -2.0 * real ( ( FFF_X[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P[n3[mm]-1][mm] ) + P[n3[pp]-1][pp] * conj( FFF_X[ii][jj][k_phi] * GA[n3[mm]-1][mm] ) ) * exp_DRQ ) * sc_rho * weight_factor ; /* [nm^5] */
									}

									/* add values from the current (p+1)-stack to Yc_X/n */
									if ( pp > 0 )
									{
										for ( mm=0; mm<10; ++mm)
										{
											Yc_X_parallel_thread[pp-1+par->nsp][mm][ii] += dYc_X_stack[mm] ;
										}
									}
								}
							}

							if ( Xn <= 0 )
							{
								/* clear dYc_n_stack */
								for ( mm=0; mm<10; ++mm)
								{
									dYc_n_stack[mm] = 0.0 ;
								}

								for ( pp=0; pp<par->nms; ++pp)
								{
									/* pp = 0,1,2,3,...,(par->nms - 1) for the 1st,2nd,3rd,4th.,...par->nms th particle in a stack */

									/* respect correct index for phi valid for all av_mode */
									if ( ( par->av_mode == 0 ) || ( par->av_mode == 2 ) || ( par->av_mode == 3 ) ) { k_phi = kk ; }
									else if ( par->av_mode == 1 ) { k_phi = ind_r3[pp] ; }

									/* scattering contributions without interparticle interferences */
									dYc_n_stack[0] += cabs2( FFF_n[ii][jj][k_phi] * GA[n3[pp]-1][pp] ) * weight_factor ; /* [nm^2] */
									dYc_n_stack[1] += cabs2( P_osl[n3[pp]-1][pp] ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
									dYc_n_stack[2] += cabs2( P_isl[n3[pp]-1][pp] ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
									dYc_n_stack[3] += cabs2( P[n3[pp]-1][pp] ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */

									dYc_n_stack[4] += 2.0 * real ( P_osl[n3[pp]-1][pp] * conj( P_isl[n3[pp]-1][pp] ) ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
									dYc_n_stack[5] += 2.0 * real ( P_osl[n3[pp]-1][pp] * conj( P[n3[pp]-1][pp] ) ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
									dYc_n_stack[6] += 2.0 * real ( P_isl[n3[pp]-1][pp] * conj( P[n3[pp]-1][pp] ) ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */

									dYc_n_stack[7] += -2.0 * real ( FFF_n[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P_osl[n3[pp]-1][pp] ) ) * sc_sld * weight_factor ; /* [nm^4] */
									dYc_n_stack[8] += -2.0 * real ( FFF_n[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P_isl[n3[pp]-1][pp] ) ) * sc_sld * weight_factor ; /* [nm^4] */ 
									dYc_n_stack[9] += -2.0 * real ( FFF_n[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P[n3[pp]-1][pp] ) ) * sc_sld * weight_factor ; /* [nm^4] */

									/* scattering contributions from interparticle interferences, applies only from the 2nd particle on (p>0) */
									for ( mm=0; mm<pp; ++mm)
									{
										exp_DRQ = exp_RQ[pp] * conj( exp_RQ[mm] ) ;
										dYc_n_stack[0] += 2.0 * cabs2( FFF_n[ii][jj][k_phi] ) * real ( GA[n3[pp]-1][pp] * conj( GA[n3[mm]-1][mm] ) * exp_DRQ ) * weight_factor ; /* [nm^2] */

										dYc_n_stack[1] += 2.0 * real ( P_osl[n3[pp]-1][pp] * conj( P_osl[n3[mm]-1][mm] ) * exp_DRQ ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
										dYc_n_stack[2] += 2.0 * real ( P_isl[n3[pp]-1][pp] * conj( P_isl[n3[mm]-1][mm] ) * exp_DRQ ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
										dYc_n_stack[3] += 2.0 * real ( P[n3[pp]-1][pp] * conj( P[n3[mm]-1][mm] ) * exp_DRQ ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */

										dYc_n_stack[4] += 2.0 * real ( ( P_osl[n3[pp]-1][pp] * conj( P_isl[n3[mm]-1][mm] ) + P_isl[n3[pp]-1][pp] * conj( P_osl[n3[mm]-1][mm] ) ) * exp_DRQ ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
										dYc_n_stack[5] += 2.0 * real ( ( P_osl[n3[pp]-1][pp] * conj( P[n3[mm]-1][mm] ) + P[n3[pp]-1][pp] * conj( P_osl[n3[mm]-1][mm] ) ) * exp_DRQ ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */
										dYc_n_stack[6] += 2.0 * real ( ( P_isl[n3[pp]-1][pp] * conj( P[n3[mm]-1][mm] ) + P[n3[pp]-1][pp] * conj( P_isl[n3[mm]-1][mm] ) ) * exp_DRQ ) * pow( sc_sld, 2) * weight_factor ; /* [nm^6] */

										dYc_n_stack[7] += -2.0 * real ( ( FFF_n[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P_osl[n3[mm]-1][mm] ) + P_osl[n3[pp]-1][pp] * conj( FFF_n[ii][jj][k_phi] * GA[n3[mm]-1][mm] ) ) * exp_DRQ ) * sc_sld * weight_factor ; /* [nm^4] */
										dYc_n_stack[8] += -2.0 * real ( ( FFF_n[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P_isl[n3[mm]-1][mm] ) + P_isl[n3[pp]-1][pp] * conj( FFF_n[ii][jj][k_phi] * GA[n3[mm]-1][mm] ) ) * exp_DRQ ) * sc_sld * weight_factor ; /* [nm^4] */
										dYc_n_stack[9] += -2.0 * real ( ( FFF_n[ii][jj][k_phi] * GA[n3[pp]-1][pp] * conj( P[n3[mm]-1][mm] ) + P[n3[pp]-1][pp] * conj( FFF_n[ii][jj][k_phi] * GA[n3[mm]-1][mm] ) ) * exp_DRQ ) * sc_sld * weight_factor ; /* [nm^4] */
									}

									/* add values from the current (p+1)-stack to Yc_X/n */
									if ( pp > 0 ) 
									{
										for ( mm=0; mm<10; ++mm)
										{
											Yc_n_parallel_thread[pp-1+par->nsp][mm][ii] += dYc_n_stack[mm] ;
										}
									}
								}
							}
						/*end FOURTH loop over the s range */
						}
					/* end THIRD loop over polar angle */
					}
				/* end SECOND.2 (parallelized) loop over the azimuthal angle range */
				}

				/******************************/
				/* Finalizing parallel region */
				/******************************/

				/* write data from different threads to dS_X/n and Y_X/n */
				if ( Xn >= 0 )
				{
					for ( mm=0; mm<nspsp; ++mm)
					{
						if ( mm<par->nsp )
						{
							for ( pp=0; pp<par->nms; ++pp)
							{
								for ( qq=0; qq<rho_multi.n; ++qq)
								{
									for ( ii=0; ii<par->np; ++ii)
									{
										#pragma omp critical
										{
											dS_X[mm][pp][qq][ii] += dS_X_parallel_thread[mm][pp][qq][ii] ;
										}
									}
								}
							}
						}
						else
						{
							for ( qq=0; qq<rho_multi.n; ++qq)
							{
								for ( ii=0; ii<par->np; ++ii)
								{
									#pragma omp critical
									{
										dS_X[mm][0][qq][ii] += dS_X_parallel_thread[mm][0][qq][ii] ;
									}
								}
							}
						}
					}
				}
				if ( Xn <= 0 )
				{
					for ( mm=0; mm<nspsp; ++mm)
					{
						if ( mm<par->nsp )
						{
							for ( pp=0; pp<par->nms; ++pp)
							{
								for ( qq=0; qq<sld_multi.n; ++qq)
								{
									for ( ii=0; ii<par->np; ++ii)
									{
										#pragma omp critical
										{
											dS_n[mm][pp][qq][ii] += dS_n_parallel_thread[mm][pp][qq][ii] ;
										}
									}
								}
							}
						}
						else
						{
							for ( qq=0; qq<sld_multi.n; ++qq)
							{
								for ( ii=0; ii<par->np; ++ii)
								{
									#pragma omp critical
									{
										dS_n[mm][0][qq][ii] += dS_n_parallel_thread[mm][0][qq][ii] ;
									}
								}
							}
						}
					}
				}

				for ( mm=0; mm<nspsp; ++mm) 
				{
					for ( pp=0; pp<10; ++pp)
					{ 
						for ( ii=0; ii<par->np; ++ii)
						{
							if ( Xn >= 0 )
							{
								#pragma omp critical
								{
									Yc_X[mm][pp][ii] += Yc_X_parallel_thread[mm][pp][ii] ;
								}
							}
							if ( Xn <= 0 )
							{
								#pragma omp critical
								{
									Yc_n[mm][pp][ii] += Yc_n_parallel_thread[mm][pp][ii] ;
								}
							}
						}
					}
				}



				/* free dynamic arrays allocated within the parallel region
				   GA, P, P_isl, P_osl, dE_X, dE_n, exp_N3Qv3, exp_N3Qw3, exp_N3Qwt3, cdummy41, cdummy42, cdummy43, exp_RQ,
				   dS_X_parallel_thread, dS_n_parallel_thread, Yc_X_parallel_thread, Yc_n_parallel_thread
				*/
				for ( mm=0; mm<par->nsp; ++mm)
				{
					free(GA[mm]) ;
					free(P[mm]) ;
					free(P_isl[mm]) ;
					free(P_osl[mm]) ;
				}
				free(GA) ;
				free(P) ;
				free(P_isl) ;
				free(P_osl) ;

				/* free dE_X/n which are a double complex matrices of size ( par->nsp x par->nms x rho/sld_multi.n ) */
				if ( Xn >= 0 )
				{
					for ( mm=0; mm<par->nsp; ++mm) 
					{
						for ( pp=0; pp<par->nms; ++pp) { free(dE_X[mm][pp]) ; }
						free(dE_X[mm]) ;
					}
					free(dE_X);
				}
				if ( Xn <= 0 )
				{
					for ( mm=0; mm<par->nsp; ++mm) 
					{
						for ( pp=0; pp<par->nms; ++pp) { free(dE_n[mm][pp]) ; }
						free(dE_n[mm]) ;
					}
					free(dE_n);
				}

				free(exp_N3Qv3) ;
				free(exp_N3Qw3) ;
				free(exp_N3Qwt3) ;
				free(cdummy41) ;
				free(cdummy42) ;
				free(cdummy43) ;

				free(exp_RQ) ;

				/* free dS_X/n_thread_parallel */
				if ( Xn >= 0 )
				{
					for ( mm=0; mm<nspsp; ++mm)
					{
						if ( mm<par->nsp )
						{
							for ( pp=0; pp<par->nms; ++pp)
							{
								for ( qq=0; qq<rho_multi.n; ++qq) { free(dS_X_parallel_thread[mm][pp][qq]) ; }
								free(dS_X_parallel_thread[mm][pp]) ;
							}
						}
						else
						{
							for ( qq=0; qq<rho_multi.n; ++qq) { free(dS_X_parallel_thread[mm][0][qq]) ; }
							free(dS_X_parallel_thread[mm][0]) ;
						}
						free(dS_X_parallel_thread[mm]) ;
					}
					free(dS_X_parallel_thread) ;
				}
				if ( Xn <= 0 )
				{
					for ( mm=0; mm<nspsp; ++mm)
					{
						if ( mm<par->nsp )
						{
							for ( pp=0; pp<par->nms; ++pp)
							{
								for ( qq=0; qq<sld_multi.n; ++qq) { free(dS_n_parallel_thread[mm][pp][qq]) ; }
								free(dS_n_parallel_thread[mm][pp]) ;
							}
						}
						else
						{
							for ( qq=0; qq<sld_multi.n; ++qq) { free(dS_n_parallel_thread[mm][0][qq]) ; }
							free(dS_n_parallel_thread[mm][0]) ;
						}
						free(dS_n_parallel_thread[mm]) ;
					}
					free(dS_n_parallel_thread) ;
				}

				/* free Yc_X/n_parallel_thread */
				if ( Xn >= 0 )
				{
					for ( mm=0; mm<nspsp; ++mm) 
					{
						for ( pp=0; pp<10; ++pp) { free(Yc_X_parallel_thread[mm][pp]) ; }
						free(Yc_X_parallel_thread[mm]) ;
					}
					free(Yc_X_parallel_thread);
				}
				if ( Xn <= 0 )
				{
					for ( mm=0; mm<nspsp; ++mm) 
					{
						for ( pp=0; pp<10; ++pp) { free(Yc_n_parallel_thread[mm][pp]) ; }
						free(Yc_n_parallel_thread[mm]) ;
					}
					free(Yc_n_parallel_thread);
				}
			} /* end pragma parallel */


			/* write now data from dS_X/n[i][m][p][q] to S_X/n[i][m][q]
			   and similarly with dB_n[m][p][q] to B_n[m][q] */
			/* for normalizations see comments further below, essentially par->nms corrections for single particles takes place via V_irr ! */
			if  ( Xn >= 0 )
			{
				for ( unsigned int m=0; m<nspsp; ++m)
				{
					if ( m < par->nsp )
					{
						/* single particles, sum over index p */
						for ( unsigned int p=0; p<par->nms; ++p)
						{
							for ( unsigned int q=0; q<rho_multi.n; ++q)
							{
								for ( unsigned int i=0; i<par->np; ++i )
									S_X[m][q][i] += dS_X[m][p][q][i] ;
							}
						}
					}
					else
					{
						/* stacks, no p summation */
						for ( unsigned int q=0; q<rho_multi.n; ++q)
						{
							for ( unsigned int i=0; i<par->np; ++i )
								S_X[m][q][i] += dS_X[m][0][q][i] ;
						}
					}
				}
			}
			if  ( Xn <= 0 )
			{
				for ( unsigned int m=0; m<nspsp; ++m)
				{
					if ( m < par->nsp )
					{
						/* single particles, sum over index p */
						for ( unsigned int p=0; p<par->nms; ++p)
						{
							for ( unsigned int q=0; q<sld_multi.n; ++q)
							{
								for ( unsigned int i=0; i<par->np; ++i )
									S_n[m][q][i] += dS_n[m][p][q][i] ;
							}
						}
					}
					else
					{
						/* stacks, no p summation */
						for ( unsigned int q=0; q<sld_multi.n; ++q)
						{
							for ( unsigned int i=0; i<par->np; ++i )
								S_n[m][q][i] += dS_n[m][0][q][i] ;
						}
					}
				}
			}


			/* if FNoMem_flag was set write FFF after the first c_rep=0 is done */
			if ( Fw_flag && FNoMem_flag && c_rep == 0 )
			{
				if ( log_flag ) { fprintf( logfile, "write_structure_amplitude()\n") ; fflush (logfile) ; }
				if ( write_structure_amplitude() )
				{
					if ( log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }
				}
				else
				{
					if ( log_flag ) { fprintf( logfile, "failed\n\n") ; fflush (logfile) ; }
				}
			}


			if ( orav_flag && FNoMem_flag && c_rep == 0 )
			{
				if ( Xn >= 0 )
				{
					if ( log_flag ) { fprintf( logfile, "compute_orientational_averaged_structure_factor( &FFF_X[0], 'X') for X-ray\n") ; fflush (logfile) ; }
					compute_orientational_averaged_structure_factor( &FFF_X[0], 'X') ;
					if ( log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }

					if ( log_flag ) { fprintf( logfile, "compute_orientational_averaged_structure_amplitude( &FFF_X[0], 'X') for X-ray\n") ; fflush (logfile) ; }
					compute_orientational_averaged_structure_amplitude( &FFF_X[0], 'X') ;
					if ( log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }
				}

				if ( Xn <= 0 )
				{
					if ( log_flag ) { fprintf( logfile, "compute_orientational_averaged_structure_factor( &FFF_n[0], 'n') for neutrons\n") ; fflush (logfile) ; }
					compute_orientational_averaged_structure_factor( &FFF_n[0], 'n') ;
					if ( log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }

					if ( log_flag ) { fprintf( logfile, "compute_orientational_averaged_structure_amplitude( &FFF_n[0], 'n') for neutrons\n") ; fflush (logfile) ; }
					compute_orientational_averaged_structure_amplitude( &FFF_n[0], 'n') ;
					if ( log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }
				}
			}


			/* Normalization for (d)S_X/n for:
				single particles output (each c_rep)
				intermediate step output (determined by the -z flag)
				final output (if c_rep = par->nr)

			   Normalization is done in three ways:
				-MULTIPLICITY accounts for number of Powder Average steps, orientational average over all orientations
				-for single particles there you have additionally the par->nms multiplicity, this is already respected in V_irr,
				 the particles volumes to that (d)S_X/n is normalized
				-SUM_WEIGHT_CORR is a correction factor, 1 when no sin-weighting is used, 2/pi, when sin-weighting was applied

			  Incoherent contribution in case of neutron scattering is dealt similarly, normalization takes place by V_irr,
			  but MULT_INC = 1 cause there is no need of powder average ! Furthermore no need of weighting factors.
			*/



			/* write dS, dB_n to files */
			if (singlefilesoutput_flag)
			{
				for ( unsigned int m=0; m<nspsp; ++m)
				{
					if ( m<par->nsp ) { idummy = par->nms ; }
					else { idummy = 1 ; }

					MULTIPLICITY[m] = 1 ;
					MULTIPLICITY[m] *= nav_ac ;
					MULTIPLICITY[m] *= nav_pc ;
					/* note that each of the dS_X/n[m][k][:] files will be written individually, that means there is no summation over k for single particles */

					sprintf ( sdummy, "%s\n%ssinglestep output at (c_rep+1) = %d\n", presp, presp, c_rep+1) ;
					sprintf ( sdummy2, "%sMULTIPLICITY = %d ( nav_ac = %d, nav_pc = %d)\n", presp, MULTIPLICITY[m], nav_ac, nav_pc) ;
					strcat( sdummy, sdummy2) ;

					for (int p=0; p<idummy; ++p)
					{
						if  ( Xn >= 0 )
						{
							/* single particles / stacks */
							if ( m < par->nsp ) { sprintf( oname, "%s_X_sp_%02d_%02d.%05d", mfname, m+1, p+1, c_rep+1) ; }
							else { sprintf( oname, "%s_X_st_%02d.%05d", mfname, m - par->nsp + 2, c_rep+1) ; }
							write_data_in_file( oname, presp, 'X', c_rep, MULTIPLICITY[m], sdummy, single_sum_weight_factor, VN[m][p], VN_isl[m][p], VN_osl[m][p], VN_dm[m][p], VN_irr[m][p], var_name(dS_X), dS_X[m][p]) ;  
						}

						if  ( Xn <= 0 )
						{
							/* single particles / stacks */
							if ( m < par->nsp ) { sprintf (oname, "%s_n_sp_%02d_%02d.%05d", mfname, m+1, p+1, c_rep+1) ; }
							else { sprintf (oname, "%s_n_st_%02d.%05d", mfname, m - par->nsp + 2, c_rep+1) ; }
							write_data_in_file( oname, presp, 'n', c_rep, MULTIPLICITY[m], sdummy, single_sum_weight_factor, VN[m][p], VN_isl[m][p], VN_osl[m][p], VN_dm[m][p], VN_irr[m][p], var_name(dS_n), dS_n[m][p], var_name(dB_n), dB_n[m][p]) ;
						}
					}
				}
			}


			/* write S_X/n, B_n to files at intermediate time steps if the flag is set */
			if (intermediatesteps_flag)
			{
				/* do the output every imsteps of (c_rep+1), even if it is the last repeat.
				   This allows later an easier batch processing of the files in order to 
				   show a dependence of results on the amount of performed repeat steps. */
				if ( ( (c_rep+1) % imsteps == 0 ) ) /* && ( (c_rep+1) != par->nr ) ) */
				{
					for ( unsigned int m=0; m<nspsp; ++m)
					{
						MULTIPLICITY[m] = 1 ;
						MULTIPLICITY[m] *= nav_ac ;
						MULTIPLICITY[m] *= nav_pc ;
						sprintf ( sdummy2, "%sMULTIPLICITY = %d ( nav_ac = %d, nav_pc = %d)\n", presp, MULTIPLICITY[m], nav_ac, nav_pc) ;

						sprintf( sdummy, "%s\n%sintermediate step output at (c_rep+1) = %d\n", presp, presp, c_rep+1) ;
						strcat( sdummy, sdummy2) ;

						if ( Xn >= 0 )
						{
							/* single particles / stacks */
							if ( m < par->nsp ) { sprintf( oname, "%s_X_sp_%02d.%05d_im", mfname, m+1, c_rep+1) ; }
							else { sprintf( oname, "%s_X_st_%02d.%05d_im", mfname, m - par->nsp + 2, c_rep+1) ; }
							write_data_in_file( oname, presp, 'X', c_rep, MULTIPLICITY[m], sdummy, sum_weight_factor, Vtot_cry[m], Vtot_isl[m], Vtot_osl[m], Vtot_dm[m], Vtot_irr[m], var_name(S_X), S_X[m]) ;

							if ( m < par->nsp ) { sprintf( oname, "%s_Y_X_sp_%02d.%05d_im", mfname, m+1, c_rep+1) ; }
							else { sprintf( oname, "%s_Y_X_st_%02d.%05d_im", mfname, m - par->nsp + 2, c_rep+1) ; }
							write_Y( oname, presp, 'X', c_rep, MULTIPLICITY[m], sdummy, sum_weight_factor, Vtot_cry[m], Vtot_isl[m], Vtot_osl[m], Vtot_dm[m], Vtot_irr[m], var_name(Yc_X), Yc_X[m]) ;
						}

						if ( Xn <= 0 )
						{
							/* single particles / stacks */
							if ( m < par->nsp ) { sprintf( oname, "%s_n_sp_%02d.%05d_im", mfname, m+1, c_rep+1) ; }
							else { sprintf( oname, "%s_n_st_%02d.%05d_im", mfname, m - par->nsp + 2, c_rep+1) ; }
							write_data_in_file( oname, presp, 'n', c_rep, MULTIPLICITY[m], sdummy, sum_weight_factor, Vtot_cry[m], Vtot_isl[m], Vtot_osl[m], Vtot_dm[m], Vtot_irr[m], var_name(S_n), S_n[m], var_name(B_n), B_n[m]) ;

							if ( m < par->nsp ) { sprintf( oname, "%s_Y_n_sp_%02d.%05d_im", mfname, m+1, c_rep+1) ; }
							else { sprintf( oname, "%s_Y_n_st_%02d.%05d_im", mfname, m - par->nsp + 2, c_rep+1) ; }
							write_Y( oname, presp, 'n', c_rep, MULTIPLICITY[m], sdummy, sum_weight_factor, Vtot_cry[m], Vtot_isl[m], Vtot_osl[m], Vtot_dm[m], Vtot_irr[m], var_name(Yc_n), Yc_n[m], var_name(Yi_n), Yi_n[m]) ;
						}
					}
				}
			}


			if ( time_flag > 1 ) { fprintf( stdout ,"PID %d | R %4d | %-7s | T %5.1f %% | %-30s\n", pid, c_rep, "Done", pr, rtstr ) ; fflush( stdout ) ; }
			if ( time_flag > 1 ) { time(&now) ; }
			if (log_flag )
			{
				if ( time_flag > 1 ) { fprintf( logfile, "%sc_rep: %d done (%ld seconds)\n\n", "\t", c_rep, (long int)(difftime( now, last)) ) ; }
				else { fprintf( logfile, "%sc_rep: %d done\n\n", "\t", c_rep) ; }
			}
			if ( time_flag > 1 ) { time(&last) ; }
			
			/* At the the end of the repeats write S_X/n, B_n */
			/* Provide an option to extend the repeats */
			if ( ( c_rep + 1 ) == par->nr )
			{
				/* Write S to files at the end of the repeats. 
				   If the repeats are extended the files will be overwritten
				   and MULTIPLICITY will be increased */

				for ( unsigned int m=0; m<nspsp; ++m)
				{
					MULTIPLICITY[m] = 1 ;
					MULTIPLICITY[m] *= nav_ac ;
					MULTIPLICITY[m] *= nav_pc ;
					sprintf ( sdummy2, "%sMULTIPLICITY = %d ( nav_ac = %d, nav_pc = %d)\n", presp, MULTIPLICITY[m], nav_ac, nav_pc) ;

					sprintf( sdummy, "%s\n%sFinal step output at (c_rep+1) = %d\n", presp, presp, c_rep+1) ;
					strcat( sdummy, sdummy2) ;

					if ( Xn >= 0 )
					{
						/* single particles / stacks */
						if ( m < par->nsp ) { sprintf( oname, "%s_X_sp_%02d.dat", mfname, m+1) ; }
						else { sprintf( oname, "%s_X_st_%02d.dat", mfname, m - par->nsp + 2) ; }
						write_data_in_file( oname, presp, 'X', c_rep, MULTIPLICITY[m], sdummy, sum_weight_factor, Vtot_cry[m], Vtot_isl[m], Vtot_osl[m], Vtot_dm[m], Vtot_irr[m], var_name(S_X), S_X[m]) ; 

						if ( m < par->nsp ) { sprintf( oname, "%s_Y_X_sp_%02d.dat", mfname, m+1) ; }
						else { sprintf( oname, "%s_Y_X_st_%02d.dat", mfname, m - par->nsp + 2) ; }
						write_Y( oname, presp, 'X', c_rep, MULTIPLICITY[m], sdummy, sum_weight_factor, Vtot_cry[m], Vtot_isl[m], Vtot_osl[m], Vtot_dm[m], Vtot_irr[m], var_name(Yc_X), Yc_X[m]) ;
					}

					if ( Xn <= 0 )
					{
						/* single particles / stacks */
						if ( m < par->nsp ) { sprintf( oname, "%s_n_sp_%02d.dat", mfname, m+1 ) ; }
						else { sprintf( oname, "%s_n_st_%02d.dat", mfname, m - par->nsp + 2) ; }
						write_data_in_file( oname, presp, 'n', c_rep, MULTIPLICITY[m], sdummy, sum_weight_factor, Vtot_cry[m], Vtot_isl[m], Vtot_osl[m], Vtot_dm[m], Vtot_irr[m], var_name(S_n), S_n[m], var_name(B_n), B_n[m]) ;

						if ( m < par->nsp ) { sprintf( oname, "%s_Y_n_sp_%02d.dat", mfname, m+1) ; }
						else { sprintf( oname, "%s_Y_n_st_%02d.dat", mfname, m - par->nsp + 2) ; }
						write_Y( oname, presp, 'n', c_rep, MULTIPLICITY[m], sdummy, sum_weight_factor, Vtot_cry[m], Vtot_isl[m], Vtot_osl[m], Vtot_dm[m], Vtot_irr[m], var_name(Yc_n), Yc_n[m], var_name(Yi_n), Yi_n[m]) ;
					}
				}


				/* -e flag, ignored in -silent mode */
				if (extension_flag)
				{
					fprintf( stdout, "The data has been written from the last (c_rep+1)=%d repeats. \n", c_rep+1) ;
					fprintf( stdout, "Do you want to proceed with repeats (y) or finish (n)? :\t");
					while (true)
					{
						scanf("%s",extensiondecision);
						if (!strcmp(extensiondecision,"y"))
						{
							fprintf( stdout, "Type a number of additional repeats:\t");
							while (true)
							{
								scanf("%s",extensionrepeats);
								repeatsincrease = (unsigned int) strtol (extensionrepeats, NULL, 10) ;
								if (repeatsincrease>0)
								{
									par->nr += repeatsincrease ;
									fprintf( stdout, "Increasing par->nr by %d up to %d \n", repeatsincrease, par->nr );
									break;
								}
								else
								{
									fprintf( stdout, "The input must be a postive integer!\n");
									fprintf( stdout, "Type a number of additional repeats:\t ");
								}
							}
							break;
						}
						else if (!strcmp(extensiondecision,"n"))
						{
							break;
						}
						else
						{
							fprintf( stdout, "Type 'y' or 'n':\t");
						}
					}
				}
			}

		/* end FIRST loop for over repeats c_rep */
		}

		if ( time_flag > 1 ) { fprintf( stdout ,"PID %d | R %4d | %-7s | %-42s\n", pid, par->nr, "Done", "Done") ; fflush( stdout ) ; }
		if ( log_flag ) { fprintf( logfile, "Done\n\n") ; }

		/********************/
		/* END OF AVERAGING */
		/********************/


		/* deallocate memory from stackcpp and then all class members */
		if (log_flag ) { fprintf( logfile, "Deallocating memory\n") ; }

		free(n1);
		free(n2);
		free(n3);
		free(r3);
		free(ind_r3);

		for ( unsigned int i=0; i<par->nms; ++i)
		{
			free(v1[i]);
			free(v2[i]);
			free(v3[i]);
			free(w1[i]);
			free(w2[i]);
			free(wt1[i]);
			free(wt2[i]);
		}
		free(v1);
		free(v2);
		free(v3);
		free(w1);
		free(w2);
		free(wt1);
		free(wt2);

		free(kf1);
		free(kf2);
		free(kf3);
		free(okf1);
		free(okf2);
		free(okf3);

		for ( unsigned int i=0; i<par->nsp; ++i)
		{
			for ( unsigned int j=0; j<par->nms; ++j)
			{
				free(w3[i][j]) ;
				free(wt3[i][j]) ;
			}
			free(w3[i]);
			free(wt3[i]);
		}
		free(w3);
		free(wt3);

		if ( Xn >= 0 )
		{
			free(drerho_dm_osl) ;
			free(drerho_osl_isl) ;
			free(drerho_isl) ;
		}
		if ( Xn <= 0 )
		{
			free(dsld_dm_osl) ;
			free(dsld_osl_isl) ;
			free(dsld_isl) ;
		}

		for ( unsigned int i=0; i<par->nsp; ++i)
		{
			free(V[i]) ;
			free(V_isl[i]) ;
			free(V_osl[i]) ;
		}
		free(V) ;
		free(V_isl) ;
		free(V_osl) ;

		for ( unsigned int i=0; i<nspsp; ++i)
		{
			free(VN[i]) ;
			free(VN_isl[i]) ;
			free(VN_osl[i]) ;
			free(VN_dm[i]) ;
			free(VN_irr[i]) ;
		}
		free(VN) ;
		free(VN_isl) ;
		free(VN_osl) ;
		free(VN_dm) ;
		free(VN_irr) ;

		free(Vtot_cry) ;
		free(Vtot_isl) ;
		free(Vtot_osl) ;
		free(Vtot_irr) ;
		free(Vtot_dm) ;

		for ( unsigned int i=0; i<par->nms; ++i) { free(R[i]); }
		free(R);

		free(d1);
		free(d2);
		free(D);

		free(MULTIPLICITY) ;

		/* the following arrays have been already deallocated within the parallel region
		   GA, P, P_isl, P_osl, dE_X, dE_n, exp_N3Qv3, exp_N3Qw3, exp_N3Qwt3, cdummy41, cdummy42, cdummy43, exp_RQ,
		   dS_X_parallel_thread, dS_n_parallel_thread, Yc_X_parallel_thread, Yc_n_parallel_thread
		*/

		/* free S_X/n (nspsp x rho/sld_multi.n x np ), 
		   dS_X/n ( nspsp x nms/1 x rho/sld_multi.n x np ),
		   and similarly for B_n, dB_n, dC_n */
		if ( Xn >= 0 )
		{
			for ( unsigned int j=0; j<nspsp; ++j)
			{
				for ( unsigned int q=0; q<rho_multi.n; ++q) { free(S_X[j][q]) ; }
				free(S_X[j]) ; 

				if ( j<par->nsp )
				{
					for ( unsigned int k=0; k<par->nms; ++k)
					{
						for ( unsigned int q=0; q<rho_multi.n; ++q) { free(dS_X[j][k][q]) ; }
						free(dS_X[j][k]) ;
					}
				}
				else
				{
					for ( unsigned int q=0; q<rho_multi.n; ++q) { free(dS_X[j][0][q]) ; }
					free(dS_X[j][0]) ;
				}
				free(dS_X[j]) ;
			}
			free(S_X) ;
			free(dS_X) ;
		}
		if ( Xn <= 0 )
		{
			for ( unsigned int j=0; j<nspsp; ++j)
			{
				for ( unsigned int q=0; q<sld_multi.n; ++q) { free(S_n[j][q]) ; }
				free(S_n[j]) ;

				if ( j<par->nsp )
				{
					for ( unsigned int k=0; k<par->nms; ++k)
					{
						for ( unsigned int q=0; q<sld_multi.n; ++q) { free(dS_n[j][k][q]) ; }
						free(dS_n[j][k]) ;
					}
				}
				else
				{
					for ( unsigned int q=0; q<sld_multi.n; ++q) { free(dS_n[j][0][q]) ; }
					free(dS_n[j][0]) ;
				}
				free(dS_n[j]) ;
			}
			free(S_n) ;
			free(dS_n) ;

			for ( unsigned int j=0; j<nspsp; ++j)
			{
				if ( j < par->nsp ) { for ( unsigned int k=0; k<par->nms; ++k) { free(dB_n[j][k]) ; } }
				else { free(dB_n[j][0]) ; }
				free(dB_n[j]) ;
				free(B_n[j]) ; 
			}
			free(dB_n) ;
			free(B_n) ;
		}

		/* free dYc(i)_X/n which are a double complex matrices of size ( par->nsp x 10(4) x par->np ) */
		if ( Xn >= 0 )
		{
			for ( unsigned int i=0; i<nspsp; ++i) 
			{
				for ( int j=0; j<10; ++j) { free(Yc_X[i][j]) ; }
				free(Yc_X[i]) ;
			}
			free(Yc_X);
		}
		if ( Xn <= 0 )
		{
			for ( unsigned int i=0; i<nspsp; ++i)
			{
				for (int j=0; j<10; ++j) { free(Yc_n[i][j]) ; }
				free(Yc_n[i]) ;
				free(Yi_n[i]) ;
			}
			free(Yc_n);
			free(Yi_n);
		}

		/* deallocate GSL uniform and discrete PRNG if par->td2==2 */
		if ( par->td2 == 2 )
		{
			gsl_ran_discrete_free(d_rand) ;
			gsl_rng_free(GSL_RNG) ;
		}

		free_class_members() ;

		if (log_flag ) { fprintf( logfile, "done\n\n") ; }

		/* end of deallocating memory */


		/* write terminating time to logfile */
		if (log_flag )
		{
			struct tm * timeinfo ;
			time(&now) ;

			timeinfo = localtime(&now) ;
			fprintf( logfile, "Program terminated at %s \n", asctime(timeinfo) );
		}

		return (0) ;
	}

	/* free-function, reserved for variables that are class members */
	void free_class_members()
	{
		free(atom) ;

		if ( Xn >= 0 )
		{
			free(rho_isl) ;
			free(rho_osl) ;
			free(rho_dm) ;
		}
		if ( Xn <= 0 )
		{
			free(sld_isl) ;
			free(sld_osl) ;
			free(sld_dm) ;
	
			free(icsd_isl) ;
			free(icsd_osl) ;
			free(icsd_dm) ;
		}
		if ( rho_multi.multiq )
		{
			free(rho_multi.XRAY_RHO_INLAY_file) ;
			free(rho_multi.XRAY_RHO_OUTLAY_file) ;
			free(rho_multi.XRAY_RHO_DM_file) ;
		}
		if ( sld_multi.multiq )
		{
			free(sld_multi.NEUT_SLD_INLAY_file) ;
			free(sld_multi.NEUT_SLD_OUTLAY_file) ;
			free(sld_multi.NEUT_SLD_DM_file) ;
		}

		/* deallocate FFF_X, FFF_n, f0table and eQvj */
		free_structure_amplitude() ;

		if ( Xn >= 0)
		{
			for (int i=0; i<=numf0; ++i) { free(f0table[i]) ; }

			free(f0table) ;
		}
		for ( unsigned int j=0; j<nav_pc; ++j)
		{
			if ( ( par->av_mode == 0 ) || ( par->av_mode == 2 ) || ( par->av_mode == 3 ) )
			{
				for ( unsigned int k=0; k<nav_ac; ++k) { free(eQvj[j][k]) ; }
			}
			else if ( par->av_mode == 1 )
			{
				for ( unsigned int k=0; k<nav_r3; ++k) { free(eQvj[j][k]) ; }
			}
			free(eQvj[j]) ;
		}
		free(eQvj) ;

		/* free databases */
		if ( Xn <= 0 ) { free_neutron_scattering_lengths_and_cross_sections() ; }
		free_standard_atomic_weight_and_relative_atomic_mass() ;
		free_cif_dictionary() ;
		free_cif_file() ;
	}


	void eval_cmd (int carg, char **varg)
	{
		int idummy ;
		double ddummy ;
		char sdummy[1024] ;

		FILE* outf ;

		ecarg = &carg ;
		evarg = varg ;

		/* allocate par struct */
		par = (struct stack_parameters *) calloc( 1, sizeof (struct stack_parameters)) ;


		/* SETUP SOME DEFAULT VALUES THAT MIGHT BE OVERWRITTEN IN THE INPUT */

		log_flag = false ;
		openmp_flag = false ;
		
		singlefilesoutput_flag = false ;
		intermediatesteps_flag = false ;
		extension_flag = false ;

		/* -sym cif */
		sym_cif_mode = false ;

		/* -sym ccp4 */
		sym_ccp4_mode = false ;
		def_ccp4_bash_file = "ccp4_bash.sh" ;
		strcpy( ccp4_bash_file, def_ccp4_bash_file) ;

		def_cif_pdb_file = "cif.pdb" ;
		strcpy( cif_pdb_file, def_cif_pdb_file) ;

		def_pdb_file = "pdbset.pdb" ;
		strcpy( pdb_file, def_pdb_file) ;

		/* -sym pdb */
		sym_pdb_mode = false ;

		symop_userdef = false ;

		/* base transformation */
		bt_flag = false ;

		/* write cif-file for Jmol */
		jmol_cif_flag = false ;

		/* perform orientational average of structure amplitude and factor */
		orav_flag = false ;

		/* par-file given by numbers */
		par_num_def = false ;

		/* read/write structure amplitude FFF */
		Fr_flag = false ;
		num_Fr_flag = 0 ;
		Fw_flag = false ;
		FNoMem_flag = false ;

		/* testing purposes */
		test_XnEquiv_flag = false ;
		test_I_at_s_equal_zero_flag = false ;

		/* log random values for n1,n2,n3,... */
		distr_flag = false ;

		/* create atomic positions in the pcr file format for FullProf Suite */
		pcr_flag = false ;

		/* apply sin theta weight factor or not */
		weight_flag = false ;

		init_n_rand_userdef = false ;
		init_d_rand_userdef = false ;
		init_MTRand_userdef = false ;

		/* default concentrations */
		conc_userdef = false ;
		conc_dm = 1.0 ;
		conc_cry = 0.01 ;

		/* -mem -time flags */
		#ifdef __FreeBSD__
		mem_flag = 0 ;
		time_flag = 0 ;
		#elif __linux__
		mem_flag = 1 ;
		time_flag = 2 ;
		#endif

		/* X/n flag */
		Xn = 0 ;
		bool xlock = false ;
		bool nlock = false ;

		/* flags for parameters that MUST be defined by the user */
		bool av_isdef = false ;
		bool cif_isdef = false ;
		bool func_isdef = false ;
		bool mfname_isdef = false ;
		bool par_isdef = false ;

		cis_isdef = false ;
		nis_isdef = false ;

		/* function */
		bool func_found = false ;
		int func_num ; /* Number of selected function */


		/* print help (usage) information if only one input argument occurs */
		if (carg < 2)
		{
			usage(varg[0]) ;
			exit(1) ;
		}


		/* allow to overwrite existing output */
		overwrite_flag = false ;

		/* silent mode */
		silent_flag = false ;

		/* default input/output folders */
		def_ifo = "in/" ;
		strcpy( ifo, def_ifo) ;
		def_ofo = "out/" ;
		strcpy( ofo, def_ofo) ;

		/* check for -ofo and -ow and +o at the beginning */
		for ( int i=1; i<carg; ++i)
		{
			if ( !strcmp(varg[i], "-l") )
			{
				/* -l, generates a logfile of the name [mfname].log */
				log_flag = true ;
			}
			if ( !strcmp(varg[i], "+o") )
			{
				/* +o mfname, set user provided master-filename, e.g. +o test */
				if ( ++i<carg )
				{
					if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
					strcpy( mfname, varg[i]) ;
				}
				else { XNDIFF_ERROR(1) ; }

				mfname_isdef = true ;
			}
			else if ( !strcmp(varg[i], "-ofo") )
			{
				/* -ofo, change default output folder */
				if ( ++i<carg ) 
				{
					if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
					strcpy( ofo, varg[i]) ;
				}
				else { XNDIFF_ERROR(1) ; }
			}
			else if ( !strcmp(varg[i], "-ow") )
			{
				overwrite_flag = true ;
			}
			else if ( !strcmp(varg[i], "-silent") )
			{
				overwrite_flag = true ;
				silent_flag = true ;
			}
		}
		/* log_flag, mfname, ofo, overwrite_flag, silent_flag are now already well defined */
		if ( !mfname_isdef) { XNDIFF_ERROR(78) ; }


		if ( !silent_flag ) { fprintf( stdout, "XNDiff --- A calculation program for X-ray and neutron powder diffraction patterns\n") ; }
		if ( !silent_flag ) { fprintf( stdout, "Version %s    %s by %s \n", XNDIFF_VERSION, XNDIFF_DATE, XNDIFF_AUTHOR) ; }


		time_t now;
		struct tm * timeinfo;
		time(&now) ;
		timeinfo = localtime(&now) ;
		strcpy( startingtime, asctime(timeinfo) ) ;


		/* check/create ifo/ofo */
		check_ifo_ofo() ;

		/* overwrite protection, if no -ow or -silent flag was explicitely provided */
		overwrite_protection() ;


		/* start log file, logfile is the global FILE-pointer to the log file */
		if ( log_flag )
		{
			/* http://www.cplusplus.com/reference/clibrary/cstdio/sprintf/ */
			/* create filename for the logfileas <mfname>.log */
			sprintf( sdummy, "%s%s.log", ofo, mfname) ;
			if ( ( logfile = fopen( sdummy, "w")) == NULL) { XNDIFF_ERROR(13) ; }

			write_call_arguments( logfile, "") ;
		}

		/* read option flags from cmd
		   Examples:
			-l +o test +cif betaPPP.cif cif_core.dic -sym ccp4 ccp4config ccp4symop.txt -distr -mem 0 -time 2 -init_n_rand -10112 +av 0 0.0 90.0 0.5 0.0 360.0 15.0 -sin +par XNDiff004.par +f stackcpp 120.0 20.0 100.0 20.0 2.522 0.382 0.0 20.0 0.0 20.0 33.2 1.5 0 1 0 0.001 0.5 50 10 0 0 1 0 0 1 10 5
			-l +o test +cif betaPPP.cif cif_core.dic -sym pdb betaPPP.pdb -distr -mem 0 -time 2 -init_n_rand -10112 +av 2 0.0 90.0 0.0 360.0 32400 +par 0.8 1.3 310.0 350.0 333.0 -0.5 1.2 6.34 0.0 0.0 0.0 +f stackcpp 120.0 20.0 100.0 20.0 2.522 0.382 0.0 20.0 0.0 20.0 33.2 1.5 0 1 0 0.001 0.5 50 10 0 0 1 0 0 1 10 5
			-l +o test +cif betaPPP.cif cif_core.dic -sym symop.txt -distr -init_n_rand -10112 -init_MTRand 100 +av 3 0.0 90.0 0.0 360.0 32400 +par 0.8 1.3 310.0 350.0 333.0 0.04 0.95 -X -conc par +f stackcpp 120.0 20.0 100.0 20.0 2.522 0.382 0.0 20.0 0.0 20.0 33.2 1.5 0 1 0 0.001 0.5 50 10 0 0 1 0 0 1 10 5
			-l +o test +cif betaPPP.cif cif_core.dic -sym symop.txt -distr -init_n_rand -10112 -init_MTRand 100 +av 3 0.0 90.0 0.0 360.0 32400 +par 0.8 1.3 rhosld002.dat 1 rhosld002.dat 2 rhosld002.dat 3 -X +f stackcpp 120.0 20.0 100.0 20.0 2.522 0.382 0.0 20.0 0.0 20.0 33.2 1.5 0 1 0 0.001 0.5 50 10 0 0 1 0 0 1 10 5
			-l +o test -openmp -ow +cif betaPPP.cif cif_core.dic -distr -orav -Fr -pcr -jmolcif -init_n_rand -10112 +av 1 0.0 90.0 0.5 0.0 60.0 -90 90.0 181 +par 0.8 1.3 -0.5 1.2 6.34 0.0 0.0 0.0 -n -conc def +f stackcpp 120.0 20.0 100.0 20.0 2.522 0.382 0.0 20.0 0.0 20.0 33.2 1.5 0 1 0 0.001 0.5 500 100 0 0 1 0 0 1 10 5 -s -z 10 -e
			-l +o test -ofo bt_twin_jmol_test +cif betaPPP.cif cif_core.dic +av 2 0.0 90.0 0.0 360.0 32400 +par XNDiff004.par -bt bt_twin -jmolcif +f stackcpp 120.0 20.0 100.0 20.0 2.522 0.382 0.0 20.0 0.0 20.0 33.2 1.5 0 1 0 0.001 0.5 50 10 0 0 1 0 0 1 10 5
			-l +o test +cif betaPPP-d98.cif cif_core.dic -sym symop.txt +av 0 0.0 90.0 0.5 0.0 360.0 15.0 -sin -init_n_rand -10112 -distr +par 0.8 1.3 310.0 350.0 333.0 rhosld002.dat 1 rhosld002.dat 2 rhosld002.dat 3 rhosld002.dat 4 rhosld002.dat 5 rhosld002.dat 6 +f stackcpp 120.0 20.0 100.0 20.0 2.522 0.382 0.0 20.0 0.0 20.0 25.0 1.5 0 1 0 0.001 0.6 60 5 0 0 1 0 0 1 10 5 -ow -jmolcif -z 2 -e -s -pcr test.pcr
			-l +o test +cif betaPPP.cif cif_core.dic -sym pdb betaPPP.pdb +av 0 0.0 90.0 0.5 0.0 360.0 1.0 -sin -init_n_rand -10112 -conc def -distr +par 1.0 0.5 310.0 350.0 333.0 -0.4939 1.1995 6.3664 0.44911 0.44911 0.01102 +f stackcpp 120.0 20.0 100.0 20.0 0.0 0.0 0.0 20.0 0.0 20.0 2.01877 0.362735 0 1 3 0.001 0.6 600 100 0 0 1 0 0 1 7 3
			-l +o debug_distr_010_N100_Nsp5_Nst5 +cif betaPPP.cif cif_core.dic -sym pdb betaPPP.pdb +av 0 0.0 90.0 10.0 0.0 360.0 10.0 -sin -init_n_rand -10112 -conc def -distr +par 1.0 0.5 310.0 350.0 333.0 -0.4939 1.1995 6.3664 0.44911 0.44911 0.01102 +f stackcpp 120.0 20.0 100.0 20.0 2.05425 0.246221 0.0 20.0 0.0 20.0 30.0 1.5 0 1 0 0.001 0.6 10 100 0 0 1 0 0 1 5 5 -openmp
			-l +o debug_distr_020_N100_Nsp5_Nst5 +cif betaPPP.cif cif_core.dic -sym pdb betaPPP.pdb -cis 0.2 0.4 0.3 0.1 0.05 +av 0 0.0 90.0 10.0 0.0 360.0 10.0 -sin -init_n_rand -10112 -conc def -distr +par 1.0 0.5 310.0 350.0 333.0 -0.4939 1.1995 6.3664 0.44911 0.44911 0.01102 +f stackcpp 120.0 20.0 100.0 20.0 0.0 0.0 0.0 20.0 0.0 20.0 30.0 1.5 0 2 0 0.001 0.6 10 100 0 0 1 0 0 1 5 5 -openmp
			-l +o debug_distr_023_N100_Nsp5_Nst5 +cif betaPPP.cif cif_core.dic -sym pdb betaPPP.pdb -cis 0.2 0.4 0.3 0.1 0.05 +av 0 0.0 90.0 10.0 0.0 360.0 10.0 -sin -init_n_rand -10112 -conc def -distr +par 1.0 0.5 310.0 350.0 333.0 -0.4939 1.1995 6.3664 0.44911 0.44911 0.01102 +f stackcpp 120.0 20.0 100.0 20.0 0.0 0.0 0.0 20.0 0.0 20.0 0.0 0.0 0 2 3 0.001 0.6 10 100 0 0 1 0 0 1 5 5 -openmp
			-l +o debug_distr_002_N100_Nsp5_Nst5 +cif betaPPP.cif cif_core.dic -sym pdb betaPPP.pdb -cis 0.2 0.4 0.3 0.1 0.05 +av 0 0.0 90.0 10.0 0.0 360.0 10.0 -sin -init_n_rand -10112 -conc def -distr +par 1.0 0.5 310.0 350.0 333.0 -0.4939 1.1995 6.3664 0.44911 0.44911 0.01102 +f stackcpp 120.0 20.0 100.0 20.0 12.0 3.0 0.0 20.0 0.0 20.0 0.0 0.0 0 0 2 0.001 0.6 10 100 0 0 1 0 0 1 5 5 -openmp
			-l +o debug_distr_020_N100_Nsp5_Nst5_v2 +cif betaPPP.cif cif_core.dic -sym pdb betaPPP.pdb -cis test.cis +av 0 0.0 90.0 10.0 0.0 360.0 10.0 -sin -init_n_rand -10112 -conc def -distr +par 1.0 0.5 310.0 350.0 333.0 -0.4939 1.1995 6.3664 0.44911 0.44911 0.01102 +f stackcpp 120.0 20.0 100.0 20.0 0.0 0.0 0.0 20.0 0.0 20.0 30.0 1.5 0 2 0 0.001 0.6 10 100 0 0 1 0 0 1 5 5 -openmp
			-l +o vbasis +cif betaPPP.cif cif_core.dic -sym symop_vbasis.txt -jmolcif +av 0 0.0 90.0 0.5 0.0 360.0 15.0 -sin +par 0.8 1.3 310.0 350.0 333.0 -0.5 1.2 6.34 0.0 0.0 0.0 +f stackcpp 120.0 20.0 100.0 20.0 2.522 0.382 0.0 20.0 0.0 20.0 15.5 0.25 0 1 0 0.001 0.6 60 1 0 0 1 0 0 1 10 5 -openmp
		*/

		if ( log_flag ) { fprintf( logfile, "Evaluate command line arguments\n") ; fflush (logfile) ; }
		for ( int i=1; i<carg; ++i)
		{
			/* allow only flags of the kind +par -sin -crazy44 +4u
			   numbers and text without +/- at the front will cause the program to terminate
			*/
			if ( !silent_flag ) { fprintf( stdout, "%s\n", varg[i]) ; fflush(stdout) ; }
			if ( log_flag ) { fprintf( logfile, "\t%s\n", varg[i]) ; fflush (logfile) ; }

			if ( ( (*varg[i] == '-') || (*varg[i] == '+') ) && ( is_numeric(varg[i],false) == false ) )
			{
				switch (varg[i][1])
				{
					/* +av */
					case 'a':
						if ( !strcmp(varg[i], "+av") )
						{
							/* av_mode mode for orientational averaging (powder average) */

							par->av_mode = (int) strtol (varg[++i], NULL, 10) ;
							if ( par->av_mode > MAX_AV_MODE ) { XNDIFF_ERROR(17) ; }

							if ( par->av_mode == 0 )
							{
								/* e.g. 0 0.0 90.0 0.5 0.0 360.0 15.0 */
								/* use a grid of polar and azimuthal angles */
								/* av_pol[0] av_pol[1] av_pol[2] range and step of polar angle vs G [degree] e.g. 0.0 90.0 0.5
								   av_azi[0] av_azi[1] av_azi[2] range and step of azimuthal angle vs. d1 [degree] e.g. 0.0 360.0 360.0 
								*/
								for ( unsigned int j=0; j<3; ++j) { par->av_pol[j] = strtod (varg[++i], NULL) ; }
								for ( unsigned int j=0; j<3; ++j) { par->av_azi[j] = strtod (varg[++i], NULL) ; }
							}
							else if ( par->av_mode == 1 )
							{
								/* e.g. 1 0.0 90.0 0.5 0.0 60.0 -90 90.0 181 */
								/* theta grid and random r3 rotation */
								/* r31 r32: mean r3 rotation and std. dev. crystals in stack on Gs (vs d1) [degree]
								   av_r3[0] av_r3[1] n_r3: range of the r3 angles av_r3[0]...av_r3[1] and numbers of grid points
								   [The randomly chosen r3 angles will be rounded to the grid specified by av_r3[0] av_r3[1] n_r3 !]
								*/
								/* by default set av_azi to {0.0 360.0 360.0} */
								for ( unsigned int j=0; j<3; ++j) { par->av_pol[j] = strtod (varg[++i], NULL) ; }
								par->r31 = strtod (varg[++i], NULL) ;
								par->r32 = strtod (varg[++i], NULL) ;
								for ( unsigned int j=0; j<2; ++j) { par->av_r3[j] = strtod (varg[++i], NULL) ; }
								par->n_r3 = (int) strtol (varg[++i], NULL, 10) ;
								
								par->av_azi[0] = 0.0 ;
								par->av_azi[1] = 360.0 ;
								par->av_azi[2] = 360.0 ;
							}
							else if ( ( par->av_mode == 2 ) || ( par->av_mode == 3 ) )
							{
								/* e.g. 2 0.0 90.0 0.0 360.0 4320 
									3 0.0 90.0 0.0 360.0 32400 
								   -sin option should not be used in MC mode !!!
								*/
								/* specify range of polar and azimuthal angles and number of MC points */
								/* MC mode */
								for ( unsigned int j=0; j<2; ++j) { par->av_pol[j] = strtod (varg[++i], NULL) ; }
								for ( unsigned int j=0; j<2; ++j) { par->av_azi[j] = strtod (varg[++i], NULL) ; }
								par->n_MC = (int) strtol (varg[++i], NULL, 10) ;
							}

							av_isdef = true ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -bt */
					case 'b':
						if ( !strcmp(varg[i], "-bt") )
						{
							/* e.g. -bt bt_twin or -bt bt_double */
							if ( ++i<carg )
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								strcpy( bt, varg[i]) ;
								bt_flag = true ;
							}
							else { XNDIFF_ERROR(1) ; }
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* +cif -conc -cis */
					case 'c':
						if ( !strcmp(varg[i], "+cif") )
						{
							/* +cif cif-file cif-dic, read the cif-file and the 
							   the cif-dictionary, e.g. +cif betaPPP.cif cif_core.dic */

							/* cif-file */
							if ( ++i<carg )
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								strcpy( cif_file, varg[i]) ;
							}
							else { XNDIFF_ERROR(1) ; }
							/* cif-dic */
							if ( ++i<carg )
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								strcpy( cif_dic, varg[i]) ;
							}
							else { XNDIFF_ERROR(1) ; }

							cif_isdef = true ;
						}
						else if ( !strcmp(varg[i], "-cis") )
						{
							/* read the c_i's if par->td==2 directly from the command line or row by row from an ascii file
							   e.g. -cis 0.2 0.4 0.3 0.1 0.05
							   e.g. -cis PPP_180x360_P02_02_ST600_100_N100.cis
							*/
							/* text or numbers */
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) )
								{
									while ( i < carg )
									{
										if ( is_numeric(varg[i]) ) { cis.push_back( strtod( varg[i], NULL) ) ; ++i ; }
										else { break ; }
									}
									--i ;
									cis_isdef = true ;
								}
								else
								{
									/* filename or flags */
									if ( *varg[i] != '-' && *varg[i] != '+')
									{
										/* read cis from file */
										read_column_from_ascii_file( string(varg[i]), ' ', 1, cis) ;
										cis_isdef = true ;
									}
									else
									{
										/* cis_isdef = false ; */
										break ;
									}
								}
							}
							else { XNDIFF_ERROR(1) ; }
						}
						else if ( !strcmp(varg[i], "-conc") )
						{
							/* might have been read already for +par -> reading twice is no problem */
							if ( ++i<carg )
							{
								if ( !strcmp( varg[i], "par") ) { conc_userdef = true ; }
								else if ( !strcmp( varg[i], "def") ) { conc_userdef = false ; }
								else { XNDIFF_ERROR(48) ; }
							}
							else { XNDIFF_ERROR(1) ; }
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -distr */
					case 'd':
						if ( !strcmp(varg[i], "-distr") )
						{
							/* -distr, table the generated random numbers for n1. n2, n3, D, d1, d2 in a file [mfname].distr */
							distr_flag = 1 ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -e */
					case 'e':
						if ( !strcmp(varg[i], "-e") )
						{
							/* -e, allows a proceeding on averaging after finishing the repeats. 
							   A question will show up after finishing the repeats whether to extend them 
							   and by which amount of repeats */
							extension_flag = true ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* +f */
					case 'f':
						if ( !strcmp(varg[i], "+f") )
						{
							if ( !( ++i < carg ) ) { XNDIFF_ERROR(1) ; }

							
							for ( unsigned int l=0; l<NOF; ++l)
							{
								/* test varg[i] on func[l].name */
								if (strcmp (func_XPPSA[l].name, varg[i]) == 0)
								{
									func_found = true ;
									func_num = l ;
									switch (func_num)
									{
										case 0: /* stackcpp function
											
											n11 n12 mean particle size along a1 and its standard deviation [nm] e.g. 120.0 20.0 (distribution according to td1)
											n21 n22 mean particle size along a2 and its standard deviation [nm] e.g. 100.0 20.0 (distribution according to td1)
											n31 n32 mean particle size along a3 and its standard deviation [nm] e.g. 2.522 0.382 (distribution according to td2)
											d11 d12 mean displacement and std. dev. of crystals in stack perpendicular to Gs e.g. 0.0 20.0 (distribution according to td1)
											d21 d22 mean displacement and std. dev. of crystals in stack perpendicular to Gs and d1 e.g. 0.0 20.0 (distribution according to td1)
											D1 D2 mean repeat distance along Gs and its standard deviation [nm] e.g. 35.0 1.5 (distribution according to td1)
											td1 type of distribution for n1, n2, d1, d2, D and r3 (+av 1 ...); e.g. 0
											-0 norm.
											-1 log. norm.
											td2 type of distribution for n3
											-0 norm.
											-1 log. norm.
											-2 discrete (requires -cis option)
											stackmode mode for particle stacks e.g. 0
											-0 use td1-distribution for d1,d2,D (kind of default for triglyceride dispersions with S100)
											-1 same as 0; additionally check for collisions with DMIN
											-2 ignore distribution for d1,d2,D; crystals lie directly on top of each other ignoring stabilizer layer (e.g. for DNA or DMPC dispersions)
											-3 same as 2; but including two stabilizer layers
											rs1 rs2 np s-range and number of points for calculation [1/nm] e.g. 0.001 0.5 1000
											nr number of repeats with averaging e.g. 1000
											G(h,k,l) hkl-values representing the direction of the scattering vektors e.g. 0 0 1 
											Gs(hs,ks,ls) hsksls-values of stacking orientation of the crystals e.g. 0 0 1 
											nsp maximal number of layers in single particles e.g. 10
											nms maximal number of particles in stacks e.g. 5

											+f func_XPPSA n11 n12 n21 n22 n31 n32 d11 d12 d21 d22 D1 D2 td1 td2 stackmode rs1 rs2 np nr h k l hs ks ls nsp nms

											*/

											if ( log_flag ) { fprintf( logfile, "\t\tFunction \"stackcpp\" selected (nr.: %d)\n", func_num) ; }
						
											/* check the right number of parameters as defined by func_XPPSA[l].nop */
											if ( i+(int)func_XPPSA[l].nop >= carg ) { XNDIFF_ERROR(5) ; }
											/* check if they are all numeric */
											for ( unsigned int j=(unsigned int)i+1; j<(unsigned int)i+1+func_XPPSA[l].nop; ++j) { if ( !is_numeric(varg[j]) ) { XNDIFF_ERROR(36) ; } }
						
											par->n11 = strtod (varg[++i], NULL) ;
											par->n12 = strtod (varg[++i], NULL) ;
											par->n21 = strtod (varg[++i], NULL) ;
											par->n22 = strtod (varg[++i], NULL) ;
											par->n31 = strtod (varg[++i], NULL) ;
											par->n32 = strtod (varg[++i], NULL) ;

											par->d11 = strtod (varg[++i], NULL) ;
											par->d12 = strtod (varg[++i], NULL) ;
											par->d21 = strtod (varg[++i], NULL) ;
											par->d22 = strtod (varg[++i], NULL) ;
											par->D1 = strtod (varg[++i], NULL) ;
											par->D2 = strtod (varg[++i], NULL) ;

											par->td1 = (int) strtol (varg[++i], NULL, 10) ;
											if (par->td1 < 0 || par->td1 > 1) { XNDIFF_ERROR(85) ; }
											par->td2 = (int) strtol (varg[++i], NULL, 10) ;
											if (par->td2 < 0 || par->td2 > 2) { XNDIFF_ERROR(86) ; }

											par->stackmode = (int) strtol (varg[++i], NULL, 10) ;
											if (par->stackmode < 0 || par->stackmode > 4) { XNDIFF_ERROR(87) ; }

											par->rs1 = strtod (varg[++i], NULL) ;
											par->rs2 = strtod (varg[++i], NULL) ;
											if (par->rs1 > par->rs2)
											{
												ddummy = par->rs1 ;
												par->rs1 = par->rs2 ;
												par->rs2 = ddummy ;
											}
											par->np = (unsigned int) strtol (varg[++i], NULL, 10) ;

											par->nr = (unsigned int) strtol (varg[++i], NULL, 10) ;

											par->h = strtod (varg[++i], NULL) ;
											par->k = strtod (varg[++i], NULL) ;
											par->l = strtod (varg[++i], NULL) ;
											par->hs = strtod (varg[++i], NULL) ;
											par->ks = strtod (varg[++i], NULL) ;
											par->ls = strtod (varg[++i], NULL) ;

											par->nsp = (unsigned int) strtol (varg[++i], NULL, 10) ;
											par->nms = (unsigned int) strtol (varg[++i], NULL, 10) ;
											break ;
									}
								}
							}

							if ( func_found == false ) { XNDIFF_ERROR(4) ; }

							func_isdef = true ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -Fr -Fw -FNoMem */
					case 'F':
						if ( !strcmp(varg[i], "-Fr") )
						{
							/* -Fr Fr_file[0] (-Fr Fr_file[1]), reads saved structure amplitude file(s) into the memory.
							   Checks will be applied later to the file(s) in order to check if the same setup
							   was used and of course which one is for neutrons and which for X-rays */
							Fr_flag = true ;
							if (++i<carg && num_Fr_flag < 2)
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								strcpy( Fr_file[num_Fr_flag], varg[i]) ;
								num_Fr_flag += 1 ;
							}
							else { XNDIFF_ERROR(1) ; }
						}
						else if ( !strcmp(varg[i], "-Fw") )
						{
							/* -Fw, write structure amplitudes from memory to files [mfname]_FFF_X.dat and [mfname]_FFF_n.dat */
							Fw_flag = true ;
						}
						else if ( !strcmp(varg[i], "-FNoMem") )
						{
							/* -FNoMem, FFF will not be computed only once at the beginnig and stored in memory, it will be computed each time in stackcpp */
							/* by this following changes:
							   -disable FFF computation in compute_structure_amplitude() but allow allocation of FFF
							   -free_structure_amplitude() must not be changed
							   -allow compute_orientational_averaged_structure_factor() only after the first c_rep=0 step, when FFF is filled
							   -allow write_structure_amplitude() only after the first c_rep=0 step, when FFF is filled
							   -fully disable read_structure_amplitude()
							*/
							FNoMem_flag = true ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -h, show help from usage function */
					case 'h':
						if ( !strcmp(varg[i], "-h") )
						{
							/* -h, show help from usage function */
							usage (varg[0]) ;
							exit(0) ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -ifo -init_n_rand -init_d_rand -init_MTRand */
					case 'i':
						/* -ifo, change default input folder */
						if ( !strcmp(varg[i], "-ifo") )
						{
							if ( ++i<carg ) 
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								strcpy( ifo, varg[i]) ;
							}
							else { XNDIFF_ERROR(1) ; }
						}
						else if ( !strcmp(varg[i], "-init_d_rand") )
						{
							/* allows initialization of d_rand (n3) with a user supplied long unsigned int number instead of using time */
							/* by this the same random numbers can be generated */
							if ( ++i<carg )
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								idummy = (int) strtol( varg[i], NULL, 10) ;
								if ( idummy < 0 ) { XNDIFF_ERROR(63) ; }
								init_d_rand = (unsigned int)idummy ;
							}
							else { XNDIFF_ERROR(1) ; }
							init_d_rand_userdef = true ;
						}
						else if ( !strcmp(varg[i], "-init_n_rand") )
						{
							/* allows initialization of n_rand (n1,n2,n3,d1,d2,D,r3) with a user supplied negative int number instead of using time */
							/* by this the same random numbers can be generated */
							if ( ++i<carg )
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								init_n_rand = (int) strtol( varg[i], NULL, 10) ;
								if ( init_n_rand >= 0 ) { XNDIFF_ERROR(62) ; }
							}
							else { XNDIFF_ERROR(1) ; }
							init_n_rand_userdef = true ;
						}
						else if ( !strcmp(varg[i], "-init_MTRand") )
						{
							/* allows initialization of MTRand (MC mode Powder Average) with a user supplied unsigned int number instead of using time */
							/* by this the same random numbers can be generated */
							if ( ++i<carg )
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								idummy = (int) strtol( varg[i], NULL, 10) ;
								if ( idummy < 0 ) { XNDIFF_ERROR(63) ; }
								init_MTRand = (unsigned int)idummy ;
							}
							else { XNDIFF_ERROR(1) ; }
							init_MTRand_userdef = true ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -jmolcif */
					case 'j':
						if ( !strcmp(varg[i], "-jmolcif") )
						{
							/* -jmolcif, write cif-file for molecular viewers like Jmol,PyMol,Rasmol,... */
							jmol_cif_flag = true ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -l */
					case 'l':
						if ( !strcmp(varg[i], "-l") )
						{
							/* -l, generates a logfile of the name [mfname].log */
							/* has been read already at beginning -> reading twice is no problem */
							log_flag = true ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -mem */
					case 'm':
						if ( !strcmp(varg[i], "-mem") )
						{
							/* -mem controls the protocolling of used memory */
							if ( ++i<carg ) 
							{
								if ( ( strlen(varg[i]) == 1 ) && isdigit(varg[i][0]) )
								{
									mem_flag = (unsigned int) strtol( varg[i], NULL, 10) ;
									if ( mem_flag > 2 ) { XNDIFF_ERROR(46) ; }
								}
								else { XNDIFF_ERROR(46) ; }
							}
							else { XNDIFF_ERROR(1) ; }
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -n */
					case 'n':
						if ( !strcmp(varg[i], "-n") )
						{
							/* might has been read already at +par -> don't read again */
							if ( !nlock )
							{
								/* -n, restricts computation on neutron scattering only. 
								(note that using together with the -X flag simultanesously
								gives the default case of X-ray AND neutron scattering)
								multiples of -n flags will be ignored by xlock
								*/
								Xn += -1 ;
								nlock = true ;
							}
						}
						else if ( !strcmp(varg[i], "-nis") )
						{
							/* read the n_i's if par->td==2 directly from the command line or row by row from an ascii file
							   e.g. -nis 0.2 0.4 0.3 0.1 0.05
							   e.g. -nis PPP_180x360_P02_02_ST600_100_N100.nis
							*/
							/* text or numbers */
							if ( ++i<carg )
							{
								if ( is_numeric(varg[i]) )
								{
									while ( i < carg )
									{
										if ( is_numeric(varg[i]) ) { nis.push_back( strtod( varg[i], NULL) ) ; ++i ; }
										else { break ; }
									}
									--i ;
									nis_isdef = true ;
								}
								else
								{
									/* filename or flags */
									if ( *varg[i] != '-' && *varg[i] != '+')
									{
										/* read nis from file */
										read_column_from_ascii_file( string(varg[i]), ' ', 1, nis) ;
										nis_isdef = true ;
									}
									else
									{
										/* nis_isdef = false ; */
										break ;
									}
								}
							}
							else { XNDIFF_ERROR(1) ; }
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* +o -orav -ofo -ow -openmp */
					case 'o':
						if ( !strcmp(varg[i], "+o") )
						{
							/* +o mfname, set user provided master-filename, e.g. +o test */
							/* has been read already at beginning -> reading twice is no problem */
							if ( ++i<carg )
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								strcpy( mfname, varg[i]) ;
							}
							else { XNDIFF_ERROR(1) ; }

							mfname_isdef = true ;
						}
						else if ( !strcmp(varg[i], "-openmp") )
						{
							/* -openmp, use OpenMP parallelization for FFF-computation and orientational average */
							openmp_flag = true ;
						}
						else if ( !strcmp(varg[i], "-orav") )
						{
							/* -orav, compute and export the orientational average of the structure factor |F|^2 */
							orav_flag = true ;
						}
						else if ( !strcmp(varg[i], "-ofo") )
						{
							/* -ofo, change default output folder */
							/* has been read already at beginning -> reading twice is no problem */
							if ( ++i<carg ) 
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								strcpy( ofo, varg[i]) ;
							}
							else { XNDIFF_ERROR(1) ; }
						}
						else if ( !strcmp(varg[i], "-ow") )
						{
							/* has been read already at beginning -> reading twice is no problem */
							overwrite_flag = true ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* +par -pcr */
					case 'p':
						if ( !strcmp(varg[i], "+par") )
						{
							/* +par <parameter-file>, reads parameter-file or
							   +par <d_isl d_osl ...> reads parameters directly from command line

							   single numerical values as well as multi-files for the contrasts are accepted,
							   however there are no checks on consistency (only single only multi)
							*/
							if ( ++i<carg )
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }

								if ( !is_numeric(varg[i]) )
								{
									/* parameter file provided */
									strcpy( par_file, varg[i]) ;
								}
								else
								{
									/* numerical values provided, write parameter file
									   that can be read later again by read_par_file()

									   an existing parameter file in ofo will be overwritten !

									   -thicknesses d_isl d_osl  must be always provided
									   -if -X option only <rho> = < rho_isl rho_osl rho_dm > must be provided
									   -if -n option only <sld> = < sld_isl sld_osl sld_dm > and <icsd> = < icsd_isl icsd_osl icsd_dm > must be provided
									   -if -conc par option also <conc> = < conc_cry conc_dm > must be provided 
									   e.g. X+n (default or -X -n), userdefined conc via -conc par
									       +par < d_isl d_osl <rho> <sld> <icsd> <conc> >
									   e.g. X only (-X), default concentration (default or via -conc def)
									       +par < d_isl d_osl <rho> >
									   e.g. n only (-n), default concentration (default or via -conc def)
									       +par < d_isl d_osl <sld> <icsd> >
									*/

									/* if -X or -n options were not yet read, check for them first
									   if they were already read, Xn will not change due to the lock-variables
									*/
									for ( int j=0; j<carg; ++j)
									{
										if ( !strcmp( varg[j], "-X") && !xlock)
										{
											Xn += 1 ;
											xlock = true ;
											break ;
										}
									}
									for ( int j=0; j<carg; ++j)
									{
										if ( !strcmp( varg[j], "-n") && !nlock)
										{
											Xn += -1 ;
											nlock = true ;
											break ;
										}
									}
									/* Xn is now well-defined */

									/* if -conc option was not yet read, check first */
									for ( int j=0; j<carg; ++j)
									{
										if ( !strcmp(varg[j], "-conc") )
										{
											if ( ++j<carg )
											{
												if ( !strcmp( varg[j], "par") ) { conc_userdef = true ; }
												else if ( !strcmp( varg[j], "def") ) { conc_userdef = false ; }
												else { XNDIFF_ERROR(48) ; }
											}
											else { XNDIFF_ERROR(1) ; }
											break;
										}
									}
									/* conc_userdef is now set */

									/* reset i to varg[i]="+par" */
									--i;
									sprintf( sdummy, "%s.par", mfname) ;
									strcpy( par_file, sdummy) ;
									sprintf( sdummy, "%s%s", ofo, par_file) ;
									if ( log_flag ) { fprintf( logfile, "\t\tWriting par-file %s\n", par_file) ; }
									if ( ( outf = fopen( sdummy, "w")) != NULL)
									{
										fprintf( outf, "# XNDiff par-file generated by XNDiff\n") ;
										write_call_arguments( outf, "# ") ;
										fprintf( outf, "\n# thicknesses of the inner and outer stabilizer layer [nm]\n") ;
										write_par_entry( outf, carg, varg, par_keyword.THICKNESS_INLAY, i, false) ;
										write_par_entry( outf, carg, varg, par_keyword.THICKNESS_OUTLAY, i, false) ;
										if ( Xn >=0 )
										{
											/* X-ray-only related parameters */
											fprintf( outf, "\n# electron densities of the inner and outer stabilizer layer [electrons/nm^3]\n") ;
											write_par_entry( outf, carg, varg, par_keyword.XRAY_RHO_INLAY, i) ;
											write_par_entry( outf, carg, varg, par_keyword.XRAY_RHO_OUTLAY, i) ;
											fprintf( outf, "\n# electron density of dispersion medium [electrons/nm^3]\n") ;
											write_par_entry( outf, carg, varg, par_keyword.XRAY_RHO_DM, i) ;
										}
										if ( Xn <=0 )
										{
											/* neutrons-only related parameters */
											fprintf( outf, "\n# neutron SLD of the inner and outer stabilizer layer [1e-6 A^(-2)]\n") ;
											write_par_entry( outf, carg, varg, par_keyword.NEUT_SLD_INLAY, i) ;
											write_par_entry( outf, carg, varg, par_keyword.NEUT_SLD_OUTLAY, i) ;
											fprintf( outf, "\n# neutron SLD of the dispersion medium [1e-6 A^(-2)]\n") ;
											write_par_entry( outf, carg, varg, par_keyword.NEUT_SLD_DM, i) ;
											fprintf( outf, "\n# neutron ICSD of the inner and outer stabilizer layer [1/cm]\n") ;
											write_par_entry( outf, carg, varg, par_keyword.NEUT_ICSD_INLAY, i) ;
											write_par_entry( outf, carg, varg, par_keyword.NEUT_ICSD_OUTLAY, i) ;
											fprintf( outf, "\n# neutron ICSD of the dispersion medium [1/cm]\n") ;
											write_par_entry( outf, carg, varg, par_keyword.NEUT_ICSD_DM, i) ;
										}
										if ( conc_userdef )
										{
											fprintf( outf, "\n# concentration of dispersion medium and crystals in the nanosuspension [vol%%]\n") ;
											write_par_entry( outf, carg, varg, par_keyword.CONC_CRY, i, false) ;
											write_par_entry( outf, carg, varg, par_keyword.CONC_DM, i, false) ;
										}
									}
									else
									{
										XNDIFF_ERROR(70) ;
									}
									fclose(outf) ;
									if ( log_flag ) { fprintf( logfile, "\t\tdone\n") ; }

									par_num_def = true ;
								}
							}
							else { XNDIFF_ERROR(1) ; }

							par_isdef = true ;
						}
						else if ( !strcmp(varg[i], "-pcr") )
						{
							/* -pcr pcr-file, generates a pcr-file with crystallographic 
							    information for use later in the FullProf Suite */
							if ( ++i<carg )
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								strcpy( pcr_file, varg[i]) ;
								pcr_flag = true ;
							}
							else { XNDIFF_ERROR(1) ; }
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -s -silent -sym -sin */
					case 's':
						if ( !strcmp(varg[i], "-s") )
						{
							/* -s, outputs scattering pattern at each single repeat step -> lots of output ! */
							singlefilesoutput_flag = true ;
						}
						else if ( !strcmp(varg[i], "-sin") )
						{
							/* -sin applies sin weighting factor in the orientational average  */
							weight_flag = true ;
						}
						else if ( !strcmp(varg[i], "-silent") )
						{
							/* has been read already at beginning -> reading twice is no problem */
							/* silent mode, prevent output in terminal
							   already set is 
							   overwrite_flag = true ;
							   set later
							   extension_flag = false ; ( potential -e flag will be ignored )
							   and restrict time_flag <=1
							   time_flag =< 1 ;
							*/
							silent_flag = true ;
						}
						else if ( !strcmp(varg[i], "-sym") )
						{
							/* -sym <ccp4 ccp4configfile <ccp4symopfile> > when using the CCP4,
							    ccp4config must be provided, symmetry operations file e.g. ccp4symop.txt is optional, 
							    if not provided, the CCP4 pdbset routine will use symmetry operations given in syminfo.lib
							   -sym <symopfile> will use symmetry operations given in symop.txt, 
							    if not provided, XNDiff will use the symmetry operations (symop, cenop) given in syminfo.lib
							   -sym <cif> will use the symmetry operations (if available) given in the cif-file (given by +cif),
							    these are in principle equivalent to those in the syminfo.lib (they incorporate both, symop and cenop)
							*/
							if ( (i+1) < carg )
							{
								if ( *varg[i+1] != '-' && *varg[i+1] != '+' )
								{
									++i;
									if ( !strcmp( varg[i], "ccp4") )
									{
										/* use CCP4 to generate symmetry operations  */
										sym_ccp4_mode = true ;

										/* read obligatory configuration file for the CCP4, e.g. ccp4config */
										if ( (i+1) < carg )
										{
											if ( *varg[i+1] != '-' && *varg[i+1] != '+' )
											{
												if ( strlen(varg[++i]) > 1023 ) { XNDIFF_ERROR(6) ; }
												strcpy( ccp4_config_file, varg[i]) ;
											}
											else { XNDIFF_ERROR(66) ; }
										}
										else { XNDIFF_ERROR(66) ; }

										/* use optionally userdefined symmetry operation, e.g. ccp4symop.txt */
										if ( (i+1) < carg )
										{
											if ( *varg[i+1] != '-' && *varg[i+1] != '+' )
											{
												/* read optional symmetry operation files for the CCP4 */
												symop_userdef = true ;

												if ( strlen(varg[++i]) > 1023 ) { XNDIFF_ERROR(6) ; }
												strcpy( symop_file, varg[i]) ;
											}
										}
										/* else { use default symmetry operations for pdbset in CCP4 by HM or Hall symbol with syminfo.lib } */
									}
									else if ( !strcmp(varg[i], "cif") )
									{
										/* read symmetry operations from cif-file if available */
										sym_cif_mode = true ;
									}
									else if ( !strcmp(varg[i], "pdb") )
									{
										/* pdb pdb-file, allows bypassing the CCP4/syminfo.lib/etc by directly 
										   providing the atom coordinates in a pdb-file, e.g. pdb betaPPP.pdb. 
										*/
										if ( (i+1) < carg )
										{
											if ( *varg[i+1] != '-' && *varg[i+1] != '+' )
											{
												if ( strlen(varg[++i]) > 1023 ) { XNDIFF_ERROR(6) ; }
												strcpy( pdb_file, varg[i]) ;
												sym_pdb_mode = true ;
											}
											else { XNDIFF_ERROR(67) ;}
										}
										else { XNDIFF_ERROR(67) ; }
									}
									else
									{
										/* use user-defined symmetry operations, e.g. from symop.txt */
										symop_userdef = true ;
										if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
										strcpy( symop_file, varg[i]) ;
									}
								}
								/* else { use symop by HM or Hall symbol with syminfo.lib } */
							}
							/* else { use symop by HM or Hall symbol with syminfo.lib } */
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -time -test_XnEquiv */
					case 't':
						if ( !strcmp(varg[i], "-time") )
						{
							/* -time controls the protocolling of elapsed time
							
							    time_flag = 0 -> protocol starting and terminating time of program in logfile
							    time_flag = 1 -> protocol additionally writing time in files and
							                     protocol additionally the time consumption in logfile for calls of :
							                     compute_eQvj()
							                     calculate_orientational_averaged_structure_amplitude()
							                     calculate_orientational_averaged_structure_factor()
							                     compute_structure_amplitude()
							                     compute_f0table_and_add_f0ind_to_atom()
							                     read_structure_amplitude()
							                     write_structure_amplitude()
							    time_flag = 2 -> print additionally progressbar on stdout for call of :
							    (default)        compute_structure_amplitude()
							                     stackcpp()
							                     print additionally time for a new stack in stackcpp() in logfile
							*/
							if ( ++i<carg ) 
							{
								if ( ( strlen(varg[i]) == 1 ) && isdigit(varg[i][0]) )
								{
									time_flag = (unsigned int) strtol( varg[i], NULL, 10) ;
									if ( time_flag > 2 ) { XNDIFF_ERROR(47) ; }
								}
								else { XNDIFF_ERROR(47) ; }
							}
							else { XNDIFF_ERROR(1) ; }
						}
						else if ( !strcmp( varg[i], "-test_XnEquiv") )
						{
							test_XnEquiv_flag = true ;
						}
						else if ( !strcmp( varg[i], "-test_I_at_s_equal_zero") )
						{
							test_I_at_s_equal_zero_flag = true ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -X */
					case 'X':
						if ( !strcmp(varg[i], "-X"))
						{
							/* might has been read already at +par -> don't read again */
							if ( !xlock)
							{
								/* -X, restricts computation on X-ray scattering only. 
								(note that using together with the -n flag simultanesously
								gives the default case of X-ray AND neutron scattering)
								multiples of -X flags will be ignored by xlock
								*/
								Xn += 1 ;
								xlock = true ;
							}
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					/* -z */
					case 'z':
						if ( !strcmp(varg[i], "-z") )
						{
							/* -z imsteps, outputs the added up scattering patterns at intermediate steps.
							   ( at every mod(imsteps,repeats) ), e.g. -z 100 */
							if ( ++i<carg )
							{
								if ( strlen(varg[i]) > 1023 ) { XNDIFF_ERROR(6) ; }
								imsteps = (unsigned int) strtol (varg[i], NULL, 10) ;
								if ( !(imsteps>0) ) { XNDIFF_ERROR(20) ; }
							}
							else { XNDIFF_ERROR(1) ; }
							intermediatesteps_flag = true ;
						}
						else { XNDIFF_ERROR(3) ; }
						break ;
					default:
						/* unknown options */
						XNDIFF_ERROR(3) ;
						break ;
				}
			}
			else
			{
				/* unknown input */
				XNDIFF_ERROR(3) ;
				break ;
			}
		}
		if ( log_flag ) { fprintf( logfile, "done\n\n") ; fflush (logfile) ; }

		/* check if all essential parameters (+...) are defined and look for incompatibilities in the input (options), incompatibilities check not yet implemented !!! */
		if ( !av_isdef || !cif_isdef || !func_isdef || !mfname_isdef || !par_isdef ) { XNDIFF_ERROR(78) ; }

		/* check for consistency with cis and nis when par->td2 = 2 */
		if ( par->td2 == 2 )
		{
			// either -nis or -cis must be provided 
			if ( cis_isdef == nis_isdef ) { XNDIFF_ERROR(12) ; }

			// check length
			if ( cis_isdef == true && cis.size() != par->nsp ) { XNDIFF_ERROR(41) ; }
			if ( nis_isdef == true && nis.size() != par->nsp ) { XNDIFF_ERROR(41) ; }
		}

		if ( silent_flag ) 
		{ 
			/* overwrite_flag = true ; ... already set */
			extension_flag = false ; 
			if ( time_flag > 1 ) { time_flag = 1 ; }
		}

		int err ;

		/* call stackcpp or other functions (not implemented) */
		switch (func_num)
		{
			case 0:
				/* call function stackcpp and if errors occur print them via XNDIFF_ERROR(err) */ 
				if ((err = stackcpp()) != 0)
				{
					XNDIFF_ERROR(err) ;
				}
				break ;
		}
		/* if logging was activated close now the log file */
		if ( log_flag ) { fclose(logfile) ; }
	}


	/* index i will we returned via reference */
	void write_par_entry(FILE *outf, int carg, char** varg, const char *keyword, int &i, bool multi=true)
	{
		if ( ++i < carg )
		{
			if ( is_numeric(varg[i]) )
			{
				/* single numeric value  */
				fprintf( outf, "%s %s\n", keyword, varg[i]) ;
			}
			else
			{
				if ( (i+1) < carg )
				{
					if ( is_unsigned_int(varg[i+1]) && multi==true )
					{ 
						/* multi-file and column number */
						fprintf( outf, "%s %s %s\n", keyword, varg[i], varg[i+1]) ;
						++i;
					}
					else
					{
						XNDIFF_ERROR(71) ;
					}
				}
				else { XNDIFF_ERROR(1) ; }
			}
		}
		else { XNDIFF_ERROR(1) ; }
	}


	/* check if already a logfile with the same mfname exists */
	void overwrite_protection()
	{
		if ( overwrite_flag == false )
		{
			FILE *inpf ;
			char sdummy[1024] ;
			char decision[1024] ;

			sprintf( sdummy, "%s%s.log", ofo, mfname) ;
			if ( ( inpf = fopen( sdummy, "r")) != NULL) 
			{
				fprintf( stdout, "Data with the same mfname %s already exists. Overwrite (y/n): ", mfname) ;
				while (true)
				{
					scanf( "%s", decision) ;
					if ( !strcmp( decision, "y") )
					{
						fprintf( stdout, "Data will be overwritten.\n");
						overwrite_flag = true ;
						break ;
					}
					else if ( !strcmp( decision, "n") )
					{
						XNDIFF_ERROR (45) ;
					}
					else
					{
						fprintf( stdout, "Type 'y' or 'n': ") ;
					}
				}
			}
		}
	}


	/* This function returns a uniform distributed random number in the range (0,1) */
	/* Call with inp < 0 for initialization, and after that with any number         */
	double i_rand ( long int *inp )
	{
		int i ;
		long int k ;
		static long int i1 = 0 ;
		static long int i2[R6] ;
		double ddummy1 ;

		if ( *inp <= 0 || !i1 )
		{
			if ( *inp > -1 )
				*inp = 1 ;
			else
				*inp = - ( *inp ) ;
			for ( i=R6+9; i>=0; --i )
			{
				k = *inp / R4 ;
				*inp = R1* ( *inp-k*R4 ) - k*R5 ;
				if ( *inp < 0 )
					*inp += R2 ;
				if ( i < R6 )
					i2[i] = *inp ;
			}
			i1 = i2[0] ;
		}
		k = *inp / R4 ;
		*inp = R1 * ( *inp-k*R4 ) - k*R5 ;
		if ( *inp < 0 )
			*inp += R2 ;
		i = i1 / R7 ;
		i1 = i2[i] ;
		i2[i] = *inp ;
		if ( ( ddummy1 = R3*i1 ) > R9 )
			return ( R9 ) ;
		else
			return ( ddummy1 ) ;
	}

	void usage (char *name)
	{
		/* OBSOLETE */
		fprintf( stdout, "\n") ;
		fprintf( stdout, "usage: %s [-l] [-h] [-info <function name>] [-o <master file name>] <function {par} {par} ...>\n", name) ;
		fprintf( stdout, "  -l, +l:    turn logging on (-l) or off (+l)\n") ;
		fprintf( stdout, "  available functions:\n") ; 
		fprintf( stdout, "    1) stack(n11,n12,n21,n22,n31,n32,n1,n2,D1,D2,rs1,rs2,np,td,ns,nr,h,k,l,hs,ks,ls,d11,d12,\n") ;
		fprintf( stdout, "             d21,d22,r11,r12,r21,r22,r31,r32,s)\n") ;
		fprintf( stdout, "         n1: Mean particle size along a1 and its standard deviation [nm] (eg: 10.0 0.5)\n") ;
		fprintf( stdout, "         n2: Mean particle size along a2 and its standard deviation [nm] (eg: 10.0 0.5)\n") ;
		fprintf( stdout, "         n3: Mean particle size along a3 and its standard deviation [nm] (eg: 10.0 0.5)\n") ;
		fprintf( stdout, "          n: Number of particles per stack and its standard deviation (eg: 10.0 2.0)\n") ;
		fprintf( stdout, "          D: Mean repeat distance along Gs and its standard deviation [nm] (eg: 20.0 5.0\n") ;
		fprintf( stdout, "         rs: S-range for calculation [1/nm] (z.B. 0.01 0.5)\n") ;
		fprintf( stdout, "         np: Number of data points in s-range\n") ;
		fprintf( stdout, "         td: Type of distribution (0: norm., 1: log. norm.)\n") ;
		fprintf( stdout, "         ns: Number of stacks\n") ;
		fprintf( stdout, "         nr: Number of repeats with averaging\n") ;
		fprintf( stdout, "          G: (h,k,l) values representing the direction of the scattering vektors (z.B. 0 0 1)\n") ;
		fprintf( stdout, "         Gs: (hs,ks,ls) values of stacking orientation of the crystals (z.B. 0.0 0.0 1.0)\n") ;
		fprintf( stdout, "         d1: Mean displacement and std. dev. of crystals in stack perpendicular to Gs (\"d1\")\n") ;
		fprintf( stdout, "         d2: Mean displacement and std. dev. of crystals in stack perpendicular to Gs and \"d1\"\n") ;
		fprintf( stdout, "         r1: Mean rotation and std. dev. of crystals, polar angle (vs Gs) [degree]\n") ;
		fprintf( stdout, "         r2: Mean rotation and std. dev. of crystals, azimuthal angle (vs \"d1\") [degree]\n") ;
		fprintf( stdout, "         r3: Mean rotation and std. dev. of crystals in stack on Gs (vs \"d1\") [degree]\n") ;
		fprintf( stdout, "         av: Average stacks orientations: <mode> <polar angle> <azimuthal angle>\n") ;
		fprintf( stdout, "                      mode:            0: No averaging\n") ;
		fprintf( stdout, "                                       1: Add all diffraction patterns with G(hkl)=2*pi*s\n") ;
		fprintf( stdout, "                                       2: Add all diffraction patterns of specified angles\n") ;
		fprintf( stdout, "                                       3: Add diffraction patterns of mode 1 and 2 both\n") ;
		fprintf( stdout, "                                      10: Special mode 0: Output structure factors instead of S\n") ;
		fprintf( stdout, "                                      11: Special mode 1: As mode 2 but do not perform dispersion medium correction\n") ;
		fprintf( stdout, "                      polar angle:     Range and step of polar angle vs G [degree]\n") ;
		fprintf( stdout, "                                          (eg: for whole range 0 180 0.1)\n") ;
		fprintf( stdout, "                      azimuthal angle: Range and step of azimuthal angle vs. \"d1\" [degree]\n") ;
		fprintf( stdout, "                                          (eg: for whole range 0 360 0.1)\n") ;
		fprintf( stdout, "          s: File name of structure data file\n") ;
		fprintf( stdout, "             File format: <space group symbol> a b c [nm] alpha beta gamma [degree]\n") ;
		fprintf( stdout, "                          <electron density of dispersion medium in electrons per nm^3>\n") ;
		fprintf( stdout, "                          <total irradiated sample volume in nm^3>\n") ;
		fprintf( stdout, "                          <atom1 label> <type> x y z B0 b11 b22 b33 b12 b13 b23\n") ;
		fprintf( stdout, "                          <atom2 label> <type> x y z B0 b11 b22 b33 b12 b13 b23\n") ;
		fprintf( stdout, "                          ...\n") ;
		fprintf( stdout, "\n") ;
	}


	void XNDIFF_ERROR(int ErrorCode)
	{
		const char *MP_ERR[] =
		{
			"No error",                                                /* 000 */
			"Missing entry after command line option",                 /* 001 */
			"Syntax error",                                            /* 002 */
			"Unknown command line option",                             /* 003 */
			"Unknown function name for func_XPPSA or func_bt",         /* 004 */
			"stackcpp - Number of input parameters does not match",    /* 005 */
			"String constant of command line parameter too long (max 1023 signs)",      /* 006 */
			"Stack-function: Unknown atom label",                      /* 007 */
			"stackcpp-function: Could not open par file",      	   /* 008 */
			"Stack-function: Syntax error in structure data file",     /* 009 */
			"Stack-function: Vector G is zero",                        /* 010 */
			"Stack-function: Vector Gs is zero",                       /* 011 */
			"For td=2 either -cis or -nis must be provided",           /* 012 */
			"Unable to open log file in cwd ",                         /* 013 */
			"Stack-function: Unable to open output file",              /* 014 */
			"Stack-function: Overlap of crystals occurs continuously", /* 015 */
			"Stack-function: Unknown space group symbol",              /* 016 */
			"stackcpp-function: Average mode flag must be between 0 and MAX_AV_MODE",   /* 017 */
			"stackcpp-function: Unable to write to the *.distr or *.pol_azi_ang file", /* 018 */
			"stackcpp-function: Unable to open one of the output files", /* 019 */
			"stackcpp-function: The number of intermediate steps for data output should be positive", /* 020 */
			"stackcpp-function: Could not open cif-dictionary", /* 021 */
			"stackcpp-function: Could not open cif-file", /* 022 */
			"stackcpp-function: Error while running the bash-script", /* 023 */
			"stackcpp-function: Could not open symmetry operations file", /* 024 */
			"stackcpp-function: Could not open par file", /* 025 */
			"stackcpp-function: Could not open CCP4-configuration file", /* 026 */
			"stackcpp-function: Could not open bash-script", /* 027 */
			"stackcpp-function: Could not open cif_pdb_file", /* 028 */
			"stackcpp-function: Could not open pdb_file", /* 029 */
			"stackcpp-function: Could not find an entry in the cif-file", /* 030 */
			"stackcpp-function: Could not open 'DeBe_NeutronNews.dat'", /* 031 */
			"stackcpp-function: Could not open 'atomweights.dat'", /* 032 */
			"stackcpp-function: Could not open one of the Fr files", /* 033 */
			"stackcpp-function: Could not open pcr-file", /* 034 */
			"stackcpp-function: Dimension mismatch für multiple parameters in par-file", /* 035 */
			"stackcpp-function: non-numeric parameter found.", /* 036 */
			"stackcpp-function: Difference in array lengths found for rho/sld", /* 037 */
			"stackcpp-function: Could not open rho/sld_multi file", /* 038 */
			"stackcpp-function: Missing column in rho/sld_multi file", /* 039 */
			"stackcpp-function: B_name=NULL or/and B=NULL in write_header_in_file()", /* 040 */
			"Length of cis or nis must be the same as nsp for td2=2", /* 041 */
			"stackcpp-function: B_name=NULL or/and B=NULL in write_Y()", /* 042 */
			"stackcpp-function: Could not create input folder (ifo)", /* 043 */
			"stackcpp-function: Could not create output folder (ofo)", /* 044 */
			"stackcpp-function: Aborted by user.", /* 045 */
			"stackcpp-function: Parameter for -mem option must be an integer in the range [0,2].", /* 046 */
			"stackcpp-function: Parameter for -time option must be an integer in the range [0,2].", /* 047 */
			"stackcpp-function: Unknown parameter for option -conc. Use either par or def.", /* 048 */
			"stackcpp-function: Missing entry THICKNESS_INLAY in par-file.", /* 049 */
			"stackcpp-function: Missing entry THICKNESS_OUTLAY in par-file.", /* 050 */
			"stackcpp-function: Missing entry XRAY_RHO_INLAY in par-file.", /* 051 */
			"stackcpp-function: Missing entry XRAY_RHO_OUTLAY in par-file.", /* 052 */
			"stackcpp-function: Missing entry XRAY_RHO_DM in par-file.", /* 053 */
			"stackcpp-function: Missing entry NEUT_SLD_INLAY in par-file.", /* 054 */
			"stackcpp-function: Missing entry NEUT_SLD_OUTLAY in par-file.", /* 055 */
			"stackcpp-function: Missing entry NEUT_SLD_DM in par-file.", /* 056 */
			"stackcpp-function: Missing entry ICSD_SLD_INLAY in par-file.", /* 057 */
			"stackcpp-function: Missing entry ICSD_SLD_OUTLAY in par-file.", /* 058 */
			"stackcpp-function: Missing entry ICSD_SLD_DM in par-file.", /* 059 */
			"stackcpp-function: Missing entry CONC_CRY in par-file.", /* 060 */
			"stackcpp-function: Missing entry CONC_DM in par-file.", /* 061 */
			"stackcpp-function: Argument for -init_n_rand flag must be a negative integer.", /* 062 */
			"stackcpp-function: Argument for -init_d_rand and -init_MTRand flag must be an unsigned integer.", /* 063 */
			"Could not open-write Jmol-cif-file", /* 064 */
			"noa is not a multiple of the noa in the cif-file. Cannot derive labels/types/bonds.", /* 065 */
			"-sym ccp4 <ccp4config-file> Missing configuration file for CCP4 <ccp4config-file>.", /* 066 */
			"-sym pdb <pdb-file> Missing pdb-file <pdb-file>.", /* 067 */
			"-sym cif: Missing entry _symmetry_equiv_pos_as_xyz in cif-file.", /* 068 */
			"-sym : space group information in cif-file missing or corresponding entry not found in syminfo.lib.", /* 069 */
			"+par : could not write par-file.", /* 070 */
			"+par : Missing or wrong parameter input.", /* 071 */
			"get_current_memory(): error in writing the bash-script memory_<OS>_<PID>.sh.", /* 072 */
			"get_current_memory(): cannot change the executable rights of the bash-script memory_<OS>_<PID>.sh.", /* 073 */
			"get_current_memory(): error in reading the output of bash-script memory_<OS>_<PID>.sh.", /* 074 */
			"Could not open syminfo.lib.", /* 075 */
			"Error in formula interpreter, incomplete expression.", /* 076 */
			"+par : Too much arguments.", /* 077 */
			"Essential paramters are missing in the input arguments list.", /* 078 */
			"Error in read_single_numeric_data_from_cfe(). The data type is not numeric.", /* 079 */
			"Error in read_single_text_data_from_cfe(). The data type is not text.", /* 080 */
			"_atom_site_fract coordinates >10 will lead to PDB-format errors when writing pdb-file.", /* 081 */
			"test_XnEquiv works only for Xn=0 (X-ray AND neutron scattering).", /* 082 */
			"test_XnEquiv should not be applied for multiple rho/sld/icsd's.", /* 083 */
			"Expected \"symop\" clause in symop-file.", /* 084 */
			"Type of distribution for td1 must be either 0,1.", /* 085 */
			"Type of distribution for td2 must be either 0,1,2.", /* 086 */
			"Stack mode must be between 0-4." /* 087 */
		} ;

		/* before exit, write and close the logfile */
		if ( log_flag && logfile != NULL )
		{
			fflush(logfile);
			fclose(logfile);
		}

		if (ErrorCode == 0)
		{
			return ;
		}

		fprintf( stderr, "XNDIFF - ERROR %d: %s\n", ErrorCode, MP_ERR[ErrorCode]) ;
		exit(1) ;
	}


	private:
	/* This function returns a normal (type = 0), log normal (type = 1)           	*/
	/* distributed random number with mean m and standard deviation s             	*/
	/* for initialization use type = -1                                           	*/
	/* http://de.wikipedia.org/wiki/Logarithmische_Normalverteilung 		*/
	/* http://de.wikipedia.org/wiki/Normalverteilung 				*/
	double n_rand ( int type, double m, double s, int const_init = 0 )
	{
		struct timeval tv ;
		static long int init ;
		static int st1 = 0 ;
		static double st2 ;
		double f, r, d1, d2 ;
	
		if ( type == -1 )
		{
			if ( const_init == 0 )
			{
				gettimeofday ( &tv, NULL ) ;
				init = - ( tv.tv_usec & 32767 ) ;
				if ( log_flag ) { fprintf( logfile, "\tRandom Number Generator Initialization with random seed: %ld\n", init ) ; }
			}
			else
			{
				init = const_init ;
				if ( log_flag ) { fprintf( logfile, "\tRandom number Generator Initialization with user-defined seed: %ld\n", init ) ; }
			}
			i_rand ( &init ) ;
			return ( 0.0 ) ;
		}
		switch ( type )
		{
			case 0:
				if ( st1 == 0 )
				{
					do
					{
						d1 = 2.0 * i_rand ( &init ) - 1.0 ;
						d2 = 2.0 * i_rand ( &init ) - 1.0 ;
						r = SQUARE ( d1 ) +SQUARE ( d2 ) ;
					}
					while ( r >= 1.0 || r == 0.0 ) ;
					f = sqrt ( -2.0*log ( r ) /r ) ;
					st2 = d1*f ;
					st1 = 1 ;
					return ( d2*f*s+m ) ;
				}
				else
				{
					st1 = 0 ;
					return ( st2*s+m ) ;
				}
				break ;
			case 1:
				if ( st1 == 0 )
				{
					do
					{
						d1 = 2.0 * i_rand ( &init ) - 1.0 ;
						d2 = 2.0 * i_rand ( &init ) - 1.0 ;
						r = SQUARE ( d1 ) +SQUARE ( d2 ) ;
					}
					while ( r >= 1.0 || r == 0.0 ) ;
					f = sqrt ( -2.0*log ( r ) /r ) ;
					st2 = d1*f ;
					st1 = 1 ;
					return ( exp ( d2*f*s+m ) ) ;
				}
				else
				{
					st1 = 0 ;
					return ( exp ( st2*s+m ) ) ;
				}
				break ;
		}
		return ( 0.0 ) ;
	}


	/* This function returns the non-dispersive part of the atomic scattering factor of "atom" at s=2sin(theta)/lambda [1/nm] 
	   http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/
	   http://ftp.esrf.eu/pub/scisoft/xop2.3/DabaxFiles/f0_InterTables.dat */
	double get_atomic_scattering_factor ( char *aatom, double s )
	{
		/* check on isotopic labels and remove them only the ionic composition is relevant for X-ray scattering */
		/* e.g. 2H -> H or 2H1- -> H1- */		
		int i;
		int strlenaatom = strlen(aatom);
		char *strptr1, *strptr2;

		strptr1 = aatom ;
		strptr2 = strptr1 + strlenaatom; /* points to terminating  0 */
		while ( isdigit(*strptr1) ) 
			++strptr1;

		i = (int)(strptr2-strptr1) ;
		char* entry = (char *) calloc( i+1, sizeof(char) ) ;
		memmove (entry, strptr1, i) ;

		double f ;
		static struct atomic_scattering_factors
		{
			const char *atom ;
			double par[9] ; /* a1  a2  a3  a4  c  b1  b2  b3  b4 */
		} asf[] =
		{ 	{"H"    , {0.4930020,0.3229120,0.1401910,4.0810000e-02,3.0380001E-03,10.51090,26.12570,3.142360,57.79970}},
			{"H."   , {0.4899180,0.2620030,0.1967670,4.9879000e-02,1.3050000E-03,20.65930,7.740390,49.55190,2.201590}},
			{"H1-"  , {0.8976610,0.5656160,0.4158150,0.1169730,2.3890000e-03,53.13680,15.18700,186.5760,3.567090}},
			{"He"   , {0.8734000,0.6309000,0.3112000,0.1780000,6.3999998e-03,9.103700,3.356800,22.92760,0.9821000}},
			{"Li"   , {1.128200,0.7508000,0.6175000,0.4653000,3.7700001e-02,3.954600,1.052400,85.39050,168.2610}},
			{"Li1+" , {0.6968000,0.7888000,0.3414000,0.1563000,1.6700000e-02,4.623700,1.955700,0.6316000,10.09530}},
			{"Be"   , {1.591900,1.127800,0.5391000,0.7029000,3.8500000e-02,43.64270,1.862300,103.4830,0.5420000}},
			{"Be2+" , {6.260300,0.8849000,0.7993000,0.1647000,-6.109200,2.7000001e-03,0.8313000,2.275800,5.114600}},
			{"B"    , {2.054500,1.332600,1.097900,0.7068000,-0.1932000,23.21850,1.021000,60.34980,0.1403000}},
			{"C"    , {2.310000,1.020000,1.588600,0.8650000,0.2156000,20.84390,10.20750,0.5687000,51.65120}},
			{"C."   , {2.260690,1.561650,1.050750,0.8392590,0.2869770,22.69070,0.6566650,9.756180,55.59490}},
			{"N"    , {12.21260,3.132200,2.012500,1.166300,-11.52900,5.7000001e-03,9.893300,28.99750,0.5826000}},
			{"O"    , {3.048500,2.286800,1.546300,0.8670000,0.2508000,13.27710,5.701100,0.3239000,32.90890}},
			{"O1-"  , {4.191600,1.639690,1.526730,-20.30700,21.94120,12.85730,4.172360,47.01790,-1.4040000e-02}},
			{"O2-." , {4.758000,3.637000,0.000000,0.000000,1.594000,7.831000,30.05000,0.000000,0.000000}},
			{"F"    , {3.539200,2.641200,1.517000,1.024300,0.2776000,10.28250,4.294400,0.2615000,26.14760}},
			{"F1-"  , {3.632200,3.510570,1.260640,0.9407060,0.6533960,5.277560,14.73530,0.4422580,47.34370}},
			{"Ne"   , {3.955300,3.112500,1.454600,1.125100,0.3515000,8.404200,3.426200,0.2306000,21.71840}},
			{"Na"   , {4.762600,3.173600,1.267400,1.112800,0.6760000,3.285000,8.842200,0.3136000,129.4240}},
			{"Na1+" , {3.256500,3.936200,1.399800,1.003200,0.4040000,2.667100,6.115300,0.2001000,14.03900}},
			{"Mg"   , {5.420400,2.173500,1.226900,2.307300,0.8584000,2.827500,79.26110,0.3808000,7.193700}},
			{"Mg2+" , {3.498800,3.837800,1.328400,0.8497000,0.4853000,2.167600,4.754200,0.1850000,10.14110}},
			{"Al"   , {6.420200,1.900200,1.593600,1.964600,1.115100,3.038700,0.7426000,31.54720,85.08860}},
			{"Al3+" , {4.174480,3.387600,1.202960,0.5281370,0.7067860,1.938160,4.145530,0.2287530,8.285240}},
			{"Si"   , {6.291500,3.035300,1.989100,1.541000,1.140700,2.438600,32.33370,0.6785000,81.69370}},
			{"Si."  , {5.662690,3.071640,2.624460,1.393200,1.247070,2.665200,38.66340,0.9169460,93.54580}},
			{"Si4+" , {4.439180,3.203450,1.194530,0.4165300,0.7462970,1.641670,3.437570,0.2149000,6.653650}},
			{"S"    , {6.905300,5.203400,1.437900,1.586300,0.8669000,1.467900,22.21510,0.2536000,56.17200}},
			{"P"    , {6.434500,4.179100,1.780000,1.490800,1.114900,1.906700,27.15700,0.5260000,68.16450}},
			{"Cl"   , {11.46040,7.196400,6.255600,1.645500,-9.557400,1.0400000E-02,1.166200,18.51940,47.77840}},
			{"Cl1-" , {18.29150,7.208400,6.533700,2.338600,-16.37800,6.6000000E-03,1.171700,19.54240,60.44860}},
			{"Ar"   , {7.484500,6.772300,0.6539000,1.644200,1.444500,0.9072000,14.84070,43.89830,33.39290}},
			{"K"    , {8.218600,7.439800,1.051900,0.8659000,1.422800,12.79490,0.7748000,213.1870,41.68410}},
			{"K1+"  , {7.957800,7.491700,6.359000,1.191500,-4.997800,12.63310,0.7674000,-2.0000001E-03,31.91280}},
			{"Ca"   , {8.626600,7.387300,1.589900,1.021100,1.375100,10.44210,0.6599000,85.74840,178.4370}},
			{"Ca2+" , {15.63480,7.951800,8.437200,0.8537000,-14.87500,-7.4000000E-03,0.6089000,10.31160,25.99050}},
			{"Sc"   , {9.189000,7.367900,1.640900,1.468000,1.332900,9.021300,0.5729000,136.1080,51.35310}},
			{"Sc3+" , {13.40080,8.027300,1.659430,1.579360,-6.666700,0.2985400,7.962900,-0.2860400,16.06620}},
			{"Ti"   , {9.759500,7.355800,1.699100,1.902100,1.280700,7.850800,0.5000000,35.63380,116.1050}},
			{"Ti2+" , {9.114230,7.621740,2.279300,8.7898999E-02,0.8971550,7.524300,0.4575850,19.53610,61.65580}},
			{"Ti3+" , {17.73440,8.738160,5.256910,1.921340,-14.65200,0.2206100,7.047160,-0.1576200,15.97680}},
			{"Ti4+" , {19.51140,8.234730,2.013410,1.520800,-13.28000,0.1788470,6.670180,-0.2926300,12.94640}},
			{"V"    , {10.29710,7.351100,2.070300,2.057100,1.219900,6.865700,0.4385000,26.89380,102.4780}},
			{"V2+"  , {10.10600,7.354100,2.288400,2.2299999E-02,1.229800,6.881800,0.4409000,20.30040,115.1220}},
			{"V3+"  , {9.431410,7.741900,2.153430,1.6865000E-02,0.6565650,6.395350,0.3833490,15.19080,63.96900}},
			{"V5+"  , {15.68870,8.142080,2.030810,-9.576000,1.714300,0.6790030,5.401350,9.972780,0.9404640}},
			{"Cr"   , {10.64060,7.353700,3.324000,1.492200,1.183200,6.103800,0.3920000,20.26260,98.73990}},
			{"Cr2+" , {9.540340,7.750900,3.582740,0.5091070,0.6168980,5.660780,0.3442610,13.30750,32.42240}},
			{"Cr3+" , {9.680900,7.811360,2.876030,0.1135750,0.5182750,5.594630,0.3343930,12.82880,32.87610}},
			{"Mn"   , {11.28190,7.357300,3.019300,2.244100,1.089600,5.340900,0.3432000,17.86740,83.75430}},
			{"Mn2+" , {10.80610,7.362000,3.526800,0.2184000,1.087400,5.279600,0.3435000,14.34300,41.32350}},
			{"Mn3+" , {9.845210,7.871940,3.565310,0.3236130,0.3939740,4.917970,0.2943930,10.81710,24.12810}},
			{"Mn4+" , {9.962530,7.970570,2.760670,5.4446999E-02,0.2518770,4.848500,0.2833030,10.48520,27.57300}},
			{"Fe"   , {11.76950,7.357300,3.522200,2.304500,1.036900,4.761100,0.3072000,15.35350,76.88050}},
			{"Fe+2" , {11.04240,7.374000,4.134600,0.4399000,1.009700,4.653800,0.3053000,12.05460,31.28090}},
			{"Fe3+" , {11.17640,7.386300,3.394800,7.2400004E-02,0.9707000,4.614700,0.3005000,11.67290,38.55660}},
			{"Co"   , {12.28410,7.340900,4.003400,2.348800,1.011800,4.279100,0.2784000,13.53590,71.16920}},
			{"Co2+" , {11.22960,7.388300,4.739300,0.7108000,0.9324000,4.123100,0.2726000,10.24430,25.64660}},
			{"Co3+" , {10.33800,7.881730,4.767950,0.7255910,0.2866670,3.909690,0.2386680,8.355830,18.34910}},
			{"Ni"   , {12.83760,7.292000,4.443800,2.380000,1.034100,3.878500,0.2565000,12.17630,66.34210}},
			{"Ni2+" , {11.41660,7.400500,5.344200,0.9773000,0.8614000,3.676600,0.2449000,8.873000,22.16260}},
			{"Ni3+" , {10.78060,7.758680,5.227460,0.8471140,0.3860440,3.547700,0.2231400,7.644680,16.96730}},
			{"Cu"   , {13.33800,7.167600,5.615800,1.673500,1.191000,3.582800,0.2470000,11.39660,64.81260}},
			{"Cu1+" , {11.94750,7.357300,6.245500,1.557800,0.8900000,3.366900,0.2274000,8.662500,25.84870}},
			{"Cu2+" , {11.81680,7.111810,5.781350,1.145230,1.144310,3.374840,0.2440780,7.987600,19.89700}},
			{"Zn"   , {14.07430,7.031800,5.165200,2.410000,1.304100,3.265500,0.2333000,10.31630,58.70970}},
			{"Zn2+" , {11.97190,7.386200,6.466800,1.394000,0.7807000,2.994600,0.2031000,7.082600,18.09950}},
			{"Ga"   , {15.23540,6.700600,4.359100,2.962300,1.718900,3.066900,0.2412000,10.78050,61.41350}},
			{"Ga3+" , {12.69200,6.698830,6.066920,1.006600,1.535450,2.812620,0.2278900,6.364410,14.41220}},
			{"Ge"   , {16.08160,6.374700,3.706800,3.683000,2.131300,2.850900,0.2516000,11.44680,54.76250}},
			{"Ge4+" , {12.91720,6.700030,6.067910,0.8590410,1.455720,2.537180,0.2058550,5.479130,11.60300}},
			{"As"   , {16.67230,6.070100,3.431300,4.277900,2.531000,2.634500,0.2647000,12.94790,47.79720}},
			{"Se"   , {17.00060,5.819600,3.973100,4.354300,2.840900,2.409800,0.2726000,15.23720,43.81630}},
			{"Br"   , {17.17890,5.235800,5.637700,3.985100,2.955700,2.172300,16.57960,0.2609000,41.43280}},
			{"Br1-" , {17.17180,6.333800,5.575400,3.727200,3.177600,2.205900,19.33450,0.2871000,58.15350}},
			{"Kr"   , {17.35550,6.728600,5.549300,3.537500,2.825000,1.938400,16.56230,0.2261000,39.39720}},
			{"Rb"   , {17.17840,9.643500,5.139900,1.529200,3.487300,1.788800,17.31510,0.2748000,164.9340}},
			{"Rb1+" , {17.58160,7.659800,5.898100,2.781700,2.078200,1.713900,14.79570,0.1603000,31.20870}},
			{"Sr"   , {17.56630,9.818400,5.422000,2.669400,2.506400,1.556400,14.09880,0.1664000,132.3760}},
			{"Sr2+" , {18.08740,8.137300,2.565400,-34.19300,41.40250,1.490700,12.69630,24.56510,-1.3800000E-02}},
			{"Y"    , {17.77600,10.29460,5.726290,3.265880,1.912130,1.402900,12.80060,0.1255990,104.3540}},
			{"Y3+"  , {17.92680,9.153100,1.767950,-33.10800,40.26020,1.354170,11.21450,22.65990,-1.3190000E-02}},
			{"Zr"   , {17.87650,10.94800,5.417320,3.657210,2.069290,1.276180,11.91600,0.1176220,87.66270}},
			{"Zr4+" , {18.16680,10.05620,1.011180,-2.647900,9.414540,1.214800,10.14830,21.60540,-0.1027600}},
			{"Nb"   , {17.61420,12.01440,4.041830,3.533460,3.755910,1.188650,11.76600,0.2047850,69.79570}},
			{"Nb3+" , {19.88120,18.06530,11.01770,1.947150,-12.91200,1.9175000E-02,1.133050,10.16210,28.33890}},
			{"Nb5+" , {17.91630,13.34170,10.79900,0.3379050,-6.393400,1.124460,2.8781001E-02,9.282060,25.72280}},
			{"Mo"   , {3.702500,17.23560,12.88760,3.742900,4.387500,0.2772000,1.095800,11.00400,61.65840}},
			{"Mo3+" , {21.16640,18.20170,11.74230,2.309510,-14.42100,1.4734000E-02,1.030310,9.536590,26.63070}},
			{"Mo5+" , {21.01490,18.09920,11.46320,0.7406250,-14.31600,1.4345000E-02,1.022380,8.788090,23.34520}},
			{"Mo6+" , {17.88710,11.17500,6.578910,0.0000000E+00,0.3449410,1.036490,8.480610,5.8881000E-02,0.0000000E+00}},
			{"Tc"   , {19.13010,11.09480,4.649010,2.712630,5.404280,0.8641320,8.144870,21.57070,86.84720}},
			{"Ru"   , {19.26740,12.91820,4.863370,1.567560,5.378740,0.8085200,8.434670,24.79970,94.29280}},
			{"Ru3+" , {18.56380,13.28850,9.326020,3.009640,-3.189200,0.8473290,8.371640,1.7662000E-02,22.88700}},
			{"Ru+4" , {18.50030,13.17870,4.713040,2.185350,1.423570,0.8445820,8.125340,3.6495000E-02,20.85040}},
			{"Rh"   , {19.29570,14.35010,4.734250,1.289180,5.328000,0.7515360,8.217580,25.87490,98.60620}},
			{"Rh3+" , {18.87850,14.12590,3.325150,-6.198900,11.86780,0.7642520,7.844380,21.24870,-1.0360000E-02}},
			{"Rh4+" , {18.85450,13.98060,2.534640,-5.652600,11.28350,0.7608250,7.624360,19.33170,-1.0200000E-02}},
			{"Pd"   , {19.33190,15.50170,5.295370,0.6058440,5.265930,0.6986550,7.989290,25.20520,76.89860}},
			{"Pd2+" , {19.17010,15.20960,4.322340,0.0000000E+00,5.291600,0.6962190,7.555730,22.50570,0.0000000E+00}},
			{"Pd4+" , {19.24930,14.79000,2.892890,-7.949200,13.01740,0.6838390,7.148330,17.91440,5.1270002E-03}},
			{"Ag"   , {19.28080,16.68850,4.804500,1.046300,5.179000,0.6446000,7.472600,24.66050,99.81560}},
			{"Ag1+" , {19.18120,15.97190,5.274750,0.3575340,5.215720,0.6461790,7.191230,21.73260,66.11470}},
			{"Ag2+" , {19.16430,16.24560,4.370900,0.0000000E+00,5.214040,0.6456430,7.185440,21.40720,0.0000000E+00}},
			{"Cd"   , {19.22140,17.64440,4.461000,1.602900,5.069400,0.5946000,6.908900,24.70080,87.48250}},
			{"Cd2+" , {19.15140,17.25350,4.471280,0.0000000E+00,5.119370,0.5979220,6.806390,20.25210,0.0000000E+00}},
			{"In"   , {19.16240,18.55960,4.294800,2.039600,4.939100,0.5476000,6.377600,25.84990,92.80290}},
			{"In3+" , {19.10450,18.11080,3.788970,0.0000000E+00,4.996350,0.5515220,6.324700,17.35950,0.0000000E+00}},
			{"Sn"   , {19.18890,19.10050,4.458500,2.466300,4.782100,5.830300,0.5031000,26.89090,83.95710}},
			{"Sn2+" , {19.10940,19.05480,4.564800,0.4870000,4.786100,0.5036000,5.837800,23.37520,62.20610}},
			{"Sn4+" , {18.93330,19.71310,3.418200,1.9300001E-02,3.918200,5.764000,0.4655000,14.00490,-0.7583000}},
			{"Sb"   , {19.64180,19.04550,5.037100,2.682700,4.590900,5.303400,0.4607000,27.90740,75.28250}},
			{"Sb3+" , {18.97550,18.93300,5.107890,0.2887530,4.696260,0.4671960,5.221260,19.59020,55.51130}},
			{"Sb5+" , {19.86850,19.03020,2.412530,0.0000000E+00,4.692630,5.448530,0.4679730,14.12590,0.0000000E+00}},
			{"Te"   , {19.96440,19.01380,6.144870,2.523900,4.352000,4.817420,0.4208850,28.52840,70.84030}},
			{"I"    , {20.14720,18.99490,7.513800,2.273500,4.071200,4.347000,0.3814000,27.76600,66.87760}},
			{"I1-"  , {20.23320,18.99700,7.806900,2.886800,4.071400,4.357900,0.3815000,29.52590,84.93040}},
			{"Xe"   , {20.29330,19.02980,8.976700,1.990000,3.711800,3.928200,0.3440000,26.46590,64.26580}},
			{"Cs"   , {20.38920,19.10620,10.66200,1.495300,3.335200,3.569000,0.3107000,24.38790,213.9040}},
			{"Cs1+" , {20.35240,19.12780,10.28210,0.9615000,3.279100,3.552000,0.3086000,23.71280,59.45650}},
			{"Ba"   , {20.33610,19.29700,10.88800,2.695900,2.773100,3.216000,0.2756000,20.20730,167.2020}},
			{"Ba2+" , {20.18070,19.11360,10.90540,0.7736340,3.029020,3.213670,0.2833100,20.05580,51.74600}},
			{"La"   , {20.57800,19.59900,11.37270,3.287190,2.146780,2.948170,0.2444750,18.77260,133.1240}},
			{"La3+" , {20.24890,19.37630,11.63230,0.3360480,2.408600,2.920700,0.2506980,17.82110,54.94530}},
			{"Ce"   , {21.16710,19.76950,11.85130,3.330490,1.862640,2.812190,0.2268360,17.60830,127.1130}},
			{"Ce3+" , {20.80360,19.55900,11.93690,0.6123760,2.090130,2.776910,0.2315400,16.54080,43.16920}},
			{"Ce4+" , {20.32350,19.81860,12.12330,0.1445830,1.591800,2.659410,0.2188500,15.79920,62.23550}},
			{"Pr"   , {22.04400,19.66970,12.38560,2.824280,2.058300,2.773930,0.2220870,16.76690,143.6440}},
			{"Pr3+" , {21.37270,19.74910,12.13290,0.9751800,1.771320,2.645200,0.2142990,15.32300,36.40650}},
			{"Pr4+" , {20.94130,20.05390,12.46680,0.2966890,1.242850,2.544670,0.2024810,14.81370,45.46430}},
			{"Nd"   , {22.68450,19.68470,12.77400,2.851370,1.984860,2.662480,0.2106280,15.88500,137.9030}},
			{"Nd3+" , {21.96100,19.93390,12.12000,1.510310,1.475880,2.527220,0.1992370,14.17830,30.87170}},
			{"Pm"   , {23.34050,19.60950,13.12350,2.875160,2.028760,2.562700,0.2020880,15.10090,132.7210}},
			{"Pm3+" , {22.55270,20.11080,12.06710,2.074920,1.194990,2.417400,0.1857690,13.12750,27.44910}},
			{"Sm"   , {24.00420,19.42580,13.43960,2.896040,2.209630,2.472740,0.1964510,14.39960,128.0070}},
			{"Sm3+" , {23.15040,20.25990,11.92020,2.714880,0.9545860,2.316410,0.1740810,12.15710,24.82420}},
			{"Eu"   , {24.62740,19.08860,13.76030,2.922700,2.574500,2.387900,0.1942000,13.75460,123.1740}},
			{"Eu2+" , {24.00630,19.95040,11.80340,3.872430,1.363890,2.277830,0.1735300,11.60960,26.51560}},
			{"Eu3+" , {23.74970,20.37450,11.85090,3.265030,0.7593440,2.222580,0.1639400,11.31100,22.99660}},
			{"Gd"   , {25.07090,19.07980,13.85180,3.545450,2.419600,2.253410,0.1819510,12.93310,101.3980}},
			{"Gd3+" , {24.34660,20.42080,11.87080,3.714900,0.6450890,2.135530,0.1555250,10.57820,21.70290}},
			{"Tb"   , {25.89760,18.21850,14.31670,2.953540,3.583240,2.242560,0.1961430,12.66480,115.3620}},
			{"Tb3+" , {24.95590,20.32710,12.24710,3.773000,0.6919670,2.056010,0.1495250,10.04990,21.27730}},
			{"Dy"   , {26.50700,17.63830,14.55960,2.965770,4.297280,2.180200,0.2021720,12.18990,111.8740}},
			{"Dy3+" , {25.53950,20.28610,11.98120,4.500730,0.6896900,1.980400,0.1433840,9.349720,19.58100}},
			{"Ho"   , {26.90490,17.29400,14.55830,3.638370,4.567960,2.070510,0.1979400,11.44070,92.65660}},
			{"Ho3+" , {26.12960,20.09940,11.97880,4.936760,0.8527950,1.910720,0.1393580,8.800180,18.59080}},
			{"Er"   , {27.65630,16.42850,14.97790,2.982330,5.920460,2.073560,0.2235450,11.36040,105.7030}},
			{"Er3+" , {26.72200,19.77480,12.15060,5.173790,1.176130,1.846590,0.1372900,8.362250,17.89740}},
			{"Tm"   , {28.18190,15.88510,15.15420,2.987060,6.756210,2.028590,0.2388490,10.99750,102.9610}},
			{"Tm3+" , {27.30830,19.33200,12.33390,5.383480,1.639290,1.787110,0.1369740,7.967780,17.29220}},
			{"Yb"   , {28.66410,15.43450,15.30870,2.989630,7.566720,1.988900,0.2571190,10.66470,100.4170}},
			{"Yb2+" , {28.12090,17.68170,13.33350,5.146570,3.709830,1.785030,0.1599700,8.183040,20.39000}},
			{"Yb3+" , {27.89170,18.76140,12.60720,5.476470,2.260010,1.732720,0.1387900,7.644120,16.81530}},
			{"Lu"   , {28.94760,15.22080,15.10000,3.716010,7.976280,1.901820,9.985190,0.2610330,84.32980}},
			{"Lu3+" , {28.46280,18.12100,12.84290,5.594150,2.975730,1.682160,0.1422920,7.337270,16.35350}},
			{"Hf"   , {29.14400,15.17260,14.75860,4.300130,8.581540,1.832620,9.599900,0.2751160,72.02900}},
			{"Hf4+" , {28.81310,18.46010,12.72850,5.599270,2.396990,1.591360,0.1289030,6.762320,14.03660}},
			{"Ta"   , {29.20240,15.22930,14.51350,4.764920,9.243540,1.773330,9.370460,0.2959770,63.36440}},
			{"Ta5+" , {29.15870,18.84070,12.82680,5.386950,1.785550,1.507110,0.1167410,6.315240,12.42440}},
			{"W"    , {29.08180,15.43000,14.43270,5.119820,9.887500,1.720290,9.225900,0.3217030,57.05600}},
			{"W6+"  , {29.49360,19.37630,13.05440,5.064120,1.010740,1.427550,0.1046210,5.936670,11.19720}},
			{"Re"   , {28.76210,15.71890,14.55640,5.441740,10.47200,1.671910,9.092270,0.3505000,52.08610}},
			{"Os"   , {28.18940,16.15500,14.93050,5.675890,11.00050,1.629030,8.979480,0.3826610,48.16470}},
			{"Os4+" , {30.41900,15.26370,14.74580,5.067950,6.498040,1.371130,6.847060,0.1651910,18.00300}},
			{"Ir"   , {27.30490,16.72960,15.61150,5.833770,11.47220,1.592790,8.865530,0.4179160,45.00110}},
			{"Ir3+" , {30.41560,15.86200,13.61450,5.820080,8.279030,1.343230,7.109090,0.2046330,20.32540}},
			{"Ir4+" , {30.70580,15.55120,14.23260,5.536720,6.968240,1.309230,6.719830,0.1672520,17.49110}},
			{"Pt"   , {27.00590,17.76390,15.71310,5.783700,11.68830,1.512930,8.811740,0.4245930,38.61030}},
			{"Pt2+" , {29.84290,16.72240,13.21530,6.352340,9.853290,1.329270,7.389790,0.2632970,22.94260}},
			{"Pt4+" , {30.96120,15.98290,13.73480,5.920340,7.395340,1.248130,6.608340,0.1686400,16.93920}},
			{"Au"   , {16.88190,18.59130,25.55820,5.860000,12.06580,0.4611000,8.621600,1.482600,36.39560}},
			{"Au1+" , {28.01090,17.82040,14.33590,6.580770,11.22990,1.353210,7.739500,0.3567520,26.40430}},
			{"Au3+" , {30.68860,16.90290,12.78010,6.523540,9.096800,1.219900,6.828720,0.2128670,18.65900}},
			{"Hg"   , {20.68090,19.04170,21.65750,5.967600,12.60890,0.5450000,8.448400,1.572900,38.32460}},
			{"Hg1+" , {25.08530,18.49730,16.88830,6.482160,12.02050,1.395070,7.651050,0.4433780,28.22620}},
			{"Hg2+" , {29.56410,18.06000,12.83740,6.899120,10.62680,1.211520,7.056390,0.2847380,20.74820}},
			{"Tl"   , {27.54460,19.15840,15.53800,5.525930,13.17460,0.6551500,8.707510,1.963470,45.81490}},
			{"Tl1+" , {21.39850,20.47230,18.74780,6.828470,12.52580,1.471100,0.5173940,7.434630,28.84820}},
			{"Tl3+" , {30.86950,18.38410,11.93280,7.005740,9.802700,1.100800,6.538520,0.2190740,17.21140}},
			{"Pb"   , {31.06170,13.06370,18.44200,5.969600,13.41180,0.6902000,2.357600,8.618000,47.25790}},
			{"Pb2+" , {21.78860,19.56820,19.14060,7.011070,12.47340,1.336600,0.4883830,6.772700,23.81320}},
			{"Pb4+" , {32.12440,18.80030,12.01750,6.968860,8.084280,1.005660,6.109260,0.1470410,14.71400}},
			{"Bi"   , {33.36890,12.95100,16.58770,6.469200,13.57820,0.7040000,2.923800,8.793700,48.00930}},
			{"Bi3+" , {21.80530,19.50260,19.10530,7.102950,12.47110,1.235600,6.241490,0.4699990,20.31850}},
			{"Bi5+" , {33.53640,25.09460,19.24970,6.915550,-6.799400,0.9165400,3.9042000E-02,5.714140,12.82850}},
			{"Po"   , {34.67260,15.47330,13.11380,7.025880,13.67700,0.7009990,3.550780,9.556420,47.00450}},
			{"At"   , {35.31630,19.02110,9.498870,7.425180,13.71080,0.6858700,3.974580,11.38240,45.47150}},
			{"Rn"   , {35.56310,21.28160,8.003700,7.443300,13.69050,0.6631000,4.069100,14.04220,44.24730}},
			{"Fr"   , {35.92990,23.05470,12.14390,2.112530,13.72470,0.6464530,4.176190,23.10520,150.6450}},
			{"Ra"   , {35.76300,22.90640,12.47390,3.210970,13.62110,0.6163410,3.871350,19.98870,142.3250}},
			{"Ra2+" , {35.21500,21.67000,7.913420,7.650780,13.54310,0.6049090,3.576700,12.60100,29.84360}},
			{"Ac"   , {35.65970,23.10320,12.59770,4.086550,13.52660,0.5890920,3.651550,18.59900,117.0200}},
			{"Ac3+" , {35.17360,22.11120,8.192160,7.055450,13.46370,0.5796890,3.414370,12.91870,25.94430}},
			{"Th"   , {35.56450,23.42190,12.74730,4.807030,13.43140,0.5633590,3.462040,17.83090,99.17220}},
			{"Th4+" , {35.10070,22.44180,9.785540,5.294440,13.37600,0.5550540,3.244980,13.46610,23.95330}},
			{"Pa"   , {35.88470,23.29480,14.18910,4.172870,13.42870,0.5477510,3.415190,16.92350,105.2510}},
			{"U"    , {36.02280,23.41280,14.94910,4.188000,13.39660,0.5293000,3.325300,16.09270,100.6130}},
			{"U3+"  , {35.57470,22.52590,12.21650,5.370730,13.30920,0.5204800,3.122930,12.71480,26.33940}},
			{"U4+"  , {35.37150,22.53260,12.02910,4.798400,13.26710,0.5165980,3.050530,12.57230,23.45820}},
			{"U6+"  , {34.85090,22.75840,14.00990,1.214570,13.16650,0.5070790,2.890300,13.17670,25.20170}},
			{"Np"   , {36.18740,23.59640,15.64020,4.185500,13.35730,0.5119290,3.253960,15.36220,97.49080}},
			{"Np3+" , {35.01360,22.72860,14.38840,1.756690,13.11300,0.4898100,2.810990,12.33000,22.65810}},
			{"Np4+" , {36.52540,23.80830,16.77070,3.479470,13.38120,0.4993840,3.263710,14.94550,105.9800}},
			{"Np6+" , {35.70740,22.61300,12.98980,5.432270,13.25440,0.5023220,3.038070,12.14490,25.49280}},
			{"Pu"   , {35.51030,22.57870,12.77660,4.921590,13.21160,0.4986260,2.966270,11.94840,22.75020}},
			{"Pu3+" , {35.84000,22.71690,13.58070,5.660160,13.19910,0.4849380,2.961180,11.53310,24.39920}},
			{"Pu4+" , {35.64930,22.64600,13.35950,5.188310,13.15550,0.4814220,2.890200,11.31600,21.83010}},
			{"Pu6+" , {35.17360,22.71810,14.76350,2.286780,13.05820,0.4732040,2.738480,11.55300,20.93030}},
			{"Am"   , {36.67060,24.09920,17.34150,3.493310,13.35920,0.4836290,3.206470,14.31360,102.2730}},
			{"Cm"   , {36.64880,24.40960,17.39900,4.216650,13.28870,0.4651540,3.089970,13.43460,88.48340}},
			{"Bk"   , {36.78810,24.77360,17.89190,4.232840,13.27540,0.4510180,3.046190,12.89460,86.00300}},
			{"Cf"   , {36.91850,25.19950,18.33170,4.243910,13.26740,0.4375330,3.007750,12.40440,83.78810}},
			{NULL   , {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}
		} ;
	

		/* search for the right sort of atom, index will be i */
		i = 0 ;
		while ( ( asf[i].atom != NULL ) && ( strcmp ( asf[i].atom, entry ) !=0 ) )
			++i ;
		if ( asf[i].atom == NULL )
		{
			XNDIFF_ERROR(7) ;
		}
		
		/* http://www.ruppweb.org/xray/comp/scatfac.htm */	
		/* f(s) = f0(s) + f'(s) + i*f''(s) */
		/* here only the nondispersive part f0(s) is computed, anamalous scattering contributions are neglected */
		/* f0[k] = c + [\sum_ {i=0}^3 a_i*exp(-b_i*(k^2)) ] */
		/* k = sin(theta) / lambda = s / 2 , factor 10 due to [1/nm] to [1/A] conversion */		
		s = SQUARE ( s / 20.0 ) ;
		f = asf[i].par[4] ;
		for (int j=0; j<4; ++j )
			f += asf[i].par[j] * exp ( -asf[i].par[j+5]*s ) ;

		free(entry);

		return ( f ) ;
	}


	/* ggt returns greatest common divisor of a[0], a[1], ..., a[num-1] */
	/* num must be greater than 0 */
	/* for zero vector a of length num, the greatest common divisor is defined as 1 */
	int ggt ( int *a, int num )
	{
		int i, j ;
		int min ;
		int ggt_flag ;
	
		min = a[0] ;
		for ( i=1; i<num; ++i )
		{
			if ( a[i] < min )
				min = a[i] ;
		}
		for ( i=min; i>1; --i )
		{
			ggt_flag = 1 ;
			for ( j=0; j<num; ++j )
			{
				if ( i % a[j] != 0 )
				{
					ggt_flag = 0 ;
					break ;
				}
			}
			if ( ggt_flag == 1 )
				return ( i ) ;
		}
		return ( 1 ) ;
	}
	
	/* returns greatest common divisor of three integers, using the ggt function above */
	/* GGT(0,0,0)=1 */
	int GGT ( int a, int b, int c )
	{
		int d[3] ;
		d[0] = a ;
		d[1] = b ;
		d[2] = c ;
		return ( ggt ( d, 3 ) ) ;
	}

	
	/* function delallspc deletes all spaces in the input string str */
	/* http://www.cs.bu.edu/teaching/cpp/string/array-vs-ptr/ */
	void delallspc (char *str)
	{
		char dummy[strlen(str)+1] ;

		for (; *str; ++str)
		{
			if (*str == ' ') 
			{
				strcpy( dummy, str+1) ;
				strcpy( str, dummy) ;
				--str ;
			}
		}
	}

	/* scalar vector product in R3*/
	double inline sp (double *a, double *b) { return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]) ; } 
	/* vector product in R3*/
	double inline *vp (double *erg, double *a, double *b)     
	{ 
		erg[0] = (a[1]*b[2]-a[2]*b[1]) ; 
		erg[1] = (a[2]*b[0]-a[0]*b[2]) ; 
		erg[2] = (a[0]*b[1]-a[1]*b[0]) ; 
		return erg ;
	}
	/*determinant in R3*/
	double inline det (double *a, double *b, double *c)   
	{ return (a[0]*b[1]*c[2]+a[1]*b[2]*c[0]+a[2]*b[0]*c[1]-a[2]*b[1]*c[0]-a[0]*b[2]*c[1]-a[1]*b[0]*c[2]) ; }
	/* absolut value */
	double inline betrag (double *a) { return (sqrt(SQUARE(a[0])+SQUARE(a[1])+SQUARE(a[2]))) ; }
	/* scalar multiplication with double scalar and double vector */
	double inline *csp (double *erg, double c, double *b)
	{
		erg[0] = c*b[0] ;
		erg[1] = c*b[1] ;
		erg[2] = c*b[2] ;
		return erg ;
	}
	/* vector addition */
	double inline *vadd (double *a, double *b)
	{
		a[0] = a[0]+b[0] ;
		a[1] = a[1]+b[1] ;
		a[2] = a[2]+b[2] ;
		return a ;
	}
	/* vector subtraction */
	double inline *vsub (double *a, double *b)
	{
		a[0] = a[0]-b[0] ;
		a[1] = a[1]-b[1] ;
		a[2] = a[2]-b[2] ;
		return a ;
	}  
	/* square of double numbers */
	double inline SQUARE (double a)
	{
		return a*a;
	}
	/* absolute squared value of a complex number */
	/* http://www.cplusplus.com/reference/std/complex/ */
	double inline cabs2(dcmplx a)
	{
		return pow( abs(a), 2) ;
	}
};



int main (int carg, char **varg)
{
	XNDiff job=XNDiff();
	job.eval_cmd(carg,varg);
	
	return (0) ;
}
