#include "../Moka/utilities.h"
#include "../Moka/cosmology.h"
#include "../Moka/halo.h"
#include "../Moka/nfwHalo.h"
#include "../Moka/nfwLens.h"
#include "../Moka/hernq_Halo.h"
#include "../Moka/hernq_Lens.h"
#include "../Moka/jaffe_Halo.h"
#include "../Moka/jaffe_Lens.h"
#include "../Moka/sisHalo.h"
#include "../Moka/sisLens.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <CCfits/CCfits>
#include "time.h"
#include <math.h>
#define READ_PLC_BINARY

extern long idhalo;
extern double zhalo;

extern const gsl_rng_type * Th;
extern gsl_rng * rh; // host halo concentration
extern const gsl_rng_type * Thfof;
extern gsl_rng * rhfof; // host halo concentration fof
extern int lenses;

// Pinocchio PLC DATA Structure
struct DATA
{

  unsigned long long int group_id;
  double true_z;
  double pos[3];
  double vel[3];
  double group_mass;
  double theta;
  double phi;
  double pv;
  double obs_z;

};

// fof properties
struct PinocchioPLC{
  int nhaloes;                 // pinocchio PLC file structure
  std:: vector<long> id;
  std:: vector<double> mass,redshift,x_c,y_c,z_c,vx,vy,vz,theta,phi,pec_vel,obs_redshift;
  std:: vector<double> concentration,ra,dec;
};

// plane list
struct PlaneList{
  std:: vector<int> plane;
  std:: vector<double> zl,zsnap;
  std:: vector<double> DlDOWN,DlUP;
  std:: vector<int> rep,snap;
};

struct gslfparams_fc{
  double vmax,rmax,h;
};

// parameters in the input file
struct InputParams{
  int nzs;               // number of source redshifts and number of cutRadius
  std:: vector<double> zs;     // redshifts of the sources
  std:: vector<double> z1;     // lower redshift in case lenses=1
  std:: string pathcats;       // path fof and subs cat
  double omegam,omegal,h0,wq;  // cosmological parameters
  int nx,ny;                   // number of pixels in x and y directions
  double fx,fy;                // size of the field of view in x and y
  std:: string cMrelation;     // c-M relation to be used
  std:: string filsigma;       // file contaning lm s relation used for c-M Zhao
  double sigmalnC;             // scatter in the c-M relation if a c-M model is used
  std:: string halodef;        // vir use FOF haloes
  std:: string PinocchioFile;  // PLC Pinocchio file or Begogna CAT
  // other varialbles in the structure
  double xmin,xmax,ymin,ymax;
  std:: vector<double> x,y;
  double fxmin,fxmax,fymin,fymax;
};

inline bool exists_file (const std::string& name);

void readPinocchio(struct PinocchioPLC *pinPLC, struct InputParams p);

void readInput(struct InputParams *p, std::string name);

double gslf_fc(double lc, void *ip);

double getC(double vmax,double rmax,double h);

double getCNeto(double m, double z);

void showStart();

void writeFits(std:: string filename,std:: valarray<float> f, int npix, int npixy,double zs,
	       double fxrad,double fyrad, double ra_c, double dec_c);

void locateHalo(struct InputParams p, double ra, double dec, double dr, int &imin, int &imax, int &jmin, int &jmax);
