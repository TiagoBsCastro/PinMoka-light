#include "util.h"

using namespace CCfits;

// Pinocchio PLC DATA Structure reading opp
std::istream & operator>>(std::istream &input, DATA &Data)
{
  input.read((char *)&Data, sizeof(Data));
  return input;
};

inline bool exists_file (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

// read PLC Pinocchio file
void readPinocchio(struct PinocchioPLC *pinPLC, struct InputParams p){
  int di0=0;
  std:: ifstream ifilinp(p.PinocchioFile.c_str());
  if(ifilinp.is_open()){
    int dummy;
    DATA data;
    int n=0;
    double ra0P=0;
    double dec0P=0;
    while( ifilinp.read( (char *)&dummy, sizeof(dummy) ) ){
       ifilinp >> data;
       ifilinp.read((char *)&dummy, sizeof(dummy));
       /*std:: cout << data.group_id << "   " <<  data.group_mass << "   " <<  data.true_z << "   " <<  data.pos[0] << "   "
       <<  data.pos[1] << "   " <<  data.pos[2] << "   " <<  data.vel[0] << "   " <<  data.vel[1] << "   "
       <<  data.vel[2] << "   " <<  data.theta << "   " <<  data.phi << "   " <<  data.pv << "   "  << data.obs_z << std:: endl;
       */
      if(data.group_id!=idhalo || fabs(data.true_z-zhalo)>1e-4){
        pinPLC->id.push_back(data.group_id);
        pinPLC->mass.push_back( data.group_mass );
        pinPLC->redshift.push_back( data.true_z );
        pinPLC->x_c.push_back( data.pos[0] );
        pinPLC->y_c.push_back( data.pos[1] );
        pinPLC->z_c.push_back( data.pos[2] );
        pinPLC->vx.push_back( data.vel[0] );
        pinPLC->vy.push_back( data.vel[1] );
        pinPLC->vz.push_back( data.vel[2] );
        pinPLC->theta.push_back( data.theta );
        pinPLC->phi.push_back( data.phi );
        pinPLC->pec_vel.push_back( data.pv );
        pinPLC->obs_redshift.push_back( data.obs_z );
        double ri = sqrt(data.pos[0]*data.pos[0]+data.pos[1]*data.pos[1]+data.pos[2]*data.pos[2]);
        double xc = ri*sin((90-data.theta)*M_PI/180.)*cos(data.phi*M_PI/180.);
        double yc = ri*sin((90-data.theta)*M_PI/180.)*sin(data.phi*M_PI/180.);
        double zc = ri*cos((90-data.theta)*M_PI/180.);
        double dec = asin(xc/ri);
        double ra  = atan(yc/zc);
        pinPLC->ra.push_back( ra );
        pinPLC->dec.push_back( dec );
        n++;
      }else{
        if(idhalo>0){
          // has been read from the file compute its quantities
          double ri = sqrt(data.pos[0]*data.pos[0]+data.pos[1]*data.pos[1]+data.pos[2]*data.pos[2]);
          double xc = ri*sin((90-data.theta)*M_PI/180.)*cos(data.phi*M_PI/180.);
          double yc = ri*sin((90-data.theta)*M_PI/180.)*sin(data.phi*M_PI/180.);
          double zc = ri*cos((90-data.theta)*M_PI/180.);
          dec0P = asin(xc/ri);
          ra0P  = atan(yc/zc);
          // check if we can build such field of view without wrapping the boundaries of the cone...
          std:: cout << " center in " << ra0P << "  " << dec0P << std:: endl;
        }
      }
    }
    pinPLC->nhaloes=n;
    if(di0==0){
      // need to recenter ra and dec if idhalo>0
      for(int i=0;i<pinPLC->nhaloes;i++){
        pinPLC->ra[i] -= ra0P;
        pinPLC->dec[i] -= dec0P;
      }
    }else{
      std:: cout << " wrapping cone boundary conditions not implemented yet! " << std:: endl;
      exit(1);
      // need to recenter ra and dec if idhalo>0 and wrap period boundaries conditions of the cone
    }
    std:: cout << " nhaloes in Pinocchio PLC file = " << n << std:: endl;
    std:: cout << "  " << std:: endl;
    //
    //std:: vector<long> id;
    //std:: vector<double> mass,redshift,x_c,y_c,z_c,vx,vy,vz,theta,phi,pec_vel,obs_redshift;
  }else{
    std:: cout << " the PLC file " << p.PinocchioFile << " does not exist " << std:: endl;
    std:: cout << " I will STOP here!!! " << std:: endl;
    exit(1);
  }
}

// read the input file "weaklensingMOKA.ini"
void readInput(struct InputParams *p, std::string name){
  // read fileinput
  std:: ifstream fin (name.c_str());
  if(fin.is_open());
  else{
    std:: cout << " weaklensingMOKA.ini file does not exist where you are running the code " << std:: endl;
    std:: cout << " I will STOP here!!! " << std:: endl;
    exit(1);
  }
  std:: string str;
  fin >> str;
  fin >> p->omegam;          //  1. omega matter
  fin >> str;
  fin >> p->omegal;          //  2. omega lambda
  fin >> str;
  fin >> p->h0;              //  3. hubble constant in unit of 100
  fin >> str;
  fin >> p->wq;              //  4. dark energy equation of state parameter
  fin >> str;
  fin >> p->nzs;             //  5. number of source redshifts to consider
  if(p->nzs<-1){
    // build differencial convergence maps
    lenses = 1;
    std:: cout << " number of zs " << p->nzs << std:: endl;
    std:: cout << " I will build differencial convergence maps " << std:: endl;
    p->nzs = -p->nzs;
  }else{
    if(fabs(p->nzs)<1e-4){
      std:: cout << " number of zs " << p->nzs << std:: endl;
      std:: cout << " check this out, I will STOP here!!! " << std:: endl;
      exit(1);
    }
    lenses = 0;
  }
  fin >> str;
  p->z1.resize(p->nzs);
  p->zs.resize(p->nzs);
  for(int i=0;i<p->nzs;i++){
    fin >> p->zs[i];              //  5. source redshift (1:nzs)
  }
  fin >> str;
  fin >> p->cutR;            //  6. radius at which cut the density profile in 1/Rvir
  fin >> str;
  fin >> p->nx;              // 7. number of pixels in x
  fin >> str;
  fin >> p->ny;              // 8. number of pixels in y
  fin >> str;
  fin >> p->fx;              // 9. size of the field of view in x
  fin >> str;
  fin >> p->fy;              // 10. size of the field of view in y
  fin >> str;
  fin >> p->filsigma;        // 11. file with lm s relation to be adopted
  fin >> str;
  fin >> p->cMrelation;      // 12. c-M relation 
  fin >> str;
  fin >> p->sigmalnC;        // 13. scatter in the c-M relation if a c-M model is used
  fin >> str;
  fin >> p->halodef;         // 14. vir use FOF haloes
  fin >> str;
  fin >> p->simcase;         // 15. MapSim files or Pinocchio PLC file
  fin >> str;
  fin >> p->PinocchioFile;   // 16. Pinocchio PLC file, needed for Pincchio and BegognaCAT
  if(fabs(p->cutR)<1e-4){
    std:: cout << " cutR parameter too small = " << p->cutR << std:: endl;
    std:: cout << " check this out ... for now I will STOP here!!! " << std:: endl;
    exit(1);
  }

  std:: cout << "  " << std:: endl;
  std:: cout << "  ___________________________________________"  << std:: endl;
  std:: cout << "  running with the following paramters:  " << std:: endl;
  std:: cout << "  Omegam           = " << p->omegam << std:: endl;
  std:: cout << "  Omegal           = " << p->omegal << std:: endl;
  std:: cout << "  Hubble parameter = " << p->h0 << std:: endl;
  std:: cout << "  w                = " << p->wq << std:: endl;
  p->z1[0] = -1;
  if(p->nzs==1){
    std:: cout << "  source redshift  = " << p->zs[0] << std:: endl;
  }else{
    std:: cout << "  source redshifts  = ";
    for(int i=0;i<p->nzs;i++){
      std:: cout << p->zs[i] << "  ";
      if(i>0){
	if(lenses==1){
	  p->z1[i] = p->zs[i-1];
	}else{
	  p->z1[i] = 0;
	}
      }
    }
    std:: cout << "  " << std:: endl;
  }
  std:: cout << "  Rcut             = " << p->cutR << std:: endl;
  std:: cout << "  path catalogues  = " << p->pathcats << std:: endl;
  std:: cout << "  nx = " << p->nx << "  ny = " << p->ny << std:: endl;
  std:: cout << "  field of view: x = " << p->fx << "  y = " << p->fy << std:: endl;
  std:: ifstream fillms;
  fillms.open(p->filsigma.c_str());
  if(fillms.is_open()){
    std:: cout << "  file lm s rel    = " << p->filsigma << std:: endl;
  }else{
    if(p->filsigma=="sim"){
      std:: cout << " no file for the lm s relation ... I will use Vmax and Rmax from the simulation " << std:: endl;
      std:: cout << " since filsigma variable is set equal to sim " << std:: endl;
    }else{
      std:: cout << " no file for the lm s relation ... I will use Neto08+Bullock01 " << std:: endl;
    }
  }
  std:: cout << "  scatter lnC      = " << p->sigmalnC << std:: endl;
  if(p->simcase=="Pinocchio"){
    std:: cout << " I will build the map assuming "<< p->halodef <<" masses " << std:: endl;
    std:: cout <<  " Effective convergence map will be built using Pinocchio PLC output file " << std:: endl;
    std:: cout << p->PinocchioFile << std:: endl;
  }else{
    std:: cout << " wrong parameter in the input file for simcase " << p->simcase  << std:: endl;
    std:: cout << " I will STOP here!!! " << std:: endl;
    exit(1);
  }
  std:: cout << "  ___________________________________________"  << std:: endl;
  std:: cout << "  " << std:: endl;
}

double gslf_fc(double lc, void *ip){
  struct gslfparams_fc *fp  = (struct gslfparams_fc *) ip;
  double c = pow(10.,lc);
  return log(200./3.*gsl_pow_3(c)/(log(1.+c)-c/(1.+c))) - log(14.426*gsl_pow_2(fp->vmax/(fp->h*100.)/(fp->rmax/1000.)));
}

double getC(double vmax,double rmax,double h){
  int status;
  int iter = 0, max_iter = 512;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
  struct  gslfparams_fc params = {vmax,rmax,h};
  F.function = &gslf_fc;
  F.params = &params;
  double lc;
  double esp = 1e-4;
  double x_lo = -3, x_hi = 4.;
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      lc = gsl_root_fsolver_root (s);
      x_lo = gsl_root_fsolver_x_lower (s);
      x_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (x_lo, x_hi,
				       0, esp);
      //if (status == GSL_SUCCESS)
      //	printf ("Converged:\n");
      //std:: cout << lc << std:: endl;
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fsolver_free (s);
  return pow(10.,lc);
}

double getCNeto(double m, double z){
  // Bullock+01 redshift evolution
  double cm = 4.67*pow(m/1.e+14,-0.11)/(1+z);
  return cm;
}

// starting will show when the code has been run
void showStart(){
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  std:: cout << "  " << std:: endl;
  std:: cout << " ********************************************************** " << std:: endl;
  std:: cout << " *                                                        * " << std:: endl;
  std:: cout << " *       Weak lensing Light-Cones from MOKA-Haloes        * " << std:: endl;
  std:: cout << " *                                                        * " << std:: endl;
  std:: cout << " *                                                        * " << std:: endl;
  std:: cout << "  " << std:: endl;
  std:: cout << "  running on: " << asctime(timeinfo) << std:: endl;
  std:: cout << " ********************************************************** " << std:: endl;
}

void writeFits(std:: string filename,std:: valarray<float> f, int npix, int npixy,double zs,
	       double fxrad,double fyrad, double ra_c, double dec_c){
  // always even
  double xo = double(npix)/2.+1;
  double yo = double(npixy)/2.+1;
  long naxis=2;
  long naxes[2]={npix,npixy};
  std::auto_ptr<FITS> fout(new FITS(filename,FLOAT_IMG,naxis,naxes));
  std::vector<long> naxex(2);
  naxex[0]=npix;
  naxex[1]=npixy;
  PHDU *phout=&fout->pHDU();
  phout->write( 1, npix*npixy, f );
  phout->addKey ("ZSOURCE",zs, "source redshift");
  phout->addKey ("RADESYS" , "FK5" , "Coordinate system");
  phout->addKey ("CTYPE1" , "RA---TAN" , "Coordinates -- projection");
  phout->addKey ("CTYPE2" , "DEC--TAN" , "Coordinates -- projection");
  phout->addKey ("EQUINOX" , 2000 , "Epoch of the equinox");
  phout->addKey ("CRPIX1",xo,"X reference pixel");
  phout->addKey ("CRPIX2",yo,"Y reference pixel");
  phout->addKey ("CUNIT1" , "deg" , " ");
  phout->addKey ("CUNIT2" , "deg" , " ");
  double xov=ra_c;
  double yov=dec_c;
  phout->addKey ("CRVAL1",xov,"Reference longitude");
  phout->addKey ("CRVAL2",yov,"Reference latitude");
  double dltx = -fxrad/npix*180./M_PI;
  double dlty = fyrad/npixy*180./M_PI;
  phout->addKey ("CDELT1",dltx,"X scale");
  phout->addKey ("CDELT2",dlty,"Y scale");
}

double uNFW(double k, nfwHalo *ha, double r){
  double m200 = ha->getMass();
  double rs = ha->getScaleRadius();
  double m = ha->massin(r);
  double f1 = 4*M_PI*ha->getDensity()*gsl_pow_3(rs)/m;
  double cR = r/rs; // equal to c if r = Rvir
  double Si=gsl_sf_Si(k*rs);
  double Ci=gsl_sf_Ci(k*rs);
  double Sic=gsl_sf_Si((1+cR)*k*rs);
  double Cic=gsl_sf_Ci((1+cR)*k*rs);
  double uk=f1*(sin(k*rs)*(Sic-Si) - sin(cR*k*rs)/(1+cR)/k/rs + cos(k*rs)*(Cic-Ci));
  return uk;
}

void locateHalo(struct InputParams p, double ra, double dec, double dr, int &imin, int &imax, int &jmin, int &jmax){
  // ----> ( CASE 1 )
  if(ra>=p.fxmin && ra<=p.fxmax && dec>=p.fymin && dec<=p.fymax){
    p.xmin = p.xmin>=p.fxmin ? p.xmin : p.fxmin;
    p.xmax = p.xmax<=p.fxmax ? p.xmax : p.fxmax;
    p.ymin = p.ymin>=p.fymin ? p.ymin : p.fymin;
    p.ymax = p.ymax<=p.fymax ? p.ymax : p.fymax;
    if(p.xmin == p.fxmin) imin = 0;
    else{
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }
    if(p.ymin == p.fymin) jmin = 0;
    else{
      jmin = locate(p.y,p.ymin);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }
    if(p.xmax == p.fxmax) imax = p.nx-1;
    else{
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }
    if(p.ymax == p.fymax) jmax = p.ny-1;
    else{
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }
  }
  // ----> ( CASE 2 )
  if(ra < p.fxmin && dec >= p.fymin && dec <= p.fymax){
    p.xmin = p.fxmin;
    imin = 0;
    p.xmax = ra + dr;
    double dx = p.fxmin - ra;
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dx));
    p.ymin = (dec - dy)>=p.fymin ? (dec - dy) : p.fymin;
    p.ymax = (dec + dr)<=p.fymax ? (dec + dy) : p.fymax;
    if(p.ymin == p.fymin) jmin = 0;
    else{
      jmin = locate(p.y,p.ymin);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }
    if(p.ymax == p.fymax) jmax = p.ny-1;
    else{
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }
    if(p.xmax>=p.fxmin){
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }else{
      imax=0;
    }
  }
  // ----> ( CASE 3 )
  if(ra > p.fxmax && dec >= p.fymin && dec <= p.fymax){
    p.xmin = ra - dr;
    p.xmax = p.fxmax;
    imax = p.nx - 1;
    double dx = ra - p.fxmax;
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dx));
    p.ymin = (dec - dy)>=p.fymin ? (dec - dy) : p.fymin ;
    p.ymax = (dec + dy)<=p.fymax ? (dec + dy) : p.fymax ;
    if(p.ymin == p.fymin) jmin = 0;
    else{
      jmin = locate(p.y,p.ymin);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }
    if(p.ymax == p.fymax) jmax = p.ny-1;
    else{
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }
    if(p.xmin<=p.fxmax){
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }else{
      imin=p.nx-1;
    }
  }
  // ----> ( CASE 4 )
  if(dec < p.fymin && ra >= p.fxmin && ra <=p.fxmax){
    p.ymin = p.fymin;
    jmin = 0;
    p.ymax = dec + dr;
    double dy = p.ymin - dec;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dy));
    p.xmin = (ra - dx)>=p.fxmin ? (ra - dx) : p.fxmin ;
    p.xmax = (ra + dx)<=p.fxmax ? (ra + dx) : p.fxmax ;
    if(p.xmin == p.fxmin) imin = 0;
    else{
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }
    if(p.xmax == p.fxmax) imax = p.nx-1;
    else{
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }
    if(p.ymax>=p.fymin){
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }else{
      jmax = 0;
    }
  }
  // ----> ( CASE 5 )
  if(dec > p.fymax && ra >= p.fxmin && ra <=p.fxmax){
    p.ymin = dec - dr;
    p.ymax = p.fymax;
    jmax = p.ny-1;
    double dy = dec-p.fymax;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dy));
    p.xmin = (ra - dx)>=p.fxmin ? (ra - dx) : p.fxmin ;
    p.xmax = (ra + dx)<=p.fxmax ? (ra + dx) : p.fxmax ;
    if(p.xmin == p.fxmin) imin = 0;
    else{
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }
    if(p.xmax == p.fxmax) imax = p.nx-1;
    else{
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }
    if(p.ymin<=p.fymax){
      jmin = locate(p.y,p.ymin);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }else{
      jmin = p.ny-1;
    }
  }
  // ---> ( CASE 6 )
  if(ra > p.fxmax && dec > p.fymax){
    p.xmax = p.fxmax;
    p.ymax = p.fymax;
    imax = p.nx-1;
    jmax = p.ny-1;
    double dxcm = ra -p.fxmax;
    double dycm = dec-p.fymax;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
    p.xmin = (ra - dx)>=p.fxmin ? (ra - dx) : p.fxmin ;
    p.ymin = (dec - dy)>=p.fymin ? (dec - dy) : p.fymin ;
    if(p.xmin == p.fxmin) imin = 0;
    else{
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }
    if(p.ymin == p.fymin) jmin = 0;
    else{
      jmin = locate(p.y,p.ymin);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }
  }
  // --> ( CASE 7 )
  if(ra > p.fxmax && dec < p.fymin){
    p.xmax = p.fxmax;
    imax = p.nx-1;
    p.ymin = p.fymin;
    jmin = 0;
    double dycm = p.fymin-dec;
    double dxcm = ra-p.fxmax;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
    p.xmin = (ra - dx)>=p.fxmin ? (ra - dx) : p.fxmin ;
    p.ymax = (dec + dy)<=p.fymax ? (dec + dy) : p.fymax ;
    if(p.xmin == p.fxmin) imin = 0;
    else{
      imin = locate(p.x,p.xmin);
      imin=GSL_MIN( GSL_MAX( imin, 0 ), p.nx-1 );
    }
    if(p.ymax == p.fxmax) jmax = p.ny-1;
    else{
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }
  }
  // ---> ( CASE 8 )
  if(ra < p.fxmin && dec < p.fymin){
    p.xmin = p.fxmin;
    p.ymin = p.fymin;
    imin = 0;
    jmin = 0;
    double dxcm = p.fxmin-ra;
    double dycm = p.fymin-dec;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
    p.xmax = (ra + dx)<=p.fxmax ? (ra + dx) : p.fxmax ;
    p.ymax = (dec + dy)<=p.fxmax ? (dec + dy) : p.fymax ;
    if(p.xmax == p.fxmax) imax = p.nx-1;
    else{
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }
    if(p.ymax == p.fymax) jmax = p.ny-1;
    else{
      jmax = locate(p.y,p.ymax);
      jmax=GSL_MIN( GSL_MAX( jmax, 0 ), p.ny-1 );
    }
  }
  // ---> ( CASE 9 )
  if(ra < p.fxmin && dec > p.fymax){
    p.xmin = p.fxmin;
    p.ymax = p.fymax;
    imin = 0;
    jmax = p.ny-1;
    double dxcm = p.fxmin-ra;
    double dycm = dec-p.fymax;
    double dx = sqrt(gsl_pow_2(dr) - gsl_pow_2(dycm));
    double dy = sqrt(gsl_pow_2(dr) - gsl_pow_2(dxcm));
    p.xmax = (ra + dx)<=p.fxmax ? (ra + dx) : p.fxmax ;
    p.ymin = (dec - dy)>=p.fymin ? (dec - dy) : p.fymin ;
    if(p.xmax == p.fxmax) imax = p.nx-1;
    else{
      imax = locate(p.x,p.xmax);
      imax=GSL_MIN( GSL_MAX( imax, 0 ), p.nx-1 );
    }
    if(p.ymax == p.fymin) jmin = 0;
    else{
      jmin = locate(p.y,p.ymax);
      jmin=GSL_MIN( GSL_MAX( jmin, 0 ), p.ny-1 );
    }
  }
}
