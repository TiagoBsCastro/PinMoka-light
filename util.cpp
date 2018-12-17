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

#ifdef READ_PLC_BINARY
// read binary PLC Pinocchio file
  void readPinocchio(struct PinocchioPLC *pinPLC, struct InputParams p){
    std:: ifstream ifilinp(p.PinocchioFile.c_str());
    if(ifilinp.is_open()){
      int dummy;
      DATA data;
      int n=0;
      while( ifilinp.read( (char *)&dummy, sizeof(dummy) ) ){
        ifilinp >> data;
        ifilinp.read((char *)&dummy, sizeof(dummy));
        std:: cout << data.group_id << "   "  <<  data.true_z << "   " <<  data.pos[0] << "   "
        <<  data.pos[1] << "   " <<  data.pos[2] << "   " <<  data.vel[0] << "   " <<  data.vel[1] << "   "
        <<  data.vel[2] << "   " <<  data.group_mass << "   " <<  data.theta << "   " <<  data.phi << "   " <<  data.pv << "   "  << data.obs_z << std:: endl;

        if(fabs(data.true_z)>1e-4){
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
          double ra  = atan2(yc,zc);
          pinPLC->ra.push_back( ra );
          pinPLC->dec.push_back( dec );
          n++;
        }
      }
      pinPLC->nhaloes=n;
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
#else
  // read PLC ASCII Pinocchio file
  void readPinocchio(struct PinocchioPLC *pinPLC, struct InputParams p){
    std:: ifstream ifilinp;
    ifilinp.open(p.PinocchioFile.c_str());
    if(ifilinp.is_open()){
      for(int i=0;i<58;i++){
        string sbut;
        ifilinp >> sbut;
      }
      long id;
      double mass,redshift,x_c,y_c,z_c,vx,vy,vz,theta,phi,pec_vel,obs_redshift;
      int n=0;
      while(ifilinp >> id
        >> redshift
        >> x_c >> y_c >> z_c
        >> vx >> vy >> vz
        >> mass
        >> theta
        >> phi
        >> pec_vel
        >> obs_redshift)
        {
          std:: cout << id << "   " <<  mass << "   " <<  redshift << "   " <<  x_c << "   "
          <<  y_c << "   " <<  z_c << "   " <<  vx << "   " <<  vy << "   "
          <<  vz << "   " <<  theta << "   " <<  phi << "   " <<  pec_vel << "   "  << obs_redshift << std:: endl;
          if(fabs(redshift)>1e-4){
            pinPLC->id.push_back(id);
            pinPLC->mass.push_back(mass);
            pinPLC->redshift.push_back(redshift);
            pinPLC->x_c.push_back(x_c);
            pinPLC->y_c.push_back(y_c);
            pinPLC->z_c.push_back(z_c);
            pinPLC->vx.push_back(vx);
            pinPLC->vy.push_back(vy);
            pinPLC->vz.push_back(vz);
            pinPLC->theta.push_back(theta);
            pinPLC->phi.push_back(phi);
            pinPLC->pec_vel.push_back(pec_vel);
            pinPLC->obs_redshift.push_back(obs_redshift);
            double ri = sqrt(x_c*x_c + y_c*y_c + z_c*z_c);
            double xc = ri*sin((90-theta)*M_PI/180.)*cos(phi*M_PI/180.);
            double yc = ri*sin((90-theta)*M_PI/180.)*sin(phi*M_PI/180.);
            double zc = ri*cos((90-theta)*M_PI/180.);
            double dec = asin(xc/ri);
            double ra  = atan(yc/zc);
            pinPLC->ra.push_back(ra);
            pinPLC->dec.push_back(dec);
            n++;
          }
        }
        pinPLC->nhaloes=n;
        // need to recenter ra and dec if idhalo>0
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
#endif

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
  std::getline(fin, str);
  std::getline(fin, str);
  p->omegam = std::stof(str);//  1. omega matter
  std::getline(fin, str);
  std::getline(fin, str);
  p->omegal = std::stof(str); //  2. omega lambda
  std::getline(fin, str);
  std::getline(fin, str);
  p->h0 = std::stof(str);     //  3. hubble constant in unit of 100
  std::getline(fin, str);
  std::getline(fin, str);
  p->wq = std::stof(str);              //  4. dark energy equation of state parameter
  std::getline(fin, str);
  std::getline(fin, str);
  p->nzs = std::stoi(str);             //  5. number of source redshifts to consider
  std::getline(fin, str);
  if(p->nzs<-1){
    // build differencial convergence maps
    p->lenses = 1;
    std:: cout << " number of zs " << p->nzs << std:: endl;
    std:: cout << " I will build differencial convergence maps " << std:: endl;
    p->nzs = -p->nzs;
  }else{
    if(fabs(p->nzs)<1e-4){
      std:: cout << " number of zs " << p->nzs << std:: endl;
      std:: cout << " check this out, I will STOP here!!! " << std:: endl;
      exit(1);
    }
    p->lenses = 0;
  }
  p->z1.resize(p->nzs);
  p->zs.resize(p->nzs);
  for(int i=0;i<p->nzs;i++){
    std::getline(fin, str);
    p->zs[i] = std::stof(str);              //  6. source redshift (1:nzs)
  }
  std::getline(fin, str);
  std::getline(fin, str);
  p->nx = std::stoi(str);              // 7. number of pixels in x
  std::getline(fin, str);
  std::getline(fin, str);
  p->ny = std::stoi(str);              // 8. number of pixels in y
  std::getline(fin, str);
  std::getline(fin, str);
  p->fx = std::stof(str);              // 9. size of the field of view in x
  std::getline(fin, str);
  std::getline(fin, str);
  p->fy = std::stof(str);              // 10. size of the field of view in y
  std::getline(fin, str);
  std::getline(fin, p->filsigma);      // 11. file with lm s relation to be adopted
  std::getline(fin, str);
  std::getline(fin, p->cMrelation);     // 12. c-M relation
  std::getline(fin, str);
  std::getline(fin, str);
  p->sigmalnC = std::stof(str);        // 13. scatter in the c-M relation if a c-M model is used
  std::getline(fin, str);
  std::getline(fin, p->halodef);   // 14. vir use FOF haloes
  std::getline(fin, str);
  std::getline(fin, p->PinocchioFile );   // 15. Pinocchio PLC file, needed for Pincchio and BegognaCAT

  std:: cout << "  " << std:: endl;
  std:: cout << "  ___________________________________________"  << std:: endl;
  std:: cout << "  running with the following paramters:  " << std:: endl;
  std:: cout << "  Omegam           = " << p->omegam << std:: endl;
  std:: cout << "  Omegal           = " << p->omegal << std:: endl;
  std:: cout << "  Hubble parameter = " << p->h0 << std:: endl;
  std:: cout << "  w                = " << p->wq << std:: endl;
  if(p->nzs==1)
    std:: cout << "  source redshift  = " << p->zs[0] << std:: endl;
  else
    std:: cout << "  source redshifts  = ";

  for(int i=0;i<p->nzs;i++){
    std:: cout << p->zs[i] << "  ";
    if(i>0){
	    if(p->lenses==1)
	      p->z1[i] = p->zs[i-1];
	    else
	      p->z1[i] = 0;

    }else{
      p->z1[i] = 0;
    }
  }
  std:: cout << "  " << std:: endl;

  std:: cout << "  path catalogues  = " << p->pathcats << std:: endl;
  std:: cout << "  nx = " << p->nx << "  ny = " << p->ny << std:: endl;
  std:: cout << "  field of view: x = " << p->fx << "  y = " << p->fy << std:: endl;
  std:: ifstream fillms;
  fillms.open(p->filsigma.c_str());
  if(fillms.is_open()){
    std:: cout << "  file lm s rel    = " << p->filsigma << std:: endl;
  }else{
    std:: cout << " no file for the lm s relation ... I will use " << p->cMrelation << std:: endl;

  }
  std:: cout << "  scatter lnC      = " << p->sigmalnC << std:: endl;
  std:: cout << " I will build the map assuming "<< p->halodef <<" masses " << std:: endl;
  std:: cout <<  " Effective convergence map will be built using Pinocchio PLC output file " << std:: endl;
  std:: cout << p->PinocchioFile << std:: endl;
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
  std::unique_ptr<FITS> fout(new FITS(filename,FLOAT_IMG,naxis,naxes));
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
