/***********************************************************/
/*         PinMoka | check this out in the Makefile!       */
/* Constructed on top of WeakLMoka code from Carlo Giocoli */
/*                                                         */
/***********************************************************/
#include "cMZhao.h"
#include "cMBhattacharya.h"
#include "power2D.h"
#include "util.h"

const gsl_rng_type * Th;
gsl_rng * rh; // host halo concentration
const gsl_rng_type * Thfof;
gsl_rng * rhfof; // host halo concentration fof
int lenses;
double zhalo;

int main(int argc, char** argv){

  if (argc!=3){
    std::cout<< "No params and Compute Pk option!! I'll Stop here!!" <<std::endl;
    return -1;
  }

  string inifile=argv[1];
  bool computePk = bool (argv[2]);
  zhalo=0;
  time_t start;
  time (&start);

  Th = gsl_rng_default;
  Thfof = gsl_rng_default;
  rh = gsl_rng_alloc (Th);
  rhfof = gsl_rng_alloc (Thfof);
  gsl_rng_set(rh,1234);
  gsl_rng_set(rhfof,1234);
  showStart();
  struct InputParams p;
  readInput(&p, inifile);
  bool twohundredc;
  if(p.halodef=="200c")
    twohundredc=true;
  else
    twohundredc=false;
  time_t tread;
  time (&tread);
  std:: cout << " time to initialize and read the input file " << difftime(tread,start) << " sec" << std:: endl;
  //************** C O S M O L O G Y **********************//
  //**************   using MOKA lib  **********************//
  cosmology co(p.omegam,p.omegal,p.h0,p.wq);
  //*******************************************************//
  time (&tread);
  std:: cout << " time to initialize cosmological model " << difftime(tread,start) << " sec" << std:: endl;
  struct PlaneList pl;
  struct PinocchioPLC pinPLC;
  readPinocchio(&pinPLC,p);
  time (&tread);
  std:: cout << " time to read data " << difftime(tread,start) << " sec" << std:: endl;
  std:: vector<double> lm,s,ls;
  std:: ifstream fillms;
  if(p.cMrelation=="Zhao"){
    fillms.open(p.filsigma.c_str());
    if(fillms.is_open()){
      double a,b;
      while(fillms>> a >> b){
        lm.push_back(a);
        s.push_back(b);
        ls.push_back(log10(b));
      }
      iniTables(&co,lm,ls);
    }
  }
  // loop on haloes with z<max zs
  double zsmax = *max_element(p.zs.begin(),p.zs.end());
  double dsmax = co.angularDist(0.0,zsmax);
  std:: vector<double> dlens(p.nzs);
  // create the pinPLC.concentration array
  for(int i=0;i<pinPLC.nhaloes;i++){
    // compute the halo concentration
    // use the c-M relation models
    double c;
    // means have read the lm s file set in the input paramater file
    if(p.cMrelation=="Zhao"){
	    // Zhat+09 with the parameters of Giocoli+13
	    // we assign the redshift of the snap for the c-M
	    c = getCZhao(&co,pinPLC.mass[i],pinPLC.redshift[i],p.halodef); // if set use virial definition
    }else if(p.cMrelation=="Neto"){
	    //Neto+08 cM relation evolving with z as in Bullock+10
	    // we assign the redshift of the snap for the c-M
	    c = getCNeto(pinPLC.mass[i],pinPLC.redshift[i]);
    }else if(p.cMrelation=="Bhattacharya"){
      //Bhattacharya+13 cM relation
      // we assign the redshift of the snap for the c-M
      c = getCBhattacharya(&co,pinPLC.mass[i],pinPLC.redshift[i],p.halodef);
    }else{
      std::cout << "c-M relation "<< p.cMrelation << " not valid!";
      exit(1);
    }
    if(p.sigmalnC>0.){
	    // add a normal scatter!!!
      if(p.cMrelation=="Zhao"){

  	    pinPLC.concentration.push_back(gsl_ran_lognormal(rh,log(c),p.sigmalnC));

      }else if(p.cMrelation=="Neto"){

  	    pinPLC.concentration.push_back(gsl_ran_lognormal(rh,log(c),p.sigmalnC));

      }else if(p.cMrelation=="Bhattacharya"){

        pinPLC.concentration.push_back( c+gsl_ran_gaussian(rh,p.sigmalnC));
	    }

    }else{

	    pinPLC.concentration.push_back(c);

    }
  }

  std:: cout << " done with the concentrations " << std:: endl;
  time (&tread);
  std:: cout << " until now I spent " << difftime(tread,start) << " sec" << std:: endl;

  double fxrad = p.fx*M_PI/180.;
  double fyrad = p.fy*M_PI/180.;
  p.fxmin = -fxrad/2.;
  p.fxmax =  fxrad/2.;
  p.fymin = -fyrad/2.;
  p.fymax =  fyrad/2.;
  // start to build the map
  //nzs is the number of different redshift sources
  std:: valarray<std:: valarray<float> > kappa(p.nzs);
  for(int i=0;i<p.nzs;i++){
    kappa[i].resize(p.nx*p.ny);
  }
  fill_linear(p.x,p.nx,-fxrad/2.,fxrad/2.);
  fill_linear(p.y,p.ny,-fyrad/2.,fyrad/2.);
  int nh;
  double ra_c=0.;
  double dec_c=0.;
  nh = pinPLC.nhaloes;

  time_t loop;
  time (&loop);
  std:: vector<double> Radii(nh);
  std::srand( time(0) );
  for(int i=0;i<nh;i++){

    halo *ha;
    nfwHalo *nha;
    nfwLens *nLha;

    double mass, concentration, redshift;
    double zli;
    double ra,dec;

    mass = pinPLC.mass[i];
    concentration = pinPLC.concentration[i];
    zli=pinPLC.redshift[i];
    redshift=zli;
    ra=pinPLC.ra[i];
    dec=pinPLC.dec[i];

    time (&tread);

    if(mass>0 && concentration>0 && zli<zsmax){

      // Creating a pointer to an instance of a halo with mass 'mass'
      //redshift 'redshift' assuming a virial halo if 'twohundrec' == False
      // or '200c' otherwise
      ha = new halo(&co,mass,redshift,twohundredc);
      // Creating a pointer to an instance of a nfwhalo
      // with concentration 'concentration'
      nha = new nfwHalo(*ha,concentration);
      // Creating a pointer to an instance of a nfwlens
      // with the redshift of source 'zsmax'
      nLha = new nfwLens(*nha,zsmax);

      double rs = nLha->getScaleRadius();
      double Radius = nLha->getRvir();
      Radii[i] = Radius;
	    // is in unit of the scale radius
	    double Rz = Radius/rs;
	    int xi0=locate(p.x,ra);
	    xi0=GSL_MIN( GSL_MAX( xi0, 0 ), p.nx-1 );
	    int yi0=locate(p.y,dec);
	    yi0=GSL_MIN( GSL_MAX( yi0, 0 ), p.ny-1 );
	    double Dl = (co.angularDist(0.,redshift)*co.cn.lightspeed*1.e-7);
      // Computing the square box around the halo angular positions "ra" and "dec"
      // that will be filled according to nfw profile
      double dr = Radius/Dl;
	    p.xmin = ra - dr;
	    p.xmax = ra + dr;
	    p.ymin = dec - dr;
	    p.ymax = dec + dr;
	    int imin=0, imax=p.nx-1, jmin=0, jmax=p.ny-1;

	    locateHalo(p,ra,dec,dr,imin,imax,jmin,jmax);

	    if(imin<0) imin =0;
	    if(jmin<0) jmin =0;
	    for(int ij=0;ij<p.nzs;ij++){
	      dlens[ij] = (co.angularDist(zli,p.zs[ij])/co.angularDist(0.0,p.zs[ij]))/
                                      (co.angularDist(zli,zsmax)/dsmax);
	    }
	    for(int jj=jmin; jj<=jmax; jj++ )
        for(int ii=imin; ii<=imax; ii++ ){
	        double dx=(p.x[ii]-ra);
	        double dy=(p.y[jj]-dec);
	        double r=sqrt( dx*dx+dy*dy );
	        if(r<=dr){
	           for(int ij=0;ij<p.nzs;ij++){
		             // this need to be in Mpc/h physical!!!
		             if(zli>=p.z1[ij] && zli<p.zs[ij]){
		               kappa[ij][ii+p.nx*jj]+=(nLha->kappaz(Dl*r/rs,Rz)*dlens[ij]);
		             }
		         }
	         }
	      }
	    delete ha;
	    delete nha;
	    delete nLha;
    }
  }

  for(int ij=0;ij<p.nzs;ij++){
    // One should not enforce <kappa>=0!
    /*float km = kappa[ij].sum()/(p.nx*p.ny);
    for(int jj=0; jj<p.ny; jj++ ) for(int ii=0; ii<p.nx; ii++ ){

	     kappa[ij][ii+p.nx*jj] = (kappa[ij][ii+p.nx*jj]-km);
    }*/
    std:: string szid;
    std:: ostringstream oszid;
    oszid << ij;
    szid = oszid.str();
    size_t name_begin,name_end;
    name_begin=p.PinocchioFile.find_last_of("/")+11;
    name_end=p.PinocchioFile.find_last_of(".")-name_begin;
    std:: string plcsufix = p.PinocchioFile.substr(name_begin,name_end);
    std:: string filout = "!kappa_" + plcsufix + "_" + szid + ".fits";
    writeFits(filout,kappa[ij],p.nx,p.ny,p.zs[ij],fxrad,fyrad,ra_c,dec_c);

    if(computePk){

      int nb = 128;
      double *ll;
      double *Pl;
      ll=new double[nb];
      Pl=new double[nb];
      powerl(kappa[ij],kappa[ij],p.nx,p.ny,fxrad,fyrad,ll,Pl,nb);
      std:: ofstream outfile;
      outfile.open("mapPowerSpectrum_" + plcsufix + "_" + szid+ ".dat");
      for(int i=0;i<nb;i++){
        outfile << ll[i] << "  " << Pl[i] << std::endl;
      }

      outfile.close();

    }

  }
  std:: cout << " end of work ... ;-) " << std:: endl;
  return 0;
}
