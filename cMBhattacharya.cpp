#include "cMBhattacharya.h"

double getCBhattacharya(cosmology *cosmology, double M, double z,std:: string halodef){

  double a = 1./(1+z);
  double D = a*cosmology->growthFactor(a)/cosmology->growthFactor(1.0);
  double nu = 1.0/D*( 1.12*pow(M/5.0/1e13,0.3) + 0.53 );
  double alpha, beta, c;

  if(halodef.compare("vir")){
    alpha = 0.9;
    beta=-0.29;
    c=7.7;
  }else if(halodef.compare("200c")){
    alpha = 0.54;
    beta=-0.35;
    c=5.9;
  }else if(halodef.compare("200m")){
    alpha = 1.15;
    beta=-0.29;
    c=9.0;
  }else{
    std::cout<< "Halo definition "<< halodef << " not valid " <<std::endl;
    exit(1);
  }

  return c*pow(D,alpha)*pow(nu,beta);

}
