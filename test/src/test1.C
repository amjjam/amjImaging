
#include "../../include/amjPImaging.H"

using namespace amjImaging;

class TestModel: public amjValue<double>, public amjRange<double>{
public:
  TestModel(){};
  virtual ~TestModel(){};
  virtual double value(const aVec &p) const;
  virtual void range(const amjLos &los, double &start, double &stop) const;
};

#include <cmath>
double TestModel::value(const aVec &p) const{
  if(fabs(p.X())<5&&fabs(p.Y())<5&&fabs(p.Z())<5)
    return 1;
  return 0;
}

void TestModel::range(const amjLos &los, double &start, double &stop) const{
  start=length(los.p())/length(los.d())-15;
  if(start<0)
    start=0;
  stop=length(los.p())/length(los.d())+15;  
}

class TestModel2: public TestModel{
public:
  TestModel2(){};
  virtual ~TestModel2(){};
  virtual double value(const aVec &p) const;
};

double TestModel2::value(const aVec &p) const{
  return 1.0/(1.0+exp(5*(-5.0-p.X())))/(1.0+exp(5*(-5.0-p.Y())))/(1.0+exp(5*(-5.0-p.Z())))/
    (1.0+exp(5*(p.X()-5)))/(1.0+exp(5*(p.Y()-5)))/(1.0+exp(5*(p.Z()-5)));
}

class TestCamera: public amjPCamera<double>{
public:
  TestCamera(double tR, double tC, aVec p, aVec d, double roll=0):
    tR(tR),tC(tC),p(p),d(d),rX(roll){
    rY=-asin(d.Z()/sqrt(d.X()*d.X()+d.Y()*d.Y()));
    rZ=atan2(d.Y(),d.X());
  }
  ~TestCamera(){};
  amjLos getlos(int);
private:
  double tR;
  double tC;
  aVec p;
  aVec d;
  double rX;
  double rY;
  double rZ;
};

amjLos TestCamera::getlos(int i){
  int iR=i/_nC;
  int iC=i%_nC;
  double altitude=tR/2-((double)iR+0.5)*tR/(double)_nR;
  double azimuth=tC/2-((double)iC+0.5)*tC/(double)_nC;
  aVec d(cos(azimuth)*cos(altitude),sin(azimuth)*cos(altitude),sin(altitude));
  d=d.rotateX(rX);
  d=d.rotateY(rY);
  d=d.rotateZ(rZ);
  //std::cout << "(" << i << "," << iR << "," << iC << ") " <<  p << ", " << d << std::endl;
  return amjLos(p,d);
}

class TestIntegrator: public amjIntegrator<double>{
public:
  TestIntegrator(){};
  virtual ~TestIntegrator(){};
  virtual double integrate(amjLos &los, amjRange<double> &range,
			   amjValue<double> &value);
};

double TestIntegrator::integrate(amjLos &los, amjRange<double> &range,
				 amjValue<double> &value){

  aVec p=los.p(),d=los.d();
  double start,stop;

  range.range(los,start,stop);

  double stepsize=0.01;
  double v=0;
  for(double s=start;s<stop;s+=stepsize){
    v+=value.value(p+d*s);
  }
  return v;
}

int main(int argc, char *argv[]){
  MPI_Init(&argc,&argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  TestModel2 model;
  int N=100,M=100;
  amjPImage<double> pimage(N,M,MPI_COMM_WORLD);
  TestCamera camera(45.0/180*M_PI,45.0/180*M_PI,
		    aVec(-20,-15,10),aVec(2,1.5,-1),0.25);
  //TestIntegrator integrator;
  amjIntegrator_gsl_qag<double> integrator;
  integrator.seteps(0.001,0.001);
  
  camera.snap(pimage,model,model,integrator);

  if(rank==0){
    amjImage<double> image;
    pimage.getimage(image);
    FILE *fp=fopen("../dat/test1.dat","w");
    fwrite(fp,image);
    fclose(fp);
  }
  
  return 0;
}
