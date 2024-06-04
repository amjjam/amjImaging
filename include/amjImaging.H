#ifndef _AMJIMAGING_H_
#define _AMJIMAGING_H_

#include <aVec.H>
#include <vector>

namespace amjImaging{
  template<class T>
  class amjImage1{
  public:
    amjImage1(){};
    amjImage1(int N){_image.resize(N);};
    virtual ~amjImage1(){};
    template<class U>
    amjImage1<T> &operator=(const amjImage1<U> &);
    unsigned int size() const {return _image.size();};
    void resize(int N){_image.resize(N);};
    T &operator[](int i){return _image[i];}
  protected:
    std::vector<T> _image;
  };
  
  template<class T> template<class U>
  amjImage1<T> &amjImage1<T>::operator=(const amjImage1<U> &other){
    if(this==&other)
      return *this;
    _image.resize(other._image.size());
    for(int i=0;i<_image.size();i++)
      _image[i]=other._image[i];
    return *this;
  }
  
  template<class T>
  class amjImage: public amjImage1<T>{
  public:
    amjImage():amjImage1<T>(),nR(0),nC(0){};
    amjImage(int nR, int nC):amjImage1<T>(nR*nC),nR(nR),nC(nC){};
    template<class U>
    amjImage(amjImage<U> &image){resize(image);};
    virtual ~amjImage(){};
    template<class U>
    amjImage<T> &operator=(const amjImage<U> &);  
    void size(unsigned int &_nR, unsigned int &_nC) const {_nR=nR; _nC=nC;};
    void resize(unsigned int _nR, unsigned int _nC){nR=_nR; nC=_nC;
      amjImage1<T>::_image.resize(nR*nC);};
    template<class U>
    void resize(amjImage<U> &image){image.size(nR,nC);
      amjImage1<T>::_image.resize(nR*nC);};
    T *operator[](int i){return &amjImage1<T>::_image[i/nC];}
    const T *operator[](int i) const {return &amjImage1<T>::_image[i/nC];}
  protected:
    unsigned int nR,nC;
  };
  
  template<class T> template<class U>
  amjImage<T> &amjImage<T>::operator=(const amjImage<U> &other){
    if(this==&other)
      return *this;
    amjImage1<T>::operator=(other);
    nR=other.nR;
    nC=other.nC;
  }
  
  template<class T>
  class amjValue{
  public:
    amjValue(){};
    virtual ~amjValue(){};
    virtual T value(const aVec &) const =0;
  };
  
  class amjLos{
  public:
    amjLos(){};
    amjLos(aVec p, aVec d):_p(p),_d(unit(d)){};
    virtual ~amjLos(){};
    void p(const aVec &p_){_p=p_;};
    void d(const aVec &d_){_d=unit(d_);};
    aVec &p(){return _p;};
    const aVec &p() const {return _p;};
    aVec &d(){return _d;};
    const aVec &d() const {return _d;};
  private:
    aVec _p;
    aVec _d;
  };

  template<class T>
  class amjRange{
  public:
    amjRange(){};
    virtual ~amjRange(){};
    virtual void range(const amjLos &los, T &start, T &stop) const =0;
  };
  
  template<class T>
  class amjIntegrator{
  public:
    amjIntegrator(){};
    virtual ~amjIntegrator(){};
    virtual T integrate(amjLos &los, amjRange<T> &range, amjValue<T> &value)=0;
  };
  
#include <gsl/gsl_integration.h>
  
  template<class T>
  class amjIntegrator_gsl_qag: public amjIntegrator<T>{
  public:
    amjIntegrator_gsl_qag(size_t nWorkspace=1000):nWorkspace(nWorkspace){
      workspace=gsl_integration_workspace_alloc(nWorkspace);
      seteps(0.1,0.1);
    }
    virtual ~amjIntegrator_gsl_qag(){
      gsl_integration_workspace_free(workspace);
    }
    virtual T integrate(amjLos &los, amjRange<T> &range, amjValue<T> &value);
    struct params{
      aVec p;
      aVec d;
      amjValue<T> *value;
    };
    void seteps(double _epsabs, double _epsrel){epsabs=_epsabs;epsrel=_epsrel;};
  private:
    size_t nWorkspace;
    gsl_integration_workspace *workspace;
    double epsabs,epsrel;
  };
  
  template<class T>
  T amjIntegrator_gsl_qag<T>::integrate(amjLos &los, amjRange<T> &range,
					amjValue<T> &value){
    gsl_function F;
    struct params p={los.p(),unit(los.d()),&value};
    F.function=[](double x, void *pp)->double{
      params *p=(params *)pp; return p->value->value(p->p+p->d*x);
    };
    F.params=&p;
    double start,stop,sum,abserr;
    range.range(los,start,stop);
    int error;
    if((error=gsl_integration_qag(&F,start,stop,epsabs,epsrel,nWorkspace,5,
				  workspace,&sum,&abserr))>0)
      std::cout << "gsl_integration_qag return " << error << std::endl;
    return sum;
  }

  int fwrite(FILE *fp, const amjLos &los);
  int fread(FILE *fp, amjLos &los);

  template<class T>
  int fwrite(FILE *fp, const amjImage<T> &image){
    unsigned int nR,nC;
    image.size(nR,nC);
    int n=fwrite(&nR,sizeof(int),1,fp);
    n+=fwrite(&nC,sizeof(int),1,fp);
    n+=fwrite(image[0],sizeof(T),nR*nC,fp);
    return n;
  }

  template<>
  int fwrite(FILE *fp, const amjImage<amjLos> &image);

  template<class T>
  int fread(FILE *fp, amjImage<T> &image){
    int nR,nC;
    int n=fread(&nR,sizeof(int),1,fp);
    n+=fread(&nC,sizeof(int),1,fp);
    image.resize(nR,nC);
    n+=fread(image[0],sizeof(T),nR*nC,fp);
    return n;
  }

  template<>
  int fread(FILE *fp, amjImage<amjLos> &image);
}

#endif
