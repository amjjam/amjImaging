#ifndef _AMJIMAGING_H_
#define _AMJIMAGING_H_

#include <aVec.H>
#include <vector>

namespace amjImaging{
  template<class T>
  class Image1{
  public:
    Image1(){};
    Image1(int N){_image.resize(N);};
    virtual ~Image1(){};
    template<class U>
    Image1<T> &operator=(const Image1<U> &);
    unsigned int size() const {return _image.size();};
    void resize(int N){_image.resize(N);};
    T &operator[](int i){return _image[i];}
    const T &operator[](int i) const {return _image[i];}
  protected:
    std::vector<T> _image;
  };
  
  template<class T> template<class U>
  Image1<T> &Image1<T>::operator=(const Image1<U> &other){
    if(this==&other)
      return *this;
    _image.resize(other._image.size());
    for(int i=0;i<_image.size();i++)
      _image[i]=other._image[i];
    return *this;
  }
  
  template<class T>
  class Image: public Image1<T>{
  public:
    Image():Image1<T>(),nR(0),nC(0){};
    Image(int nR, int nC):Image1<T>(nR*nC),nR(nR),nC(nC){};
    template<class U>
    Image(Image<U> &image){resize(image);};
    virtual ~Image(){};
    template<class U>
    Image<T> &operator=(const Image<U> &);
    unsigned int size(){return Image1<T>::size();};
    void size(unsigned int &_nR, unsigned int &_nC) const {_nR=nR; _nC=nC;};
    void resize(unsigned int _nR, unsigned int _nC){nR=_nR; nC=_nC;
      Image1<T>::_image.resize(nR*nC);};
    template<class U>
    void resize(Image<U> &image){image.size(nR,nC);
      Image1<T>::_image.resize(nR*nC);};
    //T &operator[](int i){return &Image1<T>::_image[i/nC];}
    //const T &operator[](int i) const {return &Image1<T>::_image[i/nC];}
  protected:
    unsigned int nR,nC;
  };
  
  template<class T> template<class U>
  Image<T> &Image<T>::operator=(const Image<U> &other){
    if(this==&other)
      return *this;
    Image1<T>::operator=(other);
    nR=other.nR;
    nC=other.nC;
  }

  // This is a class template which returns a value at a vector
  // point. A model will inherit from this and define the function
  // value
  template<class T>
  class Value{
  public:
    Value(){};
    virtual ~Value(){};
    virtual T value(const aVec &) const =0;
  };
  
  // This is a class which defines a line-of sight which is useful for
  // intgration along a ray in 3D. It is a pair of vectors, a starting
  // location and a direction. The return value of d() is always a
  // normalized vector (assuming the set direction is not zero
  // length).
  class Los{
  public:
    Los(){};
    Los(aVec p, aVec d):_p(p),_d(unit(d)){};
    virtual ~Los(){};
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

  // This struct template defines a range of values. It is used to
  // store the limits of a line-of-sight integration.
  template<class T>
  struct Range{
    T start;
    T stop;
  };


  // This is a class template which returns a vector of ranges. This
  // is typically inherited by a model, such that when the function
  // ranges is called with a line-of-sight, it populates a vector of
  // struct Range with the limits of segments that must be integrated
  // along to integrate along the line-of-sight for this model.
  template<class T>
  class Ranges{
  public:
    Ranges(){};
    virtual ~Ranges(){};
    virtual void ranges(const Los &los, std::vector<struct Range<T> > &) const =0;
  };

  // This is a class template for an integrator. It Takes a
  // line-of-sight, a Ranges class, and a Value class. It integrates
  // Value along the line-of-sight over the intervals returned by
  // Ranges for that line-of-sight
  template<class T>
  class Integrator{
  public:
    Integrator(){};
    virtual ~Integrator(){};
    virtual T integrate(Los &los, Ranges<T> &range, Value<T> &value)=0;
  };
  
#include <gsl/gsl_integration.h>

  // This ia a integrator with the integrate defined as a QAG
  // algorithm from GNU scientific laboratory
  template<class T>
  class Integrator_gsl_qag: public Integrator<T>{
  public:
    Integrator_gsl_qag(size_t nWorkspace=1000):nWorkspace(nWorkspace){
      workspace=gsl_integration_workspace_alloc(nWorkspace);
      seteps(0.1,0.1);
    }
    virtual ~Integrator_gsl_qag(){
      gsl_integration_workspace_free(workspace);
    }
    virtual T integrate(Los &los, Ranges<T> &range, Value<T> &value);
    struct params{
      aVec p;
      aVec d;
      Value<T> *value;
    };
    void seteps(double _epsabs, double _epsrel){epsabs=_epsabs;epsrel=_epsrel;};
  private:
    size_t nWorkspace;
    gsl_integration_workspace *workspace;
    double epsabs,epsrel;
  };
  
  template<class T>
  T Integrator_gsl_qag<T>::integrate(Los &los, Ranges<T> &ranges,
				     Value<T> &value){
    gsl_function F;
    struct params p={los.p(),unit(los.d()),&value};
    F.function=[](double x, void *pp)->double{
      params *p=(params *)pp; return p->value->value(p->p+p->d*x);
    };
    F.params=&p;
    std::vector<struct Range<double> > r;
    ranges.ranges(los,r);
    
    int error;
    double sum=0,subsum,abserr;
    for(unsigned int i=0;i<r.size();i++){
      if((error=gsl_integration_qag(&F,r[i].start,r[i].stop,epsabs,epsrel,
				    nWorkspace,5,workspace,&subsum,&abserr))>0)
	std::cout << "gsl_integration_qag return " << error << ": "
		  << "range (" << r[i].start << "," << r[i].stop << ") on p="
		  << los.p() << ", d=" << los.d() << std::endl;
      sum+=subsum;
    }
    return sum;
  }

  int fwrite(FILE *fp, const Los &los);
  int fread(FILE *fp, Los &los);

  template<class T>
  int fwrite(FILE *fp, const Image<T> &image){
    unsigned int nR,nC;
    image.size(nR,nC);
    int n=fwrite(&nR,sizeof(int),1,fp);
    n+=fwrite(&nC,sizeof(int),1,fp);
    n+=fwrite(&image[0],sizeof(T),nR*nC,fp);
    return n;
  }

  template<>
  int fwrite(FILE *fp, const Image<Los> &image);

  template<class T>
  int fread(FILE *fp, Image<T> &image){
    int nR,nC;
    int n=fread(&nR,sizeof(int),1,fp);
    n+=fread(&nC,sizeof(int),1,fp);
    image.resize(nR,nC);
    n+=fread(&image[0],sizeof(T),nR*nC,fp);
    return n;
  }

  template<>
  int fread(FILE *fp, Image<Los> &image);
}

#endif
