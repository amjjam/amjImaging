#ifndef _AMJPIMAGING_H_
#define _AMJPIMAGING_H_

#include <vector>
#include <mpi.h>
#include "amjImaging.H"

namespace amjImaging{
  template<class T>
  class amjPImage1{
  public:
    amjPImage1(unsigned int, MPI_Comm);
    virtual ~amjPImage1(){};
    unsigned int size(){return N;}
    unsigned int localsize(){return ranksize();}
    virtual int globalindex(int localindex);
    T &operator[](int i){return _image[i];};
    void gather();
    virtual void getimage(amjImaging::amjImage1<T> &);
  protected:
    MPI_Comm _comm;
  private:
    unsigned int ranksize(int r=-1);
    unsigned int maxranksize();
    unsigned int N; // Global number of pixels
    std::vector<T> _image; // vector of local pixels
    std::vector<T> _globalimage; // vector of all pixels (after gather on rank 0)
  };
  
  template<class T>
  amjPImage1<T>::amjPImage1(unsigned int _N, MPI_Comm comm_){
    MPI_Comm_dup(comm_,&_comm);
    int comm_size,comm_rank;
    MPI_Comm_size(_comm,&comm_size);
    MPI_Comm_rank(_comm,&comm_rank);
    N=_N;
    if(comm_rank==0)
      _image.resize(maxranksize()*comm_size);
    else
      _image.resize(maxranksize());
  }
  
  template<class T>
  int amjPImage1<T>::globalindex(int localindex){
    int comm_rank;
    MPI_Comm_rank(_comm,&comm_rank);
    int j=0;
    for(int i=0;i<comm_rank;i++)
      j+=ranksize(i);
    return j+localindex;
  }
  
  template<class T>
  void amjPImage1<T>::gather(){
    int comm_size,comm_rank;
    MPI_Comm_size(_comm,&comm_size);
    MPI_Comm_rank(_comm,&comm_rank);
    MPI_Gather((void*)&_image[0],maxranksize()*sizeof(T),MPI_BYTE,
	       (void*)&_image[0],maxranksize()*sizeof(T),MPI_BYTE,0,_comm);
    if(comm_rank==0){
      for(int irank=1,i=ranksize(irank);irank<comm_rank;irank++){
	memcpy(&_image[i],&_image[irank*maxranksize()],
	       ranksize(irank)*sizeof(T));
	i+=ranksize(irank);
      }
    }
  }
  
  // This function may only be called from rank=0. It will abort if it
  // is not called from rank=0
  template<class T>
  void amjPImage1<T>::getimage(amjImaging::amjImage1<T> &image){
    int comm_rank;
    MPI_Comm_rank(_comm,&comm_rank);
    if(comm_rank!=0){
      std::cout << "rank=" << comm_rank << ", may only be called from rank=0."
		<< std::endl;
      abort();
    }
    image.resize(N);
    memcpy(&image[0],&_image[0],N*sizeof(T));
  }
  
  template<class T>
  unsigned int amjPImage1<T>::ranksize(int rank){
    int comm_size;
    MPI_Comm_size(_comm,&comm_size);
    unsigned int n=N/comm_size;
    if(rank<0)
      MPI_Comm_rank(_comm,&rank);
    if(N-n*comm_size<(unsigned int)rank)
      n++;
    return n;
  }
  
  template<class T>
  unsigned int amjPImage1<T>::maxranksize(){
    int comm_size;
    MPI_Comm_size(_comm,&comm_size);
    unsigned int n=N/(unsigned int)comm_size;
    if(n*comm_size<N)
      n++;
    return n;
  }
  
  template<class T>
  class amjPImage: public amjPImage1<T>{
  public:
    amjPImage(unsigned int nR, unsigned int nC, MPI_Comm comm):
      amjPImage1<T>(nR*nC,comm),_nR(nR),_nC(nC){};
    virtual ~amjPImage(){};
    unsigned int nR() const {return _nR;};
    unsigned int nC() const {return _nC;};
    virtual void getimage(amjImaging::amjImage<T> &);
  private:
    unsigned int _nR;
    unsigned int _nC;
  };
  
  
  // This function may only be called from rank=0. It will abort if it
  // is not called from rank=0
  template<class T>
  void amjPImage<T>::getimage(amjImaging::amjImage<T> &image){
    int comm_rank;
    MPI_Comm_rank(this->_comm,&comm_rank);
    if(comm_rank!=0){
      std::cout << "rank=" << comm_rank << ", may only be called from rank=0."
		<< std::endl;
      abort();
    }
    image.resize(_nR,_nC);
    amjPImage1<T>::getimage(image);
  }
  
  template<class T>
  class amjPCamera1{
  public:
    amjPCamera1(){};
    virtual ~amjPCamera1(){};
    void snap(amjPImage1<T> &, amjImaging::amjValue<T> &,
	      amjImaging::amjRange<T> &,
	      amjImaging::amjIntegrator<T> &);
    void snap(amjImage1<amjLos> &);
    virtual amjLos getlos(int gi)=0;
  };
  
  template<class T>
  void amjPCamera1<T>::snap(amjPImage1<T> &image, amjValue<T> &value,
			    amjRange<T> &range, amjIntegrator<T> &integrator){
    amjLos los;
    for(unsigned int i=0;i<image.localsize();i++){
      los=getlos(image.globalindex(i));
      image[i]=integrator.integrate(los,range,value);
    }
    image.gather();
  }
  
  template<class T>
  void amjPCamera1<T>::snap(amjImage1<amjLos> &image){
    for(unsigned int i=0;i<image.size();i++)
      image[i]=getlos(i);
  }
  
  template<class T>
  class amjPCamera: public amjPCamera1<T>{
  public:
    amjPCamera(){};
    virtual ~amjPCamera(){};
    virtual void snap(amjPImage<T> &, amjValue<T> &, amjRange<T> &,
		      amjIntegrator<T> &);
    virtual void snap(amjImage<amjLos> &);
  protected:
    unsigned int _nR;
    unsigned int _nC;
  };
  
  template<class T>
  void amjPCamera<T>::snap(amjPImage<T> &image, amjValue<T> &value,
			   amjRange<T> &range, amjIntegrator<T> &integrator){
    _nR=image.nR();
    _nC=image.nC();
    amjPCamera1<T>::snap(image,value,range,integrator);
  }
  
  template<class T>
  void amjPCamera<T>::snap(amjImage<amjLos> &image){
    image.size(_nR,_nC);
    amjPCamera1<T>::snap(image);
  }

}

#endif
