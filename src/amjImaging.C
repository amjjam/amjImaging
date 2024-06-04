#include "../include/amjImaging.H"

namespace amjImaging{

  int fwrite(FILE *fp, const amjLos &los){
    int n=fwrite(fp,los.p());
    n+=fwrite(fp,los.d());
    return n;
  }
  

  template<>
  int fwrite(FILE *fp, const amjImage<amjLos> &image){
    unsigned int nR,nC;
    image.size(nR,nC);
    int n=fwrite(&nR,sizeof(int),1,fp);
    n+=fwrite(&nC,sizeof(int),1,fp);
    for(unsigned int iR=0;iR<nR;iR++)
      for(unsigned int iC=0;iC<nC;iC++)
	n+=fwrite(fp,image[iR][iC]);
    return n;
  }

  int fread(FILE *fp, amjLos &los){
    int n=fread(fp,los.p());
    n+=fread(fp,los.d());
    return n;
  }
  
  template<>
  int fread(FILE *fp, amjImage<amjLos> &image){
    int nR,nC;
    int n=fread(&nR,sizeof(int),1,fp);
    n+=fread(&nC,sizeof(int),1,fp);
    image.resize(nR,nC);
    for(int iR=0;iR<nR;iR++)
      for(int iC=0;iC<nC;iC++)
	n+=fread(fp,image[iR][iC]);
    return n;
  }
  
}
