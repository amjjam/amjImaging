#include "../include/amjImaging.H"

namespace amjImaging{

  int fwrite(FILE *fp, const Los &los){
    int n=fwrite(fp,los.p());
    n+=fwrite(fp,los.d());
    return n;
  }
  

  template<>
  int fwrite(FILE *fp, const Image<Los> &image){
    unsigned int nR,nC;
    image.size(nR,nC);
    int n=fwrite(&nR,sizeof(int),1,fp);
    n+=fwrite(&nC,sizeof(int),1,fp);
    for(unsigned int i=0;i<nR*nC;i++)
      n+=fwrite(fp,image[i]);
    return n;
  }

  int fread(FILE *fp, Los &los){
    int n=fread(fp,los.p());
    n+=fread(fp,los.d());
    return n;
  }
  
  template<>
  int fread(FILE *fp, Image<Los> &image){
    int nR,nC;
    int n=fread(&nR,sizeof(int),1,fp);
    n+=fread(&nC,sizeof(int),1,fp);
    image.resize(nR,nC);
    for(unsigned int i=0;i<image.size();i++)
	n+=fread(fp,image[i]);
    return n;
  }
  
}
