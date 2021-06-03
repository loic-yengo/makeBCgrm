// STD include
#include <string.h>
#include <math.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>

using namespace std;

// Eigen include
#include "Eigen/Dense"

using namespace Eigen;

float get_M_snp(string grmNfile, int n, ofstream &fileLog){
  ifstream N_bin(grmNfile.c_str(), ios::in | ios::binary);
  if(!N_bin){
    cerr<<"\t[readGRM] Error reading file "<<grmNfile<<endl;
    exit(1);
  }
  cout << "\tReading the number of SNPs for the GRM from [" + grmNfile + "].\n";
  fileLog<< "\tReading the number of SNPs for the GRM from [" + grmNfile + "].\n";
  int size = sizeof (float);
  float f_buf   = 0.0;
  float M_snp   = 0.0;
  float m_snp   = 0.0;
  float n_pairs = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      if (!(N_bin.read((char*) &f_buf, size))) throw (("Error: the size of the [" + grmNfile + "] file is incomplete?").c_str());
      if(f_buf>M_snp){
         M_snp = f_buf;
      }
      m_snp += f_buf;
      n_pairs += 1.;
    }
  }
  N_bin.close();
  m_snp = m_snp / n_pairs;
  cout << "\tMax number of SNPs: "<<M_snp<<endl;
  cout << "\tMean number of SNPs: "<<m_snp<<endl;

  fileLog<< "\t\tMax number of SNPs: "<<M_snp<<endl;
  fileLog<< "\t\tMean number of SNPs: "<<m_snp<<endl;


  return M_snp;
}

void readGRM(string grmBinfile, MatrixXf &GRM,int n, ofstream &fileLog){
  ifstream binData(grmBinfile.c_str(), ios::in | ios::binary);
  if(!binData){
    cerr << "\t[readGRM] Error reading file "<<grmBinfile<<endl;
    exit(1);
  }
  cout << "\n\tReading the GRM from [" + grmBinfile + "]." << endl;
  fileLog << "\n\tReading the GRM from [" + grmBinfile + "]." << endl;

  int size = sizeof (float);
  float f_buf  = 0.;
  float Mx_d   = 0.;
  float Mx     = 0.;
  float Mx2    = 0.;
  float nPairs = 0.; 
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      if (!(binData.read((char*) &f_buf, size))) throw (("\tError: the size of the [" + grmBinfile + "] file is incomplete?").c_str());
      GRM(i, j) = f_buf;
      if(j<i){
        GRM(j, i) = f_buf;
        Mx       += f_buf;
        Mx2      += f_buf * f_buf;
        nPairs   += 1.;
      }else{
        Mx_d += f_buf;
      }
    }
  }
  
  Mx   = Mx  / nPairs;
  Mx2  = Mx2 / nPairs;
  Mx_d = Mx_d / n;
  float varOffDiag = Mx2-Mx*Mx;
  //float Me = 2.0/varOffDiag;
  fileLog<<"\t\tMean of diagonal elements: "<<Mx_d<<".\n";
  fileLog<<"\t\tMean of off-diagonal elements: "<<Mx<<".\n";
  fileLog<<"\t\tVariance of diagonal elements: "<<varOffDiag<<".\n";
  //fileLog<<"\t\tEstimated Me = "<<Me<<"."<<endl;
  
  binData.close();
}

int main(int argc, char *argv[]){
  
  int i;
  string sw;
  
  if(argc==1){
    cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
    exit(1);
  }
  
  // Read arguments
  sw = argv[1];
  if (sw == "--help"){
    cerr<<"\t--mgrm    : list of genetic relationship matrices (GRM)."<<endl;
    cerr<<"\t--out     : Specify prefix for output GRM."<<endl;
    //cerr<<"\t--silent  : Specify whether the different steps of the calculations should be displayed on screen."<<endl;
    exit(1);
  }else{
    if (argc == 1) {
      cerr<<"\tArguments must be specified. Type --help for more details."<<endl;
      exit(1);
    }
  }
  
  string grmBinFile;
  string grmIdFile;
  string grmNFile;
  
  string grmPrefix    = "";
  string grmOutPrefix = "";
  string mgrmFile     = "";
  
  for(i = 1; i<argc;i++){
    sw = argv[i];
    if (sw == "--mgrm"){
      mgrmFile = argv[i + 1];
    }
    if (sw == "--out"){
      grmOutPrefix = argv[i + 1];
    }
  }
  
  string grmOutBinFile_w = grmOutPrefix+".within.grm.bin";
  string grmOutIdFile_w  = grmOutPrefix+".within.grm.id";

  string grmOutBinFile_b = grmOutPrefix+".between.grm.bin";
  string grmOutIdFile_b  = grmOutPrefix+".between.grm.id";

  string logFile = grmOutPrefix+".makeBCgrm.log";
  ofstream fileLog(logFile.c_str());

  clock_t tic = clock();
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    //cout <<"\t[makeBCgrm - L. Yengo]\n";
    cout <<"\tAnalysis starts : ";
    cout << (now->tm_year + 1900) << '-'
         << (now->tm_mon + 1) << '-'
         <<  now->tm_mday << " at "
         <<  now->tm_hour <<":"<<now->tm_min<<":"<<now->tm_sec
         <<  ".\n";

  string line = "";
  ifstream mgrmStream;
  int nGRM = -1;
  int n    =  0;
  mgrmStream.open(mgrmFile.c_str());
  while(mgrmStream){
    getline(mgrmStream,grmPrefix);
    if(nGRM==-1){
      grmIdFile = grmPrefix+".grm.id";
      ifstream idStream;
      ofstream fileId_w(grmOutIdFile_w.c_str());
      ofstream fileId_b(grmOutIdFile_b.c_str());
      idStream.open(grmIdFile.c_str());
      while(idStream){
        getline(idStream,line);
        if(line!=""){
          fileId_w<<line<<endl;
          fileId_b<<line<<endl;
          n++;
        }
      }
      idStream.close();
      fileId_w.close();
      fileId_b.close();
      cout<<"\t"<<n<<" samples detected."<<endl;
      fileLog<<"\t"<<n<<" samples detected."<<endl;
    }
    nGRM++;
  }
  mgrmStream.close();
  cout<<"\t"<<nGRM<<" GRM detected."<<endl;
  fileLog<<"\t"<<nGRM<<" GRM detected."<<endl;

  float *nsnps = new float[nGRM];
  float M_snps = 0.0;

  // Allocate GRMs
  MatrixXf G_all     = MatrixXf::Zero(n,n);
  MatrixXf G_tmp     = MatrixXf::Zero(n,n);
  MatrixXf G_between = MatrixXf::Zero(n,n);
  MatrixXf G_within  = MatrixXf::Zero(n,n);
 
  initParallel();
  int nThreads = nbThreads( );
  cout<<"\tDefault #Threads used is "<<nThreads<<".\n\n";

  mgrmStream.open(mgrmFile.c_str());
  for(int iGRM=0;iGRM<nGRM;iGRM++){
    /*
    clock_t tic = clock();
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    cout <<"\tStarting reading genotypes : ";
    cout << (now->tm_year + 1900) << '-' 
         << (now->tm_mon + 1) << '-'
         <<  now->tm_mday << " at "
         <<  now->tm_hour <<":"<<now->tm_min<<":"<<now->tm_sec
         <<  ".\n";
    */
    getline(mgrmStream,grmPrefix);
    grmBinFile = grmPrefix+".grm.bin";
    readGRM(grmBinFile,G_tmp,n,fileLog);
    
    grmNFile    = grmPrefix+".grm.N.bin";
    nsnps[iGRM] = get_M_snp(grmNFile,n,fileLog);
    M_snps     += nsnps[iGRM];

    G_all      += nsnps[iGRM] * G_tmp;
    G_within   += nsnps[iGRM] * nsnps[iGRM] * G_tmp * G_tmp;
    
    /*
    time_t t2 = time(0);   // get time now
    struct tm * now2 = localtime( & t2 );
    cout <<"\tAnalysis ends: ";
    cout << (now2->tm_year + 1900) << '-' 
         << (now2->tm_mon + 1) << '-'
         <<  now2->tm_mday << " at "
         <<  now2->tm_hour <<":"<<now2->tm_min<<":"<<now2->tm_sec
         << ".\n";
    clock_t toc = clock();
    float time_elapsed = (float)(toc - tic) / CLOCKS_PER_SEC;    
    printf("\tTime elapsed: %f seconds.\n\n", time_elapsed / nThreads);
    */
  }
  mgrmStream.close();
  delete [] nsnps;   
  
  cout<<"\n\tTotal number of SNPs across all GRMs is "<<M_snps<<".\n";
  fileLog<<"\n\tTotal number of SNPs across all GRMs is "<<M_snps<<".\n";

  G_between        = G_all * G_all - G_within;
  float Gb_bar     = G_between.trace() / n;
  float Gw_bar     = G_within.trace()  / n;
  float Ga_bar     = G_all.trace()     / n;
  G_between        = (1.0 / Gb_bar) * G_between;
  G_within         = (1.0 / Gw_bar) * G_within;
  G_all            = (1.0 / Ga_bar) * G_all;

  cout<<"\tTrace of (unscaled) GRM: "<<Ga_bar<<".\n\n";
  fileLog<<"\tTrace of (unscaled) GRM: "<<Ga_bar<<".\n\n";

  //printf("\tTrace of G_all: %f.\n\n",Ga_bar);

  cout<<"\n\tG_within.\n";
  for(int j=0;j<5;j++){
    for(int k=0;k<5;k++){
      cout<<"\t"<<G_all(j,k);
    }
    cout<<endl;
  }
  cout<<endl;
  
  cout<<"\n\tG_between.\n";
  for(int j=0;j<5;j++){
    for(int k=0;k<5;k++){
      cout<<"\t"<<G_between(j,k);
    }
    cout<<endl;
  }
  cout<<endl;

  // Output grm
  int j;
  float f_buf = 0.0;
  int size = sizeof (float);
   
  fstream Aw_Bin(grmOutBinFile_w.c_str(), ios::out | ios::binary);
  if (!Aw_Bin) throw (("Error: can not open the file [" + grmOutBinFile_w + "] to write.").c_str());
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      f_buf = G_all(i, j); // within = standard GRM
      Aw_Bin.write((char*) &f_buf, size);
    }
  }
  Aw_Bin.close();
  cout << "\tWithin-Chromosome GRM of " << n << " individuals has been saved in the file [" + grmOutBinFile_w + "] (in binary format)." << endl;
  fileLog << "\tWithin-Chromosome GRM of " << n << " individuals has been saved in the file [" + grmOutBinFile_w + "] (in binary format)." << endl;  
  
  fstream Ab_Bin(grmOutBinFile_b.c_str(), ios::out | ios::binary);
  if (!Ab_Bin) throw (("Error: can not open the file [" + grmOutBinFile_b + "] to write.").c_str());
  f_buf = 0.0;
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++) {
      f_buf = G_between(i, j);
      Ab_Bin.write((char*) &f_buf, size);
    }
  }
  Ab_Bin.close();
  
  cout << "\tBetween-Chromosome GRM of " << n << " individuals has been saved in the file [" + grmOutBinFile_b + "] (in binary format)." << endl;
  fileLog << "\tBetween-Chromosome GRM of " << n << " individuals has been saved in the file [" + grmOutBinFile_b + "] (in binary format)." << endl;

  time_t t2 = time(0);   // get time now
  struct tm * now2 = localtime( & t2 );
  cout <<"\n\tAnalysis ends: ";
  cout << (now2->tm_year + 1900) << '-'
       << (now2->tm_mon + 1) << '-'
       <<  now2->tm_mday << " at "
       <<  now2->tm_hour <<":"<<now2->tm_min<<":"<<now2->tm_sec
       << ".\n";
  
  fileLog <<"\n\tAnalysis ends: ";
  fileLog << (now2->tm_year + 1900) << '-'
          << (now2->tm_mon + 1) << '-'
          <<  now2->tm_mday << " at "
          <<  now2->tm_hour <<":"<<now2->tm_min<<":"<<now2->tm_sec
          << ".\n";

  clock_t toc = clock();
  float time_elapsed = (float)(toc - tic) / CLOCKS_PER_SEC;
  float TimeElapsed = time_elapsed / nThreads;
  printf("\n\tTime elapsed: %f seconds.\n\n", time_elapsed / nThreads);
  fileLog<<"\n\tTime elapsed: "<< TimeElapsed <<" seconds.\n\n";
  

  fileLog.close();
  return EXIT_SUCCESS;
}

