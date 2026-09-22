#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include <admodel.h>
  #include <time.h>
  ofstream  ofParametersMCMC("./CheckInput/NSParasMCMC.dat");
  ofstream  ofInputData("./CheckInput/Global.dat");             //check control file input
  ofstream  ofBiologyData("./CheckInput/Biology_Data.dat");     //check biology data input
  ofstream  ofCAAData("./CheckInput/Catch_Data.dat");           //check catch data input
  ofstream  ofIAAData("./CheckInput/Index_Data.dat");           //check survey data input
  ofstream  ofGMData("./CheckInput/Growth_Matrix.dat");         //check growth matrix
  ofstream  ofOGDData("./CheckInput/BRP.dat");                  //check biological reference point data
  ofstream  ofPIRData("./CheckInput/Prior.dat");                //check prior input
  ofstream  ofPJData("./CheckInput/Projection.dat");            //check projection setting
  ofstream  ofINIData("./CheckInput/Initials.dat");             //check parameter initials input
  ofstream  ofRuntimelog("./CheckInput/Runtime.log",ios::app);  //runtime log
  ofstream  WARNING("./CheckInput/warning.log");                //warning log
  // Record the Beginning and Ending time for Program! Inherit From ASAP2 Code (NOAA TOOL BOX 2011) 
  time_t    tBeginTime;
  time_t    tEndTime;
  long      lHour;
  long      lMinute;
  long      lSecond;
  double    dElapsed_time;
  int       iBeginYear;
  int       iEndYear;
  int       itimestep;
  /////////////////////////////////////////////////////////////
  streampos tmpPirror;
  streampos tmpCurrent;
  ////define///////////////////////////////////////////////////////////
  #define   Debug_Status              1
  //#define   PI                        3.141592653897323846
  #define   PI2                       6.283185307794647692
  #define   MINPOSITIVE               0.00000000001
  #define   SIGMA2FILL                9.21044   //log(100*100 +1)
  #define   DIGITALPREC               0.000001
  // 
  #define   ICHECK(object)    ofInputData << "#" #object "\n " << object << endl;
  #define   ICHECKGD(object)  ofBiologyData << "#" #object "\n " << object << endl;
  #define   ICHECKCAA(object) ofCAAData << "#" #object "\n " << object << endl;
  #define   ICHECKIAA(object) ofIAAData << "#" #object "\n " << object << endl;
  #define   ICHECKGM(object)  ofGMData << "#" #object "\n " << object << endl;
  #define   ICHECKOGD(object) ofOGDData << "#" #object "\n " << object << endl;
  #define   ICHECKPIR(object) ofPIRData << "#" #object "\n " << object << endl;
  #define   ICHECKPJD(object) ofPJData << "#" #object "\n " << object << endl;
  #define   ICHECKINI(object) ofINIData << "#" #object "\n " << object << endl;
  //
  #define   RECORD(object)    ofRuntimelog <<  object << endl;
  //
  #define   CloseRecord()     ofRuntimelog.close(); 
  #define   CloseGM()         ofGMData.close();
  #define   CloseIAA()        ofIAAData.close();
  #define   CloseCAA()        ofCAAData.close();
  #define   CloseGD()         ofBiologyData.close();
  #define   CloseOGD()        ofOGDData.close();
  #define   ClosePIR()        ofPIRData.close();
  #define   CloseGBD()        ofInputData.close();
  #define   ClosePJD()        ofPJData.close();
  #define   CloseINI()        ofINIData.close();
  //
  #if  Debug_Status 
     ofstream  ofDebug("./CheckInput/Debug.log");  
     #define   ICHECKDEB(object) ofDebug << "#" #object "\n " << object << endl;
  #endif
  double CalLogFactorial(int iFactor)
   {
      int    i;
      double dR;
      dR=0.0;
      if(iFactor>=2)
      {
        for(i=2;i<=iFactor;i++)
          dR+=log(double(i));
      }
      return dR;
   }
   //-Ln(L)
  double CalLikelihoodConst(double dObs,int iModel)
  {
      double dR;
      //double dAdd=MINPOSITIVE;
      dR=0.0;
      switch(iModel)
      {   
        case 1://robust from lobster
           dR=0.5*log(PI2);              //log(sqrt(PI2)*dSigma);
           break;
        case 2://t-dis from lobster
           dR=-log(1.32934/sqrt(4.0*PI));
           break;
        case 3://for deviation parameters e.g. F; R; q.
           dR=0.5*log(PI2);             
           break; 
        case 4://log-normal
           if(dObs>0.0)
                dR=0.5*log(PI2)+log(dObs);
           else
                dR=0.0;
           break;
        case 5://log- normal without obs
           dR=0.5*log(PI2);
           break;
        case 6://normal 
           dR=0.5*log(PI2);
           break;
        case 7:
           dR=log(PI);
           break;
        default:
           dR=0.0;
           break;
      }
      return dR;
   }
    //-Ln(L)  
   dvariable CalLikelihoodVar(double obs,double dSigma2, const prevariable & pred, const int &iModel)
   {//   dvariable pred
       dvariable dvarR;
       dvariable  dvarTemp;
       //double dAdd=MINPOSITIVE;
       dvarR=0.0;
       switch(iModel)
       { 
          case 1:
             if(obs>0.0 && pred>0.0)
                  dvarR=-log(mfexp(-square(log(obs)-log(pred))/(2*dSigma2))+0.01)+log(dSigma2)*0.5;
             else
                  dvarR=0.0;
             break;
          case 2:
              // fattail likelihood
              if(obs>0.0 && pred>0.0)
                     dvarR= 2.5*log(1.0+square(log(obs)-log(pred))/(4.0*dSigma2));
              else
                     dvarR=0.0;
              break;
          case 3://normal //Deviation 
              dvarR= 0.5*square(pred)/dSigma2+log(dSigma2)*0.5;
              break;
          case 4://log-normal 
              if(obs>0.0 && pred>0.0)
                   dvarR=0.5*square(log(obs)-log(pred))/dSigma2+log(dSigma2)*0.5;
              else
                   dvarR=0.0;
              break;
          case 5://log-normal without obs 
              if(obs>0.0 && pred>0.0)
                    dvarR=0.5*square(log(obs)-log(pred))/dSigma2;
              else
                    dvarR=0.0;
              break;
          case 6://normal  
                    dvarR=0.5*square(obs-pred)/dSigma2+log(dSigma2)*0.5;
              break;
          case 7:
                dvarTemp=0.675*sqrt(dSigma2);
                dvarR=log(dvarTemp*(1.0+square((obs-pred)/dvarTemp)));  
                break;
          default:
              dvarR=0.0;
              break;
       }
      return dvarR;
   }
    //-Ln(L)  
   dvariable CalLikelihoodVar(double obs, const prevariable  &dSigma2, const prevariable & pred, const int &iModel)
   {//   dvariable pred
       dvariable dvarR;
       dvariable  dvarTemp;
       //double dAdd=MINPOSITIVE;
       dvarR=0.0;
       switch(iModel)
       { 
          case 1:
             if(obs>0.0 && pred>0.0)
                  dvarR=-log(mfexp(-square(log(obs)-log(pred))/(2*dSigma2))+0.01)+log(dSigma2)*0.5;
             else
                  dvarR=0.0;
             break;
          case 2:
              // fattail likelihood
              if(obs>0.0 && pred>0.0)
                     dvarR= 2.5*log(1.0+square(log(obs)-log(pred))/(4.0*dSigma2));
              else
                     dvarR=0.0;
              break;
          case 3://normal //Deviation 
              dvarR= 0.5*square(pred)/dSigma2+log(dSigma2)*0.5;
              break;
          case 4://log-normal 
              if(obs>0.0 && pred>0.0)
                   dvarR=0.5*square(log(obs)-log(pred))/dSigma2+log(dSigma2)*0.5;
              else
                   dvarR=0.0;
              break;
          case 5://log-normal without obs 
              if(obs>0.0 && pred>0.0)
                    dvarR=0.5*square(log(obs)-log(pred))/dSigma2+log(dSigma2)*0.5;
              else
                    dvarR=0.0;
              break;
          case 6://normal  
                    dvarR=0.5*square(obs-pred)/dSigma2+log(dSigma2)*0.5;
              break;
          case 7:
                dvarTemp=0.675*sqrt(dSigma2);
                dvarR=log(dvarTemp*(1.0+square((obs-pred)/dvarTemp)));  
                break;
          default:
              dvarR=0.0;
              break;
       }
      return dvarR;
   }
    dvariable CalStockRecruitment( const prevariable &dvSSB,  const prevariable &dvAlpha, const prevariable &dvBelta,const int &iModel)
    {
                 dvariable dvR;
                 dvR=0.0;
                 switch(iModel)
                 {
                     case 1:
                          dvR=dvAlpha;
                          break;
                     case 2:
                          dvR=dvAlpha*dvSSB/(dvBelta+dvSSB);
                          break;
                    case 3:
                          dvR=dvAlpha*dvSSB*mfexp(-dvBelta*dvSSB);
                          break;
                    default:
                         break;
                }
                return dvR;  
    }
     dvariable CalStockRecruitmentF( double vEnri1,double vEnri2,double vEnri3,const prevariable &dvSSB,const prevariable &dvAlpha,const prevariable &dvBelta,const dvar_vector &dvEnri,const int &iModel)
     {
                 dvariable dvR;
                 dvariable dvSum;
                 int i;
                 //int imax,imin;
                 dvSum=0.0;
                 dvSum +=dvEnri(1)*vEnri1+dvEnri(2)*vEnri2+dvEnri(3)*vEnri3;
                  switch(iModel)
                 {
                     case 1:
                          dvR=dvAlpha;
                          break;
                     case 2:
                          dvR=dvAlpha*dvSSB/(dvBelta+dvSSB);
                          break;
                     case 3:
                          dvR=dvAlpha*dvSSB*mfexp(-dvBelta*dvSSB);
                          break;
                     case 4:
                          dvR=mfexp(log(dvAlpha)+dvBelta*log(dvSSB)+dvSum);//Cushing E
                          break;
                     case 5:
                          dvR=mfexp(log(dvAlpha)+log(dvSSB)-log(dvBelta+dvSSB)+dvSum);//BH E
                          break;
                     case 6:
                          dvR=dvAlpha*dvSSB*mfexp(-dvBelta*dvSSB+dvSum);//Ricker E
                          break;
                     default:
                         break;
                 }
                 return dvR;  
     }
  //////////////////////////////////////////////////////////////////////////////////
  dvariable GrowthVariance(const prevariable &dvDeltaL,const prevariable &dvDeltaLVar, const double& dL)   
  {
      dvariable dvTemp;
      dvTemp=1.0/sqrt(PI2*dvDeltaLVar)*mfexp(-(dL-dvDeltaL)*(dL- dvDeltaL)*0.5/dvDeltaLVar);
      return dvTemp;
  }
  ///////////////////////////////////////////////////////////////////////////////////
  dvariable Integrate(const prevariable &dvDeltaL,const prevariable &dvDeltaLVar,double xMin,double xMax,int Panels,int type)
  {
    dvariable result = 0.0;
    int i;
    double interval = (xMax-xMin)/Panels;
    if(type == 0)
    {
        double xIterator = xMin+interval/2;
        for(int i=0;i<Panels;i++)
        {
            result += GrowthVariance(dvDeltaL,dvDeltaLVar,xIterator);
            xIterator += interval;
        }
        result*=interval;
        return result;
    }
    else if(type == 1)
    {
        double xIterator = xMin;
        for(int i=0;i<=Panels;i++)
        {
            result += GrowthVariance(dvDeltaL,dvDeltaLVar,xIterator);
            xIterator += interval;
        }
        result -=(GrowthVariance(dvDeltaL,dvDeltaLVar,xMin)+GrowthVariance(dvDeltaL,dvDeltaLVar,xMax))*0.5;//(0.5*(func(xMin)+func(xMax)));
        result *= interval;
        return result;
    }
    else if(type == 2)
    {
        if(Panels%2!=0)
        {
            cerr<<"Integer type Panels must be even to use Simpson's rule."<<endl;
            abort();
        }
        else
        {
            double xIterator = xMin;
            double xIt2;
            for( i=0;i<=Panels;i++)
            {
                if(i%2 == 0)
                {
                    result += 2*GrowthVariance(dvDeltaL,dvDeltaLVar,xIterator);
                }
                else
                {
                    result += 4*GrowthVariance(dvDeltaL,dvDeltaLVar,xIterator);
                }
                xIterator += interval;
                if(dvDeltaLVar<0.001)
                {
                    xIt2=xIterator+interval;
                    if(xIterator<dvDeltaL && xIt2>dvDeltaL)
                            xIterator=value(dvDeltaL);
                }
            }
            result -= (GrowthVariance(dvDeltaL,dvDeltaLVar,xMin)+GrowthVariance(dvDeltaL,dvDeltaLVar,xMax));//(func(xMin)+func(xMax));
            result *= interval/3;
        }
        return result;
    }
    return -1.0;
  }
  //////////////////////////////////////////////////////////////////
   dvariable Integrate(const prevariable &dvDeltaL,const prevariable &dvDeltaLVar,double a,double b)
  {
        dvariable T2n,I2n,Tn,In;
        dvariable sigma;
	double h;
	int n=1;
	const double eps=1e-6;
	h=b-a;
	T2n=I2n=(GrowthVariance(dvDeltaL,dvDeltaLVar,a)+GrowthVariance(dvDeltaL,dvDeltaLVar,b))*0.5*h;//h*(f(a)+f(b))/2;
	In=0;
	while(fabs(I2n-In)>=eps)
	{
		Tn=T2n;
		In=I2n;
		sigma=0.0;
		for(int k=0;k<n;k++)
		{
			double x=a+(k+0.5)*h;
			sigma+=GrowthVariance(dvDeltaL,dvDeltaLVar,x);//f(x);
		}
		T2n=(Tn+h*sigma)/2.0;
		I2n=(4*T2n-Tn)/3.0;
		n*=2;
		h/=2;
	}
	return I2n;
    }
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <NSLSAP01.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
 No_Report=0;
 Ncycle=3;
ad_comm::change_datafile_name("./InputFiles/Control.DAT");                   //read control data from file named GLOBAL.DAT
  Flagtimestep.allocate("Flagtimestep");
ICHECK(Flagtimestep);
  DiYearNum.allocate("DiYearNum");
ICHECK(DiYearNum);
  DiSeasonNum.allocate("DiSeasonNum");
ICHECK(DiSeasonNum);
  DvMonthInSeason.allocate(1,DiSeasonNum,"DvMonthInSeason");
ICHECK(DvMonthInSeason);
  DiBeginYear.allocate("DiBeginYear");
ICHECK(DiBeginYear); 
  DiCalBeginYear_Ini.allocate("DiCalBeginYear_Ini");
ICHECK(DiCalBeginYear_Ini)
  DiCalEndYear_Ini.allocate("DiCalEndYear_Ini");
ICHECK(DiCalEndYear_Ini)
  DiLikelyConstFlag.allocate("DiLikelyConstFlag");
ICHECK(DiLikelyConstFlag);
  DiCohortTrackingBeginYear.allocate("DiCohortTrackingBeginYear");
 ICHECKOGD(DiCohortTrackingBeginYear)
  DiTestValue.allocate("DiTestValue");
ICHECK(DiTestValue);
    if( DiTestValue !=999)
    {
        cout<< "Input Error in Control Data"<<endl;
        RECORD("Data Input Error, the TESTVALUE not Equal 999");
        CloseRecord();
        exit(-1);                                          // echo %errorlevel%  or BOOL GetExitCodeProcess(HANDLE hProcess,LPDWORD lpExitCode)
    }
    DiCalBeginYear= DiCalBeginYear_Ini- DiBeginYear+1;
    if(DiCalBeginYear<1)
    {
         WARNING<<"Error:  first year for assessment"<<endl;
         DiCalBeginYear=1;
    }
    DiCalEndYear=DiCalEndYear_Ini- DiBeginYear+1;
    if(DiCalEndYear>DiYearNum)
    {
       WARNING<<"Error: last year for assessment"<<endl;
       DiCalEndYear=DiYearNum;
    }
    if( DiCalEndYear<= DiCalBeginYear)
    {
        RECORD("Error  DiCalEndYear LESS Than  DiCalBeginYear");
        CloseRecord();
        exit(-12);
    }
   iBeginYear=DiCalBeginYear_Ini;
   iEndYear=DiCalEndYear_Ini;
   itimestep=Flagtimestep;
   #if Debug_Status
       ICHECK(DiCalBeginYear);
       ICHECK(DiCalEndYear);
   #endif
   RECORD("Control Data End ");
   #if Debug_Status
    if (Flagtimestep==1){cout<<"Time step: year"<<endl;}
    if (Flagtimestep==4){cout<<"Time step: season"<<endl;}
      cout<<"Data Begin Year:"<<DiBeginYear<<endl;
      cout<<"Caculation Begin Year:"<<DiCalBeginYear_Ini<<endl;
      cout<<"Caculation End Year:"<<DiCalEndYear_Ini<<endl;
      cout<<""<<endl;
   #endif
   CloseGBD(); 
cout<<"Control Data Input Completed"<<endl;
ad_comm::change_datafile_name("./InputFiles/Biology_Data.DAT");
  DiSizeBinNum.allocate("DiSizeBinNum");
ICHECKGD(DiSizeBinNum);
  DvSizeBinData.allocate(0,DiSizeBinNum,"DvSizeBinData");
ICHECKGD(DvSizeBinData);
 DiWeightAtSizeDataFlag=0;
   if(DiWeightAtSizeDataFlag)
   {
     y=DiYearNum;
     z=DiSizeBinNum+1;
   }
   else
   {
     y=0;
     z=0;
   }
  #if Debug_Status
     ICHECKGD(DiWeightAtSizeDataFlag);
  #endif
  DmWeightAtSize_Ini.allocate(1,y,1,z,"DmWeightAtSize_Ini");
 if(DiWeightAtSizeDataFlag)    ICHECKGD(DmWeightAtSize_Ini); 
   if(DiWeightAtSizeDataFlag)
   {
     y=0;
     NP=0;
   }
   else
   {
     y=DiYearNum;
     NP=2+1;
   }
  DmWeightAB.allocate(1,y,1,NP,"DmWeightAB");
if(!DiWeightAtSizeDataFlag)    ICHECKGD( DmWeightAB);          
 DiMaturityAtSizeDataFlag=1;
   if(DiMaturityAtSizeDataFlag)
   {
     y=DiYearNum;
     z=DiSizeBinNum+1;
   }
   else
   {
     y=0;
     z=0;
   }
  #if Debug_Status
     ICHECKGD(DiMaturityAtSizeDataFlag);
  #endif
  DmMaturityAtSize_Ini.allocate(1,y,1,z,"DmMaturityAtSize_Ini");
 if(DiMaturityAtSizeDataFlag)   ICHECKGD(DmMaturityAtSize_Ini);
   if(!DiMaturityAtSizeDataFlag)
   {
     y=DiYearNum;
     NP=4;
   }
   else
   {
     y=0;
     NP=0;
   }
  DmMaturityFuncPara.allocate(1,y,1,NP,"DmMaturityFuncPara");
if(!DiMaturityAtSizeDataFlag)   ICHECKGD(DmMaturityFuncPara);                                                                     
 DiNaturalMortalityFlag =0;                                      // 1:estimated M
   if(DiNaturalMortalityFlag)     
   {
        s= DiYearNum  * DiSeasonNum;   
        z= DiSizeBinNum+2; // Year Season Bin1 Bin2 Bin3 ...  
   }
   else
   {
        s=0;
        z=0;
   }
  #if Debug_Status
     ICHECKGD(DiNaturalMortalityFlag);
  #endif
  DmNaturalMAtSize_Ini.allocate(1,s,1,z,"DmNaturalMAtSize_Ini");
if( DiNaturalMortalityFlag ) ICHECKGD(DmNaturalMAtSize_Ini);
   if(!DiNaturalMortalityFlag)     
   {
       // lL1= DiSeasonNum;   
        z= DiSizeBinNum; // Year Season Bin1 Bin2 Bin3 ...  
   }
   else
   {
       //lL1=0;
       z=0;
   }
  DvNaturalMSizeWeight.allocate(1,z,"DvNaturalMSizeWeight");
 if(!DiNaturalMortalityFlag)  ICHECKGD(DvNaturalMSizeWeight);
   if(!DiNaturalMortalityFlag)     
   {
       // lL1= DiSeasonNum;   
        y= DiYearNum; // Year Season Bin1 Bin2 Bin3 ...  
   }
   else
   {
       //lL1=0;
       y=0;
   }
  DvNaturalMYearWeight.allocate(1,y,"DvNaturalMYearWeight");
 if(!DiNaturalMortalityFlag)  ICHECKGD( DvNaturalMYearWeight);
  DiRecruitPrjVectNP.allocate("DiRecruitPrjVectNP");
ICHECKGD(DiRecruitPrjVectNP);                           //  how many parameters(DiRecruitPrjVectNP) was estimated
  DvRecruitmentSeaRate.allocate(1,DiSeasonNum,"DvRecruitmentSeaRate");
ICHECKGD(DvRecruitmentSeaRate)                     //     R(s)=R*rate(s)
  DnSSBmonth.allocate("DnSSBmonth");
ICHECKGD(DnSSBmonth);
    if (Flagtimestep==1)
        DvRecruitmentSeaRate(1)=1;
    DnFracYearSSB=(DnSSBmonth-1)/12;
  
    temp=sum(DvRecruitmentSeaRate);
    if(temp>0)
       DvRecruitmentSeaRate/=temp;
    else
    {
        RECORD("Error: DvRecruitmentSeaRate less Than  0");
        CloseRecord();
        exit(-21);   
    }
    temp1=0.0;
    temp=DnFracYearSSB *sum(DvMonthInSeason);
    for(i=1;i<=DiSeasonNum;i++)
    { 
       temp1=temp1+DvMonthInSeason(i);
       if(temp<=temp1)
       {
           DiSBBSeason=i;
           DnFracSeasonSSB=(temp-(temp1-DvMonthInSeason(i)))/DvMonthInSeason(i);
           break;
       }
    }
   #if Debug_Status
      ICHECKGD(DnFracYearSSB);
      ICHECKGD(DiSBBSeason);
      ICHECKGD(DnFracSeasonSSB);
      ICHECKGD(temp);
      ICHECKGD(temp1);
   #endif
  DiRSFlag.allocate("DiRSFlag");
ICHECKGD(DiRSFlag);  
  DiFlagEvntoRdevs.allocate("DiFlagEvntoRdevs");
  DiEnvNum.allocate("DiEnvNum");
  DmEnv.allocate(1,DiYearNum,1,DiEnvNum,"DmEnv");
 DiFecundityFlag=1;                                        // ??????????? egg  
DiYearBeforeEndForRDev=0;
 DiRecruitLogDevsSDEstFlag =1;
  DiNYear1ChoiceFlag.allocate("DiNYear1ChoiceFlag");
ICHECKGD(DiNYear1ChoiceFlag);       //   0: estimate the abundance of each size bin
lL1=DiSizeBinNum;
  DvNYear1Pia.allocate(1,lL1,"DvNYear1Pia");
  DiTestValueGD.allocate("DiTestValueGD");
ICHECKGD(DiTestValueGD);
    if( DiTestValueGD!=999)
    {
        cout<< "Input Error in Biology Data Part"<<endl;
        RECORD("Data Input Error, the TESTVALUE not Equal 999");
        CloseRecord();
        exit(-22);
    }
    #if Debug_Status
        ICHECKGD(DiFecundityFlag);
        ICHECKGD(DiRecruitLogDevsSDEstFlag);
    #endif
  DmWeightAtSize.allocate(1,DiYearNum,1,DiSizeBinNum);
  DmMaturityAtSize.allocate(1,DiYearNum,1,DiSizeBinNum);
  DmFecundity.allocate(1,DiYearNum,1,DiSizeBinNum);
     
     for(i=1;i<=DiYearNum;i++)
     {
        if(!DiMaturityAtSizeDataFlag)
        {
           for(j=1;j<=DiSizeBinNum;j++)
           {
              temp= DmMaturityFuncPara(i,2)/(1.0+mfexp(- DmMaturityFuncPara(i,3)*( 0.5*(DvSizeBinData(j-1)+ DvSizeBinData(j))- DmMaturityFuncPara(i,4))));
              DmMaturityAtSize(i,j)=temp;  
           }
        }
        else
        {
           for(j=1;j<=DiSizeBinNum;j++)
              DmMaturityAtSize(i,j)=DmMaturityAtSize_Ini(i,j+1);
        }
     }
    
     ///////////////////////////////////////////
    
     for(i=1;i<=DiYearNum;i++)
     {
        if(!DiWeightAtSizeDataFlag)
        {
           for(j=1;j<=DiSizeBinNum;j++)
           {//sgwj debug
               // DmWeightAtSize(i,j)=0.5*pow(DvSizeBinData(j-1),DmWeightAB(j,3))*DmWeightAB(i,2);
                //DmWeightAtSize(i,j)+=0.5*pow(DvSizeBinData(j),DmWeightAB(j,3))*DmWeightAB(i,2);
                DmWeightAtSize(i,j)=pow((DvSizeBinData(j-1)+DvSizeBinData(j))*0.5,DmWeightAB(i,3))*exp(DmWeightAB(i,2));
                DmFecundity(i,j)=DmMaturityAtSize(i,j)*DmWeightAtSize(i,j)*0.5;
           } 
        }
        
     }
     
  //!!if(DiNYear1ChoiceFlag==1||DiNYear1ChoiceFlag==2) ICHECKOGD(DvNYear1Pia);
  
     #if Debug_Status  
          ICHECKGD(DmFecundity);
          ICHECKGD(DmMaturityAtSize);
          ICHECKGD(DmWeightAtSize);
          cout<<"Biology Data Input Completed"<<endl;
     #endif
     CloseGD();
    if (Flagtimestep==1)
    { ad_comm::change_datafile_name("./InputFiles/Year/CatchDataYear.DAT");}
    else 
    {ad_comm::change_datafile_name("./InputFiles/Season/CatchDataSeason.DAT");}
  DiFleetNum.allocate("DiFleetNum");
ICHECKCAA(DiFleetNum);
  DimCatchBiomassFlag.allocate(1,DiFleetNum,1,Flagtimestep,"DimCatchBiomassFlag");
ICHECKCAA(DimCatchBiomassFlag);
  DimCSelStartSizeBin.allocate(1,DiFleetNum,1,Flagtimestep,"DimCSelStartSizeBin");
ICHECKCAA(DimCSelStartSizeBin);
  DimCSelEndSizeBin.allocate(1,DiFleetNum,1,Flagtimestep,"DimCSelEndSizeBin");
ICHECKCAA(DimCSelEndSizeBin);
  DimCatchCompLikelihoodFlag.allocate(1,DiFleetNum,1,Flagtimestep,"DimCatchCompLikelihoodFlag");
ICHECKCAA(DimCatchCompLikelihoodFlag);
  DimCatchTotalLikelihoodFlag.allocate(1,DiFleetNum,1,Flagtimestep,"DimCatchTotalLikelihoodFlag");
ICHECKCAA(DimCatchTotalLikelihoodFlag);
  DimCPUELikelihoodFlag.allocate(1,DiFleetNum,1,Flagtimestep,"DimCPUELikelihoodFlag");
ICHECKCAA(DimCPUELikelihoodFlag);
  DmCatchCompLambda.allocate(1,DiFleetNum,1,Flagtimestep,"DmCatchCompLambda");
ICHECKCAA(DmCatchCompLambda);
  DmCatchTotalLambda.allocate(1,DiFleetNum,1,Flagtimestep,"DmCatchTotalLambda");
ICHECKCAA(DmCatchTotalLambda);
  DmCPUELambda.allocate(1,DiFleetNum,1,Flagtimestep,"DmCPUELambda");
ICHECKCAA(DmCPUELambda);        
  DiCAALength.allocate("DiCAALength");
  DmCatchAtSize_Ini.allocate(1,DiCAALength,1,3+6+DiSizeBinNum,"DmCatchAtSize_Ini");
ICHECKCAA(DmCatchAtSize_Ini);                                            //adapt to the input data format, it will be replaced by other variable/const in calculation           
  Dd4CatchAtSizeO.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum);
  Dd3CatchTotalO.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear);
  Dd3CPUEO.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear);
  Dd3CatchESSInput.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear);
  Dd3CatchTotalSigma.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear);
  Dd3CatchTotalSigma2.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear);
  Dd3CPUESigma.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear);
  Dd3CPUESigma2.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear);
  Dd4CatchPropAtSizeO.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum);
  DimFleetFirstYear.allocate(1,DiFleetNum,1,Flagtimestep);
lL1=(DiCalEndYear-DiCalBeginYear+1)*Flagtimestep;                          
  DimFleetFZeroIndex.allocate(1,DiFleetNum,1,lL1);
   // Format changed 
   Dd4CatchAtSizeO.initialize();
   Dd3CatchTotalO.initialize();
   Dd3CPUEO.initialize();
   Dd3CatchESSInput.initialize();
   Dd3CatchTotalSigma.initialize();
   Dd3CatchTotalSigma2.initialize();
   Dd3CPUESigma.initialize();
   Dd3CPUESigma2.initialize();
   Dd4CatchPropAtSizeO.initialize();
   DimFleetFirstYear.initialize();
   //DivFleetFirstSeason.initialize();
   DimFleetFZeroIndex.initialize();
   for(i=1;i<=DiCAALength;i++)
   {
       kk=int( DmCatchAtSize_Ini(i,1)+0.01)-DiBeginYear+1; //year
       if(kk<DiCalBeginYear || kk> DiCalEndYear )
          continue;
       ii= int( DmCatchAtSize_Ini(i,2)+0.01); //Season
       jj= int(DmCatchAtSize_Ini(i,3)+0.01); //Fleet
       if(ii<=0||jj<=0||ii>DiSeasonNum||jj>DiFleetNum)
         continue;
       for(j=1;j<=DiSizeBinNum;j++)
       {
           Dd4CatchAtSizeO(jj,ii,kk,j)= DmCatchAtSize_Ini(i,9+j);  
       }
       Dd3CatchESSInput(jj,ii,kk)=DmCatchAtSize_Ini(i,9);
       //
       Dd3CatchTotalO(jj,ii,kk)=DmCatchAtSize_Ini(i,4);
      /////////////////////////////////////////////////////////////////
       if(Dd3CatchTotalO(jj,ii,kk)> MINPOSITIVE)
       {
              DimFleetFZeroIndex(jj,(kk-DiCalBeginYear)*Flagtimestep+ii)=1;
              if(DimFleetFirstYear(jj,ii)==0)
              {
                    DimFleetFirstYear(jj,ii)=kk;
              }      
       }
      //////////////////////////////////////////////////////////////////////
       if(DmCatchAtSize_Ini(i,5)>MINPOSITIVE)
             Dd3CatchTotalSigma2(jj,ii,kk)=DmCatchAtSize_Ini(i,5)* DmCatchAtSize_Ini(i,5);//log(DmCatchAtSize_Ini(i,5)* DmCatchAtSize_Ini(i,5)+1.0);
       else
             Dd3CatchTotalSigma2(jj,ii,kk)=SIGMA2FILL;
       Dd3CatchTotalSigma(jj,ii,kk)= sqrt(Dd3CatchTotalSigma2(jj,ii,kk));
     
       Dd3CPUEO(jj,ii,kk)=DmCatchAtSize_Ini(i,6);
       //sgwj debug  How to deal with Flag   CPUE=Catch/Effort ? or Effort(Obs)-Effort(Pred)         
       if(DmCatchAtSize_Ini(i,8)>MINPOSITIVE)
             Dd3CPUESigma2(jj,ii,kk)=DmCatchAtSize_Ini(i,8)* DmCatchAtSize_Ini(i,8);//log(DmCatchAtSize_Ini(i,7)* DmCatchAtSize_Ini(i,7)+1.0);
       else
             Dd3CPUESigma2(jj,ii,kk)=SIGMA2FILL;
       Dd3CPUESigma(jj,ii,kk)=sqrt(Dd3CPUESigma2(jj,ii,kk));
   }
   //compute the proportion of the catch at each selectivity size 
   for(i=1;i<=DiFleetNum;i++)
   {
       for(ii=1;ii<=Flagtimestep;ii++)
       {
            for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
            {
                nTemp=0.0;
                if( Dd3CatchTotalO(i,ii,j)> MINPOSITIVE)
                {
                     for(k=DimCSelStartSizeBin(i,ii);k<=DimCSelEndSizeBin(i,ii);k++)
                              nTemp+=Dd4CatchAtSizeO(i,ii,j,k);                     //sum(d3CatchAtSizeO(i,j,ivCSelStartSizeBin(i),ivCSelEndSizeBin(i)));
                     if(nTemp> MINPOSITIVE)
                         Dd4CatchPropAtSizeO(i,ii,j)(DimCSelStartSizeBin(i,ii),DimCSelEndSizeBin(i,ii))=Dd4CatchAtSizeO(i,ii,j)(DimCSelStartSizeBin(i,ii),DimCSelEndSizeBin(i,ii))/nTemp;
                }
            }
       }
   } 
  
   #if Debug_Status
     ICHECKCAA(Dd4CatchAtSizeO);
     ICHECKCAA(Dd4CatchPropAtSizeO);
     ICHECKCAA(DimFleetFZeroIndex);
     ICHECKCAA(Dd3CatchTotalO);
     ICHECKCAA(Dd3CatchTotalSigma2);
     ICHECKCAA(DimFleetFirstYear);
     ICHECKCAA(Dd3CatchESSInput);
     ICHECKCAA(Dd3CPUEO);
     ICHECKCAA(Dd3CPUESigma2);
    //cout<<"Catch Data Input Completed"<<endl;
    // exit(-99);
   #endif
   RECORD("Catch Data Input Completed");
  
  DimAvailCPUEFlag.allocate(1,DiFleetNum,1,Flagtimestep,"DimAvailCPUEFlag");
ICHECKCAA(DimAvailCPUEFlag);
DiAvailCPUENum=sum(DimAvailCPUEFlag);
  DivAvailCPUEIndexOld.allocate(1,DiAvailCPUENum);
lL1=DiFleetNum*Flagtimestep;
  DivOldIndexAvailCPUE.allocate(1,lL1);
    if(DiAvailCPUENum>0)
       DivAvailCPUEIndexOld.initialize();
    DivOldIndexAvailCPUE.initialize();
    i=0;
    for(ii=1;ii<=DiFleetNum;ii++)
    {
        for(j=1;j<=Flagtimestep;j++)
        {
           if(DimAvailCPUEFlag(ii,j))
           {
                 i++;
                 iL1=(ii-1)*Flagtimestep+j;
                 DivAvailCPUEIndexOld(i)=iL1;
                 DivOldIndexAvailCPUE(iL1)=i;
           }
       }
    }
   #if Debug_Status
     ICHECKCAA(DiAvailCPUENum);
     if(DiAvailCPUENum>0)
         ICHECKCAA( DivAvailCPUEIndexOld);
     ICHECKCAA(DivOldIndexAvailCPUE);
   #endif
  DiCPUEBlockNum.allocate("DiCPUEBlockNum");
ICHECKCAA(DiCPUEBlockNum);                              //  e.g. season 1, fleet 1, q(1,1) was same from 1967--2009,so breakpoint is to break 1967--1985---2009
  DivCPUEQChoiceFlag_Ini.allocate(1,DiCPUEBlockNum,"DivCPUEQChoiceFlag_Ini");
ICHECKCAA(DivCPUEQChoiceFlag_Ini);
    if(DiCPUEBlockNum>0)
    {
       iL1=Flagtimestep;
       iL2=DiYearNum;
       lL1=DiFleetNum+1;
    }
    else
    {   //
        iL1=-1;
        iL2=-1;
        lL1=-1; 
       //iL1=DiSeasonNum;
       //iL2=DiYearNum;
       //lL1=DiFleetNum+1;
    }
  Dd3CPUEBlockFlag.allocate(1,iL1,1,iL2,1,lL1,"Dd3CPUEBlockFlag");
ICHECKCAA(Dd3CPUEBlockFlag);
  DivCPUEBlockAvailFlag.allocate(1,DiCPUEBlockNum);
   if(DiCPUEBlockNum>0)
      DivCPUEBlockAvailFlag=0;
   // DiFleetSelAvailBlockNum=0;
   for(i=1;i<=Flagtimestep;i++)
   {
        for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
        {
            for(k=1;k<=DiFleetNum;k++)
            {
               if(Dd3CPUEO(k,i,j)>MINPOSITIVE && DimAvailCPUEFlag(k,i)==1)
               {
                    iL1=int(Dd3CPUEBlockFlag(i,j,k+1)+0.01);
                    if(iL1<1||iL1>DiCPUEBlockNum)
                            continue;
                    if(DivCPUEBlockAvailFlag(iL1)==0)
                           DivCPUEBlockAvailFlag(iL1)=1;
               }
            }
        }
   }
ICHECKCAA(DivCPUEBlockAvailFlag);
  DiFleetSelBlockNum.allocate("DiFleetSelBlockNum");
ICHECKCAA(DiFleetSelBlockNum); 
  DivFleetSelBlockFlag.allocate(1,DiFleetSelBlockNum,"DivFleetSelBlockFlag");
ICHECKCAA(DivFleetSelBlockFlag);
  DivFleetSelParaIndex_Ini.allocate(1,DiFleetSelBlockNum);
     
     DiFleetSelParaNumIni=0;
     for(i=1;i<=DiFleetSelBlockNum;i++)
     {
            DivFleetSelParaIndex_Ini(i)=   DiFleetSelParaNumIni;
            if(DivFleetSelBlockFlag(i)==1) DiFleetSelParaNumIni+=DiSizeBinNum;
            if(DivFleetSelBlockFlag(i)==2) DiFleetSelParaNumIni+=2;
            if(DivFleetSelBlockFlag(i)==3) DiFleetSelParaNumIni+=4;
     }
    
    //ICHECKCAA(DivFleetSelBlockFlag);
  Dd3FleetSelBlockFlag.allocate(1,Flagtimestep,1,DiYearNum,1,DiFleetNum+1,"Dd3FleetSelBlockFlag");
ICHECKCAA(Dd3FleetSelBlockFlag);// record  1:DiFleetSelBlockNum     
  DivFleetSelBlockAvailFlag.allocate(1,DiFleetSelBlockNum);
   
   DivFleetSelBlockAvailFlag=0;
  // DiFleetSelAvailBlockNum=0;
   for(i=1;i<=Flagtimestep;i++)
   {
        for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
        {
            for(k=1;k<=DiFleetNum;k++)
            {
               if(Dd3CatchTotalO(k,i,j)>MINPOSITIVE)
               {
                    iL1=int(Dd3FleetSelBlockFlag(i,j,k+1)+0.01);
                    if(iL1<1||iL1>DiFleetSelBlockNum)
                            continue;
                    if(DivFleetSelBlockAvailFlag(iL1)==0)
                           DivFleetSelBlockAvailFlag(iL1)=1;
               }
            }
        }
   }
  //   DiFleetSelAvailBlockNum=sum(DivFleetSelBlockAvailFlag);
  DiTestValueCAA.allocate("DiTestValueCAA");
ICHECKCAA(DiTestValueCAA);
  #if Debug_Status
     ICHECKCAA(DivCPUEBlockAvailFlag);
     ICHECKCAA(DiFleetSelParaNumIni);
     ICHECKCAA(DivFleetSelParaIndex_Ini);
     ICHECKCAA(DivFleetSelBlockAvailFlag);
     //ICHECKCAA(DiTestValueCAA);
     //cout<<"Fishery Selecitivty Data Input Completed!"<<endl;
     //cout<<"Fishery data end"<<endl;
  #endif
   if( DiTestValueCAA!=-22122)
   {
        cout<< "Input Error in Catch Data"<<endl;
        RECORD("Data Input Error, the TESTVALUE not Equal -22122");
        CloseRecord();
        exit(-4);
    }
    RECORD("Catch Data Input End");
    CloseCAA();
    cout<<"Catch Data Input Completed"<<endl;
ad_comm::change_datafile_name("./InputFiles/Survey_Data.DAT");
  DiIndexNum.allocate("DiIndexNum");
ICHECKIAA(DiIndexNum);
  DivIndexUnitsFlag_Ini.allocate(1,DiIndexNum,"DivIndexUnitsFlag_Ini");
 ICHECKIAA(DivIndexUnitsFlag_Ini);
  DivIndexSelStartSizeBin_Ini.allocate(1,DiIndexNum,"DivIndexSelStartSizeBin_Ini");
 ICHECKIAA(DivIndexSelStartSizeBin_Ini);
  DivIndexSelEndSizeBin_Ini.allocate(1,DiIndexNum,"DivIndexSelEndSizeBin_Ini");
 ICHECKIAA(DivIndexSelEndSizeBin_Ini);
  DivAvailIndexFlag.allocate(1,DiIndexNum,"DivAvailIndexFlag");
 ICHECKIAA(DivAvailIndexFlag);
  DivIndexCompLikelihoodFlag_Ini.allocate(1,DiIndexNum,"DivIndexCompLikelihoodFlag_Ini");
ICHECKIAA(DivIndexCompLikelihoodFlag_Ini);
  DivIndexTotalLikelihoodFlag_Ini.allocate(1,DiIndexNum,"DivIndexTotalLikelihoodFlag_Ini");
 ICHECKIAA( DivIndexTotalLikelihoodFlag_Ini);
  DvIndexCompLambda_Ini.allocate(1,DiIndexNum,"DvIndexCompLambda_Ini");
ICHECKIAA(DvIndexCompLambda_Ini);//*DvIndexTotalLambda(i)
  DvIndexTotalLambda_Ini.allocate(1,DiIndexNum,"DvIndexTotalLambda_Ini");
ICHECKIAA(DvIndexTotalLambda_Ini);
  DiSurveyIndexLength.allocate("DiSurveyIndexLength");
ICHECKIAA(DiSurveyIndexLength);
  DmIndexAtSizeData_Ini.allocate(1,DiSurveyIndexLength,1,6+DiSizeBinNum,"DmIndexAtSizeData_Ini");
 ICHECKIAA(DmIndexAtSizeData_Ini);
  DiSexAtSizeLamda.allocate("DiSexAtSizeLamda");
  DiIndexFlagSexComp.allocate("DiIndexFlagSexComp");
  DmSexAtSizeData_Ini.allocate(1,DiYearNum,1,DiSizeBinNum,"DmSexAtSizeData_Ini");
  DiIndexQBlockNum.allocate("DiIndexQBlockNum");
ICHECKIAA(DiIndexQBlockNum);
  DivIndexQChoiceFlag.allocate(1,DiIndexQBlockNum,"DivIndexQChoiceFlag");
ICHECKIAA(DivIndexQChoiceFlag);
  DimIndexQBlockFlag.allocate(1,DiYearNum,1,DiIndexNum+1,"DimIndexQBlockFlag");
ICHECKIAA(DimIndexQBlockFlag)
 DiAvailIndexNum=sum(DivAvailIndexFlag);  
  DivAvailIndexQFlag.allocate(1,DiIndexQBlockNum);
  DivIndexObsNum.allocate(1,DiAvailIndexNum);
   //  if(DiIndexQBlockNum>0)//??????????????????
      DivAvailIndexQFlag.initialize();
   //  if(DiAvailIndexNum>0)//?????????????????????
      DivIndexObsNum.initialize();
   for(i=1;i<=DiSurveyIndexLength;i++)
   {
       ii=int(DmIndexAtSizeData_Ini(i,1)+0.01);//year
       j=ii-DiBeginYear+1;
       if(j< DiCalBeginYear||j>DiCalEndYear)
               continue;
       jj=int(DmIndexAtSizeData_Ini(i,2)+0.01);//indexnum
       if(jj<1 ||jj>DiIndexNum)
                continue;                  //check Index number
       if(DivAvailIndexFlag(jj) && DmIndexAtSizeData_Ini(i,4)>MINPOSITIVE ) // zero or negative value for index means not included
       {
          lL2=sum(DivAvailIndexFlag(1,jj));
          DivIndexObsNum(lL2)++;
          kk=DimIndexQBlockFlag(j,jj+1);
          if(kk>=1 && kk<=DiIndexQBlockNum )
                 DivAvailIndexQFlag(kk)=1;
       }
   }
   //DiAvailIndexQNum=sum(DivAvailIndexQFlag);
  #if Debug_Status
     ICHECKIAA(DivAvailIndexQFlag);
     ICHECKIAA(DivIndexObsNum);
     ICHECKIAA(DiAvailIndexNum);
  #endif
  DivAvailIndexMirrorOld.allocate(1,DiAvailIndexNum);
    iL1=0;
    for(i=1;i<=DiIndexNum;i++)
    {
         if(DivAvailIndexFlag(i))
         {
               iL1++;
               DivAvailIndexMirrorOld(iL1)=i;
         }     
    }
   #if Debug_Status
      ICHECKIAA(DivAvailIndexMirrorOld);
   #endif
  DimIndexIndex.allocate(1,DiAvailIndexNum,1,DivIndexObsNum);
  DmIndexTotalO.allocate(1,DiAvailIndexNum,1,DivIndexObsNum);
  DmIndexTotalSigma2.allocate(1,DiAvailIndexNum,1,DivIndexObsNum);
  DmIndexTotalSigma.allocate(1,DiAvailIndexNum,1,DivIndexObsNum);
  DivIndexSelStartSizeBin.allocate(1,DiAvailIndexNum);
  DivIndexSelEndSizeBin.allocate(1,DiAvailIndexNum);
  DimIndexESSInput.allocate(1,DiAvailIndexNum,1,DivIndexObsNum);
  DvIndexTotalLambda.allocate(1,DiAvailIndexNum);
  DvIndexCompLambda.allocate(1,DiAvailIndexNum);
  DivIndexCompLikelihoodFlag.allocate(1,DiAvailIndexNum);
  DivIndexTotalLikelihoodFlag.allocate(1,DiAvailIndexNum);
  DivIndexUnitsFlag.allocate(1,DiAvailIndexNum);
  Dd3IndexMonthFlag.allocate(1,DiAvailIndexNum,1,DivIndexObsNum,1,2);
  Dd3IndexAtSizeO.allocate(1,DiAvailIndexNum,1,DivIndexObsNum,1,DiSizeBinNum);
  Dd3IndexPropAtSizeO.allocate(1,DiAvailIndexNum,1,DivIndexObsNum,1,DiSizeBinNum);
  
   //Dd3IndexAtSizeO.initialize();
   //Dd3IndexPropAtSizeO.initialize();
   //DivAvailIndexQChoiceFlag=0;
   lL1=0;//Available Index
   for(i=1;i<=DiIndexNum;i++)
   {
       if(DivAvailIndexFlag(i))
       {
            lL1++;
           
            DivIndexSelStartSizeBin(lL1)=DivIndexSelStartSizeBin_Ini(i);
            DivIndexSelEndSizeBin(lL1)=DivIndexSelEndSizeBin_Ini(i);
            DivIndexUnitsFlag(lL1)=DivIndexUnitsFlag_Ini(i);
            //DivIndexMonthFlag(lL1)= DivIndexMonthFlag_Ini(i);//sgwj 2012-04-03
            DvIndexTotalLambda(lL1)=DvIndexTotalLambda_Ini(i);
            DvIndexCompLambda(lL1)=DvIndexCompLambda_Ini(i);
            /////////////////////////////////////////////
            DivIndexCompLikelihoodFlag(lL1)=DivIndexCompLikelihoodFlag_Ini(i);
            DivIndexTotalLikelihoodFlag(lL1)=DivIndexTotalLikelihoodFlag_Ini(i);
        //    DivIndexQChoiceFlag(lL1)=DivIndexQChoiceFlag_Ini(i);
         }
   }
   //k=0; 
   DivIndexObsNum.initialize();
   Dd3IndexMonthFlag.initialize();
   Dd3IndexAtSizeO.initialize();
   Dd3IndexPropAtSizeO.initialize();
   for(i=1;i<=DiSurveyIndexLength;i++)
   {
         ii=int(DmIndexAtSizeData_Ini(i,1)+0.01);//year
         j=ii-DiBeginYear+1;
         if(j< DiCalBeginYear||j>DiCalEndYear)
               continue;
        //////////////////////////////////////////////////
         jj=int(DmIndexAtSizeData_Ini(i,2)+0.01);//indexnum
         if(jj<1 ||jj>DiIndexNum)
                continue;                  //check Index number
         if(DivAvailIndexFlag(jj) && DmIndexAtSizeData_Ini(i,4)>MINPOSITIVE) // zero or negative value for index means not included
         {
              lL2=sum(DivAvailIndexFlag(1,jj));
              DivIndexObsNum(lL2)++;
              k= DivIndexObsNum(lL2);
              DimIndexIndex(lL2,k)=j;
              ///////////////////////////////////////////////////////////////
              if(DmIndexAtSizeData_Ini(i,3)<0)
              {//sgwj 2012-04-14
                Dd3IndexMonthFlag(lL2,k,1)=-2;
                Dd3IndexMonthFlag(lL2,k,2)=-2;    
              }
              else
              {//sgwj 2012-04-14
                   nTemp=0.0;
                   for(iTemp=1;iTemp<=DiSeasonNum;iTemp++)
                   {
                         nTemp=DvMonthInSeason(iTemp)+nTemp;
                         if(DmIndexAtSizeData_Ini(i,3)<nTemp)
                         {
                               Dd3IndexMonthFlag(lL2,k,1)=iTemp;
                               Dd3IndexMonthFlag(lL2,k,2)= (nTemp-DmIndexAtSizeData_Ini(i,3)-DvMonthInSeason(iTemp))/DvMonthInSeason(iTemp);
                               break;
                         }
                   }     
              }
              ///////////////////////////////////////////////////////////////
              DmIndexTotalO(lL2,k)= DmIndexAtSizeData_Ini(i,4);
              if(DmIndexAtSizeData_Ini(i,5)>MINPOSITIVE)
                    DmIndexTotalSigma2(lL2,k)=DmIndexAtSizeData_Ini(i,5)* DmIndexAtSizeData_Ini(i,5);//log(DmIndexAtSizeData_Ini(i,4)* DmIndexAtSizeData_Ini(i,4)+1.0);
              else
                    DmIndexTotalSigma2(lL2,k)=SIGMA2FILL;
                        
             DmIndexTotalSigma(lL2,k)=sqrt( DmIndexTotalSigma2(lL2,k));
             DimIndexESSInput(lL2,k)=int(DmIndexAtSizeData_Ini(i,6)+0.001);
             nTemp=0.0;
             for(ii=1;ii<=DiSizeBinNum;ii++)
             {
               Dd3IndexAtSizeO(lL2,k,ii)= DmIndexAtSizeData_Ini(i,6+ii);
             }
             nTemp=sum(Dd3IndexAtSizeO(lL2,k)(DivIndexSelStartSizeBin(lL2),DivIndexSelEndSizeBin(lL2)));
             if(nTemp>0.0)
             {
                    for(iTemp=DivIndexSelStartSizeBin(lL2);iTemp<=DivIndexSelEndSizeBin(lL2);iTemp++)
                         Dd3IndexPropAtSizeO(lL2,k,iTemp)=Dd3IndexAtSizeO(lL2,k,iTemp)/nTemp;  //??????????????  k --->j
             }
         }
   }
  #if Debug_Status
        ICHECKIAA( DivIndexSelStartSizeBin);
        ICHECKIAA( DivIndexSelEndSizeBin);
       // ICHECKIAA( DivIndexMonthFlag);
        ICHECKIAA( Dd3IndexMonthFlag);
        ICHECKIAA( DivIndexUnitsFlag);
        ICHECKIAA( DvIndexTotalLambda);
        ICHECKIAA( DvIndexCompLambda);
        //////////////////
        ICHECKIAA( Dd3IndexAtSizeO);
        ICHECKIAA( DimIndexIndex);
        ICHECKIAA( DmIndexTotalO)
        ICHECKIAA( DmIndexTotalSigma);
        ICHECKIAA( DimIndexESSInput);
        ICHECKIAA( Dd3IndexPropAtSizeO);
     // exit(-99);
  #endif
  DivIndexFleetSelChoiceFlag_Ini.allocate(1,DiIndexNum,"DivIndexFleetSelChoiceFlag_Ini");
 ICHECK(DivIndexFleetSelChoiceFlag_Ini)
  DiIndexSelBlockNum.allocate("DiIndexSelBlockNum");
 ICHECKIAA(DiIndexSelBlockNum);//!!DiIndexSelBlockNum=DiIndexNum;              //                 ->   ICHECKIAA(DiIndexSelBlockNum)
  DivIndexSelFlag_Ini.allocate(1,DiIndexSelBlockNum,"DivIndexSelFlag_Ini");
 ICHECKIAA(DivIndexSelFlag_Ini);
  DivIndexSelParaIndex_Ini.allocate(1,DiIndexSelBlockNum);
  DimIndexSelBlockFlag.allocate(1,DiYearNum,1,DiIndexNum+1,"DimIndexSelBlockFlag");
ICHECKIAA(DimIndexSelBlockFlag);                                                                    //->  init_imatrix     DimIndexSelBlockFlag(1,DiYearNum,1,DiIndexNum+1)
  DivIndexSelBlockAvailFlag.allocate(1,DiIndexSelBlockNum);
  /*for(i=1;i<=DiYearNum;i++)                                    //sgwj 2012-04-14  delete 
  {       
      DimIndexSelBlockFlag(i,1)=i+DiBeginYear;                //sgwj 2012-04-14  delete 
      for(j=2;j<=DiIndexNum+1;j++)                           //sgwj 2012-04-14  delete 
         DimIndexSelBlockFlag(i,j)=j-1;                     //sgwj 2012-04-14  delete 
  }*/                                                        //sgwj 2012-04-14  delete 
  DivIndexSelBlockAvailFlag=0;
  
  for(j=1;j<=DiAvailIndexNum;j++)
  {
    for(i=1;i<=DivIndexObsNum(j);i++)
    {
          k=DimIndexIndex(j,i);
          ii=DivAvailIndexMirrorOld(j);
          if(DivIndexFleetSelChoiceFlag_Ini(ii)<0)
          {//not used the fleet selectivity
             jj=DimIndexSelBlockFlag(k,ii+1);
             if(jj>=1&&jj<=DiIndexSelBlockNum)
                 DivIndexSelBlockAvailFlag(jj)=1;
          }
     }
  }
  DiIndexSelParaNumIni=0;
  for (i=1;i<=DiIndexSelBlockNum;i++)
  {
       DivIndexSelParaIndex_Ini(i)=DiIndexSelParaNumIni;
       if( DivIndexSelFlag_Ini(i)==1) DiIndexSelParaNumIni+=DiSizeBinNum;
       if( DivIndexSelFlag_Ini(i)==2) DiIndexSelParaNumIni+=2;
       if( DivIndexSelFlag_Ini(i)==3) DiIndexSelParaNumIni+=4;
  }
  #if Debug_Status
        //ICHECKIAA(DivIndexFleetSelChoiceFlag_Ini)
        ICHECKIAA(DiIndexSelParaNumIni);
        ICHECKIAA(DivIndexSelParaIndex_Ini);
        ICHECKIAA(DivIndexSelBlockAvailFlag);
      
  #endif
  DivAvailIndexSelIndexOld.allocate(1,DiAvailIndexNum);
  DivIndexFleetSelChoiceFlag.allocate(1,DiAvailIndexNum);
   
  lL1=0;
  for (i=1;i<=DiIndexNum;i++)
  {
      if(DivAvailIndexFlag(i))
      {
         lL1++;
         DivIndexFleetSelChoiceFlag(lL1)=DivIndexFleetSelChoiceFlag_Ini(i);
         DivAvailIndexSelIndexOld(lL1)=i;    //which Survey selectivity block to be used.(At present,each survey have only one selectivity) 
      }
  }
  // DiIndexSelParaNumP= DiIndexSelParaNum;
  DiTestValueIAA.allocate("DiTestValueIAA");
ICHECKIAA(DiTestValueIAA);
    #if Debug_Status
        if(DiAvailIndexNum>0)
        {
           ICHECKIAA(DivIndexFleetSelChoiceFlag);
           ICHECKIAA(DivAvailIndexSelIndexOld);
        }
        cout<<"Survey Data Input Completed"<<endl;
    #endif
    if( DiTestValueIAA!=-22122)
    {
        cout<< "Input Error in Survey Index At Age  Part"<<endl;
        RECORD("Data Input Error, the TESTVALUE not Equal -22122");
        CloseRecord();
        exit(-5);
    }
    RECORD("Survey index Data End!");
    CloseIAA();
    if (Flagtimestep==1)
    { ad_comm::change_datafile_name("./InputFiles/Year/BRP_Data_Year.DAT");}
    else 
    {ad_comm::change_datafile_name("./InputFiles/Season/BRP_Data_Season.DAT");}
  DnFmax.allocate("DnFmax");
ICHECKOGD(DnFmax)
  DiBPRSelectivityFlag.allocate("DiBPRSelectivityFlag");
ICHECKOGD( DiBPRSelectivityFlag)
  //if( DiBPRSelectivityFlag) 
  //{ 
      lL2=DiSizeBinNum; 
      lL1=Flagtimestep;
  //}
  //else
  // {
  //     lL2=0;
  //     lL1=0;
  // }
  DmBPRSelectivityIni.allocate(1,lL1,1,lL2,"DmBPRSelectivityIni");
if( DiBPRSelectivityFlag==-1) ICHECKOGD(DmBPRSelectivityIni);
  DiBPRPeriod.allocate("DiBPRPeriod");
ICHECKOGD( DiBPRPeriod);
  DiBPRCalYear.allocate("DiBPRCalYear");
ICHECKOGD(DiBPRCalYear)
  DvFRInSeason.allocate(1,Flagtimestep,"DvFRInSeason");
 ICHECKOGD(DvFRInSeason)
  DivBPRBlockGrowthMatrix.allocate(1,Flagtimestep,"DivBPRBlockGrowthMatrix");
 ICHECKOGD(DivBPRBlockGrowthMatrix)
  DiTestValueOGD.allocate("DiTestValueOGD");
 ICHECKOGD(DiTestValueOGD)
    if(DiBPRCalYear >200)
         DiBPRCalYear =DiBPRCalYear -DiBeginYear+1;
    if(DiBPRCalYear<DiCalBeginYear)
          DiBPRCalYear=DiCalBeginYear;
    if(DiBPRCalYear>DiCalEndYear)
          DiBPRCalYear=DiCalEndYear;
    nTemp=sum(DvFRInSeason);
    if(nTemp>0.0)
         DvFRInSeason=DvFRInSeason/nTemp;
    else
    {
        cout<< "Error in Other Guess Data "<<endl;
        RECORD("Data Input Error, the DvFRInSeason not larger than 0.0");
        CloseRecord();
        exit(-61);
    }
    //ICHECKOGD(DvFRInSeason);
    if(DiTestValueOGD!=-22122)
    {
        cout<< "Error in Other Guess Data "<<endl;
        RECORD("Data Input Error, the TESTVALUE not Equal -22122");
        CloseRecord();
        exit(-66);
    }
    
     #if Debug_Status
        ICHECKOGD(DvFRInSeason);
        ICHECKOGD(DiBPRCalYear);
        cout<<"Biology Reference Point Data Input Completed"<<endl;
     #endif 
    if (Flagtimestep==1)
    { ad_comm::change_datafile_name("./InputFiles/Year/GrowthMatrix_Year.DAT");}
    else 
    {ad_comm::change_datafile_name("./InputFiles/Season/GrowthMatrix_Season.DAT");}
  DiGrowthMatrixFlag.allocate("DiGrowthMatrixFlag");
ICHECKGM(DiGrowthMatrixFlag);
  DiGrowthMatrixBlockNum.allocate("DiGrowthMatrixBlockNum");
ICHECKGM(DiGrowthMatrixBlockNum);
   if (1)
      {
       L1=1;
      }
   else
      {
       L1=1;
      }
  DimGrowthMatrixBlockFlag.allocate(1,DiYearNum,1,Flagtimestep+1,"DimGrowthMatrixBlockFlag");
ICHECKGM( DimGrowthMatrixBlockFlag);
  DvGrowthTimeAsYear.allocate(1,DiSeasonNum,"DvGrowthTimeAsYear");
ICHECKGM(DvGrowthTimeAsYear);//alpha
   /* for(i=1;i<=DiSeasonNum;i++)
    {//delete sgwj 2012-04-14
         if(i==1)
           DvGrowthTimeAsYear(i)=0.01;
         else if(i==2)
           DvGrowthTimeAsYear[i]=0.30;
         else if(i==3)
           DvGrowthTimeAsYear[i]=0.30;
         else if(i==4)
           DvGrowthTimeAsYear[i]=0.39;
         else
           DvGrowthTimeAsYear[i]=0.0;
    }*/
     if(!DiGrowthMatrixFlag)
     {
            lL1= DiGrowthMatrixBlockNum ;
            lL2=DiSizeBinNum;
            ii=DiSizeBinNum;
     }
        else
      {
           lL1=0;
           lL2=0;
           ii=0;
      }
  Dd3GrowthMatrix_Ini.allocate(1,lL1,1,lL2,1,ii,"Dd3GrowthMatrix_Ini");
ICHECKGM(Dd3GrowthMatrix_Ini);                                                                  //                     3  0.0  0.0 0.3 0.5  0.2
  iTestValueG.allocate("iTestValueG");
 ICHECKGM(iTestValueG)
    if( iTestValueG !=-22122)
    {
        cout<< "Error in Input Data "<<endl;
        RECORD("Data Input Error, the TESTVALUE not Equal -22122");
        CloseRecord();
        exit(-7);
    }
   #if Debug_Status
       cout<<"Growth Matrix Input Completed"<<endl;
   #endif
   CloseGM();
    if (Flagtimestep==1)
    { ad_comm::change_datafile_name("./InputFiles/Year/Prior_Year.DAT");}
    else 
    {ad_comm::change_datafile_name("./InputFiles/Season/Prior_Season.DAT");}
   if(!DiNaturalMortalityFlag)     
   {
        lL1= Flagtimestep;   
        lL2= 7; //  
   }
   else
   {
         lL1=0;
         lL2=0;
   }
  DmLorenzenA_Ini.allocate(1,lL1,1,7,"DmLorenzenA_Ini");
 if(!DiNaturalMortalityFlag)  ICHECKPIR(DmLorenzenA_Ini);
  DvLorenzenAValue.allocate(1,lL1);
  DvLorenzenALo.allocate(1,lL1);
  DvLorenzenAHi.allocate(1,lL1);
  DivLorenzenAPh.allocate(1,lL1);
  DvLorenzenASigma2.allocate(1,lL1);
  DvLorenzenASigma.allocate(1,lL1);
  DvLorenzenALambda.allocate(1,lL1);
  DivLorenzenALikelihoodFlag.allocate(1,lL1);
  DmLorenzenB_Ini.allocate(1,lL1,1,lL2,"DmLorenzenB_Ini");
 if(!DiNaturalMortalityFlag)  ICHECKGD(DmLorenzenB_Ini);
  DvLorenzenBValue.allocate(1,lL1);
  DvLorenzenBLo.allocate(1,lL1);
  DvLorenzenBHi.allocate(1,lL1);
  DivLorenzenBPh.allocate(1,lL1);
  DvLorenzenBSigma2.allocate(1,lL1);
  DvLorenzenBSigma.allocate(1,lL1);
  DvLorenzenBLambda.allocate(1,lL1);
  DivLorenzenBLikelihoodFlag.allocate(1,lL1);
   if(!DiNaturalMortalityFlag)     
   {//2011-12-15
       for(i=1;i<=lL1;i++)
       {////////////////////////////////A/////////////////////////////////////
          DvLorenzenAValue(i)=DmLorenzenA_Ini(i,1);
          DvLorenzenALo(i)= DmLorenzenA_Ini(i,2);
          DvLorenzenAHi(i)= DmLorenzenA_Ini(i,3);
          if(DmLorenzenA_Ini(i,4)>0.0)
              DivLorenzenAPh(i)= int(DmLorenzenA_Ini(i,4)+0.01);
          else
              DivLorenzenAPh(i)= int(DmLorenzenA_Ini(i,4)-0.01);
          if(DmLorenzenA_Ini(i,5)>MINPOSITIVE)
               DvLorenzenASigma2(i)=DmLorenzenA_Ini(i,5)*DmLorenzenA_Ini(i,5);//log(DmLorenzenA_Ini(i,5)*DmLorenzenA_Ini(i,5)+1.0);
          else
              DvLorenzenASigma2(i)=SIGMA2FILL;
          DvLorenzenASigma(i)=sqrt(DvLorenzenASigma2(i));
          DvLorenzenALambda(i)=DmLorenzenA_Ini(i,6);
          DivLorenzenALikelihoodFlag(i)=int(DmLorenzenA_Ini(i,7)+0.01);
          //////////////////////////B////////////////////////////////////////////
          DvLorenzenBValue(i)=DmLorenzenB_Ini(i,1);
          DvLorenzenBLo(i)= DmLorenzenB_Ini(i,2);
          DvLorenzenBHi(i)= DmLorenzenB_Ini(i,3);
          if(DmLorenzenB_Ini(i,4)>0.0)
               DivLorenzenBPh(i)= int(DmLorenzenB_Ini(i,4)+0.01);
          else
               DivLorenzenBPh(i)= int(DmLorenzenB_Ini(i,4)-0.01);
    
         if(DmLorenzenB_Ini(i,5)>MINPOSITIVE)
               DvLorenzenBSigma2(i)=log(DmLorenzenB_Ini(i,5)*DmLorenzenB_Ini(i,5)+1.0);
          else
              DvLorenzenBSigma2(i)=SIGMA2FILL;
          DvLorenzenBSigma(i)=sqrt( DvLorenzenBSigma2(i));
          DvLorenzenBLambda(i)=DmLorenzenB_Ini(i,6);
          DivLorenzenBLikelihoodFlag(i)=int(DmLorenzenB_Ini(i,7)+0.01);
       }    
   } 
   #if Debug_Status
      if(!DiNaturalMortalityFlag)
      {
        ICHECKPIR(DvLorenzenAValue);
        ICHECKPIR(DvLorenzenALo);
        ICHECKPIR(DvLorenzenAHi);
        ICHECKPIR(DivLorenzenAPh);
        ICHECKPIR(DvLorenzenASigma);
        ICHECKPIR(DvLorenzenALambda);
        ICHECKPIR(DivLorenzenALikelihoodFlag);
        ICHECKPIR(DvLorenzenBValue);
        ICHECKPIR(DvLorenzenBLo);
        ICHECKPIR(DvLorenzenBHi);
        ICHECKPIR(DivLorenzenBPh);
        ICHECKPIR(DvLorenzenBSigma);
        ICHECKPIR(DvLorenzenBLambda);
        ICHECKPIR(DivLorenzenBLikelihoodFlag);
      }
   #endif
lL1=2;
  DmRSPara_Ini.allocate(1,lL1,1,7,"DmRSPara_Ini");
ICHECKPIR(DmRSPara_Ini);                           //  belta: initial value/lower boundary/upper boundary/ phase/sd/ lambda/likelihood flag
  DvRSParaValue.allocate(1,lL1);
  DvRSParaLo.allocate(1,lL1);
  DvRSParaHi.allocate(1,lL1);
  DivRSParaPh.allocate(1,lL1);
  DvRSParaSigma2.allocate(1,lL1);
  DvRSParaSigma.allocate(1,lL1);
  DvRSParaLambda.allocate(1,lL1);
  DivRSParaLikelihoodFlag.allocate(1,lL1);
   lL1=2;
   for(i=1;i<=lL1;i++)
   {
        if(DmRSPara_Ini(i,1)>MINPOSITIVE)
              DvRSParaValue(i)=log(DmRSPara_Ini(i,1));
        else
              DvRSParaValue(i)=log(2752.0);//sgwj debug
        if(DmRSPara_Ini(i,2)>MINPOSITIVE)
                    DvRSParaLo(i)   =log(DmRSPara_Ini(i,2));
        else
                    DvRSParaLo(i)   =-15.0;
        
        if(DmRSPara_Ini(i,3)>MINPOSITIVE)
               DvRSParaHi(i)   =log(DmRSPara_Ini(i,3));
        else
               DvRSParaHi(i)   =log(20000.0); //sgwj debug 
        if(DiRSFlag ==1&& i==2)
             DivRSParaPh(i)  =-2;
        else
             DivRSParaPh(i)  =DmRSPara_Ini(i,4)>MINPOSITIVE ? int(DmRSPara_Ini(i,4)+0.01) :int(DmRSPara_Ini(i,4)-0.01);
        if(DmRSPara_Ini(i,5)> MINPOSITIVE)
             DvRSParaSigma2(i) = DmRSPara_Ini(i,5)*DmRSPara_Ini(i,5);//log( DmRSPara_Ini(i,5)*DmRSPara_Ini(i,5)+1.0);
        else
            DvRSParaSigma2(i) = SIGMA2FILL;
        DvRSParaSigma(i)=sqrt(DvRSParaSigma2(i));
        DvRSParaLambda(i) = DmRSPara_Ini(i,6);
        DivRSParaLikelihoodFlag(i)=int(DmRSPara_Ini(i,7)+0.01); 
   }
   #if Debug_Status
     ICHECKPIR(DvRSParaValue);
     ICHECKPIR(DvRSParaLo);
     ICHECKPIR(DvRSParaHi);
     ICHECKPIR(DivRSParaPh);
     ICHECKPIR(DvRSParaSigma2);
     ICHECKPIR(DvRSParaLambda);
     ICHECKPIR(DivRSParaLikelihoodFlag);
   #endif
  DmRSEnv_Ini.allocate(1,DiEnvNum,1,7,"DmRSEnv_Ini");
ICHECKPIR(DmRSPara_Ini);                           //  belta: initial value/lower boundary/upper boundary/ phase/sd/ lambda/likelihood flag
  DvRSEnvValue.allocate(1,DiEnvNum);
  DvRSEnvLo.allocate(1,DiEnvNum);
  DvRSEnvHi.allocate(1,DiEnvNum);
  DivRSEnvPh.allocate(1,DiEnvNum);
  DvRSEnvSigma2.allocate(1,DiEnvNum);
  DvRSEnvSigma.allocate(1,DiEnvNum);
  DvRSEnvLambda.allocate(1,DiEnvNum);
  DivRSEnvLikelihoodFlag.allocate(1,DiEnvNum);
   lL1=DiEnvNum;
   for(i=1;i<=lL1;i++)
   {
        if(DmRSEnv_Ini(i,1)>MINPOSITIVE)
              DvRSEnvValue(i)=log(DmRSEnv_Ini(i,1));
        else
              DvRSEnvValue(i)=log(2752.0);//sgwj debug
        if(DmRSEnv_Ini(i,2)>MINPOSITIVE)
                    DvRSEnvLo(i)   =log(DmRSEnv_Ini(i,2));
        else
                    DvRSEnvLo(i)   =-15.0;
        
        if(DmRSEnv_Ini(i,3)>MINPOSITIVE)
               DvRSEnvHi(i)   =log(DmRSEnv_Ini(i,3));
        else
               DvRSEnvHi(i)   =log(20000.0); //sgwj debug 
        if(DiRSFlag!=4&&i!=5&&i!=6)
             DivRSEnvPh(i)=-2;
        else
             DivRSEnvPh(i)  =DmRSEnv_Ini(i,4)>MINPOSITIVE ? int(DmRSEnv_Ini(i,4)+0.01) :int(DmRSEnv_Ini(i,4)-0.01);
        if(DmRSEnv_Ini(i,5)> MINPOSITIVE)
             DvRSEnvSigma2(i) = DmRSEnv_Ini(i,5)*DmRSEnv_Ini(i,5);//log( DmRSPara_Ini(i,5)*DmRSPara_Ini(i,5)+1.0);
        else
            DvRSEnvSigma2(i) = SIGMA2FILL;
        DvRSEnvSigma(i)=sqrt(DvRSEnvSigma2(i));
        DvRSEnvLambda(i) = DmRSEnv_Ini(i,6);
        DivRSEnvLikelihoodFlag(i)=int(DmRSEnv_Ini(i,7)+0.01); 
   }
   #if Debug_Status
     ICHECKPIR(DvRSEnvValue);
     ICHECKPIR(DvRSEnvLo);
     ICHECKPIR(DvRSEnvHi);
     ICHECKPIR(DivRSEnvPh);
     ICHECKPIR(DvRSEnvSigma2);
     ICHECKPIR(DvRSEnvLambda);
     ICHECKPIR(DivRSEnvLikelihoodFlag);
   #endif
  DvRecruitmentRh_Ini.allocate(1,7,"DvRecruitmentRh_Ini");
ICHECKPIR(DvRecruitmentRh_Ini)                    //    RDev(t)=sqrt( nRecruitmentR) *Rdev(t-1) -sqrt(1- nRecruitmentR) eps(t)                                          //     initial value/lower boundary/upper boundary/ phase/CV/ lambda/likelihood flag
  DvRecruitLogDevs_Ini.allocate(1,7,"DvRecruitLogDevs_Ini");
ICHECKPIR(DvRecruitLogDevs_Ini);
   DiRecruitLogDevsPhase =DvRecruitLogDevs_Ini(4)>MINPOSITIVE?int(DvRecruitLogDevs_Ini(4)+0.01):int(DvRecruitLogDevs_Ini(4)-0.01);
   DnRecruitLogDevsLambda=DvRecruitLogDevs_Ini(6);
   DiRecrDevLikelihoodFlag=int(DvRecruitLogDevs_Ini(7)+0.01);
   //for(i=DiCalBeginYear;i<=DiCalEndYear;i++)
   {     
         if(DvRecruitLogDevs_Ini(5)>MINPOSITIVE)
             DnRecruitLogDevsSigma2=DvRecruitLogDevs_Ini(5)*DvRecruitLogDevs_Ini(5);//log(DvRecruitLogDevs_Ini(5)*DvRecruitLogDevs_Ini(5)+1.0);
         else
             DnRecruitLogDevsSigma2=SIGMA2FILL;
          DnRecruitLogDevsSigma=sqrt(DnRecruitLogDevsSigma2);
   }
   //DiRecrDevLikelihoodFlag=int(DvRecruitLogDevs_Ini(7)+0.01)
  #if Debug_Status
     ICHECKPIR(DiRecruitLogDevsPhase);
     ICHECKPIR(DnRecruitLogDevsSigma);
     ICHECKPIR(DnRecruitLogDevsLambda);
     ICHECKPIR(DiRecrDevLikelihoodFlag);
  #endif
  DvRecruitLogDevsSD_Ini.allocate(1,7,"DvRecruitLogDevsSD_Ini");
ICHECKPIR(DvRecruitLogDevsSD_Ini);
  DmCPUEBiomassQE1_Ini.allocate(1,DiCPUEBlockNum,1,7,"DmCPUEBiomassQE1_Ini");
ICHECKPIR(DmCPUEBiomassQE1_Ini);
  DivCPUEQBiomassE1Ph.allocate(1,DiCPUEBlockNum);
  DvCPUEQBiomassE1Lo.allocate(1,DiCPUEBlockNum);
  DvCPUEQBiomassE1Hi.allocate(1,DiCPUEBlockNum);
  DvCPUEQBiomassE1Value.allocate(1,DiCPUEBlockNum);
  DvCPUEQBiomassE1Sigma2.allocate(1,DiCPUEBlockNum);
  DvCPUEQBiomassE1Sigma.allocate(1,DiCPUEBlockNum);
  DvCPUEQBiomassE1Lambda.allocate(1,DiCPUEBlockNum);
  DivCPUEQBiomassE1LikelihoodFlag.allocate(1,DiCPUEBlockNum);
  
   ii=0;
   for(i=1;i<=DiCPUEBlockNum;i++)
   {
                            DvCPUEQBiomassE1Value(i)=DmCPUEBiomassQE1_Ini(i,1);
                            DvCPUEQBiomassE1Lo(i)   =DmCPUEBiomassQE1_Ini(i,2);
                            DvCPUEQBiomassE1Hi(i)   =DmCPUEBiomassQE1_Ini(i,3);
                            if(DivCPUEBlockAvailFlag(i))
                                DivCPUEQBiomassE1Ph(i)=DmCPUEBiomassQE1_Ini(i,4)>0.0?int(DmCPUEBiomassQE1_Ini(i,4)+0.1):int(DmCPUEBiomassQE1_Ini(i,4)-0.1);
                            else
                                DivCPUEQBiomassE1Ph(i)=-2;
                            if(DmCPUEBiomassQE1_Ini(i,5)>MINPOSITIVE)
                                DvCPUEQBiomassE1Sigma2(i)=DmCPUEBiomassQE1_Ini(i,5)*DmCPUEBiomassQE1_Ini(i,5);//log(DmCPUEBiomassQE1_Ini(i,5)*DmCPUEBiomassQE1_Ini(i,5)+1.0);
                            else
                                DvCPUEQBiomassE1Sigma2(i)=SIGMA2FILL;
                            DvCPUEQBiomassE1Sigma(i)=sqrt( DvCPUEQBiomassE1Sigma2(i));
                            DvCPUEQBiomassE1Lambda(i)=DmCPUEBiomassQE1_Ini(i,6);
                            DivCPUEQBiomassE1LikelihoodFlag(i)=int(DmCPUEBiomassQE1_Ini(i,7)+0.1);
   }
   
    #if Debug_Status
      ICHECKPIR(DvCPUEQBiomassE1Value);
      ICHECKPIR(DvCPUEQBiomassE1Lo);
      ICHECKPIR(DvCPUEQBiomassE1Hi);
      ICHECKPIR(DivCPUEQBiomassE1Ph);
      ICHECKPIR(DvCPUEQBiomassE1Sigma);
      ICHECKPIR(DvCPUEQBiomassE1Lambda);
      ICHECKPIR(DivCPUEQBiomassE1LikelihoodFlag);
      //exit(-99);
    #endif 
  DmFleetSelParaData_Ini.allocate(1,DiFleetSelParaNumIni,1,7,"DmFleetSelParaData_Ini");
 ICHECKPIR( DmFleetSelParaData_Ini);
  DvFleetSelPara_Ini.allocate(1,DiFleetSelParaNumIni);
  DvFleetSelParaBLo.allocate(1,DiFleetSelParaNumIni);
  DvFleetSelParaBHi.allocate(1,DiFleetSelParaNumIni);
  DivFleetSelParaPh.allocate(1,DiFleetSelParaNumIni);
  DvFleetSelParaLambda.allocate(1,DiFleetSelParaNumIni);
  DvFleetSelParaSigma2.allocate(1,DiFleetSelParaNumIni);
  DvFleetSelParaSigma.allocate(1,DiFleetSelParaNumIni);
  DivFleetSelParaLikelihoodFlag.allocate(1,DiFleetSelParaNumIni);
  // lL2=1;
   for(i=1;i<=DiFleetSelBlockNum;i++)
   {
       lL1=DivFleetSelParaIndex_Ini(i)+1;//(i-1)*(DiSizeBinNum+2+4)+1;
       switch( DivFleetSelBlockFlag(i))
       {
         case 1:
           for(j=1;j<=DiSizeBinNum;j++)
           {
              DvFleetSelPara_Ini(lL1)= DmFleetSelParaData_Ini(lL1,1);
              DvFleetSelParaBLo(lL1)=DmFleetSelParaData_Ini(lL1,2);
              DvFleetSelParaBHi(lL1)=DmFleetSelParaData_Ini(lL1,3);
              if(!DivFleetSelBlockAvailFlag(i))
                  DivFleetSelParaPh(lL1)=-2;
              else
                  DivFleetSelParaPh(lL1)= DmFleetSelParaData_Ini(lL1,4)> 0.0 ? int( DmFleetSelParaData_Ini(lL1,4)+0.01): int( DmFleetSelParaData_Ini(lL1,4)-0.01);
              if( DmFleetSelParaData_Ini(lL1,5)> MINPOSITIVE)
                  DvFleetSelParaSigma2(lL1)=DmFleetSelParaData_Ini(lL1,5)* DmFleetSelParaData_Ini(lL1,5);//log( DmFleetSelParaData_Ini(lL1,5)* DmFleetSelParaData_Ini(lL1,5)+1.0);
              else
                  DvFleetSelParaSigma2(lL1)=SIGMA2FILL;
              DvFleetSelParaSigma(lL1)=sqrt(DvFleetSelParaSigma2(lL1));
              DvFleetSelParaLambda(lL1)= DmFleetSelParaData_Ini(lL1,6); 
              DivFleetSelParaLikelihoodFlag(lL1)=int(DmFleetSelParaData_Ini(lL1,7)+0.1); 
            //  lL2++;
              lL1++;
           }
           break;
         case 2:
          // lL1=lL1+DiSizeBinNum;
           for(j=1;j<=2;j++)
           {
              DvFleetSelPara_Ini(lL1)= DmFleetSelParaData_Ini(lL1,1);
              DvFleetSelParaBLo(lL1)=DmFleetSelParaData_Ini(lL1,2);
              DvFleetSelParaBHi(lL1)=DmFleetSelParaData_Ini(lL1,3);
              if(!DivFleetSelBlockAvailFlag(i))
                  DivFleetSelParaPh(lL1)=-2;
              else
                  DivFleetSelParaPh(lL1)= DmFleetSelParaData_Ini(lL1,4)>0.0 ? int( DmFleetSelParaData_Ini(lL1,4)+0.01) : int(DmFleetSelParaData_Ini(lL1,4) -0.01);
             
              if( DmFleetSelParaData_Ini(lL1,5)> MINPOSITIVE)
                  DvFleetSelParaSigma2(lL1)= DmFleetSelParaData_Ini(lL1,5)* DmFleetSelParaData_Ini(lL1,5);//log( DmFleetSelParaData_Ini(lL1,5)* DmFleetSelParaData_Ini(lL1,5)+1.0);
              else
                  DvFleetSelParaSigma2(lL1)=SIGMA2FILL;
              DvFleetSelParaSigma(lL1)=sqrt(DvFleetSelParaSigma2(lL1));
   
              DvFleetSelParaLambda(lL1)= DmFleetSelParaData_Ini(lL1,6); 
              DivFleetSelParaLikelihoodFlag(lL1)=int(DmFleetSelParaData_Ini(lL1,7)+0.1); //DvFleetSelParaLikeliFlag(lL2)=DmFleetSelParaData_Ini(lL1,7); 
             // lL2++;
              lL1++;
           }
           break;
          case 3:
          // lL1=lL1+DiSizeBinNum+2;
           for(j=1;j<=4;j++)
           {
              DvFleetSelPara_Ini(lL1)= DmFleetSelParaData_Ini(lL1,1);
              DvFleetSelParaBLo(lL1)= DmFleetSelParaData_Ini(lL1,2);
              DvFleetSelParaBHi(lL1)= DmFleetSelParaData_Ini(lL1,3);
              if(!DivFleetSelBlockAvailFlag(i))
                  DivFleetSelParaPh(lL1)=-2;
              else
                  DivFleetSelParaPh(lL1)=DmFleetSelParaData_Ini(lL1,4)> 0.0 ? int(DmFleetSelParaData_Ini(lL1,4)+0.01) :int(DmFleetSelParaData_Ini(lL1,4)-0.01);
             
              if( DmFleetSelParaData_Ini(lL1,5)> MINPOSITIVE)
                  DvFleetSelParaSigma2(lL1)=DmFleetSelParaData_Ini(lL1,5)*DmFleetSelParaData_Ini(lL1,5);//log(DmFleetSelParaData_Ini(lL1,5)*DmFleetSelParaData_Ini(lL1,5)+1.0);
              else
                  DvFleetSelParaSigma2(lL1)=SIGMA2FILL;
              DvFleetSelParaSigma(lL1)=sqrt(DvFleetSelParaSigma2(lL1));
              DvFleetSelParaLambda(lL1)= DmFleetSelParaData_Ini(lL1,6); 
              DivFleetSelParaLikelihoodFlag(lL1)=int(DmFleetSelParaData_Ini(lL1,7)+0.1); //DvFleetSelParaLikeliFlag(lL2)=DmFleetSelParaData_Ini(lL1,7); 
            //  lL2++;
              lL1++;
           }
           break;
       }
    }
  #if Debug_Status
     ICHECKPIR(DvFleetSelPara_Ini);
     ICHECKPIR(DvFleetSelParaBLo);
     ICHECKPIR(DvFleetSelParaBHi);
     //ICHECKDEB(DivFleetSelParaPhase);
     ICHECKPIR(DvFleetSelParaSigma);
     ICHECKPIR(DvFleetSelParaLambda);
  #endif
   RECORD("Fishery selectivity");
   //
  DmIndexQE1Block_Ini.allocate(1,DiIndexQBlockNum,1,7,"DmIndexQE1Block_Ini");
ICHECKPIR(DmIndexQE1Block_Ini)
  DvQBiomassE1Value.allocate(1,DiIndexQBlockNum);
  DvQBiomassE1Lo.allocate(1,DiIndexQBlockNum);
  DvQBiomassE1Hi.allocate(1,DiIndexQBlockNum);
  DivQBiomassE1Ph.allocate(1,DiIndexQBlockNum);
  DvQBiomassE1Sigma2.allocate(1,DiIndexQBlockNum);
  DvQBiomassE1Sigma.allocate(1,DiIndexQBlockNum);
  DvQBiomassE1Lambda.allocate(1,DiIndexQBlockNum);
  DivQBiomassE1LikelihoodFlag.allocate(1,DiIndexQBlockNum);
   for(i=1;i<=DiIndexQBlockNum;i++)
   {
        //  j=DivAvailIndexQMirrorOld(i);
          DvQBiomassE1Value(i)   =  DmIndexQE1Block_Ini(i,1);
          DvQBiomassE1Lo(i)      =  DmIndexQE1Block_Ini(i,2);
          DvQBiomassE1Hi(i)      =  DmIndexQE1Block_Ini(i,3);
          
          if( DivAvailIndexQFlag(i))
                DivQBiomassE1Ph(i)     =  DmIndexQE1Block_Ini(i,4)>MINPOSITIVE ? int(DmIndexQE1Block_Ini(i,4)+0.01) :int(DmIndexQE1Block_Ini(i,4)-0.01);
          else
                DivQBiomassE1Ph(i)     =-2;
          if( DmIndexQE1Block_Ini(i,5)>MINPOSITIVE)
                   DvQBiomassE1Sigma2(i)  = DmIndexQE1Block_Ini(i,5)*DmIndexQE1Block_Ini(i,5);//log( DmIndexQE1Block_Ini(i,5)*DmIndexQE1Block_Ini(i,5)+1.0);
          else
                   DvQBiomassE1Sigma2(i) = SIGMA2FILL;
         
          DvQBiomassE1Sigma(i)   =  sqrt(DvQBiomassE1Sigma2(i));
          DvQBiomassE1Lambda(i)  =  DmIndexQE1Block_Ini(i,6); 
          DivQBiomassE1LikelihoodFlag(i)  = int( DmIndexQE1Block_Ini(i,7)+0.1); 
          //DivAvailIndexQChoiceFlag
  }
  #if Debug_Status
      ICHECKPIR(DvQBiomassE1Value);
      ICHECKPIR(DvQBiomassE1Lo);
      ICHECKPIR(DvQBiomassE1Hi);
      ICHECKPIR(DivQBiomassE1Ph);
      ICHECKPIR(DvQBiomassE1Sigma);
      ICHECKPIR(DvQBiomassE1Lambda);
      ICHECKPIR(DivQBiomassE1LikelihoodFlag);
  #endif
  DmIndexSelParaData_Ini.allocate(1,DiIndexSelParaNumIni,1,7,"DmIndexSelParaData_Ini");
 ICHECKPIR( DmIndexSelParaData_Ini);
  DvIndexSelPara_Ini.allocate(1,DiIndexSelParaNumIni);
  DvIndexSelParaBLo.allocate(1,DiIndexSelParaNumIni);
  DvIndexSelParaBHi.allocate(1,DiIndexSelParaNumIni);
  DivIndexSelParaPh.allocate(1,DiIndexSelParaNumIni);
  DvIndexSelParaLambda.allocate(1,DiIndexSelParaNumIni);
  DvIndexSelParaSigma2.allocate(1,DiIndexSelParaNumIni);
  DvIndexSelParaSigma.allocate(1,DiIndexSelParaNumIni);
  DivIndexSelParaLikelihoodFlag.allocate(1,DiIndexSelParaNumIni);
  DivEstIndexPropFlag.allocate(1,DiAvailIndexNum);
   DivEstIndexPropFlag.initialize();
   DvIndexSelParaLambda.initialize(); 
   
  // DivIndexSelParaPhase=-2;//  
   for(i=1;i<=DiIndexSelParaNumIni;i++)
   {//DivIndexSelBlockAvailFlag
                              lL1=lL2=i; 
                              DvIndexSelPara_Ini(lL2)=  DmIndexSelParaData_Ini(lL1,1);
                              DvIndexSelParaBLo(lL2) =  DmIndexSelParaData_Ini(lL1,2);
                              DvIndexSelParaBHi(lL2) =  DmIndexSelParaData_Ini(lL1,3);
                              DivIndexSelParaPh(lL2)=DmIndexSelParaData_Ini(lL1,4) > MINPOSITIVE? int( DmIndexSelParaData_Ini(lL1,4)+0.01) :int( DmIndexSelParaData_Ini(lL1,4)-0.01);
                              //if(DivIndexSelParaPh(lL2)>=1)
                              //       DivEstIndexPropFlag(k)=1;
                              if(DmIndexSelParaData_Ini(lL1,5)> MINPOSITIVE)
                                   DvIndexSelParaSigma2(lL2)=DmIndexSelParaData_Ini(lL1,5)* DmIndexSelParaData_Ini(lL1,5);//log( DmIndexSelParaData_Ini(lL1,5)* DmIndexSelParaData_Ini(lL1,5)+1.0);
                              else
                                   DvIndexSelParaSigma2(lL2)=SIGMA2FILL;
                              DvIndexSelParaSigma(lL2)=sqrt(DvIndexSelParaSigma2(lL2));
                              DvIndexSelParaLambda(lL2)=  DmIndexSelParaData_Ini(lL1,6); 
                              DivIndexSelParaLikelihoodFlag(lL2)=int( DmIndexSelParaData_Ini(lL1,7)+0.01);
                              
             
   }
   
   for(i=1;i<=DiIndexSelBlockNum;i++)
   {//DivIndexSelBlockAvailFlag
        if(!DivIndexSelBlockAvailFlag(i))
        {//DivIndexSelFlag_Ini
             lL1= DivIndexSelParaIndex_Ini(i);
             //cout<<"i"<<i<<"lL1:"<<lL1<<endl;
             switch(DivIndexSelFlag_Ini(i))
             {
               case 1:
                  for(k=1;k<=DiSizeBinNum;k++)
                         DivIndexSelParaPh(lL1+k)=-2;
                  break;
               case 2:
                   DivIndexSelParaPh(lL1+1)=-2;
                   DivIndexSelParaPh(lL1+2)=-2;
                   break;
               case 3:
                   DivIndexSelParaPh(lL1+1)=-2;
                   DivIndexSelParaPh(lL1+2)=-2;
                   DivIndexSelParaPh(lL1+3)=-2;
                   DivIndexSelParaPh(lL1+4)=-2;
                   break;
             }
        }
   }
   ii=0;
   for(i=1;i<=DiIndexNum;i++)
   {
        if(DivAvailIndexFlag(i))
        {
             ii++;
             for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
             {
                jj=DimIndexSelBlockFlag(j,i+1);      //block;
                lL1= DivIndexSelParaIndex_Ini(jj);
                switch(DivIndexSelFlag_Ini(i))
                {
                  case 1:
                     for(k=1;k<=DiSizeBinNum;k++)
                     {
                         if(DivIndexSelParaPh(lL1+k)>0)
					;
                               DivEstIndexPropFlag(ii)=1;
                     }
                     break;
                  case 2:
                    for(k=1;k<=2;k++)
                    {
                         if(DivIndexSelParaPh(lL1+k)>0)
					;
                               DivEstIndexPropFlag(ii)=1;
                    }
                    break;
                  case 3:
                   for(k=1;k<=4;k++)
                    {
                         if(DivIndexSelParaPh(lL1+k)>0)
					;
                               DivEstIndexPropFlag(ii)=1;
                    }
                   break;
                }
                if(DivEstIndexPropFlag(ii)==1)
                   break;
             }
        }
   }
   #if Debug_Status
      ICHECKPIR( DvIndexSelPara_Ini);
      ICHECKPIR( DvIndexSelParaBLo);
      ICHECKPIR( DvIndexSelParaBHi);
      //
      ICHECKPIR( DvIndexSelParaSigma);
      ICHECKPIR( DvIndexSelParaLambda);
      ICHECKPIR( DivEstIndexPropFlag);
   #endif
   lL2=7;
   if( DiNYear1ChoiceFlag==0)
   {
       lL1=DiSizeBinNum;
       DiNYear1EstParaNum =DiSizeBinNum;
    }
   if( DiNYear1ChoiceFlag==1)
   {
       lL1=1;
       DiNYear1EstParaNum =1;
    }
   if( DiNYear1ChoiceFlag==2)
   {
       lL1=1;
       DiNYear1EstParaNum =1;
    }
   if( DiNYear1ChoiceFlag==7)
   {
       lL1=9;
       DiNYear1EstParaNum =9;
    }
   if( DiNYear1ChoiceFlag==3)
   {
       lL1=3;
       DiNYear1EstParaNum =3;
    }
   if( DiNYear1ChoiceFlag==4)
   {
       lL1=3;
       DiNYear1EstParaNum =3;
    }
   if( DiNYear1ChoiceFlag==5)
   {
       lL1=3;
       DiNYear1EstParaNum =3;
    }
   if( DiNYear1ChoiceFlag==6)
   {
       lL1=3;
       DiNYear1EstParaNum =3;
    }
   //else
   //{
   //    DiNYear1EstParaNum =3;
   //    lL1=3;
   //}
  //cout<<"DiNYear1ChoiceFlag"<<DiNYear1ChoiceFlag<<endl;
  //cout<<"DiNYear1EstParaNum:"<<lL1<<DiNYear1EstParaNum<<endl;
  //exit(-2);
  DmNYear1NPia_Ini.allocate(1,lL1,1,lL2,"DmNYear1NPia_Ini");
ICHECKPIR(DmNYear1NPia_Ini);
  DvNYear1ParaInital.allocate(1,DiNYear1EstParaNum);
  DvNYear1ParaLo.allocate(1,DiNYear1EstParaNum);
  DvNYear1ParaHi.allocate(1,DiNYear1EstParaNum);
  DivNYear1ParaPh.allocate(1,DiNYear1EstParaNum);
  DvNYear1ParaSigma2.allocate(1,DiNYear1EstParaNum);
  DvNYear1ParaSigma.allocate(1,DiNYear1EstParaNum);
  DvNYear1ParaLambda.allocate(1,DiNYear1EstParaNum);
  DivNYear1ParaLikelihoodFlag.allocate(1,DiNYear1EstParaNum);
  
    for( i=1;i<=DiNYear1EstParaNum;i++)
    {//DimFleetFirstYear(1,DiFleetNum,1,DiSeasonNum) 
               DvNYear1ParaInital(i) =  DmNYear1NPia_Ini(i,1);
               DvNYear1ParaLo(i)     =  DmNYear1NPia_Ini(i,2);
               DvNYear1ParaHi(i)     =  DmNYear1NPia_Ini(i,3);
               if(DiNYear1ChoiceFlag==1||DiNYear1ChoiceFlag==2)
               {
                  if(i==1)
                     DivNYear1ParaPh(i)    =  DmNYear1NPia_Ini(i,4)> MINPOSITIVE ? int( DmNYear1NPia_Ini(i,4)+0.001):int( DmNYear1NPia_Ini(i,4)-0.01);
                  else
                     DivNYear1ParaPh(i)    =  -2;
               }
               else
               {
                    
                    DivNYear1ParaPh(i)    =  DmNYear1NPia_Ini(i,4)> MINPOSITIVE ? int( DmNYear1NPia_Ini(i,4)+0.001):int( DmNYear1NPia_Ini(i,4)-0.01);
               }
               
               if( DmNYear1NPia_Ini(i,5)>MINPOSITIVE)
                  DvNYear1ParaSigma2(i)=DmNYear1NPia_Ini(i,5)* DmNYear1NPia_Ini(i,5);//log( DmNYear1NPia_Ini(i,5)* DmNYear1NPia_Ini(i,5)+1.0);
               else
                  DvNYear1ParaSigma2(i) =SIGMA2FILL;
               
               DvNYear1ParaSigma(i)    =sqrt(DvNYear1ParaSigma2(i));
               DvNYear1ParaLambda(i)   =DmNYear1NPia_Ini(i,6);
               DivNYear1ParaLikelihoodFlag(i)=int( DmNYear1NPia_Ini(i,7)+0.01);
    }
   ////////////////////////////////////////////////////////
  #if Debug_Status
   ICHECKPIR(DvNYear1ParaInital);
   ICHECKPIR(DvNYear1ParaLo);
   ICHECKPIR(DvNYear1ParaHi);
   ICHECKPIR(DivNYear1ParaPh);
   ICHECKPIR(DvNYear1ParaSigma);
   ICHECKPIR(DvNYear1ParaLambda);
   ICHECKPIR(DivNYear1ParaLikelihoodFlag);
  #endif
  ////////////////////////////////////////////////////
lL1=DiFleetNum*Flagtimestep;                    // initial value/lo/up/phase/sd/lambda/likelihood flag
  DmFYear1_Ini.allocate(1,lL1,1,7,"DmFYear1_Ini");
ICHECKPIR(DmFYear1_Ini);
  DvLogFYear1Sigma2.allocate(1,lL1);
  DvLogFYear1Sigma.allocate(1,lL1);
  DvLogFYear1Ini.allocate(1,lL1);
  DvLogFYear1Lambda.allocate(1,lL1);
  DivLogFYear1LikelihoodFlag.allocate(1,lL1);
  DvLogFYear1Lo.allocate(1,lL1);
  DvLogFYear1Hi.allocate(1,lL1);
  DivLogFYear1Ph.allocate(1,lL1);
   for(i=1;i<=DiFleetNum;i++)
   {//  DimFleetFirstYear(1,DiFleetNum,1,DiSeasonNum) 
     for(k=1;k<=Flagtimestep;k++)
     {
        iL1=(i-1)*Flagtimestep+k;
        if(DmFYear1_Ini(iL1,1)>0.0)
             DvLogFYear1Ini(iL1)=log(DmFYear1_Ini(iL1,1));//???????????
        else
             DvLogFYear1Ini(iL1)=log(0.27/4.0);//??????????
       if(DmFYear1_Ini(iL1,2)>0.0)
               DvLogFYear1Lo(iL1)=log(DmFYear1_Ini(iL1,2));//?????????????
       else
               DvLogFYear1Lo(iL1)=-25.0;
       if(DmFYear1_Ini(iL1,3)>0.0)
             DvLogFYear1Hi(iL1)=log(DmFYear1_Ini(iL1,3));//????????????
       else
             DvLogFYear1Hi(iL1)=log(DnFmax);
        //sgwj1 debug  CV -log scale
        if(DimFleetFirstYear(i,k)==0)
        {
             DivLogFYear1Ph(iL1)=-2;
        }
        else
            DivLogFYear1Ph(iL1)=DmFYear1_Ini(iL1,4)>0.0 ? int(DmFYear1_Ini(iL1,4)+0.01): int(DmFYear1_Ini(iL1,4)-0.01);
        if(DmFYear1_Ini(iL1,5)>0.0)
            DvLogFYear1Sigma2(iL1)= DmFYear1_Ini(iL1,5)*DmFYear1_Ini(iL1,5);//log(DmFYear1_Ini(iL1,5)*DmFYear1_Ini(iL1,5)+1.0);
        else
            DvLogFYear1Sigma2(iL1)=SIGMA2FILL;
        DvLogFYear1Sigma(iL1)=sqrt(DvLogFYear1Sigma2(iL1));
        DvLogFYear1Lambda(iL1)=DmFYear1_Ini(iL1,6);
        DivLogFYear1LikelihoodFlag(iL1)=int(DmFYear1_Ini(iL1,7)+0.01);
     }
   }
  #if Debug_Status
   ICHECKPIR(DvLogFYear1Ini);
   ICHECKPIR(DvLogFYear1Lo);
   ICHECKPIR(DvLogFYear1Hi);
   ICHECKPIR(DivLogFYear1Ph);
   ICHECKPIR(DvLogFYear1Sigma);
   ICHECKPIR(DvLogFYear1Lambda);
   ICHECKPIR(DivLogFYear1LikelihoodFlag);
  #endif
lL1=DiFleetNum*Flagtimestep;  
  DmFLogDevs_Ini.allocate(1,lL1,1,7,"DmFLogDevs_Ini");
ICHECKPIR(DmFLogDevs_Ini);
  DvFLogDevsSigma2.allocate(1,lL1);
  DvFLogDevsSigma.allocate(1,lL1);
  DvFLogDevsLambda.allocate(1,lL1);
  DivFLogDevsLikelihoodFlag.allocate(1,lL1);
   for(i=1;i<=lL1;i++)
   {  
           if( DmFLogDevs_Ini(i,5)>MINPOSITIVE)
                             DvFLogDevsSigma2(i)=DmFLogDevs_Ini(i,5)* DmFLogDevs_Ini(i,5);//log( DmFLogDevs_Ini(i,5)* DmFLogDevs_Ini(i,5)+1.0);
           else
                             DvFLogDevsSigma2(i)=SIGMA2FILL; 
           DvFLogDevsSigma(i)=sqrt(DvFLogDevsSigma2(i));
           DvFLogDevsLambda(i)=DmFLogDevs_Ini(i,6);
           DivFLogDevsLikelihoodFlag(i)=int(DmFLogDevs_Ini(i,7)+0.1);    
    }
  #if Debug_Status
   
   ICHECKPIR(DvFLogDevsSigma);
   ICHECKPIR(DvFLogDevsLambda);
   ICHECKPIR(DivFLogDevsLikelihoodFlag);
  
  #endif             
lL1=DiFleetNum*(DiCalEndYear-DiCalBeginYear+1)*Flagtimestep;
  DvFLogDevsValue.allocate(1,lL1);
  DvFLogDevsLo.allocate(1,lL1);
  DvFLogDevsHi.allocate(1,lL1);
  DivFLogDevsPh.allocate(1,lL1);
   for(i=1;i<=DiFleetNum;i++)
   {  
      for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
      {
          for(k=1;k<=Flagtimestep;k++)
          {
               iL1=(i-1)*(DiCalEndYear-DiCalBeginYear+1)*Flagtimestep + (j-DiCalBeginYear)*Flagtimestep +k;
               iL2=(i-1)*Flagtimestep +k;
               lL2=(j-DiCalBeginYear)*Flagtimestep+k;
               if(DimFleetFZeroIndex(i,lL2)!=0)
               {
                     DvFLogDevsValue(iL1)=0.0;      //DmFLogDevs_Ini(i,1);
                     DvFLogDevsLo(iL1)=log(MINPOSITIVE);        //DmFLogDevs_Ini(i,2);
                     DvFLogDevsHi(iL1)=log(DnFmax/MINPOSITIVE);//DmFLogDevs_Ini(i,3);
                    // DivFLogDevsPh(iL1)=DmFLogDevs_Ini(iL2,4)>MINPOSITIVE ? int(DmFLogDevs_Ini(iL2,4)+0.01) : int(DmFLogDevs_Ini(iL2,4)-0.01) ;
                     DivFLogDevsPh(iL1)=1;
               }
               else
               {
                    DvFLogDevsValue(iL1)=log(MINPOSITIVE);
                    DvFLogDevsLo(iL1)=log(MINPOSITIVE);
                    DvFLogDevsHi(iL1)=log(DnFmax/MINPOSITIVE);
                    DivFLogDevsPh(iL1)=-2;
               }
          }
       }   
    }
   
   lL1=DiFleetNum*(DiCalEndYear-DiCalBeginYear+1)*Flagtimestep;
   for(i=1;i<=DiFleetNum;i++)
   {
             for(k=1;k<=Flagtimestep;k++)
             {
                j= DimFleetFirstYear(i,k);//j==0
                if(j==0)
                 continue;
                 //  k= DivFleetFirstSeason(i);
                iL1=(i-1)*(DiCalEndYear-DiCalBeginYear+1)*Flagtimestep + (j-DiCalBeginYear)*Flagtimestep +k;
                if( iL1>=1 && iL1<=lL1)
                {
                    DvFLogDevsValue(iL1)=0.0;
                    DvFLogDevsLo(iL1)=log(MINPOSITIVE);
                    DvFLogDevsHi(iL1)=log(DnFmax/MINPOSITIVE);
                    DivFLogDevsPh(iL1)=-2;
                 }  
             }
   }     
  
  #if Debug_Status
      ICHECKPIR(DvFLogDevsValue);
      ICHECKPIR(DvFLogDevsLo);
      ICHECKPIR(DvFLogDevsHi);
      ICHECKPIR(DivFLogDevsPh);
  #endif
    lL1= 5  * DiGrowthMatrixBlockNum ;
    lL2=7;
    #if Debug_Status
    #endif
  DmVBGFParas_Ini.allocate(1,lL1,1,7,"DmVBGFParas_Ini");
 ICHECKPIR( DmVBGFParas_Ini);                           //  Linf    V    LSD     VSD    Rho(Linf,V)
 lL1=DiGrowthMatrixBlockNum;
  DvVBGFLinfValue.allocate(1,lL1);
  DvVBGFLinfLo.allocate(1,lL1);
  DvVBGFLinfHi.allocate(1,lL1);
  DivVBGFLinfPh.allocate(1,lL1);
  DvVBGFLinfSigma2.allocate(1,lL1);
  DvVBGFLinfSigma.allocate(1,lL1);
  DvVBGFLinfLambda.allocate(1,lL1);
  DivVBGFLinfLikelihoodFlag.allocate(1,lL1);
  DvVBGFVValue.allocate(1,lL1);
  DvVBGFVLo.allocate(1,lL1);
  DvVBGFVHi.allocate(1,lL1);
  DivVBGFVPh.allocate(1,lL1);
  DvVBGFVSigma2.allocate(1,lL1);
  DvVBGFVSigma.allocate(1,lL1);
  DvVBGFVLambda.allocate(1,lL1);
  DivVBGFVLikelihoodFlag.allocate(1,lL1);
  DvVBGFLSDValue.allocate(1,lL1);
  DvVBGFLSDLo.allocate(1,lL1);
  DvVBGFLSDHi.allocate(1,lL1);
  DivVBGFLSDPh.allocate(1,lL1);
  DvVBGFLSDSigma2.allocate(1,lL1);
  DvVBGFLSDSigma.allocate(1,lL1);
  DvVBGFLSDLambda.allocate(1,lL1);
  DivVBGFLSDLikelihoodFlag.allocate(1,lL1);
  DvVBGFVSDValue.allocate(1,lL1);
  DvVBGFVSDLo.allocate(1,lL1);
  DvVBGFVSDHi.allocate(1,lL1);
  DivVBGFVSDPh.allocate(1,lL1);
  DvVBGFVSDSigma2.allocate(1,lL1);
  DvVBGFVSDSigma.allocate(1,lL1);
  DvVBGFVSDLambda.allocate(1,lL1);
  DivVBGFVSDLikelihoodFlag.allocate(1,lL1);
  DvVBGFLVRhoValue.allocate(1,lL1);
  DvVBGFLVRhoLo.allocate(1,lL1);
  DvVBGFLVRhoHi.allocate(1,lL1);
  DivVBGFLVRhoPh.allocate(1,lL1);
  DvVBGFLVRhoSigma2.allocate(1,lL1);
  DvVBGFLVRhoSigma.allocate(1,lL1);
  DvVBGFLVRhoLambda.allocate(1,lL1);
  DivVBGFLVRhoLikelihoodFlag.allocate(1,lL1);
   // if(DiGrowthMatrixFlag)
   // {
      for(i=1;i<=DiGrowthMatrixBlockNum ;i++)
      {//Linf
          DvVBGFLinfValue(i)=DmVBGFParas_Ini((i-1)*5+1,1);
          DvVBGFLinfLo(i)=DmVBGFParas_Ini((i-1)*5+1,2);
          DvVBGFLinfHi(i)=DmVBGFParas_Ini((i-1)*5+1,3);
          DivVBGFLinfPh(i)=DmVBGFParas_Ini((i-1)*5+1,4)>MINPOSITIVE ?int(DmVBGFParas_Ini((i-1)*5+1,4)+0.01) :int(DmVBGFParas_Ini((i-1)*5+1,4)-0.01) ;
          if(DmVBGFParas_Ini((i-1)*5+3,5)>MINPOSITIVE)
                    DvVBGFLinfSigma2(i)=DmVBGFParas_Ini((i-1)*5+1,5)*DmVBGFParas_Ini((i-1)*5+1,5);//log(DmVBGFParas_Ini((i-1)*5+1,5)*DmVBGFParas_Ini((i-1)*5+1,5)+1.0);
          else
                    DvVBGFLinfSigma2(i)=SIGMA2FILL;
          DvVBGFLinfSigma(i)=sqrt(DvVBGFLinfSigma2(i));
          DvVBGFLinfLambda(i)=DmVBGFParas_Ini((i-1)*5+1,6);
          DivVBGFLinfLikelihoodFlag(i)=int(DmVBGFParas_Ini((i-1)*5+1,7)+0.1);
          //V
          DvVBGFVValue(i)=DmVBGFParas_Ini((i-1)*5+2,1);
          DvVBGFVLo(i)=DmVBGFParas_Ini((i-1)*5+2,2);
          DvVBGFVHi(i)=DmVBGFParas_Ini((i-1)*5+2,3);
          DivVBGFVPh(i)=DmVBGFParas_Ini((i-1)*5+2,4)>MINPOSITIVE ?int(DmVBGFParas_Ini((i-1)*5+2,4)+0.01) :int(DmVBGFParas_Ini((i-1)*5+2,4)-0.01) ;
          if(DmVBGFParas_Ini((i-1)*5+2,5)>MINPOSITIVE)
                    DvVBGFVSigma2(i)=DmVBGFParas_Ini((i-1)*5+2,5)*DmVBGFParas_Ini((i-1)*5+2,5);//log(DmVBGFParas_Ini((i-1)*5+2,5)*DmVBGFParas_Ini((i-1)*5+2,5)+1.0);
          else
                    DvVBGFVSigma2(i)=SIGMA2FILL;
          DvVBGFVSigma(i)=sqrt(DvVBGFVSigma2(i));
          DvVBGFVLambda(i)=DmVBGFParas_Ini((i-1)*5+2,6);
          DivVBGFVLikelihoodFlag(i)=int(DmVBGFParas_Ini((i-1)*5+2,7)+0.1);
          //LSD
          DvVBGFLSDValue(i)=DmVBGFParas_Ini((i-1)*5+3,1);
          DvVBGFLSDLo(i)=DmVBGFParas_Ini((i-1)*5+3,2);
          DvVBGFLSDHi(i)=DmVBGFParas_Ini((i-1)*5+3,3);
          DivVBGFLSDPh(i)=DmVBGFParas_Ini((i-1)*5+3,4)>MINPOSITIVE ?int(DmVBGFParas_Ini((i-1)*5+3,4)+0.01) :int(DmVBGFParas_Ini((i-1)*5+3,4)-0.01) ;
          if(DmVBGFParas_Ini((i-1)*5+3,5)>MINPOSITIVE)
                    DvVBGFLSDSigma2(i)=DmVBGFParas_Ini((i-1)*5+3,5)*DmVBGFParas_Ini((i-1)*5+3,5);//log(DmVBGFParas_Ini((i-1)*5+3,5)*DmVBGFParas_Ini((i-1)*5+3,5)+1.0);
          else
                    DvVBGFLSDSigma2(i)=SIGMA2FILL;
          DvVBGFLSDSigma(i)=sqrt(DvVBGFLSDSigma2(i));
          DvVBGFLSDLambda(i)=DmVBGFParas_Ini((i-1)*5+3,6);
          DivVBGFLSDLikelihoodFlag(i)=int(DmVBGFParas_Ini((i-1)*5+3,7)+0.1);
          //VSD
          DvVBGFVSDValue(i)=DmVBGFParas_Ini((i-1)*5+4,1);
          DvVBGFVSDLo(i)=DmVBGFParas_Ini((i-1)*5+4,2);
          DvVBGFVSDHi(i)=DmVBGFParas_Ini((i-1)*5+4,3);
          DivVBGFVSDPh(i)=DmVBGFParas_Ini((i-1)*5+4,4)>MINPOSITIVE ?int(DmVBGFParas_Ini((i-1)*5+4,4)+0.01) :int(DmVBGFParas_Ini((i-1)*5+4,4)-0.01) ;
          if(DmVBGFParas_Ini((i-1)*5+4,5)>MINPOSITIVE)
                    DvVBGFVSDSigma2(i)=DmVBGFParas_Ini((i-1)*5+4,5)*DmVBGFParas_Ini((i-1)*5+4,5);//log(DmVBGFParas_Ini((i-1)*5+4,5)*DmVBGFParas_Ini((i-1)*5+4,5)+1.0);
          else
                    DvVBGFVSDSigma2(i)=SIGMA2FILL;
          DvVBGFVSDSigma(i)=sqrt(DvVBGFVSDSigma2(i));
          DvVBGFVSDLambda(i)=DmVBGFParas_Ini((i-1)*5+4,6);
          DivVBGFVSDLikelihoodFlag(i)=int(DmVBGFParas_Ini((i-1)*5+4,7)+0.1);
          //LVRho
          DvVBGFLVRhoValue(i)=DmVBGFParas_Ini((i-1)*5+5,1);
          DvVBGFLVRhoLo(i)=DmVBGFParas_Ini((i-1)*5+5,2);
          DvVBGFLVRhoHi(i)=DmVBGFParas_Ini((i-1)*5+5,3);
          DivVBGFLVRhoPh(i)=DmVBGFParas_Ini((i-1)*5+5,4)>MINPOSITIVE ?int(DmVBGFParas_Ini((i-1)*5+5,4)+0.01) :int(DmVBGFParas_Ini((i-1)*5+5,4)-0.01) ;
          if(DmVBGFParas_Ini((i-1)*5+5,5)>MINPOSITIVE)
                    DvVBGFLVRhoSigma2(i)=DmVBGFParas_Ini((i-1)*5+5,5)*DmVBGFParas_Ini((i-1)*5+5,5);//log(DmVBGFParas_Ini((i-1)*5+5,5)*DmVBGFParas_Ini((i-1)*5+5,5)+1.0);
          else
                    DvVBGFLVRhoSigma2(i)=SIGMA2FILL;
          DvVBGFLVRhoSigma(i)=sqrt(DvVBGFLVRhoSigma2(i));
          DvVBGFLVRhoLambda(i)=DmVBGFParas_Ini((i-1)*5+5,6);
          DivVBGFLVRhoLikelihoodFlag(i)=int(DmVBGFParas_Ini((i-1)*5+5,7)+0.1);
          
      }
     
   // }
    #if Debug_Status
      ICHECKPIR(DvVBGFLVRhoValue);
      ICHECKPIR(DvVBGFLVRhoLo);
      ICHECKPIR(DvVBGFLVRhoHi);
      ICHECKPIR(DivVBGFLVRhoPh);
      ICHECKPIR(DvVBGFLVRhoSigma);
      ICHECKPIR(DvVBGFLVRhoLambda);
      ICHECKPIR(DivVBGFLVRhoLikelihoodFlag);
     
    #endif 
    
     //cout<<"PIRROR Data Input Over"<<endl;
  DvRPVParas_Ini.allocate(1,7,"DvRPVParas_Ini");
 ICHECKPIR(DvRPVParas_Ini);
  DvLfiftyParas_Ini.allocate(1,7,"DvLfiftyParas_Ini");
 ICHECKPIR(DvLfiftyParas_Ini);
 L1=DiCalEndYear-DiCalBeginYear+1;
  DvLfiftyLo.allocate(1,L1);
  DvLfiftyHi.allocate(1,L1);
  DivLfiftyPh.allocate(1,L1);
  DvRsexParas_Ini.allocate(1,7,"DvRsexParas_Ini");
 ICHECKPIR(DvRsexParas_Ini);
   for(i=1;i<=L1;i++)
   {   
      DvLfiftyLo(i)=DvLfiftyParas_Ini(2);
      DvLfiftyHi(i)=DvLfiftyParas_Ini(3);
      if(DvLfiftyParas_Ini(4)<0)
           DivLfiftyPh(i)=int(DvLfiftyParas_Ini(4)-0.01);
      else
           DivLfiftyPh(i)=int(DvLfiftyParas_Ini(4)+0.01);
   }
      if(DvLfiftyParas_Ini(5)<0)
           DnLfiftySigma=SIGMA2FILL;
      else
           DnLfiftySigma=DvLfiftyParas_Ini(5);
    DnLfiftySigma2=DnLfiftySigma*DnLfiftySigma;
    DnLfiftyLambda=DvLfiftyParas_Ini(6);
    DiLfiftyLikelihoodFlag=int(DvLfiftyParas_Ini(7)+0.01);
   DnRsexLo=DvRsexParas_Ini(2);
   DnRsexHi=DvRsexParas_Ini(3);
   
   if(DvRsexParas_Ini(4)<0)
       DiRsexPh=int(DvRsexParas_Ini(4)-0.01);
   else
       DiRsexPh=int(DvRsexParas_Ini(4)+0.01);
   if(DvRsexParas_Ini(5)<0)
       DnRsexSigma=SIGMA2FILL;
   else
       DnRsexSigma=DvRsexParas_Ini(5);
   DnRsexSigma2=DnRsexSigma*DnRsexSigma;
   DnRsexLambda=DvRsexParas_Ini(6);
   DiRsexLikelihoodFlag=int(DvRsexParas_Ini(7)+0.01);
  iTestValuePIR.allocate("iTestValuePIR");
 ICHECKPIR(iTestValuePIR)
    if( iTestValuePIR !=-22122)
    {
        cout<< "Error in Prior Input Data "<<endl;
        RECORD("Data Input Error, the TESTVALUE not Equal -22122 in Prior Data!");
        CloseRecord();
        exit(-8);
    }
  DvRPVLo.allocate(1,DiSizeBinNum);
  DvRPVHi.allocate(1,DiSizeBinNum);
  DivRPVPh.allocate(1,DiSizeBinNum);
   for(i=1;i<=DiSizeBinNum;i++)
   {
       DvRPVLo(i)=DvRPVParas_Ini(2);
       DvRPVHi(i)=DvRPVParas_Ini(3);
       if(DvRPVParas_Ini(4)<0)
           DivRPVPh(i)=int(DvRPVParas_Ini(4)-0.01);
       else
           DivRPVPh(i)=int(DvRPVParas_Ini(4)+0.01);
   }
   DivRPVPh(1)=-2;
   for(i=DiRecruitPrjVectNP+1;i<=DiSizeBinNum;i++)
           DivRPVPh(i)=-2;
   if(DvRPVParas_Ini(5)<0)
      DnRPVSigma=SIGMA2FILL;
   else
      DnRPVSigma=DvRPVParas_Ini(5);
    DnRPVSigma2=DnRPVSigma*DnRPVSigma;
    DnRPVLambda=DvRPVParas_Ini(6);
    DiRPVLikelihoodFlag=int(DvRPVParas_Ini(7)+0.01);
   #if Debug_Status
    
      ICHECKPIR(DvRPVLo);
      ICHECKPIR(DvRPVHi);
      ICHECKPIR(DivRPVPh);
      ICHECKPIR(DnRPVSigma);
      ICHECKPIR(DnRPVLambda);
      ICHECKPIR(DiRPVLikelihoodFlag);
      cout<<"Prior Data Input Complete"<<endl;
      // exit(-99);
   #endif
   //MnRPVLikely
   ClosePIR();
   RECORD("Prior Data Input Complete");
    if (Flagtimestep==1)
    { ad_comm::change_datafile_name("./InputFiles/Year/Projection_Year.DAT");}
    else 
    {ad_comm::change_datafile_name("./InputFiles/Season/Projection_Season.DAT");}
  DiPJDDoProjectFlag.allocate("DiPJDDoProjectFlag");
ICHECKPJD(DiPJDDoProjectFlag);
  DiPJDFinalYear.allocate("DiPJDFinalYear");
ICHECKPJD( DiPJDFinalYear)     
  DiPJDRecruitFlag.allocate("DiPJDRecruitFlag");
ICHECKPJD(DiPJDRecruitFlag);
   if( DiPJDFinalYear>200)
        DiPJDFinalYear= DiPJDFinalYear-DiBeginYear+1-DiCalEndYear;
   if(DiPJDFinalYear<0)
   {
      DiPJDFinalYear=0;
      DiPJDDoProjectFlag=0;
      WARNING<<"Project Year less than last year"<<endl;
      cout<<"Project Year less than last year"<<endl;
   }
  DvPJDRecruitData.allocate(1,DiPJDFinalYear,"DvPJDRecruitData");
if(DiPJDFinalYear>0) ICHECKPJD(DvPJDRecruitData);
  DiPJDSelectivityFlag.allocate("DiPJDSelectivityFlag");
ICHECKPJD(DiPJDSelectivityFlag);
  DivPJDRandomNumber.allocate(1,DiPJDFinalYear);
     //lHour;  lMinute;
     if(DiPJDRecruitFlag==3)
     {
         int iRandNum=int(tBeginTime);
         random_number_generator r(iRandNum);
         for(i=1;i<=DiPJDFinalYear;i++)
         {//DiCalBeginYear,DiCalEndYear
            
            lL2=int(randu(r)*(DiCalEndYear-DiCalBeginYear+1)+DiCalBeginYear);
            if(lL2>=DiCalBeginYear && lL2<=DiCalEndYear)
               DivPJDRandomNumber(i)=lL2;
            else
            {
                i--;
                continue;
            }
          }
       }
      lL2=DiSizeBinNum; 
      lL1=Flagtimestep;
      //cout<<iRandNum<<endl;
      //cout<<rand()<<endl;
      //cout<<DivPJDRandomNumber<<endl;
      //exit(-99);
  DmPJDSelectivityIni.allocate(1,lL1,1,lL2,"DmPJDSelectivityIni");
if(DiPJDSelectivityFlag) ICHECKPJD(DmPJDSelectivityIni);
  DiPJDDataFlag.allocate("DiPJDDataFlag");
ICHECKPJD(DiPJDDataFlag);
  DmPJDContrlData.allocate(1,DiPJDFinalYear,1,Flagtimestep,"DmPJDContrlData");
ICHECKPJD(DmPJDContrlData);
  DivPJDBlockGrowthMatrix.allocate(1,Flagtimestep,"DivPJDBlockGrowthMatrix");
ICHECKPJD(DivPJDBlockGrowthMatrix);
  DiTestValuePJD.allocate("DiTestValuePJD");
ICHECKPJD(DiTestValuePJD);
    if( DiTestValuePJD !=-22122)
    {
        cout<< "Error in Projection Input  Data "<<endl;
        RECORD("Data Input Error, the TESTVALUE not Equal -22122");
        CloseRecord();
        exit(-8);
    }
    ClosePJD();
    cout<<"Projection Data Input Complete"<<endl;
   if (Flagtimestep==1)
    { ad_comm::change_datafile_name("./InputFiles/Year/Parameters_Ini_Year.DAT");}
    else 
    {ad_comm::change_datafile_name("./InputFiles/Season/Parameters_Ini_Season.DAT");}
  DbnvFleetSelParams.allocate(1,DiFleetSelParaNumIni,"DbnvFleetSelParams");
ICHECKINI(DbnvFleetSelParams);
lL1=DiFleetNum*Flagtimestep;
  DbnvLogFYear1Season1.allocate(1,lL1,"DbnvLogFYear1Season1");
ICHECKINI(DbnvLogFYear1Season1);
lL1=DiFleetNum*DiYearNum*Flagtimestep; //!!lL1=(DiCalEndYear-DiCalBeginYear+1)*DiSeasonNum-1;
  DbnvLogFDevs.allocate(1,lL1,"DbnvLogFDevs");
ICHECKINI(DbnvLogFDevs);
  DbnvCPUEQBiomassE1.allocate(1,DiCPUEBlockNum,"DbnvCPUEQBiomassE1");
ICHECKINI(DbnvCPUEQBiomassE1);
  DbnvIndexSelParams.allocate(1,DiIndexSelParaNumIni,"DbnvIndexSelParams");
ICHECKINI(DbnvIndexSelParams);
  DbnvIndexQBiomassE1.allocate(1,DiIndexQBlockNum,"DbnvIndexQBiomassE1");
ICHECKINI(DbnvIndexQBiomassE1);
  DbvNYear1Para.allocate(1,DiNYear1EstParaNum,"DbvNYear1Para");
ICHECKINI(DbvNYear1Para);
lL1=2;
  DbnvLogRSPara.allocate(1,lL1,"DbnvLogRSPara");
ICHECKINI(DbnvLogRSPara);
  DbvEnvCoef.allocate(1,DiEnvNum,"DbvEnvCoef");
ICHECKINI(DbvEnvCoef);
  DbdvLogRecruitDevs.allocate(1,DiYearNum,"DbdvLogRecruitDevs");
ICHECKINI(DbdvLogRecruitDevs);
  DbnRecruitmentRh.allocate("DbnRecruitmentRh");
ICHECKINI(DbnRecruitmentRh);
iL1=1;
  DbnvRecruitLogDevsSD.allocate(1,iL1,"DbnvRecruitLogDevsSD");
ICHECKINI(DbnvRecruitLogDevsSD);
  DbnvLorenzenA.allocate(1,Flagtimestep,"DbnvLorenzenA");
ICHECKINI(DbnvLorenzenA);
  DbnvLorenzenB.allocate(1,Flagtimestep,"DbnvLorenzenB");
ICHECKINI(DbnvLorenzenB);
lL1=DiGrowthMatrixBlockNum;
  DbnvLinf.allocate(1,lL1,"DbnvLinf");
ICHECKINI(DbnvLinf);
  DbnvV.allocate(1,lL1,"DbnvV");
ICHECKINI(DbnvV);
  DbnvLSD.allocate(1,lL1,"DbnvLSD");
ICHECKINI(DbnvLSD);
  DbnvVSD.allocate(1,lL1,"DbnvVSD");
ICHECKINI(DbnvVSD);
  DbnvLVRho.allocate(1,lL1,"DbnvLVRho");
ICHECKINI(DbnvLVRho);
  DbnvRecruitPrjVect.allocate(1,DiSizeBinNum,"DbnvRecruitPrjVect");
ICHECKINI(DbnvRecruitPrjVect);
  DbnvLfifty.allocate(1,DiYearNum,"DbnvLfifty");
ICHECKINI(DbnvLfifty);
  DbnRsex.allocate("DbnRsex");
  DiTestValueINI.allocate("DiTestValueINI");
ICHECKINI(DiTestValueINI);
  DvEnvfitZscore.allocate(DiCalBeginYear,DiCalEndYear);
  
  temp1=0;
  temp2=0;
  for (k=DiCalBeginYear;k<=DiCalEndYear;k++)
  {
     Sum_temp=0;
     for(i=1;i<=DiEnvNum;i++)
     {
        Sum_temp=DmEnv(k,i)*DbvEnvCoef(i)+Sum_temp;
     }
     DvEnvfitZscore(k)=Sum_temp;
  }
  
  //cout<<DvEnvfitZscore<<endl;
  temp1=mean(DvEnvfitZscore);
  temp2=std_dev(DvEnvfitZscore);
  
  for (k=DiCalBeginYear;k<=DiCalEndYear;k++)
  {
  DvEnvfitZscore(k)=(DvEnvfitZscore(k)-temp1)/temp2;
  //DvEnvfitZscore(k)=DvEnvfitZscore(k)-temp1;
  }
  StdevEnvZ=std_dev(DvEnvfitZscore);
  //cout<<DvEnvfitZscore<<endl;
  //cout<<StdevEnvZ<<endl;
  //exit(99999);
    if(DiTestValueINI!=-22122)
    {
        cout<< "Error in Parameters Initialize  Data "<<endl;
        RECORD("Data Input Error, the TESTVALUE not Equal -22122");
        CloseRecord();
        exit(-88);
    }
    CloseINI();
  DvFleetSelLikely.allocate(1,DiFleetSelParaNumIni);
iTemp=DiFleetNum*DiSeasonNum;
  DvFYear1Likely.allocate(1,iTemp);
  DvFLogDevsLikely.allocate(1,DiFleetNum);
  DvCPUEQE1Likely.allocate(1,DiCPUEBlockNum);
  DvCatchTotalLikely.allocate(1,DiFleetNum);
  DvCatchPropLikely.allocate(1,DiFleetNum);
  DvCPUELikely.allocate(1,DiAvailCPUENum);
  DvIndexSelLikely.allocate(1,DiIndexSelParaNumIni);
  DvIndexLikely.allocate(1,DiAvailIndexNum);
  DvIndexPropLikely.allocate(1,DiAvailIndexNum);
  DvIndexQE1Likely.allocate(1,DiIndexQBlockNum);
  DvNYear1Likely.allocate(1,DiNYear1EstParaNum);
 lL1=2;
  DvRSParaLikely.allocate(1,lL1);
  DvLorenzenALikely.allocate(1,DiSeasonNum);
  DvLorenzenBLikely.allocate(1,DiSeasonNum);
  DvLinfLikely.allocate(1,DiGrowthMatrixBlockNum);
  DvVLikely.allocate(1,DiGrowthMatrixBlockNum);
  DvLSDLikely.allocate(1,DiGrowthMatrixBlockNum);
  DvVSDLikely.allocate(1,DiGrowthMatrixBlockNum);
  DvLVRhoLikely.allocate(1,DiGrowthMatrixBlockNum);
  DvTemp.allocate(1,DiSizeBinNum);
  Fd3GrowthMatrix.allocate(1,DiSeasonNum,1,DiSizeBinNum,1,DiSizeBinNum);
DiPenaltyLambdaFlag=0;
  dmFRInFleetSeason.allocate(1,DiSeasonNum,1,DiFleetNum);
  dvFRInSeason.allocate(1,DiSeasonNum);
  DvLfiftyLikely.allocate(1,DiYearNum);
cout<<"DATA INPUT COMPLETED"<<endl;
cout<<""<<endl;
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
 cout<<"PARAMETER DEFINE BEGIN"<<endl;
 RECORD("PARAMETERS SECTION!!");
 temp=0;
  PbnvFleetSelParams.allocate(1,DiFleetSelParaNumIni,DvFleetSelParaBLo,DvFleetSelParaBHi,DivFleetSelParaPh,"PbnvFleetSelParams");
  NP=0;
  for (i=1;i<=DiFleetSelParaNumIni;i++)
    {
     if (DivFleetSelParaPh(i)>MINPOSITIVE) 
        NP=NP+1;
    }
  temp=NP+temp;
  cout<<"Number of parameters for fishery seletivity: "<<NP<<endl;
lL1=DiFleetNum*Flagtimestep;
  PbnvLogFYear1Season1.allocate(1,lL1,DvLogFYear1Lo,DvLogFYear1Hi,DivLogFYear1Ph,"PbnvLogFYear1Season1");
  NP=0;
  for (i=1;i<=lL1;i++)
    {
     if (DivLogFYear1Ph(i)>MINPOSITIVE) 
        NP=NP+1;
    }
         temp=NP+temp;
lL1=DiFleetNum*(DiCalEndYear-DiCalBeginYear+1)*Flagtimestep; //!!lL1=(DiCalEndYear-DiCalBeginYear+1)*DiSeasonNum-1;
  PbnvLogFDevs.allocate(1,lL1,DvFLogDevsLo,DvFLogDevsHi,DivFLogDevsPh,"PbnvLogFDevs");
  NP1=0;
  for (i=1;i<=lL1;i++)
    {
     if (DivFLogDevsPh(i)>MINPOSITIVE) 
        NP1=NP1+1;
    }
        temp=NP1+temp;
  cout<<"Number of parameters for fishery mortality deviations: "<<NP+NP1<<endl;
  PbnvCPUEQBiomassE1.allocate(1,DiCPUEBlockNum,DvCPUEQBiomassE1Lo,DvCPUEQBiomassE1Hi,DivCPUEQBiomassE1Ph,"PbnvCPUEQBiomassE1");
  NP=0;
  for (i=1;i<=DiCPUEBlockNum;i++)
    {
     if (DivCPUEQBiomassE1Ph(i)>MINPOSITIVE) 
        NP=NP+1;
    }
          temp=NP+temp;
  cout<<"Number of parameters for catchability of fishery (E): "<<NP<<endl;
  PbnvIndexSelParams.allocate(1,DiIndexSelParaNumIni,DvIndexSelParaBLo,DvIndexSelParaBHi,DivIndexSelParaPh,"PbnvIndexSelParams");
  NP=0;
  for (i=1;i<=DiIndexSelParaNumIni;i++)
    {
     if (DivIndexSelParaPh(i)>MINPOSITIVE) 
        NP=NP+1;
    }
          temp=NP+temp;
  cout<<"Number of parameters for survey selectivity: "<<NP<<endl;
  PbnvIndexQBiomassE1.allocate(1,DiIndexQBlockNum,DvQBiomassE1Lo,DvQBiomassE1Hi,DivQBiomassE1Ph,"PbnvIndexQBiomassE1");
  NP=0;
  for (i=1;i<=DiIndexQBlockNum;i++)
    {
     if (DivQBiomassE1Ph(i)>MINPOSITIVE) 
        NP=NP+1;
    }
          temp=NP+temp;
  cout<<"Number of parameters for catchability of survey (E): "<<NP<<endl;
  PbvNYear1Para.allocate(1,DiNYear1EstParaNum,DvNYear1ParaLo,DvNYear1ParaHi,DivNYear1ParaPh,"PbvNYear1Para");
  NP=0;
  for (i=1;i<=DiNYear1EstParaNum;i++)
    {
     if (DivNYear1ParaPh(i)>MINPOSITIVE) 
        NP=NP+1;
    }
          temp=NP+temp;
  cout<<"Number of parameters for intial abundance: "<<NP<<endl;
lL1=2;
  PbnvLogRSPara.allocate(1,lL1,DvRSParaLo,DvRSParaHi,DivRSParaPh,"PbnvLogRSPara");
  PbvEnvCoef.allocate(1,DiEnvNum,DvRSEnvLo,DvRSEnvHi,DivRSEnvPh,"PbvEnvCoef");
  NP=0;
  for (i=1;i<=lL1;i++)
    {
     if (DivRSParaPh(i)>MINPOSITIVE) 
        NP=NP+1;
    }
          temp=NP+temp;
  NP=0;
  for (i=1;i<=DiEnvNum;i++)
    {
     if (DivRSEnvPh(i)>MINPOSITIVE) 
        NP=NP+1;
    }
          temp=NP+temp;
  cout<<"Number of parameters for RS: "<<NP<<endl;
  PbdvLogRecruitDevs.allocate(DiCalBeginYear,DiCalEndYear-DiYearBeforeEndForRDev,-5.0,5.0,DiRecruitLogDevsPhase,"PbdvLogRecruitDevs");
     if (DiRecruitLogDevsPhase>MINPOSITIVE) 
        {
         NP=DiCalEndYear-DiCalBeginYear+1;
        }
            temp=NP+temp;
     cout<<"Number of parameters for Rdev: "<<NP<<endl;
  PbdvLogRecruitDevsADD.allocate(DiCalEndYear-DiYearBeforeEndForRDev+1,DiCalEndYear,-5.0,5.0,DiRecruitLogDevsPhase,"PbdvLogRecruitDevsADD");
    if(DvRecruitmentRh_Ini(4) > MINPOSITIVE ) 
           lL1=int(DvRecruitmentRh_Ini(4)+0.01); 
    else
           lL1=int(DvRecruitmentRh_Ini(4)-0.01);
     nLoB=DvRecruitmentRh_Ini(2);
     nHiB=DvRecruitmentRh_Ini(3);
    if (lL1>MINPOSITIVE)
          {
          cout<<"Autocorrelation of Rdev will be estimated"<<endl;
          temp=1+temp;
          }
    else
          cout<<"Autocorrelation of Rdev will not be estimated and fixed at "<<PbnRecruitmentRh<<endl;
  PbnRecruitmentRh.allocate(nLoB,nHiB,lL1,"PbnRecruitmentRh");
      lL1=1; 
      if(DiRecruitLogDevsSDEstFlag)
      {
         if(DvRecruitLogDevsSD_Ini(4)>MINPOSITIVE )
             lL2=int(DvRecruitLogDevsSD_Ini(4)+0.01);
         else
             lL2=int(DvRecruitLogDevsSD_Ini(4)-0.01);
      }
      else
             lL2=-2;
      nLoB=DvRecruitLogDevsSD_Ini(2);
      nHiB=DvRecruitLogDevsSD_Ini(3);
  #if Debug_Status
       //ICHECKDEB(MbdvLogRecruitDevs);
       ICHECKDEB(lL1);
       ICHECKDEB(lL2);
       ICHECKDEB(nHiB);
       ICHECKDEB(nLoB);
  #endif
      if (lL2>MINPOSITIVE)
          {
          cout<<"SD of Rdev will be estimated"<<endl;
          temp=1+temp;
          }
    else
          cout<<"SD of Rdev will not be estimated and fixed at "<<DbnvRecruitLogDevsSD<<endl;
  PbnvRecruitLogDevsSD2.allocate(nLoB,nHiB,lL2,"PbnvRecruitLogDevsSD2");
    if(!DiNaturalMortalityFlag)
    {
         lL1=Flagtimestep;
    }
    else
    {    
         lL1=0;
    }
  PbnvLorenzenA.allocate(1,lL1,DvLorenzenALo,DvLorenzenAHi,DivLorenzenAPh,"PbnvLorenzenA");
  PbnvLorenzenB.allocate(1,lL1,DvLorenzenBLo,DvLorenzenBHi,DivLorenzenBPh,"PbnvLorenzenB");
  NP=0;
  NP1=0;
  for (i=1;i<=lL1;i++)
   {
     if (DivLorenzenAPh(i)>MINPOSITIVE) 
        {
         NP=NP+1;
        }
     if(PbnvLorenzenB(i)>MINPOSITIVE||DivLorenzenAPh(i)>MINPOSITIVE)
         NP1=NP1+1;
    }
  temp=NP+NP1+temp;
     cout<<"Number of parameters for natural mortality: "<<NP+NP1<<endl;
  if (NP1<MINPOSITIVE)
     cout<<"Lorenzen natural mortality is not used"<<endl;
  else
     cout<<"Lorenzen natural mortality is used"<<endl;
      if(DiGrowthMatrixFlag)
      {
            lL1=  DiGrowthMatrixBlockNum;
      }
      else
      {     
            lL1=DiGrowthMatrixBlockNum;
            for(i=1;i<=DiGrowthMatrixBlockNum;i++)
            {
                  DivVBGFLinfPh(i) =-2;
                  DivVBGFVPh(i)    =-2;
                  DivVBGFLSDPh(i)  =-2;
                  DivVBGFVSDPh(i)  =-2;
                  DivVBGFLVRhoPh(i)=-2;
            }
      }
  PbnvLinf.allocate(1,lL1,DvVBGFLinfLo,DvVBGFLinfHi,DivVBGFLinfPh,"PbnvLinf");
  NP=0;
  temp1=0;
  for (i=1;i<=lL1;i++)
   {
     if (DivVBGFLinfPh(i)>MINPOSITIVE) 
        {
         NP=NP+1;
        }
    }        
  temp=temp+NP;
  temp1=temp1+NP1;
  PbnvV.allocate(1,lL1,DvVBGFVLo,DvVBGFVHi,DivVBGFVPh,"PbnvV");
  NP=0;
  temp1=0;
  for (i=1;i<=lL1;i++)
   {
     if (DivVBGFVPh(i)>MINPOSITIVE) 
        {
         NP=NP+1;
        }
    }        
  temp=temp+NP;
  temp1=temp1+NP1;
  PbnvLSD.allocate(1,lL1,DvVBGFLSDLo,DvVBGFLSDHi,DivVBGFLSDPh,"PbnvLSD");
  NP=0;
  temp1=0;
  for (i=1;i<=lL1;i++)
   {
     if (DivVBGFLSDPh(i)>MINPOSITIVE) 
        {
         NP=NP+1;
        }
    }        
  temp=temp+NP;
  temp1=temp1+NP1;
  PbnvVSD.allocate(1,lL1,DvVBGFVSDLo,DvVBGFVSDHi,DivVBGFVSDPh,"PbnvVSD");
  NP=0;
  temp1=0;
  for (i=1;i<=lL1;i++)
   {
     if (DivVBGFVSDPh(i)>MINPOSITIVE) 
        {
         NP=NP+1;
        }
    }        
  temp=temp+NP;
  temp1=temp1+NP1;
  PbnvLVRho.allocate(1,lL1,DvVBGFLVRhoLo,DvVBGFLVRhoHi,DivVBGFLVRhoPh,"PbnvLVRho");
  NP=0;
  temp1=0;
  for (i=1;i<=lL1;i++)
   {
     if (DivVBGFLVRhoPh(i)>MINPOSITIVE) 
        {
         NP=NP+1;
        }
    }        
  temp=temp+NP;
  temp1=temp1+NP1;
  cout<<"Number of parameters for Growth matrix: "<<temp1<<endl;
 L1=DiCalEndYear-DiCalBeginYear+1;
  PbLfifty.allocate(1,L1,DvLfiftyLo,DvLfiftyHi,DivLfiftyPh,"PbLfifty");
  NP=0;
  for (i=1;i<=L1;i++)
   {
     if (DivLfiftyPh(i)>MINPOSITIVE) 
        {
         NP=NP+1;
        }
    }        
  temp=temp+NP;
  cout<<"Number of parameters for L50 of sex change: "<<NP<<endl;
  PbRsex.allocate(DnRsexLo,DnRsexHi,DiRsexPh,"PbRsex");
  NP=0;
    if (DiRsexPh>MINPOSITIVE) 
       NP=NP+1;        
  temp=temp+NP;
  cout<<"Number of parameters for Rsex: "<<NP<<endl;
  PbnvRecruitPrjVect.allocate(1,DiSizeBinNum,DvRPVLo,DvRPVHi,DivRPVPh,"PbnvRecruitPrjVect");
  NP=0;
  for (i=1;i<=DiSizeBinNum;i++)
   {
     if (DivRPVPh(i)>MINPOSITIVE) 
        {
         NP=NP+1;
        }
    }        
  temp=temp+NP;
  cout<<"Number of parameters for recruitment proportion-at-size: "<<NP<<endl;
  cout<<"TOTAL NUMBER OF PARAMETERS TO BE ESTIMATED: "<<temp<<endl;
  cout<<""<<endl;
  Md3M.allocate(1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md3M");
  #ifndef NO_AD_INITIALIZE
    Md3M.initialize();
  #endif
  Md4SelByFleet.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md4SelByFleet");
  #ifndef NO_AD_INITIALIZE
    Md4SelByFleet.initialize();
  #endif
  Md3Ft.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear,"Md3Ft");
  #ifndef NO_AD_INITIALIZE
    Md3Ft.initialize();
  #endif
  Md4FASByFleet.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md4FASByFleet");
  #ifndef NO_AD_INITIALIZE
    Md4FASByFleet.initialize();
  #endif
  Md3Ftot.allocate(1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md3Ftot");
  #ifndef NO_AD_INITIALIZE
    Md3Ftot.initialize();
  #endif
  Md3Z.allocate(1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md3Z");
  #ifndef NO_AD_INITIALIZE
    Md3Z.initialize();
  #endif
  Md3S.allocate(1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md3S");
  #ifndef NO_AD_INITIALIZE
    Md3S.initialize();
  #endif
  MmSSB_S.allocate(DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"MmSSB_S");
  #ifndef NO_AD_INITIALIZE
    MmSSB_S.initialize();
  #endif
  MmSSB_Y.allocate(DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"MmSSB_Y");
  #ifndef NO_AD_INITIALIZE
    MmSSB_Y.initialize();
  #endif
  MvSSB_S1.allocate(1,DiSizeBinNum,"MvSSB_S1");
  #ifndef NO_AD_INITIALIZE
    MvSSB_S1.initialize();
  #endif
  MvRecruitRS.allocate(DiCalBeginYear,DiCalEndYear+1,"MvRecruitRS");
  #ifndef NO_AD_INITIALIZE
    MvRecruitRS.initialize();
  #endif
  MvRecruitDev.allocate(DiCalBeginYear,DiCalEndYear,"MvRecruitDev");
  #ifndef NO_AD_INITIALIZE
    MvRecruitDev.initialize();
  #endif
  MnRSD.allocate("MnRSD");
  #ifndef NO_AD_INITIALIZE
  MnRSD.initialize();
  #endif
  MbdvLogRecruitDevs.allocate(DiCalBeginYear,DiCalEndYear,"MbdvLogRecruitDevs");
  #ifndef NO_AD_INITIALIZE
    MbdvLogRecruitDevs.initialize();
  #endif
  MvRecruitPrjVect.allocate(1,DiSizeBinNum,"MvRecruitPrjVect");
  #ifndef NO_AD_INITIALIZE
    MvRecruitPrjVect.initialize();
  #endif
  MbnvR00Para.allocate(1,1,"MbnvR00Para");
  #ifndef NO_AD_INITIALIZE
    MbnvR00Para.initialize();
  #endif
  MnZscore1.allocate("MnZscore1");
  #ifndef NO_AD_INITIALIZE
  MnZscore1.initialize();
  #endif
  MnZscore2.allocate("MnZscore2");
  #ifndef NO_AD_INITIALIZE
  MnZscore2.initialize();
  #endif
  MnDeltaL.allocate("MnDeltaL");
  #ifndef NO_AD_INITIALIZE
  MnDeltaL.initialize();
  #endif
  MnDeltaLVar.allocate("MnDeltaLVar");
  #ifndef NO_AD_INITIALIZE
  MnDeltaLVar.initialize();
  #endif
   if(DiGrowthMatrixFlag==2)
   {
         iL1=DiGrowthMatrixBlockNum*DiSeasonNum;
   }
   else
   {
         iL1=DiGrowthMatrixBlockNum;     
   }
  Md3GrowthMatrix.allocate(1,iL1,1,DiSizeBinNum,1,DiSizeBinNum,"Md3GrowthMatrix");
  #ifndef NO_AD_INITIALIZE
    Md3GrowthMatrix.initialize();
  #endif
  Md3NAS.allocate(1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md3NAS");
  #ifndef NO_AD_INITIALIZE
    Md3NAS.initialize();
  #endif
  MvNASPJD.allocate(1,DiSizeBinNum,"MvNASPJD");
  #ifndef NO_AD_INITIALIZE
    MvNASPJD.initialize();
  #endif
  MvNASPJD_1.allocate(1,DiSizeBinNum-DiRecruitPrjVectNP,"MvNASPJD_1");
  Md3CatchTotalP.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear,"Md3CatchTotalP");
  #ifndef NO_AD_INITIALIZE
    Md3CatchTotalP.initialize();
  #endif
  Md4CatchAtSizeP.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md4CatchAtSizeP");
  #ifndef NO_AD_INITIALIZE
    Md4CatchAtSizeP.initialize();
  #endif
  Md4CatchPropAtSizeP.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md4CatchPropAtSizeP");
  #ifndef NO_AD_INITIALIZE
    Md4CatchPropAtSizeP.initialize();
  #endif
  MvCPUEQ.allocate(1,DiCPUEBlockNum,"MvCPUEQ");
  #ifndef NO_AD_INITIALIZE
    MvCPUEQ.initialize();
  #endif
  MmCPUEPred.allocate(1,DiAvailCPUENum,DiCalBeginYear,DiCalEndYear,"MmCPUEPred");
  #ifndef NO_AD_INITIALIZE
    MmCPUEPred.initialize();
  #endif
  Md3NASBar.allocate(1,Flagtimestep,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md3NASBar");
  #ifndef NO_AD_INITIALIZE
    Md3NASBar.initialize();
  #endif
  MmBiomassExploit.allocate(1,DiAvailCPUENum,DiCalBeginYear,DiCalEndYear,"MmBiomassExploit");
  #ifndef NO_AD_INITIALIZE
    MmBiomassExploit.initialize();
  #endif
  Md3IndexSel.allocate(1,DiAvailIndexNum,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md3IndexSel");
  #ifndef NO_AD_INITIALIZE
    Md3IndexSel.initialize();
  #endif
  MmBiomassSurvey.allocate(1,DiAvailIndexNum,DiCalBeginYear,DiCalEndYear,"MmBiomassSurvey");
  #ifndef NO_AD_INITIALIZE
    MmBiomassSurvey.initialize();
  #endif
  Md3NASIndexBar.allocate(1,DiAvailIndexNum,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md3NASIndexBar");
  #ifndef NO_AD_INITIALIZE
    Md3NASIndexBar.initialize();
  #endif
  MvIndexQ.allocate(1,DiIndexQBlockNum,"MvIndexQ");
  #ifndef NO_AD_INITIALIZE
    MvIndexQ.initialize();
  #endif
  MmIndexP.allocate(1,DiAvailIndexNum,DiCalBeginYear,DiCalEndYear,"MmIndexP");
  #ifndef NO_AD_INITIALIZE
    MmIndexP.initialize();
  #endif
  Md3IndexPropAtSizeP.allocate(1,DiAvailIndexNum,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md3IndexPropAtSizeP");
  #ifndef NO_AD_INITIALIZE
    Md3IndexPropAtSizeP.initialize();
  #endif
  Md3FemalePropAtSizeP.allocate(1,1,DiCalBeginYear,DiCalEndYear,1,DiSizeBinNum,"Md3FemalePropAtSizeP");
  #ifndef NO_AD_INITIALIZE
    Md3FemalePropAtSizeP.initialize();
  #endif
  sdr_vSSB.allocate(DiCalBeginYear,DiCalEndYear,"sdr_vSSB");
  #ifndef NO_AD_INITIALIZE
    sdr_vSSB.initialize();
  #endif
  sdr_vSSB2.allocate(DiCalBeginYear,DiCalEndYear,"sdr_vSSB2");
  #ifndef NO_AD_INITIALIZE
    sdr_vSSB2.initialize();
  #endif
  sdr_nBPRMSY.allocate("sdr_nBPRMSY");
  #ifndef NO_AD_INITIALIZE
  sdr_nBPRMSY.initialize();
  #endif
  sdr_vNBiomass.allocate(DiCalBeginYear,DiCalEndYear,"sdr_vNBiomass");
  #ifndef NO_AD_INITIALIZE
    sdr_vNBiomass.initialize();
  #endif
  sdr_vRecruitmentP.allocate(DiCalBeginYear,DiCalEndYear,"sdr_vRecruitmentP");
  #ifndef NO_AD_INITIALIZE
    sdr_vRecruitmentP.initialize();
  #endif
  MnBPRF30SPR.allocate("MnBPRF30SPR");
  #ifndef NO_AD_INITIALIZE
  MnBPRF30SPR.initialize();
  #endif
  MnBPRF40SPR.allocate("MnBPRF40SPR");
  #ifndef NO_AD_INITIALIZE
  MnBPRF40SPR.initialize();
  #endif
  MnBPRFmsy.allocate("MnBPRFmsy");
  #ifndef NO_AD_INITIALIZE
  MnBPRFmsy.initialize();
  #endif
  MnBPRF01.allocate("MnBPRF01");
  #ifndef NO_AD_INITIALIZE
  MnBPRF01.initialize();
  #endif
  MnBPRFMax.allocate("MnBPRFMax");
  #ifndef NO_AD_INITIALIZE
  MnBPRFMax.initialize();
  #endif
  MmBPRSel.allocate(1,Flagtimestep,1,DiSizeBinNum,"MmBPRSel");
  #ifndef NO_AD_INITIALIZE
    MmBPRSel.initialize();
  #endif
  MnBPRYield.allocate("MnBPRYield");
  #ifndef NO_AD_INITIALIZE
  MnBPRYield.initialize();
  #endif
  MnBPRSSB.allocate("MnBPRSSB");
  #ifndef NO_AD_INITIALIZE
  MnBPRSSB.initialize();
  #endif
  MnBPRBiomassMSY.allocate("MnBPRBiomassMSY");
  #ifndef NO_AD_INITIALIZE
  MnBPRBiomassMSY.initialize();
  #endif
  MnBPRSSBMSY.allocate("MnBPRSSBMSY");
  #ifndef NO_AD_INITIALIZE
  MnBPRSSBMSY.initialize();
  #endif
  MnBPRYieldTotal.allocate("MnBPRYieldTotal");
  #ifndef NO_AD_INITIALIZE
  MnBPRYieldTotal.initialize();
  #endif
  MnBPRSSBTotal.allocate("MnBPRSSBTotal");
  #ifndef NO_AD_INITIALIZE
  MnBPRSSBTotal.initialize();
  #endif
  MvBPRNAS.allocate(1,DiSizeBinNum,"MvBPRNAS");
  #ifndef NO_AD_INITIALIZE
    MvBPRNAS.initialize();
  #endif
  MvBPRSSB_S.allocate(1,DiSizeBinNum,"MvBPRSSB_S");
  #ifndef NO_AD_INITIALIZE
    MvBPRSSB_S.initialize();
  #endif
  MmBPRF.allocate(1,Flagtimestep,1,DiSizeBinNum,"MmBPRF");
  #ifndef NO_AD_INITIALIZE
    MmBPRF.initialize();
  #endif
  MmBPRZ.allocate(1,Flagtimestep,1,DiSizeBinNum,"MmBPRZ");
  #ifndef NO_AD_INITIALIZE
    MmBPRZ.initialize();
  #endif
  MmBPRS.allocate(1,Flagtimestep,1,DiSizeBinNum,"MmBPRS");
  #ifndef NO_AD_INITIALIZE
    MmBPRS.initialize();
  #endif
  Md3FleetEffectiveSampleSize.allocate(1,DiFleetNum,1,Flagtimestep,DiCalBeginYear,DiCalEndYear,"Md3FleetEffectiveSampleSize");
  #ifndef NO_AD_INITIALIZE
    Md3FleetEffectiveSampleSize.initialize();
  #endif
  MmIndexEffectiveSampleSize.allocate(1,DiAvailIndexNum,DiCalBeginYear,DiCalEndYear,"MmIndexEffectiveSampleSize");
  #ifndef NO_AD_INITIALIZE
    MmIndexEffectiveSampleSize.initialize();
  #endif
  MvNASTemp.allocate(1,DiSizeBinNum,"MvNASTemp");
  #ifndef NO_AD_INITIALIZE
    MvNASTemp.initialize();
  #endif
  MnTemp1.allocate("MnTemp1");
  #ifndef NO_AD_INITIALIZE
  MnTemp1.initialize();
  #endif
  MnTemp2.allocate("MnTemp2");
  #ifndef NO_AD_INITIALIZE
  MnTemp2.initialize();
  #endif
  MnTemp3.allocate("MnTemp3");
  #ifndef NO_AD_INITIALIZE
  MnTemp3.initialize();
  #endif
  MvTemp1.allocate(1,DiSizeBinNum,"MvTemp1");
  #ifndef NO_AD_INITIALIZE
    MvTemp1.initialize();
  #endif
  MvPJDF.allocate(1,DiSizeBinNum,"MvPJDF");
  #ifndef NO_AD_INITIALIZE
    MvPJDF.initialize();
  #endif
  MvPJDZ.allocate(1,DiSizeBinNum,"MvPJDZ");
  #ifndef NO_AD_INITIALIZE
    MvPJDZ.initialize();
  #endif
  MvPJDS.allocate(1,DiSizeBinNum,"MvPJDS");
  #ifndef NO_AD_INITIALIZE
    MvPJDS.initialize();
  #endif
  MvPJDSSB_S.allocate(1,DiSizeBinNum,"MvPJDSSB_S");
  #ifndef NO_AD_INITIALIZE
    MvPJDSSB_S.initialize();
  #endif
  Md3PJDNAS.allocate(1,Flagtimestep,1,DiPJDFinalYear,1,DiSizeBinNum,"Md3PJDNAS");
  #ifndef NO_AD_INITIALIZE
    Md3PJDNAS.initialize();
  #endif
  MmPJDFtot.allocate(1,Flagtimestep,1,DiPJDFinalYear,"MmPJDFtot");
  #ifndef NO_AD_INITIALIZE
    MmPJDFtot.initialize();
  #endif
  MvPJDSSB.allocate(1,DiPJDFinalYear,"MvPJDSSB");
  #ifndef NO_AD_INITIALIZE
    MvPJDSSB.initialize();
  #endif
  MmPJDTotalCatch.allocate(1,Flagtimestep,1,DiPJDFinalYear,"MmPJDTotalCatch");
  #ifndef NO_AD_INITIALIZE
    MmPJDTotalCatch.initialize();
  #endif
  MmPJDSel.allocate(1,Flagtimestep,1,DiSizeBinNum,"MmPJDSel");
  #ifndef NO_AD_INITIALIZE
    MmPJDSel.initialize();
  #endif
  MvFleetSelLikely.allocate(1,DiFleetSelParaNumIni,"MvFleetSelLikely");
  #ifndef NO_AD_INITIALIZE
    MvFleetSelLikely.initialize();
  #endif
  MnFleetSelLikely.allocate("MnFleetSelLikely");
  #ifndef NO_AD_INITIALIZE
  MnFleetSelLikely.initialize();
  #endif
 iTemp=DiFleetNum*DiSeasonNum;
  MvFYear1Likely.allocate(1,iTemp,"MvFYear1Likely");
  #ifndef NO_AD_INITIALIZE
    MvFYear1Likely.initialize();
  #endif
  MnFYear1Likely.allocate("MnFYear1Likely");
  #ifndef NO_AD_INITIALIZE
  MnFYear1Likely.initialize();
  #endif
  MnFLogDevsLikely.allocate("MnFLogDevsLikely");
  #ifndef NO_AD_INITIALIZE
  MnFLogDevsLikely.initialize();
  #endif
  MvCPUEQE1Likely.allocate(1,DiCPUEBlockNum,"MvCPUEQE1Likely");
  #ifndef NO_AD_INITIALIZE
    MvCPUEQE1Likely.initialize();
  #endif
  MnCPUEQE1Likely.allocate("MnCPUEQE1Likely");
  #ifndef NO_AD_INITIALIZE
  MnCPUEQE1Likely.initialize();
  #endif
  MvCatchTotalLikely.allocate(1,DiFleetNum,"MvCatchTotalLikely");
  #ifndef NO_AD_INITIALIZE
    MvCatchTotalLikely.initialize();
  #endif
  MvCatchPropLikely.allocate(1,DiFleetNum,"MvCatchPropLikely");
  #ifndef NO_AD_INITIALIZE
    MvCatchPropLikely.initialize();
  #endif
  MnCatchLikely.allocate("MnCatchLikely");
  #ifndef NO_AD_INITIALIZE
  MnCatchLikely.initialize();
  #endif
  MvCPUELikely.allocate(1,DiAvailCPUENum,"MvCPUELikely");
  #ifndef NO_AD_INITIALIZE
    MvCPUELikely.initialize();
  #endif
  MnCPUELikely.allocate("MnCPUELikely");
  #ifndef NO_AD_INITIALIZE
  MnCPUELikely.initialize();
  #endif
  MvIndexSelLikely.allocate(1,DiIndexSelParaNumIni,"MvIndexSelLikely");
  #ifndef NO_AD_INITIALIZE
    MvIndexSelLikely.initialize();
  #endif
  MnIndexSelLikely.allocate("MnIndexSelLikely");
  #ifndef NO_AD_INITIALIZE
  MnIndexSelLikely.initialize();
  #endif
  MvIndexLikely.allocate(1,DiAvailIndexNum,"MvIndexLikely");
  #ifndef NO_AD_INITIALIZE
    MvIndexLikely.initialize();
  #endif
  MvIndexPropLikely.allocate(1,DiAvailIndexNum,"MvIndexPropLikely");
  #ifndef NO_AD_INITIALIZE
    MvIndexPropLikely.initialize();
  #endif
  MnIndexLikely.allocate("MnIndexLikely");
  #ifndef NO_AD_INITIALIZE
  MnIndexLikely.initialize();
  #endif
  MvIndexQE1Likely.allocate(1,DiIndexQBlockNum,"MvIndexQE1Likely");
  #ifndef NO_AD_INITIALIZE
    MvIndexQE1Likely.initialize();
  #endif
  MnIndexQE1Likely.allocate("MnIndexQE1Likely");
  #ifndef NO_AD_INITIALIZE
  MnIndexQE1Likely.initialize();
  #endif
  MvNYear1Likely.allocate(1,DiNYear1EstParaNum,"MvNYear1Likely");
  #ifndef NO_AD_INITIALIZE
    MvNYear1Likely.initialize();
  #endif
  MnNYear1Likely.allocate("MnNYear1Likely");
  #ifndef NO_AD_INITIALIZE
  MnNYear1Likely.initialize();
  #endif
 lL1=2;
  MvRSParaLikely.allocate(1,lL1,"MvRSParaLikely");
  #ifndef NO_AD_INITIALIZE
    MvRSParaLikely.initialize();
  #endif
  MnRSParaLikely.allocate("MnRSParaLikely");
  #ifndef NO_AD_INITIALIZE
  MnRSParaLikely.initialize();
  #endif
  MnRDevsSDLikely.allocate("MnRDevsSDLikely");
  #ifndef NO_AD_INITIALIZE
  MnRDevsSDLikely.initialize();
  #endif
  MnRecruitDevsLikely.allocate("MnRecruitDevsLikely");
  #ifndef NO_AD_INITIALIZE
  MnRecruitDevsLikely.initialize();
  #endif
  MnRecruitRhLikely.allocate("MnRecruitRhLikely");
  #ifndef NO_AD_INITIALIZE
  MnRecruitRhLikely.initialize();
  #endif
  MnEnvfitRdevLikely.allocate("MnEnvfitRdevLikely");
  #ifndef NO_AD_INITIALIZE
  MnEnvfitRdevLikely.initialize();
  #endif
  MvLorenzenALikely.allocate(1,Flagtimestep,"MvLorenzenALikely");
  #ifndef NO_AD_INITIALIZE
    MvLorenzenALikely.initialize();
  #endif
  MnLorenzenALikely.allocate("MnLorenzenALikely");
  #ifndef NO_AD_INITIALIZE
  MnLorenzenALikely.initialize();
  #endif
  MvLorenzenBLikely.allocate(1,Flagtimestep,"MvLorenzenBLikely");
  #ifndef NO_AD_INITIALIZE
    MvLorenzenBLikely.initialize();
  #endif
  MnLorenzenBLikely.allocate("MnLorenzenBLikely");
  #ifndef NO_AD_INITIALIZE
  MnLorenzenBLikely.initialize();
  #endif
  MvLinfLikely.allocate(1,DiGrowthMatrixBlockNum,"MvLinfLikely");
  #ifndef NO_AD_INITIALIZE
    MvLinfLikely.initialize();
  #endif
  MvVLikely.allocate(1,DiGrowthMatrixBlockNum,"MvVLikely");
  #ifndef NO_AD_INITIALIZE
    MvVLikely.initialize();
  #endif
  MvLSDLikely.allocate(1,DiGrowthMatrixBlockNum,"MvLSDLikely");
  #ifndef NO_AD_INITIALIZE
    MvLSDLikely.initialize();
  #endif
  MvVSDLikely.allocate(1,DiGrowthMatrixBlockNum,"MvVSDLikely");
  #ifndef NO_AD_INITIALIZE
    MvVSDLikely.initialize();
  #endif
  MvLVRhoLikely.allocate(1,DiGrowthMatrixBlockNum,"MvLVRhoLikely");
  #ifndef NO_AD_INITIALIZE
    MvLVRhoLikely.initialize();
  #endif
  MnGrowthMatrixLikely.allocate("MnGrowthMatrixLikely");
  #ifndef NO_AD_INITIALIZE
  MnGrowthMatrixLikely.initialize();
  #endif
  MvLfiftylikely.allocate(1,DiYearNum,"MvLfiftylikely");
  #ifndef NO_AD_INITIALIZE
    MvLfiftylikely.initialize();
  #endif
  MnLfiftylikely.allocate("MnLfiftylikely");
  #ifndef NO_AD_INITIALIZE
  MnLfiftylikely.initialize();
  #endif
  MnRsexlikely.allocate("MnRsexlikely");
  #ifndef NO_AD_INITIALIZE
  MnRsexlikely.initialize();
  #endif
  MnFMaxPenaltyLikely.allocate("MnFMaxPenaltyLikely");
  #ifndef NO_AD_INITIALIZE
  MnFMaxPenaltyLikely.initialize();
  #endif
  MnRPVLikely.allocate("MnRPVLikely");
  #ifndef NO_AD_INITIALIZE
  MnRPVLikely.initialize();
  #endif
  ofvTotal.allocate("ofvTotal");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  //Initialize the paramaters
  RECORD("PRELIMINARY_CALCS_SECTION");
  cout<<"Parameters initialized Begin"<<endl;
  ////////////////////////////////////////////////////////////////////////////////////
  for(i=1;i<=DiFleetSelParaNumIni;i++)
      PbnvFleetSelParams(i)=DbnvFleetSelParams(i);
  lL1=DiFleetNum*Flagtimestep;
  for(i=1;i<=lL1;i++)
      PbnvLogFYear1Season1(i)=DbnvLogFYear1Season1(i);
 
  for(i=1;i<=DiFleetNum;i++)
  {
      for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
      {
        for(k=1;k<=Flagtimestep;k++)
        {
            iL1=(i-1)*DiYearNum*Flagtimestep+(j-1)*Flagtimestep+k;
            iL2=(i-1)*(DiCalEndYear-DiCalBeginYear+1)*Flagtimestep+(j-DiCalBeginYear)*Flagtimestep+k;;
            PbnvLogFDevs(iL2)=DbnvLogFDevs(iL1);
        }
      }
  }
  
  for(i=1;i<=DiCPUEBlockNum;i++)
       PbnvCPUEQBiomassE1(i)=DbnvCPUEQBiomassE1(i);
  for(i=1;i<=DiIndexSelParaNumIni;i++)
      PbnvIndexSelParams(i)=DbnvIndexSelParams(i);
  for(i=1;i<=DiIndexQBlockNum;i++)
      PbnvIndexQBiomassE1(i)=DbnvIndexQBiomassE1(i);
 
  lL1=DiNYear1EstParaNum;
  for(i=1;i<=lL1;i++)
      PbvNYear1Para(i)=DbvNYear1Para(i);
  lL1=2;
  for(i=1;i<=lL1;i++)
       PbnvLogRSPara(i)= DbnvLogRSPara(i);
  lL1=DiEnvNum;
  for(i=1;i<=lL1;i++)
       PbvEnvCoef(i)= DbvEnvCoef(i);
  
  //DiCalEndYear-DiYearBeforeEndForRDev
  for(i=DiCalBeginYear;i<=DiCalEndYear-DiYearBeforeEndForRDev;i++)
      PbdvLogRecruitDevs(i)=DbdvLogRecruitDevs(i);
   for(i=DiCalEndYear-DiYearBeforeEndForRDev+1;i<=DiCalEndYear;i++)
      PbdvLogRecruitDevsADD(i)=DbdvLogRecruitDevs(i);
   PbnRecruitmentRh=DbnRecruitmentRh;
   PbnvRecruitLogDevsSD2=DbnvRecruitLogDevsSD(1);
   MnRSD=PbnvRecruitLogDevsSD2;
   if(!DiRecruitLogDevsSDEstFlag)
   {
         MnRSD=DnRecruitLogDevsSigma;
   
   }
  
  lL1=Flagtimestep;
  for(i=1;i<=lL1;i++)
  {
       PbnvLorenzenA(i)=DbnvLorenzenA(i);
       PbnvLorenzenB(i)=DbnvLorenzenB(i);
  }
 
  lL1=DiGrowthMatrixBlockNum;
  for(i=1;i<=lL1;i++)
  {
      PbnvLinf(i)=DbnvLinf(i);
      PbnvV(i)=DbnvV(i);
      PbnvLSD(i)=DbnvLSD(i);
      PbnvVSD(i)=DbnvVSD(i);
      PbnvLVRho(i)=DbnvLVRho(i);
  }
  lL1=DiSizeBinNum;
  for(i=1;i<=lL1;i++)
  {
       PbnvRecruitPrjVect(i)=DbnvRecruitPrjVect(i);
  }
  if(PbnvRecruitPrjVect(1)<=MINPOSITIVE)
         PbnvRecruitPrjVect(1)=0.8;
  
  lL1=DiCalEndYear-DiCalBeginYear+1;
  for(i=1;i<=lL1;i++)
  {
       PbLfifty(i)=DbnvLfifty(i);
  }
  PbRsex=DbnRsex;
  cout<<"Parameters initialized End"<<endl;
  //////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
 #if Debug_Status      
      ICHECKDEB(PbnvFleetSelParams);
      ICHECKDEB(PbnvLogFYear1Season1);
      ICHECKDEB(PbnvLogFDevs);
      ICHECKDEB(PbnvCPUEQBiomassE1);
      ICHECKDEB(PbnvIndexSelParams);
      ICHECKDEB(PbnvIndexQBiomassE1);
      ICHECKDEB(PbvNYear1Para);
      ICHECKDEB(PbnvLogRSPara);
      ICHECKDEB(PbdvLogRecruitDevs);
      ICHECKDEB(PbnRecruitmentRh);
      ICHECKDEB(PbnvRecruitLogDevsSD2);
      ICHECKDEB(PbnvLorenzenA);
      ICHECKDEB(PbnvLVRho);
      ICHECKDEB(PbnvRecruitPrjVect);
      //exit(-99);
  #endif 
  if(DiGrowthMatrixFlag)
  {
       //Md3GrowthMatrix.initialize();
       //Md3GrowthMatrix(2)=1;
       //Md3GrowthMatrix(1)=1;
      //cout<<"Growth Matrix Calculation Begin"<<endl;
      //cout<<"DiGrowthMatrixBlockNum"<<DiGrowthMatrixBlockNum<<endl;
       for(i=1;i<=DiGrowthMatrixBlockNum;i++)
       {
          //cout<<i<<endl;   
          // Md3GrowthMatrix(i).initialize();
           //cout<<Md3GrowthMatrix<<endl;
           if(DiGrowthMatrixFlag==2)
           {
                for(j=1;j<=Flagtimestep;j++)
                {
                      iL1=(i-1)*Flagtimestep+j;
                      Md3GrowthMatrix(iL1).initialize();
                      CalGrowthMatrix(PbnvLinf(i),PbnvV(i), PbnvLSD(i),PbnvVSD(i),PbnvLVRho(i),DvGrowthTimeAsYear(j),20,0,iL1);
                }
           }
           else
           {
                Md3GrowthMatrix(i).initialize();
                CalGrowthMatrix(PbnvLinf(i),PbnvV(i), PbnvLSD(i),PbnvVSD(i),PbnvLVRho(i),1.0,20,0,i);
           }
          //cout<<i<<endl; 
       }
       cout<<"Growth Matrix Calculation End"<<endl;
 #if Debug_Status
    ICHECKDEB(Md3GrowthMatrix);
    //exit(-99);
 #endif
  }
  else
  {
         Md3GrowthMatrix.initialize();
         for(i=1;i<=DiGrowthMatrixBlockNum;i++)
         {
             for(j=1;j<=DiSizeBinNum;j++)
             {
                 for(k=1;k<=DiSizeBinNum;k++)
                 {
                     Md3GrowthMatrix(i,j,k)= Dd3GrowthMatrix_Ini(i,j,k);
                 }
             }
          }
  }
   Md3M.initialize();
   if(DiNaturalMortalityFlag)
   {
       for(i=1;i<=Flagtimestep;i++)
       {
         for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
         {
               for(k=1;k<=DiSizeBinNum;k++)
               {
                     iL1=(i-1)*DiYearNum+j;
                     Md3M(i,j,k)=DmNaturalMAtSize_Ini(iL1,k+2);
               }
          }
       }      
   }
   else
   {//if(!DiNaturalMortalityFlag)
      for(i=1;i<=Flagtimestep;i++)
      {
         for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
         {
               for(k=1;k<=DiSizeBinNum;k++)
               {//DvNaturalMSizeWeight(k)*DvNaturalMYearWeight(j);
                     Md3M(i,j,k)=PbnvLorenzenA(i)*pow(DmWeightAtSize(j,k),PbnvLorenzenB(i))*DvNaturalMSizeWeight(k)*DvNaturalMYearWeight(j);
                     //DmNaturalMSizeWeight(i,k)*DmNaturalMYearWeight(i,j);
               }
          }
      }  
   }
  
  
  ///////////////////////////////////////////////////////////////////////////////////////
  // Fishery Selectivity and Fishery Mortality
  if(DiFleetSelParaNumIni>0)
      DvFleetSelLikely.initialize();   //   vector                     DvFleetSelLikely(1,DiFleetSelParaNum)
  DnFleetSelLikely=0.0;           //number                     DnFleetSelLikely
  // F(f,1,1)
  DvFYear1Likely.initialize();  //vector                     DvFYear1Likely(1,iFleetNum)
  DnFYear1Likely=0.0;          // number                     DnFYear1Likely
  // FDev(f)
  DvFLogDevsLikely.initialize(); //vector                     DvFLogDevsLikely(1,iFleetNum)  
  DnFLogDevsLikely=0.0;         //number                     DnFLogDevsLikely
  // CPUE ^E1
  if(DiCPUEBlockNum>0)
        DvCPUEQE1Likely.initialize(); // vector                     DvCPUEQE1Likely(1,DiAvailCPUEQNumP)
  DnCPUEQE1Likely=0.0;     // number                     DnCPUEQE1Likely
  
  //Catch at Size
  DvCatchTotalLikely.initialize();  // vector                     DvCatchTotalLikely(1,DiFleetNum)
  DvCatchPropLikely.initialize();  // vector                     DvCatchPropLikely(1,DiFleetNum)
  DnCatchLikely=0.0;             // number                     DnCatchLikely                                //sum
  //CPUE
  if(DiAvailCPUENum>0)
    DvCPUELikely.initialize();//vector                     DvCPUELikely(1,DiAvailCPUENum)
   DnCPUELikely=0.0;// number                     DnCPUELikely
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Survey Index
  // Survey Index Selectivity
  if(DiIndexSelParaNumIni>0)
      DvIndexSelLikely.initialize(); //vector                     DvIndexSelLikely(1,DiIndexSelParaNum)
  DnIndexSelLikely=0.0;         //number                     DnIndexSelLikely
  
  //Index at Size
  if(DiAvailIndexNum>0)
  {
      DvIndexLikely.initialize();//vector                     DvIndexLikely(1,DiAvailIndexNum)
      DvIndexPropLikely.initialize();//vector                     DvIndexPropLikely(1,DiAvailIndexNum)
  }
  DnIndexLikely=0.0;//number                     DnIndexLikely       
     
  //Index Q
  if(DiIndexQBlockNum>0)
       DvIndexQE1Likely.initialize();//vector                     DvIndexQE1Likely(1,DiAvailIndexQNum)
  DnIndexQE1Likely=0.0;
 
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // N for First Year
  DvNYear1Likely.initialize();//vector                     DvNYear1Likely(1,DiNYear1EstParaNum)
  DnNYear1Likely=0.0;//number                     DnNYear1Likely
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Recruitment R
   //R-S Para
  DvRSParaLikely.initialize();//vector                     DvRSParaLikely(1,lL1)
  DnRSParaLikely=0.0;//number                     DnRSParaLikely
  // R Dev SD
  DnRDevsSDLikely =0.0;// number                     DnRDevsSDLikely                  //DiRecruitLogDevsSDEstFlag
  // R Deviation
  DnRecruitDevsLikely =0.0;//number                     DnRecruitDevsLikely 
   //Rh
  DnRecruitRhLikely =0.0;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Natural Mortality
  // A
  DvLorenzenALikely.initialize();//vector                     DvLorenzenALikely(1,DiSeasonNum)   //DiNaturalMortalityFlag
  DnLorenzenALikely =0.0;//number                     DnLorenzenALikely               //DiNaturalMortalityFlag
  //B
  DvLorenzenBLikely.initialize();//vector                     DvLorenzenALikely(1,DiSeasonNum)   //DiNaturalMortalityFlag
  DnLorenzenBLikely =0.0;//number                     DnLorenzenALikely               //DiNaturalMortalityFlag
   //////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Growth Matrix 
  DvLinfLikely.initialize();// vector                     DvLinfLikely(1,DiGrowthMatrixBlockNum) // DiGrowthMatrixFlag
  DvVLikely.initialize();// vector                     DvVLikely(1,DiGrowthMatrixBlockNum)
  DvLSDLikely.initialize();// vector                     DvLSDLikely(1,DiGrowthMatrixBlockNum)
  DvVSDLikely.initialize();// vector                     DvVSDLikely(1,DiGrowthMatrixBlockNum)
  DvLVRhoLikely.initialize();//  vector                     DvLVRhoLikely(1,DiGrowthMatrixBlockNum)
  DnGrowthMatrixLikely=0.0;// number                     DnGrowthMatrixLikely   
  DvLfiftyLikely.initialize();//vector
  DnLfiftyLikely=0.0;//number
  
  DnRsexLikely=0.0;//number
         
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                  
  if(DiLikelyConstFlag)
  {
      //
      for( i=1;i<=DiFleetSelParaNumIni;i++)
      {//
            if(DivFleetSelParaPh(i)>0)
            {
                DvFleetSelLikely(i)=CalLikelihoodConst(DvFleetSelPara_Ini(i), DivFleetSelParaLikelihoodFlag(i))*DvFleetSelParaLambda(i);
                
                DnFleetSelLikely+= DvFleetSelLikely(i);
            }
      }
      //
      iTemp=DiFleetNum*Flagtimestep;
      for(i=1;i<=iTemp;i++)
      {
              if(DivLogFYear1Ph(i)>0)
              {
                DvFYear1Likely(i)=CalLikelihoodConst(DvLogFYear1Ini(i), DivLogFYear1LikelihoodFlag(i))*DvLogFYear1Lambda(i);
                DnFYear1Likely+=DvFYear1Likely(i);
              }    
      }
      //    
      for(i=1;i<=DiFleetNum;i++)
      {
           for(j=DiCalBeginYear;j<=DiCalEndYear;j++)  
           {
                 for(k=1;k<=Flagtimestep;k++)
                 {
                       lL1=(i-1)*(DiCalEndYear-DiCalBeginYear+1)*Flagtimestep+(j-DiCalBeginYear)*Flagtimestep+k; 
                       lL2=(i-1)*Flagtimestep+k;
                       if(DivFLogDevsPh(lL1)>0)
                       {
                             DvFLogDevsLikely(i)+=CalLikelihoodConst(0.0,3)*0;//*DvFLogDevsLambda(lL2); //???????????????
                            
                       }
                 }
            }
            DnFLogDevsLikely+=DvFLogDevsLikely(i);
      }
      // if(DiAvailCPUEQNumP>0)
      //{
            for(i=1;i<=DiCPUEBlockNum;i++)
            {
                 if( DivCPUEQBiomassE1Ph(i)>0)
                 {
                         DvCPUEQE1Likely(i)=CalLikelihoodConst(DvCPUEQBiomassE1Value(i),DivCPUEQBiomassE1LikelihoodFlag(i))*DvCPUEQBiomassE1Lambda(i); // vector                     DvCPUEQE1Likely(1,DiAvailCPUEQNumP)  
                 }
                 DnCPUEQE1Likely+=DvCPUEQE1Likely(i);// 
            }
      //}
      //cout<<"SSSS"<<endl;
      for(i=1;i<=DiFleetNum;i++)
      {
         for(j=1;j<=Flagtimestep;j++)
         {
              for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
              {
                 if(Dd3CatchTotalO(i,j,k)>MINPOSITIVE)
                 {
                     DvCatchTotalLikely(i)+= CalLikelihoodConst(Dd3CatchTotalO(i,j,k),DimCatchTotalLikelihoodFlag(i,j))*DmCatchTotalLambda(i,j);
                 }
              }
              iL1=(i-1)*Flagtimestep+j;
              ii=DivOldIndexAvailCPUE(iL1);
              if(ii<=0||ii>DiAvailCPUENum)
                continue;
              for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
              {
                 if(DimAvailCPUEFlag(i,j) &&  Dd3CPUEO(i,j,k)>MINPOSITIVE)
                 {
                     DvCPUELikely(ii)+= CalLikelihoodConst(Dd3CPUEO(i,j,k),DimCPUELikelihoodFlag(i,j))*DmCPUELambda(i,j);   
                 }
              }
              DnCPUELikely+=DvCPUELikely(ii);       
          }
           DnCatchLikely +=DvCatchTotalLikely(i);
      }
      cout<<"Constant likelihood value of total catch:"<<DnCatchLikely<<endl;
      // exit(-99);
      //cout<<"DvCatchPropLikely 1:"<<DvCatchPropLikely<<endl;
      //cout<<DimCatchCompLikelihoodFlag<<endl;
      for(i=1;i<=DiFleetNum;i++)
      {       
        // cout<<"Step -2:"<<endl;
        for(k=1;k<=Flagtimestep;k++)
        {
           //  cout<<"Step -1:"<<endl;
            if(DimCatchCompLikelihoodFlag(i,k)==1)  //Multinomial distribution 
            {
             // cout<<"Step 0:"<<endl;
              for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
              {
                            iTemp=int(Dd3CatchESSInput(i,k,j)+0.00001);
                            //cout<<"Step 1:"<<iTemp<<endl;
                            if(iTemp>0 && Dd3CatchTotalO(i,k,j)>MINPOSITIVE)
                            {
                                 nTemp=CalLogFactorial(iTemp);
                                 DvCatchPropLikely(i)-=1.0* nTemp;
                                 //cout<<"Step 2:"<<DvCatchPropLikely(i)<<endl;
                                 for(kk=DimCSelStartSizeBin(i,k);kk<=DimCSelEndSizeBin(i,k);kk++)
                                 {
                                        iL1=int(iTemp* Dd4CatchPropAtSizeO(i,k,j,kk)+0.5);
                                        nTemp=CalLogFactorial(iL1);
                                        DvCatchPropLikely(i)+=nTemp*DmCatchCompLambda(i,k);
                                 }
                                 //cout<<"Step 3:"<<DvCatchPropLikely(i)<<endl;
                           }
               }
            }
            else 
            {
                  for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
                  {
                            iTemp=int(Dd3CatchESSInput(i,k,j)+0.00001);  
                            if(iTemp>0  && Dd3CatchTotalO(i,k,j)>MINPOSITIVE)
                            {
                                 for(kk=DimCSelStartSizeBin(i,k);kk<=DimCSelEndSizeBin(i,k);kk++)
                                 {
                                        DvCatchPropLikely(i)+=log(PI2)*0.5*DmCatchCompLambda(i,k);
                                 }
                                //  DnCatchLikely+=DvCatchPropLikely(i);
                            }
                   }
             }
           }  
           DnCatchLikely+=DvCatchPropLikely(i);
       } 
       cout<<"Constant likelihood value of catch composition:"<<DnCatchLikely<<endl;
       //exit(-99);
       for(i=1;i<=DiIndexSelParaNumIni;i++)
       {
             if(DivIndexSelParaPh(i)>0)
             {
                DvIndexSelLikely(i)=CalLikelihoodConst(DvIndexSelPara_Ini(i),DivIndexSelParaLikelihoodFlag(i))*DvIndexSelParaLambda(i);
             }
             DnIndexSelLikely+= DvIndexSelLikely(i);
       }
       for(i=1;i<=DiAvailIndexNum;i++)
       {
          for(k=1;k<=DivIndexObsNum(i);k++)
          {
              
                     DvIndexLikely(i)+= CalLikelihoodConst(DmIndexTotalO(i,k),DivIndexTotalLikelihoodFlag(i))*DvIndexTotalLambda(i);                   
          }
          DnIndexLikely +=DvIndexLikely(i);
       }
       cout<<"Constant likelihood value of survey index:"<<DnIndexLikely<<endl;
       for(i=1;i<=DiAvailIndexNum;i++)
       {
              if(DivEstIndexPropFlag(i)!=1)
                    continue;
              if( DivIndexCompLikelihoodFlag(i)==1)
              {
                for(k=1;k<=DivIndexObsNum(i);k++)
                {
                      iTemp=DimIndexESSInput(i,k);
                      if(iTemp>0)
                      {
                                 nTemp=CalLogFactorial(iTemp);
                                 DvIndexPropLikely(i)-=1.0* nTemp;
                                 for(kk=DivIndexSelStartSizeBin(i);kk<=DivIndexSelEndSizeBin(i);kk++)
                                 {
                                        iL1=int(iTemp* Dd3IndexPropAtSizeO(i,k,kk)+0.5);
                                        nTemp=CalLogFactorial(iL1);
                                        DvIndexPropLikely(i)+=nTemp*DvIndexCompLambda(i);
                                 }
                                //cout<<DvIndexPropLikely(i)<<endl;
                       }
              }
           }
           else 
           {
            // for(i=1;i<=DiAvailIndexNum;i++)
              //{//
                for(k=1;k<=DivIndexObsNum(i);k++)
                {
                                 for(kk=DivIndexSelStartSizeBin(i);kk<=DivIndexSelEndSizeBin(i);kk++)
                                 {
                                        DvIndexPropLikely(i)+=log(PI2)*0.5*DvIndexCompLambda(i);
                                 }
                              //   DnIndexLikely += DvIndexPropLikely(i);
                 }
              //}
           }
           DnIndexLikely += DvIndexPropLikely(i);
         } 
         cout<<"Constant likelihood value of survey composition:"<<DnIndexLikely<<endl;
         for(i=1;i<=DiIndexQBlockNum;i++)
         {
            if(DivQBiomassE1Ph(i)>0)
            {
              DvIndexQE1Likely(i)=CalLikelihoodConst( DvQBiomassE1Value(i),DivQBiomassE1LikelihoodFlag(i))*DvQBiomassE1Lambda(i);
            }
            DnIndexQE1Likely+= DvIndexQE1Likely(i);
         }
         for(i=1;i<=DiNYear1EstParaNum;i++)
         {
                if(DivNYear1ParaPh(i)>0)
                {  
                  DvNYear1Likely(i)=CalLikelihoodConst(DvNYear1ParaInital(i),DivNYear1ParaLikelihoodFlag(i))*DvNYear1ParaLambda(i);
                  DnNYear1Likely+=DvNYear1Likely(i);   
                }
         }
         if(DiRSFlag==1) 
              lL1=1;
         else 
              lL1=2;
         for(i=1;i<=lL1;i++)
         {
                  if(DivRSParaPh(i)>0)
                  {
                      DvRSParaLikely(i)= CalLikelihoodConst( DvRSParaValue(i),DivRSParaLikelihoodFlag(i))*DvRSParaLambda(i);
                      DnRSParaLikely+=DvRSParaLikely(i);
                  }
         }
          
        if(DiRecruitLogDevsSDEstFlag)
        {
                  lL1=DvRecruitLogDevsSD_Ini(4)>0?int(DvRecruitLogDevsSD_Ini(4)+0.01):int(DvRecruitLogDevsSD_Ini(4)-0.01);
                  if(lL1>0)
                    DnRDevsSDLikely=CalLikelihoodConst(DvRecruitLogDevsSD_Ini(1),int(DvRecruitLogDevsSD_Ini(7)+0.01))*DvRecruitLogDevsSD_Ini(6);
        }
        if(DiRecruitLogDevsPhase >0)
        {
                     for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
                     {
                          DnRecruitDevsLikely += CalLikelihoodConst(0.0,3)*DnRecruitLogDevsLambda;//???????????????
                     }
         }
         if(DvRecruitmentRh_Ini(4)>0.01)
         {
                      DnRecruitRhLikely =CalLikelihoodConst(DvRecruitmentRh_Ini(1),int(DvRecruitmentRh_Ini(7)+0.01))*DvRecruitmentRh_Ini(6);
         }
         if(!DiNaturalMortalityFlag)
         {
                 for(i=1;i<=Flagtimestep;i++)
                 {
                       if(DivLorenzenAPh(i)>0)
                       {
                         DvLorenzenALikely(i)=CalLikelihoodConst(DvLorenzenAValue(i),DivLorenzenALikelihoodFlag(i))*DvLorenzenALambda(i);
                         DnLorenzenALikely+=DvLorenzenALikely(i);
                        }
                        if(DivLorenzenBPh(i)>0)
                        {
                            DvLorenzenBLikely(i)=CalLikelihoodConst(DvLorenzenBValue(i),DivLorenzenBLikelihoodFlag(i))*DvLorenzenBLambda(i);
                            DnLorenzenBLikely+=DvLorenzenBLikely(i);
                        } 
                  }
         }
         if(DiGrowthMatrixFlag)
         {
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(DivVBGFLinfPh(i)>0)
                     {
                       DvLinfLikely(i)= CalLikelihoodConst(DvVBGFLinfValue(i),DivVBGFLinfLikelihoodFlag(i))*DvVBGFLinfLambda(i);
                       DnGrowthMatrixLikely+=DvLinfLikely(i);
                     }
                }
                
                 for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(DivVBGFVPh(i)>0)
                     {
                       DvVLikely(i)= CalLikelihoodConst(DvVBGFVValue(i),DivVBGFVLikelihoodFlag(i))*DvVBGFVLambda(i);
                       DnGrowthMatrixLikely+=DvVLikely(i);
                     }
                }
                
                 for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(DivVBGFLSDPh(i)>0)
                     {
                       DvLSDLikely(i)= CalLikelihoodConst(DvVBGFLSDValue(i),DivVBGFLSDLikelihoodFlag(i))*DvVBGFLSDLambda(i);
                       DnGrowthMatrixLikely+=DvLSDLikely(i);
                     }
                }
                 for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(DivVBGFVSDPh(i)>0)
                     {
                       DvVSDLikely(i)= CalLikelihoodConst(DvVBGFVSDValue(i),DivVBGFVSDLikelihoodFlag(i))*DvVBGFVSDLambda(i);
                       DnGrowthMatrixLikely+=DvVSDLikely(i);
                     }
                }
                 for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(DivVBGFLVRhoPh(i)>0)
                     {
                       DvLVRhoLikely(i)= CalLikelihoodConst(DvVBGFLVRhoValue(i),DivVBGFLVRhoLikelihoodFlag(i))*DvVBGFLVRhoLambda(i);
                       DnGrowthMatrixLikely+=DvLVRhoLikely(i);
                     }
                }
       }           
  }//if
  ///////////////////////////////////
  #if Debug_Status
    ICHECKDEB(DnCatchLikely);
    ICHECKDEB(DnCPUELikely);
    ICHECKDEB(DnIndexLikely);
    ICHECKDEB(DnLorenzenALikely);
    ICHECKDEB(DnGrowthMatrixLikely);
    //exit(-99);
  #endif
  // deallocate()
  cout<<"Likelihood Constant Calculation completed"<<endl;
  cout<<""<<endl;
  //////////////////////////////////////////
}

void model_parameters::userfunction(void)
{
  ofvTotal =0.0;
  if (Flagtimestep==1)
    {
      CalNaturalMortality_Year();
      CalFleetSelectivity(); 
      CalFleetFishingMortality_Year(); 
      CalGrowthMatrix();
      CalfemalePropAtsize();
      //cout<<Md3GrowthMatrix<<endl;
      //exit(222);
      CalNumberAtSize_Year();
        if(DiRecruitLogDevsSDEstFlag)
          {
            if(DvRecruitLogDevsSD_Ini(4)<-3.0 && active(PbdvLogRecruitDevs))
              {
                if(max(MbdvLogRecruitDevs)>0.00000001)
                  MnRSD=std_dev(MvRecruitDev);//std_dev(MbdvLogRecruitDevs);//MvRecruitDev(i)
                else
                  MnRSD=0.01;   
              }
         else if(active(PbnvRecruitLogDevsSD2))
           {//PbnvRecruitLogDevsSD2
               MnRSD=PbnvRecruitLogDevsSD2;
           }
           #if Debug_Status
               ICHECKDEB(MbdvLogRecruitDevs);
               ICHECKDEB(MnRSD);
           #endif
          }
      CalPredictedCatchAtSize_Y();
      CalIndexSelectivity();
      //cout<<PbnvIndexQBiomassE1<<endl;
      CalIndexQAndComp_Y();
      //CalfemalePropAtsize();
      CalObjectiveFunc_Year();
      //cout<<sdr_vRecruitmentP<<endl;
      //cout<<Md3NAS<<endl;
      //exit(22);
      if (sd_phase())
     {
      GetBPRSel();  
      CalBiologyRP();
      //cout<<MmBPRSel<<endl;
      //cout<<DmWeightAtSize(4,22)<<endl;
      //cout<<Md4SelByFleet<<endl;
      for (k=DiRecruitPrjVectNP+1;k<=DiSizeBinNum;k++){
      	MvNASPJD_1(k-DiRecruitPrjVectNP)=MvNASPJD(k);
      }
     }
     //exit(214);
    }
  else
    {
      if(DiNaturalMortalityFlag)
        CalNaturalMortality();
      CalFleetSelectivity();
      CalFleetFishingMortality();
      //cout<<Md4FASByFleet<<endl;
      CalGrowthMatrix();
      //cout<<Md3GrowthMatrix<<endl;
      CalfemalePropAtsize();
      CalNumberAtSize();
      //cout<<Md3NAS<<endl;
      //cout<<sdr_vSSB<<endl;
      //cout<<sdr_vRecruitmentP<<endl;
        if(DiRecruitLogDevsSDEstFlag)
          {
            if(DvRecruitLogDevsSD_Ini(4)<-3.0 && active(PbdvLogRecruitDevs))
              {
                if(max(MbdvLogRecruitDevs)>0.00000001)
                  MnRSD=std_dev(MvRecruitDev);//std_dev(MbdvLogRecruitDevs);//MvRecruitDev(i)
                else
                  MnRSD=0.01;  
              }
        else if(active(PbnvRecruitLogDevsSD2))
          {//PbnvRecruitLogDevsSD2
               MnRSD=PbnvRecruitLogDevsSD2;
          }
          #if Debug_Status
              ICHECKDEB(MbdvLogRecruitDevs);
              ICHECKDEB(MnRSD);
          #endif
         }
      CalPredictedCatchAtSize();
      CalPredCPUE();
      CalIndexSelectivity();
      CalIndexQAndComp();
      //cout<<MmIndexP<<endl;
      //exit(100);
      //CalfemalePropAtsize;
      CalObjectiveFunc();
      //exit(2222);
    if (sd_phase())
     {
      //GetBPRSel();  
      //CalBiologyRP();
      for (k=DiRecruitPrjVectNP+1;k<=DiSizeBinNum;k++){
      	MvNASPJD_1(k-DiRecruitPrjVectNP)=MvNASPJD(k);
      }
      //MvNASPJD_1=MvNASPJD;
     }
     }
}

void model_parameters::CalFleetSelectivity(void)
{
  dvariable dvAlpha1;
  dvariable dvBelta1;
  dvariable dvAlpha2;
  dvariable dvBelta2;
  dvariable dvTemp1;
  dvariable dvTemp2;
  dvariable dvMax;
  dvar_vector  dvvP(1,6);
  dvar_vector  dvvSel(1,DiSizeBinNum);
  int iStartSelBin;
  double dBinwidth;
  int iMode;
  Md4SelByFleet.initialize();
  for(i=1;i<=DiFleetNum;i++)
  {
      for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
      {
           for(k=1;k<=Flagtimestep;k++)
           {
              ii=int(Dd3FleetSelBlockFlag(k,j,i+1)+0.01);              //get selectivity block
              if(ii<1||ii>DiFleetSelBlockNum)
                    continue;
             //iL1=DivFleetSelBlockOldMirrorAvail(ii);               //get available block
              if(Dd3CatchTotalO(i,k,j)>MINPOSITIVE)
              {            
                   lL2=DivFleetSelParaIndex_Ini(ii);                 //get location of parameter
                   switch(DivFleetSelBlockFlag(ii))             // whether the selectivity was same as fleet no. x x in ivFleetSelBlockFlag(lL1)
                   {
                     case 1:
                        for(kk=1;kk<=DiSizeBinNum;kk++)
                                Md4SelByFleet(i,k,j,kk)= PbnvFleetSelParams(lL2+kk);
                        break;
                     case 2:
                       dvAlpha1= PbnvFleetSelParams(lL2+1);
                       dvBelta1= PbnvFleetSelParams(lL2+2);//dvBelta1=1.0/ PbnvFleetSelParams(lL2+2);
                       for(kk=1;kk<=DiSizeBinNum;kk++)
                       {
                             nTemp=(DvSizeBinData(kk)+DvSizeBinData(kk-1))*0.5;
                             Md4SelByFleet(i,k,j,kk) =1.0/(1.0+mfexp((dvAlpha1-nTemp)* dvBelta1));
                       }
                       dvMax=max( Md4SelByFleet(i,k,j));
                       Md4SelByFleet(i,k,j)/= dvMax;
                       break;
                     case 3:
                         dvAlpha1= PbnvFleetSelParams(lL2+1);
                         dvBelta1= PbnvFleetSelParams(lL2+2);//dvBelta1=1.0/ PbnvFleetSelParams(lL2+2);
                         dvAlpha2= PbnvFleetSelParams(lL2+3);
                         dvBelta2= PbnvFleetSelParams(lL2+4);// dvBelta2=1.0/ PbnvFleetSelParams(lL2+4);
                         for(kk=1;kk<=DiSizeBinNum;kk++)
                         {     
                                  nTemp=(DvSizeBinData(kk)+DvSizeBinData(kk-1))*0.5;
                                  dvTemp1 =1.0/(1.0+mfexp((dvAlpha1-nTemp)* dvBelta1));
                                  dvTemp2 =1.0-1.0/(1.0+mfexp((dvAlpha2-nTemp)* dvBelta2));
                                  Md4SelByFleet(i,k,j,kk)=dvTemp1*dvTemp2;
                         }
                         dvMax=max( Md4SelByFleet(i,k,j));
                         Md4SelByFleet(i,k,j)/= dvMax;
                         break;   
                        case 4:
                         iMode=0;
                         dBinwidth=1.0;
                         iStartSelBin=1.0;
                         dvvP(1)= PbnvFleetSelParams(lL2+1);
                         dvvP(2)= PbnvFleetSelParams(lL2+2);//dvBelta1=1.0/ PbnvFleetSelParams(lL2+2);
                         dvvP(3)= PbnvFleetSelParams(lL2+3);
                         dvvP(4)= PbnvFleetSelParams(lL2+4);// dvBelta2=1.0/ PbnvFleetSelParams(lL2+4);
                         GetSS3DNSelectivity(dvvP,dvvSel,iStartSelBin,dBinwidth,iMode);
                         for(kk=1;kk<=DiSizeBinNum;kk++)
                         {     
                                  Md4SelByFleet(i,k,j,kk)=dvvSel(kk);
                         }
                         dvMax=max( Md4SelByFleet(i,k,j));
                         Md4SelByFleet(i,k,j)/= dvMax;
                         break;    
                        case 5:
                         iMode=1;
                         dBinwidth=1.0;
                         iStartSelBin=1.0;
                         dvvP(1)= PbnvFleetSelParams(lL2+1);
                         dvvP(2)= PbnvFleetSelParams(lL2+2);//dvBelta1=1.0/ PbnvFleetSelParams(lL2+2);
                         dvvP(3)= PbnvFleetSelParams(lL2+3);
                         dvvP(4)= PbnvFleetSelParams(lL2+4);// dvBelta2=1.0/ PbnvFleetSelParams(lL2+4);
                         GetSS3DNSelectivity(dvvP,dvvSel,iStartSelBin,dBinwidth,iMode);
                         for(kk=1;kk<=DiSizeBinNum;kk++)
                         {     
                                  Md4SelByFleet(i,k,j,kk)=dvvSel(kk);
                         }
                         dvMax=max( Md4SelByFleet(i,k,j));
                         Md4SelByFleet(i,k,j)/= dvMax;
                         break;                  
                     }
                  }
               }
     }
  }
 #if Debug_Status
 // ICHECKDEB(Md4SelByFleet);
 #endif
  return ;
}

void model_parameters::GetSS3DNSelectivity(const dvar_vector dvvP,dvar_vector &dvvSel,int iStartSelBin, double dBinwidth,int iMode)
{
  dvariable  dvPeak1;
  dvariable  dvPeak2;
  dvariable  dvUpSelx;
  dvariable  dvDownSelx;
  int        i;
  dvariable  dvTemp1;
  double     dLastBinSize;
  dLastBinSize=(DvSizeBinData(DiSizeBinNum) +DvSizeBinData(DiSizeBinNum-1))*0.5;
  if(iMode==0)
     dBinwidth=0.0;
  dvPeak1=dvvP(1);
  dvPeak2=dvPeak1+dBinwidth+(0.99*dLastBinSize-dvPeak1-dBinwidth)/(1.0+mfexp(-dvvP(2)));
  dvUpSelx=mfexp(dvvP(3));
  dvDownSelx=mfexp(dvvP(4));
  if(iMode==0)
  {
    for(i=1;i<=DiSizeBinNum;i++)
    {
        dLastBinSize=DvSizeBinData(i);
        if(dLastBinSize<dvPeak1)
        {
           dvvSel(i)=mfexp(-square(dLastBinSize-dvPeak1)/dvUpSelx);
        }
        else if(dLastBinSize< dvPeak2)
        {
             dvvSel(i)=1.0;
        }
        else
        {
             dvvSel(i)=mfexp(-square(dLastBinSize-dvPeak2)/dvDownSelx);
        }
    }
  }
  else if(iMode==1)
  {
       dvariable  dvFinal;
       dvariable  dvInit;
       int        j1;
       int        j2;
       dvariable  dvPoint1;
       dvariable  dvPoint2;
       dvariable  dvT1min;
       dvariable  dvT2min;
       dvariable  dvT1;
       dvariable  dvT2;
       dvariable  dvJoin1;
       dvariable  dvJoin2;
       dvariable  dvAsc;
       dvariable  dvDsc;
       dvFinal=dvvP(6);
       dvInit=dvvP(5);
       if(dvInit<-1000.0)
       {
            j1=-1001-int(value(dvInit)); //????????????
            if(j1>=1)
              dvvSel(1,j1)=1.0e-06;
            else
              j1=0;
       }
       else
       {
             j1=iStartSelBin-1; 
             if(dvInit>-999.)
             {
                 dvPoint1=1.0/(1.0+mfexp(-dvInit));
                 dLastBinSize=DvSizeBinData(iStartSelBin);
                 dvT1min=mfexp(-(square(dLastBinSize-dvPeak1)/dvUpSelx)); 
             }
       }
       if(dvFinal<-1000.0)
       {
           j2=-1000-int(value(dvFinal));
       }
       else
       {
          j2=DiSizeBinNum;
          if(dvFinal>-999)
          {
              dvPoint2=1.0/(1.0+mfexp(-dvFinal));
              dLastBinSize=DvSizeBinData(j2);
              dvT2min=mfexp(-(square(dLastBinSize-dvPeak2)/dvDownSelx));  // fxn at last bin
           }
       }
       for(j=j1+1;j<=j2;j++)
       {
              dLastBinSize=DvSizeBinData(j);
              dvT1=dLastBinSize-dvPeak1;  
              dvT2=dLastBinSize-dvPeak2;
              dvJoin1=1.0/(1.0+mfexp(-(20.*dvT1/(1.0+fabs(dvT1)))));  
              dvJoin2=1.0/(1.0+mfexp(-(20.*dvT2/(1.0+fabs(dvT2)))));
              if(dvInit>-999.0)
              {
                 dvAsc=dvPoint1+(1.0-dvPoint1)*(mfexp(-square(dvT1)/dvUpSelx)-dvT1min)/(1.0-dvT1min);
              }
              else
              { 
                  dvAsc=mfexp(-square(dvT1)/dvUpSelx);
              }
              if(dvFinal>-999.0)
              {
                 dvDsc=1.0+(dvPoint2-1.0)*(mfexp(-square(dvT2)/dvDownSelx)-1.0)/(dvT2min-1.0);
              }
              else
              {
                   dvDsc=mfexp(-square(dvT2)/dvDownSelx);
              }
              dvvSel(j)=dvAsc*(1.0-dvJoin1)+dvJoin1*(1.0-dvJoin2+dvDsc*dvJoin2);
       }
       if(iStartSelBin>1 && dvInit>=-1000.)
       {
              for(j=1;j<=iStartSelBin-1;j++)
              {
                  dLastBinSize=DvSizeBinData(j);
                  dvvSel(j)=square(dLastBinSize/DvSizeBinData(iStartSelBin))*dvvSel(iStartSelBin);
              }
       }
       if(j2<DiSizeBinNum) 
       {
               dvvSel(j2+1,DiSizeBinNum)=dvvSel(j2);
       }
  }
}

void model_parameters::CalGrowthMatrix(const dvariable &dvLinf,const dvariable dvK,const dvariable &dvLSigma,const dvariable &dvKSigma,const dvariable &dvRho,double dAlpha,int iStepX,int iStepY,const int iBlock)
{
  int i,j,k;
  dvariable dvEMinusK;
  dvariable dvEMinusK2;
  dvariable dv1EMinusK;
  dvariable dv1EMinusK2;
  dvariable dvTemp2;
  dvariable fun_dvIntegral;
  double    dLCurrent;
  double    dLDelta;
  #if Debug_Status
  #endif
  dvEMinusK=mfexp(-dAlpha*dvK);
  dvEMinusK2=dvEMinusK*dvEMinusK;
  dv1EMinusK=1.0-dvEMinusK;
  dv1EMinusK2=dv1EMinusK*dv1EMinusK;
  //cout<<DvSizeBinData<<endl;
  for(i=1;i<=DiSizeBinNum;i++)
  {
       if(iStepY<=1)
       {
          dLDelta=DvSizeBinData(i)- DvSizeBinData(i-1);//increment
          dLCurrent= DvSizeBinData(i-1)+ dLDelta*0.5;//mid-point of a given size class
          MnDeltaL=(dvLinf- dLCurrent)*dv1EMinusK+dLCurrent;
          //MnDeltaLVar=dvLSigma*dvLSigma* dv1EMinusK2+(dvLinf-dLCurrent)*(dvLinf-dLCurrent)*dvKSigma*dvKSigma*dvEMinusK2*dAlpha*dAlpha;
          //MnDeltaLVar+= 2*dvRho*dvLSigma*dvKSigma*dv1EMinusK*(dvLinf- dLCurrent)*dvEMinusK*dAlpha;
          //cout<<MnDeltaL<<endl;
          //cout<<MnDeltaLVar<<endl;
          for(k=i;k<=DiSizeBinNum;k++)
          {
             dLDelta=DvSizeBinData(k)- DvSizeBinData(k-1);//increment
             dLCurrent= DvSizeBinData(k-1)+ dLDelta*0.5;
             MnDeltaLVar=dvLSigma*dvLSigma* dv1EMinusK2+(dvLinf-dLCurrent)*(dvLinf-dLCurrent)*dvKSigma*dvKSigma*dvEMinusK2*dAlpha*dAlpha;
             MnDeltaLVar+= 2*dvRho*dvLSigma*dvKSigma*dv1EMinusK*(dvLinf- dLCurrent)*dvEMinusK*dAlpha;
             MnZscore1=((DvSizeBinData(k-1))-MnDeltaL)/sqrt(MnDeltaLVar);
             MnZscore2=((DvSizeBinData(k))-MnDeltaL)/sqrt(MnDeltaLVar);
             //fun_dvIntegral=Integrate(MnDeltaL,MnDeltaLVar,DvSizeBinData(k-1), DvSizeBinData(k),iStepX,2);
             //Md3GrowthMatrix(iBlock,i,k)= fun_dvIntegral;
             Md3GrowthMatrix(iBlock,i,k)=cumd_norm(MnZscore2)-cumd_norm(MnZscore1);
          }
       }
       else
       {
          dLDelta=(DvSizeBinData(i)- DvSizeBinData(i-1))/iStepY;
          for(j=1;j<=iStepY;j++)
          {
             dLCurrent= DvSizeBinData(i-1)+ dLDelta*(j-0.5);
             MnDeltaL=(dvLinf- dLCurrent)*dv1EMinusK+dLCurrent;
             MnDeltaLVar=dvLSigma*dvLSigma* dv1EMinusK2+(dvLinf- dLCurrent)*(dvLinf- dLCurrent)*dvKSigma*dvKSigma* dvEMinusK2*dAlpha*dAlpha;
             MnDeltaLVar+= 2*dvRho*dvLSigma*dvKSigma*dv1EMinusK*(dvLinf- dLCurrent)*dvEMinusK*dAlpha;
             for(k=i;k<=DiSizeBinNum;k++)
             {
                fun_dvIntegral=fun_dvIntegral=Integrate(MnDeltaL,MnDeltaLVar,DvSizeBinData(k-1), DvSizeBinData(k),iStepX,1);
                Md3GrowthMatrix(iBlock,i,k)+= fun_dvIntegral;
             }
          }    
       }
       dvTemp2=sum( Md3GrowthMatrix(iBlock,i));
       if(dvTemp2>0)
            Md3GrowthMatrix(iBlock,i)/=dvTemp2;
  }
}

void model_parameters::CalFleetFishingMortality(void)
{
  dvariable   logTemp;
  dvariable   logSTemp;
  ////////////////////////////////
  Md3Ftot.initialize();
  Md4FASByFleet.initialize();
  Md3Ft.initialize();
  Md3Z.initialize();
  Md3S.initialize();
  MmSSB_S.initialize();
  MvSSB_S1.initialize();
  for(i=1;i<=DiFleetNum;i++)
  {
      for(k=1;k<=DiSeasonNum;k++)
      {
           logTemp=0.0;
           iL1=DimFleetFirstYear(i,k);     // 
           if(iL1<=0)
             continue; 
           for(j=iL1;j<=DiCalEndYear;j++)
           {
               if(j==iL1)
               {
                    logTemp=PbnvLogFYear1Season1((i-1)*DiSeasonNum+k);
                    logSTemp=logTemp;
               }
               else
               {
                  iL2=(i-1)*(DiCalEndYear-DiCalBeginYear+1)*DiSeasonNum + (j-DiCalBeginYear)*DiSeasonNum +k;
                  lL2=(j-DiCalBeginYear)*DiSeasonNum +k;
                  if(DimFleetFZeroIndex(i,lL2)!=0) //if(DivFLogDevsPh(iL2)>=0)
                  {
                    logTemp=logSTemp;
                    logTemp= logTemp+ PbnvLogFDevs(iL2);
                    logSTemp=logTemp;
                  }
                  else
                  {
                    logTemp=-27.0;
                  }
               }
               if(logTemp>-26.0)
                   Md3Ft(i,k,j)= mfexp(logTemp);
               for(kk=1;kk<=DiSizeBinNum;kk++)
               {//DimCSelStartSizeBin(i,k)   DimCSelEndSizeBin(i,k)
                     Md4FASByFleet(i,k,j,kk)=Md3Ft(i,k,j)*Md4SelByFleet(i,k,j,kk);
                     Md3Ftot(k,j,kk)+=Md4FASByFleet(i,k,j,kk);
               }
          }
      }
  }
  ///////////////////////////////////////////////////
  for(i=1;i<=DiSeasonNum;i++)
  {
        for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
        {
             for(k=1;k<=DiSizeBinNum;k++)
             {
                 Md3Z(i,j,k)= Md3Ftot(i,j,k)+Md3M(i,j,k);
             }
        }      
  }
  Md3S=mfexp(-1.0*Md3Z);//?-1.0
  //// DiSBBSeason   DnFracSeasonSSB
  for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
  {
             for(k=1;k<=DiSizeBinNum;k++)
             {
                 MmSSB_S(j,k)=mfexp(-1.0*DnFracSeasonSSB*Md3Z(DiSBBSeason,j,k));
             }
  }  
  for(k=1;k<=DiSizeBinNum;k++)
             {
                 MvSSB_S1(k)=mfexp(-1.0*DnFracSeasonSSB*Md3Z(1,1,k));
             }
  //MvSSB_S1=mfexp(-1.0*DnFracYearSSB*Md3Z(1,1));
  //cout<<Md3Z(1,1)<<endl;
 #if Debug_Status
   ICHECKDEB(Md3Ftot);
   ICHECKDEB(Md3Z);
   ICHECKDEB(Md4FASByFleet);
   ICHECKDEB(Md3Ft);
 #endif
  return ;
}

void model_parameters::CalFleetFishingMortality_Year(void)
{
  dvariable   logTemp;
  dvariable   logSTemp;
  ////////////////////////////////
  Md3Ftot.initialize();
  Md4FASByFleet.initialize();
  Md3Ft.initialize();
  Md3Z.initialize();
  Md3S.initialize();
  MmSSB_Y.initialize();
  for(i=1;i<=DiFleetNum;i++)
  {
      for(k=1;k<=Flagtimestep;k++)
      {
           logTemp=0.0;
           iL1=DimFleetFirstYear(i,k);     // 
           if(iL1<=0)
             continue; 
           for(j=iL1;j<=DiCalEndYear;j++)
           {
               if(j==iL1)
               {
                    logTemp=PbnvLogFYear1Season1(i);
                    logSTemp=logTemp;
               }
               else
               {
                  iL2=(i-1)*(DiCalEndYear-DiCalBeginYear+1) + j;
                  if(DimFleetFZeroIndex(i,j)!=0) //if(DivFLogDevsPh(iL2)>=0)
                  {
                    logTemp=logSTemp;
                    logTemp= logTemp+ PbnvLogFDevs(iL2);
                    logSTemp=logTemp;
                    //cout<<logSTemp<<endl;
                  }
                  else
                  {
                    logTemp=-27.0;
                  }
               }
               if(logTemp>-26.0)
                   Md3Ft(i,k,j)= mfexp(logTemp);
                   //cout<<Md3Ft<<endl;
               for(kk=1;kk<=DiSizeBinNum;kk++)
               {//DimCSelStartSizeBin(i,k)   DimCSelEndSizeBin(i,k)
                     Md4FASByFleet(i,k,j,kk)=Md3Ft(i,k,j)*Md4SelByFleet(i,k,j,kk);
                     Md3Ftot(k,j,kk)+=Md4FASByFleet(i,k,j,kk);
               }
          }
      }
  }
  ///////////////////////////////////////////////////
  for(i=1;i<=Flagtimestep;i++)
  {
        for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
        {
             for(k=1;k<=DiSizeBinNum;k++)
             {
                 Md3Z(i,j,k)= Md3Ftot(i,j,k)+Md3M(i,j,k);
             }
        }      
  }
  Md3S=mfexp(-1.0*Md3Z);//?-1.0
  //// DiSBBSeason   DnFracSeasonSSB
  for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
  {
             for(k=1;k<=DiSizeBinNum;k++)
             {
                 MmSSB_Y(j,k)=mfexp(-1.0*DnFracYearSSB*Md3Z(1,j,k));
             }
  }     
  //mSSB_S=mfexp(-1.0*nFracYearSSB*mZ);
 #if Debug_Status
   ICHECKDEB(Md3Ftot);
   ICHECKDEB(Md3Z);
   ICHECKDEB(Md4FASByFleet);
   ICHECKDEB(Md3Ft);
 #endif
  return ;
}

void model_parameters::GetNYear1Season1(void)
{
       dvariable dvTempN;
       dvariable dvTempSigma;
       dvariable dvTempU;
       dvariable dvTempSigma1;
       dvariable dvTempU1;
       dvariable dvTempSigma2;
       dvariable dvTempU2;
       dvariable dvpi1;
       dvariable dvpi2;
       dvariable sump;
       double    dL;
       int       i;
       MvNASTemp.initialize();
       //cout<<DiNYear1ChoiceFlag<<endl;
       //exit(-25);
       //cout<<"DiNYear1ChoiceFlag:"<<DiNYear1ChoiceFlag<<endl;
       switch(DiNYear1ChoiceFlag)
       {
          case 0:
                for(i=1;i<=DiSizeBinNum;i++)
                        MvNASTemp(i)= PbvNYear1Para(i);
                break;
          case 1:
                dvTempN= PbvNYear1Para(1);
                for(i=1;i<=DiSizeBinNum;i++)
                       MvNASTemp(i)=DvNYear1Pia(i)* dvTempN;
                break;
          case 2:
               dvTempN= PbvNYear1Para(1);
               for(i=1;i<=DiSizeBinNum;i++)
                       MvNASTemp(i)=mfexp(DvNYear1Pia(i))/(1.0+mfexp(DvNYear1Pia(i)))* dvTempN;
               break;
          case 3://log-normal distribution 
               dvTempN=PbvNYear1Para(1);
               dvTempU=PbvNYear1Para(2);
               dvTempSigma=PbvNYear1Para(3);
               for(i=1;i<=DiSizeBinNum;i++)
               {
                       dL=0.5*(DvSizeBinData(i-1)+ DvSizeBinData(i));
                       MnTemp1=1.0/(sqrt(PI2)*dvTempSigma*dL)*mfexp(-(log(dL)-dvTempU)*(log(dL)-dvTempU)/( dvTempSigma* dvTempSigma*2.0));
                       MvNASTemp(i)=MnTemp1*dvTempN;
                       //cout<< dL<<"\t"<<dvTempU<<"\t"<<dvTempSigma<<"\t"<<MnTemp1<<"\t"<<MvNASTemp(i)<<endl;
               }
               //exit(-99);
               break;
          case 4://log-normal distribution
               dvTempN=PbvNYear1Para(1);
               dvTempU=PbvNYear1Para(2);
               dvTempSigma=PbvNYear1Para(3);
               for(i=1;i<=DiSizeBinNum;i++)
               {
                       dL=0.5*(DvSizeBinData(i-1)+ DvSizeBinData(i));
                       MnTemp1=1.0/(sqrt(PI2)*dvTempSigma*dL)*mfexp(-(log(dL)-dvTempU)*(log(dL)-dvTempU)/( dvTempSigma* dvTempSigma*2.0));
                       MvNASTemp(i)=mfexp( MnTemp1)/(1.0+mfexp( MnTemp1))* dvTempN;
               }
               break;
           case 5://normal distribution 
               dvTempN=PbvNYear1Para(1);
               dvTempU=PbvNYear1Para(2);
               dvTempSigma=PbvNYear1Para(3);
               for(i=1;i<=DiSizeBinNum;i++)
               {
                       dL=0.5*(DvSizeBinData(i-1)+ DvSizeBinData(i));
                       MnTemp1=1.0/(sqrt(PI2)*dvTempSigma)*mfexp(-(dL-dvTempU)*(dL-dvTempU)/( dvTempSigma* dvTempSigma*2.0));
                       MvNASTemp(i)=MnTemp1*dvTempN;
               }
               break;
           case 6://normal distribution 
               dvTempN=	PbvNYear1Para(1);
               dvTempU=	PbvNYear1Para(2);
               dvTempSigma= PbvNYear1Para(3);
               for(i=1;i<=DiSizeBinNum;i++)
               {
                       dL=0.5*(DvSizeBinData(i-1)+ DvSizeBinData(i));
                       MnTemp1=1.0/(sqrt(PI2)*dvTempSigma)*mfexp(-(dL-dvTempU)*(dL-dvTempU)/( dvTempSigma* dvTempSigma*2.0));
                       MvNASTemp(i)=mfexp( MnTemp1)/(1.0+mfexp( MnTemp1))* dvTempN;
               }
               break;
           case 7://mixtrue distributions
               dvTempN=	PbvNYear1Para(1);
               dvTempU=	PbvNYear1Para(2);
               dvTempSigma= PbvNYear1Para(3);
               dvTempU1=PbvNYear1Para(4);
               dvTempSigma1= PbvNYear1Para(5);
               dvTempU2=PbvNYear1Para(6);
               dvTempSigma2= PbvNYear1Para(7);
               dvpi1=PbvNYear1Para(8);
               dvpi2=PbvNYear1Para(9);
               for(i=1;i<=DiSizeBinNum;i++)
               {
                       dL=0.5*(DvSizeBinData(i-1)+ DvSizeBinData(i));
                       MnTemp1=1.0/(sqrt(PI2)*dvTempSigma)*mfexp(-(dL-dvTempU)*(dL-dvTempU)/( dvTempSigma* dvTempSigma*2.0));
                       MnTemp2=1.0/(sqrt(PI2)*dvTempSigma1)*mfexp(-(dL-dvTempU1)*(dL-dvTempU1)/( dvTempSigma1* dvTempSigma1*2.0));
                       MnTemp3=1.0/(sqrt(PI2)*dvTempSigma2)*mfexp(-(dL-dvTempU2)*(dL-dvTempU2)/( dvTempSigma2* dvTempSigma2*2.0));
                       MvNASTemp(i)=(dvpi1*MnTemp1+dvpi2*MnTemp2+(1-dvpi1-dvpi2)*MnTemp3);
               }
               sump=sum(MvNASTemp);
               for(i=1;i<=DiSizeBinNum;i++)
               {
                       MvNASTemp(i)= (MvNASTemp(i)/sump)*dvTempN;
               }
               break;
       }
  #if Debug_Status
   //ICHECKDEB(MvNASTemp);
 #endif
       return;
}

void model_parameters::CalGrowthMatrix(void)
{
  /*int k;
  if(DiGrowthMatrixFlag)
  {
      if(active(PbnvLinf(1))||active(PbnvV(1))||active(PbnvLSD(1))||active(PbnvVSD(1)))
      {
        //DiGrowthMatrixFirstCalFlag=0;
        Md3GrowthMatrix.initialize();
        for(k=1;k<=DiGrowthMatrixBlockNum;k++)
             CalGrowthMatrix(PbnvLinf(k),PbnvV(k), PbnvLSD(k),PbnvVSD(k),PbnvLVRho(k),1.0,30,0,k);
      }
 #if Debug_Status
  //ICHECKDEB(Md3GrowthMatrix);
 #endif
  }*/
  int k,iL1,j;
  if(DiGrowthMatrixFlag)
  {
        for(k=1;k<=DiGrowthMatrixBlockNum;k++)
        {
            if(active(PbnvLinf(k))||active(PbnvV(k))||active(PbnvLSD(k))||active(PbnvVSD(k)))
            {
                   if(DiGrowthMatrixFlag==2)
                   {
                         for(j=1;j<=DiSeasonNum;j++)
                         {
                              iL1=(k-1)*DiSeasonNum+j;
                              Md3GrowthMatrix(iL1).initialize();
                              CalGrowthMatrix(PbnvLinf(k),PbnvV(k), PbnvLSD(k),PbnvVSD(k),PbnvLVRho(k),DvGrowthTimeAsYear(j),20,0,iL1);
                         }
                    }
                    else
                    {
                      Md3GrowthMatrix(k).initialize();
                      CalGrowthMatrix(PbnvLinf(k),PbnvV(k), PbnvLSD(k),PbnvVSD(k),PbnvLVRho(k),1.0,30,0,k);
                    }
            }
        }
 #if Debug_Status
  //ICHECKDEB(Md3GrowthMatrix);
 #endif
  }
}

void model_parameters::CalNumberAtSize(void)
{
  // R-S
  dvariable    dvAlpha,dvBeta,dvR00;
  if(DiRSFlag ==1)
  {
        dvAlpha=mfexp(PbnvLogRSPara(1));
        dvBeta=0.0;
  }
  else
  {
        dvAlpha=mfexp(PbnvLogRSPara(1));
        dvBeta=mfexp(PbnvLogRSPara(2)); 
  }
  ///////////////////////////////////////////////
  Md3NAS.initialize();
  sdr_vNBiomass.initialize();
  MvRecruitDev.initialize();
  MvRecruitRS.initialize();
  sdr_vSSB.initialize();
  ////////////////////////////////////////////////////////////
  MvRecruitPrjVect= PbnvRecruitPrjVect/sum(PbnvRecruitPrjVect);
  /////////////////////////////////N(1,DiCalBeginYear,k)////////////////////////////////////
  GetNYear1Season1();
  for (i=1;i<=DiSizeBinNum;i++)  
             Md3NAS(1,DiCalBeginYear,i)=MvNASTemp(i);
  //cout<<Md3NAS(1,1)<<endl;
  ////////////////////////////////////////////////////////////////////////////////
  sdr_vSSB(1)=elem_prod(Md3NAS(1,1),MvSSB_S1)*DmFecundity(1); 
  sdr_vSSB2(1)=elem_prod(elem_prod(Md3NAS(1,1),MvSSB_S1),DmWeightAtSize(1))*Md3FemalePropAtSizeP(1,1)*0.5;
  //cout<<sdr_vSSB<<endl;
  //cout<<sdr_vSSB2<<endl;
  MvRecruitDev(DiCalBeginYear)=PbdvLogRecruitDevs(DiCalBeginYear);
  if(DvLfiftyParas_Ini(4)<=0&DvRsexParas_Ini(4)<=0)
     {
       dvR00=CalStockRecruitmentF(DmEnv(1,1),DmEnv(1,2),DmEnv(1,3),sdr_vSSB(1),dvAlpha,dvBeta,PbvEnvCoef,DiRSFlag)*mfexp(MvRecruitDev(DiCalBeginYear)); // MbnvR00Para(1)=MvNASTemp*DmFecundity(DiCalEndYear);
     }
     else
     {       dvR00=CalStockRecruitmentF(DmEnv(1,1),DmEnv(1,2),DmEnv(1,3),sdr_vSSB2(1),dvAlpha,dvBeta,PbvEnvCoef,DiRSFlag)*mfexp(MvRecruitDev(DiCalBeginYear)); // MbnvR00Para(1)=MvNASTemp*DmFecundity(DiCalEndYear);
     }
  /////////////////////////////////////////////////////////////////////
  MvRecruitRS(DiCalBeginYear)=dvR00;
  //cout<<MvRecruitRS<<endl;
  ///////////////////////////////////////Population dynamics//////////////////////////////////////////////
  for(i=DiCalBeginYear;i<=DiCalEndYear;i++)
  {
     if(i>DiCalBeginYear)
     {//because first year R=0 so PbdvLogRecruitDevs begin DiCalBeginYear+1
          if(i<=DiCalEndYear-DiYearBeforeEndForRDev)
          {
              MvRecruitDev(i)=sqrt(PbnRecruitmentRh+0.0000000001)*MvRecruitDev(i-1)+sqrt(1.0000000001-PbnRecruitmentRh)*PbdvLogRecruitDevs(i);
              MbdvLogRecruitDevs(i)=PbdvLogRecruitDevs(i);
            // cout<<PbdvLogRecruitDevs(i)<<endl;
          }
          else
          {
              MvRecruitDev(i)=sqrt(PbnRecruitmentRh+0.0000000001)*MvRecruitDev(i-1)+sqrt(1.0000000001-PbnRecruitmentRh)*PbdvLogRecruitDevsADD(i);
              MbdvLogRecruitDevs(i)=PbdvLogRecruitDevsADD(i);
          }
          sdr_vRecruitmentP(i)=MvRecruitRS(i)*mfexp(MvRecruitDev(i));//-0.5*MnRSD*MnRSD); // Convergence is difficult for sd calculated from MvRecruitDev when MnRSD is included 
         //cout<<"R SD is over"<<endl;
     } 
     else
         sdr_vRecruitmentP(i)=MvRecruitRS(i);//*mfexp(PbdvLogRecruitDevs(i));
     //cout<<MvRecruitDev<<endl;
     for(j=1;j<=DiSeasonNum;j++)
     {
             for (k=1;k<=DiSizeBinNum;k++)  
             {//Recruitment    
                   if(i>=DiCalBeginYear)   //MvNASTemp(i);
                        Md3NAS(j,i,k)= Md3NAS(j,i,k) +sdr_vRecruitmentP(i)*DvRecruitmentSeaRate(j)*MvRecruitPrjVect(k);
                   sdr_vNBiomass(i)+=Md3NAS(j,i,k)*DmWeightAtSize(i,k);
             }
             if(j==DiSBBSeason)
             {//SSB  
                sdr_vSSB(i)=elem_prod(Md3NAS(j,i),MmSSB_S(i))*DmFecundity(i);   //sgwj debug
                //cout<<MmSSB_S(1)<<endl;
                //cout<<Md3NAS(1,1)<<endl;
                //cout<<MmSSB_Y(i)<<endl;
                sdr_vSSB2(i)=elem_prod(elem_prod(Md3NAS(j,i),MmSSB_S(i)),DmWeightAtSize(i))*Md3FemalePropAtSizeP(1,i)*0.5;
                if(DvLfiftyParas_Ini(4)<=0&DvRsexParas_Ini(4)<=0)
                {
                MvRecruitRS(i+1)=CalStockRecruitmentF(DmEnv(i,1),DmEnv(i,2),DmEnv(i,3),sdr_vSSB(i),dvAlpha,dvBeta,PbvEnvCoef,DiRSFlag);
                }
                else
                {
                MvRecruitRS(i+1)=CalStockRecruitmentF(DmEnv(i,1),DmEnv(i,2),DmEnv(i,3),sdr_vSSB2(i),dvAlpha,dvBeta,PbvEnvCoef,DiRSFlag);
                }
                //cout<<sdr_vSSB2<<endl;
                //MvRecruitRS(i+1)=CalStockRecruitment(sdr_vSSB(i),dvAlpha,dvBeta,DiRSFlag);
             }
             //Die by Catch and natural mortality 
             MvNASTemp=elem_prod(Md3NAS(j,i),Md3S(j,i));//die        N(m+1,t)=(N(m,t)*exp(-Z)*G)
             //cout<<MvNASTemp<<endl;
             iL1=DimGrowthMatrixBlockFlag(i,j+1);
             if(DiGrowthMatrixFlag==2)
                  iL1=(iL1-1)*DiSeasonNum+j;
             //if(iL1<1||iL1>DiGrowthMatrixBlockNum||(DiGrowthMatrixFlag==2&&(iL1<1||iL1>DiGrowthMatrixBlockNum*DiSeasonNum)))
             if((DiGrowthMatrixFlag!=2 && (iL1<1||iL1>DiGrowthMatrixBlockNum))||(DiGrowthMatrixFlag==2&&(iL1<1||iL1>DiGrowthMatrixBlockNum*DiSeasonNum)))
             {
                  WARNING<<"No Growth Matrix was used"<<endl;
                  if(j<DiSeasonNum)
                        Md3NAS(j+1,i)= MvNASTemp;           //Growth Next Season
                  else
                 {
                   if(i+1<=DiCalEndYear)
                         Md3NAS(1,i+1)= MvNASTemp;           //Growth Next Year
                    else
                         MvNASPJD     = MvNASTemp; 
                  }
             }
             else
             {
               if(j<DiSeasonNum)
                     Md3NAS(j+1,i)= MvNASTemp*Md3GrowthMatrix(iL1);           //Growth Next Season
               else
               {
                  if(i+1<=DiCalEndYear)
                     Md3NAS(1,i+1)= MvNASTemp*Md3GrowthMatrix(iL1);           //Growth Next Year
                  else
                     MvNASPJD     = MvNASTemp*Md3GrowthMatrix(iL1); 
               }
            }
            //cout<<Md3NAS<<endl;
            //exit(110);
       }
   }
   //////////////////////////////////////
 #if Debug_Status
    //ICHECKDEB(MvRecruitRS);
    //ICHECKDEB(sdr_vRecruitmentP);
    ICHECKDEB(Md3NAS);
    //ICHECKDEB(Md3NAS(2));
    //ICHECKDEB(Md3NAS(3));
    //ICHECKDEB(Md3NAS(4));
 #endif
  ///////////////////////End/////////////////////////////////////////////////////////////////////
  return;
}

void model_parameters::CalNumberAtSize_Year(void)
{
  // R-S
  dvariable    dvAlpha,dvBeta,dvR00;
  if(DiRSFlag ==1)
  {
        dvAlpha=mfexp(PbnvLogRSPara(1));
        dvBeta=0.0;
  }
  else
  {
        dvAlpha=mfexp(PbnvLogRSPara(1));
        dvBeta=mfexp(PbnvLogRSPara(2)); 
  }
  ///////////////////////////////////////////////
  Md3NAS.initialize();
  sdr_vNBiomass.initialize();
  MvRecruitDev.initialize();
  MvRecruitRS.initialize();
  sdr_vSSB.initialize();
  sdr_vSSB2.initialize();
  ////////////////////////////////////////////////////////////
  MvRecruitPrjVect= PbnvRecruitPrjVect/sum(PbnvRecruitPrjVect);
  /////////////////////////////////N(1,DiCalBeginYear,k)////////////////////////////////////
  GetNYear1Season1();
  for (i=1;i<=DiSizeBinNum;i++)  
             Md3NAS(1,DiCalBeginYear,i)=MvNASTemp(i);
  //cout<<Md3NAS<<endl;
  sdr_vSSB(1)=elem_prod(Md3NAS(1,1),MmSSB_Y(1))*DmFecundity(1); 
  sdr_vSSB2(1)=elem_prod(elem_prod(Md3NAS(1,1),MmSSB_Y(1)),DmWeightAtSize(1))*Md3FemalePropAtSizeP(1,1)*0.5;
  ////////////////////////////////////////////////////////////////////////////////
     MvRecruitDev(DiCalBeginYear)=PbdvLogRecruitDevs(DiCalBeginYear);
     if(DvLfiftyParas_Ini(4)<=0&DvRsexParas_Ini(4)<=0)
     {
       dvR00=CalStockRecruitmentF(DmEnv(1,1),DmEnv(1,2),DmEnv(1,3),sdr_vSSB(1),dvAlpha,dvBeta,PbvEnvCoef,DiRSFlag)*mfexp(MvRecruitDev(DiCalBeginYear)); // MbnvR00Para(1)=MvNASTemp*DmFecundity(DiCalEndYear);
     }
     else
     {
       dvR00=CalStockRecruitmentF(DmEnv(1,1),DmEnv(1,2),DmEnv(1,3),sdr_vSSB2(1),dvAlpha,dvBeta,PbvEnvCoef,DiRSFlag)*mfexp(MvRecruitDev(DiCalBeginYear)); // MbnvR00Para(1)=MvNASTemp*DmFecundity(DiCalEndYear);
     }
  /////////////////////////////////////////////////////////////////////
  MvRecruitRS(DiCalBeginYear)=dvR00;
  //cout<<MvRecruitRS<<endl;
  ///////////////////////////////////////Population dynamics//////////////////////////////////////////////
  for(i=DiCalBeginYear;i<=DiCalEndYear;i++)
  {
     if(i>DiCalBeginYear)
     {//because first year R=0 so PbdvLogRecruitDevs begin DiCalBeginYear+1
          if(i<=DiCalEndYear-DiYearBeforeEndForRDev)
          {
              MvRecruitDev(i)=sqrt(PbnRecruitmentRh+0.0000000001)*MvRecruitDev(i-1)+sqrt(1.0000000001-PbnRecruitmentRh)*PbdvLogRecruitDevs(i);
              MbdvLogRecruitDevs(i)=PbdvLogRecruitDevs(i);
            // cout<<PbdvLogRecruitDevs(i)<<endl;
          }
          else
          {
              MvRecruitDev(i)=sqrt(PbnRecruitmentRh+0.0000000001)*MvRecruitDev(i-1)+sqrt(1.0000000001-PbnRecruitmentRh)*PbdvLogRecruitDevsADD(i);
              MbdvLogRecruitDevs(i)=PbdvLogRecruitDevsADD(i);
          }
          sdr_vRecruitmentP(i)=MvRecruitRS(i)*mfexp(MvRecruitDev(i));//-0.5*MnRSD*MnRSD); // Convergence is difficult for sd calculated from MvRecruitDev when MnRSD is included 
         //cout<<"R SD is over"<<endl;
     } 
     else
         sdr_vRecruitmentP(i)=MvRecruitRS(i);//*mfexp(PbdvLogRecruitDevs(i));
     for(j=1;j<=Flagtimestep;j++)
     {
             for (k=1;k<=DiSizeBinNum;k++)  
             {//Recruitment    
                   if(i>=DiCalBeginYear)   //MvNASTemp(i);
                        Md3NAS(j,i,k)= Md3NAS(j,i,k) +sdr_vRecruitmentP(i)*DvRecruitmentSeaRate(j)*MvRecruitPrjVect(k);
                   sdr_vNBiomass(i)+=Md3NAS(j,i,k)*DmWeightAtSize(i,k);
             }
             //cout<<DnFracYearSSB<<endl;
             //cout<<MmSSB_Y<<endl;
                sdr_vSSB(i)=elem_prod(Md3NAS(j,i),MmSSB_Y(i))*DmFecundity(i);   //sgwj debug
                sdr_vSSB2(i)=elem_prod(elem_prod(Md3NAS(j,i),MmSSB_Y(i)),DmWeightAtSize(i))*Md3FemalePropAtSizeP(1,i)*0.5; 
                if(DvLfiftyParas_Ini(4)<=0&DvRsexParas_Ini(4)<=0)
                {
                MvRecruitRS(i+1)=CalStockRecruitmentF(DmEnv(i,1),DmEnv(i,2),DmEnv(i,3),sdr_vSSB(i),dvAlpha,dvBeta,PbvEnvCoef,DiRSFlag);
                }
                else
                {
                MvRecruitRS(i+1)=CalStockRecruitmentF(DmEnv(i,1),DmEnv(i,2),DmEnv(i,3),sdr_vSSB2(i),dvAlpha,dvBeta,PbvEnvCoef,DiRSFlag);
                }
             //Die by Catch and natural mortality 
             MvNASTemp=elem_prod(Md3NAS(j,i),Md3S(j,i));//die        N(m+1,t)=(N(m,t)*exp(-Z)*G)
             iL1=DimGrowthMatrixBlockFlag(i,2);
             if(DiGrowthMatrixFlag==2)
             {
               cout<<"Growth matrix flag cannot be 2 for time step of year"<<endl;
               exit(1);
             }
             if(i+1<=DiCalEndYear)
                Md3NAS(j,i+1)= MvNASTemp*Md3GrowthMatrix(iL1);           //Growth Next Year
             else
                MvNASPJD     = MvNASTemp*Md3GrowthMatrix(iL1);                
       }
   }
   //////////////////////////////////////
 #if Debug_Status
    //ICHECKDEB(MvRecruitRS);
    //ICHECKDEB(sdr_vRecruitmentP);
    ICHECKDEB(Md3NAS);
    //ICHECKDEB(Md3NAS(2));
    //ICHECKDEB(Md3NAS(3));
    //ICHECKDEB(Md3NAS(4));
 #endif
  ///////////////////////End/////////////////////////////////////////////////////////////////////
  return;
}

void model_parameters::CalPredictedCatchAtSize(void)
{
  dvariable  dvTotal;
  Md4CatchPropAtSizeP.initialize();
  Md4CatchAtSizeP.initialize();
  Md3CatchTotalP.initialize();
  ///////////////////////////////////////////////////////////////////////////////
  for(i=1;i<=DiFleetNum;i++)
  {
    for(j=1;j<=DiSeasonNum;j++)
    {
        Md4CatchAtSizeP(i,j)=elem_prod(elem_div(Md4FASByFleet(i,j),Md3Z(j)),elem_prod(1.0-Md3S(j),Md3NAS(j)));
    }
  }
  //cout<<"Step 1"<<endl;
 //////////////////////////////////////////////////////////////////////////////////////////
  for (i=1;i<=DiFleetNum;i++)
  {
      for(j=1;j<=DiSeasonNum;j++)
      {
         for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
         {
            dvTotal=sum( Md4CatchAtSizeP(i,j,k)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j)));
            //cout<<"Step 3:"<<dvTotal<<endl;
            if ( dvTotal>0.0)
            {
               Md4CatchPropAtSizeP(i,j,k)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j))= Md4CatchAtSizeP(i,j,k)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j))/dvTotal;
               if(DimCatchBiomassFlag(i,j)==1)// in weight
                   Md3CatchTotalP(i,j,k)= Md4CatchAtSizeP(i,j,k)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j))*DmWeightAtSize(k)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j));
               else                 // in Number 
                   Md3CatchTotalP(i,j,k)= sum(Md4CatchAtSizeP(i,j,k)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j)));
            }
        }
      }
  }
 //cout<<Md4CatchPropAtSizeP<<endl;
 #if Debug_Status
   //ICHECKDEB(DmWeightAtSize)
   //ICHECKDEB(Md3CatchTotalP);
   //ICHECKDEB(Md4CatchAtSizeP(1,1));
   //ICHECKDEB(Md4CatchAtSizeP(1,2));
   //ICHECKDEB(Md4CatchAtSizeP(1,3));
   //ICHECKDEB(Md4CatchAtSizeP(1,4));
 #endif
   //cout<<"Step 2"<<endl;
}

void model_parameters::CalPredictedCatchAtSize_Y(void)
{
  dvariable  dvTotal;
  Md4CatchPropAtSizeP.initialize();
  Md4CatchAtSizeP.initialize();
  Md3CatchTotalP.initialize();
  ///////////////////////////////////////////////////////////////////////////////
  //cout<<Md4FASByFleet(2,1)<<endl;
  //cout<<Md3Z(1)<<endl;
  //cout<<Md3S(1)<<endl;
  //cout<<Md3NAS(1)<<endl;
  for(i=1;i<=DiFleetNum;i++)
  {
    Md4CatchAtSizeP(i,1)=elem_prod(elem_div(Md4FASByFleet(i,1),Md3Z(1)),elem_prod(1.0-Md3S(1),Md3NAS(1)));
  }
  //cout<<"Step 1"<<endl;
 //////////////////////////////////////////////////////////////////////////////////////////
  for (i=1;i<=DiFleetNum;i++)
  {
         for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
         {
            dvTotal=sum( Md4CatchAtSizeP(i,1,k)(DimCSelStartSizeBin(i,1),DimCSelEndSizeBin(i,1)));
            //cout<<"Step 3:"<<dvTotal<<endl;
            if ( dvTotal>0.0)
            {
               Md4CatchPropAtSizeP(i,1,k)(DimCSelStartSizeBin(i,1),DimCSelEndSizeBin(i,1))= Md4CatchAtSizeP(i,1,k)(DimCSelStartSizeBin(i,1),DimCSelEndSizeBin(i,1))/dvTotal;
               if(DimCatchBiomassFlag(i,1)==1)// in weight
                   Md3CatchTotalP(i,1,k)= Md4CatchAtSizeP(i,1,k)(DimCSelStartSizeBin(i,1),DimCSelEndSizeBin(i,1))*DmWeightAtSize(k)(DimCSelStartSizeBin(i,1),DimCSelEndSizeBin(i,1));
               else                 // in Number 
                   Md3CatchTotalP(i,1,k)= sum(Md4CatchAtSizeP(i,1,k)(DimCSelStartSizeBin(i,1),DimCSelEndSizeBin(i,1)));
            }
         }
  }
 //cout<<Md4CatchPropAtSizeP<<endl;
 #if Debug_Status
   //ICHECKDEB(DmWeightAtSize)
   //ICHECKDEB(Md3CatchTotalP);
   ICHECKDEB(Md4CatchAtSizeP);
   //ICHECKDEB(Md4CatchAtSizeP(1,2));
   //ICHECKDEB(Md4CatchAtSizeP(1,3));
   //ICHECKDEB(Md4CatchAtSizeP(1,4));
 #endif
   //cout<<"Step 2"<<endl;
}

void model_parameters::CalPredictedCatchAtSize2(void)
{
  dvariable  dvTotal;
  Md4CatchPropAtSizeP.initialize();
  Md4CatchAtSizeP.initialize();
  Md3CatchTotalP.initialize();
  ///////////////////////////////////////////////////////////////////////////////
  for(i=1;i<=DiFleetNum;i++)
  {
    for(j=1;j<=DiSeasonNum;j++)
    {
        Md4CatchAtSizeP(i,j)=elem_prod(elem_div(Md4FASByFleet(i,j),Md3Z(j)),elem_prod(1.0-Md3S(j),Md3NAS(j)));
    }
  }
 //////////////////////////////////////////////////////////////////////////////////////////
  for (i=1;i<=DiFleetNum;i++)
  {
      for(j=1;j<=DiSeasonNum;j++)
      {
         for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
         {
            if (DimCSelStartSizeBin(i,j)>DvSizeBinData(0))
            {
              dvTotal=sum( Md4CatchAtSizeP(i,j,k)(DimCSelStartSizeBin(i,j)-DvSizeBinData(0),DimCSelEndSizeBin(i,j)-DvSizeBinData(0)));
                if ( dvTotal>0.0)
                {
                    Md4CatchPropAtSizeP(i,j,k)(DimCSelStartSizeBin(i,j)-DvSizeBinData(0),DimCSelEndSizeBin(i,j)-DvSizeBinData(0))= Md4CatchAtSizeP(i,j,k)(DimCSelStartSizeBin(i,j)-DvSizeBinData(0),DimCSelEndSizeBin(i,j)-DvSizeBinData(0))/dvTotal;
                    if(DimCatchBiomassFlag(i,j)==1)// in weight
                        Md3CatchTotalP(i,j,k)= Md4CatchAtSizeP(i,j,k)(DimCSelStartSizeBin(i,j)-DvSizeBinData(0),DimCSelEndSizeBin(i,j))*DmWeightAtSize(k)(DimCSelStartSizeBin(i,j)-DvSizeBinData(0),DimCSelEndSizeBin(i,j)-DvSizeBinData(0));
                    else                 // in Number 
                        Md3CatchTotalP(i,j,k)= sum(Md4CatchAtSizeP(i,j,k)(DimCSelStartSizeBin(i,j)-DvSizeBinData(0),DimCSelEndSizeBin(i,j)-DvSizeBinData(0)));
                }
            }
            else
            {
              dvTotal=sum( Md4CatchAtSizeP(i,j,k));
                if ( dvTotal>0.0)
                {
                    Md4CatchPropAtSizeP(i,j,k)= Md4CatchAtSizeP(i,j,k)/dvTotal;
                    if(DimCatchBiomassFlag(i,j)==1)// in weight
                        Md3CatchTotalP(i,j,k)= Md4CatchAtSizeP(i,j,k)*DmWeightAtSize(k);
                    else                 // in Number 
                         Md3CatchTotalP(i,j,k)= sum(Md4CatchAtSizeP(i,j,k));
                }
            }
         }
      }
  }
  cout<<Md4CatchPropAtSizeP<<endl;
  //cout<<Md3CatchTotalP<<endl;
  //cout<<DimCSelStartSizeBin<<endl;
  //cout<<DimCSelEndSizeBin<<endl;
 #if Debug_Status
   ICHECKDEB(DmWeightAtSize)
   ICHECKDEB(Md3CatchTotalP);
   ICHECKDEB(Md4CatchAtSizeP);
 #endif
   //cout<<"Step 2"<<endl;
}

void model_parameters::CalPredCPUE(void)
{
  ivector ivCount(1,DiCPUEBlockNum);
  if(DiCPUEBlockNum<=0||DiAvailCPUENum<=0)
    return;
  ////////////////////////////////////////
  MmCPUEPred.initialize();//   Md3CPUEPred(1,DiAvailCPUENum,1,DiSeasonNum,DiCalBeginYear,DiCalEndYear)
  Md3NASBar.initialize();
  MmBiomassExploit.initialize();
  MvCPUEQ.initialize();
  //////////////////////////////////
  ivCount=0;
  for(i=1;i<=DiFleetNum;i++)
  {
    for(j=1;j<=DiSeasonNum;j++)
    {
          jj=(i-1)*DiSeasonNum+j;
          lL1=DivOldIndexAvailCPUE(jj);//
          if(lL1<1 ||lL1> DiAvailCPUENum)
              continue;
         for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
         {
            ii=Dd3CPUEBlockFlag(j,k,i+1);
            if(ii<1||ii>DiCPUEBlockNum)
                  continue;
            if(!DivCPUEBlockAvailFlag(ii))  //  i=DivCPUEQIndexFleet(ii);        // Real Fleet
                   continue;
             Md3NASBar(j,k)=elem_prod(elem_div(1.0-Md3S(j,k),Md3Z(j,k)),Md3NAS(j,k));//,mCatchWeightAtSize(j);
             MvTemp1(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j))=elem_prod(Md3NASBar(j,k)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j)),Md4FASByFleet(i,j,k)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j)));
             if(DimCatchBiomassFlag(i,j)==1)//debug sgwj 
                MmBiomassExploit(lL1,k)=MvTemp1(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j))*DmWeightAtSize(j)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j));
             else
                MmBiomassExploit(lL1,k)=sum(MvTemp1(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j)));
            switch(DivCPUEQChoiceFlag_Ini(ii))
            {
                case 1:
                     if(Dd3CPUEO(i,j,k)>MINPOSITIVE && MmBiomassExploit(lL1,k)> MINPOSITIVE)
                     {
                         if(active(PbnvCPUEQBiomassE1(ii)))
                              MvCPUEQ(ii)+=log(Dd3CPUEO(i,j,k)/pow(MmBiomassExploit(lL1,k),PbnvCPUEQBiomassE1(ii)));
                         else
                              MvCPUEQ(ii)+=log(Dd3CPUEO(i,j,k)/MmBiomassExploit(lL1,k));
                         ivCount(ii)=ivCount(ii)+1;
                     }
                     break;
                case 2:
                   if(Dd3CPUEO(i,j,k)>MINPOSITIVE && MmBiomassExploit(lL1,k)> MINPOSITIVE)
                   {
                       if(active(PbnvCPUEQBiomassE1(ii)))
                                 MvCPUEQ(ii)+=Dd3CPUEO(i,j,k)/pow(MmBiomassExploit(lL1,k),PbnvCPUEQBiomassE1(ii));
                       else
                                 MvCPUEQ(ii)+=Dd3CPUEO(i,j,k)/ MmBiomassExploit(lL1,k);
                       ivCount(ii)=ivCount(ii)+1;
                   }
                   break;
                case 3:
                   break;
                default:
                   break;
           }
      }
    }
  }
  for(ii=1;ii<=DiCPUEBlockNum;ii++)
  {
      switch(DivCPUEQChoiceFlag_Ini(ii))
      {
           case 1:
               if(ivCount(ii)>0)
                    MvCPUEQ(ii)=mfexp(MvCPUEQ(ii)/ivCount(ii));
               break;
           case 2:
               if(ivCount(ii)>0)
                   MvCPUEQ(ii)=MvCPUEQ(ii)/ivCount(ii);//  MvCPUEQ(ii)=log(MvCPUEQ(ii)/ivCount(ii));
               break;
           case 3:
               break;
          default:
              break;
       }
   }
   for(i=1;i<=DiFleetNum;i++)
   {
     for(j=1;j<=DiSeasonNum;j++)
     {
          jj=(i-1)*DiSeasonNum+j;
          lL1=DivOldIndexAvailCPUE(jj);//
          if(lL1<1 ||lL1> DiAvailCPUENum)
              continue;
          for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
          {
              ii=Dd3CPUEBlockFlag(j,k,i+1);
              if(!DivCPUEBlockAvailFlag(ii))  //  i=DivCPUEQIndexFleet(ii);        // Real Fleet
                   continue;
              if(Dd3CPUEO(i,j,k)>MINPOSITIVE && MmBiomassExploit(lL1,k)> MINPOSITIVE)
              {
                 if(active(PbnvCPUEQBiomassE1(ii)))
                        MmCPUEPred(lL1,k)=MvCPUEQ(ii)*pow(MmBiomassExploit(lL1,k),PbnvCPUEQBiomassE1(ii));
                 else  
                        MmCPUEPred(lL1,k)=MvCPUEQ(ii)*MmBiomassExploit(lL1,k);
              }
          }
      } 
   }
  return ;
}

void model_parameters::CalPredCPUE_Y(void)
{
  ivector ivCount(1,DiCPUEBlockNum);
  if(DiCPUEBlockNum<=0||DiAvailCPUENum<=0)
    return;
  ////////////////////////////////////////
  MmCPUEPred.initialize();//   Md3CPUEPred(1,DiAvailCPUENum,1,DiSeasonNum,DiCalBeginYear,DiCalEndYear)
  Md3NASBar.initialize();
  MmBiomassExploit.initialize();
  MvCPUEQ.initialize();
  //////////////////////////////////
  ivCount=0;
  for(i=1;i<=DiFleetNum;i++)
  {
    for(j=1;j<=DiSeasonNum;j++)
    {
          jj=(i-1)*DiSeasonNum+j;
          lL1=DivOldIndexAvailCPUE(jj);//
          if(lL1<1 ||lL1> DiAvailCPUENum)
              continue;
         for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
         {
            ii=Dd3CPUEBlockFlag(j,k,i+1);
            if(ii<1||ii>DiCPUEBlockNum)
                  continue;
            if(!DivCPUEBlockAvailFlag(ii))  //  i=DivCPUEQIndexFleet(ii);        // Real Fleet
                   continue;
             Md3NASBar(j,k)=elem_prod(elem_div(1.0-Md3S(j,k),Md3Z(j,k)),Md3NAS(j,k));//,mCatchWeightAtSize(j);
             MvTemp1(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j))=elem_prod(Md3NASBar(j,k)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j)),Md4FASByFleet(i,j,k)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j)));
             if(DimCatchBiomassFlag(i,j)==1)//debug sgwj 
                MmBiomassExploit(lL1,k)=MvTemp1(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j))*DmWeightAtSize(j)(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j));
             else
                MmBiomassExploit(lL1,k)=sum(MvTemp1(DimCSelStartSizeBin(i,j),DimCSelEndSizeBin(i,j)));
            switch(DivCPUEQChoiceFlag_Ini(ii))
            {
                case 1:
                     if(Dd3CPUEO(i,j,k)>MINPOSITIVE && MmBiomassExploit(lL1,k)> MINPOSITIVE)
                     {
                         if(active(PbnvCPUEQBiomassE1(ii)))
                              MvCPUEQ(ii)+=log(Dd3CPUEO(i,j,k)/pow(MmBiomassExploit(lL1,k),PbnvCPUEQBiomassE1(ii)));
                         else
                              MvCPUEQ(ii)+=log(Dd3CPUEO(i,j,k)/MmBiomassExploit(lL1,k));
                         ivCount(ii)=ivCount(ii)+1;
                     }
                     break;
                case 2:
                   if(Dd3CPUEO(i,j,k)>MINPOSITIVE && MmBiomassExploit(lL1,k)> MINPOSITIVE)
                   {
                       if(active(PbnvCPUEQBiomassE1(ii)))
                                 MvCPUEQ(ii)+=Dd3CPUEO(i,j,k)/pow(MmBiomassExploit(lL1,k),PbnvCPUEQBiomassE1(ii));
                       else
                                 MvCPUEQ(ii)+=Dd3CPUEO(i,j,k)/ MmBiomassExploit(lL1,k);
                       ivCount(ii)=ivCount(ii)+1;
                   }
                   break;
                case 3:
                   break;
                default:
                   break;
           }
      }
    }
  }
  for(ii=1;ii<=DiCPUEBlockNum;ii++)
  {
      switch(DivCPUEQChoiceFlag_Ini(ii))
      {
           case 1:
               if(ivCount(ii)>0)
                    MvCPUEQ(ii)=mfexp(MvCPUEQ(ii)/ivCount(ii));
               break;
           case 2:
               if(ivCount(ii)>0)
                   MvCPUEQ(ii)=MvCPUEQ(ii)/ivCount(ii);//  MvCPUEQ(ii)=log(MvCPUEQ(ii)/ivCount(ii));
               break;
           case 3:
               break;
          default:
              break;
       }
   }
   for(i=1;i<=DiFleetNum;i++)
   {
     for(j=1;j<=DiSeasonNum;j++)
     {
          jj=(i-1)*DiSeasonNum+j;
          lL1=DivOldIndexAvailCPUE(jj);//
          if(lL1<1 ||lL1> DiAvailCPUENum)
              continue;
          for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
          {
              ii=Dd3CPUEBlockFlag(j,k,i+1);
              if(!DivCPUEBlockAvailFlag(ii))  //  i=DivCPUEQIndexFleet(ii);        // Real Fleet
                   continue;
              if(Dd3CPUEO(i,j,k)>MINPOSITIVE && MmBiomassExploit(lL1,k)> MINPOSITIVE)
              {
                 if(active(PbnvCPUEQBiomassE1(ii)))
                        MmCPUEPred(lL1,k)=MvCPUEQ(ii)*pow(MmBiomassExploit(lL1,k),PbnvCPUEQBiomassE1(ii));
                 else  
                        MmCPUEPred(lL1,k)=MvCPUEQ(ii)*MmBiomassExploit(lL1,k);
              }
          }
      } 
   }
  return ;
}

void model_parameters::CalIndexSelectivity1(void)
{
  dvariable dvAlpha1;
  dvariable dvBelta1;
  dvariable dvAlpha2;
  dvariable dvBelta2;
  dvariable dvTemp1;
  dvariable dvTemp2;
  dvariable dvMax;
  // d3IndexSel
  for(i=1;i<=DiAvailIndexNum;i++)
  {
      for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
      {
           ii=DivAvailIndexSelIndexOld(i);       //old index
           jj=DimIndexSelBlockFlag(j,ii+1);      //block;
           if(jj<1||jj>DiIndexSelBlockNum)
                    continue;
           lL2=DivIndexSelParaIndex_Ini(jj);                 //get location of parameter
           if(DivIndexFleetSelChoiceFlag(i)<0)
           {
              switch(DivIndexSelFlag_Ini(jj))
              {
               case 1:
                  for(k=1;k<=DiSizeBinNum;k++)
                       Md3IndexSel(i,j,k)= PbnvIndexSelParams(lL2+k);
                  break;
               case 2:
                   dvAlpha1= PbnvIndexSelParams(lL2+1);
                   dvBelta1= PbnvIndexSelParams(lL2+2);//dvBelta1=1.0/ PbnvIndexSelParams(lL2+2);
                   //nTemp=
                   for(k=1;k<=DiSizeBinNum;k++)
                   {
                      nTemp=(DvSizeBinData(k)+DvSizeBinData(k-1))*0.5;
                      Md3IndexSel(i,j,k) =1.0/(1.0+mfexp((dvAlpha1-nTemp)* dvBelta1));
                   }
                   dvMax=max( Md3IndexSel(i,j));
                   Md3IndexSel(i,j)/= dvMax;
                   break;
               case 3:
                   dvAlpha1= PbnvIndexSelParams(lL2+1);
                   dvBelta1= PbnvIndexSelParams(lL2+2); //dvBelta1=1.0/ PbnvIndexSelParams(lL2+2);
                   dvAlpha2= PbnvIndexSelParams(lL2+3);
                   dvBelta2= PbnvIndexSelParams(lL2+4);//dvBelta2=1.0/ PbnvIndexSelParams(lL2+4);
                   for(k=1;k<=DiSizeBinNum;k++)
                   {
                      nTemp=(DvSizeBinData(k)+DvSizeBinData(k-1))*0.5;
                      dvTemp1 =1.0/(1.0+mfexp((dvAlpha1- nTemp)* dvBelta1));
                      dvTemp2 =1.0-1.0/(1.0+mfexp((dvAlpha2- nTemp)* dvBelta2));
                      Md3IndexSel(i,j,k)=dvTemp1*dvTemp2;
                   }
                   dvMax=max( Md3IndexSel(i,j));
                   Md3IndexSel(i,j)/= dvMax;
                   break;
           }
        }
        else
        {
              lL2=DivIndexFleetSelChoiceFlag(i);//get fleet
              /* if(DivIndexMonthFlag(i)<0)  //sgwj 2012-04-03
                 lL1=-1;
              else
              {
                   nTemp=0.0;
                   for(ii=1;ii<=DiSeasonNum;ii++)
                   {
                       nTemp=DvMonthInSeason(ii)+nTemp;
                       if(DivIndexMonthFlag(i)<nTemp)
                       {
                           lL1=ii;
                           break;
                       }
                   }
              }*/
              for(k=1;k<=DiSizeBinNum;k++)
              {
                 if(Dd3IndexMonthFlag(i,k,1)<0)
                        lL1=int(Dd3IndexMonthFlag(i,k,1)-0.001);
                 else
                        lL1=int(Dd3IndexMonthFlag(i,k,1)+0.001);
                 if(lL1<0)
                        Md3IndexSel(i,j,k)= 0.25*(Md4SelByFleet(lL2,1,j,k)+Md4SelByFleet(lL2,2,j,k)+Md4SelByFleet(lL2,3,j,k)+Md4SelByFleet(lL2,4,j,k));
                 else
                        Md3IndexSel(i,j,k)= Md4SelByFleet(lL2,lL1,j,k);
              }
         }
     }
  }
 #if Debug_Status
   //ICHECKDEB(Md3IndexSel);
 #endif
 // return;
}

void model_parameters::CalfemalePropAtsize(void)
{
  dvariable dvAlpha1;
  dvariable dvBelta1;
  int jj;
  for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
        {
          dvAlpha1=PbLfifty(j);
          dvBelta1=PbRsex;
          for(k=1;k<=DiSizeBinNum;k++)
             {
               nTemp=(DvSizeBinData(k)+DvSizeBinData(k-1))*0.5;
               Md3FemalePropAtSizeP(1,j,k)=1/(1+mfexp(-2*log(3)*(nTemp-dvAlpha1)/dvBelta1));
             }
        }
 #if Debug_Status
   ICHECKDEB(Md3FemalePropAtSizeP);
 #endif
}

void model_parameters::CalIndexSelectivity(void)
{
  dvariable dvAlpha1;
  dvariable dvBelta1;
  dvariable dvAlpha2;
  dvariable dvBelta2;
  dvariable dvTemp1;
  dvariable dvTemp2;
  dvariable dvMax;
  int jj;
  // d3IndexSel
  for(i=1;i<=DiAvailIndexNum;i++)
  {
    if(DivIndexFleetSelChoiceFlag(i)<0)
    {
        for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
        {
           ii=DivAvailIndexSelIndexOld(i);       //old index
           jj=DimIndexSelBlockFlag(j,ii+1);      //block;
           if(jj<1||jj>DiIndexSelBlockNum)
                    continue;
           lL2=DivIndexSelParaIndex_Ini(jj);                 //get location of parameter
           switch(DivIndexSelFlag_Ini(jj))
           {
               case 1:
                  for(k=1;k<=DiSizeBinNum;k++)
                       Md3IndexSel(i,j,k)= PbnvIndexSelParams(lL2+k);
                  break;
               case 2:
                   dvAlpha1= PbnvIndexSelParams(lL2+1);
                   dvBelta1= PbnvIndexSelParams(lL2+2);//dvBelta1=1.0/ PbnvIndexSelParams(lL2+2);
                   //nTemp=
                   for(k=1;k<=DiSizeBinNum;k++)
                   {
                      nTemp=(DvSizeBinData(k)+DvSizeBinData(k-1))*0.5;
                      Md3IndexSel(i,j,k) =1.0/(1.0+mfexp((dvAlpha1-nTemp)* dvBelta1));
                   }
                   dvMax=max( Md3IndexSel(i,j));
                   Md3IndexSel(i,j)/= dvMax;
                   break;
               case 3:
                   dvAlpha1= PbnvIndexSelParams(lL2+1);
                   dvBelta1= PbnvIndexSelParams(lL2+2); //dvBelta1=1.0/ PbnvIndexSelParams(lL2+2);
                   dvAlpha2= PbnvIndexSelParams(lL2+3);
                   dvBelta2= PbnvIndexSelParams(lL2+4);//dvBelta2=1.0/ PbnvIndexSelParams(lL2+4);
                   for(k=1;k<=DiSizeBinNum;k++)
                   {
                      nTemp=(DvSizeBinData(k)+DvSizeBinData(k-1))*0.5;
                      dvTemp1 =1.0/(1.0+mfexp((dvAlpha1- nTemp)* dvBelta1));
                      dvTemp2 =1.0-1.0/(1.0+mfexp((dvAlpha2- nTemp)* dvBelta2));
                      Md3IndexSel(i,j,k)=dvTemp1*dvTemp2;
                   }
                   dvMax=max( Md3IndexSel(i,j));
                   Md3IndexSel(i,j)/= dvMax;
                   break;
           }
        }
     }
     else
     {
           lL2=DivIndexFleetSelChoiceFlag(i);//get fleet
           for(jj=1;jj<=DivIndexObsNum(i);jj++)
           {
              j= DimIndexIndex(i,jj);//Year
              if(Dd3IndexMonthFlag(i,jj,1)<0)
                  lL1=int(Dd3IndexMonthFlag(i,jj,1)-0.001);
              else
                  lL1=int(Dd3IndexMonthFlag(i,jj,1)+0.001);
              for(k=1;k<=DiSizeBinNum;k++)
              {
                // cout<<"Fleet::"<<lL2<<"Season:"<<lL1<<endl;
                 if(lL1<0)
                        Md3IndexSel(i,j,k)= 0.25*(Md4SelByFleet(lL2,1,j,k)+Md4SelByFleet(lL2,2,j,k)+Md4SelByFleet(lL2,3,j,k)+Md4SelByFleet(lL2,4,j,k));
                 else
                        Md3IndexSel(i,j,k)= Md4SelByFleet(lL2,lL1,j,k);
              }
          }
     }
  }
 #if Debug_Status
   ICHECKDEB(Md3IndexSel);
 #endif
 // return;
}

void model_parameters::CalIndexQAndComp(void)
{
  ivector ivCount(1,DiIndexQBlockNum);
  if(DiIndexQBlockNum<=0||DiAvailIndexNum<=0)
    return;
  ivCount=0;
  ////////////////////////////////////////////////////
  MvIndexQ.initialize();
  Md3NASIndexBar.initialize();
  MmBiomassSurvey.initialize();
  Md3IndexPropAtSizeP.initialize();
  MmIndexP.initialize();
  /////////////////////////////////////////////////
  //d3IndexPropAtSizeP=0.0;
  for(i=1;i<=DiAvailIndexNum;i++)
  {
      ii=DivAvailIndexMirrorOld(i); // old survey index or real index by user input
      /////////////////////////////////////////////////////////////////
      for(j=1;j<=DivIndexObsNum(i);j++)
      {
          if(Dd3IndexMonthFlag(i,j,1)<0)
             lL1=int(Dd3IndexMonthFlag(i,j,1)-0.001);//sgwj 2012-04-03
          else
             lL1=int(Dd3IndexMonthFlag(i,j,1)+0.001);//sgwj 2012-04-03
          nTemp=Dd3IndexMonthFlag(i,j,2);//sgwj 2012-04-03
          k= DimIndexIndex(i,j);//Year
          if(lL1<0)
          {
            for(kk=1;kk<=DiSizeBinNum;kk++)
            {
              Md3NASIndexBar(i,k,kk) =(1.0-Md3S(1,k,kk))/Md3Z(1,k,kk)*Md3NAS(1,k,kk)*1.0/DiSeasonNum; 
              for(lL2=2;lL2<= DiSeasonNum;lL2++)
              {
                    Md3NASIndexBar(i,k,kk)+=(1.0-Md3S(lL2,k,kk))/Md3Z(lL2,k,kk)*Md3NAS(lL2,k,kk)*1.0/DiSeasonNum; 
              }
            }
          }
          else
              Md3NASIndexBar(i,k)=elem_prod(mfexp(nTemp*Md3Z(lL1,k)),Md3NAS(lL1,k));
         //////////////////////////////
          Md3NASIndexBar(i,k)=elem_prod(Md3NASIndexBar(i,k),Md3IndexSel(i,k));// NSurvey(k,t)
         //
          if(DivIndexUnitsFlag(i)==1)//
                MvTemp1=elem_prod( Md3NASIndexBar(i,k),DmWeightAtSize(k));
          else
                MvTemp1=Md3NASIndexBar(i,k);
         //
          MmBiomassSurvey(i,k)=sum(MvTemp1(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i)));//NSurvey(t)
          ////////////////////////////////////////////////////////////////////////////////////////////////////////
          iL1=DimIndexQBlockFlag(k,ii+1);//index old q
          //
          switch(DivIndexQChoiceFlag(iL1))
          {
           case 1:
                  if(active(PbnvIndexQBiomassE1(iL1)))//active(vQBiomassE1(i))
                         MvIndexQ(iL1)+=log(DmIndexTotalO(i,j)/pow( MmBiomassSurvey(i,k),PbnvIndexQBiomassE1(iL1)));
                  else
                         MvIndexQ(iL1)+=log(DmIndexTotalO(i,j)/ MmBiomassSurvey(i,k));
                 ivCount(iL1)++;
               break;
           case 2:
                 if(active(PbnvIndexQBiomassE1(iL1)))
                       MvIndexQ(iL1)+=DmIndexTotalO(i,j)/pow( MmBiomassSurvey(i,k),PbnvIndexQBiomassE1(iL1));
                 else
                       MvIndexQ(iL1)+=DmIndexTotalO(i,j)/ MmBiomassSurvey(i,k);
                 ivCount(iL1)++;
               break;
           default:
               break;
         }
      }
    }
    for(iL2=1;iL2<=DiIndexQBlockNum;iL2++)
    {
          switch(DivIndexQChoiceFlag(iL2))
         {
            case 1:
              if(ivCount(iL2)>0)
                    MvIndexQ(iL2)=mfexp(MvIndexQ(iL2)/ ivCount(iL2));
              else
                    MvIndexQ(iL2)=0.0;
              break;
            case 2:
                 if(ivCount(iL2)>0)
                    MvIndexQ(iL2)=MvIndexQ(iL2)/ ivCount(iL2);
              break;
            default:
              break;
        }     
    }
    for(i=1;i<=DiAvailIndexNum;i++)
    {
      ii=DivAvailIndexMirrorOld(i);
      for(j=1;j<=DivIndexObsNum(i);j++)
      {
         k= DimIndexIndex(i,j);//Year
         iL2=DimIndexQBlockFlag(k,ii+1);//index old q
         //
         if(active(PbnvIndexQBiomassE1(iL2)))
                MmIndexP(i,k) =  MvIndexQ(iL2)*pow(MmBiomassSurvey(i,k),PbnvIndexQBiomassE1(iL2));
         else
                MmIndexP(i,k)  = MvIndexQ(iL2)*MmBiomassSurvey(i,k);
         if(DivEstIndexPropFlag(i)!=1)
                     continue;
          //
         MnTemp1= sum(Md3NASIndexBar(i,k)(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i)));
         if(MnTemp1>0.0)
            Md3IndexPropAtSizeP(i,k)(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i))= Md3NASIndexBar(i,k)(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i))/MnTemp1;
      }
   }
   //cout<<DmIndexTotalO<<endl;
 #if Debug_Status
   ICHECKDEB(MmIndexP);
   ICHECKDEB(MvIndexQ);
   ICHECKDEB(Md3IndexPropAtSizeP);
 #endif
  return ;
}

void model_parameters::CalIndexQAndComp_Y(void)
{
  ivector ivCount(1,DiIndexQBlockNum);
  if(DiIndexQBlockNum<=0||DiAvailIndexNum<=0) //do nothing if the number of index blocks or the number of available index less than or equal to zero
    return;
  ivCount=0;
  MvIndexQ.initialize();
  Md3NASIndexBar.initialize();
  MmBiomassSurvey.initialize();
  Md3IndexPropAtSizeP.initialize();
  MmIndexP.initialize();
  for(i=1;i<=DiAvailIndexNum;i++)
    {
      ii=DivAvailIndexMirrorOld(i); // index number 1,2....
      for(j=1;j<=DivIndexObsNum(i);j++) // number of obs for a given index (e.g.,number of years if each year has an index)
        {
          if(Dd3IndexMonthFlag(i,j,1)<0) //season
             lL1=int(Dd3IndexMonthFlag(i,j,1)-0.001);
          else
             lL1=int(Dd3IndexMonthFlag(i,j,1)+0.001);
          temp=Dd3IndexMonthFlag(i,j,1);
          nTemp=Dd3IndexMonthFlag(i,j,2);//year fraction
          k= DimIndexIndex(i,j);//Year
          if(lL1<0)
          {
            for(kk=1;kk<=DiSizeBinNum;kk++)
            {
              Md3NASIndexBar(i,k,kk) =(1.0-Md3S(1,k,kk))/Md3Z(1,k,kk)*Md3NAS(1,k,kk); 
            }
          }
          else
              temp=(temp-1)*3/12+nTemp/3;
              Md3NASIndexBar(i,k)=elem_prod(mfexp(-temp*Md3Z(1,k)),Md3NAS(1,k));
              //cout<<temp<<endl;
          Md3NASIndexBar(i,k)=elem_prod(Md3NASIndexBar(i,k),Md3IndexSel(i,k));// NSurvey(k,t)
          if(DivIndexUnitsFlag(i)==1)//
                MvTemp1=elem_prod( Md3NASIndexBar(i,k),DmWeightAtSize(k));
          else
                MvTemp1=Md3NASIndexBar(i,k);
         //
          MmBiomassSurvey(i,k)=sum(MvTemp1(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i)));//NSurvey(t)
          ////////////////////////////////////////////////////////////////////////////////////////////////////////
          iL1=DimIndexQBlockFlag(k,ii+1);//index old q
          //cout<<DivIndexQChoiceFlag<<endl;
          switch(DivIndexQChoiceFlag(iL1))
          {
           case 1:
                  if(active(PbnvIndexQBiomassE1(iL1)))//active(vQBiomassE1(i))
                         MvIndexQ(iL1)+=log(DmIndexTotalO(i,j)/pow( MmBiomassSurvey(i,k),PbnvIndexQBiomassE1(iL1)));
                  else
                         MvIndexQ(iL1)+=log(DmIndexTotalO(i,j)/ MmBiomassSurvey(i,k));
                 ivCount(iL1)++;
               break;
           case 2:
                 if(active(PbnvIndexQBiomassE1(iL1)))
                       MvIndexQ(iL1)+=DmIndexTotalO(i,j)/pow( MmBiomassSurvey(i,k),PbnvIndexQBiomassE1(iL1));
                 else
                       MvIndexQ(iL1)+=DmIndexTotalO(i,j)/ MmBiomassSurvey(i,k);
                 ivCount(iL1)++;
               break;
           default:
               break;
         }
      }
    }
    for(iL2=1;iL2<=DiIndexQBlockNum;iL2++)
    {
          switch(DivIndexQChoiceFlag(iL2))
         {
            case 1:
              if(ivCount(iL2)>0)
                    MvIndexQ(iL2)=mfexp(MvIndexQ(iL2)/ ivCount(iL2));
              else
                    MvIndexQ(iL2)=0.0;
              break;
            case 2:
                 if(ivCount(iL2)>0)
                    MvIndexQ(iL2)=MvIndexQ(iL2)/ ivCount(iL2);
              break;
            default:
              break;
        }     
    }
    for(i=1;i<=DiAvailIndexNum;i++)
    {
      ii=DivAvailIndexMirrorOld(i);
      for(j=1;j<=DivIndexObsNum(i);j++)
      {
         k= DimIndexIndex(i,j);//Year
         iL2=DimIndexQBlockFlag(k,ii+1);//index old q
         //
         if(active(PbnvIndexQBiomassE1(iL2)))
                MmIndexP(i,k) =  MvIndexQ(iL2)*pow(MmBiomassSurvey(i,k),PbnvIndexQBiomassE1(iL2));
         else
                MmIndexP(i,k)  = MvIndexQ(iL2)*MmBiomassSurvey(i,k);
         if(DivEstIndexPropFlag(i)!=1)
                     continue;
          //
         MnTemp1= sum(Md3NASIndexBar(i,k)(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i)));
         if(MnTemp1>0.0)
         {
            Md3IndexPropAtSizeP(i,k)(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i))= Md3NASIndexBar(i,k)(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i))/MnTemp1;
            //Md3FemalePropAtSizeP(i,k)
         }
      }
   }
 #if Debug_Status
   ICHECKDEB(MmIndexP);
   ICHECKDEB(MvIndexQ);
   ICHECKDEB(Md3IndexPropAtSizeP);
   ICHECKDEB(MmBiomassSurvey);
 #endif
  return;
}

void model_parameters::CalObjectiveFunc(void)
{
      //CalLikelihoodVar(double obs,double dSigma2, const prevariable & pred, const int &iModel)
      ofvTotal=0.0;
      if(DiFleetSelParaNumIni>0)
          MvFleetSelLikely=0.0;
      MnFleetSelLikely=0.0;
      for( i=1;i<=DiFleetSelParaNumIni;i++)
      {//
            if(active(PbnvFleetSelParams(i)))
            {
                MvFleetSelLikely(i)=CalLikelihoodVar(DvFleetSelPara_Ini(i),DvFleetSelParaSigma2(i),PbnvFleetSelParams(i), DivFleetSelParaLikelihoodFlag(i))*DvFleetSelParaLambda(i);
                MnFleetSelLikely+= MvFleetSelLikely(i);
            }
      }
      MnFleetSelLikely+=DnFleetSelLikely;
      ofvTotal+=MnFleetSelLikely;
      //
      MvFYear1Likely=0.0;
      MnFYear1Likely=0.0;
      iTemp=DiFleetNum*DiSeasonNum;
      for(i=1;i<=iTemp;i++)
      {
              if(active(PbnvLogFYear1Season1(i)))
              {
                MvFYear1Likely(i)=CalLikelihoodVar(DvLogFYear1Ini(i),DvLogFYear1Sigma2(i), PbnvLogFYear1Season1(i),DivLogFYear1LikelihoodFlag(i))*DvLogFYear1Lambda(i);
                MnFYear1Likely+=MvFYear1Likely(i);
              }    
      }
      MnFYear1Likely+= DnFYear1Likely;
      ofvTotal+=MnFYear1Likely;
      //
      MnFLogDevsLikely=0.0;
      //MvFLogDevsLikely.initialize();
      for(i=1;i<=DiFleetNum;i++)
      {
           for(j=DiCalBeginYear;j<=DiCalEndYear;j++)  
           {
                 for(k=1;k<=DiSeasonNum;k++)
                 {
                       lL1=(i-1)*(DiCalEndYear-DiCalBeginYear+1)*DiSeasonNum+(j-DiCalBeginYear)*DiSeasonNum+k; 
                       lL2=(i-1)*DiSeasonNum+k;
                       if(active(PbnvLogFDevs(lL1)))//if(active(PbnvLogFDevs))//if(active(PbnvLogFDevs(lL1)))
                       {// DiFleetNum==>DiFleetNum*DiSeasonNum
                            // MvFLogDevsLikely(i)+=CalLikelihoodVar(0.0,DvFLogDevsSigma2(lL2),PbnvLogFDevs(lL1),3)*DvFLogDevsLambda(lL2); //CalLikelihoodConst(0.0,3)*DvFLogDevsLambda(i); 
                             MnFLogDevsLikely+=CalLikelihoodVar(0.0,DvFLogDevsSigma2(lL2),PbnvLogFDevs(lL1),3)*DvFLogDevsLambda(lL2); 
                       }
                 }
            }
            // MnFLogDevsLikely+=MvFLogDevsLikely(i);
     }
     MnFLogDevsLikely+=DnFLogDevsLikely;
     ofvTotal+=MnFLogDevsLikely;
      //
     MnCPUEQE1Likely=0.0;
     if(DiCPUEBlockNum>0)
     {
        MvCPUEQE1Likely=0.0;
     }
     for(i=1;i<=DiCPUEBlockNum;i++)
     {
                 if(active(PbnvCPUEQBiomassE1(i)))
                 {
                         MvCPUEQE1Likely(i)=CalLikelihoodVar(DvCPUEQBiomassE1Value(i),DvCPUEQBiomassE1Sigma2(i),PbnvCPUEQBiomassE1(i),DivCPUEQBiomassE1LikelihoodFlag(i))*DvCPUEQBiomassE1Lambda(i); // vector                     DvCPUEQE1Likely(1,DiAvailCPUEQNumP)
                         MnCPUEQE1Likely+=MvCPUEQE1Likely(i);// 
                 }
      }
      //cout<<"MnCPUEQE1Likely:"<<MnCPUEQE1Likely<<endl;
      MnCPUEQE1Likely+=DnCPUEQE1Likely;
      ofvTotal+=MnCPUEQE1Likely;
      //
      MnCatchLikely=0.0;
      MnCPUELikely=0.0;
      if(DiAvailCPUENum>0)
          MvCPUELikely=0.0;
      MvCatchTotalLikely=0.0;
      for(i=1;i<=DiFleetNum;i++)
      {
            for(j=1;j<=DiSeasonNum;j++)
            {
              for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
              {
                 if(Dd3CatchTotalO(i,j,k)>MINPOSITIVE)
                 {
                     MvCatchTotalLikely(i)+= CalLikelihoodVar(Dd3CatchTotalO(i,j,k),Dd3CatchTotalSigma2(i,j,k),Md3CatchTotalP(i,j,k),DimCatchTotalLikelihoodFlag(i,j))*DmCatchTotalLambda(i,j);  
                 //cout<<Dd3CatchTotalO(i,j,k)<<endl;
                 //cout<<Md3CatchTotalP(i,j,k)<<endl;
                 //cout<<MvCatchTotalLikely(i)<<endl;
                 //cout<<DimCatchTotalLikelihoodFlag(i,j)<<endl;
                 //cout<<Dd3CatchTotalSigma2(i,j,k)<<endl;
                 //cout<<DmCatchTotalLambda(i,j)<<endl;
                 }
              }
              iL1=(i-1)*DiSeasonNum+j;
              ii=DivOldIndexAvailCPUE(iL1);
              if(ii<=0||ii>DiAvailCPUENum)
                         continue;
              for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
              {
                 if(DimAvailCPUEFlag(i,j) &&  Dd3CPUEO(i,j,k)>MINPOSITIVE)
                 {    
                     MvCPUELikely(ii)+= CalLikelihoodVar(Dd3CPUEO(i,j,k),Dd3CPUESigma2(i,j,k),MmCPUEPred(ii,k),DimCPUELikelihoodFlag(i,j))*DmCPUELambda(i,j);
                 }
              }
              MnCPUELikely+=MvCPUELikely(ii);         
          }
           MnCatchLikely +=MvCatchTotalLikely(i);
      }
      cout<<"Total catch likelihood:"<<MvCatchTotalLikely<<endl;
      //MnCatchLikely+= DnCatchLikely;
     // ofvTotal+= MnCatchLikely;
      MnCPUELikely  +=DnCPUELikely;
      //cout<<"MnCatchLikely 1.5:"<<MnCatchLikely<<endl;
      ofvTotal+= MnCPUELikely;
      //if(DiCatchCompLikelihoodFlag==1)  //Multinomial distribution 
      MvCatchPropLikely=0.0;
      //MnCatchLikely=0.0;
      for(i=1;i<=DiFleetNum;i++)
      {
              for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
              {
                       for(k=1;k<=DiSeasonNum;k++)
                      {
                            iTemp=int(Dd3CatchESSInput(i,k,j)+0.00001);
                            //cout<<j<<" "<<k <<" "<<iTemp<<" ";
                            if(iTemp>0 && Dd3CatchTotalO(i,k,j)>MINPOSITIVE)
                            {
                                 iL1=DimCSelStartSizeBin(i,k);
                                 iL2=DimCSelEndSizeBin(i,k);
                                 DvTemp=Dd4CatchPropAtSizeO(i,k,j);
                                 MvTemp1=Md4CatchPropAtSizeP(i,k,j);
                                 //cout<< MvTemp1<<endl;
                                 MvCatchPropLikely(i)+=CalLikelihoodPropVar(DvTemp,MvTemp1,iTemp,iL1,iL2,DimCatchCompLikelihoodFlag(i,k))*DmCatchCompLambda(i,k);
                                 //cout<<MvCatchPropLikely(i)<<endl;
                                 //exit(-99);
                           }
                       }
              }
              MnCatchLikely+=MvCatchPropLikely(i);
      }
      cout<<"Catch composition likelihood:"<<MvCatchPropLikely<<endl;
     // cout<<Dd4CatchPropAtSizeO<<endl;
      MnCatchLikely += DnCatchLikely;
      //cout<<"MnCatchLikely 2.5:"<<MnCatchLikely<<endl;
      ofvTotal+=MnCatchLikely;
      //exit(-99);
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     if(DiIndexSelParaNumIni>0)
         MvIndexSelLikely=0.0;
     MnIndexSelLikely=0.0;
     for(i=1;i<=DiIndexSelParaNumIni;i++)
     {
             //cout<<PbnvIndexSelParams(i)<<endl;
             if(active(PbnvIndexSelParams(i))) //DivIndexSelParaPh(i)>0)
             {
                MvIndexSelLikely(i)=CalLikelihoodVar(DvIndexSelPara_Ini(i),DvIndexSelParaSigma2(i),PbnvIndexSelParams(i),DivIndexSelParaLikelihoodFlag(i))*DvIndexSelParaLambda(i);
                MnIndexSelLikely+= MvIndexSelLikely(i);
             }
     }
     MnIndexSelLikely+= DnIndexSelLikely;
     ofvTotal+=MnIndexSelLikely;
      //
     MnIndexLikely=0.0;
     if(DiAvailIndexNum>0)
           MvIndexLikely=0.0;
     //cout<<DiAvailIndexNum<<endl;
     //cout<<DmIndexTotalSigma2(1,1)<<endl;
     iL1=0;
     for(i=1;i<=DiAvailIndexNum;i++)
     {
          for(k=1;k<=DivIndexObsNum(i);k++)
          {
                     iL1++;
                     kk= DimIndexIndex(i,k);//Year
                     MvIndexLikely(i)+= CalLikelihoodVar(DmIndexTotalO(i,k),DmIndexTotalSigma2(i,k),MmIndexP(i,kk),DivIndexTotalLikelihoodFlag(i))*DvIndexTotalLambda(i);      
          }
           MnIndexLikely +=MvIndexLikely(i);
     }
     //cout<<iL1<<endl;
     //cout<<MmIndexP<<endl;
     //cout<<DmIndexTotalO<<endl;
     //cout<<DivIndexObsNum<<endl;
     cout<<"Survey index likelihood:"<<MvIndexLikely<<endl;
      // 
     if(DiAvailIndexNum>0)
           MvIndexPropLikely=0.0;
     for(i=1;i<=DiAvailIndexNum;i++)
     {
               if(DivEstIndexPropFlag(i)!=1)
                    continue;
               for(k=1;k<=DivIndexObsNum(i);k++)
               {
                      kk= DimIndexIndex(i,k);//Year
                      iTemp=DimIndexESSInput(i,k);
                      if(iTemp>0)
                      {
                                 iL1=DivIndexSelStartSizeBin(i);
                                 iL2=DivIndexSelEndSizeBin(i);
                                 DvTemp=Dd3IndexPropAtSizeO(i,k);
                                 MvTemp1=Md3IndexPropAtSizeP(i,kk);
                                 MvIndexPropLikely(i)+=CalLikelihoodPropVar(DvTemp,MvTemp1,iTemp,iL1,iL2,DivIndexCompLikelihoodFlag(i))*DvIndexCompLambda(i);
                               //  cout<<MvIndexPropLikely(i)<<endl;
                       }
              }
              MnIndexLikely+= MvIndexPropLikely(i);
      }
      cout<<"Survey composition likelihood:"<<MvIndexPropLikely<<endl;
      MnIndexLikely +=DnIndexLikely ;
      //cout<<"MnIndexLikely 2.5:"<<MnIndexLikely<<endl;
      ofvTotal+= MnIndexLikely;
      //3/29/2013
      MnLfiftylikely=0.0;
      for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
      {
          iTemp=DimIndexESSInput(DiIndexFlagSexComp,j);
          //cout<<iTemp<<endl;
                      if(iTemp>0)
                      {
                                 DvTemp=DmSexAtSizeData_Ini(j);
                                 MvTemp1=Md3FemalePropAtSizeP(1,j);
                                 MnLfiftylikely+=CalLikelihoodBinomial(DvTemp,MvTemp1,iTemp)*DiSexAtSizeLamda;
                               //  cout<<MvIndexPropLikely(i)<<endl;
                       }
       }
      cout<<MnLfiftylikely<<endl;
      ofvTotal+= MnLfiftylikely;
      //exit(-99);
        //
      MnIndexQE1Likely=0.0;
      if(DiIndexQBlockNum>0)
             MvIndexQE1Likely=0.0;
      for(i=1;i<=DiIndexQBlockNum;i++)
      {
            if(active(PbnvIndexQBiomassE1(i)))
            {
              MvIndexQE1Likely(i)=CalLikelihoodVar( DvQBiomassE1Value(i),DvQBiomassE1Sigma2(i),PbnvIndexQBiomassE1(i),DivQBiomassE1LikelihoodFlag(i))*DvQBiomassE1Lambda(i);
              MnIndexQE1Likely+= MvIndexQE1Likely(i);
            }
       }
       MnIndexQE1Likely+=DnIndexQE1Likely;
       ofvTotal+=MnIndexQE1Likely;
        //
       MvNYear1Likely=0.0;
       MnNYear1Likely=0.0;
       for(i=1;i<=DiNYear1EstParaNum;i++)
       {
                if(active(PbvNYear1Para(i)))
                {  
                  MvNYear1Likely(i)=CalLikelihoodVar(DvNYear1ParaInital(i),DvNYear1ParaSigma2(i),PbvNYear1Para(i),DivNYear1ParaLikelihoodFlag(i))*DvNYear1ParaLambda(i);
                  MnNYear1Likely+=MvNYear1Likely(i);   
                }
        }
        MnNYear1Likely+=DnNYear1Likely;
        ofvTotal+=MnNYear1Likely;
        //
        MnRSParaLikely=0.0;
        MvRSParaLikely=0.0;
        //
        if(DiRSFlag==1) 
           lL1=1; 
        else 
           lL1=2;
        for(i=1;i<=lL1;i++)
        {
                  if(active(PbnvLogRSPara(i)))
                  {
                      MvRSParaLikely(i)= CalLikelihoodVar( DvRSParaValue(i),DvRSParaSigma2(i),PbnvLogRSPara(i),DivRSParaLikelihoodFlag(i))*DvRSParaLambda(i);
                      MnRSParaLikely+=MvRSParaLikely(i);
                  }
       }
       MnRSParaLikely += DnRSParaLikely;   
       ofvTotal+=MnRSParaLikely ;
       //
       MnRDevsSDLikely=0.0;
       if(DiRecruitLogDevsSDEstFlag)
       {// lL1=DvRecruitLogDevsSD_Ini(4)>0?int(DvRecruitLogDevsSD_Ini(4)+0.01):int(DvRecruitLogDevsSD_Ini(4)-0.01);
                  if(active(PbnvRecruitLogDevsSD2))
                  {
                    if(DvRecruitLogDevsSD_Ini(5)>MINPOSITIVE)
                         nTemp=DvRecruitLogDevsSD_Ini(5)*DvRecruitLogDevsSD_Ini(5);//log(DvRecruitLogDevsSD_Ini(5)*DvRecruitLogDevsSD_Ini(5)+1.0);
                    else
                         nTemp=SIGMA2FILL;
                    MnRDevsSDLikely=CalLikelihoodVar(DvRecruitLogDevsSD_Ini(1),nTemp,PbnvRecruitLogDevsSD2,int(DvRecruitLogDevsSD_Ini(7)+0.01))*DvRecruitLogDevsSD_Ini(6);
                  }
       }
       MnRDevsSDLikely += DnRDevsSDLikely;
       ofvTotal+=MnRDevsSDLikely;//
       //
       MnRecruitDevsLikely=0.0;
       if((active(PbdvLogRecruitDevs)&DiFlagEvntoRdevs)!=1)
       {
                  if(!DiRecruitLogDevsSDEstFlag)
                        nTemp=DnRecruitLogDevsSigma;    
                  //else if(active(PbnvRecruitLogDevsSD2))
                  //      nTemp=value(PbnvRecruitLogDevsSD2);
                  else
                        nTemp=value(MnRSD);
                  nTemp=nTemp*nTemp;
                  for(k=DiCalBeginYear;k<=DiCalEndYear-DiYearBeforeEndForRDev;k++)
                  {//DiCalEndYear-DiYearBeforeEndForRDev
                        if(DvRecruitLogDevsSD_Ini(4)>0)
                             MnRecruitDevsLikely += CalLikelihoodVar(0.0,square(PbnvRecruitLogDevsSD2),PbdvLogRecruitDevs(k),3)*DnRecruitLogDevsLambda;
                        else
                             MnRecruitDevsLikely += CalLikelihoodVar(0.0,nTemp,PbdvLogRecruitDevs(k),3)*DnRecruitLogDevsLambda; 
                         //cout<<MnRSD<<" "<<MnRecruitDevsLikely<<" "<<PbdvLogRecruitDevs(k)<<endl;
                  }
                  for(k=DiCalEndYear-DiYearBeforeEndForRDev+1;k<=DiCalEndYear;k++)
                  {//DiCalEndYear-DiYearBeforeEndForRDev
                        if(DvRecruitLogDevsSD_Ini(4)>0)
                            MnRecruitDevsLikely += CalLikelihoodVar(0.0,square(PbnvRecruitLogDevsSD2),PbdvLogRecruitDevsADD(k),3)*DnRecruitLogDevsLambda;
                        else
                            MnRecruitDevsLikely += CalLikelihoodVar(0.0,nTemp,PbdvLogRecruitDevsADD(k),3)*DnRecruitLogDevsLambda;
                  }
                  //cout<<PbnvRecruitLogDevsSD2<<endl;
                  //cout<<PbdvLogRecruitDevsADD<<endl;
                  //cout<<"nTemp:"<<nTemp<<" "<<"DnRecruitLogDevsLambda"<<DnRecruitLogDevsLambda<<endl;
                  //cout<<"MnRecruitDevsLikely:"<<MnRecruitDevsLikely<<endl;
        }
        MnRecruitDevsLikely+=DnRecruitDevsLikely;
        cout<<"Recruitment likelihood:"<<MnRecruitDevsLikely<<endl;
        ofvTotal+=MnRecruitDevsLikely;//
        //exit(-99);
        MnEnvfitRdevLikely=0;
        if(DiFlagEvntoRdevs==1)
        {
                 for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
                          {
                                 MnEnvfitRdevLikely += CalLikelihoodVar(mfexp(DvEnvfitZscore(k)),StdevEnvZ,mfexp(PbdvLogRecruitDevs(k)),5);
                                 //cout<<MnEnvfitRdevLikely<<endl;
                          }
        //cout<<DvEnvfitZscore<<endl;
        //cout<<PbdvLogRecruitDevs<<endl;
        }       
        ofvTotal+=MnEnvfitRdevLikely;
        cout<<"Environment likelihood:"<<MnEnvfitRdevLikely<<endl;
        //exit(-99);
        MnRecruitRhLikely=0.0;
        if(active(PbnRecruitmentRh))
        {
          if(DvRecruitmentRh_Ini(5)>MINPOSITIVE)
                         nTemp=DvRecruitmentRh_Ini(5)*DvRecruitmentRh_Ini(5);//log(DvRecruitmentRh_Ini(5)*DvRecruitmentRh_Ini(5)+1.0);
          else
                         nTemp=SIGMA2FILL;
          MnRecruitRhLikely=CalLikelihoodVar(DvRecruitmentRh_Ini(1),nTemp,PbnRecruitmentRh,int(DvRecruitmentRh_Ini(7)+0.01))*DvRecruitmentRh_Ini(6);       
        }
        ofvTotal+= MnRecruitRhLikely;//
       //
        MnLorenzenALikely=0.0;
        if(!DiNaturalMortalityFlag)
        {   
                 MvLorenzenALikely=0.0;
                 for(i=1;i<=DiSeasonNum;i++)
                 {
                       if(active(PbnvLorenzenA(i)))
                       {
                         MvLorenzenALikely(i)=CalLikelihoodVar(DvLorenzenAValue(i),DvLorenzenASigma2(i),PbnvLorenzenA(i),DivLorenzenALikelihoodFlag(i))*DvLorenzenALambda(i);
                         MnLorenzenALikely+=MvLorenzenALikely(i);
                        }
                  }
                  MnLorenzenBLikely=0.0;
                  for(i=1;i<=DiSeasonNum;i++)
                  {
                         if(active(PbnvLorenzenB(i)))
                         {
                            MvLorenzenBLikely(i)=CalLikelihoodVar(DvLorenzenBValue(i),DvLorenzenBSigma2(i),PbnvLorenzenB(i),DivLorenzenBLikelihoodFlag(i))*DvLorenzenBLambda(i);
                            MnLorenzenBLikely+=MvLorenzenBLikely(i);
                          }
                  }
        }
        MnLorenzenALikely+=DnLorenzenALikely;
        ofvTotal+=MnLorenzenALikely;
       MnLorenzenBLikely+=DnLorenzenBLikely;
       ofvTotal+= MnLorenzenBLikely;        
        MnGrowthMatrixLikely=0.0;
        if(DiGrowthMatrixFlag)
        {
                MvLinfLikely=0.0;
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(active(PbnvLinf(i)))
                     {
                       MvLinfLikely(i)= CalLikelihoodVar(DvVBGFLinfValue(i),DvVBGFLinfSigma2(i),PbnvLinf(i),DivVBGFLinfLikelihoodFlag(i))*DvVBGFLinfLambda(i);
                       MnGrowthMatrixLikely+=MvLinfLikely(i);
                     }
                }
                MvVLikely=0.0;
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(active(PbnvV(i)))
                     {
                       MvVLikely(i)= CalLikelihoodVar(DvVBGFVValue(i),DvVBGFVSigma2(i),PbnvV(i),DivVBGFVLikelihoodFlag(i))*DvVBGFVLambda(i);
                       MnGrowthMatrixLikely+=MvVLikely(i);
                     }
                }
                MvLSDLikely=0.0;
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(active(PbnvLSD(i)))
                     {
                       MvLSDLikely(i)= CalLikelihoodVar(DvVBGFLSDValue(i),DvVBGFLSDSigma2(i),PbnvLSD(i),DivVBGFLSDLikelihoodFlag(i))*DvVBGFLSDLambda(i);
                       MnGrowthMatrixLikely+=MvLSDLikely(i);
                     }
                }
                MvVSDLikely=0.0;
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(active(PbnvVSD(i)))
                     {
                       MvVSDLikely(i)= CalLikelihoodVar(DvVBGFVSDValue(i),DvVBGFVSDSigma2(i),PbnvVSD(i),DivVBGFVSDLikelihoodFlag(i))*DvVBGFVSDLambda(i);
                       MnGrowthMatrixLikely+=MvVSDLikely(i);
                     }
                }
                MvLVRhoLikely=0.0;
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(active(PbnvLVRho(i)))
                     {
                       MvLVRhoLikely(i)= CalLikelihoodVar(DvVBGFLVRhoValue(i),DvVBGFLVRhoSigma2(i),PbnvLVRho(i),DivVBGFLVRhoLikelihoodFlag(i))*DvVBGFLVRhoLambda(i);
                       MnGrowthMatrixLikely+=MvLVRhoLikely(i);
                     }
                }
    }       
    MnGrowthMatrixLikely +=DnGrowthMatrixLikely;  
    ofvTotal+= MnGrowthMatrixLikely;
    //FMax
    MnFMaxPenaltyLikely=0.0;
    for (i=1;i<=DiFleetNum;i++)
    {
       for(k=1;k<=DiSeasonNum;k++)
       {
          for (j=DiCalBeginYear;j<=DiCalEndYear;j++)
          {
             MnTemp1=Md3Ft(i,k,j);
             if(MnTemp1>DnFmax)
               MnFMaxPenaltyLikely+=1000.*(MnTemp1-DnFmax)*(MnTemp1-DnFmax);
          }
      }
    }
    if(!last_phase())
    {
         if(DnRecruitLogDevsLambda<0.000000001)
          {
                     DnRecruitLogDevsLambda=0.05;
                     DiPenaltyLambdaFlag=1;
          }
         /*for(i=DiCalBeginYear;i<=DiCalEndYear;i++)
         {
              if(sdr_vRecruitmentP(i)<0.01*MvRecruitRS(i)||sdr_vRecruitmentP(i)>100*MvRecruitRS(i))
              {
                  //cout<<"Penalty:"<<i<<" "<<sdr_vRecruitmentP(i)<<" "<<MvRecruitRS(i)<<endl;
                  MnTemp1=sdr_vRecruitmentP(i)-MvRecruitRS(i);//sum(PbdvLogRecruitDevs(DiCalBeginYear+1,DiCalEndYear-DiYearBeforeEndForRDev-1));
                  MnFMaxPenaltyLikely+=MnTemp1*MnTemp1*0.1;
              }
         }//*/
    }
    else
    {
            if(DiPenaltyLambdaFlag)
            {
                 DnRecruitLogDevsLambda=0.0;
                 DiPenaltyLambdaFlag=0;
            }
           /* for(i=DiCalBeginYear+1;i<=DiCalEndYear;i++)
            {
                 if(sdr_vRecruitmentP(i)<0.05*sdr_vRecruitmentP(1))
                 {
                     MnTemp1=sdr_vRecruitmentP(i)-sdr_vRecruitmentP(1);
                     MnFMaxPenaltyLikely+=MnTemp1*MnTemp1*0.1;
                 }
           }*/
    }
    cout<<"Penalty:"<<MnFMaxPenaltyLikely<<endl;
    //cout<<"PbdvLogRecruitDevs:"<<PbdvLogRecruitDevs<<endl;
    ofvTotal+= MnFMaxPenaltyLikely;
    MnRPVLikely=0.0;// PbnvRecruitPrjVect(i)=DbnvRecruitPrjVect(i);
    for(i=1;i<=DiSizeBinNum;i++)
    {
        if(active(PbnvRecruitPrjVect(i)))
        {
           MnLorenzenALikely+=CalLikelihoodVar(DbnvRecruitPrjVect(i),DnRPVSigma2,PbnvRecruitPrjVect(i),DiRPVLikelihoodFlag)*DnRPVLambda;
        }
    }
    //cout<<"MnRPVLikely:"<<MnRPVLikely<<endl;
    /////////////////////////////////////////////////////////
    ofvTotal+=MnRPVLikely;
    cout<<"Objective Function Value:"<< ofvTotal<<endl;
    //exit(70);
}

void model_parameters::CalObjectiveFunc_Year(void)
{
      //CalLikelihoodVar(double obs,double dSigma2, const prevariable & pred, const int &iModel)
      ofvTotal=0.0;
      if(DiFleetSelParaNumIni>0)
          MvFleetSelLikely=0.0;
      MnFleetSelLikely=0.0;
      for( i=1;i<=DiFleetSelParaNumIni;i++)
      {//
            if(active(PbnvFleetSelParams(i)))
            {
                MvFleetSelLikely(i)=CalLikelihoodVar(DvFleetSelPara_Ini(i),DvFleetSelParaSigma2(i),PbnvFleetSelParams(i), DivFleetSelParaLikelihoodFlag(i))*DvFleetSelParaLambda(i);
                MnFleetSelLikely+= MvFleetSelLikely(i);
            }
      }
      MnFleetSelLikely+=DnFleetSelLikely;
      ofvTotal+=MnFleetSelLikely;
      //
      MvFYear1Likely=0.0;
      MnFYear1Likely=0.0;
      iTemp=DiFleetNum*Flagtimestep;
      for(i=1;i<=iTemp;i++)
      {
              if(active(PbnvLogFYear1Season1(i)))
              {
                MvFYear1Likely(i)=CalLikelihoodVar(DvLogFYear1Ini(i),DvLogFYear1Sigma2(i), PbnvLogFYear1Season1(i),DivLogFYear1LikelihoodFlag(i))*DvLogFYear1Lambda(i);
                MnFYear1Likely+=MvFYear1Likely(i);
              }    
      }
      MnFYear1Likely+= DnFYear1Likely;
      ofvTotal+=MnFYear1Likely;
      //
      MnFLogDevsLikely=0.0;
      //MvFLogDevsLikely.initialize();
      for(i=1;i<=DiFleetNum;i++)
      {
           for(j=DiCalBeginYear;j<=DiCalEndYear;j++)  
           {
                 for(k=1;k<=Flagtimestep;k++)
                 {
                       lL1=(i-1)*(DiCalEndYear-DiCalBeginYear+1)*Flagtimestep+j; 
                       lL2=(i-1)*Flagtimestep+k;
                       if(active(PbnvLogFDevs(lL1)))//if(active(PbnvLogFDevs))//if(active(PbnvLogFDevs(lL1)))
                       {// DiFleetNum==>DiFleetNum*DiSeasonNum
                            // MvFLogDevsLikely(i)+=CalLikelihoodVar(0.0,DvFLogDevsSigma2(lL2),PbnvLogFDevs(lL1),3)*DvFLogDevsLambda(lL2); //CalLikelihoodConst(0.0,3)*DvFLogDevsLambda(i); 
                             MnFLogDevsLikely+=CalLikelihoodVar(0.0,DvFLogDevsSigma2(lL2),PbnvLogFDevs(lL1),3)*DvFLogDevsLambda(lL2); 
                       }
                 }
            }
            // MnFLogDevsLikely+=MvFLogDevsLikely(i);
     }
     MnFLogDevsLikely+=DnFLogDevsLikely;
     ofvTotal+=MnFLogDevsLikely;
      //
     MnCPUEQE1Likely=0.0;
     if(DiCPUEBlockNum>0)
     {
        MvCPUEQE1Likely=0.0;
     }
     for(i=1;i<=DiCPUEBlockNum;i++)
     {
                 if(active(PbnvCPUEQBiomassE1(i)))
                 {
                         MvCPUEQE1Likely(i)=CalLikelihoodVar(DvCPUEQBiomassE1Value(i),DvCPUEQBiomassE1Sigma2(i),PbnvCPUEQBiomassE1(i),DivCPUEQBiomassE1LikelihoodFlag(i))*DvCPUEQBiomassE1Lambda(i); // vector                     DvCPUEQE1Likely(1,DiAvailCPUEQNumP)
                         MnCPUEQE1Likely+=MvCPUEQE1Likely(i);// 
                 }
      }
      //cout<<"MnCPUEQE1Likely:"<<MnCPUEQE1Likely<<endl;
      MnCPUEQE1Likely+=DnCPUEQE1Likely;
      ofvTotal+=MnCPUEQE1Likely;
      //
      MnCatchLikely=0.0;
      MnCPUELikely=0.0;
      if(DiAvailCPUENum>0)
          MvCPUELikely=0.0;
      MvCatchTotalLikely=0.0;
      for(i=1;i<=DiFleetNum;i++)
      {
            for(j=1;j<=Flagtimestep;j++)
            {
              for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
              {
                 if(Dd3CatchTotalO(i,j,k)>MINPOSITIVE)
                 {
                     MvCatchTotalLikely(i)+= CalLikelihoodVar(Dd3CatchTotalO(i,j,k),Dd3CatchTotalSigma2(i,j,k),Md3CatchTotalP(i,j,k),DimCatchTotalLikelihoodFlag(i,j))*DmCatchTotalLambda(i,j);  
                 }
              }
              iL1=i;
              ii=DivOldIndexAvailCPUE(iL1);
              if(ii<=0||ii>DiAvailCPUENum)
                         continue;
              for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
              {
                 if(DimAvailCPUEFlag(i,j) &&  Dd3CPUEO(i,j,k)>MINPOSITIVE)
                 {    
                     MvCPUELikely(ii)+= CalLikelihoodVar(Dd3CPUEO(i,j,k),Dd3CPUESigma2(i,j,k),MmCPUEPred(ii,k),DimCPUELikelihoodFlag(i,j))*DmCPUELambda(i,j);
                 }
              }
              MnCPUELikely+=MvCPUELikely(ii);         
          }
           MnCatchLikely +=MvCatchTotalLikely(i);
      }
      //cout<<"Total catch likelihood:"<<MnCatchLikely<<endl;
      //MnCatchLikely+= DnCatchLikely;
     // ofvTotal+= MnCatchLikely;
      MnCPUELikely  +=DnCPUELikely;
      //cout<<"MnCatchLikely 1.5:"<<MnCatchLikely<<endl;
      ofvTotal+= MnCPUELikely;
      //if(DiCatchCompLikelihoodFlag==1)  //Multinomial distribution 
      MvCatchPropLikely=0.0;
      //MnCatchLikely=0.0;
      for(i=1;i<=DiFleetNum;i++)
      {
              for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
              {
                       for(k=1;k<=Flagtimestep;k++)
                      {
                            iTemp=int(Dd3CatchESSInput(i,k,j)+0.00001);
                            //cout<<j<<" "<<k <<" "<<iTemp<<" ";
                            if(iTemp>0 && Dd3CatchTotalO(i,k,j)>MINPOSITIVE)
                            {
                                 iL1=DimCSelStartSizeBin(i,k);
                                 iL2=DimCSelEndSizeBin(i,k);
                                 DvTemp=Dd4CatchPropAtSizeO(i,k,j);
                                 MvTemp1=Md4CatchPropAtSizeP(i,k,j);
                                 //cout<< MvTemp1<<endl;
                                 MvCatchPropLikely(i)+=CalLikelihoodPropVar(DvTemp,MvTemp1,iTemp,iL1,iL2,DimCatchCompLikelihoodFlag(i,k))*DmCatchCompLambda(i,k);
                                 //cout<<MvCatchPropLikely(i)<<endl;
                                 //exit(-99);
                           }
                       }
              }
              MnCatchLikely+=MvCatchPropLikely(i);
      }
      //cout<<"Catch composition likelihood:"<<MnCatchLikely<<endl;
     // cout<<Dd4CatchPropAtSizeO<<endl;
      MnCatchLikely += DnCatchLikely;
      //cout<<"MnCatchLikely 2.5:"<<MnCatchLikely<<endl;
      ofvTotal+=MnCatchLikely;
      //exit(-99);
       /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     if(DiIndexSelParaNumIni>0)
         MvIndexSelLikely=0.0;
     MnIndexSelLikely=0.0;
     for(i=1;i<=DiIndexSelParaNumIni;i++)
     {
             if(active(PbnvIndexSelParams(i))) //DivIndexSelParaPh(i)>0)
             {
                MvIndexSelLikely(i)=CalLikelihoodVar(DvIndexSelPara_Ini(i),DvIndexSelParaSigma2(i),PbnvIndexSelParams(i),DivIndexSelParaLikelihoodFlag(i))*DvIndexSelParaLambda(i);
                MnIndexSelLikely+= MvFleetSelLikely(i);
             }
     }
     MnIndexSelLikely+= DnIndexSelLikely;
     ofvTotal+=MnIndexSelLikely;
      //
     MnIndexLikely=0.0;
     if(DiAvailIndexNum>0)
           MvIndexLikely=0.0;
     iL1=0;
     for(i=1;i<=DiAvailIndexNum;i++)
     {
          for(k=1;k<=DivIndexObsNum(i);k++)
          {
                     iL1++;
                     kk= DimIndexIndex(i,k);//Year
                     MvIndexLikely(i)+= CalLikelihoodVar(DmIndexTotalO(i,k),DmIndexTotalSigma2(i,k),MmIndexP(i,kk),DivIndexTotalLikelihoodFlag(i))*DvIndexTotalLambda(i);      
          }
           MnIndexLikely +=MvIndexLikely(i);
     }
     //cout<<iL1<<endl;
     //cout<<MmIndexP<<endl;
     //cout<<DmIndexTotalO<<endl;
     //cout<<DivIndexObsNum<<endl;
     //cout<<"Survey index likelihood:"<<MnIndexLikely<<endl;
      // 
     if(DiAvailIndexNum>0)
           MvIndexPropLikely=0.0;
     for(i=1;i<=DiAvailIndexNum;i++)
     {
               if(DivEstIndexPropFlag(i)!=1)
                    continue;
               for(k=1;k<=DivIndexObsNum(i);k++)
               {
                      kk= DimIndexIndex(i,k);//Year
                      iTemp=DimIndexESSInput(i,k);
                      if(iTemp>0)
                      {
                                 iL1=DivIndexSelStartSizeBin(i);
                                 iL2=DivIndexSelEndSizeBin(i);
                                 DvTemp=Dd3IndexPropAtSizeO(i,k);
                                 MvTemp1=Md3IndexPropAtSizeP(i,kk);
                                 MvIndexPropLikely(i)+=CalLikelihoodPropVar(DvTemp,MvTemp1,iTemp,iL1,iL2,DivIndexCompLikelihoodFlag(i))*DvIndexCompLambda(i);
                               //  cout<<MvIndexPropLikely(i)<<endl;
                       }
              }
              MnIndexLikely+= MvIndexPropLikely(i);
      }
      //cout<<"Survey composition likelihood:"<<MnIndexLikely<<endl;
      MnIndexLikely +=DnIndexLikely ;
      //cout<<"MnIndexLikely 2.5:"<<MnIndexLikely<<endl;
      ofvTotal+= MnIndexLikely;
      //exit(-99);
      //3/28/2013
      MnLfiftylikely=0.0;
      for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
      {
          iTemp=DimIndexESSInput(DiIndexFlagSexComp,j);
          //cout<<iTemp<<endl;
                      if(iTemp>0)
                      {
                                 DvTemp=DmSexAtSizeData_Ini(j);
                                 MvTemp1=Md3FemalePropAtSizeP(1,j);
                                 MnLfiftylikely+=CalLikelihoodBinomial(DvTemp,MvTemp1,iTemp)*DiSexAtSizeLamda;
                               //  cout<<MvIndexPropLikely(i)<<endl;
                       }
       }
      //cout<<MnLfiftylikely<<endl;
      ofvTotal+= MnLfiftylikely;
        //
      MnIndexQE1Likely=0.0;
      if(DiIndexQBlockNum>0)
             MvIndexQE1Likely=0.0;
      for(i=1;i<=DiIndexQBlockNum;i++)
      {
            if(active(PbnvIndexQBiomassE1(i)))
            {
              MvIndexQE1Likely(i)=CalLikelihoodVar( DvQBiomassE1Value(i),DvQBiomassE1Sigma2(i),PbnvIndexQBiomassE1(i),DivQBiomassE1LikelihoodFlag(i))*DvQBiomassE1Lambda(i);
              MnIndexQE1Likely+= MvIndexQE1Likely(i);
            }
       }
       MnIndexQE1Likely+=DnIndexQE1Likely;
       ofvTotal+=MnIndexQE1Likely;
        //
       MvNYear1Likely=0.0;
       MnNYear1Likely=0.0;
       for(i=1;i<=DiNYear1EstParaNum;i++)
       {
                if(active(PbvNYear1Para(i)))
                {  
                  MvNYear1Likely(i)=CalLikelihoodVar(DvNYear1ParaInital(i),DvNYear1ParaSigma2(i),PbvNYear1Para(i),DivNYear1ParaLikelihoodFlag(i))*DvNYear1ParaLambda(i);
                  MnNYear1Likely+=MvNYear1Likely(i);   
                }
        }
        MnNYear1Likely+=DnNYear1Likely;
        ofvTotal+=MnNYear1Likely;
        //
        MnRSParaLikely=0.0;
        MvRSParaLikely=0.0;
        //
        if(DiRSFlag==1) 
           lL1=1; 
        else 
           lL1=2;
        for(i=1;i<=lL1;i++)
        {
                  if(active(PbnvLogRSPara(i)))
                  {
                      MvRSParaLikely(i)= CalLikelihoodVar( DvRSParaValue(i),DvRSParaSigma2(i),PbnvLogRSPara(i),DivRSParaLikelihoodFlag(i))*DvRSParaLambda(i);
                      MnRSParaLikely+=MvRSParaLikely(i);
                  }
       }
       MnRSParaLikely += DnRSParaLikely;   
       ofvTotal+=MnRSParaLikely ;
       //
       MnRDevsSDLikely=0.0;
       if(DiRecruitLogDevsSDEstFlag)
       {// lL1=DvRecruitLogDevsSD_Ini(4)>0?int(DvRecruitLogDevsSD_Ini(4)+0.01):int(DvRecruitLogDevsSD_Ini(4)-0.01);
                  if(active(PbnvRecruitLogDevsSD2))
                  {
                    if(DvRecruitLogDevsSD_Ini(5)>MINPOSITIVE)
                         nTemp=DvRecruitLogDevsSD_Ini(5)*DvRecruitLogDevsSD_Ini(5);//log(DvRecruitLogDevsSD_Ini(5)*DvRecruitLogDevsSD_Ini(5)+1.0);
                    else
                         nTemp=SIGMA2FILL;
                    MnRDevsSDLikely=CalLikelihoodVar(DvRecruitLogDevsSD_Ini(1),nTemp,PbnvRecruitLogDevsSD2,int(DvRecruitLogDevsSD_Ini(7)+0.01))*DvRecruitLogDevsSD_Ini(6);
                  }
       }
       MnRDevsSDLikely += DnRDevsSDLikely;
       ofvTotal+=MnRDevsSDLikely;//
       //
       MnRecruitDevsLikely=0.0;
       if((active(PbdvLogRecruitDevs)&DiFlagEvntoRdevs)!=1)
       {
                  if(!DiRecruitLogDevsSDEstFlag)
                        nTemp=DnRecruitLogDevsSigma;    
                  //else if(active(PbnvRecruitLogDevsSD2))
                  //      nTemp=value(PbnvRecruitLogDevsSD2);
                  else
                        nTemp=value(MnRSD);
                  nTemp=nTemp*nTemp;
                  for(k=DiCalBeginYear;k<=DiCalEndYear-DiYearBeforeEndForRDev;k++)
                  {//DiCalEndYear-DiYearBeforeEndForRDev
                        if(DvRecruitLogDevsSD_Ini(4)>0)
                             MnRecruitDevsLikely += CalLikelihoodVar(0.0,square(PbnvRecruitLogDevsSD2),PbdvLogRecruitDevs(k),3)*DnRecruitLogDevsLambda;
                        else
                             MnRecruitDevsLikely += CalLikelihoodVar(0.0,nTemp,PbdvLogRecruitDevs(k),3)*DnRecruitLogDevsLambda; 
                         //cout<<MnRSD<<" "<<MnRecruitDevsLikely<<" "<<PbdvLogRecruitDevs(k)<<endl;
                  }
                  for(k=DiCalEndYear-DiYearBeforeEndForRDev+1;k<=DiCalEndYear;k++)
                  {//DiCalEndYear-DiYearBeforeEndForRDev
                        if(DvRecruitLogDevsSD_Ini(4)>0)
                            MnRecruitDevsLikely += CalLikelihoodVar(0.0,square(PbnvRecruitLogDevsSD2),PbdvLogRecruitDevsADD(k),3)*DnRecruitLogDevsLambda;
                        else
                            MnRecruitDevsLikely += CalLikelihoodVar(0.0,nTemp,PbdvLogRecruitDevsADD(k),3)*DnRecruitLogDevsLambda;
                  }
                  //cout<<PbdvLogRecruitDevs<<endl;
                  //cout<<PbdvLogRecruitDevsADD<<endl;
                  //cout<<"nTemp:"<<nTemp<<" "<<"DnRecruitLogDevsLambda"<<DnRecruitLogDevsLambda<<endl;
                  //cout<<"MnRecruitDevsLikely:"<<MnRecruitDevsLikely<<endl;
        }
        MnRecruitDevsLikely+=DnRecruitDevsLikely;
        cout<<"Recruitment likelihood:"<<MnRecruitDevsLikely<<endl;
        ofvTotal+=MnRecruitDevsLikely;//
        //exit(-99);
        MnEnvfitRdevLikely=0;
        if(DiFlagEvntoRdevs==1)
        {
                 for(k=DiCalBeginYear;k<=DiCalEndYear;k++)
                          {
                                 MnEnvfitRdevLikely += CalLikelihoodVar(mfexp(DvEnvfitZscore(k)),StdevEnvZ,mfexp(PbdvLogRecruitDevs(k)),4);
                          }
        //cout<<StdevEnvZ<<endl;
        }       
        cout<<DiFlagEvntoRdevs<<endl;
        cout<<PbdvLogRecruitDevs<<endl;
        //cout<<DvEnvfitZscore<<endl;
        //cout<<MnEnvfitRdevLikely<<endl;
        //exit(789);
        ofvTotal+=MnEnvfitRdevLikely;
        MnRecruitRhLikely=0.0;
        if(active(PbnRecruitmentRh))
        {
          if(DvRecruitmentRh_Ini(5)>MINPOSITIVE)
                         nTemp=DvRecruitmentRh_Ini(5)*DvRecruitmentRh_Ini(5);//log(DvRecruitmentRh_Ini(5)*DvRecruitmentRh_Ini(5)+1.0);
          else
                         nTemp=SIGMA2FILL;
          MnRecruitRhLikely=CalLikelihoodVar(DvRecruitmentRh_Ini(1),nTemp,PbnRecruitmentRh,int(DvRecruitmentRh_Ini(7)+0.01))*DvRecruitmentRh_Ini(6);       
        }
        ofvTotal+= MnRecruitRhLikely;//
       //
        MnLorenzenALikely=0.0;
        if(!DiNaturalMortalityFlag)
        {   
                 MvLorenzenALikely=0.0;
                 for(i=1;i<=Flagtimestep;i++)
                 {
                       if(active(PbnvLorenzenA(i)))
                       {
                         MvLorenzenALikely(i)=CalLikelihoodVar(DvLorenzenAValue(i),DvLorenzenASigma2(i),PbnvLorenzenA(i),DivLorenzenALikelihoodFlag(i))*DvLorenzenALambda(i);
                         MnLorenzenALikely+=MvLorenzenALikely(i);
                        }
                  }
                  MnLorenzenBLikely=0.0;
                  for(i=1;i<=Flagtimestep;i++)
                  {
                         if(active(PbnvLorenzenB(i)))
                         {
                            MvLorenzenBLikely(i)=CalLikelihoodVar(DvLorenzenBValue(i),DvLorenzenBSigma2(i),PbnvLorenzenB(i),DivLorenzenBLikelihoodFlag(i))*DvLorenzenBLambda(i);
                            MnLorenzenBLikely+=MvLorenzenBLikely(i);
                          }
                  }
        }
        MnLorenzenALikely+=DnLorenzenALikely;
        ofvTotal+=MnLorenzenALikely;
       MnLorenzenBLikely+=DnLorenzenBLikely;
       ofvTotal+= MnLorenzenBLikely;        
        MnGrowthMatrixLikely=0.0;
        if(DiGrowthMatrixFlag)
        {
                MvLinfLikely=0.0;
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(active(PbnvLinf(i)))
                     {
                       MvLinfLikely(i)= CalLikelihoodVar(DvVBGFLinfValue(i),DvVBGFLinfSigma2(i),PbnvLinf(i),DivVBGFLinfLikelihoodFlag(i))*DvVBGFLinfLambda(i);
                       MnGrowthMatrixLikely+=MvLinfLikely(i);
                     }
                }
                MvVLikely=0.0;
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(active(PbnvV(i)))
                     {
                       MvVLikely(i)= CalLikelihoodVar(DvVBGFVValue(i),DvVBGFVSigma2(i),PbnvV(i),DivVBGFVLikelihoodFlag(i))*DvVBGFVLambda(i);
                       MnGrowthMatrixLikely+=MvVLikely(i);
                     }
                }
                MvLSDLikely=0.0;
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(active(PbnvLSD(i)))
                     {
                       MvLSDLikely(i)= CalLikelihoodVar(DvVBGFLSDValue(i),DvVBGFLSDSigma2(i),PbnvLSD(i),DivVBGFLSDLikelihoodFlag(i))*DvVBGFLSDLambda(i);
                       MnGrowthMatrixLikely+=MvLSDLikely(i);
                     }
                }
                MvVSDLikely=0.0;
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(active(PbnvVSD(i)))
                     {
                       MvVSDLikely(i)= CalLikelihoodVar(DvVBGFVSDValue(i),DvVBGFVSDSigma2(i),PbnvVSD(i),DivVBGFVSDLikelihoodFlag(i))*DvVBGFVSDLambda(i);
                       MnGrowthMatrixLikely+=MvVSDLikely(i);
                     }
                }
                MvLVRhoLikely=0.0;
                for(i=1;i<=DiGrowthMatrixBlockNum;i++)
                {
                     if(active(PbnvLVRho(i)))
                     {
                       MvLVRhoLikely(i)= CalLikelihoodVar(DvVBGFLVRhoValue(i),DvVBGFLVRhoSigma2(i),PbnvLVRho(i),DivVBGFLVRhoLikelihoodFlag(i))*DvVBGFLVRhoLambda(i);
                       MnGrowthMatrixLikely+=MvLVRhoLikely(i);
                     }
                }
    }       
    MnGrowthMatrixLikely +=DnGrowthMatrixLikely;  
    ofvTotal+= MnGrowthMatrixLikely;
    //FMax
    MnFMaxPenaltyLikely=0.0;
    for (i=1;i<=DiFleetNum;i++)
    {
       for(k=1;k<=Flagtimestep;k++)
       {
          for (j=DiCalBeginYear;j<=DiCalEndYear;j++)
          {
             MnTemp1=Md3Ft(i,k,j);
             if(MnTemp1>DnFmax)
               MnFMaxPenaltyLikely+=1000.*(MnTemp1-DnFmax)*(MnTemp1-DnFmax);
          }
      }
    }
    if(!last_phase())
    {
         if(DnRecruitLogDevsLambda<0.000000001)
          {
                     DnRecruitLogDevsLambda=0.05;
                     DiPenaltyLambdaFlag=1;
          }
         /*for(i=DiCalBeginYear;i<=DiCalEndYear;i++)
         {
              if(sdr_vRecruitmentP(i)<0.01*MvRecruitRS(i)||sdr_vRecruitmentP(i)>100*MvRecruitRS(i))
              {
                  //cout<<"Penalty:"<<i<<" "<<sdr_vRecruitmentP(i)<<" "<<MvRecruitRS(i)<<endl;
                  MnTemp1=sdr_vRecruitmentP(i)-MvRecruitRS(i);//sum(PbdvLogRecruitDevs(DiCalBeginYear+1,DiCalEndYear-DiYearBeforeEndForRDev-1));
                  MnFMaxPenaltyLikely+=MnTemp1*MnTemp1*0.1;
              }
         }//*/
    }
    else
    {
            if(DiPenaltyLambdaFlag)
            {
                 DnRecruitLogDevsLambda=0.0;
                 DiPenaltyLambdaFlag=0;
            }
           /* for(i=DiCalBeginYear+1;i<=DiCalEndYear;i++)
            {
                 if(sdr_vRecruitmentP(i)<0.05*sdr_vRecruitmentP(1))
                 {
                     MnTemp1=sdr_vRecruitmentP(i)-sdr_vRecruitmentP(1);
                     MnFMaxPenaltyLikely+=MnTemp1*MnTemp1*0.1;
                 }
           }*/
    }
    //cout<<"Penalty:"<<MnFMaxPenaltyLikely<<endl;
    //cout<<"PbdvLogRecruitDevs:"<<PbdvLogRecruitDevs<<endl
    ofvTotal+= MnFMaxPenaltyLikely;
    MnRPVLikely=0.0;// PbnvRecruitPrjVect(i)=DbnvRecruitPrjVect(i);
    for(i=1;i<=DiSizeBinNum;i++)
    {
        if(active(PbnvRecruitPrjVect(i)))
        {
           MnLorenzenALikely+=CalLikelihoodVar(DbnvRecruitPrjVect(i),DnRPVSigma2,PbnvRecruitPrjVect(i),DiRPVLikelihoodFlag)*DnRPVLambda;
        }
    }
    //cout<<"MnRPVLikely:"<<MnRPVLikely<<endl;
    /////////////////////////////////////////////////////////
    ofvTotal+=MnRPVLikely;
    //cout<<"Objective Function Value:"<< ofvTotal<<endl;
}

dvariable model_parameters::CalLikelihoodPropVar2(const dvector dvrObs,const dvar_vector& dvvPred, const dvariable &dvESS,const int &iB,const int &iE)
{
           dvariable dvR;
           dvariable dvChi;
           dvariable dvTau; 
           dvariable dvVariance;
           dvariable dvResids;
           int iNBN;
           int k;
           iNBN=iE-iB+1;
           dvR=0.0;
           dvTau=1.0/dvESS; 
           for(k=iB;k<=iE;k++)
           {
                         dvChi=dvvPred(k)*(1.0-dvvPred(k));
		         dvVariance=(dvChi+0.1/iNBN)*dvTau*dvTau;
                         dvR -=log(1.0/sqrt(PI2*dvVariance));
		         dvResids=dvrObs(k)-dvvPred(k);
		         dvR -=log(mfexp(-0.5*square(dvResids)/dvVariance)+0.01);
            }                    
             return dvR;
}

dvariable model_parameters::CalLikelihoodPropVar2(const dvector dvrObs,const dvar_vector& dvvPred, const int &iESS,const int &iB,const int &iE)
{
           dvariable dvR;
           dvariable dvChi;
           dvariable dvVariance;
           dvariable dvResids;
           int iNBN;
           double dTau; 
           int k;
           iNBN=iE-iB+1;
           dvR.initialize();
           dTau=1.0/iESS; 
           for(k=iB;k<=iE;k++)
           {
                         dvChi=dvvPred(k)*(1.0-dvvPred(k));
		         dvVariance=(dvChi+0.1/iNBN)*dTau*dTau;
                         dvR +=log(dvVariance)*0.5;//dvR +=log(sqrt(dvVariance));//dvR -=log(1.0/sqrt(PI2*dvVariance));
		         dvResids=dvrObs(k)-dvvPred(k);
		         dvR -=log(mfexp(-0.5*square(dvResids)/dvVariance)+0.01);
            }                    
             return dvR;
}

dvariable model_parameters::CalLikelihoodBinomial(const dvector dvrObs,const dvar_vector& dvvPred, const int &iESS)
{
           dvariable dvR;   
           int k;       
           dvR.initialize();
           for(k=1;k<=DiSizeBinNum;k++)
           {       
                      dvR -=iESS*dvrObs(k)*log(dvvPred(k)+ MINPOSITIVE)+iESS*(1-dvrObs(k))*log(1-dvvPred(k)+ MINPOSITIVE);//??????????????????  log(0) 
                      //dvR -=iESS*dvrObs(k)*log(dvvPred(k)/dvrObs(k));
           }          
           return dvR;
}

dvariable model_parameters::CalLikelihoodPropVar1(const dvector dvrObs,const dvar_vector& dvvPred, const int &iESS, const int &iB, const int &iE)
{
           dvariable dvR;   
           int k;       
           dvR.initialize();
           for(k=iB;k<=iE;k++)
           {       
                      dvR -=iESS*dvrObs(k)*log(dvvPred(k)+ MINPOSITIVE);//??????????????????  log(0) 
                      //dvR -=iESS*dvrObs(k)*log(dvvPred(k)/dvrObs(k));
           }          
           return dvR;
}

dvariable model_parameters::CalLikelihoodPropVar(const dvector dvrObs,const dvar_vector& dvvPred, const int &iESS, const int &iB, const int &iE,int iModel)
{
         dvariable dvR;
         dvR=0.0;
         switch(iModel)
         {
             case 1:
                 dvR=CalLikelihoodPropVar1(dvrObs,dvvPred,iESS,iB,iE);
                 break;
             case 2:
                 dvR=CalLikelihoodPropVar2(dvrObs,dvvPred,iESS,iB,iE);
                 break;
             default:
                break;
         }
         return dvR;  
}

dvariable model_parameters::CalLikelihoodPropVar0(const dvector dvrObs,const dvar_vector& dvvPred, const int &iESS, const int &iB, const int &iE,int iModel)
{
           dvariable dvR;          
           double dTemp;
           int k;
           dvR.initialize();
           dTemp=0.0;
           switch(iModel)
           {
              case 1:
                  dTemp= CalLikelihoodPropVar0(dvrObs,iESS,iB,iE,0);
              case 0:
                 dvR=dTemp;
                 for(k=iB;k<=iE;k++)
                 {       
                      dvR -=iESS*dvrObs(k)*log(dvvPred(k)+ MINPOSITIVE);//??????????????????  log(0)
                 }     
           }     
           return dvR;
}

double model_parameters::CalLikelihoodPropVar0(const dvector &dvrObs, const int &iESS, const int &iB, const int &iE, const int &iModel)
{
           double dR;          
           //double dTemp;
           int k;
           dR=0.0;
           dR-=CalLogFactorial(iESS);
           for(k=iB;k<=iE;k++)
           {        
                 iTemp=int(iESS*dvrObs(k)+0.5);
                 dR+=CalLogFactorial(iTemp);
           }          
           return dR;
}

void model_parameters::CalNaturalMortality(void)
{
  if(!DiNaturalMortalityFlag)
  {
      for(i=1;i<=DiSeasonNum;i++)
      {
         if(active(PbnvLorenzenA(i))||active(PbnvLorenzenB(i)))
         {
             Md3M(i).initialize();
             for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
            {
               for(k=1;k<=DiSizeBinNum;k++)
               {
                     Md3M(i,j,k)=PbnvLorenzenA(i)*pow(DmWeightAtSize(j,k),PbnvLorenzenB(i))*DvNaturalMSizeWeight(k)*DvNaturalMYearWeight(j);
               }
            }
        }
      }
   }
   /*else
   {//
       for(i=1;i<=DiSeasonNum;i++)
       {
         for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
         {
               for(k=1;k<=DiSizeBinNum;k++)
               {
                     iL1=(i-1)*DiYearNum+j;
                     Md3M(i,j,k)=DmNaturalMAtSize_Ini(iL1,k+2);
               }
          }
       }      
   }*/
 #if Debug_Status
 // ICHECKDEB(Md3M);
 #endif
  return;
}

void model_parameters::CalNaturalMortality_Year(void)
{
  if(!DiNaturalMortalityFlag)
  {
         if(active(PbnvLorenzenA(1))||active(PbnvLorenzenB(1)))
         {
             Md3M(1).initialize();
             for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
            {
               for(k=1;k<=DiSizeBinNum;k++)
               {     
                     Md3M(1,j,k)=PbnvLorenzenA(1)*pow(DmWeightAtSize(j,k),PbnvLorenzenB(1))*DvNaturalMSizeWeight(k)*DvNaturalMYearWeight(j);
               }
            }
        }
   }
  return;
}

void model_parameters::GetBPRSel(void)
{
   //MvTemp1.initialize();
   int i,k;
   dvariable   dvSelMax;
   MmBPRSel.initialize();
   if(DiBPRSelectivityFlag==-1)
   {
     for(k=1;k<=Flagtimestep;k++)
        MmBPRSel(k)=DmBPRSelectivityIni(k);
   }
   else
   {
      if(DiBPRSelectivityFlag==0)
      {
         for(k=1;k<=Flagtimestep;k++)
         {
           // MvTemp1.initialize();
           for(i=1;i<=DiFleetNum;i++)
           {   
              MmBPRSel(k)+=Md4FASByFleet(i,k,DiBPRCalYear);
           }
           dvSelMax=max(MmBPRSel(k));
           if(dvSelMax>MINPOSITIVE)
                 MmBPRSel(k)= MmBPRSel(k)/dvSelMax;
           //MmBPRSel(k)= MmBPRSel(k)/max(MmBPRSel(k));
         }
      }
      else if(DiBPRSelectivityFlag>=1 && DiBPRSelectivityFlag<=DiFleetNum)
      {
            for(k=1;k<=Flagtimestep;k++)
            {
              MmBPRSel(k)=Md4SelByFleet(DiBPRSelectivityFlag,k,DiBPRCalYear);
            }
      }
      else if(DiBPRSelectivityFlag>DiFleetNum)
      {
          dmFRInFleetSeason.initialize();
          dvFRInSeason.initialize();
          for(k=1;k<=Flagtimestep ;k++)
          {
             for(i=1;i<=DiFleetNum;i++)
             {
                  dmFRInFleetSeason(k,i)=value(Md3Ft(i,k,DiBPRCalYear));
             }
             dvFRInSeason(k)=sum(dmFRInFleetSeason(k));
             dmFRInFleetSeason(k)/=dvFRInSeason(k);
          }
          DvFRInSeason=dvFRInSeason/sum(dvFRInSeason);
         // ofRuntimelog<<"DvFRInSeason"<<DvFRInSeason<<endl;
         // ofRuntimelog<<"dmFRInFleetSeason"<<dmFRInFleetSeason<<endl;
     }
  }
}

void model_parameters::GetSPRAndYPR(const dvariable &dvFtot,int iWeightFlag,int iModelFlag)
{
   double      dR;
   dvariable   dvF;
   //dmatrix     dmFRInFleetSeason(1,DiSeasonNum,1,DiFleetNum);
   //dvector     dvFRInSeason(1,DiSeasonNum);
   int k,i,kk,iL1;
   //DiBPRBlockGrowthMatrix DvFRInSeason
   MnBPRYieldTotal.initialize();
   MnBPRSSBTotal.initialize();
   MvBPRNAS.initialize();
   dR=1.0;//unit recruitment because we assume the balance was reached!
   //
   for(k=1;k<=DiSeasonNum;k++)
   {
       dvF=dvFtot*DvFRInSeason(k);
       if(DiBPRSelectivityFlag<=DiFleetNum)
       {
             MmBPRF(k)=dvF*MmBPRSel(k);
       }
       else
       {     
             MmBPRF(k).initialize();
             for(i=1;i<=DiFleetNum;i++)
             {
                 MmBPRF(k)+=dvF* dmFRInFleetSeason(k,i)*Md4SelByFleet(i,k,DiBPRCalYear);
             }  
       }
       MmBPRZ(k)=MmBPRF(k)+Md3M(k,DiBPRCalYear);
       MmBPRS(k)=mfexp(-MmBPRZ(k));
       if(k==DiSBBSeason)
           MvBPRSSB_S=mfexp(-MmBPRZ(k)*DnFracSeasonSSB);
   }
   //
   for(i=1;i<=DiBPRPeriod;i++)
   {
      for(k=1;k<=DiSeasonNum;k++)
      {
          if( iModelFlag)//all cohorts provide yield and ssb in one year 
          {
                   MvBPRNAS+= dR*DvRecruitmentSeaRate(k)*MvRecruitPrjVect;
          }
          else
          {//each cohort can provide yield and ssb
               if(i==1)
                    MvBPRNAS+= dR*DvRecruitmentSeaRate(k)*MvRecruitPrjVect;
               CalCatchInBPR(DiBPRCalYear,k,0.0,iWeightFlag);
               MnBPRYieldTotal+= MnBPRYield;
               if(k==DiSBBSeason)
                   MnBPRSSBTotal  +=MnBPRSSB;
          }
          //ofRuntimelog<<"Yield "<<MnBPRYieldTotal<<"  "<<MnBPRYield<<endl;
          //ofRuntimelog<<"NAS: "<<MvBPRNAS<<endl;
          MvNASTemp=elem_prod(MvBPRNAS,MmBPRS(k));                                    //die        N(t)=(N(t-1)*exp(-Z)*G)
          if(DiNaturalMortalityFlag==2)
          {
                iL1= DivBPRBlockGrowthMatrix(k);
                iL1=(iL1-1)*DiSeasonNum+k;
          }
          else
          {
              iL1=DivBPRBlockGrowthMatrix(k);
          }
          MvBPRNAS= MvNASTemp*Md3GrowthMatrix(iL1);           //Growth 
      }
  }
  if( iModelFlag)
  {
       for(k=1;k<=DiSeasonNum;k++)
       {
           MvBPRNAS+= dR*DvRecruitmentSeaRate(k)*MvRecruitPrjVect;
           if(iWeightFlag)
           {
                if(k==DiSBBSeason)
                         MnBPRSSBTotal=elem_prod(MvBPRNAS,MvBPRSSB_S)*DmFecundity(DiBPRCalYear); 
                MnBPRYieldTotal+=elem_prod(elem_div(MmBPRF(k),MmBPRZ(k)),elem_prod(MvBPRNAS,1.0-MmBPRS(k)))*DmWeightAtSize(DiBPRCalYear);
                //ofRuntimelog<<MnBPRYield*2062<<endl;
           }
           else
           {
                if(k==DiSBBSeason)
                         MnBPRSSBTotal=sum(elem_prod(MvBPRNAS,MvBPRSSB_S)); 
                MnBPRYieldTotal+=sum(elem_prod(elem_div( MmBPRF(k),MmBPRZ(k)),elem_prod(MvBPRNAS,1.0-MmBPRS(k)))); 
           }//dvFtot
            if(dvFtot<0.133 && dvFtot>0.132)
              ofRuntimelog<<"F:\t"<<dvFtot<<" FnBPRYieldTotal\t:"<< MnBPRYieldTotal<<"FnBPRSSBTotal:\t"<<MnBPRSSBTotal<<endl;
          MvNASTemp=elem_prod(MvBPRNAS,MmBPRS(k));                           //die        N(t)=(N(t-1)*exp(-Z)*G)
          if(DiNaturalMortalityFlag==2)
          {
                iL1= DivBPRBlockGrowthMatrix(k);
                iL1=(iL1-1)*DiSeasonNum+k;
          }
          else
          {
                iL1=DivBPRBlockGrowthMatrix(k);
          }
          MvBPRNAS= MvNASTemp*Md3GrowthMatrix(iL1);  
       }
  }
  return ;
}

void model_parameters::GetSPRAndYPR_Y(const dvariable &dvFtot,int iWeightFlag,int iModelFlag)
{
   double      dR;
   dvariable   dvF;
   //dmatrix     dmFRInFleetSeason(1,DiSeasonNum,1,DiFleetNum);
   //dvector     dvFRInSeason(1,DiSeasonNum);
   int k,i,kk,iL1;
   //DiBPRBlockGrowthMatrix DvFRInSeason
   MnBPRYieldTotal.initialize();
   MnBPRSSBTotal.initialize();
   MvBPRNAS.initialize();
   dR=1.0;//unit recruitment because we assume the balance was reached!
   //
   for(k=1;k<=Flagtimestep;k++)
   {
       dvF=dvFtot*DvFRInSeason(k);
       if(DiBPRSelectivityFlag<=DiFleetNum)
       {
             MmBPRF(k)=dvF*MmBPRSel(k);
       }
       else
       {     
             MmBPRF(k).initialize();
             for(i=1;i<=DiFleetNum;i++)
             {
                 MmBPRF(k)+=dvF* dmFRInFleetSeason(k,i)*Md4SelByFleet(i,k,DiBPRCalYear);
             }  
       }
       MmBPRZ(k)=MmBPRF(k)+Md3M(k,DiBPRCalYear);
       MmBPRS(k)=mfexp(-MmBPRZ(k));
       MvBPRSSB_S=mfexp(-MmBPRZ(k)*DnFracYearSSB);
   }
   //cout<<MvBPRSSB_S<<endl;
   //
   for(i=1;i<=DiBPRPeriod;i++)
   {
      for(k=1;k<=Flagtimestep;k++)
      {
          if( iModelFlag)//all cohorts provide yield and ssb in one year 
          {
                   MvBPRNAS+= dR*DvRecruitmentSeaRate(k)*MvRecruitPrjVect;
          }
          else
          {//each cohort can provide yield and ssb
               if(i==1)
                    MvBPRNAS+= dR*DvRecruitmentSeaRate(k)*MvRecruitPrjVect;
               CalCatchInBPR(DiBPRCalYear,k,0.0,iWeightFlag);
               MnBPRYieldTotal+= MnBPRYield;
               if(k==DiSBBSeason)
                   MnBPRSSBTotal  +=MnBPRSSB;
          }
          //cout<<MvBPRNAS<<endl;
          //ofRuntimelog<<"Yield "<<MnBPRYieldTotal<<"  "<<MnBPRYield<<endl;
          //ofRuntimelog<<"NAS: "<<MvBPRNAS<<endl;
          MvNASTemp=elem_prod(MvBPRNAS,MmBPRS(k));                                    //die        N(t)=(N(t-1)*exp(-Z)*G)
          //cout<<MvNASTemp<<endl;
          if(DiNaturalMortalityFlag==2)
          {
                iL1= DivBPRBlockGrowthMatrix(k);
                iL1=(iL1-1)*DiSeasonNum+k;
          }
          else
          {
              iL1=DivBPRBlockGrowthMatrix(k);
          }
          MvBPRNAS= MvNASTemp*Md3GrowthMatrix(iL1);           //Growth 
          //cout<<MvBPRNAS<<endl;
      }
  }
  //cout<<Md3GrowthMatrix(iL1)<<endl;
  if( iModelFlag)
  {
       for(k=1;k<=Flagtimestep;k++)
       {
           MvBPRNAS+= dR*DvRecruitmentSeaRate(k)*MvRecruitPrjVect;
           if(iWeightFlag)
           {
                MnBPRSSBTotal=elem_prod(MvBPRNAS,MvBPRSSB_S)*DmFecundity(DiBPRCalYear); 
                MnBPRYieldTotal+=elem_prod(elem_div(MmBPRF(k),MmBPRZ(k)),elem_prod(MvBPRNAS,1.0-MmBPRS(k)))*DmWeightAtSize(DiBPRCalYear);
                //cout<<DmFecundity(4,DiBPRCalYear)<<endl;
                //cout<<elem_prod(MvBPRNAS,MvBPRSSB_S)<<endl;
                //ofRuntimelog<<MnBPRYield*2062<<endl;
           }
           else
           {
                MnBPRSSBTotal=sum(elem_prod(MvBPRNAS,MvBPRSSB_S)); 
                MnBPRYieldTotal+=sum(elem_prod(elem_div( MmBPRF(k),MmBPRZ(k)),elem_prod(MvBPRNAS,1.0-MmBPRS(k)))); 
           }//dvFtot
            if(dvFtot<0.133 && dvFtot>0.132)
              ofRuntimelog<<"F:\t"<<dvFtot<<" FnBPRYieldTotal\t:"<< MnBPRYieldTotal<<"FnBPRSSBTotal:\t"<<MnBPRSSBTotal<<endl;
          MvNASTemp=elem_prod(MvBPRNAS,MmBPRS(k));                           //die        N(t)=(N(t-1)*exp(-Z)*G)
          if(DiNaturalMortalityFlag==2)
          {
                iL1= DivBPRBlockGrowthMatrix(k);
                iL1=(iL1-1)*DiSeasonNum+k;
          }
          else
          {
                iL1=DivBPRBlockGrowthMatrix(k);
          }
          MvBPRNAS= MvNASTemp*Md3GrowthMatrix(iL1);  
       }
  }
  return ;
}

void model_parameters::CalCatchInBPR(int iYear,int iSeason,const double &dR,const int iFlag )
{
    int k;
    MnBPRSSB.initialize();
    MnBPRYield.initialize();
    //for (k=1;k<=DiSizeBinNum;k++)  
    //      MvBPRNAS(k)+= dvR*DvRecruitmentSeaRate(iSeason)*MvRecruitPrjVect(k);
    //MvBPRNAS+= dR*DvRecruitmentSeaRate(k)*MvRecruitPrjVect;
    if(iFlag)
    {
                if(iSeason==DiSBBSeason)
                        MnBPRSSB=elem_prod(MvBPRNAS,MvBPRSSB_S)*DmFecundity(iYear); 
                MnBPRYield=elem_prod(elem_div(MmBPRF(iSeason),MmBPRZ(iSeason)),elem_prod(MvBPRNAS,1.0-MmBPRS(iSeason)))*DmWeightAtSize(iYear);
                //ofRuntimelog<<MnBPRYield*2062<<endl;
    }
    else
    {
                if(iSeason==DiSBBSeason)
                        MnBPRSSB=sum(elem_prod(MvBPRNAS,MvBPRSSB_S)); 
                MnBPRYield=sum(elem_prod(elem_div( MmBPRF(iSeason),MmBPRZ(iSeason)),elem_prod(MvBPRNAS,1.0-MmBPRS(iSeason)))); 
    }
}

dvariable model_parameters::MSY_Sub(dvariable dF,dvariable dDelta,int iWeightFlag,int iModelFlag,int iRS)
{
   dvariable dY1, dY2,dR;
   dvariable dalpha,dbeta;
   if(iRS==1)
   {
        dalpha=mfexp(PbnvLogRSPara(1));
        dbeta=0.0;
   }
   else
   {
      dalpha=mfexp(PbnvLogRSPara(1));
      dbeta=mfexp(PbnvLogRSPara(2));
   }
   //MnTemp1=dF;
   if (Flagtimestep==4)
       GetSPRAndYPR(dF, iWeightFlag, iModelFlag);
   if (Flagtimestep==1)
       GetSPRAndYPR_Y(dF, iWeightFlag, iModelFlag);
  // pt=log(value( nBPRYieldTotal)+1.0);
   dY1=CalBPRCatch(dalpha,dbeta,MnBPRSSBTotal,MnBPRYieldTotal,DiRSFlag);
   ofRuntimelog<<"DY1:"<<dY1<<" "<<MnTemp1<<" "<<dalpha<<" "<<dbeta <<endl;
  // if(iRS!=3)
  // {
 //     dR=value(nAlpha*nBPRSSBTotal-nBelta)/value(nBPRSSBTotal);
 //     dY1=dR*value( nBPRYieldTotal);
 //     
 //  }
 //  else
 //  {
 //        dR=value(log(nAlpha)+log(nBPRSSBTotal))/value(nBelta);
 //        dR/=value(nBPRSSBTotal);
 //        dY1=dR*value( nBPRYieldTotal);
 //  }
   MnTemp1=dF+ dDelta;
   if (Flagtimestep==4)
       GetSPRAndYPR(MnTemp1, iWeightFlag, iModelFlag);
   if (Flagtimestep==1)
       GetSPRAndYPR_Y(MnTemp1, iWeightFlag, iModelFlag);
   dY2=CalBPRCatch(dalpha,dbeta,MnBPRSSBTotal,MnBPRYieldTotal,DiRSFlag);
    ofRuntimelog<<"DY2:"<<dY2<<" "<<MnTemp1<<" "<<dalpha<<" "<<dbeta <<endl;
   //ofRuntimelog<<"DY2:"<<dY2<<" "<<MnTemp1<<endl;
   // if(iRS!=3)
  // {
  //    dR=value(nAlpha*nBPRSSBTotal-nBelta)/value(nBPRSSBTotal);
  //    dY2=dR*value( nBPRYieldTotal);
  // }
  // else
  // {
  //       dR=value(log(nAlpha)+log(nBPRSSBTotal))/value(nBelta);
  //       dR/=value(nBPRSSBTotal);
  //       dY2=dR*value( nBPRYieldTotal);
  // }
  // dpt=log(value( nBPRYieldTotal)+1.0);
  // nTemp1=nAlpha*nBPRSSBTotal-nBelta;
   dR=(dY2-dY1);///dDelta;//dY2-dY1;//(dY2-dY1)/dDelta;
   return dR;
}

dvariable model_parameters::MSY_Bisection(dvariable dF1, dvariable dF2, dvariable dDelta,int iWeightFlag,int iModelFlag,int iRS)
{
   dvariable bf,fmid;
   int j;
   dvariable  xmid,rtb,dx;
   bf=MSY_Sub(dF1,dDelta/3.0, iWeightFlag, iModelFlag,iRS);
   fmid=MSY_Sub(dF2,dDelta/3.0, iWeightFlag, iModelFlag,iRS);
   ofRuntimelog<<"bf:\t"<<bf<<endl;
   ofRuntimelog<<"fmid:\t"<<fmid<<endl;
   if(fabs(bf)<DIGITALPREC)
      return dF1;
   if(fabs(fmid)<DIGITALPREC)
      return dF2;
   if (bf*fmid > 0.0)
   {
            ofRuntimelog<<"Root must be bracketed for bisection in rtbis of MSY"<<endl;
   }
   else
   {
             rtb = bf < 0.0 ? (dx=dF1-dF2,dF2): (dx=dF2-dF1,dF1); //Orient the search so that f>0
             for (j=1;j<=20;j++) 
             { //lies at x+dx.
                     dx = dx *0.5;
                     xmid=rtb+dx;
                     fmid=MSY_Sub(xmid,dDelta/3.0, iWeightFlag, iModelFlag,iRS); //Bisection loop.
                     //cout<<"MSYtestF = "<<xmid<<" Dyield = "<<fmid<<endl;
                     if(fabs(fmid)<DIGITALPREC)
                     {
                            rtb=xmid;
                            return rtb;
                     }
                     if (fabs(dx) < dDelta) 
                           return rtb;
                     if (fmid > 0.0) 
                            rtb=xmid;
             }
               ofRuntimelog<<"Too many bisections in rtbis"<<endl;
   }
   return 0.0;
}

void model_parameters::GetSPRAndYPR(double &FnBPRSSBTotal,double &FnBPRYieldTotal,double dvFtot,int iWeightFlag,int iModelFlag)
{
   double   dR;
   double   dvF;
   dmatrix  FmBPRF(1,DiSeasonNum,1,DiSizeBinNum);
   dmatrix  FmBPRZ(1,DiSeasonNum,1,DiSizeBinNum);
   dmatrix  FmBPRS(1,DiSeasonNum,1,DiSizeBinNum);
   dvector  FvBPRSSB_S(1,DiSizeBinNum);
   //d3array 
   //d3array  Fd3GrowthMatrix(1,DiSeasonNum,1,DiSizeBinNum,1,DiSizeBinNum);
   dvector  FvBPRNAS(1,DiSizeBinNum);
   dvector  FvNASTemp(1,DiSizeBinNum);
   //dmatrix  dmFRInFleetSeason(1,DiSeasonNum,1,DiFleetNum);
   //dvector  dvFRInSeason(1,DiSeasonNum);
   int k,i,kk,iL;
   FnBPRYieldTotal=0.0;
   FnBPRSSBTotal=0.0;
   dR=1.0;//unit recruitment because we assume the balance was reached!
   //dvector  dvFRInSeason(1,DiSeasonNum)
   for(k=1;k<=DiSeasonNum;k++)
   {
       dvF=dvFtot*DvFRInSeason(k);
       if(DiNaturalMortalityFlag==2)
       {
                iL= DivBPRBlockGrowthMatrix(k);
                iL=(iL-1)*DiSeasonNum+k;
       }
       else
       {
                iL=DivBPRBlockGrowthMatrix(k);
       }
       FmBPRF(k).initialize();
       for(kk=1;kk<=DiSizeBinNum;kk++)
       {
          if(DiBPRSelectivityFlag<=DiFleetNum)
                 FmBPRF(k,kk)=dvF*value(MmBPRSel(k,kk));
          else
          {
             for(i=1;i<=DiFleetNum;i++)
             {
                 FmBPRF(k,kk)+=dvF* dmFRInFleetSeason(k,i)*value(Md4SelByFleet(i,k,DiBPRCalYear,kk));
             }  
          }
           FmBPRZ(k,kk)=FmBPRF(k,kk)+value(Md3M(k,DiBPRCalYear,kk));
           for(i=1;i<=DiSizeBinNum;i++)
           {
              Fd3GrowthMatrix(k,kk,i)= value(Md3GrowthMatrix(iL,kk,i));
           }
       }
       FmBPRS(k)=mfexp(-FmBPRZ(k));
       if(k==DiSBBSeason)
               FvBPRSSB_S=mfexp(-FmBPRZ(k)*DnFracSeasonSSB);
   }
   //
   for(i=1;i<=DiBPRPeriod;i++)
   {
      for(k=1;k<=DiSeasonNum;k++)
      {
          if( iModelFlag)//all cohorts provide yield and ssb in one year 
          {
                  //cout<<"Slow Here!"<<endl;
                 // for (kk=1;kk<=DiSizeBinNum;kk++)  
                 //        MvBPRNAS(kk)+= dR*DvRecruitmentSeaRate(k)*MvRecruitPrjVect(kk);
                   FvBPRNAS+= dR*DvRecruitmentSeaRate(k)*value(MvRecruitPrjVect);
                // CalCatchInBPR(DiBPRCalYear,k,dvR,iWeightFlag);
          }
          else
          {//each cohort can provide yield and ssb
               if(i==1)
                  FvBPRNAS+= dR*DvRecruitmentSeaRate(k)*value(MvRecruitPrjVect);
               if(iWeightFlag)
               {
                    if(k==DiSBBSeason)
                        FnBPRSSBTotal+=elem_prod(FvBPRNAS,FvBPRSSB_S)*DmFecundity(DiBPRCalYear); 
                    FnBPRYieldTotal+=elem_prod(elem_div(FmBPRF(k),FmBPRZ(k)),elem_prod(FvBPRNAS,1.0-FmBPRS(k)))*DmWeightAtSize(DiBPRCalYear);
               }
               else
               {
                    if(k==DiSBBSeason)
                        FnBPRSSBTotal+=sum(elem_prod(FvBPRNAS,FvBPRSSB_S)); 
                    FnBPRYieldTotal+=elem_prod(elem_div(FmBPRF(k),FmBPRZ(k)),elem_prod(FvBPRNAS,1.0-FmBPRS(k)))*DmWeightAtSize(DiBPRCalYear);
               }
          }
          //ofRuntimelog<<"Yield "<<MnBPRYieldTotal<<"  "<<MnBPRYield<<endl;
          //ofRuntimelog<<"NAS: "<<MvBPRNAS<<endl;
          FvNASTemp=elem_prod(FvBPRNAS,FmBPRS(k));                                    //die        N(t)=(N(t-1)*exp(-Z)*G)
          FvBPRNAS= FvNASTemp*Fd3GrowthMatrix(k);           //Growth 
      }
  }
  if( iModelFlag)
  {
      // if(dvFtot>1.95||dvFtot<0.001)
      //     ofRuntimelog<<"F:\t"<<dvFtot<<"FvBPRNAS:"<<endl<<FvBPRNAS<<endl;
       for(k=1;k<=DiSeasonNum;k++)
       {
           FvBPRNAS+= dR*DvRecruitmentSeaRate(k)*value(MvRecruitPrjVect);
           if(iWeightFlag)
           {
                    if(k==DiSBBSeason)
                        FnBPRSSBTotal+=elem_prod(FvBPRNAS,FvBPRSSB_S)*DmFecundity(DiBPRCalYear); 
                    FnBPRYieldTotal+=elem_prod(elem_div(FmBPRF(k),FmBPRZ(k)),elem_prod(FvBPRNAS,1.0-FmBPRS(k)))*DmWeightAtSize(DiBPRCalYear);
           }
           else
           {
                    if(k==DiSBBSeason)
                        FnBPRSSBTotal+=sum(elem_prod(FvBPRNAS,FvBPRSSB_S)); 
                    FnBPRYieldTotal+=elem_prod(elem_div(FmBPRF(k),FmBPRZ(k)),elem_prod(FvBPRNAS,1.0-FmBPRS(k)))*DmWeightAtSize(DiBPRCalYear);
           }
          if(dvFtot<0.133 && dvFtot>0.132)
              ofRuntimelog<<"F:\t"<<dvFtot<<" FnBPRYieldTotal\t:"<< FnBPRYieldTotal<<"FnBPRSSBTotal:\t"<<FnBPRSSBTotal<<endl;
          FvNASTemp=elem_prod(FvBPRNAS,FmBPRS(k));                           //die        N(t)=(N(t-1)*exp(-Z)*G)
          FvBPRNAS= FvNASTemp*Fd3GrowthMatrix(k);   
       }
  }
  //cout<<"FnBPRYieldTotal BBB:"<<FnBPRYieldTotal<<"FnBPRSSBTotal"<<FnBPRSSBTotal<<endl;
  return ;
}

void model_parameters::GetSPRAndYPR_Y(double &FnBPRSSBTotal,double &FnBPRYieldTotal,double dvFtot,int iWeightFlag,int iModelFlag)
{
   double   dR;
   double   dvF;
   dmatrix  FmBPRF(1,Flagtimestep,1,DiSizeBinNum);
   dmatrix  FmBPRZ(1,Flagtimestep,1,DiSizeBinNum);
   dmatrix  FmBPRS(1,Flagtimestep,1,DiSizeBinNum);
   dvector  FvBPRSSB_S(1,DiSizeBinNum);
   //d3array 
   //d3array  Fd3GrowthMatrix(1,DiSeasonNum,1,DiSizeBinNum,1,DiSizeBinNum);
   dvector  FvBPRNAS(1,DiSizeBinNum);
   dvector  FvNASTemp(1,DiSizeBinNum);
   //dmatrix  dmFRInFleetSeason(1,DiSeasonNum,1,DiFleetNum);
   //dvector  dvFRInSeason(1,DiSeasonNum);
   int k,i,kk,iL;
   FnBPRYieldTotal=0.0;
   FnBPRSSBTotal=0.0;
   dR=1.0;//unit recruitment because we assume the balance was reached!
   //dvector  dvFRInSeason(1,DiSeasonNum)
   for(k=1;k<=Flagtimestep;k++)
   {
       dvF=dvFtot*DvFRInSeason(k);
       if(DiNaturalMortalityFlag==2)
       {
                iL= DivBPRBlockGrowthMatrix(k);
                iL=(iL-1)*DiSeasonNum+k;
       }
       else
       {
                iL=DivBPRBlockGrowthMatrix(k);
       }
       FmBPRF(k).initialize();
       for(kk=1;kk<=DiSizeBinNum;kk++)
       {
          if(DiBPRSelectivityFlag<=DiFleetNum)
                 FmBPRF(k,kk)=dvF*value(MmBPRSel(k,kk));
          else
          {
             for(i=1;i<=DiFleetNum;i++)
             {
                 FmBPRF(k,kk)+=dvF* dmFRInFleetSeason(k,i)*value(Md4SelByFleet(i,k,DiBPRCalYear,kk));
             }  
          }
           FmBPRZ(k,kk)=FmBPRF(k,kk)+value(Md3M(k,DiBPRCalYear,kk));
           for(i=1;i<=DiSizeBinNum;i++)
           {
              Fd3GrowthMatrix(k,kk,i)= value(Md3GrowthMatrix(iL,kk,i));
           }
       }
       FmBPRS(k)=mfexp(-FmBPRZ(k));
       FvBPRSSB_S=mfexp(-FmBPRZ(k)*DnFracYearSSB);
   }
   //cout<<FmBPRS<<endl;
   //cout<<FvBPRSSB_S<<endl;
   //
   for(i=1;i<=DiBPRPeriod;i++)
   {
      for(k=1;k<=Flagtimestep;k++)
      {
          if( iModelFlag)//all cohorts provide yield and ssb in one year 
          {
                  //cout<<"Slow Here!"<<endl;
                 // for (kk=1;kk<=DiSizeBinNum;kk++)  
                 //        MvBPRNAS(kk)+= dR*DvRecruitmentSeaRate(k)*MvRecruitPrjVect(kk);
                   FvBPRNAS+= dR*DvRecruitmentSeaRate(k)*value(MvRecruitPrjVect);
                // CalCatchInBPR(DiBPRCalYear,k,dvR,iWeightFlag);
          }
          else
          {//each cohort can provide yield and ssb
               if(i==1)
                  FvBPRNAS+= dR*DvRecruitmentSeaRate(k)*value(MvRecruitPrjVect);
               if(iWeightFlag)
               {
                    if(k==DiSBBSeason)
                        FnBPRSSBTotal+=elem_prod(FvBPRNAS,FvBPRSSB_S)*DmFecundity(DiBPRCalYear); 
                    FnBPRYieldTotal+=elem_prod(elem_div(FmBPRF(k),FmBPRZ(k)),elem_prod(FvBPRNAS,1.0-FmBPRS(k)))*DmWeightAtSize(DiBPRCalYear);
               }
               else
               {
                    if(k==DiSBBSeason)
                        FnBPRSSBTotal+=sum(elem_prod(FvBPRNAS,FvBPRSSB_S)); 
                    FnBPRYieldTotal+=elem_prod(elem_div(FmBPRF(k),FmBPRZ(k)),elem_prod(FvBPRNAS,1.0-FmBPRS(k)))*DmWeightAtSize(DiBPRCalYear);
               }
          }
          //ofRuntimelog<<"Yield "<<MnBPRYieldTotal<<"  "<<MnBPRYield<<endl;
          //ofRuntimelog<<"NAS: "<<MvBPRNAS<<endl;
          FvNASTemp=elem_prod(FvBPRNAS,FmBPRS(k));                                    //die        N(t)=(N(t-1)*exp(-Z)*G)
          FvBPRNAS= FvNASTemp*Fd3GrowthMatrix(k);           //Growth 
      }
  }
  //cout<<FvBPRNAS<<endl;
  if( iModelFlag)
  {
      // if(dvFtot>1.95||dvFtot<0.001)
      //     ofRuntimelog<<"F:\t"<<dvFtot<<"FvBPRNAS:"<<endl<<FvBPRNAS<<endl;
       for(k=1;k<=Flagtimestep;k++)
       {
           FvBPRNAS+= dR*DvRecruitmentSeaRate(k)*value(MvRecruitPrjVect);
           if(iWeightFlag)
           {
                    FnBPRSSBTotal+=elem_prod(FvBPRNAS,FvBPRSSB_S)*DmFecundity(DiBPRCalYear); 
                    FnBPRYieldTotal+=elem_prod(elem_div(FmBPRF(k),FmBPRZ(k)),elem_prod(FvBPRNAS,1.0-FmBPRS(k)))*DmWeightAtSize(DiBPRCalYear);
           }
           else
           {
                    FnBPRSSBTotal+=sum(elem_prod(FvBPRNAS,FvBPRSSB_S)); 
                    FnBPRYieldTotal+=elem_prod(elem_div(FmBPRF(k),FmBPRZ(k)),elem_prod(FvBPRNAS,1.0-FmBPRS(k)))*DmWeightAtSize(DiBPRCalYear);
           }
          if(dvFtot<0.133 && dvFtot>0.132)
              ofRuntimelog<<"F:\t"<<dvFtot<<" FnBPRYieldTotal\t:"<< FnBPRYieldTotal<<"FnBPRSSBTotal:\t"<<FnBPRSSBTotal<<endl;
          FvNASTemp=elem_prod(FvBPRNAS,FmBPRS(k));                           //die        N(t)=(N(t-1)*exp(-Z)*G)
          FvBPRNAS= FvNASTemp*Fd3GrowthMatrix(k);   
       }
  }
  //cout<<"FnBPRYieldTotal BBB:"<<FnBPRYieldTotal<<"FnBPRSSBTotal"<<FnBPRSSBTotal<<endl;
  return ;
}

double model_parameters::MSY_Sub(double dF,double dDelta,int iWeightFlag,int iModelFlag,int iRS)
{
   double dY1, dY2,dR;
   double dalpha,dbeta;
   double FnBPRSSBTotal,FnBPRYieldTotal;
   if(iRS==1)
   {
        dalpha=mfexp(value(PbnvLogRSPara(1)));
        dbeta=0.0;
   }
   else
   {
      dalpha=mfexp(value(PbnvLogRSPara(1)));
      dbeta=mfexp(value(PbnvLogRSPara(2)));
   }
   FnBPRSSBTotal=FnBPRYieldTotal=0.0;
   //MnTemp1=dF;
   if (Flagtimestep==4)
       GetSPRAndYPR(FnBPRSSBTotal,FnBPRYieldTotal,dF, iWeightFlag, iModelFlag);
   if (Flagtimestep==1)
       GetSPRAndYPR_Y(FnBPRSSBTotal,FnBPRYieldTotal,dF, iWeightFlag, iModelFlag);
  // pt=log(value( nBPRYieldTotal)+1.0);
   dY1=CalBPRCatch(dalpha,dbeta,FnBPRSSBTotal,FnBPRYieldTotal,DiRSFlag);
  // cout<<"FnBPRYieldTotal 1:"<<FnBPRYieldTotal<<"FnBPRSSBTotal"<<FnBPRSSBTotal<<endl;
    ofRuntimelog<<"DY1:"<<dY1<<"\tdF\t"<<dF<<"\t"<<dalpha<<"\t"<<dbeta <<"SPR:\t"<<FnBPRSSBTotal<<"YPR:\t"<<FnBPRYieldTotal<<endl;
    if(dY1<0)
       dY1=FnBPRYieldTotal*dalpha;
  // if(iRS!=3)
  // {
 //     dR=value(nAlpha*nBPRSSBTotal-nBelta)/value(nBPRSSBTotal);
 //     dY1=dR*value( nBPRYieldTotal);
 //     
 //  }
 //  else
 //  {
 //        dR=value(log(nAlpha)+log(nBPRSSBTotal))/value(nBelta);
 //        dR/=value(nBPRSSBTotal);
 //        dY1=dR*value( nBPRYieldTotal);
 //  }
   dF=dF+ dDelta;
   if (Flagtimestep==4)
       GetSPRAndYPR(FnBPRSSBTotal,FnBPRYieldTotal,dF, iWeightFlag, iModelFlag);
   if (Flagtimestep==1)
       GetSPRAndYPR_Y(FnBPRSSBTotal,FnBPRYieldTotal,dF, iWeightFlag, iModelFlag);
   dY2=CalBPRCatch(dalpha,dbeta,FnBPRSSBTotal,FnBPRYieldTotal,DiRSFlag);
   //cout<<"FnBPRYieldTotal 2:"<<FnBPRYieldTotal<<"FnBPRSSBTotal"<<FnBPRSSBTotal<<endl;
   //dY2=CalBPRCatch(dalpha,dbeta,FnBPRSSBTotal,FnBPRYieldTotal,DiRSFlag);
    ofRuntimelog<<"DY2:"<<dY2<<"\tdF"<<dF<<" "<<dalpha<<" "<<dbeta <<"SPR:\t"<<FnBPRSSBTotal<<"YPR:\t"<<FnBPRYieldTotal<<endl;
   //ofRuntimelog<<"DY2:"<<dY2<<" "<<MnTemp1<<endl;
   // if(iRS!=3)
  // {
  //    dR=value(nAlpha*nBPRSSBTotal-nBelta)/value(nBPRSSBTotal);
  //    dY2=dR*value( nBPRYieldTotal);
  // }
  // else
  // {
  //       dR=value(log(nAlpha)+log(nBPRSSBTotal))/value(nBelta);
  //       dR/=value(nBPRSSBTotal);
  //       dY2=dR*value( nBPRYieldTotal);
  // }
  // dpt=log(value( nBPRYieldTotal)+1.0);
  // nTemp1=nAlpha*nBPRSSBTotal-nBelta;
    if(dY2<0)
       dY2=FnBPRYieldTotal*dalpha;
   dR=(dY2-dY1);///dDelta;//dY2-dY1;//(dY2-dY1)/dDelta;
   return dR;
}

double model_parameters::MSY_Bisection(double dF1, double dF2, double dDelta,int iWeightFlag,int iModelFlag,int iRS)
{
   double bf,fmid;
   int j;
   double  xmid,rtb,dx;
   bf=MSY_Sub(dF1,dDelta/3.0, iWeightFlag, iModelFlag,iRS);
   fmid=MSY_Sub(dF2,dDelta/3.0, iWeightFlag, iModelFlag,iRS);
   ofRuntimelog<<"bf:\t"<<bf<<endl;
   ofRuntimelog<<"fmid:\t"<<fmid<<endl;
   if(fabs(bf)<DIGITALPREC)
      return dF1;
   if(fabs(fmid)<DIGITALPREC)
      return dF2;
   if (bf*fmid > 0.0)
   {
            ofRuntimelog<<"Root must be bracketed for bisection in rtbis of MSY Double"<<endl;
   }
   else
   {
             rtb = bf < 0.0 ? (dx=dF1-dF2,dF2): (dx=dF2-dF1,dF1); //Orient the search so that f>0
             for (j=1;j<=20;j++) 
             { //lies at x+dx.
                     dx = dx *0.5;
                     xmid=rtb+dx;
                     fmid=MSY_Sub(xmid,dDelta/3.0, iWeightFlag, iModelFlag,iRS); //Bisection loop.
                     //cout<<"MSYtestF = "<<xmid<<" Dyield = "<<fmid<<endl;
                     if(fabs(fmid)<DIGITALPREC)
                     {
                            rtb=xmid;
                            return rtb;
                     }
                     if (fabs(dx) < dDelta) 
                           return rtb;
                     if (fmid > 0.0) 
                            rtb=xmid;
             }
               ofRuntimelog<<"Too many bisections in rtbis"<<endl;
   }
   return 0.0;
}

double model_parameters::SPR_Bisection(double dF1,double dF2, double dDelta, double dTargetegg,int iWeightFlag,int iModelFlag)
{
    double  dx, bf, fmid, xmid, rtb;
    double FnBPRSSBTotal,FnBPRYieldTotal;
    int j;
    //MnTemp1.initialize();
    //MnTemp1=dF1;
    if (Flagtimestep==4)
        GetSPRAndYPR(FnBPRSSBTotal,FnBPRYieldTotal,dF1, iWeightFlag, iModelFlag);
    if (Flagtimestep==1)
        GetSPRAndYPR_Y(FnBPRSSBTotal,FnBPRYieldTotal,dF1, iWeightFlag, iModelFlag);
    bf=dTargetegg-FnBPRSSBTotal;//value(MnBPRSSBTotal);
    //MnTemp1=dF2;
    if (Flagtimestep==4)
        GetSPRAndYPR(FnBPRSSBTotal,FnBPRYieldTotal,dF2, iWeightFlag, iModelFlag);
    if (Flagtimestep==1)
        GetSPRAndYPR_Y(FnBPRSSBTotal,FnBPRYieldTotal,dF2, iWeightFlag, iModelFlag);
    fmid=dTargetegg-FnBPRSSBTotal;
    if(fabs(bf)<DIGITALPREC)
      return dF1;
    if(fabs(fmid)<DIGITALPREC)
      return dF2;
    if (bf*fmid > 0.0)
    {
             ofRuntimelog<<"Root must be bracketed for bisection in rtbis of F10%"<<endl;
    }
    else
    {
            rtb = bf < 0.0 ? (dx=dF2-dF1,dF1) : (dx=dF1-dF2,dF2); //Orient the search so that f>0
            for (j=1;j<=20;j++) { //lies at x+dx.
                 dx = dx *0.5;
                 xmid=rtb+dx;
                 //MnTemp1=xmid;
                 if (Flagtimestep==4)
                     GetSPRAndYPR(FnBPRSSBTotal,FnBPRYieldTotal,xmid, iWeightFlag, iModelFlag);
                 if (Flagtimestep==1)
                     GetSPRAndYPR_Y(FnBPRSSBTotal,FnBPRYieldTotal,xmid, iWeightFlag, iModelFlag);
                 fmid=dTargetegg-FnBPRSSBTotal;//value(MnBPRSSBTotal); //Bisection loop.
                 if(fabs(fmid)<DIGITALPREC)
                 {
                            rtb=xmid;
                            return rtb;
                 }
                 if (fabs(dx) < dDelta) 
                           return rtb;
                 if(fmid < 0.0) 
                           rtb=xmid;
           }
           ofRuntimelog<<"Too many bisections in rtbis"<<endl;
    }
     return 0.0;
}

double model_parameters::F01_Bisection(double dF1,double dF2, double dDelta, double dTargetegg,int iWeightFlag,int iModelFlag)
{
   double  dSlopeOrigin;
   double  dA,dB,dC;
   double  dS1,dS2;
   double  dR;
   double  dSlopeTarget;
   double  dnTemp1;
   double  dSlopeCurrent;
   double FnBPRSSBTotal,FnBPRYieldTotal;
   int     i;
   //MnTemp1.initialize();
   //MnTemp1=dDelta;
   if (Flagtimestep==4)
       GetSPRAndYPR(FnBPRSSBTotal,FnBPRYieldTotal,dDelta, iWeightFlag, iModelFlag);
   if (Flagtimestep==1)
       GetSPRAndYPR_Y(FnBPRSSBTotal,FnBPRYieldTotal,dDelta, iWeightFlag, iModelFlag);
   dSlopeOrigin=FnBPRYieldTotal/dDelta;
   dA=dF1;
   dB=dF2;
   dSlopeTarget= dTargetegg*dSlopeOrigin;
   for (i=1;i<=20;i++)
   {
      dC=(dA+dB)/2.0;
      dnTemp1=dC+dDelta;
      if (Flagtimestep==4)
          GetSPRAndYPR(FnBPRSSBTotal,FnBPRYieldTotal,dnTemp1, iWeightFlag, iModelFlag);
      if (Flagtimestep==1)
          GetSPRAndYPR_Y(FnBPRSSBTotal,FnBPRYieldTotal,dnTemp1, iWeightFlag, iModelFlag);
      dS1= FnBPRYieldTotal;
      dnTemp1=dC;
      if (Flagtimestep==4)
          GetSPRAndYPR(FnBPRSSBTotal,FnBPRYieldTotal,dnTemp1, iWeightFlag, iModelFlag);
      if (Flagtimestep==1)
          GetSPRAndYPR_Y(FnBPRSSBTotal,FnBPRYieldTotal,dnTemp1, iWeightFlag, iModelFlag);
      dS2=FnBPRYieldTotal;
      dSlopeCurrent=(dS1-dS2)/dDelta;
      if (dSlopeCurrent <  dSlopeTarget)
           dB=dC;
      else
           dA=dC;
    }
    dR=dC;
    return dR;
}

void model_parameters::CalBiologyRP(void)
{
   double  dFinterval=5e-03;
   //dvariable  dvFMin,dvFMax;
   int     iWeightFlag, iModelFlag;
   dvariable  dalpha,dBeta;
   double  dFMin,dFMax;
   double  FnBPRSSBTotal,FnBPRYieldTotal;
   dalpha=0.0;
   dBeta=0.0;
   //dvFMin=0.0;
   //dvFMax=0.0;
   sdr_nBPRMSY.initialize();
   MnBPRFmsy.initialize();
   MnBPRF30SPR.initialize();
   MnBPRF40SPR.initialize();
   MnBPRF01.initialize();
   MnBPRFMax.initialize();
   MnBPRSSBMSY.initialize();
   if(DiRSFlag==1)
   {
        dalpha=mfexp(PbnvLogRSPara(1));
        dBeta=0.0;
        //sdr_nBPRMSY=mfexp(PbnvLogRSPara(1));
   }
   else
   {
      dalpha=mfexp(PbnvLogRSPara(1));
      dBeta=mfexp(PbnvLogRSPara(2));
   }
   iWeightFlag=1;
   iModelFlag=1;
   dFMin=0.0;
   dFMax= DnFmax;
   //ofRuntimelog<<" "<<DvFRInSeason<<endl;
   ///
   dFMin =MSY_Bisection(dFMin,dFMax, dFinterval,iWeightFlag,iModelFlag,DiRSFlag);
   MnBPRFmsy=dFMin;
   cout<<"F_msy:"<<dFMin<<endl;
   //MnTemp1=MnBPRFmsy;
   //GetSPRAndYPR(FnBPRSSBTotal,FnBPRYieldTotal,dFMin, iWeightFlag, iModelFlag);
   if (Flagtimestep==4)
       GetSPRAndYPR(MnBPRFmsy,iWeightFlag,iModelFlag);
   if (Flagtimestep==1)
       GetSPRAndYPR_Y(MnBPRFmsy,iWeightFlag,iModelFlag);
   cout<<"YPR_fmsy:"<<MnBPRYieldTotal<<endl;
   cout<<"SPR_fmsy:"<<MnBPRSSBTotal<<endl;
   MnBPRBiomassMSY=MvBPRNAS*DmWeightAtSize(DiBPRCalYear);//biomass produce MSY
   //cout<<MvBPRNAS<<endl;
   cout<<"Bmsy:"<<MnBPRBiomassMSY<<endl;
   // cout<<"Slow"<<FnBPRSSBTotal<<endl;
   //cout<<"Slow 0:"<<endl;
   sdr_nBPRMSY=CalBPRCatch(dalpha,dBeta,MnBPRSSBTotal,MnBPRYieldTotal,DiRSFlag);//MSY
   //cout<<dalpha<<endl;
   cout<<"MSY:"<<sdr_nBPRMSY<<endl;
   MnBPRBiomassMSY*=sdr_nBPRMSY/MnBPRYieldTotal;//???
   MnBPRSSBMSY=MnBPRSSBTotal*sdr_nBPRMSY/MnBPRYieldTotal;
   cout<<"SSBmsy:"<<MnBPRSSBMSY<<endl;
   //cout<<"Slow 1:"<<endl;
   //exit(-99);
   //cout<<FnBPRSSBTotal<<endl;
   //cout<<FnBPRYieldTotal<<endl;
   if (Flagtimestep==4)
       GetSPRAndYPR(FnBPRSSBTotal,FnBPRYieldTotal,0.0, iWeightFlag, iModelFlag);
   if (Flagtimestep==1)
       GetSPRAndYPR_Y(FnBPRSSBTotal,FnBPRYieldTotal,0.0, iWeightFlag, iModelFlag);
   nTemp=FnBPRSSBTotal;//value(MnBPRSSBTotal);
   cout<<"SPR_f0:"<<nTemp<<endl;
   MnBPRF30SPR=SPR_Bisection(0.0, DnFmax, dFinterval, nTemp*0.3,iWeightFlag,iModelFlag);//search F30
   cout<<"F30%:"<<MnBPRF30SPR<<endl;
   MnBPRF40SPR=SPR_Bisection(0.0, DnFmax, (dFinterval), nTemp*0.4,iWeightFlag,iModelFlag);//search F40
   cout<<"F40%:"<<MnBPRF40SPR<<endl;
   /////////////////////////////////////////////////////////////////////////////////////////////////////
   MnBPRF01=F01_Bisection(0.0, DnFmax, (dFinterval),0.1,iWeightFlag,iModelFlag);
   cout<<"F0.1:"<<MnBPRF01<<endl;
   MnBPRFMax=F01_Bisection(0.0, DnFmax, (dFinterval),0.0,iWeightFlag,iModelFlag);
   cout<<"Fmax:"<<MnBPRFMax<<endl;
   //*/
}

dvariable model_parameters::CalBPRCatch(const dvariable &nAlpha,const dvariable &nBelta,const dvariable &dvSPR,const dvariable &dvYPR,int iRS)
{
    dvariable dvR;
    dvariable dvTemp;
    dvR=0.0;
    switch(iRS)
    {
       case 1:
              dvR=nAlpha*dvYPR;
              break;
       case 2:
              dvTemp=(nAlpha*dvSPR-nBelta)/dvSPR;
              dvR=  dvTemp*dvYPR;
              break;
       case 3:
              dvTemp=(log(nAlpha)+log(dvSPR))/nBelta;
              dvTemp/=dvSPR;
              dvR =  dvTemp* dvYPR;
              break;
   } 
  return dvR;
}

double model_parameters::CalBPRCatch(const dvariable &nAlpha,const dvariable &nBelta,const double  &dSPR, const double &dYPR, const int iRS)
{
    double dR;
    double dTemp;
    dR=0.0;
    switch(iRS)
    {
            case 1:
                dR=value(nAlpha)*dYPR;
                break;
            case 2:
                dTemp=value(nAlpha*dSPR-nBelta)/dSPR;
                dR=  dTemp*dYPR;
                break;
            case 3:
                dTemp=value((log(nAlpha)+log(dSPR))/nBelta);
                dTemp/=dSPR;
                dR =  dTemp* dYPR;
                break;
    } 
    return dR;
}

dvariable model_parameters::CalBPRCatch(const dvariable &nAlpha,const dvariable &nBelta,const double  &dSPR, const double &dYPR, const int iRS,int iFlag)
{
    dvariable  dR;
    dvariable  dTemp;
    dR=0.0;
    switch(iRS)
    {
            case 1:
                dR=nAlpha*dYPR;
                break;
            case 2:
                dTemp=(nAlpha*dSPR-nBelta)/dSPR;
                dR=  dTemp*dYPR;
                break;
            case 3:
                dTemp=((log(nAlpha)+log(dSPR))/nBelta);
                dTemp/=dSPR;
                dR =  dTemp* dYPR;
                break;
    } 
    return dR;
}

double model_parameters::CalBPRCatch(double nAlpha,double nBelta,const double  &dSPR, const double &dYPR, const int iRS)
{
    double dR;
    double dTemp;
    dR=0.0;
    switch(iRS)
    {
            case 1:
                dR=nAlpha*dYPR;
                break;
            case 2:
                dTemp=(nAlpha*dSPR-nBelta)/dSPR;
                dR=  dTemp*dYPR;
                break;
            case 3:
                dTemp=((log(nAlpha)+log(dSPR))/nBelta);
                dTemp/=dSPR;
                dR =  dTemp* dYPR;
                break;
    } 
    return dR;
}

void model_parameters::WritePAR(void)
{
  int i,j,k,ii,lL2;
  dvariable dvAlpha1;
  dvariable dvBelta1;
  dvariable dvAlpha2;
  dvariable dvBelta2;  
  dvariable dvarR;
  ofstream  ofPAROutput("./PAROUTPUT.dat");
   for(ii=1;ii<=DiFleetSelBlockNum;ii++)
   {     
                   ofPAROutput<<"#Catch Selectivity Block:"<<ii<<endl; 
                   lL2=DivFleetSelParaIndex_Ini(ii);                 //get location of parameter
                   switch(DivFleetSelBlockFlag(ii))             // whether the selectivity was same as fleet no. x x in ivFleetSelBlockFlag(lL1)
                   {
                     case 1:
                        for(kk=1;kk<=DiSizeBinNum;kk++)
                                ofPAROutput<<PbnvFleetSelParams(lL2+kk)<<"\t";
                        ofPAROutput<<endl;
                        break;
                     case 2:
                       dvAlpha1= PbnvFleetSelParams(lL2+1);
                       dvBelta1= PbnvFleetSelParams(lL2+2);//dvBelta1=1.0/ PbnvFleetSelParams(lL2+2);
                       ofPAROutput<< dvAlpha1<<"\t"<< dvBelta1<<endl;
                       break;
                     case 3:
                         dvAlpha1= PbnvFleetSelParams(lL2+1);
                         dvBelta1= PbnvFleetSelParams(lL2+2);//dvBelta1=1.0/ PbnvFleetSelParams(lL2+2);
                         dvAlpha2= PbnvFleetSelParams(lL2+3);
                         dvBelta2= PbnvFleetSelParams(lL2+4);// dvBelta2=1.0/ PbnvFleetSelParams(lL2+4);
                         ofPAROutput<< dvAlpha1<<"\t"<< dvBelta1<<"\t"<< dvAlpha2<<"\t"<< dvBelta2<<endl;
                         break;                  
                    }
  }
  ///////////////////////////
  for(i=1;i<=DiFleetNum;i++)
  {
      ofPAROutput<<"#Fishing mortality in first year for each season"<<endl;
      for(k=1;k<=Flagtimestep;k++)
      {
        ofPAROutput<<PbnvLogFYear1Season1((i-1)*DiSeasonNum+k)<<"\t";
      }
      ofPAROutput<<endl;
  }
  /////////////////////////////////////////////
  for(i=1;i<=DiFleetNum;i++)
  {
      ofPAROutput<<"#Fleet:"<<i<<endl;
      for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
      {
        for(k=1;k<=Flagtimestep;k++)
        {
            iL2=(i-1)*(DiCalEndYear-DiCalBeginYear+1)*DiSeasonNum+(j-DiCalBeginYear)*DiSeasonNum+k;;
            ofPAROutput<<PbnvLogFDevs(iL2)<<"\t";
        }
         ofPAROutput<<endl;
      }
  }
  ofPAROutput<<"#CPUE Block Q E1"<<endl;
  for(i=1;i<=DiCPUEBlockNum;i++)
     ofPAROutput<<PbnvCPUEQBiomassE1(i)<<"\t";
  ofPAROutput<<endl;
  //
  for(jj=1;jj<=DiIndexSelBlockNum;jj++)
  {
              lL2=DivIndexSelParaIndex_Ini(jj);   
              ofPAROutput<<"#Index Selectivity Block:"<<jj<<endl; 
              switch(DivIndexSelFlag_Ini(jj))
              {
               case 1:
                  for(k=1;k<=DiSizeBinNum;k++)
                       ofPAROutput<<PbnvIndexSelParams(lL2+k)<<endl;
                  ofPAROutput<<endl;
                  break;
               case 2:
                   dvAlpha1= PbnvIndexSelParams(lL2+1);
                   dvBelta1= PbnvIndexSelParams(lL2+2);//dvBelta1=1.0/ PbnvIndexSelParams(lL2+2);
                   ofPAROutput<<dvAlpha1<<"\t"<<dvBelta1<<endl;
                   //nTemp=
                   break;
               case 3:
                   dvAlpha1= PbnvIndexSelParams(lL2+1);
                   dvBelta1= PbnvIndexSelParams(lL2+2); //dvBelta1=1.0/ PbnvIndexSelParams(lL2+2);
                   dvAlpha2= PbnvIndexSelParams(lL2+3);
                   dvBelta2= PbnvIndexSelParams(lL2+4);//dvBelta2=1.0/ PbnvIndexSelParams(lL2+4);
                   ofPAROutput<<dvAlpha1<<"\t"<<dvBelta1<<"\t"<<dvAlpha2<<"\t"<<dvBelta2<<endl;
                   break;
           }
   }
  //
   ofPAROutput<<"#Survey Index Block Q E1"<<endl;
   for(i=1;i<=DiIndexQBlockNum;i++)
      ofPAROutput<<PbnvIndexQBiomassE1(i)<<"\t";
   ofPAROutput<<endl;
  //
   ofPAROutput<<"#Initial abundance parameters3"<<endl;
   for(i=1;i<=DiNYear1EstParaNum;i++)
      ofPAROutput<<PbvNYear1Para(i)<<"\t";
   ofPAROutput<<endl;
  //
  ofPAROutput<<"#Recruitment-Spawnning Relationship Parameters"<<endl;
   for(i=1;i<=2;i++)
      ofPAROutput<<PbnvLogRSPara(i)<<"\t";
   ofPAROutput<<endl;
  //
   ofPAROutput<<"#Recruitment deviation log scale"<<endl;
   for(i=DiCalBeginYear;i<=DiCalEndYear-DiYearBeforeEndForRDev;i++)
      ofPAROutput<<PbdvLogRecruitDevs(i)<<"\t";
   for(i=DiCalEndYear-DiYearBeforeEndForRDev+1;i<=DiCalEndYear;i++)
      ofPAROutput<< PbdvLogRecruitDevsADD(i)<<"\t";
   ofPAROutput<<endl;
  //ofRuntimelog
  //ofRuntimelog<<"#Recruitment autocorrelation Coefficient"<<endl;
  //ofRuntimelog<<PbnvR00Para<<endl;
  ofPAROutput<<"#Recruitment autocorrelation Coefficient"<<endl;
  ofPAROutput<<PbnRecruitmentRh<<endl;
  ofPAROutput<<"##Recruitment deviation (log scale) SD "<<endl;
  //ofPAROutput<<PbnvRecruitLogDevsSD2<<endl;
  ofPAROutput<<MnRSD<<endl;
  //ofPAROutput<<"##Recruitment R00"<<endl;
  //ofPAROutput<<PbnvR00Para<<endl;
   ofPAROutput<<"#Abundance for season one at first year"<<endl;
   GetNYear1Season1();
   for(i=1;i<=DiSizeBinNum;i++)
   {
      dvarR= MvNASTemp(i) +sdr_vRecruitmentP(DiCalBeginYear)*DvRecruitmentSeaRate(1)*MvRecruitPrjVect(i);
      ofPAROutput<<MvNASTemp(i)<<"\t"<<sdr_vRecruitmentP(DiCalBeginYear)<<"\t"<<DvRecruitmentSeaRate(1)<<"\t"<<MvRecruitPrjVect(i)<<"\t"<<dvarR<<endl;
   }
   ofPAROutput<<endl;
  ofPAROutput<<"#Natural Mortality Coefficient A"<<endl;
  for(i=1;i<=Flagtimestep;i++)
     ofPAROutput<<PbnvLorenzenA(i)<<"\t";
  ofPAROutput<<endl;
  ofPAROutput<<"#Natural Mortality Coefficient B"<<endl;
  for(i=1;i<=Flagtimestep;i++)
     ofPAROutput<<PbnvLorenzenB(i)<<"\t";
  ofPAROutput<<endl;
  ofPAROutput<<"#Recruitment Project vector for size was estimated "<<endl;
  for(i=1;i<=DiSizeBinNum;i++)
    ofPAROutput<<PbnvRecruitPrjVect(i)<<"\t";
  ofPAROutput<<endl;
  ofPAROutput.close();
}

void model_parameters::WriteCohortTrack(int iRecruitmentY)
{
    ofstream  ofCohortTrackOutput("./DebugFile/CohortTrack.dat");
    dvar_vector    dvAbundance(1,DiSizeBinNum);
    int j,i,k,iL1;
    dvAbundance.initialize();
   /////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    if(iRecruitmentY>=DiCalBeginYear &&  iRecruitmentY<=DiCalEndYear)
    {
       ////////////////////////////////////////////////////////
       ofCohortTrackOutput<<"Year\t"<<"Season\t";
       for (k=1;k<=DiSizeBinNum-1;k++)
                    ofCohortTrackOutput<<"SizeBin"<<k<<"\t";
       ofCohortTrackOutput<<"SizeBin"<<k<<endl;
       //////////////////////////////////////////////////////////
       for(i=iRecruitmentY;i<=DiCalEndYear;i++)
       {
          for(j=1;j<=Flagtimestep;j++)
          {  
             if(i==iRecruitmentY)
             {
                 for (k=1;k<=DiSizeBinNum;k++)
                 {
                    dvAbundance(k)= dvAbundance(k)+sdr_vRecruitmentP(i)*DvRecruitmentSeaRate(j)*MvRecruitPrjVect(k);
                 }
             }
            ////////////////////////////////////////////////////
             ofCohortTrackOutput<<i+DiBeginYear<<"\t"<<j<<"\t";
             for (k=1;k<=DiSizeBinNum-1;k++)
                    ofCohortTrackOutput<<dvAbundance(k)<<"\t";
             ofCohortTrackOutput<<dvAbundance(k)<<endl;
            //////////////////////////////////////////////////////////
             MvNASTemp=elem_prod(dvAbundance,Md3S(j,i));//die 
             iL1=DimGrowthMatrixBlockFlag(i,j+1);
             if(DiGrowthMatrixFlag==2)
                  iL1=(iL1-1)*DiSeasonNum+j;
             if((DiGrowthMatrixFlag!=2 && (iL1<1||iL1>DiGrowthMatrixBlockNum))||(DiGrowthMatrixFlag==2&&(iL1<1||iL1>DiGrowthMatrixBlockNum*DiSeasonNum)))
             {
                  WARNING<<"No Growth Matrix was used in WriteCohortTrack"<<endl;
                  dvAbundance= MvNASTemp;                   //Growth Next Season
             }
             else
             {
                 dvAbundance= MvNASTemp*Md3GrowthMatrix(iL1);           //Growth Next Season
             } 
         }
      }
    }
    else
    {
        ofCohortTrackOutput<<"Error: the input year less than DiCalBeginYear or larger than  DiCalEndYear"<<endl;
    }
    ofCohortTrackOutput.close();
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  if(last_phase())
  {
        GetBPRSel();
        CalBiologyRP();
        if(DiCohortTrackingBeginYear>0)
        {
             iTemp= DiCohortTrackingBeginYear-DiBeginYear+1;
             WriteCohortTrack(iTemp);
        } 
  }
  //report << "#Norther_Shrimp_length_structure_model_01" << endl;
  //report << "#Beginning_time_for_run:" << ctime(&tBeginTime) << endl;
  // Likelihood
  report << "Objective_function_value" << endl << ofvTotal << endl;
  report<<"BeginYear"<<endl<<iBeginYear<<endl;
  report<<"EndYear"<<endl<<iEndYear<<endl;
  report<<"N_year"<<endl<<iEndYear-iBeginYear+1<<endl;
  report<<"Time_step"<<endl<<itimestep<<endl;
  report<<"Size_bins"<<endl<<DvSizeBinData<<endl;
  report<<"Weight_length"<<endl<<DmWeightAtSize<<endl;
  report<<"Maturity_length"<<endl<<DmMaturityAtSize<<endl;
  report<<"Growth"<<endl<<Md3GrowthMatrix<<endl;
  //for(i=1;i<=(iEndYear-iBeginYear+1)*itimestep;i++)
  //   {
  //    report<<"Growth"<<i<<endl<<Md3GrowthMatrix(i)<<endl;
  //   }
  report<<"Growth_ratio"<<endl<<DvGrowthTimeAsYear<<endl;
  report<<"Recrut_Season_ratio"<<endl<<DvRecruitmentSeaRate<<endl;
  report<<"Recrut_Sizebin_ratio"<<endl<<MvRecruitPrjVect<<endl;
   for(i=1;i<=DiFleetNum;i++)
   {
     for(k=1;k<=Flagtimestep;k++)
     {
         Md3FleetEffectiveSampleSize(i,k)=0.0;
         for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
         {
            MnTemp1=norm2(Md4CatchPropAtSizeP(i,k,j)(DimCSelStartSizeBin(i,k),DimCSelEndSizeBin(i,k))-Dd4CatchPropAtSizeO(i,k,j)(DimCSelStartSizeBin(i,k),DimCSelEndSizeBin(i,k)));
           if(MnTemp1>0.0)
           {
             Md3FleetEffectiveSampleSize(i,k,j)=Md4CatchPropAtSizeP(i,k,j)(DimCSelStartSizeBin(i,k),DimCSelEndSizeBin(i,k))*(1.0-Md4CatchPropAtSizeP(i,k,j)(DimCSelStartSizeBin(i,k),DimCSelEndSizeBin(i,k)));
             Md3FleetEffectiveSampleSize(i,k,j)/=MnTemp1;
           }
         }
      }
   }
   report<<"BiomassORNum"<<endl<<DimCatchBiomassFlag<<endl;
   //Catch(fleet*season*year)
   report<<"Catch_Obs"<<endl<<Dd3CatchTotalO<<endl;
   report<<"Catch_Pred"<<endl<<Md3CatchTotalP<<endl;
   // Catch ESS (fleet*season*year)
   report<<"Catch_ESSinput"<<endl<<Dd3CatchESSInput<<endl;
   report<<"Catch_ESSpred"<<endl<<Md3FleetEffectiveSampleSize<<endl;
   // Catch Composition
   report<<"Catch_Comp_Obs"<<endl<<Dd4CatchPropAtSizeO<<endl;
   report<<"Catch_Comp_Pred"<<endl<<Md4CatchPropAtSizeP<<endl;
  // Index Obs and Pred  
   if(DiAvailIndexNum>0)
       MmIndexEffectiveSampleSize=0.0;
   for(i=1;i<=DiAvailIndexNum;i++)
   {
        for(k=1;k<=DivIndexObsNum(i);k++)//                           for(j=DiCalBeginYear;j<=DiCalEndYear;j++)
       {
          j= DimIndexIndex(i,k);
          MnTemp1=norm2(Md3IndexPropAtSizeP(i,j)(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i))-Dd3IndexPropAtSizeO(i,k)(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i)));
          if(MnTemp1>0.0)
          {
             MmIndexEffectiveSampleSize(i,j)=Md3IndexPropAtSizeP(i,j)(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i))*(1.0-Md3IndexPropAtSizeP(i,j)(DivIndexSelStartSizeBin(i),DivIndexSelEndSizeBin(i)));
             MmIndexEffectiveSampleSize(i,j)/=MnTemp1;
          }
       }
   }
   report<<"Survey_Index_Obs"<<endl<<DmIndexTotalO<<endl;
   report<<"Survey_Index_Pred"<<endl<<MmIndexP<<endl;
   //ESS
   report<<"Survey_ESSinput"<<endl<<DimIndexESSInput<<endl;
   report<<"Survey_ESSpred"<<endl<<MmIndexEffectiveSampleSize<<endl;
   report<<"DimIndexIndex"<<endl<<DimIndexIndex<<endl;
   //Index Comp
   report<<"Survey_Comp_Obs"<<endl<<Dd3IndexPropAtSizeO<<endl;
   report<<"Survey_Comp_Pred"<<endl<<Md3IndexPropAtSizeP<<endl;
   //CPUE Q
   report<<"Block_CPUE_Catchability_Pred"<<endl<<MvCPUEQ<<endl;
   //INDEX Q
   report<<"Block_Survey_Catchability_Pred"<<endl<<MvIndexQ<<endl;
   //Fleet Selectivity
   report<<"Fleet_selectivity"<<endl<<Md4SelByFleet<<endl;
   //Index Selectivity
   report<<"Index_selectivity"<<endl<<Md3IndexSel<<endl;
   //Ft by Fleet
   report<<"Fishing_Mortality"<<endl<<Md3Ft<<endl;
   //M
   report<<"Natural_Mortality"<<endl<<Md3M<<endl;
   //report<<"#recruitment_and_spawning_stock_relatioship "<<endl;
   if(DiRSFlag!=1)
   {
      if(DiRSFlag==2)
         {
             report<<"B-H_Model_Alpha"<<endl<<mfexp(PbnvLogRSPara(1))<<endl;
             report<<"B-H_Model_Beta"<<endl<<mfexp(PbnvLogRSPara(2))<<endl;
         }
      if(DiRSFlag==3)
         {
             report<<"Ricker_Model_Alpha"<<endl<<mfexp(PbnvLogRSPara(1))<<endl;
             report<<"Ricker_Model_Beta"<<endl<<mfexp(PbnvLogRSPara(2))<<endl;
         }
      if(DiRSFlag==4)
         {
             report<<"CushingE_Model_Alpha"<<endl<<mfexp(PbnvLogRSPara(1))<<endl;
             report<<"CushingE_Model_Beta"<<endl<<mfexp(PbnvLogRSPara(2))<<endl;
             report<<"CushingE_Model_Env"<<endl<<mfexp(PbvEnvCoef)<<endl;
         }
      if(DiRSFlag==5)
         {
             report<<"BH-E_Model_Alpha"<<endl<<mfexp(PbnvLogRSPara(1))<<endl;
             report<<"BH-E_Model_Beta"<<endl<<mfexp(PbnvLogRSPara(2))<<endl;
             report<<"BH-E_Model_Env"<<endl<<mfexp(PbvEnvCoef)<<endl;
         }
      if(DiRSFlag==6)
         {
             report<<"Ricker-E_Model_Alpha"<<endl<<mfexp(PbnvLogRSPara(1))<<endl;
             report<<"Ricker-E_Model_Beta"<<endl<<mfexp(PbnvLogRSPara(2))<<endl;
             report<<"Ricker-E_Model_Env"<<endl<<mfexp(PbvEnvCoef)<<endl;
         }
   }
   else
      report<<"Mean_Recruitment:"<<endl<<mfexp(PbnvLogRSPara(1))<<endl;
   //Recruitment
   report<<"recruitment_Pred"<<endl<<sdr_vRecruitmentP<<endl;
   report<<"recruitment_log_Dev"<<endl<<PbdvLogRecruitDevs<<endl;
   report<<"recruitment_log_DevSD2"<<endl<<PbnvRecruitLogDevsSD2<<endl;   
   report<<"Environment Signal"<<endl<<DvEnvfitZscore<<endl;
   //SSB 
   report<<"Spawning_stock_Biomass"<<endl<<sdr_vSSB2<<endl;
   report<<"Spawning_stock_Biomass_input"<<endl<<sdr_vSSB<<endl;   
   report<<"SexAtSizeLamda"<<endl<<DvLfiftyParas_Ini(4)<<endl;
   //N
   report<<"Abundance_at_Size"<<endl<<Md3NAS<<endl;
   report<<"Abundance_proj"<<endl<<MvNASPJD<<endl;
   //
   report<<"FemaleProp_at_Size_obs"<<endl<<DmSexAtSizeData_Ini<<endl;
   report<<"FemaleProp_at_Size"<<endl<<Md3FemalePropAtSizeP<<endl;   
   report<<"Lfifty"<<endl<<PbLfifty<<endl;
  ////////////////////////////////////////////////////
  //BRP
  //report<<"#Biology_Reference_Point"<<endl;
  report<<"Fmax"<<endl<<MnBPRFMax<<endl;
  report<<"F0.1"<<endl<<MnBPRF01<<endl;
  report<<"F30SPR"<<endl<<MnBPRF30SPR<<endl;
  report<<"F40SPR"<<endl<<MnBPRF40SPR<<endl;
  report<<"FMSY"<<endl<<MnBPRFmsy<<endl;
  report<<"MSY"<<endl<<sdr_nBPRMSY<<endl;
  report<<"Bmsy"<<endl<<MnBPRBiomassMSY<<endl;
  report<<"SSBmsy"<<endl<<MnBPRSSBMSY<<endl;
  report<<"YPR"<<endl<<MnBPRYieldTotal<<endl;
  report<<"SPR"<<endl<<MnBPRSSBTotal<<endl;
  //////////////////////////////////////////////////////////////
  report<<"likelihood_total"<<endl<<ofvTotal<<endl;
  report<<"likelihood_tcatch"<<endl<<MvCatchTotalLikely<<endl;
  report<<"likelihood_pcatch"<<endl<<MvCatchPropLikely<<endl;
  report<<"likelihood_index"<<endl<<MvIndexLikely<<endl;
  report<<"likelihood_pindex"<<endl<<MvIndexPropLikely<<endl;
  report<<"likelihood_Recruit"<<endl<<MnRecruitDevsLikely<<endl;
  report<<"likelihood_pfemal"<<endl<<MnLfiftylikely<<endl;
  ////////////////////////////////////////////////////////////
  //CHECK NUMBER
  //gmax 
  // If this is small (usually 1e-3 or less) then convergence is
  // likely, if large (>1) then convergence is unlikely
  report<<"#objective_function_value::pobjfun->gmax:"<<endl<< objective_function_value::pobjfun->gmax<<endl;
  report<<"#Check Number:"<<endl<<-221<<endl;
  //End!!!!!!!!!!!!!!!!!!!
  // */
  //WriteMCMC
   if(objective_function_value::pobjfun->gmax>1.0  && last_phase() )
   {
       cout<<"gmax is large than  1.0, convergence is unlikely, exit!!"<<endl;
       exit(-99);
   }
}

void model_parameters::set_runtime(void)
{
  dvector temp("{0.01, 0.001, 0.0001}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
  dvector temp1("{5000, 10000, 10000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
}

void model_parameters::final_calcs()
{
  // Inherit From ASAP2 Code (NOAA TOOL BOX 2011) 
  // this code is based on the Widow Rockfish model (from Erik H. Williams, NMFS-Santa Cruz, now Beaufort)
  time(&tEndTime);
  dElapsed_time = difftime(tEndTime,tBeginTime);
  lHour = long(dElapsed_time)/3600;
  lMinute = long(dElapsed_time)%3600/60;
  lSecond = (long(dElapsed_time)%3600)%60;
  //print on screen 
  //cout<<endl<<endl<<"Beginning Time: "<<ctime(&tBeginTime);
  cout<<"Finish Time: "<<ctime(&tEndTime);
  cout<<"Elapsed Time: ";
  cout<<lHour<<" Hours, "<<lMinute<<" minutes, "<<lSecond<<" seconds."<<endl<<endl<<endl;
  //record in log file
  ofRuntimelog<<endl<<endl<<"Beginning Time: "<<ctime(&tBeginTime);
  ofRuntimelog<<endl<<endl<<"Beginning Time: "<<ctime(&tBeginTime);
  ofRuntimelog<<"End Time: "<<ctime(&tEndTime);
  ofRuntimelog<<"The Elapse Time: ";
  ofRuntimelog<<lHour<<" Hours, "<<lMinute<<" minutes, "<<lSecond<<" seconds."<<endl<<endl<<endl;
  ofRuntimelog.close();
  char cCommandLine[200];
  sprintf(cCommandLine,"copy NSLSAP01.rep NSLSAP01_%d_%d_%d.rep",iBeginYear,iEndYear,itimestep);
  system(cCommandLine);
 #if  Debug_Status 
     ofDebug.close();
 #endif
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize=20000000;      //5,000,000 bytes of memory for variable objects.
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000);  
  gradient_structure::set_MAX_NVAR_OFFSET(50000);  
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);
  time(&tBeginTime); //    Get Beginning Time
  //Record the beginning time 
  ofRuntimelog<< endl << "Beginning time : " << ctime(&tBeginTime) << endl;
  //print on screen
  cout<<endl<<" NSLSAP Beginning Run:  "<< ctime(&tBeginTime) << endl;
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
