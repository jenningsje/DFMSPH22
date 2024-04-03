#include "dfmsph_def.h"

//************* changed 2019 04 ********************************************
int main()
{
	code_end=10;
	NUM = FUN_Number();
    FUN_INP();  //Reading and checking the input information
	AP13=pow(AP,0.3333333);AT13=pow(AT,0.3333333);APT13=AP13+AT13;
	Hom= C2b=C3b=BZ=-10;BZ=ZT*ZP/APT13;
	if(code_end<0)return(1);  //code_end<0  in the case of the wrong initialization of any variable or function
	Elab=ECM*APT/AT;
	for(is=1;is<ikGup;is++){rhoMh[is]=-200;rhoPh[is]=-200;}
	if(key_vNN<2)
	{
		if(key_exc==1)FUN_HYEX();	//Initializes the coefficients for the exchange part of the nuclear interaction
		if(key_UC>0) FUN_HCYEX();   //Initializes the coefficients for the exchange part of the Coulomb interaction
	}
	Uex=0.;
	UCex=0.;

	printf("\n\n <<BARR <<<<<<<< RPT0=%6.2f >>>>>>>>>>>\n\n",RPT0);
	
	//****** former BARR(RCCstart,RCCfin,RCCstep);
	//********Searches for the barrier as the maximum of the function the UNCD(R)
	printf("\n <<mai-63 BARR>>>  RCCstart=%7.3f RCCfin=%7.3f RCCstep=%7.3f",RCCstart,RCCfin,RCCstep);
    printf("\n     rat   RCCstart    RCM   iter     UN        UC        Utot");
	for(RCM=RCCstart,ii=0;RCM>RCCfin;RCM-=RCCstep,ii++)
	{
		Entry=-6;FUN_CheckIndex(0,ii,irmax-1);
		if(key_vNN==2)
	    {
			UN=FUN_UNMF(0.,k_up);UC=FUN_UC00(0.,k_up);UNCD=UC+UN;
    	}
		if(key_vNN>2)UNCD=FUN_DFPDEL(RCM);
		
		if(key_vNN<2)	
		{
			for(iter=0;iter<iter_up/1;iter++)  //The iterative procedure for calculating UNCD
			{
				iterdim[ii]=iter;
				if(key_exc==1)
					{UNCD=FUN_DFPFIN(RCM);ratiter=fabs((Umem-UDFP)/(Umem+UDFP-alittle));}
				if(key_exc==0){UNCD=FUN_DFPDEL(RCM);ratiter=alittle;}
				if(ratiter<eps_iter && iter>0)break;
			}
		}
		irup=ii;
		UCdim[ii]=UC;UCDdim[ii]=UCD;
		UNdim[ii]=UN;
		Utotdim[ii]=UNCD;
		RCMdim[ii]=RCM;
			
		if(ii==(ii/50)*50)printf("\n%9.2e %7.3f %7.3f %3ld  %12.4e %10.2e %12.4e ",ratiter,RCCstart,RCM,iter,UN,UC,UNCD);
	}
	
	//******** ** ** Finding the maximum of the function Utot(RCM) *******************
	
	for(RCM=RCCstart,ii=2;RCM>RCCfin;RCM-=RCCstep,ii++)
	{
		if(Utotdim[ii]>Utotdim[ii-1] && Utotdim[ii]>Utotdim[ii+1])
		{BfusDFPsph=Utotdim[ii];r_bar=RCM;i_bar=ii;break;}  
	}
	
	
	if(ZP>0)Hom=FUN_Hom();
	if(irup>irmax)irup=irmax;

	if(chi2WS>0){printf("\n........ Calculating chi2WSmin WS-approximation ..fast ................ \n");
		FUN_WS_DFP();}  //Fits the UN near the barrier by a Woods-Saxon profile
	if(chi2GK>0){printf("\n..Wait calculating chi2GKmin GK-approximation .this lasts about a minute.... \n");
		FUN_GK_DFP();}  //Fits the UN near the barrier by a Gross-Kalinowski profile
	
	printf("\n\n <<mai-39 <<<<<<<< r_bar=%7.3f >>>q_bar=%7.3f>> Ub=%7.2f >>>>>>\n\n",
	r_bar,r_bar/RPT0,BfusDFPsph);

	if((r_bar>(RCCfin-RCCstep) && r_bar<(RCCfin+RCCstep)) || r_bar==RCCstart)
		{printf("\n<<<<< Interval (%4.2lf, %4.2lf) doesn't contain the barrier! >>>>>",RCCstart,RCCfin);code_end=-44;}

	if(code_end<0){printf("\n<<<<<<<<<<<<<<<<<<<<<<<\n<<<<<  ABNORMAL TERMINATION  >>>>> code_end=%3ld\n\n",code_end);
		getchar();}
	else
	{ 
		printf("\n<< NUM=%5d <<<<<<<<<<  SUCCESSFUL TERMINATION  >>>>>>>>>>>>>>>\n",NUM);
		printf("\n<<<<<<<<<<<<<<<<<<<< See output files >>>>>>>>>>>>>>>>>>>>");
		printf("	\n <out_vNN_.c>  <OUT_U_R.c>  <out_dfmsph22.c> <out_inp_dfm22.c> \n");
	}
	FUN_DFPOUT();  //Prints the output information to the file out_dfmsph.out
	getchar();
	return(1);
}


/***************************** WS_DFP *********************************/
/******** approx the DFP potential near the barrier by a WS profile ***********/
void FUN_WS_DFP()
{
	
	chi2WSmin=100;rWSmem=10;VWSmem=-10000;aWSmem=1;
	for(rWS=rWSmin;rWS<rWSmax;rWS+=rWSstep)
 	{
		RP0WS=rWS*AP13;RT0WS=rWS*AT13;
 		for(VWS=VWSmin;VWS<VWSmax;VWS+=VWSstep)
   		{
	  		for(aWS=aWSmin;aWS<aWSmax;aWS+=aWSstep)
			{
		   		i1=0;
				for(ii=0;ii<irup;ii+=deliR)
				{
					if(RCMdim[ii]>fRBfin*r_bar && RCMdim[ii]<fRBstart*r_bar)
					{
				    	UN=VWS/(1.+exp((RCMdim[ii]-rWS*AP13-rWS*AT13)/aWS));
				     	epsVN=(UNdim[ii]-UN)/(UNdim[ii]+UN);
				     	chi2WS+=epsVN*epsVN;i1+=1;
					}
			  	}
		  		chi2WS/=i1;
		   		if(chi2WS<chi2WSmin){chi2WSmin=chi2WS;rWSmem=rWS;VWSmem=VWS;aWSmem=aWS;}
			}
		}
	}
	i1chi2=i1;
	return;
}


/************** 2012 07 21 *** GK_DFP *********************************/
/******** approx the DFP potential by Gross ln *************************/
void FUN_GK_DFP()
{
	double r0appr,Rappr,aappr,chi2appr,epsVNappr,VN,Aappr[3];

	for(ii=0;ii<5;ii++)Amemdim[ii]=90.;chi2GKmin=90;

	for(r0appr=rGKmin;r0appr<rGKmax;r0appr+=rGKstep)
 	{
		Rappr=r0appr*APT13;
		for(Aappr[0]=Amindim[0];Aappr[0]<Amaxdim[0];Aappr[0]+=Astepdim[0])
   		{
			for(aappr=aGKmin;aappr<aGKmax;aappr+=aGKstep)
			{
				for(Aappr[1]=Amindim[1];Aappr[1]<Amaxdim[1];Aappr[1]+=Astepdim[1])
				{
					for(Aappr[2]=Amindim[2];Aappr[2]<Amaxdim[2];Aappr[2]+=Astepdim[2])
					{
						i1=0;
						for(ii=0;ii<irup;ii+=deliR)
						{
							if(RCMdim[ii]>fRBfin*r_bar && RCMdim[ii]<fRBstart*r_bar)
							{
								dR=	RCMdim[ii]-Rappr;
								VN=-Aappr[0]-Aappr[1]*dR-Aappr[2]*pow(dR,2.);
								VN*=log(1.+exp(-1.*dR/aappr));
								epsVNappr=(UNdim[ii]-VN)/(UNdim[ii]+VN);
								chi2GK+=epsVNappr*epsVNappr;i1+=1;
							}
						}
						chi2GK/=i1;
						if(chi2GK<chi2GKmin){chi2GKmin=chi2GK;Amemdim[3]=r0appr;
							for(ii=0;ii<3;ii++)Amemdim[ii]=Aappr[ii];
							Amemdim[4]=aappr;}
					}
				}
			}
		}
	}

	printf("\n i1chi2 r0appr   aapr   A1appr   A2appr    A3appr      chi2GKmin");
	printf("\n %3ld   %6.3f   %6.3f  %6.3f  %6.3f    %6.3f       %9.2e",
		   i1chi2, Amemdim[3],Amemdim[4],Amemdim[0],Amemdim[1],Amemdim[2],chi2GKmin);
	i1chi2=i1; 
	return;
}