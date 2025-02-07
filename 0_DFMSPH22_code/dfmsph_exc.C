#include "dfmsph_def.h"
/****************** Exchange parts of the DFM potential *****************************/

/****************************** KFERMI **********************************************/
/********************** The Fermi momentum ******************************************/
double FUN_KFERMI(double r, long iNUC)
{
	static double kF;
	rho=FUN_RHOtable(r,iNUC,0);
	kF=pow(1.5*pi*pi*rho,0.6666666);
	absGRADrho=sqrt(DrhoDr*DrhoDr);
	LAPLACErho=D2rhoDr2+2.*DrhoDr/r;
	kF+=5.*CS*absGRADrho*absGRADrho/rho/rho/3.;
	if(r>0.2*(RP0*(-iNUC+1)+RT0*iNUC)) kF+=5.*LAPLACErho/rho/36.;
	if(kF<0.){printf("\n<<<<<<<< U_EXC: wrong kF=%9.2e >>>>>>>>>\n\n",kF);getchar();code_end=-51;return(1.);}
	kF=sqrt(kF);
	return(kF);
}

/************************** HEX*****************************************************/
double FUN_HEX(double t, double s, long iNUC)
{
	static long ih;
	static double jB0_t_rr,z3,kF,integ,rr,h,z;
	h=0;integ=0;
	for(ih=1;ih<ikGup;ih++)
	{
		rr=kMdim[ih]*rup/k_up;
		kF=FUN_KFERMI(rr,iNUC);
		z=kF*s;z3=z*z*z+alittle;
		jB0_t_rr=1.;if(fabs(t*rr)>alittle)jB0_t_rr=sin(t*rr)/(t*rr);
		IntMh[ih]=rho*3.*(sin(z)-z*cos(z))/z3*jB0_t_rr*rr*rr;
		rhoMh[ih]=rho;
	    rr=kPdim[ih]*rup/k_up;
		kF=FUN_KFERMI(rr,iNUC);
		z=kF*s;z3=z*z*z+alittle;
		jB0_t_rr=1.;if(fabs(t*rr)>alittle)jB0_t_rr=sin(t*rr)/(t*rr);
		IntPh[ih]=rho*3.*(sin(z)-z*cos(z))/z3*jB0_t_rr*rr*rr;
		rhoPh[ih]=rho;
		integ+=WW[ih]*(IntPh[ih]+IntMh[ih]);
   	}
	integ*=(rup-0.)/2.;
	h=integ*4.*pi;
	return(h);
}

/*********************************YEX****************************************************/
double FUN_YEX(double t, double s, long iNUC)
{
	static long iy;
	static double y,IntP,IntM,integ;
	y=0;integ=0;
	for(iy=1;iy<ikGup;iy++)
	{
    	if(rhoMh[iy]<0. || rhoPh[iy]<0)
		{
			printf("\n<<<<<<<<<<<<<<<<<<<too early to call YEX, \n<<<<<<<<<<<<<<<IntMh and IntPh are not yet defined command?\n");
			code_end=-123;getchar();return(-10.);
		}
		IntM=IntMh[iy]*exp(-beDD*rhoMh[iy]);
		IntP=IntPh[iy]*exp(-beDD*rhoPh[iy]);
		integ+=WW[iy]*(IntP+IntM);
   	}
	integ*=(rup-0.)/2.;
	y=integ*4.*pi;
	return(y);
}

/*********************************ZEX****************************************************/
double FUN_ZEX(double t, double s, long iNUC, long m)
{
	static long iz;
	static double z,IntP,IntM,integ;
	z=0;integ=0;
	for(iz=1;iz<ikGup;iz++)
	{
    	if(rhoMh[iz]<0. || rhoPh[iz]<0)
		{
			printf("\n<<<<<<<<<<<<<<<<<<<too early to call ZEX, \n<<<<<<<<<<<<<<<IntMh and IntPh are not yet defined command?\n");
			code_end=-124;getchar();return(-10.);
		}
		IntM=IntMh[iz]*rhoMh[iz];
		IntP=IntPh[iz]*rhoPh[iz];
		if(m==2){
			IntM*=rhoMh[iz];
			IntP*=rhoPh[iz];}
		integ+=WW[iz]*(IntP+IntM);
   }
	integ*=(rup-0.)/2.;
	z=integ*4.*pi;
	return(z);
}

/*** Initializing of the coefficients HEX, YEX and ZEX for the exchange part of the nuclear interaction ***/
void FUN_HYEX()
{
	long it,is;
	for(iNUC=0;iNUC<2;iNUC++)
 	{
  		tup=k_up;
		rup=sup=Crup*RT0;
  		  	
		for(it=1;it<ikGup;it++)
  		{
   			tttM[it]=kMdim[it]*tup/k_up;
   			for(is=1;is<ikGup;is++)
			{
				sssM[is]=kMdim[is]*sup/k_up;
	  		 	hdimMM[it][is][iNUC]=FUN_HEX(tttM[it],sssM[is],iNUC);
	  		 	if(fabs(alDD)>alittle)ydimMM[it][is][iNUC]=FUN_YEX(tttM[it],sssM[is],iNUC);
	  			if(fabs( gDD)>alittle)
				{
					zdimMM[it][is][iNUC][0]=FUN_ZEX(tttM[it],sssM[is],iNUC,1);
	  				zdimMM[it][is][iNUC][1]=FUN_ZEX(tttM[it],sssM[is],iNUC,2);
				}
	  	 		sssP[is]=kPdim[is]*sup/k_up;
	  			hdimMP[it][is][iNUC]=FUN_HEX(tttM[it],sssP[is],iNUC);
	  	 		if(fabs(alDD)>alittle)ydimMP[it][is][iNUC]=FUN_YEX(tttM[it],sssP[is],iNUC);
	  			if(fabs( gDD)>alittle)
				{
					zdimMP[it][is][iNUC][0]=FUN_ZEX(tttM[it],sssP[is],iNUC,1);
	  				zdimMP[it][is][iNUC][1]=FUN_ZEX(tttM[it],sssP[is],iNUC,2);
				}
			}
			tttP[it]=kPdim[it]*tup/k_up;
			for(is=1;is<ikGup;is++)
			{
	 			sssM[is]=kMdim[is]*sup/k_up;
	   			hdimPM[it][is][iNUC]=FUN_HEX(tttP[it],sssM[is],iNUC);
	   			if(fabs(alDD)>alittle)ydimPM[it][is][iNUC]=FUN_YEX(tttP[it],sssM[is],iNUC);
				if(fabs(gDD)>alittle)
				{
					zdimPM[it][is][iNUC][0]=FUN_ZEX(tttP[it],sssM[is],iNUC,1);
					zdimPM[it][is][iNUC][1]=FUN_ZEX(tttP[it],sssM[is],iNUC,2);
				}
   				sssP[is]=kPdim[is]*sup/k_up;
				hdimPP[it][is][iNUC]=FUN_HEX(tttP[it],sssP[is],iNUC);
				if(fabs(alDD)>alittle)ydimPP[it][is][iNUC]=FUN_YEX(tttP[it],sssP[is],iNUC);
				if(fabs(gDD)>alittle)
   				{
	  				zdimPP[it][is][iNUC][0]=FUN_ZEX(tttP[it],sssP[is],iNUC,1);
	   				zdimPP[it][is][iNUC][1]=FUN_ZEX(tttP[it],sssP[is],iNUC,2);
				}
    		}
		}
 	}
	return;
}

/*** Initializing of the coefficients HEX for the exchange part of the Coulomb interaction ***/
void FUN_HCYEX()
{
	long it,is;
	for(iNUC=0;iNUC<2;iNUC++)
 	{
  		tup=k_up;
		rup=sup=Crup*RT0;

		//diffusP[0]=diffusCP; diffusT[0]=diffusCT;
		typ=0;
		for(it=1;it<ikGup;it++)
  		{
   			tttM[it]=kMdim[it]*tup/k_up;
   			for(is=1;is<ikGup;is++)
			{
				sssM[is]=kMdim[is]*sup/k_up;
	  		 	hCdimMM[it][is][iNUC]=FUN_HEX(tttM[it],sssM[is],iNUC);
	  	 		sssP[is]=kPdim[is]*sup/k_up;
	  			hCdimMP[it][is][iNUC]=FUN_HEX(tttM[it],sssP[is],iNUC);
			}
			tttP[it]=kPdim[it]*tup/k_up;
			for(is=1;is<ikGup;is++)
			{
	 			sssM[is]=kMdim[is]*sup/k_up;
	   			hCdimPM[it][is][iNUC]=FUN_HEX(tttP[it],sssM[is],iNUC);
   				sssP[is]=kPdim[is]*sup/k_up;
				hCdimPP[it][is][iNUC]=FUN_HEX(tttP[it],sssP[is],iNUC);
    		}
		}
		typ=1;
 	}
	return;
}

/**************** The exchange part of the NN interactions *****************/
double FUN_VEX(double s)
{
	static double vex;
	vex=0;if(fabs(s)<alittle)return(alittle);
	else vex=Aex1*exp(-aex1*s)/aex1/s+Aex2*exp(-aex2*s)/aex2/s+Aex3*exp(-aex3*s)/aex3/s;
	return(vex);
}

/************** The exchange part of the Coulomb interactions ***************/
/************ between two points of the nuclei separated by s ***************/
double FUN_VCEX(double s)
{
	static double vCex;
	vCex=0.; if(fabs(s)<alittle)return(alittle);
	else vCex=1./s;
	return(vCex);
}

/************* GEX ******** The multiplier G ********************************/
double FUN_GEX(double R, double s)
{
	long itG;
	static double jB0_R_tt,integ,IntP,IntM,tt,G;
	integ=0;
	if(priz!=-1 && priz!=1){printf("\n<<<<<<<<<<<<bad priz=%3ld>>>>>>>>>",priz);code_end=-124;return(1e4);}
	Entry=-111;FUN_CheckIndex(1,is,ikGup-1);
	for(itG=1;itG<ikGup;itG++)
	{
    	Entry=-222;FUN_CheckIndex(1,itG,ikGup-1);
	  	tt=tttM[itG];
		if(priz==-1)for(iNUC=0;iNUC<2;iNUC++)
		{
			hdim[itG][is][iNUC]=hdimMM[itG][is][iNUC];
		 	ydim[itG][is][iNUC]=ydimMM[itG][is][iNUC];
			zdim[itG][is][iNUC][0]=zdimMM[itG][is][iNUC][0];
			zdim[itG][is][iNUC][1]=zdimMM[itG][is][iNUC][1];
		}
		if(priz==1)for(iNUC=0;iNUC<2;iNUC++)
		{
			hdim[itG][is][iNUC]=hdimMP[itG][is][iNUC];
		 	ydim[itG][is][iNUC]=ydimMP[itG][is][iNUC];
			zdim[itG][is][iNUC][0]=zdimMP[itG][is][iNUC][0];
			zdim[itG][is][iNUC][1]=zdimMP[itG][is][iNUC][1];
		}
		jB0_R_tt=1.;if(fabs(tt*R)>alittle)jB0_R_tt=sin(tt*R)/(tt*R);
    	IntM=(hdim[itG][is][0]*hdim[itG][is][1]
			+alDD*ydim[itG][is][0]*ydim[itG][is][1]
			-gDD *(zdim[itG][is][0][1]*zdim[itG][is][1][0]+zdim[itG][is][1][1]*zdim[itG][is][0][0])
			)*jB0_R_tt*tt*tt;
	  	tt=tttP[itG];
		if(priz==-1)for(iNUC=0;iNUC<2;iNUC++)
		{
			hdim[itG][is][iNUC]=hdimPM[itG][is][iNUC];
		 	ydim[itG][is][iNUC]=ydimPM[itG][is][iNUC];
			zdim[itG][is][iNUC][0]=zdimPM[itG][is][iNUC][0];
			zdim[itG][is][iNUC][1]=zdimPM[itG][is][iNUC][1];
		}
		if(priz==1)for(iNUC=0;iNUC<2;iNUC++)
		{
			hdim[itG][is][iNUC]=hdimPP[itG][is][iNUC];
		 	ydim[itG][is][iNUC]=ydimPP[itG][is][iNUC];
			zdim[itG][is][iNUC][0]=zdimPP[itG][is][iNUC][0];
			zdim[itG][is][iNUC][1]=zdimPP[itG][is][iNUC][1];
		}
		jB0_R_tt=1.;if(fabs(tt*R)>alittle)jB0_R_tt=sin(tt*R)/(tt*R);
		IntP=(hdim[itG][is][0]*hdim[itG][is][1]
			+alDD*ydim[itG][is][0]*ydim[itG][is][1]
			-gDD *(zdim[itG][is][0][1]*zdim[itG][is][1][0]+zdim[itG][is][1][1]*zdim[itG][is][0][0])
			)*jB0_R_tt*tt*tt;
		integ+=WW[itG]*(IntP+IntM);
   	}
	integ*=(tup-0.)/2.;
	G=CDD/2./pi/pi*integ;
	return(G);
}

/************ GCEX ******** The multiplier GC *******************************/
double FUN_GCEX(double R, double s)
{
	long itG;
	static double jB0_R_tt,integ,IntP,IntM,tt,G;
	integ=0;
	if(priz!=-1 && priz!=1){printf("\n<<<<<<<<<<<<bad priz=%3ld>>>>>>>>>>>>>",priz);code_end=-134;return(1e4);}
	Entry=-111;FUN_CheckIndex(1,is,ikGup-1);
	for(itG=1;itG<ikGup;itG++)
	{
    	Entry=-223;FUN_CheckIndex(1,itG,ikGup-1);
	  	tt=tttM[itG];
		if(priz==-1)for(iNUC=0;iNUC<2;iNUC++)
			hCdim[itG][is][iNUC]=hCdimMM[itG][is][iNUC];
		if(priz==1)for(iNUC=0;iNUC<2;iNUC++)
			hCdim[itG][is][iNUC]=hCdimMP[itG][is][iNUC];
		jB0_R_tt=1.;if(fabs(tt*R)>alittle)jB0_R_tt=sin(tt*R)/(tt*R);
    	IntM=hCdim[itG][is][0]*hCdim[itG][is][1]*jB0_R_tt*tt*tt;
	  	tt=tttP[itG];
		if(priz==-1)for(iNUC=0;iNUC<2;iNUC++)
			hCdim[itG][is][iNUC]=hCdimPM[itG][is][iNUC];
		if(priz==1)for(iNUC=0;iNUC<2;iNUC++)
			hCdim[itG][is][iNUC]=hCdimPP[itG][is][iNUC];
 		jB0_R_tt=1.;if(fabs(tt*R)>alittle)jB0_R_tt=sin(tt*R)/(tt*R);
		IntP=hCdim[itG][is][0]*hCdim[itG][is][1]*jB0_R_tt*tt*tt;
		integ+=WW[itG]*(IntP+IntM);
   	}
	integ*=(tup-0.)/2.;
	G=1./2./pi/pi*integ;
	return(G);
}

/************* The exchange part of the nuclear interaction *****************/
double FUN_UEX(double El, double R)
{
	static double argB,jB0_k_R_ss,k_R,integ,IntP,IntM,ss,vex;
	k_R=sqrt(2./m_red)*m_nucleon/hbar*sqrt(fabs(ECM-UDFP));
	integ=0;
	for(is=1;is<ikGup;is++)
	{
		ss=sssM[is];priz=-1;
		jB0_k_R_ss=1.;argB=k_R*ss;
		if(fabs(argB)>alittle)
		{
			if(ECM>UDFP)jB0_k_R_ss=sin(argB)/(argB);
   			else
   			{
   				if(argB<500.) jB0_k_R_ss=(exp(argB)-exp(-argB))/2./argB;
   				else
   				{
   					code_end=-125;printf("\n<<<<< code_end=%3ld <<<< bad argument for Bessel function = %9.2e \n",
   						                                code_end,argB);
   					printf("\nUex=%9.2e ECM=%9.2e UDFP=%9.2e RCM=%9.2e m_red=%9.2e\n vex=%9.2e k_R=%9.2e ss=%9.2e\n>>>command?\n",
   				               Uex,ECM,UDFP,RCM,m_red,vex,k_R,ss);getchar();return(-10.);
   				}
   			}
   		}
  	  	GexMdim[is]=FUN_GEX(R,ss);vex=FUN_VEX(ss);
  		IntM=GexMdim[is]*jB0_k_R_ss*vex*ss*ss;
  		ss=sssP[is];priz=1;
  		jB0_k_R_ss=1.;argB=k_R*ss;
  		if(fabs(argB)>alittle)
		{
			if(ECM>UDFP)jB0_k_R_ss=sin(argB)/(argB);
  			else jB0_k_R_ss=(exp(argB)-exp(-argB))/2./argB;
		}
  	  	GexPdim[is]=FUN_GEX(R,ss);
  		IntP=GexPdim[is]*jB0_k_R_ss*FUN_VEX(ss)*ss*ss;
  		integ+=WW[is]*(IntP+IntM);
	}
	integ*=(sup-0.)/2.;

	integ=4.*pi*(1.-CEl*El/(AP+alittle))*integ;
	return(integ);
}

/************** The exchange part of the Coulomb interactions ***************/
double FUN_UCEX(double El, double R)
{
	static double argB,jB0_k_R_ss,k_R,integ,IntP,IntM,ss,vCex;
	k_R=sqrt(2./m_red)*m_nucleon/hbar*sqrt(fabs(ECM-UDFP));
	integ=0;
	for(is=1;is<ikGup;is++)
	{
		ss=sssM[is];priz=-1;
		jB0_k_R_ss=1.;argB=k_R*ss;
		if(fabs(argB)>alittle)
		{
			if(ECM>UDFP)jB0_k_R_ss=sin(argB)/(argB);
   			else
   			{
   				if(argB<500.) jB0_k_R_ss=(exp(argB)-exp(-argB))/2./argB;
   				else
   				{
   					code_end=-126;printf("\n<<<<< code_end=%3ld <<<< bad argument for Bessel function = %9.2e \n",
   						                                code_end,argB);
   					printf("\n  UCex=%9.2e ECM=%9.2e UDFP=%9.2e RCM=%9.2e m_red=%9.2e   >>>command?\n",
   				               UCex,ECM,UDFP,RCM,m_red);getchar();return(-10.);
   				}
   			}
   		}
  	  	GCexMdim[is]=FUN_GCEX(R,ss);vCex=FUN_VCEX(ss);
  		IntM=GCexMdim[is]*jB0_k_R_ss*vCex*ss*ss;
  		ss=sssP[is];priz=1;
  		jB0_k_R_ss=1.;argB=k_R*ss;
  		if(fabs(argB)>alittle)
		{
			if(ECM>UDFP)jB0_k_R_ss=sin(argB)/(argB);
  			else jB0_k_R_ss=(exp(argB)-exp(-argB))/2./argB;
		}
  	  	GCexPdim[is]=FUN_GCEX(R,ss);
  		IntP=GCexPdim[is]*jB0_k_R_ss*FUN_VCEX(ss)*ss*ss;
  		integ+=WW[is]*(IntP+IntM);
	}
	integ*=(sup-0.)/2.;

	integ=8.*pi*a_coul*(2./pi)*ZP*ZT/AP/AT*integ;
	return(integ);
}