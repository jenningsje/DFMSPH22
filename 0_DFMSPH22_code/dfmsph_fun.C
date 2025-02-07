#include "dfmsph_def.h"
/************** Set of functions in use ******************************************/

/*********************************************************************************/
/**********double FUN_CubeInterpD2gDx2(double x, double *ydim, double *gdim)**********/
/*********************************************************************************/
double FUN_CubeInterpD2gDx2(double x, double *ydim, double *gdim)
{
	double D2gDx2;
	D2gDx2=	 gdim[0]*(x-ydim[3])/(ydim[0]-ydim[1])/(ydim[0]-ydim[2])/(ydim[0]-ydim[3])
			+gdim[0]*(x-ydim[2])/(ydim[0]-ydim[1])/(ydim[0]-ydim[2])/(ydim[0]-ydim[3])
			+gdim[0]*(x-ydim[3])/(ydim[0]-ydim[1])/(ydim[0]-ydim[2])/(ydim[0]-ydim[3])
			+gdim[0]*(x-ydim[1])/(ydim[0]-ydim[1])/(ydim[0]-ydim[2])/(ydim[0]-ydim[3])
			+gdim[0]*(x-ydim[2])/(ydim[0]-ydim[1])/(ydim[0]-ydim[2])/(ydim[0]-ydim[3])
			+gdim[0]*(x-ydim[1])/(ydim[0]-ydim[1])/(ydim[0]-ydim[2])/(ydim[0]-ydim[3])
			+gdim[1]*(x-ydim[3])/(ydim[1]-ydim[0])/(ydim[1]-ydim[2])/(ydim[1]-ydim[3])
			+gdim[1]*(x-ydim[2])/(ydim[1]-ydim[0])/(ydim[1]-ydim[2])/(ydim[1]-ydim[3])
			+gdim[1]*(x-ydim[3])/(ydim[1]-ydim[0])/(ydim[1]-ydim[2])/(ydim[1]-ydim[3])
			+gdim[1]*(x-ydim[0])/(ydim[1]-ydim[0])/(ydim[1]-ydim[2])/(ydim[1]-ydim[3])
			+gdim[1]*(x-ydim[2])/(ydim[1]-ydim[0])/(ydim[1]-ydim[2])/(ydim[1]-ydim[3])
			+gdim[1]*(x-ydim[0])/(ydim[1]-ydim[0])/(ydim[1]-ydim[2])/(ydim[1]-ydim[3])
			+gdim[2]*(x-ydim[3])/(ydim[2]-ydim[0])/(ydim[2]-ydim[1])/(ydim[2]-ydim[3])
			+gdim[2]*(x-ydim[1])/(ydim[2]-ydim[0])/(ydim[2]-ydim[1])/(ydim[2]-ydim[3])
			+gdim[2]*(x-ydim[3])/(ydim[2]-ydim[0])/(ydim[2]-ydim[1])/(ydim[2]-ydim[3])
			+gdim[2]*(x-ydim[0])/(ydim[2]-ydim[0])/(ydim[2]-ydim[1])/(ydim[2]-ydim[3])
			+gdim[2]*(x-ydim[1])/(ydim[2]-ydim[0])/(ydim[2]-ydim[1])/(ydim[2]-ydim[3])
			+gdim[2]*(x-ydim[0])/(ydim[2]-ydim[0])/(ydim[2]-ydim[1])/(ydim[2]-ydim[3])
			+gdim[3]*(x-ydim[2])/(ydim[3]-ydim[0])/(ydim[3]-ydim[1])/(ydim[3]-ydim[2])
			+gdim[3]*(x-ydim[1])/(ydim[3]-ydim[0])/(ydim[3]-ydim[1])/(ydim[3]-ydim[2])
			+gdim[3]*(x-ydim[2])/(ydim[3]-ydim[0])/(ydim[3]-ydim[1])/(ydim[3]-ydim[2])
			+gdim[3]*(x-ydim[0])/(ydim[3]-ydim[0])/(ydim[3]-ydim[1])/(ydim[3]-ydim[2])
			+gdim[3]*(x-ydim[1])/(ydim[3]-ydim[0])/(ydim[3]-ydim[1])/(ydim[3]-ydim[2])
			+gdim[3]*(x-ydim[0])/(ydim[3]-ydim[0])/(ydim[3]-ydim[1])/(ydim[3]-ydim[2]);
	
	return(D2gDx2);
}

/************************************************************************************/
/*************double CubeInterpDgDx(double x, double *ydim, double *gdim)************/
/************************************************************************************/
double FUN_CubeInterpDgDx(double x, double *ydim, double *gdim)
{
	double DgDx;
	DgDx=gdim[0]*(x-ydim[2])*(x-ydim[3])/(ydim[0]-ydim[1])/(ydim[0]-ydim[2])/(ydim[0]-ydim[3])
		+gdim[0]*(x-ydim[1])*(x-ydim[3])/(ydim[0]-ydim[1])/(ydim[0]-ydim[2])/(ydim[0]-ydim[3])
	    +gdim[0]*(x-ydim[1])*(x-ydim[2])/(ydim[0]-ydim[1])/(ydim[0]-ydim[2])/(ydim[0]-ydim[3])
    	
		+gdim[1]*(x-ydim[2])*(x-ydim[3])/(ydim[1]-ydim[0])/(ydim[1]-ydim[2])/(ydim[1]-ydim[3])
		+gdim[1]*(x-ydim[0])*(x-ydim[3])/(ydim[1]-ydim[0])/(ydim[1]-ydim[2])/(ydim[1]-ydim[3])
		+gdim[1]*(x-ydim[0])*(x-ydim[2])/(ydim[1]-ydim[0])/(ydim[1]-ydim[2])/(ydim[1]-ydim[3])

	
		+gdim[2]*(x-ydim[1])*(x-ydim[3])/(ydim[2]-ydim[0])/(ydim[2]-ydim[1])/(ydim[2]-ydim[3])
		+gdim[2]*(x-ydim[0])*(x-ydim[3])/(ydim[2]-ydim[0])/(ydim[2]-ydim[1])/(ydim[2]-ydim[3])
		+gdim[2]*(x-ydim[0])*(x-ydim[1])/(ydim[2]-ydim[0])/(ydim[2]-ydim[1])/(ydim[2]-ydim[3])

	
		+gdim[3]*(x-ydim[1])*(x-ydim[2])/(ydim[3]-ydim[0])/(ydim[3]-ydim[1])/(ydim[3]-ydim[2])
		+gdim[3]*(x-ydim[0])*(x-ydim[2])/(ydim[3]-ydim[0])/(ydim[3]-ydim[1])/(ydim[3]-ydim[2])
		+gdim[3]*(x-ydim[0])*(x-ydim[1])/(ydim[3]-ydim[0])/(ydim[3]-ydim[1])/(ydim[3]-ydim[2]);

	return(DgDx);
}

/*********************************************************************************/
/*************double CubeInterp(double x, double *ydim, double *gdim)*************/
/*********************************************************************************/
double FUN_CubeInterp(double x, double *ydim, double *gdim)
{
	double g;
	g=gdim[0]*(x-ydim[1])*(x-ydim[2])*(x-ydim[3])/(ydim[0]-ydim[1])/(ydim[0]-ydim[2])/(ydim[0]-ydim[3])
	+gdim[1]*(x-ydim[0])*(x-ydim[2])*(x-ydim[3])/(ydim[1]-ydim[0])/(ydim[1]-ydim[2])/(ydim[1]-ydim[3])
	+gdim[2]*(x-ydim[0])*(x-ydim[1])*(x-ydim[3])/(ydim[2]-ydim[0])/(ydim[2]-ydim[1])/(ydim[2]-ydim[3])
	+gdim[3]*(x-ydim[0])*(x-ydim[1])*(x-ydim[2])/(ydim[3]-ydim[0])/(ydim[3]-ydim[1])/(ydim[3]-ydim[2]);
	
	return(g);
}

/*********************CheckIndex***************************************************/
/* Checks values of an array index to make impossible overrun the array dimension */
long FUN_CheckIndex(long down, long index, long up)
{
	if(index<down)
  	{
   		printf("\n=====CheckIndex====Entry=%ld=index=%ld=smaller than down=%ld==\n",Entry,index,down);
   		index=down;code_end=Entry;
  	}
	if(index>up)
  	{
   		printf("\n=====CheckIndex====Entry=%ld=index=%ld==larger than up=%ld=\n",Entry,index,up);
   		index=up;code_end=Entry;
  	}
	return(index);
}

/***************spherical Bessel function of the first kind by Rayleigh formula******************/
double FUN_j0BRay(double x)
{
	double y;
	y=1;
	if(fabs(x)>alittle)y=sin(x)/x;return(y);}

/**first derivative of spherical Bessel function of the first kind by Rayleigh formula***********/
double FUN_Dj0BRayDx(double x)
{
	double y;
	y=0;
	if(fabs(x)>alittle)y=cos(x)/x-sin(x)/x/x;
	return(y);}

/**second derivative of spherical Bessel function of the first kind by Rayleigh formula**********/
double FUN_D2j0BRayDx2(double x)
{
	double y,x2;
	y=-1./3.;x2=x*x;
	if(fabs(x)>alittle)y=sin(x)/x*(2./x2-1.)-2.*cos(x)/x2;
	return(y);
}

/**** Coefficients for integration by the Gauss method **************************************/
void FUN_GAUSS()
{
	XX20[ 0]=0.;  //20 knots in the integration by the Gauss method
	XX20[ 1]=0.076526521133497,
	XX20[ 2]=0.227785851141645,
	XX20[ 3]=0.373706088715420,
	XX20[ 4]=0.510867001950827,
	XX20[ 5]=0.636053680726515,
	XX20[ 6]=0.746331906460151,
	XX20[ 7]=0.839116971822219,
	XX20[ 8]=0.912234428251326,
	XX20[ 9]=0.963971927277914,
	XX20[10]=0.993128599185095;
	WW20[ 0]=0.,  //Corresponding weight functions
	WW20[ 1]=0.152753387130726,
	WW20[ 2]=0.149172986472604,
	WW20[ 3]=0.142096109318382,
	WW20[ 4]=0.131688638449177,
	WW20[ 5]=0.118194531961518,
	WW20[ 6]=0.101930119817240,
	WW20[ 7]=0.083276741576705,
	WW20[ 8]=0.062672048334109,
	WW20[ 9]=0.040601429800387,
	WW20[10]=0.017614007139152;

	XX40[ 0]=0.,  //40 knots in the integration by the Gauss method
	XX40[ 1]=0.038772417506051,
	XX40[ 2]=0.116084070675255,
	XX40[ 3]=0.192697580701372,
	XX40[ 4]=0.268152185007254,
	XX40[ 5]=0.341994090825758,
	XX40[ 6]=0.413779204371605,
	XX40[ 7]=0.483075801686179,
	XX40[ 8]=0.549467125095128,
	XX40[ 9]=0.612553889667980,
	XX40[10]=0.671956684614180,
	XX40[11]=0.727318255189927,
	XX40[12]=0.778305651426519,
	XX40[13]=0.824612230833312,
	XX40[14]=0.865959503212260,
	XX40[15]=0.902098806968874,
	XX40[16]=0.932812808278677,
	XX40[17]=0.957916819213792,
	XX40[18]=0.977259949983774,
	XX40[19]=0.990726238699457,
	XX40[20]=0.998237709710559;
	WW40[ 0]=0.,  //Corresponding weight functions
	WW40[ 1]=0.077505947978425,
	WW40[ 2]=0.077039818164248,
	WW40[ 3]=0.076110361900626,
	WW40[ 4]=0.074723169057968,
	WW40[ 5]=0.072886582395804,
	WW40[ 6]=0.070611647391287,
	WW40[ 7]=0.067912045815234,
	WW40[ 8]=0.064804013456601,
	WW40[ 9]=0.061306242492929,
	WW40[10]=0.057439769099392,
	WW40[11]=0.053227846983936,
	WW40[12]=0.048695807635072,
	WW40[13]=0.043870908185673,
	WW40[14]=0.038782167974472,
	WW40[15]=0.033460195282548,
	WW40[16]=0.027937006980023,
	WW40[17]=0.022245849194167,
	WW40[18]=0.016421058381908,
	WW40[19]=0.010498284531153,
	WW40[20]=0.004521277098533;

	XX96[ 0]=0.,  //96 knots in the integration by the Gauss method
	XX96[ 1]=0.016276744849602969579,
	XX96[ 2]=0.048812985136049731112,
	XX96[ 3]=0.081297495464425558994,
	XX96[ 4]=0.113695850110665920911,
	XX96[ 5]=0.145973714654896941989,
	XX96[ 6]=0.178096882367618602759,
	XX96[ 7]=0.210031310460567203603,
	XX96[ 8]=0.241743156163840012328,
	XX96[ 9]=0.273198812591049141487,
	XX96[10]=0.304364944354496353024,
	XX96[11]=0.335208522892625422616,
	XX96[12]=0.365696861472313635031,
	XX96[13]=0.395797649828908603285,
	XX96[14]=0.425478988407300545365,
	XX96[15]=0.454709422167743008636,
	XX96[16]=0.483457973920596359768,
	XX96[17]=0.511694177154667673586,
	XX96[18]=0.539388108324357436227,
	XX96[19]=0.566510418561397168404,
	XX96[20]=0.593032364777572080684,
	XX96[21]=0.618925840125468570386,
	XX96[22]=0.644163403784967106798,
	XX96[23]=0.668718310043916153953,
	XX96[24]=0.692564536642171561344,
	XX96[25]=0.715676812348967626225,
	XX96[26]=0.738030643744400132851,
	XX96[27]=0.759602341176647498703,
	XX96[28]=0.780369043867433217604,
	XX96[29]=0.800308744139140817229,
	XX96[30]=0.819400310737931675539,
	XX96[31]=0.837623511228187121494,
	XX96[32]=0.854959033434601455463,
	XX96[33]=0.871388505909296502874,
	XX96[34]=0.886894517402420416057,
	XX96[35]=0.901460635315852341319,
	XX96[36]=0.915071423120898074206,
	XX96[37]=0.927712456722308690965,
	XX96[38]=0.939370339752755216932,
	XX96[39]=0.950032717784437635756,
	XX96[40]=0.959688291448742539300,
	XX96[41]=0.968326828463264212174,
	XX96[42]=0.975939174585136466453,
	XX96[43]=0.982517263563014677447,
	XX96[44]=0.988054126329623799481,
	XX96[45]=0.992543900323762624572,
	XX96[46]=0.995981842987209290650,
	XX96[47]=0.998364375863181677724,
	XX96[48]=0.999689503883230766828;

	WW96[ 0]=0.,  //Corresponding weight functions
	WW96[ 1]=0.032550614492363166242,
	WW96[ 2]=0.032516118713868835987,
	WW96[ 3]=0.032447163714064269364,
	WW96[ 4]=0.032343822568575928429,
	WW96[ 5]=0.032206204794030250669,
	WW96[ 6]=0.032034456231992663218,
	WW96[ 7]=0.031828758894411006535,
	WW96[ 8]=0.031589330770727168558,
	WW96[ 9]=0.031316425596861355813,
	WW96[10]=0.031010332586313837423,
	WW96[11]=0.030671376123669149014,
	WW96[12]=0.030299915420827593794,
	WW96[13]=0.029896344136328385984,
	WW96[14]=0.029461089958167905970,
	WW96[15]=0.028994614150555236543,
	WW96[16]=0.028497411065085385646,
	WW96[17]=0.027970007616848334440,
	WW96[18]=0.027412962726029242823,
	WW96[19]=0.026826866725591762198,
	WW96[20]=0.026212340735672413913,
	WW96[21]=0.025570036005349361499,
	WW96[22]=0.024900633222483610288,
	WW96[23]=0.024204841792364691282,
	WW96[24]=0.023483399085926219842,
	WW96[25]=0.022737069658329374001,
	WW96[26]=0.021966644438744349195,
	WW96[27]=0.021172939892191298988,
	WW96[28]=0.020356797154333324595,
	WW96[29]=0.019519081140145022410,
	WW96[30]=0.018660679627411467385,
	WW96[31]=0.017782502316045260838,
	WW96[32]=0.016885479864245172450,
	WW96[33]=0.015970562902562291381,
	WW96[34]=0.015038721026994938006,
	WW96[35]=0.014090941772314860916,
	WW96[36]=0.013128229566961572637,
	WW96[37]=0.012151604671088319635,
	WW96[38]=0.011162102099838498591,
	WW96[39]=0.010160770535008415758,
	WW96[40]=0.009148671230783386633,
	WW96[41]=0.008126876925698759217,
	WW96[42]=0.007096470791153865269,
	WW96[43]=0.006058545504235961683,
	WW96[44]=0.005014202742927517693,
	WW96[45]=0.003964554338444686674,
	WW96[46]=0.002910731817934946408,
	WW96[47]=0.001853960788946921732,
	WW96[48]=0.000796792065552012429;
	return;
}

/****************************Number**********************************/
/*************** Calculates the number of the run *******************/
 int FUN_Number()
{
	FILE *numbeR;
	if( (numbeR=fopen("number.c","r")) == NULL ) NUM=1;
	else
	{
  		fscanf(numbeR,"%d",&NUM);
  		fclose(numbeR);
	}
	numbeR=fopen("number.c","w");
	fprintf(numbeR,"%6d",NUM+1);
	fclose(numbeR);
	return(NUM);
}