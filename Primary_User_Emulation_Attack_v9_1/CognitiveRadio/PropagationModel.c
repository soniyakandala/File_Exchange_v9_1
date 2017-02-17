/************************************************************************************
 * Copyright (C) 2016                                                               *
 * TETCOS, Bangalore. India                                                         *
 *                                                                                  *
 * Tetcos owns the intellectual property rights in the Product and its content.     *
 * The copying, redistribution, reselling or publication of any or all of the       *
 * Product or its content without express prior written consent of Tetcos is        *
 * prohibited. Ownership and / or any other right relating to the software and all  *
 * intellectual property rights therein shall remain at all times with Tetcos.      *
 *                                                                                  *
 * Author:    Shashi Kant Suman                                                     *
 *                                                                                  *
 * ---------------------------------------------------------------------------------*/
#include "main.h"

int fn_NetSim_Propagation_CalculatePathLoss(double dTxPower,
	double dDistance,
	double dFrequency,
	double dPathLossExponent,
	double* dReceivedPower)
{
	double fpi=3.1417f;		
	double nGainTxr=0;		
	double nGainRver=0;	
	int nDefaultDistance=1; 
	double fA1,fWaveLength=0.0; 
	double fA1dB, fA2dB;
	double fTXPowerdBm;

	fTXPowerdBm=10 * log10(dTxPower);
	if(dPathLossExponent==0.0 || dDistance==0.0)
	{
		*dReceivedPower = fTXPowerdBm;
		return 0;
	}
	// get the gain of the Transmitter
	nGainTxr=0;
	// get the gain of the receiver
	nGainRver=0;
	// Calculate Lambda
	fWaveLength=(double)(300.0/(dFrequency * 1.0));//pathloss
	// Calculate  (4*3.14*do)
	fA1=fWaveLength/(4*(double) fpi * nDefaultDistance );
	//Calculate  20log10[ Lambda/(4*3.14*do) ]
	fA1dB =  20 * log10(fA1);
	//Calculate  10 * n *log10 (do/d)
	fA2dB =  10 * dPathLossExponent * log10(nDefaultDistance /dDistance);
	//Calculate  [Pt]  +  [Gt]  +  [Gr]  +  20log10[ Lemda/(4*3.14*do) ] + ( 10 * n *log10 (do/d) )
	*dReceivedPower = fTXPowerdBm + nGainTxr + nGainRver + fA1dB + fA2dB;
	return 0;
}

int fn_NetSim_Propagation_CalculateShadowLoss(unsigned long* ulSeed1,
	unsigned long* ulSeed2,
	double* dReceivedPower,
	double dStandardDeviation)
{
	double ldGaussianDeviate=0.0;
	double ldRandomNumber = 0.0;
	static int nIset = 0;
	static double fGset=0;
	double fFac,fRsq,fV1,fV2;	
	double d_sign1;
	if(dStandardDeviation == 0.0)
		return 0;
	if(nIset==0)
	{
		do
		{		
			// call the random number generator function
			ldRandomNumber = fn_NetSim_Utilities_GenerateRandomNo(ulSeed1,ulSeed2)/NETSIM_RAND_MAX;		
			fV1=(double)(2.0*ldRandomNumber-1.0);

			ldRandomNumber = fn_NetSim_Utilities_GenerateRandomNo(ulSeed1,ulSeed2)/NETSIM_RAND_MAX;			
			fV2=(double)(2.0*ldRandomNumber-1.0);

			fRsq=fV1*fV1+fV2*fV2;
		}while(fRsq>=1.0 || fRsq==0.0);
		fFac=(double)(sqrt(-2.0*log(fRsq)/fRsq));
		fGset=fV1*fFac;	
		nIset=1;
		ldGaussianDeviate = fV2*fFac;
	}
	else
	{
		nIset=0;
		ldGaussianDeviate = fGset;
	}
	d_sign1 = fn_NetSim_Utilities_GenerateRandomNo(ulSeed1,ulSeed2)/NETSIM_RAND_MAX;
	
	// This is done to ensure there is constructive and destructive shadowing
	if(d_sign1 <= 0.5) 
	{
		// Assign the Received power due to shadowloss.	
		*dReceivedPower -= ldGaussianDeviate * sqrt(dStandardDeviation); 
	}
	else
	{	
		// Assign the Received power due to shadowloss.	
		*dReceivedPower += ldGaussianDeviate * sqrt(dStandardDeviation);
	}

	return 0;
}

/*********************************************************************************************
                    This Function is used to calculate the Fading Loss(Rayleigh Fading)
**********************************************************************************************/
int fn_NetSim_Propagation_CalculateFadingLoss(unsigned long *ulSeed1, unsigned long *ulSeed2,double *dReceivedPower,int fm1)
{
	double sign;
	double dFadingPower;
	
	if(fm1==1)// if rayleigh fading is set
	{
		dFadingPower = fn_NetSim_Utilities_GenerateRandomNo(ulSeed1,ulSeed2)/NETSIM_RAND_MAX;
		if(dFadingPower<0 ||dFadingPower>1) // Sanity check -- this condition should never occur
		{
			dFadingPower=0;
			fnNetSimError("Random Number generated while calculating fading loss is not in the range [0,1]");
			return 0;
		}
		else
		{
			dFadingPower = log(dFadingPower);
			dFadingPower = 10*log10(-dFadingPower); //Rayleigh Fading is an exponential distribution with mean 1. So, we take an inverse exponential funtion to calculate Fading loss.
			sign = fn_NetSim_Utilities_GenerateRandomNo(ulSeed1,ulSeed2)/NETSIM_RAND_MAX;
			if(sign<0.5)
			{
				dFadingPower = -1*dFadingPower; //this is done to ensure constructive/destructive fading
			}

		}
		
	}
	else
	{
		dFadingPower = 0;
	}
	*dReceivedPower += dFadingPower;
	return 0;
}

#include "802_22.h"

int fn_NetSim_CalculateReceivedPower(NETSIM_ID TX,NETSIM_ID RX)
{
	double txpower;
	double power;
	double dDistance = fn_NetSim_Utilities_CalculateDistance(DEVICE_POSITION(TX),DEVICE_POSITION(RX));
	NetSim_LINKS* link=DEVICE_PHYLAYER(TX,1)->pstruNetSimLinks;
	double frequency;

	if(TX==RX)
	{
		dReceivedPower[TX-1][RX-1]=0;
		return 1;
	}
	if(DEVICE_TYPE(TX)==BASESTATION)
	{
		txpower=((BS_PHY*)DEVICE_PHYVAR(TX,1))->dTXPower;
		frequency = ((BS_PHY*)DEVICE_PHYVAR(TX,1))->pstruOpratingChannel->dLowerFrequency;
	}
	else
	{
		txpower=((CPE_PHY*)DEVICE_PHYVAR(TX,1))->dTXPower;
		frequency = ((CPE_PHY*)DEVICE_PHYVAR(TX,1))->pstruOperatingChannel->dLowerFrequency;
	}
	//calculate the path loss
	fn_NetSim_Propagation_CalculatePathLoss(txpower,
		dDistance,
		frequency,
		link->puniMedProp.pstruWirelessLink.dPathLossExponent,
		&power);

	//Calculate the shadow loss
	fn_NetSim_Propagation_CalculateShadowLoss(
		&NETWORK->ppstruDeviceList[TX-1]->ulSeed[0],
		&NETWORK->ppstruDeviceList[TX-1]->ulSeed[1],
		&power,
		link->puniMedProp.pstruWirelessLink.dStandardDeviation);

	fprintf(stdout,"\n%d\t%d\t%lf\n",TX,RX,power);
	dReceivedPower[TX-1][RX-1]=power;
	return 1;
}

int fn_NetSim_CR_CalulateReceivedPower()
{
	NETSIM_ID i,j;
	printf("CR Protocol, Calculating Received power...\n");
	fprintf(stdout,"TX\tRX\tPower\n");
	dReceivedPower=(double**)calloc(NETWORK->nDeviceCount,sizeof* dReceivedPower);
	for(i=0;i<NETWORK->nDeviceCount;i++)
	{
		if(DEVICE_MACLAYER(i+1,1)->nMacProtocolId == MAC_PROTOCOL_IEEE802_22)
		{
			dReceivedPower[i]=(double*)calloc(NETWORK->nDeviceCount,sizeof* dReceivedPower[i]);
			for(j=0;j<NETWORK->nDeviceCount;j++)
			{
				if(DEVICE_MACLAYER(j+1,1)->nMacProtocolId == MAC_PROTOCOL_IEEE802_22)
				{
					fn_NetSim_CalculateReceivedPower(i+1,j+1);
				}
			}
		}
	}
	return 0;
}