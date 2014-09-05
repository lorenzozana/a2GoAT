#include "PPi0Example.h"

PPi0Example::PPi0Example()
{ 
    GHistBGSub::InitCuts(-20, 15, -100, -40);
    GHistBGSub::AddRandCut(35, 95);
    
    SetTarget(115417.416); // Mass of the Target Tin mass AT REST: Sn-120: M=111688.180MeV, Sn124: M=115417.416MeV, Pb208: M=193728.955MeV, Ni58: M=53966.385
        
    time 	= new GH1("time", 	"time", 	1400, -700, 700);
    time_cut 	= new GH1("time_cut", 	"time_cut", 	1400, -700, 700);

    time_2g 	= new GH1("time_2g",	"time_2g", 	1400, -700, 700);
    time_2g_cut = new GH1("time_2g_cut","time_2g_cut", 	1400, -700, 700);

    IM 		= new GH1("IM", 	"IM", 		400,   0, 400);
    IM_2g 	= new GH1("IM_2g", 	"IM_2g", 	400,   0, 400);
  
    MMass		= new GH1("MMass", 	"MMass", 	 	1200,   100000., 120000.);     
    MMass_2g	= new GH1("MMass_2g", 	"MMass_2g", 	1200, 100000., 120000.); 

    MMom		= new GH1("MMom", 	"MMom;q(fm^{-1})", 	 	400, 0., 2.0);     
    MMom_2g	= new GH1("MMom_2g", 	"MMom_2g;q(fm^{-1})", 	400,   0, 2.0); 

    MMomDan		= new GH1("MMomDan", 	"MMom;q(fm^{-1})", 	 	400, 0., 2.0);     
    MMomDan_2g	= new GH1("MMomDan_2g", 	"MMom_2g;q(fm^{-1})", 	400,   0, 2.0); 

    DeltaE_CM_Dan_2g	= new GH1("DeltaE_CM_Dan_2g", 	"DeltaE_{CM} 2 #gamma ; MeV", 	400, -60., 60.); 

    int bin_q=300;
    int bin_e=23;
    char Title[60];
    char Title2[160];
    for (int i=0; i<bin_q; i++) {
      for (int j=0; j<bin_e; j++) {
	sprintf(Title,"pimissen_q_%d_%.3f_%.3f",j,i*0.005,(i*0.005+0.005));
	sprintf(Title2,"DeltaE_{CM} 2 #gamma for q_{bin}=%i and E_{bin}=%i; DeltaE_{CM}(MeV)",i,j);
	DeltaE_Missmom_BeamE[i*bin_e+j] =  new GH1(Title,Title2,100,-60,60);
      }
    }

}

PPi0Example::~PPi0Example()
{
}

Bool_t	PPi0Example::Start()
{
    if(!IsGoATFile())
    {
        cout << "ERROR: Input File is not a GoAT file." << endl;
        return kFALSE;
    }
    SetAsPhysicsFile();

    TraverseValidEvents();

	return kTRUE;
}

void	PPi0Example::ProcessEvent()
{
	// fill time diff (tagger - pi0), all pi0
	FillTime(*pi0,time);
	FillTimeCut(*pi0,time_cut);
	
	// fill missing mass, all pi0
	FillMissingMass(*pi0,MMass);

	// fill missing momentum, all pi0
	FillMissingMomentum(*pi0,MMom);

	// fill missing momentum calculated using Dan routine, all pi0
	FillMissingMomentumDan(*pi0,MMomDan);
	
	
	// fill invariant mass, all pi0
	FillMass(*pi0,IM);
		
	// Some neutral decays
    for (Int_t i = 0; i < pi0->GetNParticles(); i++)
	{		
        // Fill MMass for 2 photon decay
        if ((pi0->GetNSubParticles(i) == 2) & (pi0->GetNSubPhotons(i) == 2))
        {
		// fill time diff (tagger - pi0), this pi0
		FillTime(*pi0,i,time_2g);
		FillTimeCut(*pi0,i,time_2g_cut);
			
		// fill missing mass, this pi0
            	FillMissingMass(*pi0,i,MMass_2g);

		// fill missing momentum, this pi0
            	FillMissingMomentum(*pi0,i,MMom_2g);

		// fill missing momentum calculated using Dan routine, this pi0
            	FillMissingMomentumDan(*pi0,i,MMomDan_2g);

		FillDeltaE_Missmom(*pi0,DeltaE_Missmom_BeamE);
            
		// fill invariant mass, this pi0
        	FillMass(*pi0,i,IM_2g);

		FillDeltaE(*pi0,i,DeltaE_CM_Dan_2g);
        }

	}
	
}

void	PPi0Example::ProcessScalerRead()
{

    //time.ScalerReadCorrection(5);
}
