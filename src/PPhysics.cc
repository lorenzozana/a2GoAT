#ifndef __CINT__

#include "PPhysics.h"

PPhysics::PPhysics() 
{ 
}

PPhysics::~PPhysics()
{
}

Bool_t	PPhysics::Init(const Char_t *configfile)
{
	
	return kTRUE;
}

void	PPhysics::Reconstruct()
{
}

void PPhysics::FillMissingMomentum(const GTreeParticle& tree, GH1* gHist)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			FillMissingMomentum(tree, i, j, gHist);
		}
	}
}

void PPhysics::FillMissingMomentum(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
        FillMissingMomentum(tree, particle_index, i, gHist);
	}
}

void PPhysics::FillMissingMomentum(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist)
{
    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

	// Fill GH1
	gHist->Fill(missingp4.Rho()/197.3,time);					

}

void PPhysics::FillDeltaE(const GTreeMeson& tree, GH1* gHist)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			FillDeltaE(tree, i, j, gHist);
		}
	}
}

void PPhysics::FillDeltaE(const GTreeMeson& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
        FillDeltaE(tree, particle_index, i, gHist);
	}
}

void PPhysics::FillDeltaE(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist)
{
    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    
    

    // Fill GH1
    gHist->Fill(CalcDeltaEDan(tree,particle_index,tagger_index),time);					

}

void PPhysics::FillDeltaE_Missmom(const GTreeMeson& tree, GH1** gHist)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
		  FillDeltaE_Missmom(tree, i, j, gHist);
		}
	}
}

void PPhysics::FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index, GH1** gHist)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
	  FillDeltaE_Missmom(tree, particle_index, i, gHist);
	}
}

void PPhysics::FillDeltaE_Missmom(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index, GH1** gHist)
{
  double qmin=0.0; //min and max in fm^-1
  double qmax=1.5;
  int nqbin = 150;
  int qbin;
  int nEbin=23;
  int Ebin=-1;
  double Ebin_v[] = {135,140,145,150,160,170,180,190,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,580};
    
    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    Double_t q=CalcMissingMomentumDan(tree,particle_index,tagger_index);
    Double_t E_beam=CalcBeamE(tagger_index);
    Double_t DeltaE=CalcDeltaEDan(tree,particle_index,tagger_index);

    for (int i=0; i<nEbin; i++) {
      if (E_beam >= Ebin_v[i] && E_beam<Ebin_v[i+1]) Ebin=i; 
    }
    qbin = int(TMath::Floor(q));
    if ( Ebin != -1 && qbin >0 && qbin< nqbin ) {
      gHist[qbin*nEbin+Ebin]->Fill(DeltaE,time);
    }
}


Double_t PPhysics::CalcMissingMomentum(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.Rho()/197.3;
}

void PPhysics::FillMissingMomentumDan(const GTreeParticle& tree, GH1* gHist)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			FillMissingMomentumDan(tree, i, j, gHist);
		}
	}
}

void PPhysics::FillMissingMomentumDan(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
        FillMissingMomentumDan(tree, particle_index, i, gHist);
	}
}

void PPhysics::FillMissingMomentumDan(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist)
{
    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    
    particle	= tree.Particle(particle_index);
    beam 		= TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(tagger_index),tagger->GetPhotonBeam_E(tagger_index));
  
  
    double theta_pi0 = particle.Theta(); // theta pi0
    double mpi0 = 134.9766;
    double costheta = TMath::Cos(theta_pi0);
    double beta = beam.E()/(beam.E()+target.M());
    double Egamma_c = beam.E()*TMath::Sqrt( (1-beta)/(1+beta) );
    double gamma_l = 1./(TMath::Sqrt(1-beta*beta));
    double minv = 2*beam.E()*target.M() + target.M2();
    double Epi_c = TMath::Sqrt(minv)/2. + (mpi0*mpi0 - target.M2())/(2.*TMath::Sqrt(minv));
    double Epi = (Epi_c + TMath::Sqrt( Epi_c*Epi_c - ( 1 - beta*beta*costheta*costheta )*( gamma_l*gamma_l*beta*beta*mpi0*mpi0*costheta*costheta + Epi_c*Epi_c ) ) ) / ( gamma_l*(1-beta*beta*costheta*costheta) );

    double qsq = (Egamma_c - Epi_c)*(Egamma_c - Epi_c) + 2.*beam.E()*(Epi - TMath::Sqrt(Epi*Epi - mpi0*mpi0)*costheta) - mpi0*mpi0;

    Double_t q = TMath::Sqrt(qsq)/197.3;

	// Fill GH1
	gHist->Fill(q,time);					

}


Double_t PPhysics::CalcBeamE(Int_t tagger_index){
  return tagger->GetPhotonBeam_E(tagger_index);
}


Double_t PPhysics::CalcMissingMomentumDan(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
  particle	= tree.Particle(particle_index);
  beam 		= TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(tagger_index),tagger->GetPhotonBeam_E(tagger_index));
  
  
  double theta_pi0 = particle.Theta(); // theta pi0
  double mpi0 = 134.9766;
  double costheta = TMath::Cos(theta_pi0);
  double beta = beam.E()/(beam.E()+target.M());
  double Egamma_c = beam.E()*TMath::Sqrt( (1-beta)/(1+beta) );
  double gamma_l = 1./(TMath::Sqrt(1-beta*beta));
  double minv = 2*beam.E()*target.M() + target.M2();
  double Epi_c = TMath::Sqrt(minv)/2. + (mpi0*mpi0 - target.M2())/(2.*TMath::Sqrt(minv));
  double Epi = (Epi_c + TMath::Sqrt( Epi_c*Epi_c - ( 1 - beta*beta*costheta*costheta )*( gamma_l*gamma_l*beta*beta*mpi0*mpi0*costheta*costheta + Epi_c*Epi_c ) ) ) / ( gamma_l*(1-beta*beta*costheta*costheta) );

  double qsq = (Egamma_c - Epi_c)*(Egamma_c - Epi_c) + 2.*beam.E()*(Epi - TMath::Sqrt(Epi*Epi - mpi0*mpi0)*costheta) - mpi0*mpi0;

  Double_t q = TMath::Sqrt(qsq)/197.3;
    
  return q;
  
}

Double_t PPhysics::CalcDeltaEDan(const GTreeMeson& tree, Int_t particle_index, Int_t tagger_index)
{
  particle	= tree.Particle(particle_index);
  beam 		= TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(tagger_index),tagger->GetPhotonBeam_E(tagger_index));
 
  int n_sub_part;
  int n_sub_phot;
  Double_t  E1, E2;
  Double_t E_diff=10000; // Value out of range in case one does not have just 2 photons

  // number of meson is ever 1 on the tree
  n_sub_part = tree.GetNSubParticles(0);
  n_sub_phot = tree.GetNSubPhotons(0);
  TLorentzVector gamma[n_sub_phot];
  if ( n_sub_phot ==2) {
    for (int i=0; i< n_sub_phot; i++) {
      gamma[i] = tree.SubPhotons(0,i);
    }

    E1 = gamma[0].E(); // photon 1 energy
    E2 = gamma[1].E(); // photon 2 energy
    Double_t beta = (beam.E()/(beam.E() + target.M())); 
    Double_t lorentz_gamma =  1/(sqrt(1 - beta*beta));
    Double_t Xform = (E1 - E2)/(E1 + E2);
    Double_t psi = gamma[0].Vect().Angle(gamma[1].Vect());
    Double_t costheta1 = gamma[0].CosTheta();
    Double_t costheta2 = gamma[1].CosTheta();
    Double_t mpi0 = 134.9766;
    Double_t M = target.M();
    Double_t Egamma=beam.E();
    E_diff = lorentz_gamma*((sqrt(2*mpi0*mpi0/((1-Xform*Xform)*(1-cos(psi))))) 
			    - ( beta*(E1*costheta1 + E2*costheta2)) ) 
      - ( (2*Egamma*M + mpi0*mpi0)/(2*sqrt(2*Egamma*M + M*M)) );
  }
  return E_diff;
  
}


void PPhysics::FillMissingMass(const GTreeParticle& tree, GH1* gHist)
{
	for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
			FillMissingMass(tree, i, j, gHist);
		}
	}
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
    for (Int_t i = 0; i < tagger->GetNTagged(); i++)
	{
        FillMissingMass(tree, particle_index, i, gHist);
	}
}

void PPhysics::FillMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index, GH1* gHist)
{
    // calc particle time diff
    time = tagger->GetTagged_t(tagger_index) - tree.GetTime(particle_index);
    
    // calc missing p4
    missingp4 = CalcMissingP4(tree, particle_index,tagger_index);

	// Fill GH1
	gHist->Fill(missingp4.M(),time);					

}



Double_t PPhysics::CalcMissingMass(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree, particle_index, tagger_index);

	return missingp4.M();
}

Double_t PPhysics::CalcMissingEnergy(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    missingp4 	= CalcMissingP4(tree,particle_index, tagger_index);

	return missingp4.T();
}

TLorentzVector PPhysics::CalcMissingP4(const GTreeParticle& tree, Int_t particle_index, Int_t tagger_index)
{
    particle	= tree.Particle(particle_index);
    beam 		= TLorentzVector(0.,0.,tagger->GetPhotonBeam_E(tagger_index),tagger->GetPhotonBeam_E(tagger_index));
	missingp4 	= beam + target - particle;						

	return missingp4;
}

void PPhysics::FillTime(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
            time = tagger->GetTagged_t(j) - tree.GetTime(i);
			gHist->Fill(time);
		}
	}
}

void PPhysics::FillTime(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
	for (Int_t j = 0; j < tagger->GetNTagged(); j++)
	{
		time = tagger->GetTagged_t(j) - tree.GetTime(particle_index);
		gHist->Fill(time);
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
        for (Int_t j = 0; j < tagger->GetNTagged(); j++)
		{
            		time = tagger->GetTagged_t(j) - tree.GetTime(i);
			if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
		}
	}
}

void PPhysics::FillTimeCut(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
	for (Int_t j = 0; j < tagger->GetNTagged(); j++)
	{
		time = tagger->GetTagged_t(j) - tree.GetTime(particle_index);
		if((GHistBGSub::IsPrompt(time)) || (GHistBGSub::IsRandom(time))) gHist->Fill(time);
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, GH1* gHist)
{
    for (Int_t i = 0; i < tree.GetNParticles(); i++)
	{
		gHist->Fill(tree.Particle(i).M());
	}
}

void PPhysics::FillMass(const GTreeParticle& tree, Int_t particle_index, GH1* gHist)
{
	gHist->Fill(tree.Particle(particle_index).M());
}

Bool_t 	PPhysics::Write()
{
	return kTRUE;
}
#endif
