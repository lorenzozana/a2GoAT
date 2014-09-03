#ifndef __PPi0Example_h__
#define __PPi0Example_h__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string> 

#include "GTreeManager.h"
#include "PPhysics.h"
#include "TH2F.h"

class	PPi0Example  : public PPhysics
{
private:
    GH1*	time;
    GH1*	time_cut;
    GH1*	time_2g;      
    GH1*	time_2g_cut;   
     
    GH1*	IM;
    GH1*	IM_2g;
    
    GH1*	MMass;
    GH1*	MMass_2g; 

    GH1*	MMom;
    GH1*	MMom_2g; 

    GH1*	MMomDan;
    GH1*	MMomDan_2g; 

    GH1*	DeltaE_CM_Dan_2g; 
    GH1*       DeltaE_Missmom_BeamE[3500];
protected:
    virtual Bool_t  Start();

    virtual void    ProcessEvent();
    virtual void	ProcessScalerRead();
			
public:
    PPi0Example();
    virtual ~PPi0Example();

    //virtual Bool_t	Init(const char* configfile);

};
#endif
