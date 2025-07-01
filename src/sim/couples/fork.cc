// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fork.h"
#include "fork_prop.h"
#include "meca.h"
#include <iostream>
#include <fstream>
#include <random>


//m and s are the mean and std dev of the NATURAL LOG of a
//distribution of angles (by Akanksha) in RADIANS
double m = -2.29;
double s = 1.16; //these are from Akanksha's angle data

std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<> d(m, s);

Fork::Fork(ForkProp const* p, Vector const& w)
: Couple(p, w)
{
#if ( DIM == 2 )
    sine = prop()->rest_dir.YY;
#endif
	//sinus = 0;
	equilibrium_angle = exp(d(gen));
}


Fork::~Fork()
{
}


void Fork::setInteractions(Meca& meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();
    
    meca.addLink(pt1, pt2, prop()->stiffness);
    
#if ( DIM == 2 )
    if ( prop()->flip )
    {
        Vector2 dir = prop()->rest_dir;
        // flip the angle to match the current configuration of the bond
        sine = std::copysign(dir.YY, cross(pt1.diff(), pt2.diff()));
        dir.YY = sine;
        meca.addTorque(pt1, pt2, dir, prop()->angular_stiffness);
    }
    else
    {
        meca.addTorque(pt1, pt2, prop()->rest_dir, prop()->angular_stiffness);
    }
    //meca.addTorquePoliti(pt1, pt2, dir, prop()->angular_stiffness);
#elif ( DIM == 3 )
	
	//set the resting angle?
	//prop()->rest_dir.XX = cos(equilibrium_angle);
    //prop()->rest_dir.YY = sin(equilibrium_angle);
	
    meca.addTorque(pt1, pt2, prop()->rest_dir, prop()->angular_stiffness);
#endif
}

