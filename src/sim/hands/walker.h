// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef WALKER_H
#define WALKER_H

#include "digit.h"
#include "walker_prop.h"


/// A Hand that move with discrete steps of fixed size along a Fiber
/**
 The Walker moves like a Motor, but it makes discrete jumps along the fiber.
 
 The size of the step size is set by @ref DigitPar digit:step_size.
 The number of steps in one `time_step` is a stochastic integer  
 with a Poisson distribution.
 
 As defined in Hand, detachment increases exponentially with force.

 The Digit will use the Lattice and will not step into a site,
 if this site already has a Hand bound with a overlapping `footprint`.
 
 See Examples and the @ref WalkerPar.
 @ingroup HandGroup 

 Note: Conventional kinesin can make backward steps under heavy load:\n
 Carter, N. & Cross, R. Mechanics of the kinesin step. Nature 435, 308–312 (2005).
 http://dx.doi.org/doi:10.1038/nature03528
 
 According to this work, a real kinesin stalls because the probability of making
 a forward step, which is decreased by load, becomes equal to the probability of
 making a backward step.

 With antagonistic force however, the Walker will stall in a immobile configuration,
 because it can only make forward steps, wereas a stalled kinesin keeps moving
 back-and-forth over a few lattice sites.
 Thus one should be careful when using Walker to model Kinesin.
*/
class Walker : public Digit
{
private:

    /// disabled default constructor
    Walker();

public:
    
    /// Property
    WalkerProp const* prop() const { return static_cast<WalkerProp const*>(Hand::prop); }

    /// constructor
    Walker(WalkerProp const*, HandMonitor*);
    
    
    /// attach and update variables
    void attach(FiberSite const&);
    
    /// simulate when `this` is attached but not under load
    void stepUnloaded();
    
    /// simulate when `this` is attached and under load
    void stepLoaded(Vector const& force);
    
};

#endif

