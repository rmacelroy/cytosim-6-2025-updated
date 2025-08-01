// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
#ifndef RESCUER_PROP_H
#define RESCUER_PROP_H

#include "hand_prop.h"


/// Additional Property for Rescuer
/**
 @ingroup Properties
 */
class RescuerProp : public HandProp
{
    friend class Rescuer;
    
public:
    
    /**
     @defgroup RescuerPar Parameters of Rescuer
     @ingroup Parameters
     Inherits @ref HandPar.
     @{
     */
    
    /// probability of rescuing a Fiber
    /**
     This parameter is used when the shrinking end of the fiber is reaching the position of the Rescuer.
     
     The parameter determines between two possible outcome:
     - the fiber nearest end state is set to STATE_GREEN.
     - the fiber continues to shrink, and the Rescuer detaches
     .
     */    
    real rescue_chance;
    
    /// @}
    
public:
    
    /// constructor
    RescuerProp(const std::string& n) : HandProp(n)  { clear(); }
    
    /// destructor
    ~RescuerProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new RescuerProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

