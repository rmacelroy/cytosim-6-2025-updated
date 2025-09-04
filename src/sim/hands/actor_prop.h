// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ACTOR_PROP_H
#define ACTOR_PROP_H

#include "hand_prop.h"


/// Additional Property for Actor
/**
 @ingroup Properties
 */
class ActorProp : public HandProp
{
    friend class Actor;
    
public:
    
    /**
     @defgroup ActorPar Parameters of Actor
     @ingroup Parameters
     Inherits @ref HandPar.
     @{
     */
    
    FiberEnd track_end;
    
    /// if true, bind only to growing fiber ends
    bool bind_only_growing_end;
    
    /// rate of event
    real rescue_chance;
    
    /// @}

public:
    
    /// constructor
    ActorProp(const std::string& n) : HandProp(n)  { clear(); }
    
    /// destructor
    ~ActorProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor*) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new ActorProp(*this); }
    
    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

