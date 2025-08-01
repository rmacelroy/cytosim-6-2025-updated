// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FORK_PROP_H
#define FORK_PROP_H

#include "vector2.h"
#include "couple_prop.h"


/// Additional Property for Fork
/**
 @ingroup Properties
*/
class ForkProp : public CoupleProp
{
    friend class Fork;
    
public:
    
    /**
     @defgroup ForkPar Parameters of Fork
     @ingroup Parameters
     Inherits @ref CouplePar.
     @{
     */
    
    /// Resting angle in radian (set as `torque[0]`)
    real rest_angle;
    
    /// Stiffness of the angular link, in Torque per radians (pN.um/radian) (set as `torque[1]`)
    real angular_stiffness;
    
    /// Allow the angle to flip in 2D
    bool flip;
    
    /// @}
    
    /// derived variable: [cos(angle), sin(angle)]
    Vector2 rest_dir;

public:
    
    /// constructor
    ForkProp(const std::string& n) : CoupleProp(n)  { clear(); }
    
    /// destructor
    ~ForkProp() { }
    
    /// return a Hand with this property
    Couple * newCouple() const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new ForkProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

