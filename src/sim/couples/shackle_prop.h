// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SHACKLE_PROP_H
#define SHACKLE_PROP_H

#include "couple_prop.h"


/// Additional Property for Shackle
/**
 @ingroup Properties
*/
class ShackleProp : public CoupleProp
{
    friend class Shackle;
    
public:
    
    /**
     @defgroup ShacklePar Parameters of Shackle
     @ingroup Parameters
     Inherits @ref CouplePar.
     @{
     */
    
    
    /// @}

public:
    
    /// constructor
    ShackleProp(const std::string& n) : CoupleProp(n)  { clear(); }
    
    /// destructor
    ~ShackleProp() { }
    
    /// return a Hand with this property
    Couple * newCouple() const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new ShackleProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

