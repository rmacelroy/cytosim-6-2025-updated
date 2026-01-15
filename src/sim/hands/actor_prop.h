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
	
   /// mobility coefficient (default=0)
   /*
    The displacement of a slider is proportional to the projected force component:

        force_parallel_component = force_vector . direction_of_fiber
        speed = mobility * force_parallel_component
        abscissa += time_step * speed
    
    The mobility has unit of speed per force = um . s^-1 . pN^-1
    
    In principle, the ‘mobility’ is related to the diffusion constant of the molecule
    along the filament lattice by Einstein’s relation: `diffusion = mobility * kT`.
    So measuring the 1D diffusion constant of the molecules bound to a filament
    may provide an estimate of the mobility.
    */
   real movability;
   
   /// unidimensional diffusion coefficient while bound to a Fiber
   real line_diffusion;

   /// stiffness used to calculate the mobility associated with passive movements
   /**
    If this parameter is set, `mobility` is calculated as:

        mobility = 1.0 / ( stiffness * time_step(sim) );

    Corresponding to the maximum mobility possible, since the Hand will
    move in one time-step to the position where the force originates.
    
    The stiffness needs to be set from the link stiffness `S`,
    and considering the mobility of the other link.
    
    Generally you may use:
    -# stiffness = S if the Hand belongs to a Single, or if it is the only Slider in a Couple
    -# stiffness = 2 * S if the Hand is part of a Couple, made of two Slider
    .
    
    By default, this parameter is unset, and `mobility` is used unmodified.
    */
   real stiffness;
   
   /// @}

private:
    
    /// derived variable:
    real line_diffusion_dt, movability_dt;
	
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
	
    /// perform additional tests for the validity of parameters, given the elasticity
    void checkStiffness(real stiff, real len, real mul, real kT) const;
    
    /// return a carbon copy of object
    Property* clone() const { return new ActorProp(*this); }
    
    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

