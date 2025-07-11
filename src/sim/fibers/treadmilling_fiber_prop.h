// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef TREADMILLING_FIBER_PROP
#define TREADMILLING_FIBER_PROP

#include "fiber_prop.h"


/// additional Property for TreadmillingFiber
/**
 @ingroup Properties
 Assembly is a continuous process, and can occur at both ends, typically:
 
      delta_length = growing_speed * time_step
 
 */
class TreadmillingFiberProp : public FiberProp
{
    friend class TreadmillingFiber;
    
public:
    
    /**
     @defgroup TreadmillingFiberPar Parameters of TreadmillingFiber
     @ingroup Parameters
     Inherits @ref FiberPar.
     @{
     */
    
    /// see @ref TreadmillingFiber

    /// Characteristic force for polymer assembly
    real growing_force[2];
    
    /// Speed of assembly
    real growing_speed[2];
    
    /// Speed of disassembly
    real shrinking_speed[2];
    
    /// @}
    
private:
    
    /// growing_speed * time_step
    real growing_speed_dt[2];
    
    /// 1.0 / growing_force
    real growing_force_inv[2];
    
    /// shrinking_speed * time_step
    real shrinking_speed_dt[2];
    
public:
    
    /// constructor
    TreadmillingFiberProp(const std::string& n) : FiberProp(n) { clear(); }

    /// destructor
    ~TreadmillingFiberProp() { }
    
    /// return a Fiber with this property
    Fiber* newFiber() const;
    
    /// set default values
    void clear();
       
    /// set using a Glossary
    void read(Glossary&);
   
    /// check and derive parameter values
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new TreadmillingFiberProp(*this); }

    /// write
    void write_values(std::ostream&) const;

};

#endif

