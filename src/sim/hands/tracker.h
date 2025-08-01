// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef TRACKER_H
#define TRACKER_H

#include "hand.h"
#include "tracker_prop.h"

/// A Hand that only binds near the end of a Fiber
/**
 The Tracker is a Hand, and thus can bind and unbind from fibers,
 but it can only bind at a certain distance from the end of the fiber.
  
 if @ref TrackerPar "bind_only_growing_end" is set, 
 binding will occur only if the fiber end is growing
 
 The parameter @ref TrackerPar "track_end" can
 specify that the Tracker will always be attached at the end of the filament.
 
 See Examples and the @ref TrackerPar.
 @ingroup HandGroup
 */
class Tracker : public Hand
{
private:
    
    /// disabled default constructor
    Tracker();

public:
    
    /// Property
    TrackerProp const* prop() const { return static_cast<TrackerProp const*>(Hand::prop); }

    /// constructor
    Tracker(TrackerProp const*, HandMonitor* h);
    
    
    /// attach the hand at the position described by site, or return false
    bool attachmentAllowed(FiberSite&) const;
    
    /// simulate when `this` is attached but not under load
    void stepUnloaded();
    
    /// simulate when `this` is attached and under load
    void stepLoaded(Vector const& force);
    
};

#endif

