// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "exceptions.h"
#include "glossary.h"
#include "simul_prop.h"
#include "actor_prop.h"
#include "actor.h"


Hand * ActorProp::newHand(HandMonitor* m) const
{
    return new Actor(this, m);
}


void ActorProp::clear()
{
    HandProp::clear();
    
    rescue_chance = 0;
    track_end             = NO_END;
    bind_only_growing_end = false;
}


void ActorProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    glos.set(rescue_chance, "rescue_prob", "rescue_probability");
    glos.set(rescue_chance, "rescue_chance");
    
    glos.set(track_end,  "track_end", {{"off",       NO_END},
#if BACKWARD_COMPATIBILITY < 100
                                       {"none",      NO_END},
#endif
                                       {"plus_end",  PLUS_END},
                                       {"minus_end", MINUS_END},
                                       {"both_ends", BOTH_ENDS}});

    glos.set(bind_only_growing_end, "bind_only_growing_end");
}


void ActorProp::complete(Simul const& sim)
{
    HandProp::complete(sim);
    
    if ( rescue_chance < 0 )
        throw InvalidParameter("rescuer:rescue_chance must be >= 0");
}


void ActorProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "rescue_chance", rescue_chance);
    write_value(os, "track_end",             track_end);
    write_value(os, "bind_only_growing_end", bind_only_growing_end);
}

