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
    movability = 0;
    stiffness = -1;
    line_diffusion = 0;
}


void ActorProp::read(Glossary& glos)
{
    HandProp::read(glos);
    
    glos.set(rescue_chance, "rescue_prob", "rescue_probability");
    glos.set(rescue_chance, "rescue_chance");
    
    glos.set(movability, "mobility");
    glos.set(stiffness, "stiffness");
    glos.set(line_diffusion, "diffusion");
}


void ActorProp::complete(Simul const& sim)
{
    HandProp::complete(sim);
    
    if ( rescue_chance < 0 )
        throw InvalidParameter("rescuer:rescue_chance must be >= 0");
	
    if ( movability < 0 )
        throw InvalidParameter(name()+"mobility must be >= 0");
    
    if ( primed(sim) && movability <= 0 )
        std::clog << "WARNING: slider `" << name() << "' will not slide because mobility=0\n";

    if ( line_diffusion < 0 )
        throw InvalidParameter(name()+":diffusion must be >= 0");
    
    // use Einstein's relation to get a mobility:
    if ( line_diffusion == 0 )
    {
        line_diffusion = movability * boltzmann(sim);
        Cytosim::log(name(), ":diffusion <--- ", line_diffusion, "\n");
    }
    
    /*
     This is for unidimensional diffusion along the filaments, and we want:
     var(dx) = 2 D time_step, given that we use dx = diffusion_dt * RNG.sreal()
     Since `sreal()` is uniformly distributed in [-1, 1], its variance is 1/3,
     and we need `diffusion_dt^2 = 6 D time_step`
     */
    line_diffusion_dt = std::sqrt(6.0 * line_diffusion * time_step(sim));

    /*
     Explicit
     */
    
    movability_dt = time_step(sim) * movability;
    
    if ( stiffness > 0 )
    {
        /*
         We devise here an implicit integration approach, assuming:
         - that all other elements of the simulation are static
         - that the link is Hookean of zero resting length:
         force = stiffness * offset
         However, this is true only if the Slider is part of a Single or a plain Couple.
         This does not hold in particular for any of the non-zero resting length Couple or Single.
         J. Ward found that in this case, the numerical precision is not improved compared to
         the explicit integration above.
         */
        std::clog << "         slider:mobility explicit = " << movability_dt;
        movability_dt = -std::expm1( - movability_dt * stiffness ) / stiffness;
        std::clog << "   implicit = " << movability_dt << '\n';
    }
}


void ActorProp::checkStiffness(real stiff, real len, real mul, real kT) const
{
    HandProp::checkStiffness(stiff, len, mul, kT);
    
    /*
     Estimate numerical stability from mobility and stiffness
     */
    real e = movability_dt * stiff * mul;
    if ( e > 2.0 + REAL_EPSILON )
    {
        std::ostringstream oss;
        oss << "Actor `" << name() << "' is unstable since:\n";
        oss << PREF << "time_step * mobility * stiffness = " << e << '\n';
        oss << PREF << "-> reduce mobility or time_step\n";
        throw InvalidParameter(oss.str());
        //std::clog << oss.str();
    }
    else
        Cytosim::log("   Actor `", name(), "' stability = ", e, "\n");
}


void ActorProp::write_values(std::ostream& os) const
{
    HandProp::write_values(os);
    write_value(os, "rescue_chance", rescue_chance);
    write_value(os, "mobility",  movability);
    write_value(os, "stiffness", stiffness);
    write_value(os, "diffusion", line_diffusion);
}

