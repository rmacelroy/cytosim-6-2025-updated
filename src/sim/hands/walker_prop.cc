// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "dim.h"
#include "cymdef.h"
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "walker_prop.h"
#include "walker.h"
#include "simul_part.h"


Hand * WalkerProp::newHand(HandMonitor* m) const
{
    return new Walker(this, m);
}


void WalkerProp::clear()
{
    DigitProp::clear();

    stall_force       = 0;
    unloaded_speed    = 0;
    unbinding_chance  = 0;
    walking_rate_dt   = 0;
    var_rate_dt       = 0;
}


void WalkerProp::read(Glossary& glos)
{
    DigitProp::read(glos);
    
    glos.set(stall_force,    "stall_force", "force");
    glos.set(unloaded_speed, "unloaded_speed", "speed");
#if BACKWARD_COMPATIBILITY < 100
    if ( glos.set(unloaded_speed, "max_speed") )
        Cytosim::warn("'max_speed' is deprecated: use 'unloaded_speed'\n");
#endif
    glos.set(unbinding_chance, "unbinding_chance", 0, "unbinding", 2);

    if ( glos.has_key("dangling_chance") )
        Cytosim::warn ("use `hold_growing_end` instead of `dangling_chance`\n");
    
#if BACKWARD_COMPATIBILITY < 100
    if ( glos.set(hold_growing_end, 2, "hold_fiber") )
        Cytosim::warn("use hand:hold_growing_end instead of hand:hold_fiber\n");
#endif
}


void WalkerProp::complete(Simul const& sim)
{
    DigitProp::complete(sim);
   
    if ( primed(sim) && stall_force <= 0 )
        throw InvalidParameter("walker:stall_force must be > 0");
    
    if ( unbinding_chance < 0 )
        throw InvalidParameter("walker:unbinding_chance must be >= 0");

    if ( unbinding_chance > 1 )
        throw InvalidParameter("walker:unbinding_chance must be <= 1");
    
    walking_rate_dt = time_step(sim) * abs_real(unloaded_speed) / step_size;
    var_rate_dt = std::copysign(walking_rate_dt/stall_force, unloaded_speed);
}


void WalkerProp::checkStiffness(real stiff, real len, real mul, real kT) const
{
    DigitProp::checkStiffness(stiff, len, mul, kT);

#if ( 0 )
    /*
     Compare mobility with stiffness: this can induce instability
     */
    real ef = abs_real(set_speed_dt) * stiff * mul / stall_force;
    if ( unloaded_speed  &&  ef > 0.5 )
    {
        Cytosim::warn("simulating `", name(), "' may fail as:\n",
                      PREF, "time_step * stiffness * unloaded_speed / stall_force = ", ef, '\n',
                      PREF, "-> reduce time_step (really!)\n";
        //throw InvalidParameter(oss.str());
    }
    
    /*
     Compare the energy in a link due to the equipartition theorem
     to the maximum force that the motor can sustain before detaching:
     1/2 kT * DIM  <<  1/2 stiffness x^2 ~ 1/2 force^2 / stiffness;
     */
    if ( std::sqrt( DIM * kT * stiff ) > stall_force )
    {
        Cytosim::warn("The stall force of `", name(), "' is too small:\n",
                      PREF, "DIM * kT * stiffness > stall_force\n",
                      PREF, "-> reduce stiffness or increase stall_force\n");
    }
    
    /*
     Compare the force created by traveling during the time 1/unbinding_rate,
     and compare to stall_force. This is limit the efficiency of the motor.
     */
    ef = abs_real( stiff * unloaded_speed / ( unbinding_rate * stall_force ));
    if ( unbinding_rate && unloaded_speed  &&  ef < 1 )
    {
        Cytosim::warn("Walker `", name(), "' is inefficient since\n", PREF,
                      "stiffness * unloaded_speed / unbinding_rate << stall_force",
                      "; ratio = ", ef, "\n");
    }
    
    
    /*
     Compare the force reached in one step with the stall force
     */
    if ( abs_real( step_size * stiff ) > 0.5 * stall_force )
        Cytosim::warn("attention:  stiffness * digit:step > stall_force / 2\n");
#endif
}


void WalkerProp::write_values(std::ostream& os) const
{
    DigitProp::write_values(os);
    write_value(os, "stall_force",      stall_force);
    write_value(os, "unloaded_speed",   unloaded_speed);
    write_value(os, "unbinding_chance", unbinding_chance);
}

