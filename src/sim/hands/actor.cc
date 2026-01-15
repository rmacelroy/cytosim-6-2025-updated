// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "actor.h"
#include "actor_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul_part.h"
#include "hand_monitor.h"

//------------------------------------------------------------------------------

Actor::Actor(ActorProp const* p, HandMonitor* h)
: Hand(p,h)
{
    //throw InvalidParameter("the actor class is unfinished");
}

//------------------------------------------------------------------------------


void Actor::handleDisassemblyM()
{
    Fiber * fib = modifiableFiber();

    if ( hAbs < fib->abscissaM() && fib->freshAssemblyM() < 0 )
    {
        if ( RNG.test(prop()->rescue_chance) )
        {
            //revert the last disassembly step:
            fib->undoGrowM();
            // induce rescue:
            fib->setEndStateM(STATE_GREEN);
        }
        else
            detach();
    }
}

void Rescuer::handleDisassemblyP()
{
    Fiber * fib = modifiableFiber();

    /* increase MT length to cover position of Hand after checking that increment
    would be positive, since another Hand might have done the same job already. */
    if ( hAbs > fib->abscissaP() && fib->freshAssemblyP() < 0 )
    {
        if ( RNG.test(prop()->rescue_chance) )
        {
            //revert the last disassembly step:
            fib->undoGrowP();
            // induce rescue:
            fib->setEndStateP(STATE_GREEN);
            //assert_true(hAbs <= fib->abscissaP())
        }
        else
            detach();
    }
}

void Actor::stepUnloaded()
{
    assert_true( attached() );
	
    real a = hAbs + prop()->line_diffusion_dt * RNG.sreal();
    
    const real M = hFiber->abscissaM();
    const real P = hFiber->abscissaP();
    
    if ( a <= M )
    {
        if ( RNG.test_not(prop()->hold_growing_end[1]) )
            return detach();
        a = M;
    }
    
    if ( a >= P )
    {
        if ( RNG.test_not(prop()->hold_growing_end[0]) )
            return detach();
        a = P;
    }

    moveTo(a);
    
}


void Actor::stepLoaded(Vector const& force)
{
    assert_true( attached() );
    
    real load = dot(force, dirFiber());
    real diff = prop()->line_diffusion_dt * RNG.sreal();
    real a = hAbs + diff + prop()->movability_dt * load;
    
    const real M = hFiber->abscissaM();
    const real P = hFiber->abscissaP();
    
    if ( a <= M )
    {
        if ( RNG.test_not(prop()->hold_growing_end[1]) )
            return detach();
        a = M;
    }
    
    if ( a >= P )
    {
        if ( RNG.test_not(prop()->hold_growing_end[0]) )
            return detach();
        a = P;
    }

    moveTo(a);

}

