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

bool Actor::attachmentAllowed(FiberSite& sit) const
{
    if ( !Hand::attachmentAllowed(sit) )
        return false;

    // check if fiber end in growing:
    if ( prop()->bind_only_growing_end && !sit.fiber()->isGrowing(sit.nearestEnd()) )
        return false;
    
    return true;
}


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
    
    switch ( prop()->track_end )
    {
        case NO_END:
            break;
            
        case PLUS_END:
            relocateP();
            break;
            
        case MINUS_END:
            relocateM();
            break;
            
        case BOTH_ENDS:
            moveToEnd(nearestEnd());
            break;
            
        default:
            throw InvalidParameter("invalid value of tracker:track_end");
    }
}


void Actor::stepLoaded(Vector const& force)
{
    assert_true( attached() );

    switch ( prop()->track_end )
    {
        case NO_END:
            break;
            
        case PLUS_END:
            relocateP();
            break;
            
        case MINUS_END:
            relocateM();
            break;
            
        case BOTH_ENDS:
            moveToEnd(nearestEnd());
            break;
            
        default:
            throw InvalidParameter("invalid value of tracker:track_end");
    }
}

