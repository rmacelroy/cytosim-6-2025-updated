# Filament Assembly Dynamics

This document provides an overview of available models covering filament growth and shrinkage. 
They are accessible by the `activity` parameter, and implemented in separate classes, derived 
from the `fiber` class. 

Generally, names of parameters are similar between the different classes: `growing_speed`, `growing_force`, `max_length`, `total_polymer`.

Note that Vertices may be added automatically to a filament as it grow, and similarly removed upon shrinking.
All points on a filament remain at any time regularly distributed and equidistant. The number of points is always such that the distance separating adjacent points is as close as possible to the `segmentation` parameter.

# 1. Growing filament

Use the filament's activity `grow` to make it grow smoothly, up to a maximum length:

    set fiber microtubule
    {
        ...
        activity = grow
        growing_speed = 1
        growing_force = 3
        ...
    }

 The basic parameters are:

 * `growing_speed`, the base assembly rate in um/s.
 * `growing_force`, the characteristic force of polymerization in pN.

Positive values of the rate correspond to assembly, and negative values to disassembly. 
The value of `growing_speed` must be specified, but by default `growing_force = inf`. 
Thus by default, growth will occur at constant speed, but if `growing_force` is specified, the growth will depend on force (see below).

### Initial state

For all type of dynamic filaments, one needs to specify the initial state in the 'new' command where they are created, for example:

    new microtubule
    {
        length = 1
        plus_end = grow
    }

The exact possible values for `plus_end` would depends on the class, but in general `static`, `grow` and `shrink` should be supported. The ends are static by default.

### Defining a maximum length

Additional parameters can be defined:

    set fiber microtubule
    {
        ...
        max_length = 8
        ...
    }

In this case, the length of a fiber will not exceed `max_length` (which by default is infinite).

### Assembly can depend on available monomers
    
To simulate a finite pool of monomer, define the parameter `total_polymer`:

    set fiber microtubule
    {
        ...
        total_polymer = 100
        ...
    }

In this case, growth rate will be linearly dependent on the availability of polymer:

    speed = growing_speed * free_polymer
         
In this equation, `free_polymer` is a number in [0,1], representing the fraction of available monomers:

    free_polymer = 1.0 - sum(all_fiber_length) / total_polymer

This parameter does not need to be defined, and by default, `total_polymer = inf` and `free_polymer = 1`.

### Assembly can depend on force

You can make assembly exponentially dependent on force.  
Disassembly will always occurs at the specified rate.
For assembly, only the component of the force parallel to the direction of the fiber at the end is taken into account:
 
     force = dot(force_vector, fiber_direction);
     
The projected force is negative ( antagonistic ) if it is directed against fiber assembly.

     if ( force < 0 )
         speed = growing_speed * free_polymer * exp( force / growing_force ) + growing_off_speed;
     else
         speed = growing_speed * free_polymer + growing_off_speed;

### Assembly at the minus-end

Assembly can be defined for the Minus end as well by providing a second value to these parameters:

    set fiber microtubule
    {
        ...
        growing_speed = 1, 1
        growing_force = 3, 3
        ...
    }

Here the minus-end (second column) will assembles at the same speed as the plus end (first column).

# 2. Treadmilling filament

 This is not documented yet!

# 3. Standard microtubule dynamic instability

Filaments can grow and shrink according to the classical two-states dynamic instability model as defined by T. Hill.
It is in many ways an extension of the Growing filament model above.
The parameter `activity = classic` turns this feature on:

    set fiber microtubule
    {
        ...
        activity = classic
        growing_speed    = 0.06
        shrinking_speed  = -0.15
        catastrophe_rate = 0.004,  0.04
        rescue_rate      = 0
        growing_force    = 1.67
        total_polymer    = 100
        min_length       = 0.5
        ...
    }

 The basic parameters are:

 * `growing_speed`, the base assembly rate in um/s (must be non-negative),
 * `shrinking_speed`, the base disassembly rate in um/s (must be non-positive),
 * `growing_force`, the characteristic force of polymerization in pN,
 * `castastrophe_rate`, the rate of stochastic switch from the assembly state to the disassembly state,
 * `rescue_rate`, the rate of stochastic switch from the disassembly state to the assembly state,
 * `min_length`, the length below which the filament stops shriking
 
The length of a fiber will not exceed `max_length` (which by default is infinite), and any Fiber shorter than `min_length` will be deleted.

In this model, the lengths can be exponentially distributed, and in this case the mean length is

    average_length = growing_speed / catastrophe_rate
    
For more information, check:

    Physical aspects of the growth and regulation of microtubule structures.
    Dogterom, M. & Leibler, S.; Phys Rev Lett 70, 1347–1350 (1993).

### Assembly depends on force and available monomers

Assembly is exponentially decreased by antagonistic force, and linearly dependent on the availability of polymer.  
Disassembly always occurs at the specified rate.
Only the component of the force parallel to the direction of the fiber at the end is taken into account:
 
     force = force_vector * fiber_direction;
     
The projected force is negative ( antagonistic ) if it is directed against fiber assembly.

     if ( force < 0 )
         speed = growing_speed * free_polymer * exp( force / growing_force ) + growing_off_speed;
     else
         speed = growing_speed * free_polymer + growing_off_speed;
         
In this equation, `free_polymer` is a number in [0,1], representing the fraction of free monomers.
It is defined as:

    free_polymer = 1.0 - sum(all_fiber_length) / total_polymer

### Catastrophe rate depends on growth speed

If assembly is reduced following the equations above, this also increase the catastrophe rate.
This effect was described in:

>[Measurement of the Force-Velocity Relation for Growing Microtubules](http://www.sciencemag.org/content/278/5339/856.abstract)  
>Marileen Dogterom and Bernard Yurke; Science Vol 278 pp 856-860; 1997
 
>[Dynamic instability of MTs is regulated by force](http://www.jcb.org/cgi/doi/10.1083/jcb.200301147)  
> M.Janson, M. de Dood, M. Dogterom; Journal of Cell Biology Vol 161, Nb 6, 2003
 Figure 2 C

The two values of the castrophe rate are:

1. rate for a microtubule that is freely growing
2. rate for a microtubule that is fully stalled

The second value should be set higher than the first one.
In this case, it corresponds to a delay of 25 seconds:

    set fiber microtubule
    {
        ...
        activity = classic
        catastrophe_rate = 0.004,  0.04
        ...
    }

### Setting the initial state of the ends

The ends of a dynamic filament can have three possible states:

* static = white
* growing = green
* shrinking = red

These correspond to the default color in which they are displayed. 

It is possible to set the state of the ends of a filament when it is created:

    new microtubule
    {
        ...
        end_state = green, white 
        ...
    }

The first value applies to the plus end and the second value to the minus end.

While playing with this features, press `r` to display information about the microtubules lengths and states while displaying the simulation.


# 4. Discrete Microtubule Dynamic Instability 

Microtubules can grow and shrink according to a four-state dynamic instability model:

> [A theory of microtubule catastrophes and their regulation](http://www.pnas.org/content/106/50/21173)  
> Brun L, Rupp B, Ward J, Nedelec F; PNAS 106 (50) 21173-21178; 2009

Set the parameter `activity = dynamic` to enable this model on:

    ...    
    set fiber microtubule
    {
        ...
        activity        = dynamic
        unit_length     = 0.008
        growing_speed   = 0.06
        shrinking_speed = -0.15
        hydrolysis_rate = 0.06
        growing_force   = 1.7
        ...
    }
    ...

 The basic parameters are:

 * `unit_length`, the effective size of the monomers,
 * `growing_speed`, the base assembly rate in um/s (must be non-negative),
 * `shrinking_speed`, the base disassembly rate in um/s (must be non-positive),
 * `growing_force`, the characteristic force of polymerization in pN,
 * `hydrolysis_rate`, the rate of GTP hydrolysis in the microtubule lattice.

In this model the length of the filament is always a multiple of `unit_length`.
The average growing/shrinking speed is set by `growing_speed` and `shrinking_speed`, respectively. At each time step, the increment of length is a multiple of `unit_length` and is a stochastic number following a Poisson law. Cytosim uses the Gillespie scheme to calculate these increments.

In this model, the effective catastrophe rate is set indirectly by `hydrolysis_rate`.


### Authors

Francois Nedelec (14.09.2017)

