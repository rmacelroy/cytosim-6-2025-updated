// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2023

#include <iostream>
#include <numeric>
#include <list>
#include <set>
#include "tokenizer.h"
#include "organizer.h"
#include "matrix22.h"

#include "random_pcg.h"
using namespace PCG32;


/// width of columns in formatted output, in number of characters
static int column_width = 10;

/// use this macro at the beginning of a line of comment
#define COM "\n% " << std::setw(column_width-2)

/// use this macro at the beginning of a new line of data
#define LIN '\n' << std::setw(column_width)

/// use this macro to separate between values on a line of data
#define SEP ' ' << std::setw(column_width-1)

#include "accumulator.h"

/// add white-space right of `str` to reach length 'n * column_width - p' at most
static std::string ljust(std::string const& str, size_t n, size_t p = 0)
{
    size_t a = n * column_width;
    size_t s = a - std::min(a, p+str.size());
    return str + std::string(s, ' ');
}

/// add white-space left of `str` to reach length 'n * column_width - p' at most
static std::string rjust(std::string const& str, size_t n, size_t p = 1)
{
    size_t a = n * column_width;
    size_t s = a - std::min(a, p+str.size());
    return std::string(s, ' ') + str;
}

/// repeat string DIM times with X, Y, Z endings as appropriate
static std::string repeatXYZ(std::string const& str)
{
    std::string res(rjust(str+"X", 1));
#if ( DIM > 1 )
    res += " " + rjust(str+"Y", 1);
#endif
#if ( DIM > 2 )
    res += " " + rjust(str+"Z", 1);
#endif
    return res;
}

/// remove any 's' at the end of the argument
static void remove_trailing_s(std::string & str)
{
    if ( str.size() > 2  &&  str.at(str.size()-1) == 's' )
        str.resize(str.size()-1);
}

//------------------------------------------------------------------------------

/**
 combines multiple report, if `what` has multiple instructions separated by ','
 for example:
 
     report fiber:force,fiber:length
 
 parameters can be specifed as in:
 
     report fiber:cluster{couple=1},fiber:length

 */
void Simul::poly_report(std::ostream& out, std::string what, Glossary& opt, int frm) const
{
    int ver = 2;
    opt.set(ver, "verbose");
    if (( ver & 2 ) && frm >= 0 )
        out << "% frame   " << frm << '\n';
    std::stringstream is(what);
    while ( is.good() )
    {
        std::string arg = Tokenizer::get_polysymbol(is, false);
        std::string blk = Tokenizer::get_block(is, '{');
        try {
            if ( blk.empty() )
            {
                //out << "\nSimul::report(" << arg << ")";
                mono_report(out, arg, opt, ver);
            }
            else
            {
                //out << "\nSimul::report(" << arg << ", " << blk << ")";
                Glossary glos(opt);
                glos.read(blk, 0);
                mono_report(out, arg, glos, ver);
                opt.add_reads(glos);
            }
        }
        catch( Exception & e )
        {
            out << '\n' << e.brief();
        }
        // another report instruction should be separated by a comma:
        if ( is.peek() != ',' )
            break;
        is.get();
    }
    if ( ver & 2 )
    {
        out << "\n% end report\n\n";
    }
    {
        // check for unused characters in instruction stream
        std::string str;
        std::getline(is, str);
        if ( str.size() > 0 )
            throw InvalidParameter("unexpected `" + str + "' in report string");
    }
    //opt.print_counts(std::cerr);
}


/**
 Surround the report with comments to identify start/end
 */
void Simul::mono_report(std::ostream& out, std::string const& arg, Glossary& opt, int ver) const
{
    std::streamsize p = 4;
    opt.set(p, "precision");
    opt.set(column_width, "column") || opt.set(column_width, "width");

    // adjust floating-point notation:
    out.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::streamsize sp = out.precision(p);

    if ( ver & 1 )
    {
        //out << "% start\n";
        out << "% time " << std::to_string(time()) << '\n';
    }
    if ( ver > 0 )
    {
        out << "% report " << arg << " " << opt.to_string();
        report_one(out, out, arg, opt);
    }
    else
    {
        std::ofstream nul("/dev/null");
        report_one(out, nul, arg, opt);
    }
    if ( ver & 1 )
        out << "\n% end\n\n";
    out.precision(sp);
}


/**
 Split 'arg' into `who:what` and call report_one() accordingly
 Surround the report with comments to identify start/end
 */
void Simul::report_one(std::ostream& out, std::ostream& com, std::string const& arg, Glossary& opt) const
{
    std::string who = arg, what;
    
    // split argument string into who:what
    std::string::size_type pos = arg.find(':');
    if ( pos != std::string::npos )
    {
        who  = arg.substr(0, pos);
        what = arg.substr(pos+1);
    }
    
    // allow for 's' to be present or not at the ends of words:
    remove_trailing_s(who);
    remove_trailing_s(what);
    
    //std::clog << "report("<< who << "|" << what << ")\n";
    if ( isCategory(who) )
    {
        int split = false;
        if ( opt.set(split, "split") && split )
        {
            int veb = 1;
            opt.peek(veb, "verbose");
            // generate a separate report for all classes in this category:
            PropertyList plist = properties.find_all(who);
            std::ofstream nul("/dev/null");
            if ( veb && plist.size() )
            {
                com << COM << "split:";
                for ( Property const* sel : plist )
                    com << SEP << sel->name();
            }
            for ( Property const* sel : plist )
            {
                if ( veb ) {
                    report_one(out, com, who, sel, what, opt);
                    veb = 0;
                } else
                    report_one(out, nul, who, sel, what, opt);
            }
        }
        else
        {
            // pool all objects in this category:
            report_one(out, com, who, nullptr, what, opt);
        }
    }
    else
    {
        // check if `who` is the name of a property:
        Property const* sel = nullptr;
        sel = properties.find(who);
        if ( sel )
            report_one(out, com, sel->category(), sel, what, opt);
        else
            report_one(out, com, who, nullptr, what, opt);
    }
}


/**
 WHAT            | Output
 ----------------|--------------------------------------------------------------
 `bead`          | Position of beads
 `couple`        | Summary with number of couples in each state
 `fiber`         | Length and position of the ends of fibers
 `single`        | Number and state of singles
 `solid`         | Position of center and first point of solids
 `sphere`        | Position of center and first point of spheres
 `organizer`     | Position of the center of asters and other organizers
 `field`         | Total quantity of substance in field and Lattices
 
 
 WHAT               | Output
 -------------------|----------------------------------------------------------
 `simul:time`       | Time
 `simul:inventory`  | summary list of objects
 `simul:property`   | All properties
 `simul:parameter`  | global parameters
 `simul:NAME`       | parameters for Property 'NAME'


 WHAT                    | Output
 ------------------------|------------------------------------------------------
 `fiber:position`        | Position and orientation of fibers
 `fiber:age`             | Average age of fibers
 `fiber:length`          | Average length and variance of lengths of fibers
 `fiber:distribution`    | length distribution of fiber lengths (option: `max` and `interval`)
 `fiber:dynamic`         | Number of fiber classified by Dynamic state of plus end
 `fiber:point`           | coordinates of vertices of all fibers
 `fiber:displacement`    | mean squared displacement of fibers since the last call
 `fiber:moments`         | standard deviation of vertices of all fibers
 `fiber:speckle`         | coordinates of points randomly distributed along all fibers (option: `interval`)
 `fiber:sample`          | coordinates of points newly distributed along all fibers (option: `interval`)
 `fiber:segment`         | information about lengths of segments, number of kinks
 `fiber:end`             | Positions and dynamic states of all fiber ends
 `fiber:force`           | Position of vertices and Forces acting on vertices
 `fiber:tension`         | Internal stress along fibers
 `fiber:energy`          | Fiber's elastic bending energy
 `fiber:confinement`     | Force applied by fibers on their confinement Space
 `fiber:binder`          | Positions of bridging hands along each fiber
 `fiber:lattice`         | Total quantity on fiber's lattices
 `fiber:mesh`            | Total quantity on fiber's meshes
 `fiber:intersection`    | Intersections point of fibers
 `fiber:hand`            | Position of hands attached to fibers
 `fiber:link`            | Positions of attached hands for which stiffness > 0
 `fiber:cluster`         | Clusters made of fibers connected by Couples


 WHAT                    | Output
 ------------------------|------------------------------------------------------
 `bead:all`              | Position of beads
 `bead:single`           | Number of Beads with no single attached, 1 single attached etc.
 `solid:hand`            | Number of hands and number of attached hands on Solids
 `solid:orientation`     | Solid's position and orientation vectors
 `spindle:indice`        | Two scalar indices that caracterize the organization of fibers
 `spindle:profile`       | Number of right- and left-pointing fiber as a function of position
 `single:all`            | Position and force of singles
 `single:force`          | Average and maximum force of singles
 `NAME_OF_SINGLE`        | Position and force of singles of given class
 `NAME:position`         | Position and force of singles of class NAME
 `couple:state`          | Position and state of all couples
 `NAME_OF_COUPLE`        | Position and state of couples of given class
 `couple:link`           | detailed information on doubly-attached couples
 `couple:configuration`  | number of Couples in { X, P, A, V, T } states
 `couple:force`          | Average and maximum of tension in the couple links
 `couple:histogram`      | Histogram of tension in the couple links
 `couple:active`         | Position of active couples
 `couple:anatomy`        | Composition of couples
 `NAME:position`         | Position of couples of class NAME
 `couple:hands`          | Composition of couples
 
 */
void Simul::report_one(std::ostream& out, std::ostream& com, std::string const& who,
                       Property const* sel, std::string const& what, Glossary& opt) const
{
    if ( who == "fiber" )
    {
        if ( what.empty() || what == "position" )
            return reportFibers(out, com, sel);
        if ( what == "plus_end" )
            return reportFiberEnds(out, com, PLUS_END, sel);
        if ( what == "minus_end" )
            return reportFiberEnds(out, com, MINUS_END, sel);
        if ( what == "end" )
            return reportFiberEnds(out, com, BOTH_ENDS, sel);
        if ( what == "point" )
            return reportFiberPoints(out, com, sel);
        if ( what == "displacement" )
            return reportFiberDisplacement(out, com, sel);
        if ( what == "direction" )
            return reportFiberDirections(out, com, sel);
        if ( what == "moment" )
            return reportFiberMoments(out, com);
        if ( what == "speckle" )
            return reportFiberSpeckles(out, com, opt);
        if ( what == "sample" )
            return reportFiberSamples(out, com, opt);
        if ( what == "length" )
            return reportFiberLengths(out, com, sel);
        if ( what == "mark" )
            return reportFiberMarkedLengths(out, com, sel);
        if ( what == "distribution" || what == "histogram" )
            return reportFiberLengthHistogram(out, com, opt);
        if ( what == "tension" )
            return reportFiberTension(out, com, opt);
        if ( what == "segment" )
            return reportFiberSegments(out, com);
        if ( what == "energy" )
            return reportFiberBendingEnergy(out, com);
        if ( what == "extension" )
            return reportFiberExtension(out, com);
        if ( what == "nematic" )
            return reportFiberNematic(out, com, opt);
        if ( what == "end_state" || what == "dynamic" )
        {
            std::ofstream nul("/dev/null");
            reportFiberEndState(out, com, PLUS_END, sel);
            reportFiberEndState(out, nul, MINUS_END, sel);
            return;
        }
        if ( what == "plus_state" )
            return reportFiberEndState(out, com, PLUS_END, sel);
        if ( what == "minus_state" )
            return reportFiberEndState(out, com, MINUS_END, sel);
        if ( what == "force" )
            return reportFiberForces(out, com);
        if ( what == "confine_force" )
            return reportFiberConfineForce(out, com);
        if ( what == "confinement" )
            { reportFiberConfinement(out, com); return; }
        if ( what == "cluster" )
            return reportClusters(out, com, opt);
        if ( what == "age" )
            return reportFiberAge(out, com);
        if ( what == "intersection" )
            return reportFiberIntersections(out, com, opt);
        if ( what == "hand" )
            return reportFiberHands(out, com);
        if ( what == "link" )
            return reportFiberLinks(out, com);
        if ( what == "lattice" )
            return reportFiberLattice(out, com, sel);
        if ( what == "density" )
            return reportFiberDensityTotal(out, com, false, sel);
        if ( what == "density_avg" )
            return reportFiberDensityTotal(out, com, true, sel);
        if ( what == "density_all" )
            return reportFiberDensity(out, com, false, sel);
        if ( what == "connector" )
            return reportFiberConnectors(out, com, opt);
        if ( what == "collision" )
            return reportFiberCollision(out, com, sel, opt);

        throw InvalidSyntax("I can only report fiber: position, end, minus_end, plus_end, "\
                            "point, moment, speckle, sample, segment, dynamic, length, extension,"\
                            "distribution, tension, force, cluster, age, energy, hand, link\n");
    }
    if ( who == "bead" )
    {
        if ( what == "position" || what.empty() )
            return reportBeadPosition(out, com, sel);
        if ( what == "single" )
            return reportBeadSingles(out, com);
        throw InvalidSyntax("I can only report bead: position, single\n");
    }
    if ( who == "solid" )
    {
        if ( what == "hand" )
            return reportSolidHands(out, com, sel);
        if ( what == "orientation" )
            return reportSolidOrientation(out, com, sel);
        if ( what == "position" || what.empty() )
            return reportSolidPosition(out, com, sel);
        throw InvalidSyntax("I can only report solid: hand, position\n");
    }
    if ( who == "space" )
    {
        if ( what == "force" )
            return reportSpaceForce(out, com);
        if ( what.empty() )
            return reportSpace(out, com);
        throw InvalidSyntax("I can only report `space` and space:force\n");
    }
    if ( who == "sphere" )
    {
        if ( what == "position" || what.empty() )
            return reportSpherePosition(out, com, sel);
        throw InvalidSyntax("I can only report sphere:position\n");
    }
    if ( who == "single" )
    {
        if ( what.empty() )
            return reportSingle(out, com, sel);
        if ( what == "link" )
            return reportSingleLink(out, com, sel);
        if ( what == "state" )
            return reportSingleState(out, com, sel);
        if ( what == "force" )
            return reportSingleForce(out, com, sel);
        if ( what == "position" )
            return reportSinglePosition(out, com, sel);
        throw InvalidSyntax("I can only report single: link, state, force, position\n");
    }
    if ( who == "couple" )
    {
        if ( what.empty() )
            return reportCouple(out, com, sel);
        if ( what == "list" )
            return reportCoupleList(out, com, sel);
        if ( what == "state" )
            return reportCoupleState(out, com, sel);
        if ( what == "link" )
            return reportCoupleLink(out, com, sel);
        if ( what == "configuration" )
            return reportCoupleConfiguration(out, com, sel, opt);
        if ( what == "force" )
            return reportCoupleForce(out, com, sel);
        if ( what == "histogram" )
            return reportCoupleForceHistogram(out, com, opt);
        if ( what == "active" )
            return reportCoupleActive(out, com, sel);
        if ( what == "anatomy" )
            return reportCoupleAnatomy(out, com, sel);
        throw InvalidSyntax("I can only report couple: state, link, configuration, active, force, anatomy\n");
    }
    if ( who == "organizer" )
    {
        if ( what.empty() )
            return reportOrganizer(out, com);
        throw InvalidSyntax("I can only report `organizer'\n");
    }
    if ( who == "aster" )
    {
        if ( what.empty() )
            return reportAster(out, com);
        throw InvalidSyntax("I can only report `aster'\n");
    }
    if ( who == "field" )
    {
        return reportField(out, com);
    }
    if ( who == "simul" || who.empty() )
    {
        if ( what.empty() )
            return reportSimul(out, com);
        if ( what == "time" )
            return reportTime(out);
        if ( what == "inventory" )
            return reportInventory(out);
        if ( what == "property" || what == "parameter" )
            return writeProperties(out, false);
        throw InvalidSyntax("I can only report simul: time, inventory, property\n");
    }
    if ( who == "property" )
    {
        if ( what.empty() || what=="all" )
            return writeProperties(out, false);
        Property const* p = findProperty(what);
        if ( p )
            return p->write(out);
        throw InvalidSyntax("unknown property");
    }
    if ( who == "spindle" )
    {
        if ( what == "indice" )
            return reportSpindleIndices(out, com);
        if ( what == "length" )
            return reportSpindleLength(out, com, opt);
        if ( what == "minus_end" )
            return reportMarkedFiberEnds(out, com, opt);
        if ( what == "profile" )
            return reportSpindleProfile(out, com, opt);
        if ( what == "fitnes" )
            return reportSpindleFitness(out, com, opt);
        throw InvalidSyntax("I can only report spindle: indices, profile, fitness\n");
    }
    if ( who == "network" )
    {
        if ( what == "bridge" )
            return reportNetworkBridges(out, com, opt);
        if ( what == "size" )
            return reportNetworkSize(out, com);
    }
    if ( who == "ring" )
        return reportRing(out, com);
    if ( who == "platelet" )
        return reportPlatelet(out, com);
    if ( who == "ashbya" )
        return reportAshbya(out, com);
    if ( who == "something" )
        return reportSomething(out, com);

    if ( what.empty() )
        throw InvalidSyntax("unknown report `"+who+"'\n");
    else
        throw InvalidSyntax("unknown report `"+who+":"+what+"'\n");
}

//------------------------------------------------------------------------------
#pragma mark - Fiber Aggregated Properties

/**
 Export average length and variance of lengths for each class of fiber
 */
void Simul::reportFiberAge(std::ostream& out, std::ostream& com) const
{
#if !FIBER_HAS_BIRTHTIME
    out << SEP << "birthtime information disabled at compile time";
    return;
#endif
    com << COM << ljust("class", 2, 2) << SEP << "count" << SEP << "avg_birth";
    com << SEP << "var_birth" << SEP << "avg_age" << SEP << "min_age" << SEP << "max_age";
    
    const real now = time();

    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        ObjectList objs = fibers.collect(fp);
        size_t cnt = 0;
        real avg = 0, var = 0, mn = INFINITY, mx = -INFINITY;
        fibers.infoBirthtime(objs, cnt, avg, var, mn, mx);
        out << LIN << ljust(fp->name(), 2);
        out << SEP << cnt;
        out << SEP << avg;
        out << SEP << var;
        out << SEP << now-mx;
        out << SEP << now-avg;
        out << SEP << now-mn;
    }
}


/**
Export average length and variance of length for a class of fiber
*/
static void printFiberLengths(std::ostream& out, ObjectList const& objs)
{
    size_t cnt = 0;
    std::streamsize p = out.precision();
    real avg = 0, var = 0, mn = INFINITY, mx = -INFINITY, off = 0;
    FiberSet::infoLength(objs, cnt, avg, var, mn, mx, off);
    out << SEP << cnt;
    out.precision(3);
    out << SEP << std::fixed << avg;
    out << SEP << std::fixed << var;
    out << SEP << std::fixed << mn;
    out << SEP << std::fixed << mx;
    out.precision(1);
    out << SEP << std::fixed << avg*cnt;
    out << SEP << std::fixed << off;
    out.precision(p);
}

/**
 Export average length and variance for each class of fiber
 */
void Simul::reportFiberLengths(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << ljust("class", 2, 2) << SEP << "count" << SEP << "avg_len" << SEP << "var_len";
    com << SEP << "min_len" << SEP << "max_len" << SEP << "total" << SEP << "off";
    
    if ( sel )
    {
        ObjectList objs = fibers.collect(sel);
        out << LIN << ljust(sel->name(), 2);
        printFiberLengths(out, objs);
    }
    else
    {
        for ( Property const* i : properties.find_all("fiber") )
        {
            ObjectList objs = fibers.collect(i);
            out << LIN << ljust(i->name(), 2);
            printFiberLengths(out, objs);
        }
    }
}

/**
 Export average length and variance for fibers grouped by mark
 */
void Simul::reportFiberMarkedLengths(std::ostream& out, std::ostream& com, Property const*) const
{
    com << COM << ljust("mark", 1, 2) << SEP << "count" << SEP << "avg_len" << SEP << "var_len";
    com << SEP << "min_len" << SEP << "max_len" << SEP << "total" << SEP << "off";
    
    ObjectMark sup = 0;
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
        sup = std::max(sup, fib->mark());

    for ( ObjectMark k = 0; k <= sup; ++k )
    {
        uintptr_t val = k;
        ObjectList objs = fibers.collect(match_mark, reinterpret_cast<void*>(val));
        if ( objs.size() )
        {
            out << LIN << ljust(std::to_string(val), 1);
            printFiberLengths(out, objs);
        }
    }
}


/**
 Export length histograms for each class of fiber
 */
void Simul::reportFiberLengthHistogram(std::ostream& out, std::ostream& com, Glossary & opt) const
{
    const size_t BMAX = 256;
    unsigned cnt[BMAX+1];

    real sup = 0;
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
        sup = std::max(sup, fib->length());

    size_t nbin = BMAX;
    real delta = ( sup > 2 ) ? 1 : 0.1;
    opt.set(delta, "interval", "bin");
    if ( !opt.set(nbin, "interval", 1, "bin", 1) )
        nbin = std::ceil(sup/delta);
    nbin = std::min(nbin, BMAX);
    
    if ( 1 )
    {
        com << COM << "bin " << delta << " count " << fibers.size();
        out << LIN << ljust("scale", 2);
        std::streamsize p = out.precision(2);
        for ( size_t u = 0; u <= nbin; ++u )
            out << " " << std::setw(5) << delta * ( u + 0.5 );
        out.precision(p);
    }
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        
        for ( size_t u = 0; u <= nbin; ++u )
            cnt[u] = 0;
        
        for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
        {
            if ( fib->prop == fp )
            {
                size_t u = (size_t)std::floor( fib->length() / delta );
                ++cnt[std::min(u, nbin)];
            }
        }

        out << LIN << ljust(fp->name(), 2);
        for ( size_t u = 0; u <= nbin; ++u )
            out << " " << std::setw(5) << cnt[u];
    }
}


/**
 Export number of fiber, classified according to dynamic state of one end
 */
void Simul::reportFiberEndState(std::ostream& out, std::ostream& com, FiberEnd end, Property const* sel) const
{
    std::string name;
    if ( sel )
        name = sel->name() + ":";
    name.append(end==PLUS_END ?"plus":"minus");
    
    com << COM << ljust("class", 2, 2) << SEP << "total" << SEP << "static";
    com << SEP << "green" << SEP << "yellow" << SEP << "orange" << SEP << "red";
    
    constexpr size_t MAX = 5;
    size_t cnt[MAX+1] = { 0 };
    size_t sum = 0;
    
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( !sel || sel == fib->prop )
        {
            ++sum;
            state_t x = std::min(fib->endState(end), (state_t)MAX);
            ++cnt[x];
        }
    }
    
    out << LIN << ljust(name, 2) << SEP << sum;
    for ( size_t i = 0; i < MAX; ++i )
        out << SEP << cnt[i];
}


//------------------------------------------------------------------------------
#pragma mark - Fiber conformation


void Simul::reportFiberSegments(std::ostream& out, std::ostream& com) const
{
    com << COM << ljust("class", 2, 2) << SEP << "fibers" << SEP << "joints";
    com << SEP << "kinks" << SEP << "min_seg" << SEP << "max_seg" << SEP << "err_seg";
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        
        size_t cnt, points;
        real mn = 0, mx = 0, dv = 0;
        
        ObjectList objs = fibers.collect(fp);
        fibers.infoSegments(objs, cnt, points, mn, mx, dv);
        out << LIN << ljust(fp->name(), 2);
        out << SEP << cnt;
        out << SEP << points - 2 * cnt;
        out << SEP << fibers.nbKinks(objs);
        out << SEP << std::fixed << mn;
        out << SEP << std::fixed << mx;
        out << SEP << std::fixed << dv;
    }
}


/**
 Export fiber elastic bending energy
 */
void Simul::reportFiberBendingEnergy(std::ostream& out, std::ostream& com) const
{
    com << COM << ljust("bending_energy",2,2) << SEP << "count";
    com << SEP << "sum" << SEP << "avg" << SEP << "var" << SEP << "rigidity";
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        ObjectList objs = fibers.collect(fp);
        size_t cnt = 0;
        real avg = 0, var = 0;
        fibers.infoBendingEnergy(objs, cnt, avg, var);
        if ( cnt > 0 )
        {
            out << LIN << ljust(fp->name(), 2);
            out << SEP << cnt;
            out << SEP << avg*cnt;
            out << SEP << avg;
            out << SEP << var;
            out << SEP << fp->rigidity;
        }
    }
}


void Simul::reportFiberExtension(std::ostream& out, std::ostream& com) const
{
    com << COM << ljust("end_to_end_dist",2,2) << SEP << "count";
    com << SEP << "avg" << SEP << "var" << SEP << "min" << SEP << "max";
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        size_t cnt = 0;
        real avg = 0, var = 0, in = INFINITY, ax = -INFINITY;
        FiberProp const* p = static_cast<FiberProp const*>(i);
        for ( Object const* f : fibers.collect(p) )
        {
            Fiber const* fib = Fiber::toFiber(f);
            real x = ( fib->posEndP() - fib->posEndM() ).norm();
            in = std::min(in, x);
            ax = std::max(ax, x);
            avg += x;
            var += x * x;
            ++cnt;
        }
        if ( cnt > 0 )
        {
            avg /= cnt;
            var -= square(avg)*cnt;
            if ( cnt > 1 )
                var /= real(cnt-1);
            out << LIN << ljust(p->name(), 2);
            out << SEP << cnt;
            out << SEP << avg;
            out << SEP << var;
            out << SEP << in;
            out << SEP << ax;
        }
    }
}


void Simul::reportFiberNematic(std::ostream& out, std::ostream& com, FiberProp const* fp, Space const* spc) const
{
    com << COM << ljust("class",2,2) << SEP << "time" << SEP << "count";
    com << SEP << "order" << SEP << "dirX" << SEP << "dirY" << SEP << "dirZ";
    if ( spc )
        com << SEP << "ortho" << SEP << "dirX" << SEP << "dirY" << SEP << "dirZ";
    
    ObjectList objs = fibers.collect(fp);
    if ( objs.size() > 0 )
    {
        real S = 0;
        real vec[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
        S = FiberSet::infoNematic(objs, vec);
        out << LIN << ljust(fp->name(), 2);
        out << SEP << time();
        out << SEP << objs.size();
        out << SEP << S;
        out << SEP << vec[0] << SEP << vec[1] << SEP << vec[2];
        if ( spc )
        {
            S = FiberSet::infoOrthoNematic(objs, vec, spc);
            out << SEP << S;
            out << SEP << vec[0] << SEP << vec[1] << SEP << vec[2];
        }
    }
}


void Simul::reportFiberNematic(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    std::string str;
    Space const* spc = nullptr;
    if ( opt.set(str, "space") )
        spc = findSpace(opt.value("space"));

    for ( Property const* i : properties.find_all("fiber") )
    {
        reportFiberNematic(out, com, static_cast<FiberProp const*>(i), spc);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Fiber Lattice & Density


/**
 Report quantity of substance in the fiber's Lattice
 */
void Simul::reportFiberLattice(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << ljust("class", 2, 2) << SEP << "count" << SEP << "vacant";
    com << SEP << "avg" << SEP << "avg_if" << SEP << "min" << SEP << "max";
    
    size_t cnt = 0, vac = 0;
    real sum = 0, mn = INFINITY, mx = -INFINITY;
    
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( !sel || sel == fib->prop )
            fib->infoLattice(cnt, vac, sum, mn, mx);
    }
    
    std::streamsize p = out.precision(4);
    out << LIN << ljust("fiber:lattice", 2);
    out << SEP << cnt;
    out << SEP << vac;
    out << SEP << sum / (real)cnt;
    out << SEP << sum / (real)(cnt-vac);
    out.precision(2);
    out << SEP << std::fixed << mn;
    out << SEP << std::fixed << mx;
    out.precision(p);
}


/**
 Report quantity of substance in the fiber's Density
 */
void Simul::reportFiberDensityTotal(std::ostream& out, std::ostream& com, bool density, Property const* sel) const
{
    com << COM << ljust("class", 2, 2) << SEP << "total";
    com << SEP << "avg" << SEP << "min" << SEP << "max" << SEP << "length";
    
    size_t cnt = 0, vac = 0;
    real len = 0, sum = 0, mn = INFINITY, mx = -INFINITY;
    
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( !sel || sel == fib->prop )
            fib->infoDensity(len, cnt, vac, sum, mn, mx, density);
    }
    
    std::streamsize p = out.precision(4);
    out << LIN << ljust("fiber:mesh", 2);
    out << SEP << sum;
    out << SEP << sum / (real)cnt;
    out.precision(6);
    out << SEP << std::fixed << mn;
    out << SEP << std::fixed << mx;
    out << SEP << std::setprecision(3) << len;
    out.precision(p);
}

/**
 Report quantity of substance in the fiber's Density
 */
void Simul::reportFiberDensity(std::ostream& out, std::ostream& com, bool density, Property const* sel) const
{
    com << COM << ljust("fiber", 2, 2) << SEP << "total";
    com << SEP << "avg" << SEP << "min" << SEP << "max" << SEP << "length";
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        size_t cnt = 0, vac = 0;
        real len = 0, sum = 0, mn = INFINITY, mx = -INFINITY;
        if ( !sel || sel == fib->prop )
        {
            fib->infoDensity(len, cnt, vac, sum, mn, mx, density);
            std::streamsize p = out.precision(4);
            out << LIN << ljust(fib->reference(), 2);
            out << SEP << sum;
            out << SEP << sum / (real)cnt;
            out << SEP << std::fixed << std::setprecision(6) << mn;
            out << SEP << std::fixed << std::setprecision(6) << mx;
            out << SEP << std::setprecision(3) << len;
            out.precision(p);
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Fiber Individual Properties

/**
 Export length, position and directions at center of fibers
 */
void Simul::reportFiber(std::ostream& out, Fiber const* fib) const
{
    out << LIN << fib->prop->number();
    out << SEP << fib->identity();
    out << SEP << fib->length();
    out << SEP << fib->posEnd(CENTER);
    out << SEP << fib->dirEnd(CENTER);
    out << SEP << (fib->posEndM()-fib->posEndP()).norm();
    real C = dot(fib->dirEndM(), fib->dirEndP());
#if ( DIM >= 3 )
    real S = cross(fib->dirEndM(), fib->dirEndP()).norm();
#else
    real S = cross(fib->dirEndM(), fib->dirEndP());
#endif
    out << SEP << C << SEP << S;
    out << SEP << organizers.findOrganizerID(fib);
}

    
/// qsort function comparing length of Fibers
static int compareFiberLength(Object const* A, Object const* B)
{
    real a = static_cast<Fiber const*>(A)->length();
    real b = static_cast<Fiber const*>(B)->length();
    return ( a < b ) - ( a > b ); //sort in decreasing length
    //return ( a > b ) - ( a < b );
}


/**
 This will sort fibers by decreasing length
 Export length, position and directions at center of fibers
 */
void Simul::reportFibersSorted(std::ostream& out, std::ostream& com, Property const* sel)
{
    fibers.pool_.blinksort(compareFiberLength);
    
    com << COM << "class" << SEP << "identity" << SEP << "length";
    com << SEP << repeatXYZ("pos") << SEP << repeatXYZ("dir") << SEP << "endToEnd";
    com << SEP << "cos" << SEP << "sin" << SEP << "organizer";

    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( !sel || sel == fib->prop )
            reportFiber(out, fib);
    }
}


/**
 Export length, position and directions at center of fibers
 */
void Simul::reportFibers(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity" << SEP << "length";
    com << SEP << repeatXYZ("pos") << SEP << repeatXYZ("dir") << SEP << "endToEnd";
    com << SEP << "cos" << SEP << "sin" << SEP << "organizer";
    
    // list fibers in the order of the inventory:
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( !sel || sel == fib->prop )
            reportFiber(out, fib);
    }
}


/**
 Export dynamic state, positions and directions of fiber
 Argument `end` can be { MINUS_END, PLUS_END, BOTH_ENDS }
 */
void Simul::reportFiberEnds(std::ostream& out, std::ostream& com, FiberEnd end, Property const* sel) const
{
    com << COM << "class" << SEP << "identity" << SEP << "length";
    if ( end & PLUS_END )
        com << SEP << "stateP" << SEP << repeatXYZ("posP") << SEP << repeatXYZ("dirP");
    if ( end & MINUS_END )
        com << SEP << "stateM" << SEP << repeatXYZ("posM") << SEP << repeatXYZ("dirM");
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( sel && sel != fib->prop )
            continue;
        out << LIN << fib->prop->number();
        out << SEP << fib->identity();
        out << SEP << fib->length();
        if ( end & PLUS_END )
        {
            out << SEP << fib->endStateP();
            out << SEP << fib->posEndP();
            out << SEP << fib->dirEndP();
        }
        if ( end & MINUS_END )
        {
            out << SEP << fib->endStateM();
            out << SEP << fib->posEndM();
            out << SEP << fib->dirEndM();
        }
    }
}


/**
 Export Fiber-number, position of vertices
 */
void Simul::reportFiberPoints(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "identity" << SEP << repeatXYZ("pos") << SEP << "curvature";

    // list fibers in the order of the inventory:
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( sel && sel != fib->prop )
            continue;
        com << COM << "fiber " << fib->reference() << "  " << fib->segmentation();
        
        for ( index_t p = 0; p < fib->nbPoints(); ++p )
        {
            out << LIN << fib->identity();
            out << SEP << fib->posPoint(p);
            out << SEP << fib->curvature(p);
        }
    }
}


/**
 Export positions of points taken randomly along all fibers,
 but that remain static with respect to the lattice of each fiber,
 during the life-time of this fiber.
 
 This is meant to simulate the `speckle microscopy` that is obtained
 in microcscopy with a low amount of fluorescent-monomers.
 
 The distance between the speckles follows an exponential distribution
 with an average defined by a parameter `gap` or `1/density` defined in `opt`.
 */
void Simul::reportFiberSpeckles(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    real gap = 1;
    if ( opt.set(gap, "density") )
        gap = 1.0 / gap;
    else
        opt.set(gap, "interval", "gap");
    constexpr real TINY = 0x1p-32;

    Fiber const* fib = fibers.first();
    while ( fib )
    {
        com << COM << "fiber " << fib->reference();
        
        // generate speckles below the origin of abscissa
        if ( fib->abscissaM() < 0 )
        {
            uint64_t Z = pcg32_init(fib->signature());
            real a = gap * std::log(pcg32(Z)*TINY);
            while ( a > fib->abscissaP() )
            {
                a += gap * std::log(pcg32(Z)*TINY);
            }
            while ( a >= fib->abscissaM() )
            {
                out << '\n' << fib->pos(a) << " " << a;
                a += gap * std::log(pcg32(Z)*TINY);
            }
        }
        // generate speckles above the origin of abscissa
        if ( fib->abscissaP() > 0 )
        {
            uint64_t Z = pcg32_init(~fib->signature());
            real a = -gap * std::log(pcg32(Z)*TINY);
            while ( a < fib->abscissaM() )
                a -= gap * std::log(pcg32(Z)*TINY);
            while ( a <= fib->abscissaP() )
            {
                out << '\n' << fib->pos(a) << " " << a;
                a -= gap * std::log(pcg32(Z)*TINY);
            }
        }
        
        fib = fib->next();
    }
}


/**
 Export positions of points taken randomly along all fibers,
 changing the distribution at every time.
 */
void Simul::reportFiberSamples(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    real gap = 1;
    if ( opt.set(gap, "density") )
        gap = 1.0 / gap;
    else
        opt.set(gap, "interval", "gap");
    
    FiberSiteList loc(1024);
    fibers.uniFiberSites(loc, gap);
    
    Fiber const* ofib = nullptr;
    for ( FiberSite & i : loc )
    {
        if ( ofib != i.fiber() )
        {
            com << COM << "fiber " << i.fiber()->reference();
            ofib = i.fiber();
        }
        
        out << LIN << i.pos();
    }
}


/**
 Export Mean Squared Displacement of fiber's minus ends since the last call
 to this function.
 \todo: it would be simpler to add a 'Vector oldPos;' directly in class Fiber.
 */
void Simul::reportFiberDisplacement(std::ostream& out, std::ostream& com, Property const* sel) const
{
    typedef std::map <ObjectID, Vector> fiber_map;
    static fiber_map positions;
    static double past = 0;
    
    com << COM << "delta_time nb_fibers mean_squared_displacement";
    
    real sum = 0;
    size_t cnt = 0;
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( sel && sel != fib->prop )
            continue;
        /* using minus end to avoid effects of plus-end growth,
         but there is no check that the minus-end is not growing... */
        Vector pos = fib->posEndM();
        fiber_map::iterator i = positions.find(fib->identity());
        if ( i != positions.end() )
        {
            ++cnt;
            sum += distanceSqr(pos, i->second);
            i->second = pos;
        }
        else
        {
            positions[fib->identity()] = pos;
        }
    }
    
    if ( cnt > 0 )
        out << LIN << time() - past << SEP << cnt << SEP << sum / (real)cnt;
    else
        out << LIN << time() - past << SEP << 0 << SEP << 0;
    
    past = time();
}


/** This is a hack to maintain some credibility */
static bool isSymmetricAroundAxisZ(std::string const& shape)
{
    if ( shape == "sphere" ) return true;
    if ( shape == "cylinderZ" ) return true;
    if ( shape == "ellipse" ) return true;
    if ( shape == "torus" ) return true;
    if ( shape == "ring" ) return true;
    if ( shape == "disc" ) return true;
    return false;
}

void Simul::reportFiberDirections(std::ostream& out, std::ostream& com, Property const* sel) const
{
    Space const* spc = spaces.master();
    if ( sel )
        spc = static_cast<FiberProp const*>(sel)->confine_space;
    if ( !isSymmetricAroundAxisZ(spc->prop->shape) )
        throw InvalidParameter("reportFiberDirections() cannot handle non symmetric Space");

    real sum = 0;
    Vector eZ(0, 0, 1);
    Vector2 avg(0, 0);
    Matrix22 mat(0, 0);
#if ( DIM == 3 )
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( sel && sel != fib->prop )
            continue;
        for ( index_t p = 0; p < fib->nbSegments(); ++p )
        {
            Vector pos = fib->midPoint(p);
            Vector dir = fib->dirSegment(p);
            Vector nor = spc->normalToEdge(pos);
            Vector tan = cross(eZ, nor);
            real n = normSqr(tan);
            if ( n > 0.5 )
            {
                Vector2 V(dot(dir, tan/sqrt(n)), dir.ZZ);
                V.normalize();
                mat(0,0) += V.XX * V.XX;
                mat(1,0) += V.XX * V.YY;
                mat(1,1) += V.YY * V.YY;
                avg += V;
                ++sum;
            }
        }
    }
#endif
    real X = 0, Y = 0, S = 0;
    if ( sum > 0 )
    {
        avg /= sum;
        mat *= 2.0 / sum;
        mat(0,1) = mat(1,0);
        X = std::sqrt(mat(0,0) * 0.5);
        Y = std::sqrt(mat(1,1) * 0.5);
        // subtract trace:
        mat(0,0) -= 1.0;
        mat(1,1) -= 1.0;
        //std::clog << mat << " ";
        // eigenvalue of a 2x2 traceless symmetric matrix:
        S = sqrt(square(mat(0,0))+square(mat(1,0)));
    }
    // polar order parameter:
    real M = norm(avg);
    com << COM << "n_seg" << SEP << "nematic" << SEP << "polar" << SEP << "avg_t" << SEP << "avg_z" << SEP << "mom_t" << SEP << "mom_z";
    out << LIN << sum << SEP << S << SEP << M << SEP << avg.XX << SEP << avg.YY << SEP << X << SEP << Y;
}


/**
 Export first and second-order moments of vertices for each class of Fiber
 */
void Simul::reportFiberMoments(std::ostream& out, std::ostream& com) const
{
    com << COM << ljust("class", 2, 2) << SEP << "sum";
    com << SEP << "avgX" << SEP << "avgY" << SEP << "avgZ";
    com << SEP << "varX" << SEP << "varY" << SEP << "varZ" << SEP << "var_sum";
    out << std::fixed;
    
    Accumulator acc;
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        
        acc.reset();
        
        for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
        {
            if ( fib->prop == fp )
            {
                const real w = fib->segmentation();
                acc.add(0.5*w, fib->posEndM());
                for ( index_t n = 1; n < fib->lastPoint(); ++n )
                    acc.add(w, fib->posPoint(n));
                acc.add(0.5*w, fib->posEndP());
            }
        }
        
        acc.subtract_mean();
        out << LIN << ljust(fp->name(), 2);
        acc.print(out, 0);
    }
}

/**
 Export first and second-order moments of vertices for each class of Fiber
 */
void Simul::reportNetworkSize(std::ostream& out, std::ostream& com) const
{
    com << COM << "polymer" << SEP << "surface";
    out << std::fixed;
    Accumulator acc;
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        const real w = fib->segmentation();
        acc.add(0.5*w, fib->posEndM());
        for ( index_t n = 1; n < fib->lastPoint(); ++n )
            acc.add(w, fib->posPoint(n));
        acc.add(0.5*w, fib->posEndP());
    }
    acc.subtract_mean();
    real S = 2 * M_PI * acc.total_variance();
    out << LIN << acc.total_length() << SEP << S;
}

//------------------------------------------------------------------------------
#pragma mark - Fiber forces

/**
 Export Fiber-number, position of vertices and tension in each segment
 */
void Simul::reportFiberForces(std::ostream& out, std::ostream& com) const
{
    computeForces();

    com << COM << "identity" << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force") << SEP << "tension";
    
    // list fibers in the order of the inventory:
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        com << COM << "fiber " << fib->reference();
        
        for ( index_t p = 0; p < fib->nbPoints(); ++p )
        {
            out << LIN << fib->identity();
            out << SEP << fib->posPoint(p);
            out << SEP << fib->netForce(p);
            if ( p == fib->lastPoint() )
                out << SEP << 0.0;
            else
                out << SEP << fib->tension(p);
        }
    }
}


/**
 Sum of the internal tensions from fiber segments that intersect a plane
 specified in `opt`.
 The plane is defined by <em> n.pos + a = 0 </em>

     plane = NORMAL, SCALAR

 */
void Simul::reportFiberTension(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    computeForces();
    
    com << COM << "count" << SEP << "sum_force" << SEP << "min_force" << SEP << "max_force";

    Vector n(1,0,0);
    real ten = 0, inf = 0, sup = 0;
    size_t cnt = 0;
    if ( opt.value_is("plane", 0, "all") )
    {
        // extending the comments:
        for ( int d = 1; d < DIM; ++d )
            out << SEP << "count" << SEP << "force";
        
        // plane orthogonal to X:
        fibers.infoTension(cnt, ten, inf, sup, Vector(1,0,0), 0);
        out << LIN << cnt << SEP << ten << SEP << inf << SEP << sup;
#if ( DIM > 1 )
        // plane orthogonal to Y:
        fibers.infoTension(cnt, ten, inf, sup, Vector(0,1,0), 0);
        out << SEP << cnt << SEP << ten << SEP << inf << SEP << sup;
#endif
#if ( DIM > 2 )
        // plane orthogonal to Z:
        fibers.infoTension(cnt, ten, inf, sup, Vector(0,0,1), 0);
        out << SEP << cnt << SEP << ten << SEP << inf << SEP << sup;
#endif
    }
    else if ( opt.set(n, "plane") )
    {
        real a = 0;
        opt.set(a, "plane", 1);
        com << COM << "fiber tension orthogonal to plane: (" << n << ").pos = " << -a;
        fibers.infoTension(cnt, ten, inf, sup, n, a);
        out << LIN << cnt << SEP << ten << SEP << inf << SEP << sup;
    }
    else
    {
        // if no plane is specified, sum all tension from all segments
        fibers.infoTension(cnt, ten, inf, sup);
        out << LIN << cnt << SEP << ten << SEP << inf << SEP << sup;
    }
}


/** Attention: This does not handle all cases */
static bool confinementApplies(Confinement mode, Space const* spc, Vector const& pos)
{
    switch ( mode )
    {
        case CONFINE_OFF: return false;
        case CONFINE_INSIDE: return spc->outside(pos);
        case CONFINE_OUTSIDE: return spc->inside(pos);
        case CONFINE_ON: return true;
        case CONFINE_PLUS_END: return true;
        case CONFINE_MINUS_END: return true;
        case CONFINE_BOTH_ENDS: return true;
        case CONFINE_MINUS_OUT: return spc->inside(pos);
        case CONFINE_PLUS_OUT: return spc->inside(pos);
        default: return true;
    }
    return true;
}

/** Attention: This does not handle all cases */
static bool vertexIsConfined(Confinement mode, Fiber const* fib, size_t inx)
{
    switch ( mode )
    {
        case CONFINE_OFF: return false;
        case CONFINE_INSIDE: return true;
        case CONFINE_OUTSIDE: return true;
        case CONFINE_ON: return true;
        case CONFINE_PLUS_END: return ( inx == fib->lastPoint() );
        case CONFINE_MINUS_END: return ( inx == 0 );
        case CONFINE_BOTH_ENDS: return ( inx == 0 || inx == fib->lastPoint() );
        case CONFINE_MINUS_OUT: return ( inx == 0 );
        case CONFINE_PLUS_OUT: return ( inx == fib->lastPoint() );
        default: return true;
    }
    return true;
}


/**
 Export total magnitude of force exerted by Fiber on the confinement
 */
void Simul::reportFiberConfineForce(std::ostream& out, std::ostream& com) const
{
    com << COM << "confinement forces";
    com << COM << "identity" << SEP << "vertex" << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force");
     
     // list fibers in the order of the inventory:
     for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
     {
         com << COM << "fiber " << fib->reference();
         Space const* spc = findSpace(fib->prop->confine_spec);
         const real stiff = fib->prop->confine_stiff[0];
         const Confinement mode = fib->prop->confine;

         for ( index_t p = 0; p < fib->nbPoints(); ++p )
         {
             Vector w, pos = fib->posPoint(p);
             if ( vertexIsConfined(mode, fib, p) && confinementApplies(mode, spc, pos) )
             {
                 w = spc->project(pos);
                 out << LIN << fib->identity();
                 out << SEP << p << SEP << pos;
                 out << SEP << stiff * ( w - pos );
             }
         }
     }
}

/**
 Export total magnitude of force exerted by Fiber on the confinement.
 The radial components of the forces are summed up, which is only meaningful
 in very particular systems, when the geometry is circular around the Z-axis!
 */
real Simul::reportFiberConfinement(std::ostream& out, std::ostream& com) const
{
    com << COM << "count" << SEP << repeatXYZ("force") << SEP << "radial";
    size_t cnt = 0;
    Vector sum(0,0,0);
    real   rad = 0;
    
#if ( DIM > 1 )
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        Space const* spc = findSpace(fib->prop->confine_spec);
        const real stiff = fib->prop->confine_stiff[0];
        const Confinement mode = fib->prop->confine;

        if ( !isSymmetricAroundAxisZ(spc->prop->shape) )
            throw InvalidParameter("reportFiberConfinement() cannot handle non symmetric Space");
        
        for ( index_t p = 0; p < fib->nbPoints(); ++p )
        {
            Vector w, pos = fib->posPoint(p);
            if ( vertexIsConfined(mode, fib, p) && confinementApplies(mode, spc, pos) )
            {
                ++cnt;
                w = spc->project(pos);
                // assuming the goemetry is rotational-symmetric around the Z axis:
                Vector dir = normalize(Vector(pos.XX, pos.YY, 0));
                Vector vec = stiff * ( pos - w );
                sum += vec;
                rad += dot(vec, dir);
            }
        }
    }
#endif
    out << LIN << cnt << SEP << sum << SEP << rad;
    return rad;
}


//------------------------------------------------------------------------------
#pragma mark - bound Hands per Fiber


void Simul::reportFiberHands(std::ostream& out, std::ostream& com) const
{
    com << COM << "fib_type" << SEP << "fib_id" << SEP << "class" << SEP << "abs";
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( fib->nbAttachedHands() > 0 )
        {
            com << COM << "on fiber " << fib->reference();
            fib->sortHands();
            for ( Hand const* h = fib->firstHand(); h; h = h->next() )
            {
                out << LIN << fib->prop->number();
                out << SEP << fib->identity();
                out << SEP << h->property()->number();
                out << SEP << h->abscissa();
            }
        }
    }
}


void Simul::reportFiberLinks(std::ostream& out, std::ostream& com) const
{
    com << COM << "fib_type" << SEP << "fib_id" << SEP << "class" << SEP << "abs" << SEP << "position";
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( fib->nbAttachedHands() > 0 )
        {
            com << COM << "on fiber " << fib->reference();
            fib->sortHands();
            for ( Hand const* h = fib->firstHand(); h; h = h->next() )
            {
                if ( h->linkStiffness() > 0 )
                {
                    out << LIN << fib->prop->number();
                    out << SEP << fib->identity();
                    out << SEP << h->property()->number();
                    out << SEP << h->abscissa();
                    out << SEP << h->linkFoot();
                    out << SEP << h->linkStiffness();
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Networks


void Simul::reportFiberIntersections(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    int details = 2;
    real up = 0;
    opt.set(up, "threshold");
    opt.set(details, "details");
    
    const real sup = up * up;
    
    if ( details == 2 )
    {
        com << COM << "id1" << SEP << "abs1";
        com << SEP << "id2" << SEP << "abs2" << SEP << repeatXYZ("pos");
    }
    Accumulator acc;
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        unsigned cnt = 0;
        for ( index_t ii = 0; ii < fib->nbSegments(); ++ii )
        {
            FiberSegment seg(fib, ii);
            for ( Fiber const* fox = fibers.nextID(fib); fox; fox = fibers.nextID(fox) )
            {
                for ( index_t jj = 0; jj < fox->nbSegments(); ++jj )
                {
                    FiberSegment soc(fox, jj);
                    real abs1, abs2;
                    real dis2 = seg.shortestDistanceSqr(soc, abs1, abs2);
                    if (( dis2 < sup ) && seg.within(abs1) && soc.within(abs2))
                    {
                        ++cnt;
                        Vector pos1 = seg.midPoint(abs1/seg.len());
                        //Vector pos2 = loc2.pos(abs2/loc2.len());
                        if ( details == 2 )
                        {
                            out << LIN << fib->identity();
                            out << SEP << abs1 + seg.abscissa1();
                            out << SEP << fox->identity();
                            out << SEP << abs2 + soc.abscissa1();
                            out << SEP << pos1;
                        }
                        acc.add(pos1);
                    }
                }
            }
        }
        if ( cnt && details >= 1 )
        {
            com << COM << "total" << SEP << fib->identity() << SEP << cnt;
        }
    }
    acc.subtract_mean();
    acc.print_doc(out);
    acc.print(out, 1);
}


/// accessory class to analyse the connections in a network of fibers
struct Connector
{
    real a;
    real s;
    long f;
    long g;
    long h;
    Connector(real as, long fs) { a = as; f = fs; g = -1; h = -1; s = 0; }
    Connector(real as, long fs, long gs) { a = as; f = fs; g = gs; h = -1; s = 0; }
    Connector(real as, long fs, long gs, long hs) { a = as; f = fs; g = gs; h = hs; s = 0; }
    Connector(real as, long fs, long gs, long hs, real ss) { a = as; f = fs; g = gs; h = hs; s = ss; }
    bool operator < (const Connector& r) const { return a < r.a; }
};


/**
 This is the older version
 */
void Simul::reportFiberConnectors(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    int details = 2;
    opt.set(details, "details");

    if ( details > 1 )
    {
        com << COM << "class  identity      abs1    fiber1     hand1      dist";
        com << "      abs2    fiber2     hand2      dist ...";
    }
    else
    {
        com << COM << "fiber connectors";
    }
    
    // used to calculate the size of the network from the position of connectors
    Accumulator acc;
    typedef std::vector<Connector> clist_t;
    typedef std::map<ObjectID, clist_t> map_t;
    
    map_t map;

    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        map.clear();
        // check all connecting Hands and record abscissa, depending on the fiber that is linked
        for ( Hand const* h = fib->firstHand(); h; h = h->next() )
        {
            Hand const* g = h->otherHand();
            if ( g && g->attached() )
            {
                ObjectID f2 = g->fiber()->identity();
                map[f2].push_back(Connector(h->abscissa(), h->property()->number()));
            }
        }
        if ( map.size() )
        {
            if ( details > 1 )
            {
                out << LIN << fib->prop->number() << SEP << fib->identity();
            }
            
            clist_t list;
            // average all the abscissa linking to the same fiber:
            for ( map_t::const_iterator mi = map.begin(); mi != map.end(); ++mi )
            {
                clist_t const& sec = mi->second;
                real a = 0.0;
                // number of connector of each type
                int c1 = 0, c2 = 0;
                for ( clist_t::const_iterator ci = sec.begin(); ci != sec.end(); ++ci )
                {
                    a += ci->a;
                    if ( ci->f == 1 ) ++c1;
                    if ( ci->f == 2 ) ++c2;
                 }
                a /= sec.size();
                list.push_back(Connector(a, mi->first, c1, c2, sec.size()-c1-c2));
            }
            // sort the list in increasing abscissa
            std::sort(list.begin(), list.end());
            
            clist_t::const_iterator p = list.begin();
            for ( clist_t::const_iterator c = list.begin(); c != list.end(); ++c )
            {
                if ( details > 1 )
                {
                    out << SEP << c->a;
                    out << SEP << c->f;
                    out << SEP << std::to_string(c->g)+"+"+std::to_string(c->h);
                    // calculate direct distance to previous point of intersection:
                    if ( c != p )
                        out << SEP << ( fib->pos(p->a) - fib->pos(c->a) ).norm();
                    else
                        out << SEP << 0;
                }
                p = c;
                acc.add(fib->pos(p->a));
            }
            if ( details > 0 )
            {
                com << COM << "total";
                com << SEP << fib->prop->number();
                com << SEP << fib->identity();
                com << SEP << list.size();
            }
        }
    }
    acc.subtract_mean();
    acc.print_doc(out);
    acc.print(out, 1);
}

/**
 F. Nedelec, 18/08/2017
 */
void Simul::reportNetworkBridges(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    int details = 0;
    opt.set(details, "details");

    com << COM << "length" << SEP << "speed" << SEP << "type1" << SEP << "type2";

    typedef std::vector<Connector> clist_t;
    typedef std::map<ObjectID, clist_t> map_t;
    
    map_t map;
    
    HandProp const* hp1 = findProperty<HandProp>("hand", 1);
    HandProp const* hp2 = findProperty<HandProp>("hand", 2);
    
    const real speedh1 = hp1->motorSpeed();
    const real speedh2 = hp2->motorSpeed();

    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        map.clear();
        // check all connecting Hands and record abscissa, depending on the fiber that is linked
        for ( Hand const* h = fib->firstHand(); h; h = h->next() )
        {
            Hand const* g = h->otherHand();
            if ( g && g->attached() )
            {
                ObjectID f2 = g->fiber()->identity();
                if ( h->property() == hp1 )
                    map[f2].push_back(Connector(h->abscissa(), 1));
                else if ( h->property() == hp2 )
                    map[f2].push_back(Connector(h->abscissa(), 2));
                else
                    com << COM << "report network:bridge can only handle 2 hand types";
            }
        }
        if ( map.size() )
        {
            clist_t list;
            // average all the abscissa linking to the same fiber:
            for ( map_t::const_iterator mi = map.begin(); mi != map.end(); ++mi )
            {
                clist_t const& sec = mi->second;
                // average abscissa:
                real a = 0.0;
                // number of connector of each type
                int c1 = 0, c2 = 0;
                for ( clist_t::const_iterator ci = sec.begin(); ci != sec.end(); ++ci )
                {
                    a += ci->a;
                    if ( ci->f == 1 ) ++c1;
                    if ( ci->f == 2 ) ++c2;
                }
                a /= sec.size();
                real speed = 0;
                if ( c1 > 0 && c2 > 0 )
                    speed = std::min(speedh1, speedh2);
                else if ( c1 > 0 )
                    speed = speedh1;
                else if ( c2 > 0 )
                    speed = speedh2;
                list.push_back(Connector(a, mi->first, c1, c2, speed));
            }
            // sort the list in increasing abscissa
            std::sort(list.begin(), list.end());
            
            if ( details > 0 )
            {
                com << COM << "connectors on fiber f" << fib->identity() << ":";
                // print all connector attachment positions:
                com << COM << "abscissa" << SEP << "fiber_id" << SEP << "speed" << SEP << "type";
                for ( clist_t::const_iterator c = list.begin(); c != list.end(); ++c )
                {
                    out << LIN << c->a;
                    out << SEP << c->f;
                    out << SEP << c->s;
                    out << SEP << std::to_string(c->g)+"+"+std::to_string(c->h);
                }
            }
            if ( list.size() > 1 )
            {
                // print all bridges
                if ( details > 0 )
                    com << COM << "length" << SEP << "speed" << SEP << "type1" << SEP << "type2";
#if ( 1 )
                for ( clist_t::const_iterator p = list.begin(); p != list.end(); ++p )
                for ( clist_t::const_iterator c = p+1; c != list.end(); ++c )
#else
                for ( clist_t::const_iterator p = list.begin(), c = p+1; c != list.end(); ++p, ++c )
#endif
                {
                    out << LIN << c->a - p->a;
                    out << SEP << c->s - p->s;
                    out << SEP << std::to_string(p->g)+"+"+std::to_string(p->h);
                    out << SEP << std::to_string(c->g)+"+"+std::to_string(c->h);
                }
            }
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Beads, Solid, Space


void Simul::reportTime(std::ostream& out) const
{
    out << LIN << time();
}


void Simul::reportInventory(std::ostream& out) const
{
    //com << COM << "properties:";
    //properties.write_names(out, "");
    //com << COM << "objects:";
    spaces.report(out);
    fields.report(out);
    fibers.report(out);
    spheres.report(out);
    beads.report(out);
    solids.report(out);
    singles.report(out);
    couples.report(out);
    organizers.report(out);
    events.report(out);
}

template < typename SET >
static void reportSet(std::ostream& out, SET& set, PropertyList const& properties)
{
    for ( Property const* i : properties.find_all(set.title()) )
    {
        ObjectID id = 0;
        index_t points = 0, sup = 0;
        ObjectList objs = set.collect(i);
        for ( Object * o : objs )
        {
            Mecable * mec = Simul::toMecable(o);
            if ( mec )
            {
                points += mec->nbPoints();
                sup = std::max(sup, mec->nbPoints());
                id = std::max(id, mec->identity());
            }
        }
        if ( points > 0 )
        {
            out << LIN << ljust(i->name(), 2);
            out << SEP << objs.size() << SEP << points;
            out << SEP << sup << SEP << id;
        }
    }
}

void Simul::reportSimul(std::ostream& out, std::ostream& com) const
{
    com << COM << ljust("class", 2, 2) << SEP << "objects" << SEP << "vertices";
    com << SEP << "largest" << SEP << "identity";
    reportSet(out,  fibers, properties);
    reportSet(out,  solids, properties);
    reportSet(out, spheres, properties);
    reportSet(out,   beads, properties);
}

/**
 Export position of all organizers
 */
void Simul::reportOrganizer(std::ostream& out, std::ostream& com) const
{
    com << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");

    for ( Organizer const* obj=organizers.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->property()->number();
        out << SEP << obj->identity();
        out << SEP << obj->position();
        out << SEP << obj->nbOrganized();
    }
}


/**
 Export position of Asters
 */
void Simul::reportAster(std::ostream& out, std::ostream& com) const
{
    com << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
    
    for ( Organizer const* obj=organizers.first(); obj; obj=obj->next() )
    {
        if ( obj->tag() == Organizer::ASTER_TAG )
        {
            out << LIN << obj->property()->number();
            out << SEP << obj->identity();
            out << SEP << obj->position();
        }
    }
}


/**
 Export position of Beads
 */
void Simul::reportBeadPosition(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
    
    for ( Bead const* obj=beads.first(); obj; obj=obj->next() )
    {
       if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->position();
    }
}


/**
 Export number of beads classified as a function of
 the number of grafted Single that are attached to Fibers
 */
void Simul::reportBeadSingles(std::ostream& out, std::ostream& com) const
{
    com << COM << "identity" << "amount(nb_attached_hands)";
    
    std::map<ObjectID, int> cnt;
    
    for ( Single const* i=singles.firstA(); i; i=i->next() )
    {
        Bead const* obj = Bead::toBead(i->base());
        if ( obj )
            ++cnt[ obj->identity() ];
    }

    const int max = 12;
    int nb[max] = { 0 };
    for ( Bead const* obj=beads.first(); obj; obj=obj->next() )
        ++nb[ cnt[obj->identity()] ];
    
    for ( int c = 0; c < max; ++c )
        out << " " << std::setw(3) << nb[c];
}


/**
 Export position of Solids
 */
void Simul::reportSolidPosition(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity" << SEP << repeatXYZ("cen");
    com << SEP << repeatXYZ("point1") << SEP << repeatXYZ("point2");
        
    for ( Solid const* obj=solids.first(); obj; obj=obj->next() )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->centroid();
        out << SEP << obj->posPoint(0);
        if ( obj->nbPoints() > 1 )
        out << SEP << obj->posPoint(1);
        
        if ( modulo )
        {
            Vector pos = obj->centroid();
            modulo->fold(pos);
            out << SEP << pos;
        }
    }
}

/**
 Export orientation of Solids
 */
void Simul::reportSolidOrientation(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity" << SEP << repeatXYZ("cen");
    com << SEP << repeatXYZ("dir");
        
    for ( Solid const* obj=solids.first(); obj; obj=obj->next() )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->centroid();
        out << SEP << obj->orientation();
    }
}

/**
 Export position of Solids with counts of Hands and attached Hands
 */
void Simul::reportSolidHands(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
    com << SEP << "nb_hand" << SEP << "nb_link";
        
    for ( Solid const* obj = solids.firstID(); obj; obj = solids.nextID(obj) )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        Vector pos = obj->centroid();
        if ( modulo ) modulo->fold(pos);
        out << SEP << pos;
        SingleList anchored = singles.collectWrists(obj);
        int cnt = 0;
        for ( Single const* s : anchored )
            cnt += s->attached();
        out << SEP << anchored.size() << SEP << cnt;
    }
}


/**
 Report position of Sphere
 */
void Simul::reportSpherePosition(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity";
    com << SEP << repeatXYZ("center") << SEP << repeatXYZ("point2");
        
    for ( Sphere const* obj=spheres.first(); obj; obj=obj->next() )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->posPoint(0);
        if ( obj->nbPoints() > 1 )
        out << SEP << obj->posPoint(1);
    }
}


/**
 Report something about Space (incomplete)
 */
void Simul::reportSpace(std::ostream& out, std::ostream& com) const
{
    com << COM << "class" << SEP << "identity";
    
    for ( Space const* obj=spaces.firstID(); obj; obj=spaces.nextID(obj) )
    {
        out << LIN << obj->prop->name();
        out << SEP << obj->identity();
        out << SEP << std::fixed << obj->prop->shape;
    }
}


/**
 Report force on Space
 */
void Simul::reportSpaceForce(std::ostream& out, std::ostream& com) const
{
    com << COM << "class" << SEP << "identity" << SEP << "shape";
    
    for ( Space const* obj=spaces.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->name();
        out << SEP << obj->identity();
        out << SEP << obj->prop->shape;
        obj->report(out);
    }
}


/**
 Report quantity of substance in Field
 */
void Simul::reportField(std::ostream& out, std::ostream& com) const
{
    if ( fields.size() == 0 )
        return;
    com << COM << ljust("class", 2, 2);
    com << SEP << "total" << SEP << "avg" << SEP << "min" << SEP << "max";
    
    FieldCell s, a, n, x;
    // report total substance in each Field
    for ( Field const* obj=fields.first(); obj; obj=obj->next() )
    {
        const real alpha = 1.0 / obj->cellVolume();
        obj->infoValues(s, a, n, x);
        out << LIN << ljust(obj->prop->name(), 2);
        out << SEP << s;
        out << SEP << a * alpha;
        out << SEP << n * alpha;
        out << SEP << x * alpha;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Single

void writeSingle(std::ostream& out, Single const* obj, Simul const* simul)
{
    out << LIN << obj->prop->number();
    out << SEP << obj->identity();
    out << SEP << obj->position();
    Fiber const* fib = obj->fiber();
    if ( fib )
    {
        out << SEP << obj->force();
        out << SEP << fib->identity();
        out << SEP << obj->abscissa();
        out << SEP << simul->organizers.findOrganizerID(fib);
    }
    else
    {
        out << SEP << Vector(0,0,0);
        out << SEP << "0";
        out << SEP << "nan";
        out << SEP << "0";
    }
}


/**
 Export details of Singles, possiby selecting for a certain kind
 */
void Simul::reportSingleState(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity";
    com << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force");
    com << SEP << "fiber" << SEP << "abscissa" << SEP << "aster";
    
    for ( Single const* obj = singles.firstID(); obj; obj = singles.nextID(obj) )
        if ( !sel || sel == obj->prop )
            writeSingle(out, obj, this);
}


/**
 Export details of attached Singles
 */
void Simul::reportSinglePosition(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
    com << SEP << "fiber" << SEP << "abscissa";
        
    for ( Single const* obj = singles.firstID(); obj; obj = singles.nextID(obj) )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->posFoot();
        if ( obj->attached() )
        {
            out << SEP << obj->fiber()->identity();
            out << SEP << obj->abscissa();
        }
        else
        {
            out << SEP << 0 << SEP << 0;
        }
    }
}

/**
 Export details of attached Singles
 */
void Simul::reportSingleLink(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity";
    com << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force");
    com << SEP << "fiber" << SEP << "abscissa" << SEP << "aster";
        
    for ( Single const* obj = singles.firstID(); obj; obj = singles.nextID(obj) )
    {
        if ( obj->attached()  && ( !sel || sel == obj->prop ))
            writeSingle(out, obj, this);
    }
}


/**
 Export number of Single in each state
 */
void Simul::reportSingle(std::ostream& out, std::ostream& com, Property const* sel) const
{
    constexpr PropertyID SUP = 128;
    
    size_t free[SUP+1] = { 0 };
    size_t bound[SUP+1] = { 0 };
    size_t based[SUP+1] = { 0 };
    
    for ( Single const* i = singles.firstF(); i ; i = i->next() )
    {
        assert_true(!i->attached());
        PropertyID x = std::min(i->prop->number(), SUP);
        ++free[x];
        based[x] += ( i->base() != nullptr );
    }
    
    for ( Single const* i=singles.firstA(); i ; i=i->next() )
    {
        assert_true(i->attached());
        PropertyID x = std::min(i->prop->number(), SUP);
        ++bound[x];
        based[x] += ( i->base() != nullptr );
    }
    
    com << COM << ljust("single", 2, 2);
    com << SEP << "total";
    com << SEP << "based";
    com << SEP << "free";
    com << SEP << "bound";
    
    for ( Property const* i : properties.find_all("single") )
    {
        if ( sel && i != sel )
            continue;
        out << LIN << ljust(i->name(), 2);
        PropertyID x = i->number();
        if ( x < SUP )
        {
            out << SEP << free[x] + bound[x];
            out << SEP << based[x];
            out << SEP << free[x];
            out << SEP << bound[x];
        }
        else
            out << SEP << " out-of-range ";
    }
}


/**
 Export average properties of Couples forces
 */
void Simul::reportSingleForce(std::ostream& out, std::ostream& com, Property const* sel) const
{
    constexpr PropertyID SUP = 8;
    real cnt[SUP+1] = { 0 };
    real avg[SUP+1] = { 0 };
    real sup[SUP+1] = { 0 };
    real len[SUP+1] = { 0 };

    // accumulate counts:
    for ( Single const* i=singles.firstA(); i; i=i->next() )
    {
        if ( i->hasLink() && ( !sel || sel == i->prop ))
        {
            PropertyID x = std::min(i->prop->number(), SUP);
            real f = i->force().norm();
            avg[x] += f;
            cnt[x] += 1;
            sup[x] = std::max(sup[x], f);
            len[x] = std::max(len[x], i->stretch().norm());
        }
    }
    
    com << COM << ljust("single", 2, 2) << SEP << "avg_force" << SEP << "max_force" << SEP << "max_len";
    for ( PropertyID i = 0; i < SUP; ++i )
    {
        if ( cnt[i] > 0 )
        {
            Property const* p = properties.find_or_die("single", i);
            out << LIN << ljust(p->name(), 2);
            out << SEP << avg[i] / cnt[i];
            out << SEP << sup[i];
            out << SEP << len[i];
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Couple

void Simul::reportCoupleList(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "couples" << SEP << "count";

    size_t cnt = couples.collect(match_property, sel).size();

    if ( sel )
        out << LIN << sel->name();
    else
        out << LIN << "all";

    out << SEP << cnt;
 }


void writeCouple(std::ostream& out, Couple const* obj)
{
    out << LIN << obj->prop->number();
    out << SEP << obj->identity();
    out << SEP << obj->active();
    out << SEP << obj->position();

    Fiber const* fib = obj->fiber1();
    if ( fib )
    {
        out << SEP << fib->identity();
        out << SEP << obj->abscissa1();
    }
    else
    {
        out << SEP << "0";
        out << SEP << "nan";
    }

    fib = obj->fiber2();
    if ( fib )
    {
        out << SEP << fib->identity();
        out << SEP << obj->abscissa2();
    }
    else
    {
        out << SEP << "0";
        out << SEP << "nan";
    }
}

        
/**
 Export position of Couples of a certain kind
 */
void Simul::reportCoupleState(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity" << SEP << "active" << SEP << repeatXYZ("pos");
    com << SEP << "fiber1" << SEP << "abscissa1" << SEP << "fiber2" << SEP << "abscissa2";
    
    for ( Couple const* obj=couples.firstFF(); obj ; obj=obj->next() )
        if ( !sel || sel == obj->prop )
            writeCouple(out, obj);
    
    for ( Couple const* obj=couples.firstAF(); obj ; obj=obj->next() )
        if ( !sel || sel == obj->prop )
            writeCouple(out, obj);
    
    for ( Couple const* obj=couples.firstFA(); obj ; obj=obj->next() )
        if ( !sel || sel == obj->prop )
            writeCouple(out, obj);
    
    for ( Couple const* obj=couples.firstAA(); obj ; obj=obj->next() )
        if ( !sel || sel == obj->prop )
            writeCouple(out, obj);
}


/**
 Export position of active Couples of a certain kind
 */
void Simul::reportCoupleActive(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "state" << SEP << repeatXYZ("pos");
    
    for ( Couple const* obj=couples.firstFF(); obj ; obj=obj->next() )
        if ( obj->active()  &&  obj->prop == sel )
            out << "\n 0 " << obj->position();
   
    for ( Couple const* obj=couples.firstAF(); obj ; obj=obj->next() )
        if ( obj->prop == sel )
            out << "\n 1 " << obj->position();
    
    for ( Couple const* obj=couples.firstFA(); obj ; obj=obj->next() )
        if ( obj->prop == sel )
            out << "\n 2 " << obj->position();
    
    for ( Couple const* obj=couples.firstAA(); obj ; obj=obj->next() )
        if ( obj->prop == sel )
            out << "\n 3 " << obj->position();
}


/**
 Export position and force of Couples that are bound to 2 filaments
 */
void Simul::reportCoupleLink(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "class" << SEP << "identity";
    com << SEP << "fiber1" << SEP << "abscissa1";// << SEP << repeatXYZ("pos1");
    com << SEP << "fiber2" << SEP << "abscissa2";// << SEP << repeatXYZ("pos2");
    com << SEP << "force" << SEP << "cos_angle";
        
    for ( Couple const* obj=couples.firstAA(); obj ; obj=obj->next() )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        
        out << SEP << obj->fiber1()->identity();
        out << SEP << obj->abscissa1();
        //out << SEP << obj->posHand1();

        out << SEP << obj->fiber2()->identity();
        out << SEP << obj->abscissa2();
        //out << SEP << obj->posHand2();

        out << SEP << obj->force().norm();
        out << SEP << obj->cosAngle();
    }
}


/**
 Export configuration of bridging couple, as
 P: parallel side-side links
 A: antiparallel side-side links
 X: other side-side links
 T: side-end
 V: end-end
 
 T and V are defined with respect to the `end`, at distance `threshold`,
 both can be set as parameters.
 
 by Jamie Li Rickman for
 Determinants of Polar versus Nematic Organization in Networks of Dynamic Microtubules
 and Mitotic Motors, Cell 2018
 */
void Simul::reportCoupleConfiguration(std::ostream& out, std::ostream& com, Property const* sel,
                                      Glossary& opt) const
{
    real dis = 0.010;  // max distance to end to constitute T or V link
    opt.set(dis, "distance", "threshold");
    
    size_t T[8] = { 0 };
    for ( Couple const* obj=couples.firstAA(); obj ; obj=obj->next() )
    {
        if ( !sel || sel == obj->prop )
            ++T[obj->configuration(dis)];
    }
    size_t S = T[0]+T[1]+T[2]+T[3]+T[4]+T[5]+T[6];
    
    com << COM << "couples" << SEP << "Total" << SEP << "P" << SEP << "A";
    com << SEP << "X" << SEP << "T+" << SEP << "V+" << SEP << "T-" << SEP << "V-";

    if ( sel )
        out << LIN << sel->name();
    else
        out << LIN << "all";

    out << SEP << S << SEP << T[0] << SEP << T[1] << SEP << T[2] << SEP << T[3] << SEP << T[4] << SEP << T[5] << SEP << T[6];
 }



/**
 Export average properties of Couples forces
 */
void Simul::reportCoupleForce(std::ostream& out, std::ostream& com, Property const* sel) const
{
    constexpr PropertyID SUP = 8;
    real cnt[SUP+1] = { 0 };
    real avg[SUP+1] = { 0 };
    real sup[SUP+1] = { 0 };
    real len[SUP+1] = { 0 };

    // accumulate counts:
    for ( Couple const* i=couples.firstAA(); i ; i = i->next() )
    {
        if ( !sel || sel == i->prop )
        {
            PropertyID x = std::min(i->prop->number(), SUP);
            real f = i->force().norm();
            avg[x] += f;
            cnt[x] += 1;
            sup[x] = std::max(sup[x], f);
            len[x] = std::max(len[x], i->stretch().norm());
        }
    }
        
    com << COM << ljust("couple", 2) << SEP << "count" << SEP << "avg_force" << SEP << "max_force" << SEP << "max_len";
    for ( PropertyID i = 0; i < SUP; ++i )
    {
        if ( cnt[i] > 0 )
        {
            Property const* p = properties.find_or_die("couple", i);
            out << LIN << ljust(p->name(), 2);
            out << SEP << (int)cnt[i];
            out << SEP << avg[i] / cnt[i];
            out << SEP << sup[i];
            out << SEP << len[i];
        }
    }
}


/**
 Export histogram of Couples forces
 */
void Simul::reportCoupleForceHistogram(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    const PropertyID SUP = 8;
    const unsigned MAX = 256;
    size_t cnt[SUP+1][MAX+1];

    real delta = 0.5;
    unsigned nbin = 64;
    opt.set(delta, "interval");
    opt.set(nbin, "interval", 1);
    nbin = std::min(nbin, MAX);

    // reset counts:
    for ( PropertyID ii = 0; ii < SUP; ++ii )
    for ( unsigned jj = 0; jj <= nbin; ++jj )
        cnt[ii][jj] = 0;
    
    // accumulate counts:
    for ( Couple const* i=couples.firstAA(); i ; i = i->next() )
    {
        PropertyID x = std::min(i->prop->number(), SUP);
        unsigned f = (unsigned)( i->force().norm() / delta );
        if ( f < nbin )
            ++cnt[x][f];
        else
            ++cnt[x][nbin];
    }
    
    if ( 1 )
    {
        com << COM << "force_distribution" << " (`scale` indicates the center of each bin)";
        out << LIN << ljust("scale", 2);
        std::streamsize p = out.precision(3);
        for ( size_t u = 0; u <= nbin; ++u )
            out << " " << std::setw(5) << delta * ( u + 0.5 );
        out.precision(p);
    }
    
    for ( PropertyID x = 0; x < SUP; ++x )
    {
        size_t sum = 0;
        for ( size_t jj = 0; jj <= nbin; ++jj )
            sum += cnt[x][jj];
        if ( sum )
        {
            Property const* p = properties.find_or_die("couple", x);
            out << LIN << ljust(p->name(), 2);
            for ( size_t jj = 0; jj <= nbin; ++jj )
                out << ' ' << std::setw(5) << cnt[x][jj];
        }
    }
}


/**
 Export number of Couples in each state
 */
void Simul::reportCouple(std::ostream& out, std::ostream& com, Property const* sel) const
{
    constexpr PropertyID SUP = 128;
    int act[SUP] = { 0 }, cnt[SUP][4];
    
    //reset counts:
    for ( PropertyID i = 0; i < SUP; ++i )
    {
        cnt[i][0] = 0;
        cnt[i][1] = 0;
        cnt[i][2] = 0;
        cnt[i][3] = 0;
    }
    
    for ( Couple const* i=couples.firstFF(); i ; i = i->next() )
    {
        assert_true(!i->attached1() && !i->attached2());
        PropertyID x = std::min(i->prop->number(), SUP);
        act[x] += ( i->active() );
        ++(cnt[x][0]);
    }
    
    for ( Couple const* i=couples.firstAF(); i ; i = i->next() )
    {
        assert_true(i->attached1() && !i->attached2());
        PropertyID x = std::min(i->prop->number(), SUP);
        act[x] += ( i->active() );
        ++(cnt[x][1]);
    }
    for ( Couple const* i=couples.firstFA(); i ; i = i->next() )
    {
        assert_true(!i->attached1() && i->attached2());
        PropertyID x = std::min(i->prop->number(), SUP);
        act[x] += ( i->active() );
        ++(cnt[x][2]);
    }
    
    for ( Couple const* i=couples.firstAA(); i ; i = i->next() )
    {
        assert_true(i->attached1() && i->attached2());
        PropertyID x = std::min(i->prop->number(), SUP);
        act[x] += ( i->active() );
        ++(cnt[x][3]);
    }
    
    com << COM << ljust("couple", 2, 2);
    com << SEP << "total" << SEP << "active";
    com << SEP << "FF" << SEP << "AF" << SEP << "FA"<< SEP << "AA";
    
    for ( Property const* i : properties.find_all("couple") )
    {
        if ( sel && i != sel )
            continue;
        out << LIN << ljust(i->name(), 2);
        PropertyID x = i->number();
        if ( x < SUP )
        {
            out << SEP << cnt[x][0]+cnt[x][1]+cnt[x][2]+cnt[x][3];
            out << SEP << act[x];
            for ( size_t d = 0; d < 4; ++d )
                out << SEP << cnt[x][d];
        }
        else
            out << SEP << "out-of-range";
    }
}


/**
 Export composition of Couple classes
 */
void Simul::reportCoupleAnatomy(std::ostream& out, std::ostream& com, Property const* sel) const
{
    com << COM << "hand_id" << SEP << rjust("hand_name", 2);
    
    for ( Property const* i : properties.find_all("hand") )
    {
        HandProp const* p = static_cast<HandProp const*>(i);
        out << LIN << p->number();
        out << SEP << rjust(p->name(), 2);
    }

    com << COM << "couple_id" << SEP << rjust("couple_name", 2);
    com << SEP << rjust("hand1", 2) << SEP << rjust("hand2", 2);
    
    for ( Property const* i : properties.find_all("couple") )
    {
        if ( sel && i != sel )
            continue;
        CoupleProp const* p = static_cast<CoupleProp const*>(i);
        out << LIN << p->number();
        out << SEP << rjust(p->name(), 2);
        out << SEP << rjust(p->hand1_prop->name(), 2);
        out << SEP << rjust(p->hand2_prop->name(), 2);
    }
}

//------------------------------------------------------------------------------
#pragma mark - Cluster Analysis

/**
 equalize flag() if connected by a Couple
 */
void Simul::flagClustersCouples() const
{
    for ( Couple const* X=couples.firstAA(); X ; X=X->next() )
    {
        ObjectFlag f = X->fiber1()->flag();
        ObjectFlag g = X->fiber2()->flag();
        if ( f != g )
            changeFlags(std::max(f, g), std::min(f, g));
    }
}

/**
 equalize flag() if connected by a Couple of given type
 */
void Simul::flagClustersCouples(Property const* sel) const
{
    for ( Couple const* X=couples.firstAA(); X ; X=X->next() )
    {
        if ( X->prop == sel )
        {
            ObjectFlag f = X->fiber1()->flag();
            ObjectFlag g = X->fiber2()->flag();
            if ( f != g )
                changeFlags(std::max(f, g), std::min(f, g));
        }
    }
}

/**
 equalize flag() if connected through Single of class Wrist:
 */
void Simul::flagClustersSingles() const
{
    for ( Single * X=singles.firstA(); X; X=X->next() )
    {
        Mecable const* B = X->base();
        if ( B )
        {
            ObjectFlag f = X->fiber()->flag();
            ObjectFlag g = B->flag();
            if ( f != g )
                changeFlags(std::max(f, g), std::min(f, g));
        }
    }
}


void Simul::flagClusters(bool C, bool S, bool M) const
{
    if ( ! ( C | S | M ) )
        throw InvalidSyntax("you must specify a cluster type: couple=1 or solid=1 or meca=1");

    setUniqueFlags();
    if ( C ) flagClustersCouples();
    if ( S ) flagClustersSingles();
    if ( M ) flagClustersMeca();
    flagLargestCluster(1UL);
}


/// class to store info about a Cluster
struct Cluster
{
    size_t     cnt;
    ObjectFlag flg;

    Cluster(ObjectFlag f, size_t n) : cnt(n), flg(f) {}
        
    /// Compare function
    bool operator < (Cluster const&b) const
    {
        if ( cnt == b.cnt ) return flg < b.flg;
        return cnt > b.cnt;
    }
};


/**
Set Mecable::flag() to 'f' for Mecables in the largest cluster
*/
void Simul::flagLargestCluster(ObjectFlag f) const
{
    // for large systems, it would be better to use a std:unordered_map here:
    typedef std::map<ObjectFlag, size_t> map_t;
    map_t map;

    // collect number of fibers with same 'flag' value:
    for ( Fiber* fib = fibers.first(); fib; fib = fib->next() )
        ++map[fib->flag()];
    
    size_t size = 0;
    ObjectFlag largest = 0;
    // find largest clusters in 'map':
    for ( map_t::value_type const& i : map )
    {
        if ( i.second > size )
        {
            largest = i.first;
            size = i.second;
        }
    }
    
    if ( size > 2 )
    {
        // swap 'largest' and 'f':
        for ( Fiber* fib = fibers.first(); fib; fib = fib->next() )
        {
            if ( fib->flag() == f )
                fib->flag(largest);
            else if ( fib->flag() == largest )
                fib->flag(f);
        }
    }
}


/// to sort fibers in order of identity
bool comp_fibers(Fiber const* i, Fiber const* j)
{
    return ( i->identity() < j->identity() );
}

/**
Order Clusters from 1 to max, in order of decreasing size,
and set fiber:flag() to the corresponding cluster index.
@return number of clusters
*/
size_t reportOrderedClusters(std::ostream& out, std::ostream& com, Fiber* first, size_t threshold, size_t details)
{
    typedef std::vector<Fiber*> list_t;
    typedef std::map<ObjectFlag, list_t> map_t;
    // the std::set will keep its elements ordered:
    typedef std::multiset<Cluster> set_t;
    map_t map;
    set_t clusters;

    // extract clusters in 'map' and reset fiber's flag:
    for ( Fiber* F = first; F; F = F->next() )
    {
        //std::clog << fib->reference() << " " << fib->flag() << "\n";
        map[F->flag()].push_back(F);
        F->flag(0);
    }
    
    // insert clusters with size information to get them sorted:
    for ( map_t::value_type const& i : map )
    {
        size_t s = i.second.size();
        if ( s >= threshold )
            clusters.emplace(i.first, s);
    }
    
    if ( details > 1 )
    {
        com << COM << "cluster_id" << SEP << "nb_fibers :" << SEP << "fiber_id";
        out << LIN << clusters.size() << " clusters:";
    }
    else if ( details > 0 )
    {
        com << COM << "cluster_id" << SEP << "nb_fibers";
        out << LIN << clusters.size() << " clusters:";
    }
    
    // consider clusters by decreasing size:
    ObjectFlag idx = 0;
    for ( set_t::value_type const& c : clusters )
    {
        ++idx;
        list_t & list = map[c.flg];
        
        for ( Fiber * i : list )
            i->flag(idx);

        if ( details > 0 )
        {
            std::sort(list.begin(), list.end(), comp_fibers);
            out << LIN << c.flg << SEP << c.cnt << ":";
            size_t cnt = std::min(list.size(), details);
            for ( size_t i = 0; i < cnt; ++i )
                out << " " << list[i]->identity();
        }
    }
    return idx;
}


/**
 Export clusters defined by Simul::flagClusters()
 Clusters are ordered in decreasing size.
 */
void Simul::reportClusters(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    size_t details = 128;
    bool C = false, S = false, M = false;

    opt.set(details, "details");
    opt.set(C, "couples", "couple");
    opt.set(S, "singles", "single");
    opt.set(M, "meca");
    
    flagClusters(C, S, M);
    
    com << COM << "cluster by couples " << C << " solids " << S << " meca " << M;
    reportOrderedClusters(out, com, fibers.first(), 2, details);
}


//------------------------------------------------------------------------------
#pragma mark - Platelet Microtubule Ring

/**
 Evaluates if the Fiber distribution makes a connected ring around the Z-axis
 @returns number of rings
 FJN 8.07.2017, 8.11.2018, for Blood Platelet project
 */
size_t Simul::flagRing() const
{
    flagClusters(1, 0, 0);

    typedef std::list<ObjectFlag> list_t;
    list_t ring;
    
    // include all cluster into 'ring'
    for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
    {
        ObjectFlag f = fib->flag();
        if ( f > 0 )
        {
            list_t::const_iterator i = std::find(ring.begin(), ring.end(), f);
            if ( i == ring.end() )
                ring.push_back(f);
        }
    }
    
    // rotate plane around the Z-axis and find intersecting fibers
    for ( unsigned a = 0; a < 360; ++a )
    {
        real ang = a * M_PI / 180.0;
        Vector nor( std::cos(ang), std::sin(ang), 0.0);
        Vector dir(-std::sin(ang), std::cos(ang), 0.0);
        
        list_t sec;
        for ( Fiber const* fib=fibers.first(); fib; fib=fib->next() )
        {
            for ( index_t s = 0; s < fib->nbSegments(); ++s )
            {
                // check that fiber intersect with plane:
                real abs = fib->planarIntersect(s, nor, 0);
                if ( 0 <= abs  &&  abs < 1 )
                {
                    // check that intersection is located on 'dir' side of Z-axis:
                    Vector pos = fib->midPoint(s, abs);
                    if ( dot(pos, dir) > 0 )
                    {
                        // transfer cluster if it was already present before:
                        ObjectFlag f = fib->flag();
                        list_t::iterator i = std::find(ring.begin(), ring.end(), f);
                        if ( i != ring.end() )
                        {
                            ring.erase(i);
                            sec.push_back(f);
                        }
                    }
                }
            }
        }
        ring = sec;
    }
    
    if ( ring.size() == 1 )
    {
        ObjectFlag f = *ring.begin();
        for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
        {
            if ( fib->flag() == f )
                fib->flag(1);
            else
                fib->flag(0);
        }
    }
    else if ( ring.size() > 0 )
    {
        // unflag all fibers that are not part of a ring:
        for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
        {
            if ( std::find(ring.begin(), ring.end(), fib->flag()) == ring.end() )
                fib->flag(0);
        }
    }
    else
    {
        // unflag all
        for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
            fib->flag(0);
    }
    return ring.size();
}


/**
 Calculate the length of the ring and its mean radius
 FJN 15.09.2018, for Blood Platelet project
*/
void Simul::analyzeRing(ObjectFlag flg, real& length, real& radius) const
{
    length = 0.0;
    radius = 0.0;

    Vector cen_old;
    unsigned rad_cnt = 0;
    
    // rotate plane around the Z-axis and find intersecting fibers
    for ( unsigned a = 0; a < 360; ++a )
    {
        real ang = a * M_PI / 180.0;
        Vector nor( std::cos(ang), std::sin(ang), 0);
        Vector dir(-std::sin(ang), std::cos(ang), 0);
        
        unsigned cen_cnt = 0;
        Vector cen(0,0,0);
        for ( Fiber const* fib=fibers.first(); fib; fib=fib->next() )
        {
            // only consider fiber that are part of the ring:
            if ( fib->flag() == flg )
            {
                for ( index_t s = 0; s < fib->nbSegments(); ++s )
                {
                    // check that fiber intersect with plane:
                    real abs = fib->planarIntersect(s, nor, 0);
                    if ( 0 <= abs  &&  abs < 1 )
                    {
                        Vector pos = fib->midPoint(s, abs);
                        // check that intersection is located on 'dir' side of Z-axis:
                        real H = dot(pos, dir);
                        if ( H > 0 )
                        {
                            radius += H;
                            ++rad_cnt;
                            cen += pos;
                            ++cen_cnt;
                        }
                    }
                }
            }
        }
        if ( cen_cnt > 0 )
        {
            cen /= cen_cnt;
            if ( a > 0 )
                length += ( cen - cen_old ).norm();
            cen_old = cen;
        }
    }
    radius /= rad_cnt;
}

/**
Calculate the length of the ring and its radius
*/
void Simul::reportRing(std::ostream& out, std::ostream& com) const
{
    com << COM << "nb_rings" << SEP << "length" << SEP << "radius";
    size_t nring = flagRing();
    if ( nring == 1 )
    {
        real len, rad;
        analyzeRing(1, len, rad);
        out << LIN << 1 << SEP << std::fixed << len<< SEP << rad ;
    }
    else
    {
        out << LIN << nring << SEP << 0.0 << SEP << 0.0;
    }
}


/**
 Evaluates if the Fiber distribution makes a connected ring around the Z-axis
 FJN 8.07.2017, for Blood Platelet project
 */
void Simul::reportPlatelet(std::ostream& out, std::ostream& com) const
{
    size_t nfib = 0;
    real pol = 0, var = 0, mn = INFINITY, mx = -INFINITY, off = 0;
    FiberSet::infoLength(fibers.collect(), nfib, pol, var, mn, mx, off);
    pol *= nfib;
    
    if ( nfib > 1024 )
    {
        // the calculation can take too much time with lots of fibers, so we cut it here:
        com << COM << "nb_fiber" << SEP << "polymer";
        out << LIN << nfib << SEP << std::fixed << pol;
        return;
    }
    
    computeForces();

    real t = 0, ten = 0, inf = 0, sup = 0;
    size_t c, cnt = 0;

    // rotate plane around the Z-axis and find intersecting fibers
    for ( unsigned a = 0; a < 180; a += 10 )
    {
        real ang = a * M_PI / 180.0;
        Vector dir(std::cos(ang), std::sin(ang), 0);
        fibers.infoTension(c, t, inf, sup, dir, 0);
        cnt += 2;  // every plane should intersect the ring twice
        ten += t;
    }
    ten /= (real)cnt;
    
    std::ofstream nos("/dev/null");
    real force = reportFiberConfinement(nos, com);

    real len = 0.0, rad = 0.0;
    if ( flagRing() == 1 )
        analyzeRing(1, len, rad);

    std::streamsize p = out.precision();
    if ( p < 3 ) out.precision(3);
    com << COM << "nb_fiber" << SEP << "polymer" << SEP << "tension" << SEP << "force" << SEP << "length" << SEP << "radius";
    out << LIN << nfib << SEP << std::fixed << pol << SEP << ten << SEP << force << SEP << len << SEP << rad;
    out.precision(p);
}


//------------------------------------------------------------------------------
#pragma mark - Mitotic Spindle

/**
 Export indices calculated by FiberSet::infoSpindle
 */
void Simul::reportSpindleIndices(std::ostream& out, std::ostream& com) const
{
    com << COM << "amount" << SEP << "radial" << SEP << "polar";
    real ixa, ixp;
    fibers.infoSpindle(ixa, ixp, Vector(1,0,0), 0, 10, 0.5);
    out << LIN << fibers.size();
    out << SEP << ixa;
    out << SEP << ixp;
}


/**
 Export average position of beads "condensate" located on left and right sides,
 and the distance between the two barycenters. Right/Left are defined on X-axis.
 */
void Simul::reportSpindleLength(std::ostream& out, std::ostream& com, Glossary&) const
{
    Vector axis(1,0,0);
    BeadProp * bip = findProperty<BeadProp>("bead", "condensate");
    if ( bip )
    {
        com << COM << "left_cnt" << SEP << repeatXYZ("left_");
        com << SEP << "right_cnt" << SEP << repeatXYZ("right_") << SEP << "distance";
        Vector R(0,0,0), L(0,0,0);
        size_t nR = 0, nL = 0;
        
        for ( Object const* i : beads.collect(bip) )
        {
            Vector pos = i->position();
            if ( dot(pos, axis) < 0 ) { L += pos; ++nL; } else { R += pos; ++nR; }
        }
        if ( nL ) L /= nL;
        if ( nR ) R /= nR;
        out << LIN << nL << SEP << L << SEP << nR << SEP << R << SEP << norm(R-L);
    }
}


/**
 Export number and average position of minus ends located on left and right sides,
 and the distance between the two barycenters. Right/Left are defined on X-axis.
 */
void Simul::reportMarkedFiberEnds(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    Vector axis(1,0,0);
    com << COM << ljust("mark", 1, 2) << SEP << "left_cnt" << SEP << repeatXYZ("left_");
    com << SEP << "right_cnt" << SEP << repeatXYZ("right_") << SEP << "distance";

    ObjectMark sup = 0;
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
        sup = std::max(sup, fib->mark());

    for ( ObjectMark k = 0; k <= sup; ++k )
    {
        uintptr_t val = k;
        ObjectList objs = fibers.collect(match_mark, reinterpret_cast<void*>(val));
        if ( objs.size() )
        {
            Vector R(0,0,0), L(0,0,0);
            size_t nR = 0, nL = 0;
            out << LIN << ljust(std::to_string(val), 1);
            for ( Object const* i : objs )
            {
                Vector pos = static_cast<Fiber const*>(i)->posEndM();
                if ( dot(pos, axis) < 0 ) { L += pos; ++nL; } else { R += pos; ++nR; }
            }
            if ( nL ) L /= nL;
            if ( nR ) R /= nR;
            out << SEP << nL << SEP << L << SEP << nR << SEP << R << SEP << norm(R-L);
        }
    }
}


/**
 Export number of Fibers pointing left and right,
 that intersect a plane parallel to YZ.
 The planes are distributed regularly every 0.5 um along the X-axis.
 */
void Simul::reportSpindleProfile(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    com << COM << "position" << SEP << "leftward" << SEP << "rightward" << SEP << "right-left";
    Vector n(1,0,0);
    real m = 10, interval = 0.5;
    opt.set(interval, "interval");
    int R, L;
    for ( real x = -m ; x <= m ; x += interval )
    {
        fibers.infoPlane(R, L, n, -x);
        out << LIN << x;
        out << SEP << L;
        out << SEP << R;
        out << SEP << R-L;
    }
}


/// print some coefficients calculated from the distribution of fibers
void Simul::reportSpindleFitness(std::ostream& out, std::ostream& com, Glossary& opt) const
{
    size_t cnt = 0, left = 0;
    real std = 0, dev = 0;
    com << COM << "kin_left" << SEP << "kin_right" << SEP << "kin_devX" << SEP << "kin_devYZ";
    com << SEP << "bead_left" << SEP << "bead_right" << SEP << "bead_dev";
    real half_length = 4.0;
    opt.set(half_length, "half_length");
    /// check positions of kinetochores (Solid):
    SolidProp * sop = findProperty<SolidProp>("solid", "kinetochore");
    if ( sop )
    {
        ObjectList list = solids.collect(sop);
        for ( Object const* i : list )
        {
            ++cnt;
            Solid const* sol = Solid::toSolid(i);
            Vector pos = sol->posPoint(0);
            left += ( pos.XX < 0 );
            std += square(pos.XX);
            dev += pos.normYZSqr();
        }
        std /= cnt;
        dev /= cnt;
        out << LIN << left << SEP << cnt-left << SEP << std << SEP << dev;
    }
    
    /// check position of condensate (Bead):
    cnt = 0; left = 0;
    BeadProp * bip = findProperty<BeadProp>("bead", "condensate");
    if ( bip )
    {
        ObjectList list = beads.collect(bip);
        for ( Object const* i : list )
        {
            ++cnt;
            Vector pos = i->position();
            left += ( pos.XX < 0 );
            std += square(abs_real(pos.XX) - half_length) + pos.normYZSqr();
        }
        std /= cnt;
        out << SEP << left << SEP << cnt-left << SEP << std;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Misc

/**
 l'angle entre un vecteur 1 (centre du noyau --> SPB)
 et un vecteur 2 (axe de l'hyphe; gauche --> droite = sens du flow).
 */
void Simul::reportAshbya(std::ostream& out, std::ostream& com) const
{
    com << COM << "class id point_0, vector_1, angle";
    for ( Solid const* obj=solids.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->posPoint(0);
        if ( obj->nbPoints() > 1 )
        {
            Vector vec = normalize(obj->diffPoints(0));
            Vector dir(1,0,0);
            out << SEP << vec;
            out << SEP << std::acos(dot(vec, dir));
        }
    }
}


/**
 Categorize the configuration of two microtubules, checking for possible collision
 This calculates 3 boolean values: K = catastrophe, X = crossing, Z = zippered
 1.10.2021 -- 11.2022,
 7.04.2024: added case 'B' if `fib` crosses below `fox`, and 'M' if `fox` has moved
 */
void Simul::reportFiberCollision(std::ostream& out, Fiber* fib, Fiber const* fox, const int print) const
{
    const real NO_CONTACT = 777;
    const real COS20 = 0.94;
    const real COS10 = 0.985;
    const real L1 = 1.0; // threshold used to decide on events
    const real M0 = L1 * L1 * 0.25; // threshold used to decide on 'M'

    static int decided = 0;
    static real abs = -77; // abscissa at first contact on fox
    static real len = -77; // length of fib at first contact
    static real ang = NO_CONTACT; // angle at first contact
    static real dis = INFINITY; // minimum distance reached
    static Vector hit(0,0,0); // position of first contact
    static Vector aim(0,0,0); // direction of plus-end at first contact
    static char kat = 'U'; // category
    
    bool K = 0; // catastrophe
    bool X = 0; // crossing
    bool T = 0; // tangent
    bool Z = 0; // zippering
    bool B = 0; // fib below obstacle fox
    bool M = 0; // obstacle (fox) has moved
    bool P = 0; // fib tip is not deflected
    
    if ( fib && fox )
    {
        const real sup = 3 * fib->prop->steric_radius;
        
        // K = plus-end is shrinking
        K = ( fib->endStateP() == STATE_RED );

        Vector tip = fib->posEndP();
        Vector dir = fib->dirEndP();
        // find closest point to 'fib' plus end, located on 'fox':
        real dpe, aaa = fox->projectPoint(tip, dpe);
        Vector dirFox = fox->dir(aaa);
        dis = std::min(dis, dpe);
        
        // check if obstacle has moved after contact by more than 250 nm:
        if ( abs > 0 )
            M = ( distanceSqr(hit, fox->pos(abs)) > M0 );
        
        // plus-tip of 'fib' is within distance 'sup' from 'fox':
        bool contact = ( dpe < sup*sup );
        if ( contact )
        {
            // check angle of fib's tip with respect to fox at closest point:
            real C = dot(dir, dirFox);
            real A = std::acos(C);
            // the angle and abscissa are set at first contact:
            if ( ang == NO_CONTACT )
            {
                ang = A;
                abs = aaa;
                hit = fox->pos(aaa);
                len = fib->abscissaP();
                aim = dir;
            }
            // 'zippering' implies 'being tangent' and 'significant growth' along the obstacle
            T = ( abs_real(C) > COS20 );   // tangent within 20 degrees
            // the distance zipped is measured along the obstacle, using abscissa:
            Z = ( abs_real(aaa-abs) > 2*L1 ); // has zipped over more than 2*L1
        }
        else if ( ang != NO_CONTACT )
        {
            // the plus tip is growing in the same direction as when contact was made:
            P = ( dot(aim, dir) > COS10 );

            // the plus-tip may have crossed the other filament if it is not in contact
            // consider a point back, and check if it is on opposite side of 'fox'
            real a = min_real(0, 2*len-fib->abscissaP());
            Vector bak = fib->pos(a);
            real ddd, bbb = fox->projectPoint(bak, ddd);
            // bbb = 0.5 * ( aaa + bbb );  // OLD formula before 13.04.2023
            // the abscissa on `fox` below the intersection is estimated with Thales's theorem:
            real alpha = std::sqrt(dpe) / ( std::sqrt(dpe) + std::sqrt(ddd) );
            bbb = aaa + ( bbb - aaa ) * alpha; // abscissa on 'fox'

            Vector mid = fox->pos(bbb);  // position adjacent to the intersection
            Vector axs = fox->dir(bbb);
            Torque TP = cross(axs, tip-mid);
            Torque TM = cross(axs, bak-mid);
            /**
             For X, the plus-tip should be on the other side of the fiber, at least by L1,
             And the angle should be within 20 degree of the angle at contact
            */
            X = ( dot(TP, TM) < 0 );
            
            if ( X )
            {
                // check if 'fib' is below 'fox':
                Space const* spc = fib->prop->confine_space;
                Vector prj = spc->project(mid);  // on the edge
                Vector dwn = spc->normalToEdge(mid); // directed outward
                real ccc = fib->abscissaP() - alpha * 2 * L1; // intersection on 'fib'
                // calculate distances to the edge:
                real foxZ = dot(prj-mid, dwn);
                real fibZ = dot(prj-fib->pos(ccc), dwn);
                B = ( fibZ < foxZ );
#if 0
                // adding a Couple for visual debugging
                CoupleProp * CP = static_cast<CoupleProp*>(properties.find("couple", 1));
                if ( CP && couples.size() < 8 ) {
                    Couple * C = CP->newCouple();
                    C->setPosition(mid);
                    C->hand1()->attachTo(fox, bbb);
                    C->hand2()->attachTo(fib, ccc);
                    const_cast<Simul*>(this)->couples.add(C);
                }
#endif
            }
        }

        if ( kat == 'U' )
        {
            if ( K )
            {
                // rescue microtubule if no contact was ever made:
                if ( ang == NO_CONTACT )
                {
                    fib->setEndStateP(STATE_GREEN);
                    //out << fib->reference() << " rescued\n";
                    return;
                }
                // recorded catastrophes must be at contact
                kat = ( contact ? 'K' : 'k' );
                // in any case, we can stop the simulation
                decided = 1;
            }
            else if ( M ) kat = 'M';
            else if ( X && P ) kat = B?'B':'X';
            // can be 'Z' only if 'X' is not true
            else if ( Z ) kat = 'Z';
            else if ( X ) kat = 'Y';

            // since these states are final, we can terminate the simulation
            if ( strchr("BKMXYZ", kat) )
            {
                decided = 1;
                text_ = std::string(1, kat);
            }
        }
    }

    if ( print )
    {
        out << LIN << std::fixed << std::setprecision(5) << ang;
        out << SEP << std::fixed << std::setprecision(5) << std::sqrt(dis);
        out << SEP << K << SEP << X << SEP << Z << SEP << P << SEP << T << SEP << kat << " ";
    }
    if ( print & 1 )
    {
        // reset static variables for next round:
        decided = 0;
        abs = -77;
        ang = NO_CONTACT;
        dis = INFINITY;
        hit.set(0,0,0);
        aim.set(0,0,0);
        kat = 'U';
        text_ = " ";
    }
    else if ( decided == 1 )
    {
        end_now(0);
        decided = 2;
    }
}


void Simul::reportFiberCollision(std::ostream& out, std::ostream& com, Property const* sel, Glossary& opt) const
{
#if ( DIM == 1 )
    throw InvalidParameter("fiber:collision meaningless in 1D");
#endif
    
    int print = 0;
    opt.set(print, "print");

    if ( fibers.size() > 2 )
        throw InvalidParameter("fiber:collision can only handle 2 fibers");
    
    Fiber * fib = nullptr, * fox = nullptr;
    for ( Fiber * f = fibers.first(); f; f = f->next() )
    {
        if ( f->prop == sel )
            fib = f;
        else
            fox = f;
    }
    // discard condition when fiber has reached its max length
    if ( fib && fib->length()+0.01 > fib->prop->max_length )
        return;
    
    if ( fox && fib )
        reportFiberCollision(out, fib, fox, print);
}


/**
 Export end-to-end distance of Fiber
 */
void Simul::reportSomething(std::ostream& out, std::ostream& com) const
{
    com << COM << "something";
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        Vector ee = fib->posEndP() - fib->posEndM();
        out << ee.norm() << " ";
    }
}
