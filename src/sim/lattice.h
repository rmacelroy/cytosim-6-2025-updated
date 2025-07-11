// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University

#ifndef LATTICE_H
#define LATTICE_H

#include <cmath>
#include <iostream>
#include "assert_macro.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "real.h"


/// Array of discrete sites aligned with the abscissa of a Fiber
/**
 A Lattice can have different purposes:
 - to limit the local density of attached Hand,
 - to simulate traffic-jam of motors moving on Fibers,
 - to simulate binding/unbinding of a continuous substance.
 .
 
 The first 2 points are implemented in a specialized Hand: Digit.
 
 The index of a site can be positive or negative, but always corresponds
 to the abscissa taken along the Fiber:
 
     index = (int) std::floor( abscissa / unit );
 
 The Lattice may hold different types of information at each site, 
 eg. a 'number of molecule' or the concentration of something.
 The linear density can be derived by dividing by the unit length:
 
     density = value(index) / unit();
 
 */
template <typename CELL>
class Lattice
{
public:
    
    /// type used to index the Lattice
    /** This must be a signed integer type */
    typedef int  lati_t;

    /// type stored in each Lattice cell
    typedef CELL cell_t;
    
private:
    
    /// lowest valid index (can be negative or positive)
    lati_t laInf;
    
    /// highest valid index plus one (laSup > laInf)
    lati_t laSup;
    
    /// index of cell containing the minus end
    lati_t laIndexM;
    
    /// index of cell containing the plus end
    lati_t laIndexP;
    
    
    /// Array of sites, of valid range [laInf, laSup[
    /** This pointer uses an offset to match its range : laCell = laPtr - laInf */
    cell_t * laCell;
    
    /// Pointer to allocated memory for cell (not used for access)
    cell_t * laPtr;

    /// distance between adjacent sites
    real     laUnit;
    
    /// index of first cell entirely within minus and plus ends
    lati_t   laEntry;
    
    /// index of last cell entirely within minus and plus ends
    lati_t   laFence;

#pragma mark - Allocation
    
    /// allocate for indices within [inf, sup[, and copy values from `src`
    void allocate_copy(lati_t inf, lati_t sup, cell_t src[], lati_t src_inf, lati_t src_sup)
    {
        //std::clog << this << " Lattice::allocate [" << inf << ", " << sup << "[\n";
        assert_true( inf <= sup );
        
        cell_t * ptr = new cell_t[sup-inf];
        cell_t * mem = ptr - inf;
        
        // reset new array:
        for ( lati_t s = inf; s < sup; ++s )
            mem[s] = 0;

        // transfer `src` data from given range
        if ( src )
        {
            for ( lati_t s = src_inf; s < src_sup; ++s )
            {
                assert_true( inf <= s );
                assert_true( s < sup );
                mem[s] = src[s];
            }
        }
        
        delete[] laPtr;

        laInf = inf;
        laSup = sup;
        laPtr = ptr;
        laCell = mem;
    }
    
    
    /// allocate sites for indices within [inf, sup[, conserving existing indices
    void allocate(lati_t inf, lati_t sup, lati_t margin)
    {
        assert_true( inf <= sup );
        
        // return if existing boundaries are sufficient
        if (( laCell != nullptr ) & ( laInf <= inf ) & ( sup <= laSup ))
            return;
        
        inf -= margin;
        sup += margin;
        
        if ( laCell )
        {
            // only extend the current coverage:
            inf = std::min(inf, laInf);
            sup = std::max(sup, laSup);
            //std::clog<<"Lattice "<<this<<" covers ["<<laInf<<","<<laSup<<"[\n";
        }
        
        allocate_copy(inf, sup, laCell, laInf, laSup);
    }
    
    
    /// free all claimed memory
    void deallocate()
    {
        //std::clog<<"Lattice realeased\n";
        delete[] laPtr;
        laPtr = nullptr;
        laCell = nullptr;
        laInf = 0;
        laSup = 0;
    }
    
    /// forbid access to data via automatic type conversion:
    template<typename T> cell_t   value(const T s) const;
    template<typename T> cell_t&  site(const T s);
    template<typename T> cell_t&  operator[](const T s);
    
#pragma mark - Construction
    
public:
    
    /// Constructor
    Lattice()
    {
        laInf = 0;
        laSup = 0;
        laUnit = 0;
        laPtr = nullptr;
        laCell = nullptr;
        laIndexM = 0;
        laIndexP = 0;
        laEntry = 0;
        laFence = 0;
    }
    
    /// Copy constructor
    Lattice(const Lattice & lat)
    {
        //std::clog << this << " Lattice::Lattice [" << lat.laInf << ", " << lat.laSup << "[\n";
        
        //copy member variables:
        laUnit = lat.laUnit;

        // allocate and copy lat's data:
        allocate_copy(lat.laInf, lat.laSup, lat.laCell, lat.laInf, lat.laSup);
    }
    
    /// assignment operator
    Lattice & operator = (const Lattice & lat)
    {
        laUnit = lat.laUnit;

        // allocate and copy lat's data:
        allocate_copy(lat.laInf, lat.laSup, lat.laCell, lat.laInf, lat.laSup);
        return *this;
    }

    /// Destructor
    ~Lattice() { deallocate(); }
    
    
    /// set distance betwen adjacent sites (the size of a site)
    void setUnit(real u)
    {
        assert_true( u >= REAL_EPSILON );
        laUnit = u;
    }
    
    /// true if lattice unit size was set
    bool ready() const
    {
        return ( laUnit > REAL_EPSILON );
    }

    /// change distance betwen adjacent sites
    void changeUnit(real u)
    {
        if ( u < REAL_EPSILON )
            throw InvalidParameter("lattice:unit must be > 0");
        if ( laUnit != u )
        {
            std::clog << "lattice:unit " << laUnit << " -> " << u << "\n";
            //preserve the same abscissa range, with the new lattice unit
            real s = laUnit * laInf;
            real e = laUnit * laSup + laUnit;
            laUnit = u;
            if ( laCell )
                setRange(s, e);
        }
    }
    
#pragma mark - Edges
    
    /// set cell values outside the valid range
    void markEdges(const cell_t val)
    {
        for ( lati_t i = laInf; i < laIndexM; ++i ) laCell[i] = val;
        for ( lati_t i = laIndexP; i < laSup; ++i ) laCell[i] = val;
    }
    
    /// set the range of valid abscissa
    void setRange(real a, real b)
    {
        assert_true( laUnit > REAL_EPSILON );
#if 0
        if ( !std::is_same<real, cell_t>::value && laCell )
            markEdges(0);
#endif
        laIndexM = index(a);
        laIndexP = index_sup(b) - 1;

        laEntry = index_sup(a);
        laFence = index(b) - 1;

        /* allocate with a safety margin of 8 cells */
        allocate(laIndexM, laIndexP+1, 8);
#if 0
        if ( !std::is_same<real, cell_t>::value )
            markEdges(~0);
#endif
        //std::clog << this << ":setRange(" << a << ", " << b << ") [" << laEntry << ", " << laFence << "]\n";
    }

    /// index of site containing the minus end
    lati_t indexM() const { return laIndexM; }
    
    /// index of site containing the plus end
    lati_t indexP() const { return laIndexP; }
    
    /// index of lowest cell that is entirely within minus and plus ends
    lati_t entry()  const { return laEntry; }
    
    /// index of highest cell that is entirely within minus and plus ends
    lati_t fence()  const { return laFence; }

    /// first valid index
    lati_t inf()    const { return laInf; }
    
    /// one + last valid index
    lati_t sup()    const { return laSup; }

    /// distance between adjacent sites
    real   unit()   const { return laUnit; }
    
#pragma mark - Index / Abscissa

    /// index of the site containing abscissa `a`
    lati_t index(real a)       const { return (lati_t)std::floor(a/laUnit); }
    
    /// index of the site just above the one containing abscissa `a`
    lati_t index_sup(real a)   const { return (lati_t)std::ceil(a/laUnit); }

    /// index of the site closest to abscissa `a`
    lati_t index_round(real a) const { return (lati_t)round(a/laUnit); }

    /// true if site `i` is covered by the lattice allocated range
    bool valid(lati_t i)     const { return (( laInf <= i ) & ( i < laSup )); }
    
    /// true if site `i` is not covered by the lattice allocated range
    bool invalid(lati_t i)   const { return (( i < laInf ) | ( laSup <= i )); }
    
    /// true if site `i` is completely between Minus and Plus ends
    bool betweenMP(lati_t i) const { return (( laEntry <= i ) & ( i <= laFence )); }
    
    /// return 2 if site `i` is partly or entirely below the minus end, 1 if above the plus end
    int outsideMP(lati_t i) const { return 2*( i < laEntry ) | ( laFence < i ); }
    
    /// true if site `i` is partly or entirely below the minus end
    bool belowM(lati_t i) const { return ( i < laEntry ); }
    
    /// true if site `i` is partly or entirely above the plus end
    bool aboveP(lati_t i) const { return ( laFence < i ); }

    
    /// the site of index `h` covers the abscissa range `unit * h < s < unit * ( h + 1 )`
    /**
     abscissa(h) returns the abscissa corresponding to the beginning of the site.
     The range covered by site 'h' is [ abscissa(h), abscissa(h+1) ], and the
     abscissa of the center is abscissa(h+0.5).
     */
    real abscissa(const real s) const { return s * laUnit; }

#pragma mark - data access

    /// pointer to data array
    cell_t*       data()        const  { return laCell; }
    
    /// modifiable value at index `s`, equivalent to []
    cell_t&       data(lati_t s)       { assert_true(valid(s)); return laCell[s]; }

    /// constant value at index `s`, equivalent to []
    cell_t const& data(lati_t s) const { assert_true(valid(s)); return laCell[s]; }
    
    /// reference to Site at index s
    cell_t& operator[](lati_t s)       { assert_true(valid(s)); return laCell[s]; }
    
    /// value at abscissa `a`, with convertion to site index, unlike operator []
    cell_t&       cell(real a)         { lati_t s=index(a); assert_true(valid(s)); return laCell[s]; }
    
    /// value at abscissa `a`, with convertion to site index, unlike operator []
    cell_t const& cell(real a)   const { lati_t s=index(a); assert_true(valid(s)); return laCell[s]; }

    /// set all sites to `value`
    void clear(cell_t value = 0)
    {
        for ( lati_t s = laInf; s < laSup; ++s )
            laCell[s] = value;
    }
    
    /// set gradient of values from zero to `alpha`
    void set_gradient(cell_t alpha)
    {
        cell_t val = alpha * laUnit;
        for ( lati_t s = laInf; s < laSup; ++s )
            laCell[s] = s * val;
    }

#pragma mark - Transfer
    
    /// transfer cells within `[s, e[` to *this, and reset transferred values
    void take(Lattice<CELL> & lat, lati_t is, lati_t ie, cell_t zero)
    {
        //std::clog << " Lattice::take [" << is << ", " << ie << "[\n";
        // select valid range:
        is = std::max(is, std::max(laInf, lat.laInf));
        ie = std::min(ie, std::min(laSup, lat.laSup));
        
        CELL * s = laCell + is;
        CELL *const e = laCell + ie;
        CELL * o = lat.laCell + is;

        // transfer values
        while ( s < e )
        {
            *s = *o;
            *o = zero;
            ++s;
            ++o;
        }
    }
    
    
    /// transfer values in [inf, e[
    void takeM(Lattice<CELL> & lat, const lati_t e)
    {
        take(lat, laInf, e, 0);
    }
    
    /// transfer values in [s, sup[
    void takeP(Lattice<CELL> & lat, const lati_t s)
    {
        take(lat, s, laSup, 0);
    }
    
#pragma mark - Collect

    /// sum all sites in `[s, e[`, set them to zero and return total in `res`
    template <typename SUM>
    void collectCells(SUM& res, lati_t is, lati_t ie)
    {
        res = 0;
        // limit boundaries
        CELL * s = laCell + std::max(is, laInf);
        CELL * e = laCell + std::min(ie, laSup);
        
        // collect
        while ( s < e )
        {
            res += (SUM)*s;
            *s = 0;
            ++s;
        }
    }
    
    
    /// sum all sites that are entirely below the minus end
    template <typename SUM>
    void collectM(SUM& res)
    {
        collectCells(res, laInf, laIndexM);
    }
    
    /// sum all sites that are entirely above the plus end
    template <typename SUM>
    void collectP(SUM& res)
    {
        collectCells(res, laIndexP+1, laSup);
    }
    
    /// sum of values in the entire lattice; set sites to zero
    template <typename SUM>
    void collect(SUM& res) const
    {
        collectCells(res, laInf, laSup);
    }
    
    //------------------------------------------------------------------------------
    #pragma mark - Sums

    /// sum of values for all sites in `[s, e[` and return total in `res`
    template <typename SUM>
    void sum(SUM& res, lati_t is, lati_t ie) const
    {
        // check boundaries
        CELL * s = laCell + std::max(is, laInf);
        CELL * e = laCell + std::min(ie, laSup);
        
        // collect
        res = 0;
        while ( s < e )
        {
            res += *s;
            ++s;
        }
    }
    
    /// sum of values for all sites below `e` (not-included)
    template <typename SUM>
    void sumM(SUM& res, const lati_t e) const
    {
        sum(res, laInf, e);
    }

    /// sum of values for all sites above `s` (included)
    template <typename SUM>
    void sumP(SUM& res, const lati_t s) const
    {
        sum(res, s, laSup);
    }

    
    /// sum of values in the entire lattice
    template <typename SUM>
    void sum(SUM& res) const
    {
        sum(res, laInf, laSup);
    }

    /// sum of values in the entire lattice using 'real' as accumulator
    real sum() const
    {
        real res;
        sum(res, laInf, laSup);
        return res;
    }

#pragma mark - I/O
    
    /// write data within [inf, sup[ to file
    void write_data(Outputter& out, lati_t inf, lati_t sup) const;
    
    /// write specified range to file
    void write(Outputter& out, lati_t inf, lati_t sup) const
    {
        if ( sup < inf )
            throw InvalidIO("incoherent Lattice boundaries");
        if ( laUnit < REAL_EPSILON )
            throw InvalidIO("incoherent Lattice unit value");

        out.writeInt32(inf);
        out.writeInt32(sup);
        out.writeFloat(laUnit);
        
        if ( laCell )
            write_data(out, inf, sup);
        else
        {
            out.writeUInt16(0);
            out.writeUInt8(0);
            out.writeUInt8(0);
        }
    }
    
    /// write all data to file
    void write(Outputter& out) const
    {
        write(out, laInf, laSup);
    }

    /// clear all cells and read data from file
    void read(Inputter& in)
    {
        lati_t inf = in.readInt32();
        lati_t sup = in.readInt32();
        real uni = in.readFloat();
        in.readUInt16();
        in.readUInt8();
        int nbytes = in.readUInt8();
        
        if ( sup < inf )
            throw InvalidIO("incoherent Lattice boundaries");
        if ( uni < REAL_EPSILON )
            throw InvalidIO("incoherent Lattice unit value");
        
        //only change unit length if the difference is significant
        if ( laUnit <= 0 )
            laUnit = uni;
        else if ( abs_real( laUnit - uni ) > 1e-6 )
            changeUnit(uni);
        
        allocate(inf, sup, 0);
        clear();
        
        //std::clog << "reading fiber:lattice " << nbytes << "bytes " << sup-inf << " cells\n";
        
        if ( nbytes == 1 )
        {
            for ( lati_t s = inf; s < sup; ++s )
                laCell[s] = in.readUInt8();
        }
        else if ( nbytes == 2 )
        {
            for ( lati_t s = inf; s < sup; ++s )
                laCell[s] = in.readUInt16();
        }
        else if ( nbytes == 4 )
        {
            for ( lati_t s = inf; s < sup; ++s )
                laCell[s] = in.readUInt32();
        }
        else if ( nbytes == 8 )
        {
            for ( lati_t s = inf; s < sup; ++s )
                laCell[s] = in.readUInt64();
        }
        else
        {
            for ( lati_t s = inf; s < sup; ++s )
                laCell[s] = in.readFloat();
        }
    }
    
    #pragma mark - doc

    /// quantification
    void info(size_t& cnt, size_t& vac, real& sum, real& mn, real& mx, bool density) const
    {
        const real scale = ( density ? 1.0/unit() : 1.0 );
        const auto sup = indexP();
        for ( auto i = indexM(); i <= sup; ++i )
        {
            ++cnt;
            if ( data(i) == 0 )
                ++vac;
            else
            {
                sum += data(i);
                real x = data(i) * scale;
                mn = std::min(mn, x);
                mx = std::max(mx, x);
            }
        }
    }
    
    /// printout
    void dump(FILE *f)
    {
        Outputter out(f, 0);
        write(out);
    }
    
    /// check basic requirements
    int bad()
    {
        if ( laUnit <= 0 )   return 1;
        if ( !laCell )       return 2;
        if ( laInf > laSup ) return 4;
        return 0;
    }
};

#endif
