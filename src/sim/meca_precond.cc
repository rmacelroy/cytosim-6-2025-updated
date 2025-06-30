// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "xtbsv.h"
#include "xtrsm.h"

/**
 machine_time() provides a timer used to evaluate the Preconditionner's performance
This can be disabled if you are not using auto-selection of preconditionners
 */
#if ( 1 )
static inline unsigned long machine_time()
{
    return 0;
}
#elif defined(__APPLE__)
#  include <mach/mach_time.h>
static inline unsigned long machine_time()
{
    // the units of this counter is not specified, but seems close to nanosec
    return mach_absolute_time() >> 10;
}
#else
#  include <sys/time.h>
static inline unsigned long machine_time()
{
    timespec tv;
    clock_gettime(CLOCK_MONOTONIC, &tv);
    return 1000 * (unsigned long)tv.tv_sec + tv.tv_nsec / 1000;
}
#endif

//------------------------------------------------------------------------------

/// if you like Alsatian specialities, set this to a prime number
#define CHOUCROUTE 7
#define SAUERKRAUT 7

/// leading dimension of the banded matrix used for iso symmetric blocks
constexpr index_t ISOB_KD = 2;
constexpr index_t ISOB_LDD = 3;

/*
 number of off-diagonals and leading dimension for the non-isotropic banded matrix
 We use a band preconditionner with 2*DIM off-diagonals to include the near-
 diagonal terms from the blocks that are offset by 2 from the matrix diagonal
 */
constexpr index_t BAND_NUD = 6;

// should allocate to also hold the true diagonal: BAND_LDD > BAND_NUD
constexpr index_t BAND_LDD = BAND_NUD+2;


template < typename REAL >
static void PRINT_MAT(std::string const& msg, index_t lin, index_t col, const REAL* mat, index_t ldd, int digits=0)
{
    VecPrint::full(msg, lin, col, mat, ldd, digits);
}

//------------------------------------------------------------------------------
#pragma mark - Apply Blocks

/// apply banded symmetric isotropic preconditionner block
static inline void applyPrecondIsoB(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();

#if CHOUCROUTE
    assert_true( ISOB_KD == 2 );
    alsatian_iso_xpbtrsL<DIM>(nbp, mec->pblock(), ISOB_LDD, Y);
#else
    /*
     we cannot call lapack::DPBTRS('L', bks, KD, 1, mec->pblock(), KD+1, Y, bks, &info)
     because the coordinates of the vector 'Y' are not contiguous but offset by 'DIM'.
     But calling DTBSV gets the required work done.
     */
    for ( int d = 0; d < DIM; ++d )
    {
        blas::xtbsv('L', 'N', 'N', nbp, ISOB_KD, mec->pblock(), ISOB_LDD, Y+d, DIM);
        blas::xtbsv('L', 'T', 'N', nbp, ISOB_KD, mec->pblock(), ISOB_LDD, Y+d, DIM);
    }
#endif
}


/// apply symmetric isotropic preconditionner block
static inline void applyPrecondIsoS(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();
    /*
     we cannot call lapack::DPOTRS('L', nbp, mec->pblock(), nbp, Y, DIM, &info);
     because the coordinates of the vector 'Y' are not contiguous but offset by 'DIM'.
     */
#if 0
    real * tmp = new_real(DIM*nbp);
    copy_real(DIM*nbp, Y, tmp);
    iso_xpotrsL<DIM>(nbp, mec->pblock(), nbp, tmp);
    std::clog << "\n "; VecPrint::print(nbp, tmp, 3, 100.0);
    free_real(tmp);
#endif
#if CHOUCROUTE
    alsatian_iso_xpotrsL<DIM>(nbp, mec->pblock(), nbp, Y);
#else
    iso_xpotrsL<DIM>(nbp, mec->pblock(), nbp, Y);
#endif
    //std::clog << "\nL"; VecPrint::print(nbp, Y, 3, 100.0);
}


/// apply non-symmetric but isotropic preconditionner block
static inline void applyPrecondIsoP(Mecable const* mec, real* Y)
{
    int nbp = mec->nbPoints();
#if SAUERKRAUT
    alsatian_iso_xgetrsN<DIM>(nbp, mec->pblock(), nbp, mec->pivot(), Y);
#else
    iso_xgetrsN<DIM>(nbp, mec->pblock(), nbp, mec->pivot(), Y);
#endif
}


/// apply banded symmetric preconditionner block, double precision
static inline void applyPrecondBand(Mecable const* mec, real* Y)
{
    const int bks = DIM * mec->nbPoints();
    assert_true( (int)BAND_NUD < bks );
#if SAUERKRAUT && REAL_IS_DOUBLE && USE_SIMD
    static_assert(BAND_NUD==6, "BAND_NUD should be 6");
    alsatian_xtbsvLNN6K_SSE(bks, (float*)mec->pblock(), BAND_LDD, Y);
    alsatian_xtbsvLTN6K_SSE(bks, (float*)mec->pblock(), BAND_LDD, Y);
#elif SAUERKRAUT
    alsatian_xpbtrsLK<BAND_NUD>(bks, mec->pblock(), BAND_LDD, Y);
#elif 1
    blas_xtbsvLN<'N'>(bks, BAND_NUD, mec->pblock(), BAND_LDD, Y, 1);
    blas_xtbsvLT<'N'>(bks, BAND_NUD, mec->pblock(), BAND_LDD, Y, 1);
#elif 1
    blas::xtbsv('L', 'N', 'N', bks, BAND_NUD, mec->pblock(), BAND_LDD, Y, 1);
    blas::xtbsv('L', 'T', 'N', bks, BAND_NUD, mec->pblock(), BAND_LDD, Y, 1);
#else
    int info = 0;
    lapack::xpbtrs('L', bks, BAND_NUD, 1, mec->pblock(), BAND_LDD, Y, bks, &info);
    assert_true(info==0);
#endif
}


/// apply full size symmetric preconditionner block, single precision
static inline void applyPrecondHalf(Mecable const* mec, real* Y)
{
    const int bks = DIM * mec->nbPoints();
#if SAUERKRAUT
    // assuming that diagonal terms of the preconditionner block have been inverted:
    alsatian_xpotrsL(bks, (float*)mec->pblock(), bks, Y);
#elif 1
    iso_xpotrsL<1>(bks, mec->pblock(), bks, Y);
#else
    int info = 0;
    lapack::xpotrs('L', bks, 1, mec->pblock(), bks, Y, bks, &info);
    assert_true(info==0);
#endif
}


/// apply full size non-symmetric preconditionner block, single precision
static inline void applyPrecondFull(Mecable const* mec, real* Y)
{
    const int bks = DIM * mec->nbPoints();
#if CHOUCROUTE && USE_SIMD
    // assuming that diagonal terms of the preconditionner block have been inverted:
    alsatian_xgetrsN_SSE(bks, (float*)mec->pblock(), bks, mec->pivot(), Y);
#elif CHOUCROUTE
    alsatian_xgetrsN(bks, (float*)mec->pblock(), bks, mec->pivot(), Y);
#elif 1
    // translated LAPACK's reference code:
    lapack_xgetrsN(bks, mec->pblock(), bks, mec->pivot(), Y);
#else
    // using LAPACK's library
    int info = 0;
    lapack::xgetrs('N', bks, 1, mec->pblock(), bks, mec->pivot(), Y, bks, &info);
    assert_true(info==0);
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Apply Preconditionner


/// apply preconditionner block corresponding to Mecable
static inline void applyPreconditionner(Mecable const* mec, real* Y)
{
    switch ( mec->blockType() )
    {
        case 0: break;
        case 1: applyPrecondIsoB(mec, Y); break;
        case 2: applyPrecondIsoS(mec, Y); break;
        case 3: applyPrecondIsoP(mec, Y); break;
        case 4: applyPrecondBand(mec, Y); break;
        case 5: applyPrecondHalf(mec, Y); break;
        case 6: applyPrecondFull(mec, Y); break;
        default: ABORT_NOW("unknown Mecable::blockType()"); break;
    }
}

/// apply preconditionner to entire system: Y <- Preconditionner * X
/** This can be done in parallel if the preconditionner is block-diagonal */
void Meca::precondition(const real* X, real* Y) const
{
    unsigned long rdt = machine_time();
    if ( Y != X )
        copy_real(dimension(), X, Y);
    
    #pragma omp parallel for
    for ( Mecable const* mec : mecables )
    {
        const index_t inx = DIM * mec->matIndex();
            applyPreconditionner(mec, Y+inx);
    }
    // CPU time is recorded, but only for info, so this code can be removed
    cycles_ += machine_time() - rdt;
}


/// total allocated memory size for preconditionner
size_t Meca::preconditionnerSize() const
{
    size_t res = 0;
    for ( Mecable const* mec : mecables )
        res += mec->blockLimit();
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Debug Preconditionner

[[maybe_unused]]
static void printPreconditionnerBlock(Mecable* mec, index_t sup)
{
    const index_t bks = DIM * mec->nbPoints();
    index_t S = std::min(bks, sup);
    real * blk = mec->pblock();

    PRINT_MAT("Diagonal block"+std::to_string(bks), S, S, mec->pblock(), bks);
    
    std::clog << "S ";
    char str[32];
    for ( index_t i = 0; i < S; ++i )
    {
        real sum = 0;
        for ( index_t j = 0; j < bks; ++j )
            sum += blk[j+bks*i];
        snprintf(str, sizeof(str), " %8.3f", sum);
        std::clog << str;
    }
    std::clog << "\n";
}


/**
 This version builds the diagonal block indirectly using Meca:multiply().
 This is a slow method that calls 'multiply()' n-times, where
 'n' is the size of the block.
 
 This should be used for validation only.
*/
void Meca::extractBlock(const Mecable* mec, real* res) const
{
    const index_t dim = dimension();
    const index_t bks = DIM * mec->nbPoints();
    const index_t off = DIM * mec->matIndex();
    
    assert_true( off+bks <= dim );
    real * vec = new_real(dim);
    real * tmp = new_real(dim);
    
    zero_real(dim, vec);
    //zero_real(bks*bks, res);
    
    // proceed column by column:
    for ( index_t jj = 0; jj < bks; ++jj )
    {
        vec[jj+off] = 1;
        multiply(vec, tmp);
        vec[jj+off] = 0;
        copy_real(bks, tmp+off, res+jj*bks);
    }
    
    free_real(vec);
    free_real(tmp);
}


/**
 DEBUG: compare `blk` with block extracted using extractBlock()
 */
void Meca::verifyBlock(const Mecable * mec, const real* blk)
{
    const index_t bks = DIM * mec->nbPoints();
    real * wrk = new_real(bks*bks);
    
    extractBlock(mec, wrk);
    
    real mag = blas::nrm2(bks*bks, blk);
    blas::xaxpy(bks*bks, -1.0, blk, 1, wrk, 1);
    real err = blas::nrm2(bks*bks, wrk);
 
    std::clog << "verifyBlock ";
    std::clog << std::setw(8) << mec->reference() << "  " << std::setw(4) << bks;
    std::clog << " | B - K | = " << err << ' ';

    if ( err > mag * 1e-6 )
    {
        //VecPrint::sparse(std::clog, bks, bks, wrk, bks, 3, (real)0.1);
        extractBlock(mec, wrk);
        PRINT_MAT("\nextracted", bks, bks, wrk, bks);
        PRINT_MAT("computed", bks, bks, blk, bks);
    }
    
    free_real(wrk);
}


/**
 Multiply here `blk` with the dynamic block extracted by extractBlock()
 and check that we recover the identity matrix
 */
void Meca::checkBlock(const Mecable * mec, const real* blk)
{
    const index_t bks = DIM * mec->nbPoints();
    
    if ( 0 == mec->blockType() )
        return;

    std::clog << "  checkBlock " << mec->blockType() << " ";
    std::clog << std::setw(10) << mec->reference() << " " << std::setw(6) << bks;
    
    real * wrk = new_real(bks*bks);
    real * mat = new_real(bks*bks);
    real * vec = new_real(bks);
    
    extractBlock(mec, wrk);   // wrk <- MEC_BLOCK
   
    copy_real(bks*bks, wrk, mat);  // mat <- wrk
    for ( index_t i = 0; i < bks; ++i )
    {
        applyPreconditionner(mec, mat+bks*i);
        mat[i+bks*i] -= 1;
    }
    real err = blas::nrm2(bks*bks, mat) / bks;
    std::clog << " | 1 - PM | = " << err;
    
    if ( 1 == mec->blockType() )
    {
        // use power iterations to estimate largest eigenvalue
        blas::xcopy(bks, vRHS+DIM*mec->matIndex(), 1, vec, 1);
        real eig = largest_eigenvalue(bks, blk, mec->pivot(), wrk, -1.0, vec, mat);
        std::clog << "  eigen(1-PM) = " << eig;
    }
    
    std::clog << '\n';

    if ( err > 1 )
    {
        std::cout << "\n";
        // print preconditionner block for visual inspection:
        PRINT_MAT("matrix", bks, bks, wrk, bks);
        PRINT_MAT("precond", bks, bks, blk, bks);
        PRINT_MAT("precond * matrix", bks, bks, mat, bks);
    }
    free_real(vec);
    free_real(mat);
    free_real(wrk);
}


//------------------------------------------------------------------------------
#pragma mark - Extract Diagonal Block

/**
 Build a banded block that can be used as a preconditionner:
 
     I - time_step * mobility * Rigidity
 
 This block is square, symmetric, definite positive and well-behaved
 Banded storage is used and, mat(i, j) is stored in blk[i-j+ldd*j]
 */
void Meca::getIsoBandedBlock(const Mecable * mec, real* res, index_t kd, index_t ldd) const
{
    const index_t nbp = mec->nbPoints();

    const real beta = -tau_ * mec->pointMobility();
    real jR = mec->jointRigidity();
    if ( jR != 0 )
    {
        if ( ldd != 3 )
            zero_real(ldd*nbp, res);
        setBendingRigidity<1>(res, ldd-1, nbp, beta*jR);
    }
    else
        zero_real(ldd*nbp, res);
    
    //PRINT_MAT("bending elasticity band", 3, nbp, res, ldd);

    /*
     The matrix `res` is stored in 'packed symmetric banded storage':
     usually, mat(i, j) is stored in mat[i+ldd*j] but with banded storage,
        mat(i, j) is stored in mat[i-j+ldd*j] for i > j
     We can get the correct addressing using `ldd-1` instead of `ldd`
     */

#if USE_ISO_MATRIX
    mISO.addLowerBand(beta, res, ldd-1, mec->matIndex(), nbp, 1, kd);
    if ( useFullMatrix )
#endif
    {
        mFUL.addDiagonalTrace(beta/DIM, res, ldd-1, mec->matIndex(), nbp, DIM, kd, false);
    }
    
    // add Identity matrix to band storage:
    for ( index_t i = 0; i < nbp; ++i )
        res[ldd*i] += 1;
}


/**
 Extract block of reduced dimension that does not include projection:
 
     I - time_step * point_mobility ( mISO + mFUL )
 
 The result is constructed by using efficient methods from mISO and mFUL
 The block is symmetric and can be factorized by Cholesky, which may fail
 Normal storage is used and, mat(i, j) is stored in blk[i+nbp*j]
 */
void Meca::getIsoBlock(const Mecable * mec, real* res) const
{
    const index_t nbp = mec->nbPoints();
    
    zero_real(nbp*nbp, res);
    
#if SEPARATE_RIGIDITY_TERMS
    // set the Rigidity terms:
    real jR = mec->jointRigidity();
    if ( jR != 0 )
        //addBendingRigidityLower<1>(res, nbp, mec->nbPoints(), jR);
        addBendingRigidity<1>(res, nbp, mec->nbPoints(), jR);
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, nbp, mec->matIndex(), nbp, 1);
    if ( useFullMatrix )
#endif
    {
        mFUL.addDiagonalTrace(1.0/DIM, res, nbp, mec->matIndex(), nbp, DIM, nbp, true);
    }
    
    // the projection is not called, so we scale by mobility
    const real beta = -tau_ * mec->pointMobility();

    //blas::xscal(bs*bs, beta, res, 1);
    for ( index_t n = 0; n < nbp*nbp; ++n )
        res[n] = beta * res[n];
    
    // add Identity matrix:
    for ( index_t i = 0; i < nbp*nbp; i += nbp+1 )
        res[i] += 1;
}


/**
 Get a diagonal block corresponding to an Object, which is:
 
     I - time_step * mob * ( mISO + mFUL )
 
 The result is constructed by using functions from mISO and mFUL, and then
 multiplied by the vertex mobility, but projection is not applied.
 This block is banded and symmetric, and can be factorized by Cholesky's method!
 */
void Meca::getBandedBlock(const Mecable * mec, real* res, index_t ldd, index_t rank) const
{
    const index_t nbp = mec->nbPoints();
    const index_t bks = DIM * nbp;

    zero_real(ldd*bks, res);

    // multiply by mobility coefficient, skipping projection
    const real beta = -tau_ * mec->pointMobility();
    
#if SEPARATE_RIGIDITY_TERMS
    // set the Rigidity terms:
    real jR = mec->jointRigidity();
    if ( jR != 0 )
    {
        setBendingRigidity<DIM>(res, ldd-1, nbp, beta*jR);
        //PRINT_MAT("Rigidity block ", bks, bks, res, bks);
    }
#endif
#if USE_ISO_MATRIX
    mISO.addLowerBand(beta, res, ldd-1, mec->matIndex(), nbp, 1, rank/DIM);
#endif
    // PRINT_MAT("\niso", ldd, bks, res, ldd);
    copy_lower_subspace<DIM>(bks, res, ldd-1, rank);
    // PRINT_MAT("\ncopy_subspace", ldd, bks, res, ldd);
#if USE_ISO_MATRIX
    if ( useFullMatrix )
#endif
        mFUL.addLowerBand(beta, res, ldd-1, mec->matIndex(), nbp, DIM, rank);
    
    // add Identity matrix to band storage:
    for ( index_t i = 0; i < bks; ++i )
        res[ldd*i] += 1;
}

/**
 Get a diagonal block corresponding to an Object, which is:
 
     I - time_step * mob * ( mISO + mFUL )
 
 The result is constructed by using functions from mISO and mFUL, and then
 multiplied by the vertex mobility, to approximate the dynamics.
 This block is square and symmetric, and can be factorized by Cholesky's method!
 */
void Meca::getHalfBlock(const Mecable * mec, real* res) const
{
    const index_t nbp = mec->nbPoints();
    const index_t bks = DIM * nbp;
    
    zero_real(bks*bks, res);
    
#if SEPARATE_RIGIDITY_TERMS
    // set the Rigidity terms:
    real jR = mec->jointRigidity();
    if ( jR != 0 )
    {
        addBendingRigidityLower<DIM>(res, bks, mec->nbPoints(), jR);
        //PRINT_MAT("Rigidity block", bks, bks, res, bks);
    }
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, bks, mec->matIndex(), nbp, DIM);
#endif
    copy_lower_subspace<DIM, true>(bks, res, bks);
#if USE_ISO_MATRIX
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalBlock(res, bks, mec->matIndex(), nbp, DIM);
    
    // multiply by mobility coefficient, skipping projection
    const real beta = -tau_ * mec->pointMobility();
    
    //blas::xscal(bs*bs, beta, res, 1);
    for ( index_t n = 0; n < bks*bks; ++n )
        res[n] = beta * res[n];
    
    // add Identity matrix:
    for ( index_t i = 0; i < bks*bks; i += bks+1 )
        res[i] += 1;
}


/**
 Get the total diagonal block corresponding to an Object, which is:
 
     I - time_step * P * ( mISO + mFUL + diffP )
 
 The result is constructed by using functions from mISO and mFUL
 This block is square but not symmetric!
 */
void Meca::getFullBlock(const Mecable * mec, real* res) const
{
    const index_t nbp = mec->nbPoints();
    const index_t bks = DIM * nbp;
    
    zero_real(bks*bks, res);
    
#if SEPARATE_RIGIDITY_TERMS
    // set the Rigidity terms:
    real jR = mec->jointRigidity();
    if ( jR != 0 )
    {
        addBendingRigidityLower<DIM>(res, bks, mec->nbPoints(), jR);
        //PRINT_MAT("Rigidity block", bks, bks, res, bks);
    }
#endif
#if USE_ISO_MATRIX
    mISO.addDiagonalBlock(res, bks, mec->matIndex(), nbp, DIM);
#endif
    copy_lower_subspace<DIM, true>(bks, res, bks);
#if USE_ISO_MATRIX
    if ( useFullMatrix )
#endif
        mFUL.addDiagonalBlock(res, bks, mec->matIndex(), nbp, DIM);
    
    //PRINT_MAT("mISO+mFUL block", bks, bks, res, bks);
    
#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
        mec->addProjectionDiff(res);
#endif

    // include the projection, by applying it to each column vector:
    /* This could be vectorized */
    for ( index_t i = 0; i < bks; ++i )
        mec->projectForces(res+bks*i, res+bks*i);

    // scale
    const real beta = -tau_ * mec->leftoverMobility();
    //blas::xscal(bs*bs, beta, res, 1);
    for ( index_t n = 0; n < bks*bks; ++n )
        res[n] = beta * res[n];
    
    // add Identity matrix:
    for ( index_t i = 0; i < bks*bks; i += bks+1 )
        res[i] += 1;
}

//------------------------------------------------------------------------------
#pragma mark - Compute Preconditionner Block

/**
 Compute a preconditionner block corresponding to 'mec'
 The dimension is reduced by DIM and banded with diagonal + 2 off-diagonals
 This block is usually symmetric definite positive, and is factorized by Cholesky's method
 */
void Meca::computePrecondIsoB(Mecable* mec, real* tmp)
{
    const index_t nbp = mec->nbPoints();
    mec->blockSize(DIM*nbp, ISOB_LDD*nbp, 0);
    
    /**
     Factorize banded matrix with Andre-Louis Cholesky's method
     born 15.10.1875 in Montguyon, France
     died 31.08.1918 in Bagneux, following wounds received in battle.
     */
    int info = 0;
    unsigned short bt = 0;
    if ( ISOB_LDD <= nbp )
    {
        getIsoBandedBlock(mec, mec->pblock(), ISOB_KD, ISOB_LDD);
#if 0
        PRINT_MAT("isoBand", ISOB_LDD, nbp, mec->pblock(), ISOB_LDD);
        // visual comparison with getIsoBlock()
        getFullBlock(mec, tmp);
        PRINT_MAT("fullBlock", DIM*nbp, DIM*nbp, tmp, DIM*nbp);
        //getIsoBlock(mec, tmp);
        //PRINT_MAT("isoBlock", nbp, nbp, tmp, nbp);
#endif
        // calculate Banded Cholesky factorization:
#if CHOUCROUTE
        alsatian_xpbtf2L(nbp, ISOB_KD, mec->pblock(), ISOB_LDD, &info);
#else
        lapack::xpbtf2('L', nbp, ISOB_KD, mec->pblock(), ISOB_LDD, &info);
#endif
        bt = 1;
    }
    else
    {
        getIsoBlock(mec, mec->pblock());
#if 0
        PRINT_MAT("isoBlock", nbp, nbp, mec->pblock(), nbp);
        getFullBlock(mec, tmp);
        PRINT_MAT("fullBlock", DIM*nbp, DIM*nbp, tmp, DIM*nbp);
#endif
        // calculate Cholesky factorization:
#if CHOUCROUTE
        alsatian_xpotf2L(nbp, mec->pblock(), nbp, &info);
#else
        lapack::xpotf2('L', nbp, mec->pblock(), nbp, &info);
#endif
        bt = 2;
    }

    if ( 0 == info )
    {
        mec->blockType(bt);
        //PRINT_MAT("factorized", 3, nbp, mec->pblock(), ISOB_LDD);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute Band Preconditionner block of size " << nbp << "\n";
        ++bump_;
        //computePrecondIsoS(mec);
    }
}


/**
 Compute a preconditionner block corresponding to 'mec':
 Block of dimension reduced by DIM and without projection
 This block is symmetric definite positive, and is factorized by Cholesky's method
 */
void Meca::computePrecondIsoS(Mecable* mec)
{
    const index_t nbp = mec->nbPoints();
#if 0
    const index_t bks = DIM * nbp;
    mec->blockSize(bks, bks*bks, bks);
    getFullBlock(mec, mec->pblock());
    project_matrix<DIM>(nbp, mec->pblock(), bks, mec->pblock(), nbp);
    PRINT_MAT("projected: ", nbp, nbp, mec->pblock(), nbp);
#endif

    mec->blockSize(DIM*nbp, nbp*nbp, 0);
    
    //getIsoBandedBlock(mec, mec->pblock(), ISOB_LDD);
    //PRINT_MAT("banded preconditionner ", 3, nbp, mec->pblock(), ISOB_LDD);

    getIsoBlock(mec, mec->pblock());
    
    //PRINT_MAT("isoSymm", nbp, nbp, mec->pblock(), nbp);

    // calculate Cholesky factorization:
    int info = 0;
#if CHOUCROUTE
    alsatian_xpotf2L(nbp, mec->pblock(), nbp, &info);
#else
    lapack::xpotf2('L', nbp, mec->pblock(), nbp, &info);
#endif
    
    if ( 0 == info )
    {
        mec->blockType(2);
        //checkBlock(mec, blk);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute IsoS Preconditionner block of size " << nbp << "\n";
        ++bump_;
    }
}

/**
Compute a preconditionner block corresponding to 'mec':
 Block of dimension reduced by DIM, including the projection
 The block is not symmetric and is factorized by LU decomposition
 */
void Meca::computePrecondIsoP(Mecable* mec, real* tmp)
{
    const index_t nbp = mec->nbPoints();
    int info = 0;

    const index_t bks = DIM * nbp;
    mec->blockSize(nbp, nbp*nbp, bks);
    real * blk = mec->pblock();

    getFullBlock(mec, tmp);
    // project matrix to reduce dimensionality:
    project_matrix<DIM>(nbp, tmp, bks, blk, nbp);

    //PRINT_MAT("BlockIsoP", nbp, nbp, blk, nbp);

    // calculate LU factorization:
#if SAUERKRAUT
    alsatian_xgetf2(nbp, blk, nbp, mec->pivot(), &info);
#else
    lapack::xgetf2(nbp, nbp, blk, nbp, mec->pivot(), &info);
#endif

    if ( 0 == info )
    {
        mec->blockType(3);
        //checkBlock(mec, blk);
#if SAUERKRAUT && REAL_IS_DOUBLE
        convert_to_floats(nbp*nbp, blk, (float*)blk);
#endif
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute IsoP Preconditionner block of size " << nbp << "\n";
        ++bump_;
    }
}

/**
 Compute banded symmetric preconditionner block corresponding to 'mec',
 factorized by Cholesky's method.
 */
void Meca::computePrecondBand(Mecable* mec, real* tmp)
{
    assert_true(BAND_NUD < BAND_LDD);
    const index_t bks = DIM * mec->nbPoints();
    const index_t lin = std::min(BAND_LDD, bks);
#if SAUERKRAUT && REAL_IS_DOUBLE
    mec->blockSize(bks, 4+bks*lin/2, 0);
    // use temporary memory to build matrix block:
    real * blk = tmp;
#else
    mec->blockSize(bks, bks*lin, 0);
    real * blk = mec->pblock();
#endif

    int info = 0;
    unsigned short bt = 0;

    if ( BAND_LDD < bks )
    {
        getBandedBlock(mec, blk, BAND_LDD, BAND_NUD);
#if ( 0 )
        PRINT_MAT("fullBand", BAND_NUD+1, bks, blk, BAND_LDD);
        getHalfBlock(mec, tmp);
        PRINT_MAT("halfBlock", bks, bks, tmp, bks);
        getIsoBandedBlock(mec, tmp, ISOB_KD, ISOB_LDD);
        PRINT_MAT("isoBand", ISOB_LDD, mec->nbPoints(), tmp, ISOB_LDD);
        getBandedBlock(mec, blk, BAND_LDD, BAND_NUD);
#endif
        // calculate Cholesky factorization for band storage:
#if SAUERKRAUT
        alsatian_xpbtf2L(bks, BAND_NUD, blk, BAND_LDD, &info);
        //alsatian_xpbtf2L_lapack(bks, BAND_NUD, blk, BAND_LDD, &info);
        //PRINT_MAT("factorizedBand", BAND_NUD+1, bks, blk, BAND_LDD, 2);
        //modify_matrix<DIM>(BAND_NUD+1, bks, blk, BAND_LDD);
        //PRINT_MAT("factorizedBand", BAND_NUD+1, bks, blk, BAND_LDD, 2);
#else
        lapack::xpbtf2('L', bks, BAND_NUD, blk, BAND_LDD, &info);
#endif
        bt = 4;
    }
    else
    {
        getHalfBlock(mec, blk);
        // calculate Cholesky factorization:
#if SAUERKRAUT
        alsatian_xpotf2L(bks, blk, bks, &info);
#else
        lapack::xpotf2('L', bks, blk, bks, &info);
#endif
        bt = 5;
    }
    
    if ( 0 == info )
    {
#if SAUERKRAUT && REAL_IS_DOUBLE
        convert_to_floats(bks*lin, blk, (float*)mec->pblock());
#endif
        mec->blockType(bt);
        //checkBlock(mec, blk);
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute band Preconditionner block of size " << bks << "\n";
        ++bump_;
    }
}

/**
 Compute preconditionner block corresponding to 'mec'
 This block is symmetric, and factorized by Cholesky's method.
 */
void Meca::computePrecondHalf(Mecable* mec, real* tmp)
{
    const index_t bks = DIM * mec->nbPoints();
#if SAUERKRAUT && REAL_IS_DOUBLE
    mec->blockSize(bks, 4+bks*bks/2, bks);
    // use temporary memory to build matrix block:
    real * blk = tmp;
#else
    mec->blockSize(bks, bks*bks, 0);
    real * blk = mec->pblock();
#endif
    getHalfBlock(mec, blk);
    
    //PRINT_MAT("halfBlock", bks, bks, mec->pblock(), bks);

    int info = 0;
    // calculate Cholesky factorization:
#if SAUERKRAUT
    alsatian_xpotf2L(bks, blk, bks, &info);
#else
    lapack::xpotf2('L', bks, blk, bks, &info);
#endif

    if ( 0 == info )
    {
        mec->blockType(5);
        //checkBlock(mec, blk);
        //PRINT_MAT("half", bks, bks, blk, bks);
#if SAUERKRAUT && REAL_IS_DOUBLE
        convert_to_floats(bks*bks, blk, (float*)mec->pblock());
#endif
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute half Preconditionner bloc of size " << bks << "\n";
        ++bump_;
    }
}

/**
Compute preconditionner block corresponding to 'mec'
 */
void Meca::computePrecondFull(Mecable* mec, real* tmp)
{
    const index_t bks = DIM * mec->nbPoints();
    
#if CHOUCROUTE && REAL_IS_DOUBLE
    mec->blockSize(bks, 4+bks*bks/2, bks);
    real * blk = tmp;
#else
    mec->blockSize(bks, bks*bks, bks);
    real * blk = mec->pblock();
#endif
    
    getFullBlock(mec, blk);
    //verifyBlock(mec, blk);
    //PRINT_MAT("fullBlock", bks, bks, blk, bks);
    
    // calculate LU factorization:
    int info = 0;
#if CHOUCROUTE
    alsatian_xgetf2(bks, blk, bks, mec->pivot(), &info);
#else
    lapack::xgetf2(bks, bks, blk, bks, mec->pivot(), &info);
#endif
    
    if ( 0 == info )
    {
        mec->blockType(6);
        //checkBlock(mec, blk);
        //if ( bks < 4 ) PRINT_MAT("Full", bks, bks, blk, bks);
#if CHOUCROUTE && REAL_IS_DOUBLE
        convert_to_floats(bks*bks, blk, (float*)mec->pblock());
#endif
    }
    else
    {
        mec->blockType(0);
        //std::clog << "failed to compute full Preconditionner block of size " << bks << "\n";
        ++bump_;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Compute Preconditionner

/*
 Here the preconditionner blocks are calculated,
 according to Simul::precondition
 */
void Meca::computePreconditionner()
{
    bump_ = 0;
    index_t sup = ( 1 + DIM * largestMecable() ) & ~1;
    real * tmp = new_real(4+sup*sup);

    switch( precond_ )
    {
        case 0:
            for ( Mecable * mec : mecables )
                mec->blockType(0);
            break;
        case 1:
            for ( Mecable * mec : mecables )
            {
                if ( mec->tag() == Solid::TAG )
                    computePrecondFull(mec, tmp);
                else
                    computePrecondIsoB(mec, tmp);
            }
            break;
        case 2:
            for ( Mecable * mec : mecables )
                computePrecondIsoS(mec);
            break;
        case 3:
            for ( Mecable * mec : mecables )
                computePrecondIsoP(mec, tmp);
            break;
        case 4:
            for ( Mecable * mec : mecables )
                computePrecondBand(mec, tmp);
            break;
        case 5:
            for ( Mecable * mec : mecables )
                computePrecondHalf(mec, tmp);
            break;
        case 6:
            for ( Mecable * mec : mecables )
                computePrecondFull(mec, tmp);
            break;
        case 7:
            for ( Mecable * mec : mecables )
            {
                if ( mec->tag() == Bead::TAG )
                    computePrecondIsoS(mec);
                else
                    computePrecondFull(mec, tmp);
            }
            break;
        default:
            throw InvalidParameter("unknown `precondition' value");
            break;
    }
    free_real(tmp);
}

