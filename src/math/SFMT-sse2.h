#pragma once
/**
 * @file  SFMT-sse2.h
 * @brief SIMD oriented Fast Mersenne Twister(SFMT) for Intel SSE2
 *
 * @author Mutsuo Saito (Hiroshima University)
 * @author Makoto Matsumoto (Hiroshima University)
 *
 * @note We assume LITTLE ENDIAN in this file
 *
 * Copyright (C) 2006, 2007 Mutsuo Saito, Makoto Matsumoto and Hiroshima
 * University. All rights reserved.
 *
 * The new BSD License is applied to this software, see LICENSE.txt
 */

#ifndef SFMT_SSE2_H
#define SFMT_SSE2_H

/**
 * parameters used by sse2.
 */
static const w128_t sse2_param_mask = {{SFMT_MSK1, SFMT_MSK2, SFMT_MSK3, SFMT_MSK4}};

/**
 * This function represents the recursion formula.
 * @param r an output
 * @param a a 128-bit part of the interal state array
 * @param b a 128-bit part of the interal state array
 * @param c a 128-bit part of the interal state array
 * @param d a 128-bit part of the interal state array
 */
static inline __m128i mm_recursion(__m128i x, __m128i y, __m128i z, __m128i v)
{
    y = _mm_srli_epi32(y, SFMT_SR1);
    z = _mm_srli_si128(z, SFMT_SR2);
    v = _mm_slli_epi32(v, SFMT_SL1);
    z = _mm_xor_si128(z, x);
    z = _mm_xor_si128(z, v);
    x = _mm_slli_si128(x, SFMT_SL2);
    y = _mm_and_si128(y, sse2_param_mask.xmm);
    z = _mm_xor_si128(z, x);
    return _mm_xor_si128(z, y);
}

/**
 * This function fills the internal state array with pseudorandom
 * integers.
 * @param sfmt SFMT internal state
 */
void sfmt_gen_rand_all(sfmt_t * sfmt)
{
    int i;
    __m128i r1, r2, rr;
    __m128i * state = (__m128i*)sfmt->state;
    
    r1 = state[SFMT_N - 2];
    r2 = state[SFMT_N - 1];
    for (i = 0; i < SFMT_N - SFMT_POS1; ++i) {
        rr = mm_recursion(state[i], state[i + SFMT_POS1], r1, r2);
        state[i] = rr;
        r1 = r2;
        r2 = rr;
    }
    for (i = SFMT_N - SFMT_POS1; i < SFMT_N; ++i) {
        rr = mm_recursion(state[i], state[i + (SFMT_POS1-SFMT_N)], r1, r2);
        state[i] = rr;
        r1 = r2;
        r2 = rr;
    }
}

/**
 * This function fills the user-specified array with pseudorandom
 * integers.
 * @param sfmt SFMT internal state.
 * @param array an 128-bit array to be filled by pseudorandom numbers.
 * @param size number of 128-bit pseudorandom numbers to be generated.
 */
static void gen_rand_array(sfmt_t * sfmt, w128_t * buffer, int size)
{
    int i, j;
    __m128i r1, r2, rr;
    __m128i * state = (__m128i*)sfmt->state;
    __m128i * array = (__m128i*)buffer;
    
    r1 = state[SFMT_N - 2];
    r2 = state[SFMT_N - 1];
    for (i = 0; i < SFMT_N - SFMT_POS1; i++) {
        rr = mm_recursion(state[i], state[i + SFMT_POS1], r1, r2);
        array[i] = rr;
        r1 = r2;
        r2 = rr;
    }
    for (; i < SFMT_N; i++) {
        rr = mm_recursion(state[i], array[i - (SFMT_N-SFMT_POS1)], r1, r2);
        array[i] = rr;
        r1 = r2;
        r2 = rr;
    }
    for (; i < size - SFMT_N; i++) {
        rr = mm_recursion(array[i - SFMT_N], array[i - (SFMT_N-SFMT_POS1)], r1, r2);
        array[i] = rr;
        r1 = r2;
        r2 = rr;
    }
    for (j = 0; j < 2 * SFMT_N - size; j++) {
        state[j] = array[j + size - SFMT_N];
    }
    for (; i < size; i++, j++) {
        rr = mm_recursion(array[i - SFMT_N], array[i - (SFMT_N-SFMT_POS1)], r1, r2);
        array[i] = rr;
        r1 = r2;
        r2 = rr;
        state[j] = rr;
    }
}


#endif
