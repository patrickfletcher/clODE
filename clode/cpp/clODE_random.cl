#ifndef CLODE_RANDOM_H_
#define CLODE_RANDOM_H_

#include "realtype.cl"

// A private stream of PRNs will be generated, since the whole ODE solver is per work item.
// Seed can be created as global array. final state of RNG can be returned to global

//seeding the RNG
// seed with i=get_global_id[0]? or better to pass in more random seed?
//=>make a rand_init kernel: hash global_id to create RNGstate array.
// return final state back to a global array to continue stream if desired?

//TODO: not clear how to make host code aware of rngstatetype. Only use 64bit int methods?
// -> template rngData struct on realtype, rngstatetype. Pass N_RNGSTATE as compiler define to OpenCL? Store array of structs (not just rngData.state)

#define XORSHIRO128_PLUS
//~ #define XORSHIFT128_PLUS

//Wrappers with fixed function signatures. Conditional compilation could be used to select the backend "next" function

#ifdef XORSHIRO128_PLUS

//TODO: this can be set as compiler flag to enable different RNG algorithms
#define N_RNGSTATE (2)

typedef ulong rngstatetype;
#define RNGNORM RCONST(5.421010862427522e-20)
// #if defined(CLODE_SINGLE_PRECISION)
// 	#define RNGNORM (5.4210109e-20f) // 1.0/(2^64 - 1) single precision
// #elif defined(CLODE_DOUBLE_PRECISION)
// 	#define RNGNORM (5.421010862427522e-20) // 1.0/(2^64 - 1) double precision
// #endif

static inline ulong rotl(const ulong x, int k)
{
	return (x << k) | (x >> (64 - k));
}

ulong next(__private rngstatetype s[])
{
	const ulong s0 = s[0];
	ulong s1 = s[1];
	const ulong result = s0 + s1;

	s1 ^= s0;
	s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
	s[1] = rotl(s1, 36);				   // c

	return result;
}

#endif

#ifdef XORSHIFT128_PLUS

#define N_RNGSTATE (2)

typedef ulong rngstatetype;
#define RNGNORM RCONST(5.421010862427522e-20)

ulong next(ulong s[])
{
	ulong s1 = s[0];
	const ulong s0 = s[1];
	const ulong result = s0 + s1;
	s[0] = s0;
	s1 ^= s1 << 23;							 // a
	s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
	return result;
}

#endif

//hold the RNG's state and other data in a struct
typedef struct rngData
{
	rngstatetype state[N_RNGSTATE];
	bool randnUselast;
	realtype randnLast;
} rngData;

//return uniform pseudorandom number in [0,1)
inline realtype rand(__private rngstatetype *state)
{
	rngstatetype result = next(state);
	return result * RNGNORM;
};

//return normally distributed pseudorandom number N(0,1)
//polar method, generates two at a time requiring  external storage of useLast switch and y2.... ugly! but no access to work-item private static vars..
inline realtype randn(__private rngData *rd)
{

	realtype x1, x2, w, y1;

	if (rd->randnUselast)
	{
		y1 = rd->randnLast;
		rd->randnUselast = 0;
	}
	else
	{
		do
		{
			x1 = RCONST(2.0) * rand(rd->state) - RCONST(1.0);
			x2 = RCONST(2.0) * rand(rd->state) - RCONST(1.0);
			w = x1 * x1 + x2 * x2;
		} while (w >= RCONST(1.0));

		w = sqrt((-RCONST(2.0) * log(w)) / w);
		y1 = x1 * w;
		rd->randnLast = x2 * w;
		rd->randnUselast = 1;
	}

	return y1;
};

#endif //CLODE_RANDOM_H_
