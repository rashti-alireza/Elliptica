#ifndef elliptica_system_LIB_H
#define elliptica_system_LIB_H

/* define handy macro for unsigned */
#define Uint unsigned

/* define inline calls.
// to make it portable, each inline function is defined in the common
// header file, thus one better to use static storage. also to avoid
// compiler warnings of unused functions when inline has not been defined,
// unused attribute is defined. */
#ifdef INLINE_FUNC
# define INLINE static inline
# define INLINE_WARN_UNUSED_FUNC 
#else
# define INLINE static
# define INLINE_WARN_UNUSED_FUNC __attribute__((unused))
#endif

#endif

