#ifdef TCC_TARGET_X86_64
#define TCC_LIBTCC1 "libtcc1-64.a"
#else
#define TCC_LIBTCC1 "libtcc1-32.a"
#endif
#ifdef __GNUC__
#define GCC_MAJOR __GNUC___
#define GCC_MINOR __GNUC_MINOR__
#endif
#define TCC_VERSION "0.9.27"
