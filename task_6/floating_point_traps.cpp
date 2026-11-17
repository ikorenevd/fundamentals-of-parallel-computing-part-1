#include "./floating_point_traps.h"

#ifdef ENABLE_FPE
#include <cfenv>

bool configure_floating_point_traps()
{
    const int exceptions = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW;
    if (std::feclearexcept(exceptions) != 0)
        return false;
    #if defined(__linux__)
        return ::feenableexcept(exceptions) != -1;
    #elif defined(__APPLE__)
        std::fenv_t env;
        if (std::fegetenv(&env) != 0)
            return false;
        #if defined(__aarch64__) || defined(__arm64__)
            // ARM FPCR exception-enable bits are shifted by 8 from FE_* flags.
            const auto traps = static_cast<decltype(env.__fpcr)>(exceptions) << 8;
            env.__fpcr |= traps;
            if (std::fesetenv(&env) != 0 || std::fegetenv(&env) != 0)
                return false;
            return (env.__fpcr & traps) == traps;
        #elif defined(__x86_64__) || defined(__i386__)
            // Unmask exceptions in both the x87 and SSE control registers.
            env.__control = static_cast<decltype(env.__control)>(env.__control & ~exceptions);
            env.__mxcsr &= ~(static_cast<unsigned int>(exceptions) << 7);
            if (std::fesetenv(&env) != 0 || std::fegetenv(&env) != 0)
                return false;
            return (env.__control & exceptions) == 0
                && (env.__mxcsr & (static_cast<unsigned int>(exceptions) << 7)) == 0;
        #else
            #error Unsupported macOS architecture for floating-point traps
        #endif
    #else
        #error Floating-point traps require Linux or macOS
    #endif
}
#else
bool configure_floating_point_traps()
{
    return true;
}
#endif
