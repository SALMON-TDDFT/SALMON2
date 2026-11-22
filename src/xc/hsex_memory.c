/* Process high-water mark; each diagnostic rank is launched separately. */
#include <sys/resource.h>
double hsex_peak_rss_mib(void) {
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) return -1.0;
#ifdef __APPLE__
    return (double)usage.ru_maxrss / (1024.0 * 1024.0);
#else
    return (double)usage.ru_maxrss / 1024.0;
#endif
}
