#ifndef NETWORKIT_AUXILIARY_PERF_EVENT_HPP_
#define NETWORKIT_AUXILIARY_PERF_EVENT_HPP_

#include <cstdint>

namespace Aux {
class PerfEventCountHardware {
private:
    int fd;

public:
    #ifdef __linux__
    static constexpr bool is_available = true;
    #else
    static constexpr bool is_available = false;
    #endif

    enum Event : uint64_t {
        CPU_CYCLES = 0,
        INSTRUCTIONS = 1,
        CACHE_REFERENCES = 2,
        CACHE_MISSES = 3,
        BRANCH_INSTRUCTIONS = 4,
        BRANCH_MISSES = 5,
        BUS_CYCLES = 6,
        STALLED_CYCLES_FRONTEND = 7,
        STALLED_CYCLES_BACKEND = 8,
        REF_CPU_CYCLES = 9
    };

    PerfEventCountHardware(Event event);
    PerfEventCountHardware(const PerfEventCountHardware &) = delete;
    PerfEventCountHardware(PerfEventCountHardware &&) noexcept;
    PerfEventCountHardware *operator=(const PerfEventCountHardware &) = delete;
    ~PerfEventCountHardware();

    void enable();
    void disable();
    void reset();

    uint64_t readValue();
};
} // namespace Aux

#endif
