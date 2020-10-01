#ifdef __linux__
#include <cstring>
#include <unistd.h>
#include <asm/unistd.h>
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#endif

#include <networkit/auxiliary/PerfEventCountHardware.hpp>

#include <networkit/auxiliary/Log.hpp>

#ifdef __linux__
namespace Aux {
namespace {
long perf_event_open(struct perf_event_attr *hw_event, pid_t pid, int cpu, int group_fd,
                     unsigned long flags) {
    int ret;

    ret = syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);
    return ret;
}
} // namespace
#endif

PerfEventCountHardware::PerfEventCountHardware(PerfEventCountHardware::Event event) : fd(-1) {
#ifdef __linux__
    struct perf_event_attr pe;
    memset(&pe, 0, sizeof(struct perf_event_attr));
    pe.type = PERF_TYPE_HARDWARE;
    pe.size = sizeof(struct perf_event_attr);
    pe.config = event;
    pe.disabled = 1;
    pe.exclude_kernel = 1;
    pe.exclude_hv = 1;

    fd = perf_event_open(&pe, 0, -1, -1, 0);

    if (fd == -1) {
        throw std::runtime_error("Error opening leader " + std::to_string(pe.config));
    }

    reset();
#endif
}

PerfEventCountHardware::PerfEventCountHardware(PerfEventCountHardware &&o) : fd(o.fd) {
    o.fd = -1;
}

void PerfEventCountHardware::enable() {
#ifdef __linux__
    ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
#endif
}

void PerfEventCountHardware::disable() {
#ifdef __linux__
    ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
#endif
}

void PerfEventCountHardware::reset() {
#ifdef __linux__
    ioctl(fd, PERF_EVENT_IOC_RESET, 0);
#endif
}

uint64_t PerfEventCountHardware::readValue() {
    uint64_t result = 0;
#ifdef __linux__
    read(fd, &result, sizeof(uint64_t));
#endif
    return result;
}

PerfEventCountHardware::~PerfEventCountHardware() {
#ifdef __linux__
    if (fd != -1) {
        disable();
        close(fd);
        fd = -1;
    }
#endif
}
} // namespace Aux
