#ifndef NETWORKIT_AUXILIARY_SPIN_LOCK_HPP_
#define NETWORKIT_AUXILIARY_SPIN_LOCK_HPP_

// networkit-format

#include <atomic>

namespace Aux {

class Spinlock {
public:
    void lock() {
        while (spinner.test_and_set(std::memory_order_acquire)) {
            /* spin */
        }
    }
    void unlock() { spinner.clear(std::memory_order_release); }

private:
    std::atomic_flag spinner = ATOMIC_FLAG_INIT;
};

class SpinlockArray {
public:
    SpinlockArray(size_t n) : spinlocks(n) {}
    SpinlockArray(const SpinlockArray &other) : SpinlockArray(other.size()) {}
    Spinlock &operator[](size_t idx) { return spinlocks[idx]; }
    size_t size() const { return spinlocks.size(); }

private:
    std::vector<Spinlock> spinlocks;
};

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_SPIN_LOCK_HPP_
