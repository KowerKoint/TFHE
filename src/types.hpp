#include <cmath>
#include <cstdint>

namespace TFHE {
struct TorusValue {
    uint32_t val;

    TorusValue(double d) {
        d = fmod(d, 1.);
        if (d < 0.) d = 1. - d;
        d *= (1LL << 32);
        val = d;
    }

    TorusValue(bool b) {
        if (b)
            val = 1U << 29;
        else
            val = -(1U << 29);
    }

    operator double() const {
        double d = val;
        d /= (1LL << 32);
        return d;
    }

    operator bool() const { return !(val >> 31 & 1U); }
};
}  // namespace TFHE