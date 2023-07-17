#include <random>

namespace TFHE {
class Random {
    std::random_device rd;
    std::normal_distribution<> normal_dist;
    std::uniform_real_distribution<> uniform_dist;

public:
    Random() : rd{}, normal_dist{}, uniform_dist{} {}

    double normal() { return normal_dist(rd); }

    double uniform() { return uniform_dist(rd); }
};
}  // namespace TFHE