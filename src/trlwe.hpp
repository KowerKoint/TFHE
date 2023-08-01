#include <algorithm>
#include <iostream>

#include "random.hpp"
#include "types.hpp"

namespace TFHE {
struct TRLWEParameterDefault {
    constexpr static int N = 512;
    constexpr static int K = 2;
    constexpr static double ALPHA = 2.9802322387695312e-08;  // 2^(-25)
};

template <typename TRLWEParameter = TRLWEParameterDefault>
class TRLWE {
public:
    constexpr static int N = TRLWEParameter::N;
    constexpr static int K = TRLWEParameter::K;
    constexpr static double ALPHA = TRLWEParameter::ALPHA;

private:
    Random random;
    Vector<Polynomial<TorusValue, N>, K> generate_a() {
        Vector<Polynomial<TorusValue, N>, K> a;
        for (int i = 0; i < K; i++) {
            std::generate(a[i].begin(), a[i].end(),
                [this] { return random.uniform_torus(); });
        }
        return a;
    }
    Polynomial<TorusValue, N> generate_e() {
        Polynomial<TorusValue, N> e;
        std::generate(
            e.begin(), e.end(), [this] { return random.normal() * ALPHA; });
        return e;
    }

    template <int N>
    constexpr static Vector<Polynomial<TorusValue, N>, K + 1> concat_ba(
        const Polynomial<TorusValue, N>& b,
        const Vector<Polynomial<TorusValue, N>, K>& a) {
        Vector<Polynomial<TorusValue, N>, K + 1> ba;
        ba[0] = b;
        for (int i = 0; i < K; i++) ba[1 + i] = a[i];
        return ba;
    }

    template <int N>
    constexpr static std::pair<Polynomial<TorusValue, N>,
        Vector<Polynomial<TorusValue, N>, K>>
    decompose_ba(const Vector<Polynomial<TorusValue, N>, K + 1>& ba) {
        Polynomial<TorusValue, N> b = ba[0];
        Vector<Polynomial<TorusValue, N>, K> a;
        for (int i = 0; i < K; i++) a[i] = ba[1 + i];
        return std::make_pair(b, a);
    }

public:
    TRLWE() : random{} {}

    Vector<Polynomial<bool, N>, K> generate_s() {
        Vector<Polynomial<bool, N>, K> s;
        for (int i = 0; i < K; i++) {
            s[i] = random.bit_polynomial<N>();
        }
        return s;
    }

    Vector<Polynomial<TorusValue, N>, K + 1> encrypt_binary_polynomial(
        const Polynomial<bool, N>& m, const Vector<Polynomial<bool, N>, K>& s) {
        Polynomial<TorusValue, N> m_torus;
        std::transform(m.begin(), m.end(), m_torus.begin(),
            [](bool b) { return TorusValue(b); });
        auto a = generate_a();
        auto e = generate_e();
        auto b = a.dot(s) + m_torus + e;
        auto tmp = b - a.dot(s);
        return concat_ba(b, a);
    }

    Polynomial<bool, N> decrypt_binary_polynomial(
        const Vector<Polynomial<TorusValue, N>, K + 1>& ba,
        const Vector<Polynomial<bool, N>, K>& s) {
        auto [b, a] = decompose_ba<N>(ba);
        auto tmp = b - a.dot(s);
        auto torus_vector = b - a.dot(s);
        Polynomial<bool, N> ans;
        std::transform(torus_vector.begin(), torus_vector.end(), ans.begin(),
            [](const TorusValue& val) { return (bool)val; });
        return ans;
    }
};
}  // namespace TFHE