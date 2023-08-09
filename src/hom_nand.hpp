#pragma once

#include "tlwe_key_switching.hpp"
#include "trgsw.hpp"

namespace TFHE {
template <typename TLWELv0Parameter = TLWEParameter128BitSecurity,
    typename TLWELv1Parameter = TLWELv1ParameterDefault,
    typename TRLWEParameter = TRLWEParameterDefault,
    typename TRGSWParameter = TRGSWParameterDefault,
    typename TLWEKeySwitchingParameter = TLWEKeySwitchingParameterDefault>
class HomNAND {
    constexpr static int N_LV0 = TLWELv0Parameter::N;
    constexpr static int N_LV1 = TLWELv1Parameter::N;
    constexpr static int N = TRLWEParameter::N;
    constexpr static int K = TRGSWParameter::K;
    constexpr static int BG_BIT = TRGSWParameter::BG_BIT;
    constexpr static int L = TRGSWParameter::L;
    constexpr static int BASE_BIT = TLWEKeySwitchingParameter::BASE_BIT;
    constexpr static int T = TLWEKeySwitchingParameter::T;

private:
    TLWE<TLWELv0Parameter> tlwe_lv0;
    TLWE<TLWELv1Parameter> tlwe_lv1;
    TRLWE<TRLWEParameter> trlwe;
    TRGSW<TRGSWParameter> trgsw;
    TLWEKeySwitching<decltype(tlwe_lv0), decltype(tlwe_lv1),
        TLWEKeySwitchingParameter>
        key_switching;

public:
    struct SecretKey {
        Vector<bool, N_LV0> tlwe_lv0_s;
        Vector<Polynomial<bool, N>, K> trlwe_s;
    };
    struct EvaluateKey {
        Vector<Matrix<Polynomial<TorusValue, N>, (K + 1) * L, K + 1>, N_LV0> bk;
        Vector<
            Vector<Vector<Vector<TorusValue, N_LV0 + 1>, (1 << BASE_BIT)>, T>,
            N_LV1>
            ks;
    };
    using Cipher = Vector<TorusValue, N_LV0 + 1>;

    HomNAND()
        : tlwe_lv0{},
          tlwe_lv1{},
          trlwe{},
          trgsw{trlwe},
          key_switching{tlwe_lv0} {}

    SecretKey generate_secret_key() {
        return {tlwe_lv0.generate_s(), trlwe.generate_s()};
    }

    EvaluateKey make_evaluate_key(const SecretKey& secret_key) {
        EvaluateKey evaluate_key;
        evaluate_key.bk =
            trgsw.make_bk(secret_key.tlwe_lv0_s, secret_key.trlwe_s);
        auto tlwe_lv1_s = trlwe.extract_tlwe_lv0_key(secret_key.trlwe_s);
        evaluate_key.ks = key_switching.make_ks(secret_key.tlwe_lv0_s,
            trlwe.extract_tlwe_lv0_key(secret_key.trlwe_s));
        return evaluate_key;
    }

    Cipher encrypt(bool m, const SecretKey& secret_key) {
        return tlwe_lv0.encrypt_single_binary(m, secret_key.tlwe_lv0_s);
    }

    bool decrypt(const Cipher& c, const SecretKey& secret_key) {
        return tlwe_lv0.decrypt_single_binary(c, secret_key.tlwe_lv0_s);
    }

    Cipher nand(
        const Cipher& a, const Cipher& b, const EvaluateKey& evaluate_key) {
        Cipher c;
        c[0] = TorusValue(true);
        c -= a + b;
        auto lv1 = trgsw.gate_bootstrapping_tlwe_to_tlwe(c, evaluate_key.bk);
        return key_switching.identity_key_switch(lv1, evaluate_key.ks);
    }
};
}  // namespace TFHE