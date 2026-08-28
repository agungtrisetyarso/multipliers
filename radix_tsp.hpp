#pragma once
// Compact radix encoding of permutations / TSP, following
// Trisetyarso, Yulianti, Surendro, "Two multipliers, not one" (JOCO draft).
// All definitions match Sections 4–7 of the manuscript.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace joco {

inline int ceil_log2(int n) {
    int ell = 0;
    while ((1 << ell) < n) ++ell;
    return ell;
}

struct Encoding {
    int n = 0;
    int ell = 0;
    int q = 0;          // 2^ell - n illegal codes
    int m = 0;          // n * ell bits
    int alphabet = 0;   // 2^ell
};

inline Encoding make_encoding(int n) {
    Encoding e;
    e.n = n;
    e.ell = ceil_log2(n);
    e.alphabet = 1 << e.ell;
    e.q = e.alphabet - n;
    e.m = n * e.ell;
    return e;
}

// Packed bitstring: bit r of register t lives at bit (t*ell + r), LSB first.
inline int register_value(std::uint64_t bits, int t, const Encoding& e) {
    const std::uint64_t mask = (e.ell == 64) ? ~0ull : ((1ull << e.ell) - 1ull);
    return static_cast<int>((bits >> (t * e.ell)) & mask);
}

inline std::uint64_t set_register(std::uint64_t bits, int t, int value, const Encoding& e) {
    const std::uint64_t mask = (e.ell == 64) ? ~0ull : ((1ull << e.ell) - 1ull);
    bits &= ~(mask << (t * e.ell));
    bits |= (static_cast<std::uint64_t>(value) & mask) << (t * e.ell);
    return bits;
}

struct PenaltyObjective {
    int P = 0;          // P_ill + P_col  (Lemma 4.4)
    int Pill = 0;
    int Pcol = 0;
    double C = 0.0;     // C_TSP          (Definition 4.5)
    int Z = 0;          // number of zero-cost cyclic arcs
};

inline PenaltyObjective evaluate(std::uint64_t bits, const Encoding& e,
                                 const std::vector<std::vector<double>>& d) {
    PenaltyObjective out;
    std::vector<int> val(e.n);
    std::vector<int> freq(e.alphabet, 0);
    for (int t = 0; t < e.n; ++t) {
        val[t] = register_value(bits, t, e);
        freq[val[t]]++;
        if (val[t] >= e.n) out.Pill++;
    }
    for (int u = 0; u < e.alphabet; ++u) {
        const int f = freq[u];
        if (f >= 2) out.Pcol += f * (f - 1) / 2;
    }
    out.P = out.Pill + out.Pcol;
    for (int t = 0; t < e.n; ++t) {
        const int u = val[t];
        const int v = val[(t + 1) % e.n];
        if (u < e.n && v < e.n) {
            out.C += d[u][v];
            if (d[u][v] == 0.0) out.Z++;
        } else {
            out.Z++;
        }
    }
    return out;
}

// Left-to-right repair Π_LR (Definition 7.2).
inline std::vector<int> repair_LR(std::uint64_t bits, const Encoding& e) {
    std::vector<int> raw(e.n);
    for (int t = 0; t < e.n; ++t) raw[t] = register_value(bits, t, e) % e.n;
    std::vector<int> used(e.n, 0);
    std::vector<int> pi(e.n, -1);
    for (int t = 0; t < e.n; ++t) {
        if (!used[raw[t]]) {
            pi[t] = raw[t];
            used[raw[t]] = 1;
        }
    }
    int next = 0;
    for (int t = 0; t < e.n; ++t) {
        if (pi[t] < 0) {
            while (used[next]) ++next;
            pi[t] = next;
            used[next] = 1;
        }
    }
    return pi;
}

// Legal-first repair Π_LF (Definition 7.2).
inline std::vector<int> repair_LF(std::uint64_t bits, const Encoding& e) {
    std::vector<int> used(e.n, 0);
    std::vector<int> pi(e.n, -1);
    for (int t = 0; t < e.n; ++t) {
        const int v = register_value(bits, t, e);
        if (v < e.n && !used[v]) {
            pi[t] = v;
            used[v] = 1;
        }
    }
    int next = 0;
    for (int t = 0; t < e.n; ++t) {
        if (pi[t] < 0) {
            while (used[next]) ++next;
            pi[t] = next;
            used[next] = 1;
        }
    }
    return pi;
}

inline int changed_set_size(std::uint64_t bits, const std::vector<int>& pi, const Encoding& e) {
    int a = 0;
    for (int t = 0; t < e.n; ++t)
        if (pi[t] != register_value(bits, t, e)) ++a;
    return a;
}

inline double tour_cost(const std::vector<int>& pi, const std::vector<std::vector<double>>& d) {
    const int n = static_cast<int>(pi.size());
    double c = 0;
    for (int t = 0; t < n; ++t) c += d[pi[t]][pi[(t + 1) % n]];
    return c;
}

// Distance tables of Theorem 5.5 / Examples 10.x.
inline std::vector<std::vector<double>> unit_complete(int n) {
    std::vector<std::vector<double>> d(n, std::vector<double>(n, 1.0));
    for (int i = 0; i < n; ++i) d[i][i] = 0;
    return d;
}

inline std::vector<std::vector<double>> spread_cycle(int n, double M) {
    std::vector<std::vector<double>> d(n, std::vector<double>(n, M));
    for (int i = 0; i < n; ++i) {
        d[i][i] = 0;
        d[i][(i + 1) % n] = 1;
        d[i][(i + n - 1) % n] = 1;
    }
    return d;
}

inline std::vector<std::vector<double>> path_metric(int n) {
    std::vector<std::vector<double>> d(n, std::vector<double>(n, 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) d[i][j] = std::abs(i - j);
    return d;
}

// Unit square on 4 vertices in cyclic order 0-1-2-3.
inline std::vector<std::vector<double>> unit_square() {
    const double s2 = std::sqrt(2.0);
    std::vector<std::vector<double>> d(4, std::vector<double>(4, 0));
    d[0][1] = d[1][0] = 1;
    d[1][2] = d[2][1] = 1;
    d[2][3] = d[3][2] = 1;
    d[3][0] = d[0][3] = 1;
    d[0][2] = d[2][0] = s2;
    d[1][3] = d[3][1] = s2;
    return d;
}

// Optimal tour length for the three closed-form families.
inline double c_star_unit_complete(int n) { return static_cast<double>(n); }
inline double c_star_spread(int n) { return static_cast<double>(n); }
inline double c_star_path(int n) { return 2.0 * (n - 1); }

inline double dmax_of(const std::vector<std::vector<double>>& d) {
    double m = 0;
    for (auto& row : d)
        for (double v : row) m = std::max(m, v);
    return m;
}

struct EnumStats {
    double L = std::numeric_limits<double>::infinity();
    double maxC = 0;
    double cstar = std::numeric_limits<double>::infinity();
    double UF = 0;
    int nu = std::numeric_limits<int>::max();          // min {P : C=0, P>0}
    double minC_P1 = std::numeric_limits<double>::infinity();
    double lambda_star = 0;
    double lambda_tr_LF = 0;
    double lambda_tr_LR = 0;
    int max_locality_num_LF = 0;   // max 2|A|
    int max_locality_den_LF = 1;   // corresponding P
    int max_locality_num_LR = 0;
    int max_locality_den_LR = 1;
    int n_feasible = 0;
    int cube = 0;
    // path-metric census of sharp strings (P=1, C = n-2)
    int sharp = 0;
    int sharp_LR_opt = 0, sharp_LR_mid = 0, sharp_LR_worst = 0;
    double sharp_LR_max = 0;
    int sharp_LF_opt = 0;
};

}  // namespace joco
