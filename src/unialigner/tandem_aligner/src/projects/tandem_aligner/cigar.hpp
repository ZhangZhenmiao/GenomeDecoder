//
// Created by Andrey Bzikadze on 05/19/22.
//

#pragma once

#include <iostream>
#include <vector>
#include <fstream>
#include <list>
#include "common/logging.hpp"

namespace tandem_aligner {

enum class CigarMode { M, I, D, X };

inline char cigar_mode2str(const CigarMode &fragment);

struct CigarFragment {
    int64_t length{0};
    CigarMode mode{};
};

class Cigar {
    std::list<CigarFragment> cigar_vec;

 public:
    Cigar() = default;

    [[nodiscard]] bool empty() const { return cigar_vec.empty(); }
    void extend(int64_t length, CigarMode mode);

    [[nodiscard]] size_t QueryLength() const;
    [[nodiscard]] size_t TargetLength() const;

    [[nodiscard]] double GetMutability() const;

    void AssertValidity(const std::string &target,
                        const std::string &query) const;

    friend std::ostream &operator<<(std::ostream &os, const Cigar &cigar);

    decltype(cigar_vec)::iterator begin() { return cigar_vec.begin(); }
    decltype(cigar_vec)::iterator end() { return cigar_vec.end(); }
    [[nodiscard]] decltype(cigar_vec)::const_iterator begin() const {
        return cigar_vec.begin();
    }
    [[nodiscard]] decltype(cigar_vec)::const_iterator end() const {
        return cigar_vec.end();
    }
    [[nodiscard]] decltype(cigar_vec)::const_iterator cbegin() const {
        return cigar_vec.cbegin();
    }
    [[nodiscard]] decltype(cigar_vec)::const_iterator cend() const {
        return cigar_vec.cend();
    }

    [[nodiscard]] int64_t Size() const { return cigar_vec.size(); }

    decltype(cigar_vec)::iterator
    AssignInterval(Cigar other, decltype(cigar_vec)::iterator it) {
        for (int i = 0; i < 2; ++i) {
            if (it!=cigar_vec.end()) {
                it = cigar_vec.erase(it);
            }
        }
        return cigar_vec.insert(it, other.cigar_vec.begin(),
                                other.cigar_vec.end());
    }

    decltype(cigar_vec)::iterator Erase(decltype(cigar_vec)::iterator pos) {
        return cigar_vec.erase(pos);
    }
    decltype(cigar_vec)::iterator Insert(decltype(cigar_vec)::iterator pos,
                                         const CigarFragment &value) {
        return cigar_vec.insert(pos, value);
    }

    void Summary(logging::Logger &) const;
};

std::ostream &operator<<(std::ostream &os, const Cigar &cigar);

} // namespace tandem_aligner