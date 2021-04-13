/*! \file   methods.cxx
 *   \brief  Source file for the methods routines.
 *
 *   For more details see the class documentation.
 *
 *   \date   21.06.2010
 *   \author Jan-Hendrik Koehne
 */

// #include <stdlib.h>

#include "PROPOSAL/methods.h"
#include <algorithm>
#include <string>

namespace PROPOSAL {

namespace Helper {

    std::string Centered(int width, const std::string& str, char fill)
    {
        int len = str.length();
        if (width < len) {
            return str;
        }

        int diff = width - len;
        int pad1 = diff / 2;
        int pad2 = diff - pad1;
        return std::string(pad1, fill) + str + std::string(pad2, fill);
    }

    bool case_insensitive_comp::operator()(
        std::string const& lhs_ref, std::string const& rhs_ref) const
    {
        auto lhs = std::string(lhs_ref);
        auto rhs = std::string(rhs_ref);
        transform(lhs.begin(), lhs.end(), lhs.begin(), ::tolower);
        transform(rhs.begin(), rhs.end(), rhs.begin(), ::tolower);

        return lhs < rhs;
    }

} // namespace Helper

} // namespace PROPOSAL
