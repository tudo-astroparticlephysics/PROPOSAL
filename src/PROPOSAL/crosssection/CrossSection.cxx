#include "PROPOSAL/crosssection/CrossSection.h"

namespace PROPOSAL {
template <>
size_t crosssection_hasher(
    size_t hash_diggest, std::shared_ptr<const EnergyCutSettings> const& ptr)
{
    hash_combine(hash_diggest, ptr->GetHash());
    return hash_diggest;
}
}
