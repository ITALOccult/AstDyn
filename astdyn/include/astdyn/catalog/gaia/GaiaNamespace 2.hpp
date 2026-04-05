#pragma once
// Namespace alias: users can refer to gaia types as astdyn::catalog::gaia::*
// while the implementation namespace remains ioc::gaia:: (unchanged internally).
#include <astdyn/catalog/gaia/types.h>
#include <astdyn/catalog/gaia/unified_gaia_catalog.h>

namespace astdyn {
namespace catalog {
namespace gaia = ::ioc::gaia;
}  // namespace catalog
}  // namespace astdyn
