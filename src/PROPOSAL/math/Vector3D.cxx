
#include <cmath>
#include <sstream>

#include "PROPOSAL/math/Vector3D.h"
#include "PROPOSAL/methods.h"

using namespace PROPOSAL;

void Vector3D::SetCoordinates(double val1, double val2, double val3) {
    coordinates = {val1, val2, val3};
}

bool Vector3D::operator==(const Vector3D& rhs) const {
    for (size_t i = 0; i<3; i++) {
        if (coordinates[i] != rhs.coordinates[i] )
            return false;
    }
    return true;
}

bool Vector3D::operator!=(const Vector3D& rhs) const {
    return !(*this == rhs);
}

double& Vector3D::operator[](size_t idx) {
    return coordinates.at(idx);
}

const double& Vector3D::operator[](size_t idx) const {
    return coordinates.at(idx);
}

namespace PROPOSAL {
    std::ostream &operator<<(std::ostream &os, const Vector3D &vector) {
        std::stringstream ss;
        ss << " Vector3D (" << &vector << ") ";
        os << Helper::Centered(60, ss.str()) << '\n';
        vector.print(os);
        os << Helper::Centered(60, "");
        return os;
    }
}