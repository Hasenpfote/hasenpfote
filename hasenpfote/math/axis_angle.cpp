#include <string>
#include "axis_angle.h"

namespace hasenpfote{ namespace math{

AxisAngle::AxisAngle(const AxisAngle& a)
{
    axis = a.axis;
    angle = a.angle;
}

AxisAngle::AxisAngle(const Vector3& axis, float angle)
{
    this->axis = axis;
    this->angle = angle;
}

AxisAngle& AxisAngle::operator = (const AxisAngle& a)
{
    axis = a.axis;
    angle = a.angle;
    return *this;
}

std::string AxisAngle::ToString() const
{
    return "AxisAngle{" + axis.ToString() + ", " + std::to_string(angle) + "}";
}

}}