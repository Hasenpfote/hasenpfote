/*!
* @file axis_angle.h
* @brief Axis-angle class.
* @author Hasenpfote
* @date 2016/04/20
*/
#pragma once
#include "vector3.h"

namespace hasenpfote{ namespace math{

class AxisAngle final
{
public:
    Vector3 axis;   //!< an unit vector.
    float angle;    //!< an angle in radians.

public:
/* Constructor */

    AxisAngle() = default;
    AxisAngle(const AxisAngle& a);
    AxisAngle(const Vector3& axis, float angle);

/* Destructor */

    ~AxisAngle() = default;

/* Assignment operator */

    AxisAngle& operator = (const AxisAngle& a);

/* Debug */
    friend std::ostream& operator<<(std::ostream& os, const AxisAngle& a);
};

/* Inline */

}}