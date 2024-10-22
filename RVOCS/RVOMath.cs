using Lockstep.Math;
using System;

namespace RVO
{
    public struct RVOMath
    {
        public static LFloat RVO_EPSILON = 1;//0.1 1¼Ó5¸ö0

        public static LFloat abs(LVector2 vector)
        {

            return sqrt(absSq(vector));
        }

        public static LFloat absSq(LVector2 vector)
        {
            return vector * vector;
        }

        public static LVector2 normalize(LVector2 vector)
        {
            if (vector.x() == 0 &&vector.y() == 0)
            {
                return vector;
            }
            return vector / abs(vector);
        }

        internal static LFloat det(LVector2 vector1, LVector2 vector2)
        {
            return vector1.x_ * vector2.y_ - vector1.y_ * vector2.x_;
        }

        internal static LFloat distSqPointLineSegment(LVector2 vector1, LVector2 vector2, LVector2 vector3)
        {
            LFloat r = ((vector3 - vector1) * (vector2 - vector1)) / absSq(vector2 - vector1);

            if (r < 0.0f)
            {
                return absSq(vector3 - vector1);
            }

            if (r > 1.0f)
            {
                return absSq(vector3 - vector2);
            }

            return absSq(vector3 - (vector1 + r * (vector2 - vector1)));
        }

        internal static LFloat fabs(LFloat scalar)
        {
            return Math.Abs(scalar);
        }

        internal static LFloat leftOf(LVector2 a, LVector2 b, LVector2 c)
        {
            return det(a - c, b - a);
        }

        internal static LFloat sqr(LFloat scalar)
        {
            return scalar * scalar;
        }

        internal static LFloat sqrt(LFloat scalar)
        {
            return (LFloat)Math.Sqrt((float)scalar);
        }
    }
}
