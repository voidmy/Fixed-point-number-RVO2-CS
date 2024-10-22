using Lockstep.Math;

namespace RVO
{
    public struct LVector2
    {
        internal LFloat x_;
        internal LFloat y_;

        public LVector2(LFloat x, LFloat y)
        {
            x_ = x;
            y_ = y;
        }

        public override string ToString()
        {
            return "(" + x_.ToString() + "," + y_.ToString() + ")";
        }

        public LFloat x()
        {
            return x_;
        }

        public LFloat y()
        {
            return y_;
        }

        public static LFloat operator *(LVector2 vector1, LVector2 vector2)
        {
            // operator '*' might need to be adjusted for LFloat type
            return vector1.x_ * vector2.x_ + vector1.y_ * vector2.y_;
        }

        public static LVector2 operator *(LFloat scalar, LVector2 vector)
        {
            return vector * scalar;
        }

        public static LVector2 operator *(LVector2 vector, LFloat scalar)
        {
            return new LVector2(vector.x_ * scalar, vector.y_ * scalar);
        }

        public static LVector2 operator /(LVector2 vector, LFloat scalar)
        {
            return new LVector2(vector.x_ / scalar, vector.y_ / scalar);
        }

        public static LVector2 operator +(LVector2 vector1, LVector2 vector2)
        {
            return new LVector2(vector1.x_ + vector2.x_, vector1.y_ + vector2.y_);
        }

        public static LVector2 operator -(LVector2 vector1, LVector2 vector2)
        {
            return new LVector2(vector1.x_ - vector2.x_, vector1.y_ - vector2.y_);
        }

        public static LVector2 operator -(LVector2 vector)
        {
            return new LVector2(-vector.x_, -vector.y_);
        }
    }
}
