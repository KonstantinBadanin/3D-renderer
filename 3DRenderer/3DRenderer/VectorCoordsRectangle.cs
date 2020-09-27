using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace _3DRenderer
{
    class Vector //OK
    {
        public Coords Origin { get; protected set; }
        public Coords Direction { get; protected set; }

        public Vector(Coords origin, Coords direction)
        {
            Origin = origin ?? throw new ArgumentNullException(paramName: nameof(origin));
            Direction = direction ?? throw new ArgumentNullException(paramName: nameof(direction));
        }

        public Vector()
        {
            Origin = Direction = new Coords();
        }

        public static Vector operator *(double factor, Vector a)
        //Returns multiplied vector located at the same origin.
        {
            if (a == null)
                throw new ArgumentNullException(paramName: nameof(a));
            return new Vector(a.Origin, factor * a.Direction);
        }
    }
    class Coords
    {
        public double ScalarProd(Coords b)
        {
            return X * b.X + Y * b.Y + Z * b.Z;
        }

        public string ToStringOfInts()
        {
            return Math.Round(X) + " " + Math.Round(Y) + " " + Math.Round(Z);
        }

        public double X { get; protected set; }
        public double Y { get; protected set; }
        public double Z { get; protected set; }

        public override bool Equals(object obj)
        {
            if (obj == null)
                throw new ArgumentNullException(paramName: nameof(obj));
            if (!GetType().Equals(obj.GetType()))
            {
                return false;
            }
            else
            {
                return this == (Coords)obj;
            }
        }

        public Coords VectProd(Coords arg)
        {
            return new Coords(Y * arg.Z - Z * arg.Y, Z * arg.X - X * arg.Z, X * arg.Y - Y * arg.X);
        }

        public Coords()
        {
            X = Y = Z = 0;
        }

        public double Lenght()
        {
            return Math.Sqrt(X * X + Y * Y + Z * Z);
        }

        public override int GetHashCode()
        {
            return base.GetHashCode();
        }

        public Coords(double x, double y, double z)
        {
            X = x; Y = y; Z = z;
        }

        public static Coords operator +(Coords a, Coords b)
        {
            if (a == null)
                throw new ArgumentNullException(paramName: nameof(a));
            if (b == null)
                throw new ArgumentNullException(paramName: nameof(b));
            return new Coords(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
        }

        public static Coords operator -(Coords a, Coords b)
        {
            if (a == null)
                throw new ArgumentNullException(paramName: nameof(a));
            if (b == null)
                throw new ArgumentNullException(paramName: nameof(b));
            return new Coords(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
        }

        public static Coords operator -(Coords a)
        {
            if (a == null)
                throw new ArgumentNullException(paramName: nameof(a));
            return new Coords(-a.X, -a.Y, -a.Z);
        }

        public static Coords operator *(double a, Coords b)
        {
            if (b == null)
                throw new ArgumentNullException(paramName: nameof(b));
            return new Coords(a * b.X, a * b.Y, a * b.Z);
        }

        public static Coords operator /(Coords a, double b)
        {
            if (a == null)
                throw new ArgumentNullException(paramName: nameof(a));
            return new Coords(a.X / b, a.Y / b, a.Z / b);
        }

        public static bool operator ==(Coords a, Coords b)
        {
            if ((a as object == null) || (b as object == null))
            {
                return false;
            }
            return (a.X == b.X) && (a.Y == b.Y) && (a.Z == b.Z);
        }

        public static bool operator !=(Coords a, Coords b)
        {
            if (a == null)
                throw new ArgumentNullException(paramName: nameof(a));
            if (b == null)
                throw new ArgumentNullException(paramName: nameof(b));
            return !(a == b);
        }

        public override string ToString()
        {
            return X.ToString() + " " + Y.ToString() + " " + Z.ToString() + "\n";
        }

        public double CosA(Coords b)
        {
            if (b == null)
                throw new ArgumentNullException(paramName: nameof(b));
            double eps = 0.0001;
            if ((Lenght() > eps) && (b.Lenght() > eps))
            {
                double cosA = ScalarProd(b) / (Lenght() * b.Lenght());
                return (cosA > 1) ? 1 : (cosA < -1) ? -1 : cosA;
            }
            else
                throw new ArgumentException("one or two arguments has 0 length");
        }
    }
}
