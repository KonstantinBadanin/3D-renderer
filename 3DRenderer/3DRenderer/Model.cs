using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using OpenTK;
using OpenTK.Graphics.OpenGL;
using OpenTK.Input;
using Assimp;

namespace _3DRenderer
{
    abstract class Model : IHaveNormals, IIntersectable
    {
        public abstract void Relocate(Coords arg);
        public abstract double Transparency { get; protected set; } //0-NO TRANSPARENCY, 1-FULL TRANSPARENCY
        public abstract double WarpCoef { get; protected set; }
        public abstract double Refl { get; protected set; }
        public abstract int Spec { get; protected set; }
        public abstract Color Color { get; protected set; }
        public abstract void IntersectRayModel(Coords origin, Coords direction, ref double t1, ref double t2, ref Coords normal);
        public abstract Coords GetUnitNormal(Coords hit_point, bool rayIn = true, Coords normal = null, Coords direction = null);

        protected Model(Color color, int spec, double refl, double warp)
        {
            Color = color; Spec = spec; WarpCoef = warp; Refl = refl > 1 ? 1 : refl < 0 ? 0 : refl; Transparency = 1 - Refl;
        }
    }

    class Polygonal : Model
    {
        public override void Relocate(Coords arg) { }
        public override double Transparency { get; protected set; }
        public override double WarpCoef { get; protected set; }
        public override double Refl { get; protected set; }
        public override int Spec { get; protected set; }
        public Scene Obj { get; protected set; }  //Only one object per one instance.
        //public List<Coords> Normals { get; protected set; }
        public override Color Color { get; protected set; }

        public override Coords GetUnitNormal(Coords hitPoint, bool rayIn, Coords normalUnit, Coords direction)
        {
            if (rayIn)
            {
                normalUnit = (normalUnit.CosA(direction) > 0) ? normalUnit : -normalUnit;
            }
            else
            {
                normalUnit = (normalUnit.CosA(direction) < 0) ? normalUnit : -normalUnit;
            }
            return normalUnit;
        }

        public override void IntersectRayModel(Coords origin, Coords direction, ref double t1, ref double t2, ref Coords normal)
        {
            double tMin = double.PositiveInfinity;
            foreach (Face face in Obj.Meshes.First().Faces)
            {
                Vector3D[] vertices = new Vector3D[3]
                {
                    Obj.Meshes.First().Vertices[face.Indices[0]],
                    Obj.Meshes.First().Vertices[face.Indices[1]],
                    Obj.Meshes.First().Vertices[face.Indices[2]]
                };
                Coords[] verts = new Coords[3];
                for (int i = 0; i < 3; i++)
                {
                    verts[i] = new Coords(vertices[i].X, vertices[i].Y, vertices[i].Z);
                }
                double v1 = 0.001;
                double v2 = double.PositiveInfinity;
                Plane tmp = new Plane(verts[0], verts[1], verts[2], Color.White);
                Coords proxy = new Coords();
                tmp.IntersectRayModel(origin, direction, ref v1, ref v2, ref proxy);
                if (v1 != double.PositiveInfinity)
                {
                    Coords e1 = verts[1] - verts[0];
                    Coords e2 = verts[2] - verts[0];
                    Coords hitPoint = v1 * direction + origin;
                    Coords pointMoved = hitPoint - verts[0];
                    double beta = (pointMoved.Y * e1.X - e1.Y * pointMoved.X) / (e2.Y * e1.X - e1.Y * e2.X);
                    double alpha = (pointMoved.X - e2.X * beta) / e1.X;
                    alpha = (double.IsNaN(alpha) || double.IsInfinity(alpha)) ? (pointMoved.Y - e2.Y * beta) / e1.Y : alpha;
                    alpha = (double.IsNaN(alpha) || double.IsInfinity(alpha)) ? (pointMoved.Z - e2.Z * beta) / e1.Z : alpha;
                    double sum = alpha + beta;
                    if ((alpha * beta > 0) && (sum <= 1) && (sum >= 0))
                    {
                        t1 = (hitPoint - origin).Lenght();
                        if (t1 < tMin)
                        {
                            tMin = t1;
                            normal = new Coords(tmp.A, tmp.B, tmp.C);
                            normal /= normal.Lenght();
                        }
                    }
                }
            }
            normal = (tMin == double.PositiveInfinity) ? new Coords() : normal;
            t1 = tMin;
            t2 = double.PositiveInfinity;
            return;
        }

        public List<Coords> ListOf3dVectorsToListOfCoords(List<Vector3D> lst)
        {
            List<Coords> res = new List<Coords>();
            foreach (Vector3D item in lst)
            {
                res.Add(new Coords(item.X, item.Y, item.Z));
            }
            return res;
        }

        public Polygonal(string path, Color color, int specularity, double reflection, double warp_coef) : base(color, specularity, reflection, warp_coef)
        {
            using (AssimpContext Model_raw = new AssimpContext())
            {
                //Normals = new List<Coords>();
                Obj = Model_raw.ImportFile(path);
                foreach (var item in Obj.Meshes)
                {
                    //Normals = ListOf3dVectorsToListOfCoords(item.Normals);
                    break;
                    //Protection from second mesh.
                }
            }
        }
    }

    class Sphere : Model  //Ok
    {
        public override void Relocate(Coords offset)
        {
            Center += offset ?? throw new ArgumentNullException(paramName: nameof(offset));
        }

        public override double Transparency { get; protected set; }
        public override double WarpCoef { get; protected set; }
        public override double Refl { get; protected set; }
        public override int Spec { get; protected set; }
        public Coords Center { get; protected set; }
        public double Radius { get; protected set; }
        public override Color Color { get; protected set; }

        public override void IntersectRayModel(Coords origin, Coords direction, ref double t1, ref double t2, ref Coords normal)
        {
            Coords oc = origin - Center;
            double k1 = direction.ScalarProd(direction);
            double k2 = 2 * oc.ScalarProd(direction);
            double k3 = oc.ScalarProd(oc) - Radius * Radius;
            double Discr = k2 * k2 - 4 * k1 * k3;
            if (Discr < 0)
            {
                t1 = t2 = double.PositiveInfinity;
                return;
            }
            else
            {
                t1 = (-k2 + Math.Sqrt(Discr)) / (2 * k1);
                t2 = (-k2 - Math.Sqrt(Discr)) / (2 * k1);
                return;
            }
        }

        public override Coords GetUnitNormal(Coords hit_point, bool rayIn, Coords normal, Coords direction)
        {
            _ = hit_point ?? throw new ArgumentNullException(paramName: nameof(hit_point));
            Coords k = hit_point - Center;
            k /= k.Lenght();
            return k;
        }

        public Sphere(Coords center, double radius, Color color, int specularity, double reflection, double warpCoef) : base(color, specularity, reflection, warpCoef)
        {
            Center = center ?? throw new ArgumentNullException(paramName: nameof(center)); Radius = radius;
        }
    }

    class Plane : Model
    {
        public override void Relocate(Coords arg) { }
        public override double Transparency { get; protected set; }
        public override double WarpCoef { get; protected set; } = 1;
        public override Color Color { get; protected set; }
        public override int Spec { get; protected set; }
        public override double Refl { get; protected set; }

        public override Coords GetUnitNormal(Coords hit_point, bool rayIn, Coords normal, Coords direction)
        {
            _ = hit_point ?? throw new ArgumentNullException(paramName: nameof(hit_point));
            Coords normalUnit = new Coords(A, B, C);
            normalUnit /= normalUnit.Lenght();
            return normalUnit;
        }

        public override void IntersectRayModel(Coords origin, Coords direction, ref double t1, ref double t2, ref Coords normal)
        {
            Coords coefs = new Coords(A, B, C);
            double numerator = -D - coefs.ScalarProd(origin);
            double denominator = coefs.ScalarProd(direction);
            double eps = 0.0001;
            t1 = ((denominator > eps) || (denominator < -eps)) ? numerator / denominator : double.PositiveInfinity;
            t1 = (t1 < 0) ? double.PositiveInfinity : t1;
            t2 = double.PositiveInfinity;
        }

        public double A { get; protected set; }
        public double B { get; protected set; }
        public double C { get; protected set; }
        public double D { get; protected set; }

        public Plane(double a, double b, double c, double d, Color color, int specularity = 0, double reflection = 0) : base(color, specularity, reflection, 1)
        {
            A = a; B = b; C = c; D = d;
        }

        public Plane(Coords a1, Coords a2, Coords a3, Color color, int specularity = 0, double reflection = 0) : base(color, specularity, reflection, 1)
        {
            A = (a2.Y - a1.Y) * (a3.Z - a1.Z) - (a2.Z - a1.Z) * (a3.Y - a1.Y);
            B = (a2.Z - a1.Z) * (a3.X - a1.X) - (a2.X - a1.X) * (a3.Z - a1.Z);
            C = (a2.X - a1.X) * (a3.Y - a1.Y) - (a2.Y - a1.Y) * (a3.X - a1.X);
            D = A * (-a1.X) + B * (-a1.Y) + C * (-a1.Z);
        }

        public override string ToString()
        {
            return A.ToString() + " " + B.ToString() + " " + C.ToString() + " " + D.ToString() + "\n";
        }
    }

    interface IIntersectable
    {
        void IntersectRayModel(Coords origin, Coords direction, ref double t1, ref double t2, ref Coords normal);
    }

    interface IHaveNormals
    {
        Coords GetUnitNormal(Coords hit_point, bool rayIn, Coords normal, Coords direction);
    }
}
