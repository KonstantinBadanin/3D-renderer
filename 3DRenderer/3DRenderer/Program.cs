using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Drawing;
using OpenTK;
using OpenTK.Graphics.OpenGL;
using System.IO;
using Assimp;

namespace _3DRenderer
{
    class Proxy
    {
        public int CurrentXPosition { get; protected set; }
        public int CurrentYPosition { get; protected set; }
        public List<List<Vector>> ListsOrigins { get; protected set; }
        public List<List<Coords>> ListsCorners { get; protected set; }
        public int Limit { get; protected set; }
        public double Max { get; protected set; }

        public Proxy(List<List<Vector>> origins, List<List<Coords>> corners, double max, int limit, int i, int j)
        {
            ListsOrigins = origins ?? throw new ArgumentNullException(paramName: nameof(origins));
            Max = max; Limit = limit;
            ListsCorners = corners ?? throw new ArgumentNullException(paramName: nameof(corners));
            CurrentXPosition = i; CurrentYPosition = j;
        }
    }

    class RayInfo
    {
        public Coords Origin { get; private set; }
        public Coords Direction { get; private set; }
        public int Limit { get; private set; }
        public bool In { get; private set; }

        public override string ToString()
        {
            return Origin.ToStringOfInts() + " " + In.ToString() + " " + Limit.ToString();
        }

        public RayInfo(Coords origin, Coords direction, int limit1, bool flag)
        {
            Origin = origin ?? throw new ArgumentNullException(paramName: nameof(origin));
            if (direction == null)
                throw new ArgumentNullException(paramName: nameof(direction));
            var a = direction.Lenght();
            double eps = 0.001;
            if ((a < 1 + eps) && (a > 1 - eps))
                Direction = direction;
            else
                throw new ArgumentException("argument must have length 1");
            Limit = limit1; In = flag;
        }
    }

    class ThreadStruct
    {
        public List<Coords> Lst { get; protected set; } = new List<Coords>();
        public Color Colour { get; protected set; }

        public ThreadStruct(List<Coords> lst, Color colour)
        {
            Lst = lst ?? throw new ArgumentNullException(paramName: nameof(lst));
            Colour = colour;
        }
    }

    class ThreadList
    {
        public List<ThreadStruct> Lst { get; protected set; } = new List<ThreadStruct>();
        public readonly object locker = new object();

        public ThreadList()
        {
            Lst = new List<ThreadStruct>();
            locker = new object();
        }

        public void Add(List<Coords> points, Color color)
        {
            lock (locker)
            {
                Lst.Add(new ThreadStruct(points, color));
            }
        }

        public void DrawLst()
        {
            for (int item = 0; item < Lst.Count(); item++)
            {
                GL.Begin(OpenTK.Graphics.OpenGL.PrimitiveType.Quads);
                GL.Color3(Lst[item].Colour);
                GL.Vertex3(Lst[item].Lst[0].X, Lst[item].Lst[0].Y, Lst[item].Lst[0].Z);
                GL.Vertex3(Lst[item].Lst[1].X, Lst[item].Lst[1].Y, Lst[item].Lst[1].Z);
                GL.Vertex3(Lst[item].Lst[3].X, Lst[item].Lst[3].Y, Lst[item].Lst[3].Z);
                GL.Vertex3(Lst[item].Lst[2].X, Lst[item].Lst[2].Y, Lst[item].Lst[2].Z);
                GL.End();
                Console.WriteLine(Lst[item].Lst[0].X + " " + Lst[item].Lst[0].Y + "\n");
            }
            return;
        }
    }

    enum WarpState : int
    {
        FullInnerReflection,
        WarpAndInnerReflection,
        WarpAndOuterReflection
    }

    class Program
    {
        private static double Eps { get; set; } = 0.0001;
        private static bool BlurReflections { get; set; } = false;

        private static List<Model> Scene = new List<Model>
            {
                //Addition of objects to the scene.
                //new Plane(25,10,-19,500,Color.White),
                //new Sphere(new Coords(-620, 100, 950), 100, Color.Pink, 100,0.2,1.1),
                new Sphere(new Coords(450, 50, 800), 100, Color.Red, 500,0.5,2.5),
                new Sphere(new Coords(50, 350, 900), 150, Color.Green, 500,0.3,1.1),
                new Sphere(new Coords(50, 100, 1500), 70, Color.Pink, 500,0.3,1.5),
                new Sphere(new Coords(0, -1000, 1000), 950, Color.Yellow, 25,0.6,3),
                new Polygonal(@"C:\Users\kaste\3D-renderer\3DRenderer\3DRenderer\test\cvv.fbx",Color.Blue, 500,0.8,1.5)
                //new Sphere(new Coords(500, 350, 1000), 49.9, Color.FromArgb(255,255,255), -1,0)
                //new Plane(0,0.001,0.001,-2,Color.Red,-1,1)
            };

        private static List<Light> Lights = new List<Light>
            {
                //Addition of light to the scene.
                //new Dirl(new Coords(5,2,5),0.2),         
                //new Point(new Coords(0, 0, 100), 0.4),
                new Point(new Coords(481, 672, 962), 0.4),
                //new Point(new Coords(700, 500, 25), 0.4),
                //new Point(new Coords(650, 500, 25), 0.4),
                //new Point(new Coords(600, 500, 25), 0.4),
                //new Point(new Coords(550, 500, 25), 0.4),
                new Ambient(0.4),
            };

        private static Color GetMixedColor(Color a, Color b, double reflOrTransp)
        {
            int k = (int)(a.R * (1 - reflOrTransp) + b.R * reflOrTransp);
            int l = (int)(a.G * (1 - reflOrTransp) + b.G * reflOrTransp);
            int m = (int)(a.B * (1 - reflOrTransp) + b.B * reflOrTransp);
            k = (k >= 0) ? (k <= 255) ? k : 255 : 0;
            l = (l >= 0) ? (l <= 255) ? l : 255 : 0;
            m = (m >= 0) ? (m <= 255) ? m : 255 : 0;
            return Color.FromArgb(k, l, m);
        }

        private static void WrapTraceRay(object rayData)
        {
            bool ray_in = false;
            Proxy obj = (rayData is Proxy) ? rayData as Proxy : null;
            for (int l = 0; l < obj.ListsOrigins.Count(); l++)
            {
                int i = 0; int k = 0; int j = 0;
                for (int m = 0; m < obj.ListsOrigins[l].Count(); m++)
                //Supersampling.
                {
                    Color a1 = TraceRayUnit(
                        new RayInfo(obj.ListsOrigins[l][m].Origin, obj.ListsOrigins[l][m].Direction, obj.Limit, ray_in),
                        0, double.PositiveInfinity, obj.CurrentXPosition, obj.CurrentYPosition);
                    i += a1.R; k += a1.G; j += a1.B;
                }
                _thread_pool.Add(obj.ListsCorners[l], Color.FromArgb(255, i / 4, k / 4, j / 4));
            }
            return;
        }

        private static Coords DivRay1(Coords refl_unit, Coords normal_unit, double rnd1, double rnd2)
        {
            Coords refl_proj_normal = refl_unit.CosA(normal_unit) * normal_unit;
            double i = (rnd1 / 5 + 0.8);
            double j = (rnd2 / 5 + 0.8);
            Coords refl_corrected_unit = i * (refl_unit - refl_proj_normal) + j * refl_proj_normal;
            refl_corrected_unit /= refl_corrected_unit.Lenght();
            return refl_corrected_unit;
        }

        private static Model GetClosestModelAndDistanceToIt(RayInfo ray, double tMin, double tMax, ref double closestT, ref Coords normalForPolygonal)
        {
            Model closestModel = null;
            closestT = (ray.In) ?
                ClosestIntersection(ray.Origin, ray.Direction, tMin, tMax, ref Scene, ref closestModel, ref normalForPolygonal) :
                ClosestIntersection(ray.Origin, ray.Direction, 0.0001, tMax, ref Scene, ref closestModel, ref normalForPolygonal);
            return closestModel;
        }

        private static Color TraceRayUnit(RayInfo ray, double tMin, double tMax, int it, int jt)
        {
            double closestT = 0;
            Coords normalForPolygonal = new Coords();
            Model closestModel = GetClosestModelAndDistanceToIt(ray, tMin, tMax, ref closestT, ref normalForPolygonal);
            if (closestModel == null)
            {
                return Color.Black;
            }
            Coords rayHit = ray.Origin + closestT * ray.Direction;
            Coords normal_unit = closestModel.GetUnitNormal(rayHit, ray.In, normalForPolygonal, ray.Direction);
            Color local_color = (ray.In) ? Color.Black :
                ComputeColor(ComputeLightning(
                        rayHit, normal_unit, Lights, Scene, -ray.Direction, closestModel.Spec
                        ), closestModel);
            Coords refl_unit = ReflRayUnit(-ray.Direction, normal_unit, ray.In);
            if (ray.Limit < 1)
            {
                return local_color;
            }
            if (closestModel.Transparency < Eps)
            {
                double refl = 0.8;
                Color refl_color =
                    TraceRayUnit(new RayInfo(rayHit, refl_unit, ray.Limit - 1, ray.In), 0.001, double.PositiveInfinity, it, jt);
                return GetMixedColor(local_color, refl_color, refl);
            }
            else
            {
                var warped = WarpRay(ray.Direction, normal_unit, closestModel.WarpCoef, ray.In);
                switch (warped.Item2)
                {
                    case WarpState.FullInnerReflection:
                        {
                            double refl = 0.8;
                            Color refl_color = GetReflColor(refl_unit, normal_unit, rayHit, ray.Limit, ray.In, it, jt);
                            return GetMixedColor(local_color, refl_color, refl);
                        }
                    case WarpState.WarpAndInnerReflection:
                        {
                            double refl = 0.9;
                            Color refl_color = GetReflColor(refl_unit, normal_unit, rayHit, ray.Limit, ray.In, it, jt);
                            Color warped_color =
                                TraceRayUnit(new RayInfo(rayHit, warped.Item1, ray.Limit - 1, !ray.In), 0.001, double.PositiveInfinity, it, jt);
                            return GetMixedColor(local_color, GetMixedColor(
                                warped_color, refl_color, closestModel.Refl
                                ), refl);
                        }
                    case WarpState.WarpAndOuterReflection:
                        {
                            double refl = 0.8;
                            Color refl_color = GetReflColor(refl_unit, normal_unit, rayHit, ray.Limit, ray.In, it, jt);
                            Color warped_color =
                                TraceRayUnit(new RayInfo(rayHit, warped.Item1, ray.Limit - 1, !ray.In), 0.001, double.PositiveInfinity, it, jt);
                            return GetMixedColor(local_color, GetMixedColor(
                                warped_color, refl_color, closestModel.Refl
                                ), refl);
                        }
                    default:
                        return Color.White;
                }

            }
        }

        private static Color GetReflColor(Coords refl_unit, Coords normal_unit, Coords ray_hit, int limit, bool ray_in, int it, int jt)
        {
            return BlurReflections ? GetBlurredColor(refl_unit, normal_unit, ray_hit, limit, ray_in, it, jt) :
                TraceRayUnit(new RayInfo(ray_hit, refl_unit, limit - 1, ray_in), 0.001, double.PositiveInfinity, it, jt);
        }

        private static Color GetBlurredColor(Coords refl_unit, Coords normal_unit, Coords ray_hit, int limit, bool ray_in, int it, int jt)
        {
            int a1 = 0; int a2 = 0; int a3 = 0;
            int numrays = 30;
            int op = 0;
            Random RandGen = new Random();
            while (op < numrays)
            {
                Coords refl_corrected_unit = DivRay1(refl_unit, normal_unit, RandGen.NextDouble(), RandGen.NextDouble());
                if (normal_unit.CosA(refl_corrected_unit) > 0)
                {
                    Color b = TraceRayUnit(new RayInfo(ray_hit, refl_corrected_unit, limit - 1, ray_in), 0.001, double.PositiveInfinity, it, jt);
                    a1 += b.R;
                    a2 += b.G;
                    a3 += b.B;
                    op++;
                }
            }
            a1 /= numrays; a2 /= numrays; a3 /= numrays;
            return Color.FromArgb(255, a1, a2, a3);
        }

        private static Coords BuildResultVector(Coords incomingRayUnit, double cosAlpha, Coords normalUnit, double sinGamma, bool rayIn)
        {
            double eps = 0.001;
            incomingRayUnit /= cosAlpha;
            normalUnit = rayIn ? -normalUnit : normalUnit;
            Coords border_unit = incomingRayUnit + normalUnit;
            border_unit /= border_unit.Lenght();
            border_unit = double.IsNaN(border_unit.X) ? new Coords() : border_unit;
            var k = border_unit.Lenght();
            if ((border_unit != new Coords()) && !((k < 1 + eps) && (k > 1 - eps)))
            {
                throw new Exception();
            }
            double cos_gamma = Math.Sqrt(1 - sinGamma * sinGamma);
            cos_gamma = (border_unit == new Coords()) ? 1 : cos_gamma;
            return sinGamma * border_unit - cos_gamma * normalUnit;
        }

        private static Tuple<Coords, WarpState> WarpRay(Coords incoming_ray_unit, Coords normal_unit, double coef, bool ray_in)
        {
            if (!((incoming_ray_unit.Lenght() < 1.001) && (incoming_ray_unit.Lenght() > 0.999)))
            {
                throw new Exception();
            }
            if (!((normal_unit.Lenght() < 1.001) && (normal_unit.Lenght() > 0.999)))
            {
                throw new Exception();
            }
            coef = (ray_in) ? 1 / coef : coef;
            double eps = 0.001;
            double cos_alpha = incoming_ray_unit.ScalarProd(normal_unit);
            cos_alpha = cos_alpha < 0 ? -cos_alpha : cos_alpha;
            double sin_alpha = Math.Sqrt(1 - cos_alpha * cos_alpha);
            double sin_gamma = sin_alpha / coef;
            sin_gamma = (sin_gamma > 1) ? 1 : (sin_gamma < -1) ? -1 : sin_gamma;
            if (coef < 1)
            {
                if (sin_alpha < coef)
                {
                    Coords warped_ray_unit = BuildResultVector(incoming_ray_unit, cos_alpha, normal_unit, sin_gamma, ray_in);
                    if (double.IsNaN((warped_ray_unit).X) || !((warped_ray_unit.Lenght() < 1 + eps) && (warped_ray_unit.Lenght() > 1 - eps)))
                    {
                        throw new Exception();
                    }
                    return new Tuple<Coords, WarpState>
                        (warped_ray_unit, WarpState.WarpAndInnerReflection);
                }
                else
                {
                    return new Tuple<Coords, WarpState>(new Coords(), WarpState.FullInnerReflection);
                }
            }
            else
            {
                Coords warped_ray_unit = BuildResultVector(incoming_ray_unit, cos_alpha, normal_unit, sin_gamma, ray_in);
                if (double.IsNaN((warped_ray_unit).X) || !((warped_ray_unit.Lenght() < 1 + eps) && (warped_ray_unit.Lenght() > 1 - eps)))
                {
                    throw new Exception();
                }
                return new Tuple<Coords, WarpState>
                    (warped_ray_unit, WarpState.WarpAndOuterReflection);
            }
        }

        private static Color ComputeColor(double light, Model a)
        {
            int k = (int)(a.Color.R * light);
            int l = (int)(a.Color.G * light);
            int m = (int)(a.Color.B * light);
            k = (k >= 0) ? (k <= 255) ? k : 255 : 0;
            l = (l >= 0) ? (l <= 255) ? l : 255 : 0;
            m = (m >= 0) ? (m <= 255) ? m : 255 : 0;
            return Color.FromArgb(k, l, m);
        }

        private static double ComputeLightning(Coords ray_hit, Coords normal_unit, List<Light> light, List<Model> scene, Coords direction_inverted, int specularity)
        {
            double i = 0;
            foreach (var item in light)
            {
                if (item is Ambient) i += item.Power;
                else
                {
                    Coords l = new Coords();
                    double t_max = 0;
                    if (item is Directed)
                    {
                        l = (item as Directed).Direction;
                        t_max = double.PositiveInfinity;
                    }
                    if (item is Point)
                    {
                        l = (item as Point).Position - ray_hit;
                        t_max = 1;
                    }
                    Model sh_sphere = null;
                    Coords proxy = new Coords();
                    double sh_t = ClosestIntersection(ray_hit, l, 0.001, t_max, ref scene, ref sh_sphere, ref proxy);
                    if (sh_sphere != null) continue;
                    double h = normal_unit.ScalarProd(l);
                    if (h > 0) i += item.Power * normal_unit.CosA(l);
                    if (specularity != -1)
                    {
                        Coords r = ReflRayUnit(normal_unit, l, false);
                        double rv = r.ScalarProd(direction_inverted);
                        if (rv > 0) i += item.Power * Math.Pow(r.CosA(direction_inverted), specularity);
                    }
                }
            }
            return i;
        }

        private static double ClosestIntersection(Coords origin, Coords direction, double t_min, double t_max, ref List<Model> scene, ref Model mod, ref Coords normal)
        {
            double closest_t = double.PositiveInfinity;
            foreach (var item in scene)
            {
                double t1 = -1, t2 = -1;
                if (item is Sphere)
                {
                    (item as Sphere).IntersectRayModel(origin, direction, ref t1, ref t2, ref normal);
                    if ((t1 >= t_min) && (t1 <= t_max) && (t1 < closest_t))
                    {
                        closest_t = t1; mod = item as Sphere;
                    }
                    if ((t2 >= t_min) && (t2 <= t_max) && (t2 < closest_t))
                    {
                        closest_t = t2; mod = item as Sphere;
                    }
                }
                if (item is Plane)
                {
                    (item as Plane).IntersectRayModel(origin, direction, ref t1, ref t2, ref normal);
                    if ((t1 >= t_min) && (t1 <= t_max) && (t1 < closest_t))
                    {
                        closest_t = t1; mod = item as Plane;
                    }
                }
                if (item is Polygonal)
                {
                    (item as Polygonal).IntersectRayModel(origin, direction, ref t1, ref t2, ref normal);
                    if ((t1 >= t_min) && (t1 <= t_max) && (t1 < closest_t))
                    {
                        closest_t = t1; mod = item as Polygonal;
                    }
                }
            }
            return closest_t;
        }

        private static Coords ReflRayUnit(Coords incoming_unit, Coords normal_unit, bool ray_in)
        {
            normal_unit = ray_in ? -normal_unit : normal_unit;
            double cosPhi = normal_unit.ScalarProd(incoming_unit);
            cosPhi = (cosPhi > 1) ? 1 : (cosPhi < -1) ? -1 : cosPhi;
            Coords reflected_unit = cosPhi * 2 * normal_unit - incoming_unit;
            reflected_unit /= reflected_unit.Lenght();
            return reflected_unit;
        }

        private static double GetZOfSphereOrReturnNaN(double x, double y, double radius)
        {
            return Math.Sqrt(radius * radius - x * x - y * y);
        }

        private static ThreadList _thread_pool = new ThreadList();

        private static void AddSphericalLight(Coords center, int radius, int quantity_on_diameter)
        {
            double dl = 1 / (0.524 * quantity_on_diameter * quantity_on_diameter * quantity_on_diameter);
            for (int y = -radius; y <= radius; y += 2 * radius / quantity_on_diameter)//example:500,350,1000,rad:50,quantity:20
            {
                for (int x = -radius; x <= radius; x += 2 * radius / quantity_on_diameter)
                {
                    double k = GetZOfSphereOrReturnNaN(x, y, radius);
                    if (k != double.NaN)
                    {
                        Lights.Add(new Point(new Coords(x + center.X, y + center.Y, center.Z + k), dl));
                        Lights.Add(new Point(new Coords(x + center.X, y + center.Y, center.Z - k), dl));
                    }
                }
            }
        }

        public static void Main(string[] args)
        {
            using (var g = new GameWindow())
            {
                g.Load += (sender, e) =>
                {
                    g.VSync = VSyncMode.On;
                    g.Location = new System.Drawing.Point(0, 0);
                    g.Width = 1200;
                    g.Height = 748;
                };
                g.Resize += (sender, e) => { GL.Viewport(0, 0, g.Width, g.Height); };
                g.UpdateFrame += (sender, e) => { };
                g.RenderFrame += (sender, e) =>
                {
                    GL.ClearColor(Color.SkyBlue);
                    GL.Enable(EnableCap.DepthTest);
                    GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);
                    float distance = 500;
                    Matrix4 p = Matrix4.CreateOrthographic(g.Width, g.Height, 0, distance);
                    GL.MatrixMode(MatrixMode.Projection);
                    GL.LoadMatrix(ref p);
                    Coords camera = new Coords();
                    Matrix4 modelview = Matrix4.LookAt((float)camera.X, (float)camera.Y, (float)camera.Z, 0, 0, 1, 0, 1, 0);
                    GL.MatrixMode(MatrixMode.Modelview);
                    GL.LoadMatrix(ref modelview);

                    List<List<Coords>> Vertices = new List<List<Coords>>();
                    List<List<Vector>> Rays = new List<List<Vector>>();
                    List<Thread> ThreadSeq = new List<Thread>();
                    const int x_resolution = 1200;
                    const int num_threads = 34;
                    const int y_resolution = 748;
                    const int recurs_limit = 2;
                    double Xbegin = -g.Width / 2.0;
                    double Xstep = g.Width / (double)x_resolution;
                    double Ybegin = -g.Height / 2.0;
                    double Ystep = g.Height / (double)y_resolution;
                    for (int j = 0; j < y_resolution; j++)
                    {
                        int it = -1;
                        for (int i = 0; i < x_resolution; i++)
                        {
                            it = i;
                                //if ((i > 423) && (i < 425) && (j > 1140) && (j < 1142))
                                //if ((i > 423) && (i < 425) && (j > 1132) && (j < 1134))
                                //if ((i > 473) && (i < 475) && (j > 731) && (j < 733))
                                //if ((i > 482) && (i < 484) && (j > 731) && (j < 733))
                                //if ((i > 478) && (i < 480) && (j > 788) && (j < 790))
                                //if ((i > 418) && (i < 420) && (j >740) && (j < 742))
                                //if((i==510)&&(j==375))
                                {
                                List<Coords> tmp1 = new List<Coords>
                        {
                                new Coords(Xbegin + i * Xstep, Ybegin + j * Ystep, 1),
                                new Coords(Xbegin + (i + 1) * Xstep, Ybegin + j * Ystep, 1),
                                new Coords(Xbegin + i * Xstep, Ybegin + (j + 1) * Ystep, 1),
                                new Coords(Xbegin + (i + 1) * Xstep, Ybegin + (j + 1) * Ystep, 1)
                        };
                                List<Vector> tmp2 = new List<Vector>
                        {
                                new Vector(new Coords(Xbegin + (i+0.25) * Xstep, Ybegin + (j+0.25) * Ystep, 1),new Coords(0,0,1)),
                                new Vector(new Coords(Xbegin + (i+0.75) * Xstep, Ybegin + (j+0.25) * Ystep, 1),new Coords(0,0,1)),
                                new Vector(new Coords(Xbegin + (i+0.25) * Xstep, Ybegin + (j+0.75) * Ystep, 1),new Coords(0,0,1)),
                                new Vector(new Coords(Xbegin + (i+0.75) * Xstep, Ybegin + (j+0.75) * Ystep, 1),new Coords(0,0,1)),
                        };
                                Rays.Add(tmp2);
                                Vertices.Add(tmp1);
                            }
                        }
                        List<List<Vector>> Rays1 = Rays;
                        List<List<Coords>> Vertices1 = Vertices;
                        Proxy Tmp = new Proxy(Rays1, Vertices1, double.PositiveInfinity, recurs_limit, it, j);
                            //ThreadPool.QueueUserWorkItem(WrapTraceRay, Tmp);//T
                            Thread thread = new Thread(WrapTraceRay);
                        thread.Start(Tmp);
                        ThreadSeq.Add(thread);
                        Console.WriteLine(j + "\n");
                        Rays = new List<List<Vector>>();
                        Vertices = new List<List<Coords>>();
                        if (ThreadSeq.Count >= num_threads)
                        {
                            foreach (var item in ThreadSeq) { item.Join(); Console.WriteLine("======ENDED======\n"); }
                            ThreadSeq = new List<Thread>();
                            _thread_pool.DrawLst();
                            _thread_pool = new ThreadList();
                        }
                    }
                    g.SwapBuffers();
                        //Scene[0].Relocate(new Coords(-6.2533, 0, 0));
                        //Scene[2].Relocate(new Coords(8.6666, 0, 0));
                    };
                g.Run(60, 60);
            }
            return;
        }
    }
}
