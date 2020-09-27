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
    abstract class Light
    {
        public virtual double Power { get; protected set; }
    }

    class Point : Light //OK
    {
        public Coords Position { get; protected set; }
        public override double Power { get; protected set; }
        public Point(Coords pos, double pow)
        {
            Position = pos ?? throw new ArgumentNullException(paramName: nameof(pos));
            Power = pow;
        }
    }

    class Directed : Light //OK
    {
        public Coords Direction { get; protected set; }
        public override double Power { get; protected set; }
        public Directed(Coords dir, double pow)
        {
            Direction = dir ?? throw new ArgumentNullException(paramName: nameof(dir));
            Power = pow;
        }
    }

    class Ambient : Light //OK
    {
        public override double Power { get; protected set; }
        public Ambient(double pow)
        {
            Power = pow;
        }
    }
}
