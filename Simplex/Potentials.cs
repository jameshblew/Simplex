using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using System.Diagnostics;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace Simplex
{
    #region Potential definitions
    public abstract class BasePotential
    {
        public abstract double U(int q1, int q2, double rsquared);
    }

    public class Born : BasePotential
    {
        public double c { get; set; } = 1920.0; //au

        public override double U(int q1, int q2, double rsquared)
        {
            return (q1 * q2 / Math.Sqrt(rsquared)) + (c / Math.Pow(rsquared, 4));
        }
    }

    public class BornAttractive : BasePotential
    {
        public double c { get; set; } = 1920.0; //au

        public override double U(int q1, int q2, double rsquared)
        {
            return (-Math.Abs(q1 * q2) / Math.Sqrt(rsquared)) + (c / Math.Pow(rsquared, 4));
        }
    }

    public class LennardJones : BasePotential
    {
        public double epsilon { get; set; } = 1.728;      //mHartrees
        public double requilibrium { get; set; } = 7.049; //Bohrs

        public override double U(int q1, int q2, double rsquared)
        {
            double rRatio = requilibrium * requilibrium / rsquared;
            return epsilon * Math.Pow(rRatio, 6) - 2 * epsilon * Math.Pow(rRatio, 3);
        }
    }

    public class Morse : BasePotential
    {
        public double De { get; set; } = -1.66;  //eV
        public double r0 { get; set; } = 2.47;   //Angstrom
        public double beta { get; set; } = 1.49; //Inverse Angstrom

        public override double U(int q1, int q2, double rsquared)
        {
            return De * Math.Pow(1 - Math.Exp(-beta * (rsquared - r0)), 2) - De;
        }
    }
    #endregion Potential definitions

    public class Particle : IComparable<Particle>
    {
        protected GeometryModel3D model = new GeometryModel3D();
        protected Vector3D posVector;
        protected int number;
        protected int q;

        private static MeshGeometry3D sphereMesh = new IcoSphereCreator().Create(3); //Only create the mesh once, copy afterwards
        private static int index = 0;

        public Vector3D Position
        {
            get { return posVector; }
            set
            {
                Debug.WriteLine("Changing position of Particle " + number + ".");
                model.Transform = new TranslateTransform3D(value);
                posVector = value;
            }
        }
        public int Charge
        {
            get { return q; }
            set
            {
                Debug.WriteLine("Changing charge of Particle " + number + ".");
                q = value;
                if (q == 0)
                    model.Material = new DiffuseMaterial(Brushes.Green);
                else if (q < 0)
                    model.Material = new DiffuseMaterial(Brushes.Blue);
                else if (q > 0)
                    model.Material = new DiffuseMaterial(Brushes.Red);
                else
                    throw new ArithmeticException("WTF Error: Charge is not correctly interpreted.");
            }
        }
        public bool IsStatic { get; set; }
        public double Potential { get; set; } = 0;

        public Particle(double x, double y, double z, int q)
        {
            Position = new Vector3D(x, y, z);
            Charge = q;
            number = ++index;
            Console.WriteLine("Creating Particle number " + number + ".");
        }
        
        public int CompareTo(Particle other)
        {
            if (other == null) return 1;

            return Potential.CompareTo(other.Potential);
        }

        public void Show(Model3DGroup group)
        {
            Console.WriteLine("Creating Particle model.");
            model.Geometry = sphereMesh;
            group.Children.Add(model);
        }

        public double calculatePotential(Particle other, BasePotential potential)
        {
            if (other == null) throw new NullReferenceException();
            if (this == other) return 0;

            double r2 = (Position - other.Position).LengthSquared;
            if (Math.Abs(r2) < 0.01) throw new Exception("Particles too close together!");

            return potential.U(Charge, other.Charge, r2);
        }
    }

    //Following code taken from the blog "catch 22" by Andreas Kahler.
    //http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
    //Slight modifications have been made.
    public class IcoSphereCreator
    {
        private struct TriangleIndices
        {
            public int v1;
            public int v2;
            public int v3;

            public TriangleIndices(int v1, int v2, int v3)
            {
                this.v1 = v1;
                this.v2 = v2;
                this.v3 = v3;
            }
        }

        //System.IO.StreamWriter file = new System.IO.StreamWriter(@"C:\Users\James\Desktop\mesh.txt");

        private MeshGeometry3D geometry;
        private int index;
        private Dictionary<long, int> middlePointIndexCache;

        // add vertex to mesh, fix position to be on unit sphere, return index
        private int addVertex(Point3D p)
        {
            double length = Math.Sqrt(p.X * p.X + p.Y * p.Y + p.Z * p.Z);
            geometry.Positions.Add(new Point3D(p.X / length, p.Y / length, p.Z / length));
            return index++;
        }

        // return index of point in the middle of p1 and p2
        private int getMiddlePoint(int p1, int p2)
        {
            // first check if we have it already
            bool firstIsSmaller = p1 < p2;
            long smallerIndex = firstIsSmaller ? p1 : p2;
            long greaterIndex = firstIsSmaller ? p2 : p1;
            long key = (smallerIndex << 32) + greaterIndex;

            int ret;
            if (middlePointIndexCache.TryGetValue(key, out ret))
            {
                return ret;
            }

            // not in cache, calculate it
            Point3D point1 = geometry.Positions[p1];
            Point3D point2 = geometry.Positions[p2];
            Point3D middle = new Point3D(
                (point1.X + point2.X) / 2.0,
                (point1.Y + point2.Y) / 2.0,
                (point1.Z + point2.Z) / 2.0);

            // add vertex makes sure point is on unit sphere
            int i = addVertex(middle);

            // store it, return index
            middlePointIndexCache.Add(key, i);
            return i;
        }

        public MeshGeometry3D Create(int recursionLevel)
        {
            geometry = new MeshGeometry3D();
            middlePointIndexCache = new Dictionary<long, int>();
            index = 0;


            // create 12 vertices of a icosahedron
            var t = (1.0 + Math.Sqrt(5.0)) / 2.0;

            addVertex(new Point3D(-1, t, 0));
            addVertex(new Point3D(1, t, 0));
            addVertex(new Point3D(-1, -t, 0));
            addVertex(new Point3D(1, -t, 0));

            addVertex(new Point3D(0, -1, t));
            addVertex(new Point3D(0, 1, t));
            addVertex(new Point3D(0, -1, -t));
            addVertex(new Point3D(0, 1, -t));

            addVertex(new Point3D(t, 0, -1));
            addVertex(new Point3D(t, 0, 1));
            addVertex(new Point3D(-t, 0, -1));
            addVertex(new Point3D(-t, 0, 1));


            // create 20 triangles of the icosahedron
            var faces = new List<TriangleIndices>();

            // 5 faces around point 0
            faces.Add(new TriangleIndices(0, 11, 5));
            faces.Add(new TriangleIndices(0, 5, 1));
            faces.Add(new TriangleIndices(0, 1, 7));
            faces.Add(new TriangleIndices(0, 7, 10));
            faces.Add(new TriangleIndices(0, 10, 11));

            // 5 adjacent faces 
            faces.Add(new TriangleIndices(1, 5, 9));
            faces.Add(new TriangleIndices(5, 11, 4));
            faces.Add(new TriangleIndices(11, 10, 2));
            faces.Add(new TriangleIndices(10, 7, 6));
            faces.Add(new TriangleIndices(7, 1, 8));

            // 5 faces around point 3
            faces.Add(new TriangleIndices(3, 9, 4));
            faces.Add(new TriangleIndices(3, 4, 2));
            faces.Add(new TriangleIndices(3, 2, 6));
            faces.Add(new TriangleIndices(3, 6, 8));
            faces.Add(new TriangleIndices(3, 8, 9));

            // 5 adjacent faces 
            faces.Add(new TriangleIndices(4, 9, 5));
            faces.Add(new TriangleIndices(2, 4, 11));
            faces.Add(new TriangleIndices(6, 2, 10));
            faces.Add(new TriangleIndices(8, 6, 7));
            faces.Add(new TriangleIndices(9, 8, 1));


            // refine triangles
            for (int i = 0; i < recursionLevel; i++)
            {
                var faces2 = new List<TriangleIndices>();
                foreach (var tri in faces)
                {
                    // replace triangle by 4 triangles
                    int a = getMiddlePoint(tri.v1, tri.v2);
                    int b = getMiddlePoint(tri.v2, tri.v3);
                    int c = getMiddlePoint(tri.v3, tri.v1);

                    faces2.Add(new TriangleIndices(tri.v1, a, c));
                    faces2.Add(new TriangleIndices(tri.v2, b, a));
                    faces2.Add(new TriangleIndices(tri.v3, c, b));
                    faces2.Add(new TriangleIndices(a, b, c));
                }
                faces = faces2;
            }

            // done, now add triangles to mesh
            foreach (var tri in faces)
            {
                geometry.TriangleIndices.Add(tri.v1);
                geometry.TriangleIndices.Add(tri.v2);
                geometry.TriangleIndices.Add(tri.v3);
            }

            Console.WriteLine("Generating points");
            /**
            file.WriteLine(geometry.Positions.Count);
            foreach (Point3D point in geometry.Positions)
            {
                file.WriteLine("geometry.Positions.Add(new Point3D(" + point.ToString() + "));");
            }
            file.WriteLine(geometry.TriangleIndices.Count);
            foreach (int vert in geometry.TriangleIndices)
            {
                file.WriteLine("geometry.TriangleIndices.Add(" + vert.ToString() + ");");
            }
            **/
            return geometry;
        }
    }

    /// <summary>
    /// Application logic for MainWindow.xaml
    /// </summary>
    partial class MainWindow : Window
    {
        //NM coefficients, probably not going to be alterable.
        private double alpha = 1.0;
        private double gamma = 2.0;
        private double rho = 0.5;
        private double sigma = 0.5;

        public List<Particle> NMStep(List<Particle> pList, BasePotential potential)
        {
            //TODO: Nelder Mead step process
            // exclude stationary particles

            //TODO: calculate potentials

            //sort by potential
            pList.Sort();

            //TODO: calculate centroid point for all but last
            //TODO: reflection
            //TODO: expansion
            //TODO: contraction
            //TODO: reduction

            throw new NotImplementedException("Too soon, bro.");
        }

        public List<Particle> calcPotentials(List<Particle> pList, BasePotential potential)
        {
            //TODO: calculate the potentials
            int n = pList.Count;
            double placeholder;
            double[,] U = new double[n, n];

            //Form potentials matrix (symmetric, real, zero-diagonal)
            for (int i = 0; i < n; i++)
            {
                for (int j = i; i < n; j++)
                {
                    placeholder = pList[i].calculatePotential(pList[j], potential);
                    U[i, j] = placeholder;
                    U[j, i] = placeholder;
                }
            }

            //Sum the rows (or columns) for individual center potentials
            for (int i = 0; i < n; i++)
            {
                pList[i].Potential = 0;
                for (int j = 0; j < n; j++)
                {
                    pList[i].Potential += U[i, j];
                }
            }

            return pList;
        }
    }
}
