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
    //if (Math.Abs(r2) < 0.01) throw new Exception("Particles " + one.number + " and " + two.number + " are too close together!");
    #region Potential definitions
    public abstract class BasePotential
    {
        public abstract double U(Particle one, Particle two);
        public abstract Vector3D delU(Particle one, Particle two);
    }

    public class Born : BasePotential
    {
        public double c { get; set; } = 1920.0; //au

        public override double U(Particle one, Particle two)
        {
            double rsquared = (one.Position - two.Position).LengthSquared;
            return (one.Charge * two.Charge / Math.Sqrt(rsquared)) + (c / Math.Pow(rsquared, 4));
        }

        public override Vector3D delU(Particle one, Particle two)
        {
            double rsquared = (one.Position - two.Position).LengthSquared;
            return -(one.Charge * two.Charge / Math.Pow(rsquared, 1.5) + 8 * c / Math.Pow(rsquared, 5)) * (two.Position - one.Position);
        }
    }

    public class BornAttractive : BasePotential
    {
        public double c { get; set; } = 1920.0; //au

        public override double U(Particle one, Particle two)
        {
            double rsquared = (one.Position - two.Position).LengthSquared;
            return (-Math.Abs(one.Charge * two.Charge) / Math.Sqrt(rsquared)) + (c / Math.Pow(rsquared, 4));
        }

        public override Vector3D delU(Particle one, Particle two)
        {
            double rsquared = (one.Position - two.Position).LengthSquared;
            return (Math.Abs(one.Charge * two.Charge) / Math.Pow(rsquared, 1.5) - 8 * c / Math.Pow(rsquared, 5)) * (two.Position - one.Position);
        }
    }

    public class LennardJones : BasePotential
    {
        public double epsilon { get; set; } = 1.728;      //mHartrees
        public double requilibrium { get; set; } = 7.049; //Bohrs

        public override double U(Particle one, Particle two)
        {
            double rsquared = (one.Position - two.Position).LengthSquared;
            double rRatio = requilibrium * requilibrium / rsquared;
            return epsilon * (Math.Pow(rRatio, 6) - 2 * Math.Pow(rRatio, 3));
        }

        public override Vector3D delU(Particle one, Particle two)
        {
            double rsquared = (one.Position - two.Position).LengthSquared;
            double rRatio = requilibrium * requilibrium / rsquared;
            return (-6 * epsilon / rsquared * (2 * Math.Pow(rRatio, 6) - Math.Pow(rRatio, 3))) * (two.Position - one.Position);
        }
    }

    public class Morse : BasePotential
    {
        public double De { get; set; } = -1.66;  //eV
        public double r0 { get; set; } = 2.47;   //Angstrom
        public double beta { get; set; } = 1.49; //Inverse Angstrom

        public override double U(Particle one, Particle two)
        {
            double rsquared = (one.Position - two.Position).LengthSquared;
            return De * Math.Pow(1 - Math.Exp(-beta * (rsquared - r0)), 2) - De; //TODO: check formula
        }

        public override Vector3D delU(Particle one, Particle two)
        {
            throw new NotImplementedException("Overeager bastard...these are difficult derivatives!");
        }
    }
    #endregion Potential definitions

    public class Particle : IComparable<Particle>
    {
        protected GeometryModel3D model = new GeometryModel3D();
        protected Vector3D posVector;
        protected Vector3D velVector;
        protected int number;
        protected int q;

        private static MeshGeometry3D sphereMesh = new IcoSphereCreator().Create(3); //Only create the mesh once, copy afterwards
        private static double velmult = 0.1;
        private static int index = 0;

        public Vector3D Position
        {
            get { return posVector; }
            set
            {
                Debug.WriteLine("Changing position of Particle " + number + ".");
                posVector = value;
                model.Transform = new TranslateTransform3D(posVector);
            }
        }
        public Vector3D Velocity
        {
            get { return velVector; }
            set
            {
                Console.WriteLine("Changing velocity of Particle " + number + ".");
                velVector = value;
                /*** TODO: revisit 'freeze' later
                if (velVector.LengthSquared < 0.0001)
                {
                    IsStatic = true;  //Once velocity gets low enough, freeze the particle.
                    Console.WriteLine("Velocity low, freezing Particle " + number + ".");
                }
                ***/
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
            Console.WriteLine("Creating Particle " + number + ".");
        }
        public Particle(Vector3D pos, int q)
        {
            Position = pos;
            Charge = q;
            number = ++index;
            Console.WriteLine("Creating Particle " + number + ".");
        }
        public Particle(Particle other) //Copy constructor
        {
            Position = other.Position;
            Charge = other.Charge;
            number = ++index;
            Console.WriteLine("Creating Particle " + number + " as clone of Particle " + other.number + ".");
        }

        ~Particle()
        {
            Console.WriteLine("Destroying Particle " + number + ".");
        }

        public int CompareTo(Particle other)
        {
            if (other == null) return 1;
            return Potential.CompareTo(other.Potential);
        }

        public void Show(Model3DGroup group)
        {
            Console.WriteLine("Displaying Particle " + number + "'s model.");
            if (model.Geometry == null) model.Geometry = sphereMesh;
            group.Children.Add(model);
        }

        public void Hide(Model3DGroup group)
        {
            Console.WriteLine("Hiding Particle " + number + "'s model.");
            group.Children.Remove(model);
        }

        public void updatePotential(List<Particle> others, BasePotential potential)
        {
            Potential = 0;

            foreach (Particle other in others)
            {
                if (other == null) throw new ArgumentNullException("Particle " + other.number);
                if (this == other) continue;

                Potential += potential.U(this, other);
            }
        }

        public void updateVelocity(List<Particle> others, BasePotential potential)
        {
            Velocity *= velmult; //Keep some of the prior velocity, makes optimization nonlinear

            foreach (Particle other in others)
            {
                if (other == null) throw new ArgumentNullException("Particle " + other.number);
                if (this == other) continue;

                Velocity += potential.delU(this, other);
            }
            Console.WriteLine("Particle " + number + " velocity: " + Velocity.ToString());
        }

        public void updatePosition(double stepSize)
        {
            Position += Velocity * stepSize;
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
        public void NewtonianStep(double stepSize = 1.0)
        {
            foreach (Particle part in pList)
            {
                part.updateVelocity(pList, potential);
            }
            //TODO: convergence criteria??

            pList.NewtonStep(stepSize);
        }


        //TODO: short-circuit logic for only two particles
        public void NMStep()
        {
            //NM coefficients, probably not going to be alterable.
            double alpha = 1.0;
            double gamma = 2.0;
            double rho = 0.5;
            double sigma = 0.5;

            //Calculate and sort by potential
            foreach (Particle part in pList)
            {
                part.updatePotential(pList, potential);
            }
            pList.Sort();

            //Calculate centroid point for all Particles but the last
            //TODO: exclude stationary particles (increase findCentroid argument)
            Vector3D centroid = findCentroid(1);
            Particle worst = new Particle(pList.Last());
            int last = pList.Count - 1;

            //reflection TODO: find why it's binding here
            Particle reflectedPoint = new Particle(centroid + alpha * (centroid - worst.Position), worst.Charge);
            reflectedPoint.updatePotential(pList, potential);

            if (reflectedPoint.Potential < pList[last - 1].Potential && reflectedPoint.Potential > pList[0].Potential)
            {
                pList[last].Hide(group);
                pList[last] = reflectedPoint;
                pList[last].Show(group);
                return;
            }

            //expansion
            if (reflectedPoint.Potential < pList[0].Potential)
            {
                Particle expandedPoint = new Particle(reflectedPoint.Position + gamma * (reflectedPoint.Position - centroid), worst.Charge);
                expandedPoint.updatePotential(pList, potential);

                pList[last].Hide(group);

                if (expandedPoint.Potential < reflectedPoint.Potential)
                    pList[last] = expandedPoint;
                else pList[last] = reflectedPoint;
                pList[last].Show(group);
                return;
            }

            //contraction
            Particle contractedPoint = new Particle(centroid + rho * (worst.Position - centroid), worst.Charge);
            contractedPoint.updatePotential(pList, potential);

            if (contractedPoint.Potential < worst.Potential)
            {
                pList[last].Hide(group);
                pList[last] = contractedPoint;
                pList[last].Show(group);
                return;
            }

            //reduction
            for (int i = 1; i < pList.Count; i++)
            {
                pList[i].Position = pList[0].Position + sigma * (pList[i].Position - pList[0].Position);
            }
        }

        public Vector3D findCentroid(int exclude = 1)
        {
            Vector3D center = new Vector3D(0, 0, 0);
            int n = pList.Count - exclude;
            
            for (int i = 0; i < n; i++)
            {
                center += pList[i].Position;
            }

            return (center / n);
        }
    }

    public static class ListExtensions
    {
        public static void ShowAll(this IEnumerable<Particle> pList, Model3DGroup group)
        {
            foreach (Particle part in pList) part.Show(group);
        }

        public static void NewtonStep(this IEnumerable<Particle> pList, double stepSize)
        {
            foreach (Particle part in pList) part.updatePosition(stepSize);
        }
    }
}
