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
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace Simplex
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        GeometryModel3D sphere = new GeometryModel3D(null, new DiffuseMaterial(Brushes.Red));

        public MainWindow()
        {
            InitializeComponent();

            IcoSphereCreator gen = new IcoSphereCreator();
            sphere.Geometry = gen.Create(2);
            group.Children.Add(sphere);
            GeometryModel3D othersphere = new GeometryModel3D(sphere.Geometry, new DiffuseMaterial(Brushes.Green));
            group.Children.Add(othersphere);
        }

        private void button_Click(object sender, RoutedEventArgs e)
        {
            sphere.Transform = new TranslateTransform3D(1, 1, 1);
            testWindow.InvalidateVisual();
        }
    }
}
