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
        Particle testSphere;

        public MainWindow()
        {
            InitializeComponent();

            testSphere = new Particle(-1, 0, 0, 1);
            testSphere.Show(group);
        }

        private void button_Click(object sender, RoutedEventArgs e)
        {
            Particle sp2 = new Particle(0, 0, 0, -1);
            sp2.Show(group);
            testSphere.Position = new Vector3D(1, 0, 0);
            testWindow.InvalidateVisual();
        }
    }
}
