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
        List<Particle> pList = new List<Particle>();
        BasePotential potential;

        public MainWindow()
        {
            InitializeComponent();

            pList.Add(new Particle(-2, -2, 0, -1));
            pList.Add(new Particle(2, -2, 0, 1));
            pList.Add(new Particle(0, 3, 0, 1));
            pList.ShowAll(group);

            potential = new BornAttractive();
        }

        private void button_Click(object sender, RoutedEventArgs e)
        {
            //while (true)
            {
                NMStep(pList, potential);
                //await Task.Delay(50);
                testWindow.InvalidateVisual();
            }
        }
    }
}
