using System;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Windows;
using Microsoft.Win32;
using SVD;

namespace QuasiStationar
{
  /// <summary>
  /// Interaction logic for MainWindow.xaml
  /// </summary>
  public partial class MainWindow
  {
    public MainWindow()
    {
      InitializeComponent();
      InitializeComboBox();
    }

    private const int VectorLength = 40;

    private string fileName;

    public struct SVD
    {
      public double[,] A;
      public double[,] U;
      public double[,] V;
      public double[] S;
    }

    public int[] Ranges { get; private set; }

    private readonly Collection<SVD> svd = new Collection<SVD>();
    private readonly Collection<double> inputData = new Collection<double>();

    private void ButtonRunClick(object sender, RoutedEventArgs e)
    {
      buttonRun.IsEnabled = false;
      buttonOpen.IsEnabled = false;

      progressBar.Value = 0;

      svd.Clear();

      if (fileName == null)
      {
        CreateTestFile();
        ShowFile();
      }

      var bw = new BackgroundWorker { WorkerReportsProgress = true };
      bw.DoWork += BwDoWork;
      bw.RunWorkerCompleted += BwRunWorkerCompleted;
      bw.ProgressChanged += BwProgressChanged;
      bw.RunWorkerAsync(fileName);
    }

    private void ButtonOpenClick(object sender, RoutedEventArgs e)
    {
      var openFileDialog = new OpenFileDialog
                               {
                                 Filter = "All files (*.*)|*.*",
                                 RestoreDirectory = true
                               };

      if (openFileDialog.ShowDialog().Value)
      {
        try
        {
          fileName = openFileDialog.FileName;
          ShowFile(fileName);
        }
        catch (Exception ex)
        {
          MessageBox.Show(ex.Message);
        }
      }
    }

    void BwDoWork(object sender, DoWorkEventArgs e)
    {
      var bw = sender as BackgroundWorker;

      var matrixCount = 0;

      var matrixArray = e.Argument == null ? GetMatrix() : GetMatrix((string)e.Argument);

      foreach (var item in matrixArray)
      {
        var tempSvd = new SVD { A = item };
        var currentSvd = new SingularValueDecomposition(tempSvd.A);

        tempSvd.U = currentSvd.U;
        tempSvd.V = currentSvd.V;
        tempSvd.S = currentSvd.S;

        svd.Add(tempSvd);

        //SaveMatrix("a" + matrixCount + ".txt", tempSvd.A);
        //SaveMatrix("u" + matrixCount + ".txt", tempSvd.U);
        //SaveMatrix("v" + matrixCount + ".txt", tempSvd.V);
        //SaveMatrix("s" + matrixCount + ".txt", tempSvd.S);

        matrixCount++;

        //Thread.Sleep(500);
        bw.ReportProgress((int)((double)matrixCount / (double)matrixArray.Length * 100));
      }
    }

    void BwRunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
    {
      buttonRun.IsEnabled = true;
      buttonOpen.IsEnabled = true;
      progressBar.Value = 0;

      var quasi = new QuasiStationarChecker(svd) { Epsilon = Convert.ToDouble(comboBox.SelectedItem) };
      textBox.Text = quasi.Check();
      Ranges = quasi.GetRanges();
      SaveResult();
      //MessageBox.Show(quasi.Check());
    }

    void BwProgressChanged(object sender, ProgressChangedEventArgs e)
    {
      progressBar.Value = e.ProgressPercentage;
    }

    private void CreateTestFile(string fileName = "Test.txt")
    {
      using (var stream = new StreamWriter(fileName))
      {
        var random = new Random();

        for (var i = 0; i < 10; i++)
        {
          for (var j = 0; j < 100; j++)
          {
            stream.Write(random.NextDouble().ToString().Substring(0, 5) + " ");
          }
          stream.WriteLine();
        }
      }
    }

    #region GetMatrix (separated vectors file format)
    private double[][,] GetMatrix()
    {
      const string fileName = "Test.txt";

      try
      {
        double[][,] matrix;

        using (var stream = new StreamReader(fileName))
        {
          var vectors = stream.ReadToEnd().Split('\n', '\r').Where(c => !string.IsNullOrEmpty(c)).ToArray();

          matrix = new double[vectors.Length][,];

          for (var vectorCount = 0; vectorCount < vectors.Length; vectorCount++)
          {
            var vector = vectors[vectorCount].Split(' ').Where(c => !string.IsNullOrEmpty(c)).ToArray();
            var length = vector.Length;

            matrix[vectorCount] = new double[length, length];

            for (var i = 0; i < length; i++)
            {
              var baseItem = double.Parse(vector[i]);

              for (var j = 0; j < length; j++)
              {
                var temp = double.Parse(vector[j]);
                matrix[vectorCount][i, j] = temp / baseItem;
              }
            }
          }
        }

        return matrix;
      }
      catch (Exception e)
      {
        MessageBox.Show(e.Message);
      }

      return null;
    }
    #endregion

    private double[][,] GetMatrix(string fileName)
    {
      try
      {
        double[][,] matrix;

        using (var stream = new StreamReader(fileName))
        {
          var items = stream.ReadToEnd().Split('\n', '\r', ' ').Where(c => !string.IsNullOrEmpty(c)).ToArray();

          var vectorsCount = items.Length / VectorLength / 10;

          //for (var i = 0; i < vectorsCount * VectorLength; i++)
          //{
          //  items[i] = items[i].Replace(',', '.');
          //}

          matrix = new double[vectorsCount][,];

          for (var iCount = 0; iCount < vectorsCount; iCount++)
          {
            matrix[iCount] = new double[VectorLength, VectorLength];

            for (var i = 0; i < VectorLength; i++)
            {
              var baseItem = double.Parse(items[i * iCount]);
              var delta = VectorLength * iCount;

              for (var j = 0; j < VectorLength; j++)
              {
                var temp = double.Parse(items[j + delta]);
                matrix[iCount][i, j] = temp / baseItem;
              }
            }
          }
        }

        return matrix;
      }
      catch (Exception e)
      {
        MessageBox.Show(e.Message);
      }

      return null;
    }

    private void SaveMatrix(string fileName, double[,] matrix)
    {
      try
      {
        using (var stream = new StreamWriter(fileName))
        {
          for (var i = 0; i < matrix.GetLength(0); i++)
          {
            for (var j = 0; j < matrix.GetLength(1); j++)
            {
              stream.Write(matrix[i, j] + "\t");
            }

            stream.WriteLine();
          }
        }
      }
      catch (Exception e)
      {
        MessageBox.Show(e.Message + " matrix");
      }
    }

    private void SaveMatrix(string fileName, double[] vector)
    {
      try
      {
        using (var stream = new StreamWriter(fileName))
        {
          for (var i = 0; i < vector.Length; i++)
          {
            stream.Write(vector[i] + "\t");
          }
        }
      }
      catch (Exception e)
      {
        MessageBox.Show(e.Message + " vector");
      }
    }

    private void SaveResult()
    {
      using (var sw = new StreamWriter("QuasiResult.txt"))
      {
        using (var stream = new StreamReader(fileName))
        {
          var items = stream.ReadToEnd().Split('\n', '\r', ' ').Where(c => !string.IsNullOrEmpty(c)).ToArray();

          var vectorsCount = items.Length/VectorLength/10;

          for (var i = 0; i < vectorsCount*VectorLength; i++)
          {
            items[i] = items[i].Replace(',', '.');

            var index = i/VectorLength;

            if (Ranges.Contains(index))
            {
              if (i % VectorLength == 0)
              {
                sw.WriteLine("==================");
                sw.WriteLine("======= " + index + " =======");
                sw.WriteLine("==================");
              }
            }

            sw.WriteLine(items[i]);
          }
        }
      }
    }

    private void ShowFile(string fileName = "Test.txt")
    {
      using (var stream = new StreamReader(fileName))
      {
        textBox.Text = stream.ReadToEnd();
      }
    }

    private void InitializeComboBox()
    {
      for (var iCount = 0; iCount < 1000; iCount += 25)
        comboBox.Items.Add(iCount);

      comboBox.SelectedValue = 600;
    }
  }
}
