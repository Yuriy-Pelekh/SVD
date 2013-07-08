using MathNet.Numerics.LinearAlgebra.Double;
using Microsoft.Win32;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Windows;

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

      worker = new BackgroundWorker {WorkerReportsProgress = true};
      worker.DoWork += BackgroundWorkerDoWork;
      worker.RunWorkerCompleted += BackgroundWorkerRunWorkerCompleted;
      worker.ProgressChanged += BackgroundWorkerProgressChanged;
    }

    private readonly BackgroundWorker worker;
    private short[] data;
    private int vectorLength;

    private readonly List<SVD> svd = new List<SVD>();

    private void ButtonRunClick(object sender, RoutedEventArgs e)
    {
      UpdateUIButtonsState(false);

      svd.Clear();
      
      worker.RunWorkerAsync(data);
    }

    private void UpdateUIButtonsState(bool isEnabled)
    {
      ButtonRun.IsEnabled = isEnabled;
      ButtonOpen.IsEnabled = isEnabled;
    }

    private void ButtonOpenClick(object sender, RoutedEventArgs e)
    {
      var openFileDialog = new OpenFileDialog
                             {
                               Filter = "Audio Files (.wav)|*.wav|Text Files (.txt)|*.txt|All Files (*.*)|*.*",
                               RestoreDirectory = true
                             };

      if (openFileDialog.ShowDialog() == true)
      {
        try
        {
          var fileName = openFileDialog.FileName;
          var extension = fileName.Substring(fileName.LastIndexOf(".", StringComparison.Ordinal) + 1).ToLower();

          switch (extension)
          {
            case FileTypes.WAV:
              using (var stream = openFileDialog.OpenFile())
              {
                // Skip header
                stream.Seek(44, SeekOrigin.Begin);
                data = ReadBinaryFile(stream);
              }
              ButtonRun.IsEnabled = true;
              break;
            case FileTypes.TXT:
              using (var stream = openFileDialog.OpenFile())
              {
                data = ReadTextFile(stream);
              }
              ButtonRun.IsEnabled = true;
              break;
            default:
              TextBoxScreen.Text = "Not supported file format!";
              ButtonRun.IsEnabled = false;
              return;
          }

          TextBoxScreen.Text = string.Join(" ", data);
        }
        catch (Exception ex)
        {
          MessageBox.Show(ex.Message);
        }
      }
    }

    private void BackgroundWorkerDoWork(object sender, DoWorkEventArgs e)
    {
      var matrixArray = GetMatrix((short[]) e.Argument);

      foreach (var item in matrixArray)
      {
        var tempSvd = new SVD(item);
        var currentSvd = new DenseMatrix(tempSvd.A).Svd(true);

        tempSvd.U = currentSvd.U().ToArray();
        tempSvd.V = currentSvd.VT().ToArray();
        tempSvd.S = currentSvd.S().ToArray();

        svd.Add(tempSvd);

        worker.ReportProgress((int) ((double) svd.Count/matrixArray.Length*100));
      }
    }

    private void BackgroundWorkerRunWorkerCompleted(object sender, RunWorkerCompletedEventArgs e)
    {
      var epsilon = Convert.ToDouble(ComboBoxEpsilon.SelectedItem);

      var quasi = new QuasiStationary(svd.Select(item => item.S).ToList(), epsilon);
      quasi.Split();
      TextBoxScreen.Text = quasi.QuasiIndexesString;
      SaveResult(quasi.QuasiIndexes);

      ProgressBar.Value = 0;
      UpdateUIButtonsState(true);
    }

    private void BackgroundWorkerProgressChanged(object sender, ProgressChangedEventArgs e)
    {
      ProgressBar.Value = e.ProgressPercentage;
    }

    private double[][,] GetMatrix(short[] items)
    {
      try
      {
        var vectorsCount = items.Length/vectorLength;
        var matrix = new double[vectorsCount][,];

        for (var iCount = 0; iCount < vectorsCount; iCount++)
        {
          matrix[iCount] = new double[vectorLength,vectorLength];

          for (var i = 0; i < vectorLength; i++)
          {
            var baseItem = items[i*iCount];
            var delta = vectorLength*iCount;

            for (var j = 0; j < vectorLength; j++)
            {
              var temp = items[j + delta];
              matrix[iCount][i, j] = (double) temp/baseItem;
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

    private short[] GetPoints(IList<string> hexBytes)
    {
      var points = new List<short>();

      for (var i = 0; i < hexBytes.Count; i += 2)
      {
        points.Add(Convert.ToInt16(hexBytes[i] + hexBytes[i + 1], 16));
      }

      return points.ToArray();
    }

    private short[] ReadTextFile(Stream stream)
    {
      string dataString;
      using (var streamReader = new StreamReader(stream))
      {
        dataString = streamReader.ReadToEnd();
      }
      return dataString.Split(' ').Select(short.Parse).ToArray();
    }

    private short[] ReadBinaryFile(Stream stream)
    {
      byte[] loadedBytes;
      using (var buffer = new MemoryStream())
      {
        stream.CopyTo(buffer);
        loadedBytes = buffer.GetBuffer();
      }
      var dataString = BitConverter.ToString(loadedBytes);
      return GetPoints(dataString.Split('-'));
    }

    private void SaveResult(int[] ranges)
    {
      using (var sw = new StreamWriter("QuasiResult.txt"))
      {
        var items = data;

        var vectorsCount = items.Length/vectorLength;

        for (var i = 0; i < vectorsCount*vectorLength; i++)
        {
          var segmentIndex = i/vectorLength;

          if (ranges.Contains(segmentIndex))
          {
            if (i%vectorLength == 0)
            {
              sw.WriteLine("==================");
              sw.WriteLine("======= " + segmentIndex + " =======");
              sw.WriteLine("==================");
            }
          }
          else
          {
            if (i%vectorLength == 0)
            {
              sw.WriteLine("======= " + segmentIndex + " =======");
            }
          }

          sw.WriteLine(items[i]);
        }
      }
    }

    private void InitializeComboBox()
    {
      for (var iCount = 0; iCount <= 5000; iCount += 25)
      {
        ComboBoxEpsilon.Items.Add(iCount);
      }

      ComboBoxEpsilon.SelectedIndex = 30;
    }

    private void RadioButtonChecked(object sender, RoutedEventArgs e)
    {
      if (RadioButton60.IsChecked == true)
      {
        vectorLength = 60;
      }
      else if (RadioButton120.IsChecked == true)
      {
        vectorLength = 120;
      }
      else
      {
        vectorLength = 100;
      }
    }
  }
}
