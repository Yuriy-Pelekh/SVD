using System;
using System.Collections.ObjectModel;
using System.Linq;
using SingularValueDecomposition = QuasiStationar.MainWindow.SVD;

namespace QuasiStationar
{
  class QuasiStationarChecker
  {
    public QuasiStationarChecker(Collection<SingularValueDecomposition> svd)
    {
      this.svd = svd;
    }

    private readonly Collection<SingularValueDecomposition> svd;

    private double epsilon = 600;
    public double Epsilon
    {
      get { return epsilon; }
      set { epsilon = value; }
    }

    private double Check(SingularValueDecomposition svd1, SingularValueDecomposition svd2)
    {
      return Math.Abs(svd1.S.Length != svd2.S.Length ? double.NaN : svd1.S.Select((t, i) => (t - svd2.S[i])).Sum());
    }

    private string rangesString;
    private int[] rangesArray;

    public string Check()
    {
      var result = "0";

      for (var i = 1; i < svd.Count; i++)
      {
        var res = Check(svd[i - 1], svd[i]);

        if (res < epsilon)
        {
          result += " + ";// "(" + res + ") ";
        }
        else
        {
          result += "\n";// "(" + res + ") ";
        }

        result += i;
      }

      rangesString = result;

      return result;
    }

    public int[] GetRanges()
    {
      if (rangesString == null)
        Check();

      var ranges = rangesString.Split('\n');

      rangesArray = ranges.Select(range => Convert.ToInt32(range.Split('+').Last())).ToArray();

      return rangesArray;
    }
  }
}
