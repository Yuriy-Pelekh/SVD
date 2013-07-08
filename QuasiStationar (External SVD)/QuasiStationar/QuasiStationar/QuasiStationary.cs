using System;
using System.Collections.Generic;
using System.Linq;

namespace QuasiStationar
{
  public class QuasiStationary
  {
    private readonly List<double[]> vectors;
    private readonly double epsilon;

    public int[] QuasiIndexes { get; private set; }
    public string QuasiIndexesString { get; private set; }

    public QuasiStationary(List<double[]> vectors, double epsilon)
    {
      this.vectors = vectors;
      this.epsilon = epsilon;
    }

    private double GetQuasiSegmentValue(double[] svd1, double[] svd2)
    {
      return Math.Abs(svd1.Length != svd2.Length
                        ? double.NaN
                        : Math.Sqrt(svd1.Select((t, i) => (Math.Pow(Math.Abs(t - svd2[i]), 2))).Sum()));
    }

    public List<List<int>> Split()
    {
      var result = new List<List<int>> {new List<int> {0}};

      var currentQuasiIndex = 0;

      for (var index = 1; index < vectors.Count; index++)
      {
        var quasiSegmentValue = GetQuasiSegmentValue(vectors[currentQuasiIndex], vectors[index]);

        if (!(quasiSegmentValue < epsilon))
        {
          currentQuasiIndex = index;
          result.Add(new List<int>());
        }

        result.Last().Add(index);
      }

      QuasiIndexesString = string.Join(Environment.NewLine, result.Select(s => string.Join(" + ", s)));
      QuasiIndexes = result.Select(v => v[0]).ToArray();

      return result;
    }
  }
}
