namespace QuasiStationar
{
  public class SVD
  {
    /// <summary>
    /// Source matrix.
    /// </summary>
    public double[,] A;
    
    /// <summary>
    /// SVD matrix 1.
    /// </summary>
    public double[,] U;
    
    /// <summary>
    /// SVD matrix 2.
    /// </summary>
    public double[,] V;
    
    /// <summary>
    ///  Diagonal matrix - vector.
    /// </summary>
    public double[] S;

    public SVD(double[,] a)
    {
      A = a;
    }
  }
}
