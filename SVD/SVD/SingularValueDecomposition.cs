using System;
using System.Collections.Generic;

namespace SVD
{
    public sealed class SingularValueDecomposition
    {
        #region Fields

        private int _mRow;
        private int _nCol;

        #endregion

        #region Properties

        #region Never used

        //public double Condition
        //{
        //    get { return S[0] / S[Math.Min(_mRow, _nCol) - 1]; }
        //}

        //public double TwoNorm
        //{
        //    get { return S[0]; }
        //}

        //public int Rank
        //{
        //    get
        //    {
        //        var eps = Math.Pow(2.0, -52.0);
        //        var tol = Math.Max(_mRow, _nCol) * S[0] * eps;

        //        #region Same as LiNQ-expession

        //        //int rank = 0;

        //        //for (int i = 0; i < _s.Length; i++)
        //        //{
        //        //    if (_s[i] > tol)
        //        //    {
        //        //        rank++;
        //        //    }
        //        //}

        //        //return rank;

        //        #endregion

        //        return S.Count(t => t > tol);

        //    }
        //}

        #endregion
        
        // A = U * S * V
        
        public double[,] A { get; private set; }

        public double[] S { get; private set; }

        public double[,] V { get; private set; }

        public double[,] U { get; private set; }

        public double Threshold
        {
            get { return Double.Epsilon * Math.Max(_mRow, _nCol) * S[0]; }
        }

        #endregion

        #region Constructors

        public SingularValueDecomposition(double[,] value)
        {
            Run(value, true, true);
        }

        #region Old

        //public SingularValueDecomposition(double[,] value)
        //    : this(value, true, true)
        //{
        //}

        //public SingularValueDecomposition(double[,] value, bool computeLeftSingularVectors, bool computeRightSingularVectors)
        //{
        //    double[,] a = (double[,])value.Clone();
        //    _mRow = value.GetLength(0); // rows
        //    _nCol = value.GetLength(1); // cols
        //    int nu = Math.Min(_mRow, _nCol);
        //    _s = new double[Math.Min(_mRow + 1, _nCol)];
        //    _u = new double[_mRow, nu];
        //    _v = new double[_nCol, _nCol];
        //    double[] e = new double[_nCol];
        //    double[] work = new double[_mRow];
        //    bool wantu = computeLeftSingularVectors;
        //    bool wantv = computeRightSingularVectors;

        //    // Reduce A to bidiagonal form, storing the diagonal elements in s and the super-diagonal elements in e.
        //    int nct = Math.Min(_mRow - 1, _nCol);
        //    int nrt = Math.Max(0, Math.Min(_nCol - 2, _mRow));
        //    for (int k = 0; k < Math.Max(nct, nrt); k++)
        //    {
        //        if (k < nct)
        //        {
        //            // Compute the transformation for the k-th column and place the k-th diagonal in s[k].
        //            // Compute 2-norm of k-th column without under/overflow.
        //            _s[k] = 0;
        //            for (int i = k; i < _mRow; i++)
        //            {
        //                _s[k] = Hypotenuse(_s[k], a[i, k]);
        //            }

        //            if (_s[k] != 0.0)
        //            {
        //                if (a[k, k] < 0.0)
        //                {
        //                    _s[k] = -_s[k];
        //                }

        //                for (int i = k; i < _mRow; i++)
        //                {
        //                    a[i, k] /= _s[k];
        //                }

        //                a[k, k] += 1.0;
        //            }

        //            _s[k] = -_s[k];
        //        }

        //        for (int j = k + 1; j < _nCol; j++)
        //        {
        //            if ((k < nct) & (_s[k] != 0.0))
        //            {
        //                // Apply the transformation.
        //                double t = 0;
        //                for (int i = k; i < _mRow; i++)
        //                    t += a[i, k] * a[i, j];
        //                t = -t / a[k, k];
        //                for (int i = k; i < _mRow; i++)
        //                    a[i, j] += t * a[i, k];
        //            }

        //            // Place the k-th row of A into e for the subsequent calculation of the row transformation.
        //            e[j] = a[k, j];
        //        }

        //        if (wantu & (k < nct))
        //        {
        //            // Place the transformation in U for subsequent back
        //            // multiplication.
        //            for (int i = k; i < _mRow; i++)
        //                _u[i, k] = a[i, k];
        //        }

        //        if (k < nrt)
        //        {
        //            // Compute the k-th row transformation and place the k-th super-diagonal in e[k].
        //            // Compute 2-norm without under/overflow.
        //            e[k] = 0;
        //            for (int i = k + 1; i < _nCol; i++)
        //            {
        //                e[k] = Hypotenuse(e[k], e[i]);
        //            }

        //            if (e[k] != 0.0)
        //            {
        //                if (e[k + 1] < 0.0)
        //                    e[k] = -e[k];

        //                for (int i = k + 1; i < _nCol; i++)
        //                    e[i] /= e[k];

        //                e[k + 1] += 1.0;
        //            }

        //            e[k] = -e[k];
        //            if ((k + 1 < _mRow) & (e[k] != 0.0))
        //            {
        //                // Apply the transformation.
        //                for (int i = k + 1; i < _mRow; i++)
        //                    work[i] = 0.0;

        //                for (int j = k + 1; j < _nCol; j++)
        //                    for (int i = k + 1; i < _mRow; i++)
        //                        work[i] += e[j] * a[i, j];

        //                for (int j = k + 1; j < _nCol; j++)
        //                {
        //                    double t = -e[j] / e[k + 1];
        //                    for (int i = k + 1; i < _mRow; i++)
        //                        a[i, j] += t * work[i];
        //                }
        //            }

        //            if (wantv)
        //            {
        //                // Place the transformation in V for subsequent back multiplication.
        //                for (int i = k + 1; i < _nCol; i++)
        //                    _v[i, k] = e[i];
        //            }
        //        }
        //    }

        //    // Set up the final bidiagonal matrix or order p.
        //    int p = Math.Min(_nCol, _mRow + 1);
        //    if (nct < _nCol) _s[nct] = a[nct, nct];
        //    if (_mRow < p) _s[p - 1] = 0.0;
        //    if (nrt + 1 < p) e[nrt] = a[nrt, p - 1];
        //    e[p - 1] = 0.0;

        //    // If required, generate U.
        //    if (wantu)
        //    {
        //        for (int j = nct; j < nu; j++)
        //        {
        //            for (int i = 0; i < _mRow; i++)
        //                _u[i, j] = 0.0;
        //            _u[j, j] = 1.0;
        //        }

        //        for (int k = nct - 1; k >= 0; k--)
        //        {
        //            if (_s[k] != 0.0)
        //            {
        //                for (int j = k + 1; j < nu; j++)
        //                {
        //                    double t = 0;
        //                    for (int i = k; i < _mRow; i++)
        //                        t += _u[i, k] * _u[i, j];

        //                    t = -t / _u[k, k];
        //                    for (int i = k; i < _mRow; i++)
        //                        _u[i, j] += t * _u[i, k];
        //                }

        //                for (int i = k; i < _mRow; i++)
        //                    _u[i, k] = -_u[i, k];

        //                _u[k, k] = 1.0 + _u[k, k];
        //                for (int i = 0; i < k - 1; i++)
        //                    _u[i, k] = 0.0;
        //            }
        //            else
        //            {
        //                for (int i = 0; i < _mRow; i++)
        //                    _u[i, k] = 0.0;
        //                _u[k, k] = 1.0;
        //            }
        //        }
        //    }

        //    // If required, generate V.
        //    if (wantv)
        //    {
        //        for (int k = _nCol - 1; k >= 0; k--)
        //        {
        //            if ((k < nrt) & (e[k] != 0.0))
        //            {
        //                //TODO: Check if this is a bug.
        //                //   for (int j = k + 1; j < nu; j++)
        //                // The correction would be:
        //                for (int j = k + 1; j < _nCol; j++)
        //                {
        //                    double t = 0;
        //                    for (int i = k + 1; i < _nCol; i++)
        //                        t += _v[i, k] * _v[i, j];

        //                    t = -t / _v[k + 1, k];
        //                    for (int i = k + 1; i < _nCol; i++)
        //                        _v[i, j] += t * _v[i, k];
        //                }
        //            }

        //            for (int i = 0; i < _nCol; i++)
        //                _v[i, k] = 0.0;
        //            _v[k, k] = 1.0;
        //        }
        //    }

        //    // Main iteration loop for the singular values.
        //    int pp = p - 1;
        //    int iter = 0;
        //    double eps = Math.Pow(2.0, -52.0);
        //    while (p > 0)
        //    {
        //        int k, kase;

        //        // Here is where a test for too many iterations would go.
        //        // This section of the program inspects for
        //        // negligible elements in the s and e arrays.  On
        //        // completion the variables kase and k are set as follows.
        //        // kase = 1     if s(p) and e[k-1] are negligible and k<p
        //        // kase = 2     if s(k) is negligible and k<p
        //        // kase = 3     if e[k-1] is negligible, k<p, and s(k), ..., s(p) are not negligible (qr step).
        //        // kase = 4     if e(p-1) is negligible (convergence).
        //        for (k = p - 2; k >= -1; k--)
        //        {
        //            if (k == -1)
        //                break;

        //            if (Math.Abs(e[k]) <= eps * (Math.Abs(_s[k]) + Math.Abs(_s[k + 1])))
        //            {
        //                e[k] = 0.0;
        //                break;
        //            }
        //        }

        //        if (k == p - 2)
        //        {
        //            kase = 4;
        //        }
        //        else
        //        {
        //            int ks;
        //            for (ks = p - 1; ks >= k; ks--)
        //            {
        //                if (ks == k)
        //                    break;

        //                double t = (ks != p ? Math.Abs(e[ks]) : 0.0) + (ks != k + 1 ? Math.Abs(e[ks - 1]) : 0.0);
        //                if (Math.Abs(_s[ks]) <= eps * t)
        //                {
        //                    _s[ks] = 0.0;
        //                    break;
        //                }
        //            }

        //            if (ks == k)
        //                kase = 3;
        //            else if (ks == p - 1)
        //                kase = 1;
        //            else
        //            {
        //                kase = 2;
        //                k = ks;
        //            }
        //        }

        //        k++;

        //        // Perform the task indicated by kase.
        //        switch (kase)
        //        {
        //            // Deflate negligible s(p).
        //            case 1:
        //                {
        //                    double f = e[p - 2];
        //                    e[p - 2] = 0.0;
        //                    for (int j = p - 2; j >= k; j--)
        //                    {
        //                        double t = Hypotenuse(_s[j], f);
        //                        double cs = _s[j] / t;
        //                        double sn = f / t;
        //                        _s[j] = t;
        //                        if (j != k)
        //                        {
        //                            f = -sn * e[j - 1];
        //                            e[j - 1] = cs * e[j - 1];
        //                        }

        //                        if (wantv)
        //                        {
        //                            for (int i = 0; i < _nCol; i++)
        //                            {
        //                                t = cs * _v[i, j] + sn * _v[i, p - 1];
        //                                _v[i, p - 1] = -sn * _v[i, j] + cs * _v[i, p - 1];
        //                                _v[i, j] = t;
        //                            }
        //                        }
        //                    }
        //                }
        //                break;

        //            // Split at negligible s(k).
        //            case 2:
        //                {
        //                    double f = e[k - 1];
        //                    e[k - 1] = 0.0;
        //                    for (int j = k; j < p; j++)
        //                    {
        //                        double t = Hypotenuse(_s[j], f);
        //                        double cs = _s[j] / t;
        //                        double sn = f / t;
        //                        _s[j] = t;
        //                        f = -sn * e[j];
        //                        e[j] = cs * e[j];
        //                        if (wantu)
        //                        {
        //                            for (int i = 0; i < _mRow; i++)
        //                            {
        //                                t = cs * _u[i, j] + sn * _u[i, k - 1];
        //                                _u[i, k - 1] = -sn * _u[i, j] + cs * _u[i, k - 1];
        //                                _u[i, j] = t;
        //                            }
        //                        }
        //                    }
        //                }
        //                break;

        //            // Perform one qr step.
        //            case 3:
        //                {
        //                    // Calculate the shift.
        //                    double scale = Math.Max(Math.Max(Math.Max(Math.Max(Math.Abs(_s[p - 1]), Math.Abs(_s[p - 2])), Math.Abs(e[p - 2])), Math.Abs(_s[k])), Math.Abs(e[k]));
        //                    double sp = _s[p - 1] / scale;
        //                    double spm1 = _s[p - 2] / scale;
        //                    double epm1 = e[p - 2] / scale;
        //                    double sk = _s[k] / scale;
        //                    double ek = e[k] / scale;
        //                    double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
        //                    double c = (sp * epm1) * (sp * epm1);
        //                    double shift = 0.0;
        //                    if ((b != 0.0) | (c != 0.0))
        //                    {
        //                        shift = Math.Sqrt(b * b + c);
        //                        if (b < 0.0)
        //                            shift = -shift;
        //                        shift = c / (b + shift);
        //                    }

        //                    double f = (sk + sp) * (sk - sp) + shift;
        //                    double g = sk * ek;

        //                    // Chase zeros.
        //                    for (int j = k; j < p - 1; j++)
        //                    {
        //                        double t = Hypotenuse(f, g);
        //                        double cs = f / t;
        //                        double sn = g / t;
        //                        if (j != k)
        //                            e[j - 1] = t;
        //                        f = cs * _s[j] + sn * e[j];
        //                        e[j] = cs * e[j] - sn * _s[j];
        //                        g = sn * _s[j + 1];
        //                        _s[j + 1] = cs * _s[j + 1];
        //                        if (wantv)
        //                        {
        //                            for (int i = 0; i < _nCol; i++)
        //                            {
        //                                t = cs * _v[i, j] + sn * _v[i, j + 1];
        //                                _v[i, j + 1] = -sn * _v[i, j] + cs * _v[i, j + 1];
        //                                _v[i, j] = t;
        //                            }
        //                        }

        //                        t = Hypotenuse(f, g);
        //                        cs = f / t;
        //                        sn = g / t;
        //                        _s[j] = t;
        //                        f = cs * e[j] + sn * _s[j + 1];
        //                        _s[j + 1] = -sn * e[j] + cs * _s[j + 1];
        //                        g = sn * e[j + 1];
        //                        e[j + 1] = cs * e[j + 1];
        //                        if (wantu && (j < _mRow - 1))
        //                        {
        //                            for (int i = 0; i < _mRow; i++)
        //                            {
        //                                t = cs * _u[i, j] + sn * _u[i, j + 1];
        //                                _u[i, j + 1] = -sn * _u[i, j] + cs * _u[i, j + 1];
        //                                _u[i, j] = t;
        //                            }
        //                        }
        //                    }

        //                    e[p - 2] = f;
        //                    iter = iter + 1;
        //                }
        //                break;

        //            // Convergence.
        //            case 4:
        //                {
        //                    // Make the singular values positive.
        //                    if (_s[k] <= 0.0)
        //                    {
        //                        _s[k] = (_s[k] < 0.0 ? -_s[k] : 0.0);
        //                        if (wantv)
        //                            for (int i = 0; i <= pp; i++)
        //                                _v[i, k] = -_v[i, k];
        //                    }

        //                    // Order the singular values.
        //                    while (k < pp)
        //                    {
        //                        if (_s[k] >= _s[k + 1])
        //                            break;

        //                        double t = _s[k];
        //                        _s[k] = _s[k + 1];
        //                        _s[k + 1] = t;
        //                        if (wantv && (k < _nCol - 1))
        //                            for (int i = 0; i < _nCol; i++)
        //                            {
        //                                t = _v[i, k + 1];
        //                                _v[i, k + 1] = _v[i, k];
        //                                _v[i, k] = t;
        //                            }

        //                        if (wantu && (k < _mRow - 1))
        //                            for (int i = 0; i < _mRow; i++)
        //                            {
        //                                t = _u[i, k + 1];
        //                                _u[i, k + 1] = _u[i, k];
        //                                _u[i, k] = t;
        //                            }

        //                        k++;
        //                    }

        //                    iter = 0;
        //                    p--;
        //                }
        //                break;
        //        }
        //    }
        //}

        #endregion

        #endregion

        #region Solution

        #region Never used

        //public double[,] Solve(double[,] value)
        //{
        //    // Additionally an important property is that if there does not exists a solution
        //    // when the matrix A is singular but replacing 1/Li with 0 will provide a solution
        //    // that minimizes the residue |AX -Y|. SVD finds the least squares best compromise
        //    // solution of the linear equation system. Interestingly SVD can be also used in an
        //    // over-determined system where the number of equations exceeds that of the parameters.

        //    // L is a diagonal matrix with non-negative matrix elements having the same
        //    // dimension as A, Wi ? 0. The diagonal elements of L are the singular values of matrix A.

        //    var y = value;

        //    // Create L*, which is a diagonal matrix with elements
        //    //    L*[i] = 1/L[i]  if L[i] < e, else 0, 
        //    // where e is the so-called singularity threshold.

        //    // In other words, if L[i] is zero or close to zero (smaller than e),
        //    // one must replace 1/L[i] with 0. The value of e depends on the precision
        //    // of the hardware. This method can be used to solve linear equations
        //    // systems even if the matrices are singular or close to singular.

        //    //singularity threshold
        //    var e = Threshold;

        //    var ls = new double[S.Length, S.Length];
        //    for (var i = 0; i < S.Length; i++)
        //    {
        //        if (Math.Abs(S[i]) < e)
        //        {
        //            ls[i, i] = 0.0;
        //        }
        //        else
        //        {
        //            ls[i, i] = 1.0 / S[i];
        //        }
        //    }

        //    //(V x L*) x Ut x Y
        //    var vl = new double[V.GetLength(0), ls.GetLength(1)];
        //    for (var i = 0; i < V.GetLength(0); i++)
        //    {
        //        for (var j = 0; j < ls.GetLength(1); j++)
        //        {
        //            for (var k = 0; k < ls.GetLength(0); k++)
        //            {
        //                vl[i, j] += V[i, k] * ls[k, j];
        //            }
        //        }
        //    }

        //    //(V x L* x Ut) x Y
        //    var vlu = new double[vl.GetLength(0), U.GetLength(0)];
        //    for (var i = 0; i < vl.GetLength(0); i++)
        //    {
        //        for (var j = 0; j < U.GetLength(0); j++)
        //        {
        //            for (var k = 0; k < U.GetLength(0); k++)
        //            {
        //                vlu[i, j] += vl[i, k] * U[j, k];
        //            }
        //        }
        //    }

        //    //(V x L* x Ut x Y)
        //    var x = new double[vlu.GetLength(0), y.GetLength(1)];
        //    for (var i = 0; i < vlu.GetLength(0); i++)
        //    {
        //        for (var j = 0; j < y.GetLength(1); j++)
        //        {
        //            for (var k = 0; k < vlu.GetLength(1); k++)
        //            {
        //                x[i, j] += vlu[i, k] * y[k, j];
        //            }
        //        }
        //    }

        //    return x;
        //}

        //public double[] Solve(double[] value)
        //{
        //    // Additionally an important property is that if there does not exists a solution
        //    // when the matrix A is singular but replacing 1/Li with 0 will provide a solution
        //    // that minimizes the residue |AX -Y|. SVD finds the least squares best compromise
        //    // solution of the linear equation system. Interestingly SVD can be also used in an
        //    // over-determined system where the number of equations exceeds that of the parameters.

        //    // L is a diagonal matrix with non-negative matrix elements having the same
        //    // dimension as A, Wi ? 0. The diagonal elements of L are the singular values of matrix A.

        //    //singularity threshold
        //    var t = Threshold;

        //    var y = value;

        //    // Create L*, which is a diagonal matrix with elements
        //    //    L*i = 1/Li  if Li = e, else 0, 
        //    // where e is the so-called singularity threshold.

        //    // In other words, if Li is zero or close to zero (smaller than e),
        //    // one must replace 1/Li with 0. The value of e depends on the precision
        //    // of the hardware. This method can be used to solve linear equations
        //    // systems even if the matrices are singular or close to singular.

        //    var ls = new double[S.Length, S.Length];
        //    for (var i = 0; i < S.Length; i++)
        //    {
        //        if (Math.Abs(S[i]) < t)
        //        {
        //            ls[i, i] = 0;
        //        }
        //        else
        //        {
        //            ls[i, i] = 1.0 / S[i];
        //        }
        //    }

        //    //(V x L*) x Ut x Y
        //    var vl = new double[V.GetLength(0), ls.GetLength(1)];
        //    for (var i = 0; i < V.GetLength(0); i++)
        //    {
        //        for (var j = 0; j < ls.GetLength(1); j++)
        //        {
        //            for (var k = 0; k < V.GetLength(1); k++)
        //            {
        //                vl[i, j] += V[i, k] * ls[k, j];
        //            }
        //        }
        //    }

        //    //(V x L* x Ut) x Y
        //    var vlu = new double[vl.GetLength(0), U.GetLength(0)];
        //    for (var i = 0; i < vl.GetLength(0); i++)
        //    {
        //        for (var j = 0; j < U.GetLength(0); j++)
        //        {
        //            for (var k = 0; k < vl.GetLength(1); k++)
        //            {
        //                vlu[i, j] += vl[i, k] * U[j, k];
        //            }
        //        }
        //    }

        //    //(V x L* x Ut x Y)
        //    var x = new double[y.Length];

        //    for (var i = 0; i < vlu.GetLength(0); i++)
        //    {
        //        for (var j = 0; j < y.Length; j++)
        //        {
        //            x[i] += vlu[i, j] * y[j];
        //        }
        //    }

        //    return x;
        //}

        #endregion

        public double Hypotenuse(double a, double b)
        {
            var r = 0.0;

            if (Math.Abs(a) > Math.Abs(b))
            {
                r = b / a;
                r = Math.Abs(a) * Math.Sqrt(1 + r * r);
            }
            else if (b != 0)
            {
                r = a / b;
                r = Math.Abs(b) * Math.Sqrt(1 + r * r);
            }

            return r;
        }

        #endregion

        #region Methods

        private void Run(double[,] value, bool computeLeftSingularVectors, bool computeRightSingularVectors)
        {
            #region Variables

            //var aNormalized = Normalize((double[,])value.Clone());
            //A = (double[,])aNormalized.Clone();

            A = value;

            _mRow = value.GetLength(0); // rows
            _nCol = value.GetLength(1); // cols
            var nu = Math.Min(_mRow, _nCol);
            S = new double[Math.Min(_mRow + 1, _nCol)];
            U = new double[_mRow, nu];
            V = new double[_nCol, _nCol];
            var e = new double[_nCol];
            var work = new double[_mRow];
            var wantU = computeLeftSingularVectors;
            var wantV = computeRightSingularVectors;

            #endregion

            #region Reduce A to bidiagonal form, storing the diagonal elements in s and the super-diagonal elements in e.

            int nct;
            int nrt;
            ReduceA(A, e, work, wantU, wantV, out nct, out nrt);

            #endregion

            #region Set up the final bidiagonal matrix or order p.

            // Set up the final bi-diagonal matrix or order p.
            var p = Math.Min(_nCol, _mRow + 1);

            if (nct < _nCol)
                S[nct] = A[nct, nct];

            if (_mRow < p)
                S[p - 1] = 0.0;

            if (nrt + 1 < p)
                e[nrt] = A[nrt, p - 1];

            e[p - 1] = 0.0;

            #endregion

            #region If required, generate U.

            if (wantU)
            {
                GenerateU(nu, nct);
            }

            #endregion

            #region If required, generate V.

            if (wantV)
            {
                GenerateV(e, nrt);
            }

            #endregion

            #region Main iteration loop for the singular values.

            GenerateS(e, wantU, wantV, p);

            #endregion

        }

        private void ReduceA(double[,] a, IList<double> e, IList<double> work, bool wantu, bool wantv, out int nct, out int nrt)
        {
            // Reduce A to bi-diagonal form, storing the diagonal elements in s and the super-diagonal elements in e.
            nct = Math.Min(_mRow - 1, _nCol);
            nrt = Math.Max(0, Math.Min(_nCol - 2, _mRow));
            for (var k = 0; k < Math.Max(nct, nrt); k++)
            {
                if (k < nct)
                {
                    // Compute the transformation for the k-th column and place the k-th diagonal in s[k].
                    // Compute 2-norm of k-th column without under/overflow.
                    S[k] = 0;
                    for (var i = k; i < _mRow; i++)
                    {
                        S[k] = Hypotenuse(S[k], a[i, k]);
                    }

                    if (S[k] != 0.0)
                    {
                        if (a[k, k] < 0.0)
                        {
                            S[k] = -S[k];
                        }

                        for (var i = k; i < _mRow; i++)
                        {
                            a[i, k] /= S[k];
                        }

                        a[k, k] += 1.0;
                    }

                    S[k] = -S[k];
                }

                for (var j = k + 1; j < _nCol; j++)
                {
                    if ((k < nct) & (S[k] != 0.0))
                    {
                        // Apply the transformation.
                        double t = 0;
                        for (var i = k; i < _mRow; i++)
                            t += a[i, k] * a[i, j];
                        t = -t / a[k, k];
                        for (var i = k; i < _mRow; i++)
                            a[i, j] += t * a[i, k];
                    }

                    // Place the k-th row of A into e for the subsequent calculation of the row transformation.
                    e[j] = a[k, j];
                }

                if (wantu & (k < nct))
                {
                    // Place the transformation in U for subsequent back
                    // multiplication.
                    for (var i = k; i < _mRow; i++)
                        U[i, k] = a[i, k];
                }

                if (k < nrt)
                {
                    // Compute the k-th row transformation and place the k-th super-diagonal in e[k].
                    // Compute 2-norm without under/overflow.
                    e[k] = 0;
                    for (var i = k + 1; i < _nCol; i++)
                    {
                        e[k] = Hypotenuse(e[k], e[i]);
                    }

                    if (e[k] != 0.0)
                    {
                        if (e[k + 1] < 0.0)
                            e[k] = -e[k];

                        for (var i = k + 1; i < _nCol; i++)
                            e[i] /= e[k];

                        e[k + 1] += 1.0;
                    }

                    e[k] = -e[k];
                    if ((k + 1 < _mRow) & (e[k] != 0.0))
                    {
                        // Apply the transformation.
                        for (var i = k + 1; i < _mRow; i++)
                            work[i] = 0.0;

                        for (var j = k + 1; j < _nCol; j++)
                            for (var i = k + 1; i < _mRow; i++)
                                work[i] += e[j] * a[i, j];

                        for (var j = k + 1; j < _nCol; j++)
                        {
                            var t = -e[j] / e[k + 1];
                            for (var i = k + 1; i < _mRow; i++)
                                a[i, j] += t * work[i];
                        }
                    }

                    if (wantv)
                    {
                        // Place the transformation in V for subsequent back multiplication.
                        for (var i = k + 1; i < _nCol; i++)
                            V[i, k] = e[i];
                    }
                }
            }
        }

        private void GenerateU(int nu, int nct)
        {
            for (var j = nct; j < nu; j++)
            {
                for (var i = 0; i < _mRow; i++)
                    U[i, j] = 0.0;

                U[j, j] = 1.0;
            }

            for (var k = nct - 1; k >= 0; k--)
            {
                if (S[k] != 0.0)
                {
                    for (var j = k + 1; j < nu; j++)
                    {
                        double t = 0;

                        for (var i = k; i < _mRow; i++)
                            t += U[i, k] * U[i, j];

                        t = -t / U[k, k];

                        for (var i = k; i < _mRow; i++)
                            U[i, j] += t * U[i, k];
                    }

                    for (var i = k; i < _mRow; i++)
                        U[i, k] = -U[i, k];

                    U[k, k] = 1.0 + U[k, k];

                    for (var i = 0; i < k - 1; i++)
                        U[i, k] = 0.0;
                }
                else
                {
                    for (var i = 0; i < _mRow; i++)
                        U[i, k] = 0.0;

                    U[k, k] = 1.0;
                }
            }
        }

        private void GenerateV(IList<double> e, int nrt)
        {
            for (var k = _nCol - 1; k >= 0; k--)
            {
                if ((k < nrt) & (e[k] != 0.0))
                {
                    //TODO: Check if this is a bug.
                    //   for (var j = k + 1; j < nu; j++)
                    // The correction would be:
                    for (var j = k + 1; j < _nCol; j++)
                    {
                        double t = 0;

                        for (var i = k + 1; i < _nCol; i++)
                            t += V[i, k] * V[i, j];

                        t = -t / V[k + 1, k];

                        for (var i = k + 1; i < _nCol; i++)
                            V[i, j] += t * V[i, k];
                    }
                }

                for (var i = 0; i < _nCol; i++)
                    V[i, k] = 0.0;

                V[k, k] = 1.0;
            }
        }

        private void GenerateS(IList<double> e, bool wantu, bool wantv, int p)
        {
            // Main iteration loop for the singular values.
            var pp = p - 1;
            var iter = 0;
            var eps = Math.Pow(2.0, -52.0);

            while (p > 0)
            {
                int k, kase;

                // Here is where a test for too many iterations would go.
                // This section of the program inspects for
                // negligible elements in the s and e arrays.  On
                // completion the variables kase and k are set as follows.
                // kase = 1     if s(p) and e[k-1] are negligible and k<p
                // kase = 2     if s(k) is negligible and k<p
                // kase = 3     if e[k-1] is negligible, k<p, and s(k), ..., s(p) are not negligible (qr step).
                // kase = 4     if e(p-1) is negligible (convergence).
                for (k = p - 2; k >= -1; k--)
                {
                    if (k == -1)
                        break;

                    if (Math.Abs(e[k]) <= eps * (Math.Abs(S[k]) + Math.Abs(S[k + 1])))
                    {
                        e[k] = 0.0;
                        break;
                    }
                }

                if (k == p - 2)
                {
                    kase = 4;
                }
                else
                {
                    int ks;
                    for (ks = p - 1; ks >= k; ks--)
                    {
                        if (ks == k)
                            break;

                        double t = (ks != p ? Math.Abs(e[ks]) : 0.0) + (ks != k + 1 ? Math.Abs(e[ks - 1]) : 0.0);
                        if (Math.Abs(S[ks]) <= eps * t)
                        {
                            S[ks] = 0.0;
                            break;
                        }
                    }

                    if (ks == k)
                        kase = 3;
                    else if (ks == p - 1)
                        kase = 1;
                    else
                    {
                        kase = 2;
                        k = ks;
                    }
                }

                k++;

                // Perform the task indicated by kase.
                switch (kase)
                {
                    // Deflate negligible s(p).
                    case 1:
                        {
                            var f = e[p - 2];
                            e[p - 2] = 0.0;
                            for (var j = p - 2; j >= k; j--)
                            {
                                var t = Hypotenuse(S[j], f);
                                var cs = S[j] / t;
                                var sn = f / t;

                                S[j] = t;

                                if (j != k)
                                {
                                    f = -sn * e[j - 1];
                                    e[j - 1] = cs * e[j - 1];
                                }

                                if (wantv)
                                {
                                    for (var i = 0; i < _nCol; i++)
                                    {
                                        t = cs * V[i, j] + sn * V[i, p - 1];
                                        V[i, p - 1] = -sn * V[i, j] + cs * V[i, p - 1];
                                        V[i, j] = t;
                                    }
                                }
                            }
                        }
                        break;

                    // Split at negligible s(k).
                    case 2:
                        {
                            var f = e[k - 1];
                            e[k - 1] = 0.0;
                            for (var j = k; j < p; j++)
                            {
                                var t = Hypotenuse(S[j], f);
                                var cs = S[j] / t;
                                var sn = f / t;
                                S[j] = t;
                                f = -sn * e[j];
                                e[j] = cs * e[j];

                                if (wantu)
                                {
                                    for (var i = 0; i < _mRow; i++)
                                    {
                                        t = cs * U[i, j] + sn * U[i, k - 1];
                                        U[i, k - 1] = -sn * U[i, j] + cs * U[i, k - 1];
                                        U[i, j] = t;
                                    }
                                }
                            }
                        }
                        break;

                    // Perform one qr step.
                    case 3:
                        {
                            // Calculate the shift.
                            var scale = Math.Max(Math.Max(Math.Max(Math.Max(Math.Abs(S[p - 1]), Math.Abs(S[p - 2])), Math.Abs(e[p - 2])), Math.Abs(S[k])), Math.Abs(e[k]));
                            var sp = S[p - 1] / scale;
                            var spm1 = S[p - 2] / scale;
                            var epm1 = e[p - 2] / scale;
                            var sk = S[k] / scale;
                            var ek = e[k] / scale;
                            var b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                            var c = (sp * epm1) * (sp * epm1);
                            var shift = 0.0;

                            if ((b != 0.0) | (c != 0.0))
                            {
                                shift = Math.Sqrt(b * b + c);

                                if (b < 0.0)
                                    shift = -shift;

                                shift = c / (b + shift);
                            }

                            var f = (sk + sp) * (sk - sp) + shift;
                            var g = sk * ek;

                            // Chase zeros.
                            for (var j = k; j < p - 1; j++)
                            {
                                var t = Hypotenuse(f, g);
                                var cs = f / t;
                                var sn = g / t;

                                if (j != k)
                                    e[j - 1] = t;

                                f = cs * S[j] + sn * e[j];
                                e[j] = cs * e[j] - sn * S[j];
                                g = sn * S[j + 1];
                                S[j + 1] = cs * S[j + 1];

                                if (wantv)
                                {
                                    for (var i = 0; i < _nCol; i++)
                                    {
                                        t = cs * V[i, j] + sn * V[i, j + 1];
                                        V[i, j + 1] = -sn * V[i, j] + cs * V[i, j + 1];
                                        V[i, j] = t;
                                    }
                                }

                                t = Hypotenuse(f, g);
                                cs = f / t;
                                sn = g / t;
                                S[j] = t;
                                f = cs * e[j] + sn * S[j + 1];
                                S[j + 1] = -sn * e[j] + cs * S[j + 1];
                                g = sn * e[j + 1];
                                e[j + 1] = cs * e[j + 1];

                                if (wantu && (j < _mRow - 1))
                                {
                                    for (var i = 0; i < _mRow; i++)
                                    {
                                        t = cs * U[i, j] + sn * U[i, j + 1];
                                        U[i, j + 1] = -sn * U[i, j] + cs * U[i, j + 1];
                                        U[i, j] = t;
                                    }
                                }
                            }

                            e[p - 2] = f;
                            iter = iter + 1;
                        }
                        break;

                    // Convergence.
                    case 4:
                        {
                            // Make the singular values positive.
                            if (S[k] <= 0.0)
                            {
                                S[k] = (S[k] < 0.0 ? -S[k] : 0.0);
                                if (wantv)
                                    for (var i = 0; i <= pp; i++)
                                        V[i, k] = -V[i, k];
                            }

                            // Order the singular values.
                            while (k < pp)
                            {
                                if (S[k] >= S[k + 1])
                                    break;

                                var t = S[k];
                                S[k] = S[k + 1];
                                S[k + 1] = t;

                                if (wantv && (k < _nCol - 1))
                                {
                                    for (var i = 0; i < _nCol; i++)
                                    {
                                        t = V[i, k + 1];
                                        V[i, k + 1] = V[i, k];
                                        V[i, k] = t;
                                    }
                                }

                                if (wantu && (k < _mRow - 1))
                                {
                                    for (var i = 0; i < _mRow; i++)
                                    {
                                        t = U[i, k + 1];
                                        U[i, k + 1] = U[i, k];
                                        U[i, k] = t;
                                    }
                                }

                                k++;
                            }

                            iter = 0;
                            p--;
                        }
                        break;
                }
            }
        }

        #region Normalize / Denormalize

        //private double[,] Normalize(double[,] matrix)
        //{
        //    for (var i = 0; i < matrix.GetLength(0); i++)
        //    {
        //        for (var j = 0; j < matrix.GetLength(1); j++)
        //        {
        //            matrix[i, j] /= 255.0;
        //        }
        //    }

        //    return matrix;
        //}

        //public double[,] Denormalize(double[,] matrix)
        //{
        //    for (var i = 0; i < matrix.GetLength(0); i++)
        //    {
        //        for (var j = 0; j < matrix.GetLength(1); j++)
        //        {
        //            matrix[i, j] *= 255.0;
        //        }
        //    }

        //    return matrix;
        //}

        #endregion

        #endregion
    }
}
