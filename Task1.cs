public class Task1(int n = 3, double a = 1.1, double b = 2.5, int N = 2, int L = 2)
{
    protected static Func<double, double> f = (x) => 0.5 * Math.Cos(2 * x) * Math.Exp(2 * x / 5.0) + 2.4 * Math.Sin(1.5 * x) * Math.Exp(-6 * x) + 6 * x;
    static Func<double, double> p = (x) => 1.0 / Math.Pow(x - 1.1, 2.0 / 5.0);
    protected static Func<double, double> f_full = (x) => f(x) * p(x);

    Func<double, double> IntegrandGaussError(double[] x_vals, Func<double, double> p)
    {
        return x =>
        {
            double product_sq = 1.0;
            foreach (double xi in x_vals)
            {
                product_sq *= Math.Pow(x - xi, 2);
            }
            return product_sq * p(x);
        };
    }

    static Vector<double> ComputeMoments(int n, double a, double b)
    {
        var x = MathS.Var("x");
        var integrals = Vector<double>.Build.Dense(n);
        var mu = Vector<double>.Build.Dense(n);

        for (int i = 0; i < n; i++)
        {
            var expr = x.Pow(MathS.FromString($"{i}-2/5")).Simplify().Integrate(x);
            var I = (expr.Substitute(x, b - 1.1) - expr.Substitute(x, a - 1.1)).EvalNumerical();
            integrals[i] = (double)I;
            for (int j = 0; j <= i; j++)
            {
                mu[i] += MathNet.Numerics.SpecialFunctions.Binomial(i, j) * Math.Pow(1.1, i - j) * integrals[j];
            }
        }
        return mu;
    }

    double ComputeTol(int n, double a, double b, MethodNames M)
    {
        if (M == MethodNames.NewtonCotes)
        {
            var derivative = new NumericalDerivative(n + 1, (n + 1) / 2);
            double max = 0;
            for (int i = 0; i < 100; i++)
            {
                double t = a + i * (b - a) / 999.0;
                max = Math.Abs(derivative.EvaluateDerivative(f, t, n)) > max ? Math.Abs(derivative.EvaluateDerivative(f, t, n)) : max;
            }

            double[] x_vals = MathNet.Numerics.Generate.LinearSpaced(n, 0, b - a);
            var terms = x_vals.Skip(1).Select(xi => $"(x - {xi.ToString(System.Globalization.CultureInfo.InvariantCulture)})").ToList();
            string w = string.Join(" * ", terms);
            var x = MathS.Var("x");
            var uzl = MathS.FromString(w).Simplify();
            var p = (x.Pow(MathS.FromString("3/5")) * uzl).Integrate(x);
            double factor = Math.Pow(-1, n - 1);
            double tot = 0;
            for (int i = 1; i < n; i++) { tot += factor * (double)(p.Substitute(x, x_vals[i]) - p.Substitute(x, x_vals[i - 1])).EvalNumerical(); factor *= -1; }
            double estimate = max / MathNet.Numerics.SpecialFunctions.Factorial(n) * tot;
            return estimate;
        }
        else
        {
            var derivative = new NumericalDerivative(2 * n + 1, (2 * n + 1) / 2);
            double max = 0;
            for (int i = 0; i < 100; i++)
            {
                double t = a + i * (b - a) / 999.0;
                max = Math.Abs(derivative.EvaluateDerivative(f, t, 2 * n)) > max ? Math.Abs(derivative.EvaluateDerivative(f, t, 2 * n)) : max;
            }

            var mu = ComputeMoments(2 * n, a, b);
            var X = Matrix<double>.Build.Dense(n, n);
            for (int i = 0; i < n; i++) { X.SetRow(i, mu.SubVector(i, n)); }
            var coef = X.Solve(-mu.SubVector(n, n));
            var poly = new double[n + 1];
            for (int i = 0; i < n; i++) { poly[i] = coef[i]; }
            poly[n] = 1;
            var polynomial = new Polynomial(poly);
            var x_vals = polynomial.Roots().Select(x => x.Real).OrderBy(x => x).ToArray();
            var integrand_function = IntegrandGaussError(x_vals, p);
            double tot = NewtonCotesTrapeziumRule.IntegrateComposite(integrand_function, a + 1e-8, b, 100);
            double estimate = max / MathNet.Numerics.SpecialFunctions.Factorial(2 * n) * tot;
            return estimate;
        }

    }

    protected double SmallNewtonCotes(int n, double a, double b)
    {
        var mu = ComputeMoments(n, a, b);
        var x_vals = Vector<double>.Build.Dense(n, i => a + i * (b - a) / (n - 1));
        var X = Matrix<double>.Build.Dense(n, n);
        for (int i = 0; i < n; i++) { X.SetRow(i, x_vals.Map(x => Math.Pow(x, i))); }
        var A = X.Solve(mu);
        double res1 = A.DotProduct(x_vals.Map(f));
        return res1;
    }

    protected double SmallGaussian(int n, double a, double b)
    {
        var mu = ComputeMoments(2 * n, a, b);
        var X = Matrix<double>.Build.Dense(n, n);
        for (int i = 0; i < n; i++) { X.SetRow(i, mu.SubVector(i, n)); }
        var coef = X.Solve(-mu.SubVector(n, n));
        var poly = new double[n + 1];
        for (int i = 0; i < n; i++) { poly[i] = coef[i]; }
        poly[n] = 1;
        var polynomial = new Polynomial(poly);
        var x_vals = Vector<double>.Build.DenseOfArray(polynomial.Roots().Select(x => x.Real).OrderBy(x => x).ToArray());
        var Y = Matrix<double>.Build.Dense(n, n);
        for (int i = 0; i < n; i++) { Y.SetRow(i, x_vals.Map(x => Math.Pow(x, i))); }
        var A = Y.Solve(mu.SubVector(0, n));
        double res1 = A.DotProduct(x_vals.Map(f));
        return res1;
    }

    protected double OneIteration(int n, double a, double b, Method F)
    {
        double h = (b - a) / n;
        double val = 0;
        for (int i = 0; i < n; i++)
        {
            val += F(3, a + i * h, a + (i + 1) * h);
        }
        return val;
    }

    double CompositeQuadrature(int N, int L, double a, double b, Method F, out double m, double tol = 1e-6)
    {
        int j = 3;
        var S = new double[3];
        var H = new double[3];
        for (int i = 0; i < 3; i++) { S[i] = OneIteration(N, a, b, F); H[i] = (b - a) / N; N *= L; }
        m = -Math.Log((S[2] - S[1]) / (S[1] - S[0])) / Math.Log(L);
        var s = Vector<double>.Build.DenseOfArray(new double[] { S[1] - S[0], S[2] - S[1] });
        var X = Matrix<double>.Build.DenseOfArray(new double[,]
        {
    {Math.Pow(H[0], m)*(1-1.0/Math.Pow(L, m)), Math.Pow(H[0], m+1)*(1-1.0/Math.Pow(L, m+1)) },
    {Math.Pow(H[1], m)*(1-1.0/Math.Pow(L, m)), Math.Pow(H[1], m+1)*(1-1.0/Math.Pow(L, m+1)) }
        });
        var C = X.Solve(s);
        double eps = C.DotProduct(Vector<double>.Build.DenseOfArray(new double[] { Math.Pow(H[2], m), Math.Pow(H[2], m + 1) }));
        Console.WriteLine(string.Format("Iteration: {0}, m: {1:F2}, tol: {2}", j, m, eps));
        while (Math.Abs(eps) > tol)
        {
            S[0] = S[1]; S[1] = S[2];
            H[0] = H[1]; H[1] = H[2]; H[2] = (b - a) / N;
            S[2] = OneIteration(N, a, b, F); N *= L; j++;
            m = -Math.Log((S[2] - S[1]) / (S[1] - S[0])) / Math.Log(L);
            s = Vector<double>.Build.DenseOfArray(new double[] { S[1] - S[0], S[2] - S[1] });
            X = Matrix<double>.Build.DenseOfArray(new double[,]
            {
        {Math.Pow(H[0], m)*(1-1.0/Math.Pow(L, m)), Math.Pow(H[0], m+1)*(1-1.0/Math.Pow(L, m+1)) },
        {Math.Pow(H[1], m)*(1-1.0/Math.Pow(L, m)), Math.Pow(H[1], m+1)*(1-1.0/Math.Pow(L, m+1)) }
            });
            C = X.Solve(s);
            eps = C.DotProduct(Vector<double>.Build.DenseOfArray(new double[] { Math.Pow(H[2], m), Math.Pow(H[2], m + 1) }));
            Console.WriteLine(string.Format("Iteration: {0}, m: {1:F2}, tol: {2}", j, m, eps));
        }
        return S[2];
    }

    double OptimalQuadrature(int N, int L, double a, double b, double m, Method F, double tol = 1e-6)
    {
        var S = new double[2];
        double h = (b - a) / N;
        for (int i = 0; i < 2; i++) { S[i] = OneIteration(N, a, b, F); N *= L; }
        double hopt = h * Math.Pow(tol * (1 - Math.Pow(L, -m)) / (S[1] - S[0]), 1.0 / m);
        Console.WriteLine(string.Format("Running composite quadrature with optimal step: {0:F3}", hopt));
        return CompositeQuadrature((int)Math.Ceiling((b - a) / (0.95 * hopt)), L, a, b, F, out _);
    }

    public void print()
    {
        var res = DoubleExponentialTransformation.Integrate(f_full, a + 1e-12, b, 1000);
        Console.WriteLine(string.Format("True value: {0:F6}", res));

        Console.WriteLine(string.Format("Newton-Cotes: {0:F6} +- {1:F6}\n", SmallNewtonCotes(n, a, b), ComputeTol(n, a, b, MethodNames.NewtonCotes)));

        Console.WriteLine("Running composite quadrature");
        Console.WriteLine(string.Format("{0:F6}\n", CompositeQuadrature(N, L, a, b, SmallNewtonCotes, out double m)));
        Console.WriteLine(string.Format("{0:F6}\n", OptimalQuadrature(N, L, a, b, m, SmallNewtonCotes)));

        Console.WriteLine(string.Format("Gauss: {0:F6} +- {1:F6}\n", SmallGaussian(n, a, b), ComputeTol(n, a, b, MethodNames.Gaussian)));

        Console.WriteLine("Running composite quadrature");
        Console.WriteLine(string.Format("{0:F6}\n", CompositeQuadrature(N, L, a, b, SmallGaussian, out m)));
        Console.WriteLine(string.Format("{0:F6}\n", OptimalQuadrature(N, L, a, b, m, SmallGaussian)));
    }

    protected delegate double Method(int n, double a, double b);
    enum MethodNames { NewtonCotes, Gaussian }
}

