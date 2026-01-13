public class Task2
{
    static double CubeRoot(double value)
    {
        return value < 0 ? -Math.Pow(-value, 1.0 / 3.0) : Math.Pow(value, 1.0 / 3.0);
    }
    static Vector<double> ODE(double x, Vector<double> y)
    {
        var res = Vector<double>.Build.Dense(4);
        res[0] = 2 * x * CubeRoot(y[1]) * y[3];
        res[1] = 6 * x * Math.Exp(y[2] - 1) * y[3];
        res[2] = 6 * x * y[3];
        res[3] = -2 * x * Math.Log(y[0]);
        return res;
    }
    static Vector<double> TrueODE(double x)
    {
        var res = Vector<double>.Build.Dense(4);
        res[0] = Math.Exp(Math.Sin(Math.Pow(x, 2)));
        res[1] = Math.Exp(3 * Math.Sin(Math.Pow(x, 2)));
        res[2] = 3 * Math.Sin(Math.Pow(x, 2)) + 1;
        res[3] = Math.Cos(Math.Pow(x, 2));
        return res;
    }
    static Vector<double> SecOrderIt(double x0, Vector<double> y0, double h, out Vector<double> K1, Vector<double>? K = null)
    {
        K1 = K is null ? ODE(x0, y0) : K;
        var K2 = ODE(x0 + 0.6 * h, y0 + 0.6 * h * K1);
        return y0 + h * (1.0 / 6.0 * K1 + 5.0 / 6.0 * K2);
    }

    static Vector<double> FourthOrderIt(double x0, Vector<double> y0, double h, out Vector<double> K1, Vector<double>? K = null)
    {
        K1 = K is null ? ODE(x0, y0) : K;
        var K2 = ODE(x0 + 0.5 * h, y0 + 0.5 * h * K1);
        var K3 = ODE(x0 + 0.5 * h, y0 + 0.5 * h * K2);
        var K4 = ODE(x0 + h, y0 + h * K3);
        return y0 + h * (1.0 / 6.0 * K1 + 1.0 / 3.0 * K2 + 1.0 / 3.0 * K3 + 1.0 / 6.0 * K4);
    }
    Result1 BuildTrue(double x0, double x1, double h)
    {
        int n = (int)Math.Ceiling((x1 - x0) / h);
        h = (x1 - x0) / n;
        var res = new Result1(x0, TrueODE(x0));
        for (int i = 0; i < n; i++)
        {
            var z = TrueODE(x0 + (i + 1) * h);
            res.add(x0 + (i + 1) * h, z);
        }
        return res;
    }
    Result1 BuildTrueL(List<double> x)
    {
        var res = new Result1();
        foreach (double i in x)
        {
            var z = TrueODE(i);
            res.add(i, z);
        }
        return res;
    }
    Result1 FirstAlgorithm(Method method, double x0, double x1, Vector<double> y0, double h)
    {
        int n = (int)Math.Ceiling((x1 - x0) / h);
        h = (x1 - x0) / n;
        var res = new Result1(x0, y0);
        for (int i = 0; i < n; i++)
        {
            var z = method(res.x_vals[^1], Vector<double>.Build.DenseOfArray(new double[] { res[0][^1], res[1][^1], res[2][^1], res[3][^1] }), h, out _);
            res.add(x0 + (i + 1) * h, z);
        }
        return res;
    }

    double GetInitialStep(double x0, double x1, Vector<double> y0, MethodName M, out int i, double tol = 1e-12)
    {
        Vector<double> K = ODE(x0, y0);
        int p = ((int)M + 1) * 2;
        var delta = Math.Pow(1.0 / Math.Max(Math.Abs(x0), Math.Abs(x1)), p + 1) + Math.Pow(K.L2Norm(), p + 1);
        double h = Math.Pow(tol / delta, 1.0 / (p + 1));
        if (K.L1Norm() < 1)
        {
            i = 3;
            y0 += h * ODE(x0, y0);
            var K1 = ODE(x0 + h, y0);
            delta = Math.Pow(1.0 / Math.Max(Math.Abs(x0 + h), Math.Abs(x1)), p + 1) + Math.Pow(K1.L2Norm(), p + 1);
            double h1 = Math.Pow(tol / delta, 1.0 / (p + 1));
            return Math.Min(h, h1);
        }
        else
        {
            i = 1;
            return h;
        }
    }

    Result2 SecondAlgorithm(double x0, double x1, Vector<double> y0, MethodName M, Method m, double rtol = 1e-6, double atol = 1e-12)
    {
        double x = x0;
        Vector<double> y = y0;
        double count = 0;
        int s = ((int)M + 1) * 2;
        double hmax = 0.25/*GetInitialStep(x, x1, y, M, out int i, atol)*/;
        count += 3;
        double h = hmax;
        var res = new Result2(x0, y0);
        while (x < x1)
        {
            while (true)
            {
                double h1 = h;
                double h2 = h1 / 2.0;
                Vector<double> y1 = m(x, y, h1, out Vector<double> K1);
                Vector<double> y12 = m(x, y, h2, out _, K1);
                Vector<double> y2 = m(x + h2, y12, h2, out _);
                count += 3 * s - 1;
                var eps = (y2 - y1) / (1 - Math.Pow(2, -s));
                double tol = rtol * y2.InfinityNorm() + atol;
                if (eps.InfinityNorm() > tol * Math.Pow(2, s))
                {
                    res.addh(x, h1);
                    h /= 2;
                }
                else if (tol < eps.InfinityNorm() && eps.InfinityNorm() <= tol * Math.Pow(2, s))
                {
                    res.add(x + h2, y12);
                    res.add(x + h1, y2);
                    res.addh(x, h1);
                    res.addhtrue(x, h2);
                    x += h;
                    h /= 2;
                    y = y2;
                    break;
                }
                else if (tol * Math.Pow(2, -(s + 1)) <= eps.InfinityNorm() && eps.InfinityNorm() <= tol)
                {
                    res.add(x + h1, y1);
                    res.addhtrue(x, h1);
                    x += h;
                    y = y1;
                    break;
                }
                else
                {
                    res.add(x + h1, y1);
                    res.addhtrue(x, h1);
                    x += h;
                    h = Math.Min(2 * h, hmax);
                    y = y1;
                    break;
                }
            }
        }
        res.count = count;
        return res;
    }
    void FirstPlot()
    {
        var h_vals = Vector<double>.Build.Dense(7);
        var norm1_vals = Vector<double>.Build.Dense(7);
        var norm2_vals = Vector<double>.Build.Dense(7);
        var y0 = Vector<double>.Build.DenseOfArray(new double[] { 1, 1, 1, 1 });
        var ytrue = TrueODE(5);
        for (int i = 0; i < 7; i++)
        {
            h_vals[i] = Math.Pow(2, i - 5 - 6);
            var y1 = FirstAlgorithm(SecOrderIt, 0, 5, y0, h_vals[i]).ReturnLast();
            var y2 = FirstAlgorithm(FourthOrderIt, 0, 5, y0, h_vals[i]).ReturnLast();
            
            norm1_vals[i] = (ytrue - y1).L2Norm();
            norm2_vals[i] = (ytrue - y2).L2Norm();
        }
        var h_ = h_vals.Map(x => Math.Log2(x)).ToArray();
        var norm1_ = norm1_vals.Map(x => Math.Log2(x)).ToArray();
        var norm2_ = norm2_vals.Map(x => Math.Log2(x)).ToArray();

        var plt = new Plot();
        var scatter1 = plt.Add.Scatter(h_, norm1_);
        scatter1.LegendText = "Second order RK";
        scatter1.MarkerStyle.Shape = ScottPlot.MarkerShape.FilledCircle;
        scatter1.LineStyle = ScottPlot.LineStyle.None; // Отображаем только точки

        var scatter2 = plt.Add.Scatter(h_, norm2_);
        scatter2.LegendText = "\"The\" Runge-Kutta method";
        scatter2.MarkerStyle.Shape = ScottPlot.MarkerShape.FilledCircle;
        scatter2.LineStyle = ScottPlot.LineStyle.None;

        var line1 = plt.Add.Line(h_[0], norm1_[0], h_[6], 2 * (h_[6] - h_[0]) + norm1_[0]);
        line1.Color = Color.FromARGB(0xFF0000FF); // Синий
        line1.LineStyle.Pattern = LinePattern.Dashed;
        line1.LegendText = "Slope 2";

        var line2 = plt.Add.Line(h_[0], norm2_[0], h_[6], 4 * (h_[6] - h_[0]) + norm2_[0]);
        line2.Color = Color.FromARGB(0xFFFF0000); // Красный
        line2.LineStyle.Pattern = LinePattern.Dashed;
        line2.LegendText = "Slope 4";

        plt.Title("||y_true - y|| на последнем шаге ~ h");
        plt.Axes.Bottom.Label.Text = "Шаг интегрирования";
        plt.Axes.Left.Label.Text = "Евклидова норма погрешности";
        plt.ShowLegend();

        plt.Axes.AutoScale();

        plt.SavePng("../../../plots/first_plot.png", 850, 600);
    }

    void SecondPlot()
    {
        var y0 = Vector<double>.Build.DenseOfArray(new double[] { 1, 1, 1, 1 });
        var y1 = new Vector<double>[2];
        var y2 = new Vector<double>[2];
        for (int i = 1; i < 3; i++)
        {
            y1[i - 1] = FirstAlgorithm(SecOrderIt, 0, 5, y0, Math.Pow(0.5, i)).ReturnLast();
            y2[i - 1] = FirstAlgorithm(FourthOrderIt, 0, 5, y0, Math.Pow(0.25, i)).ReturnLast();
        }
        var r1 = ((y1[1] - y1[0]) / (1 - Math.Pow(2, -2))).L2Norm();
        var r2 = ((y2[1] - y2[0]) / (1 - Math.Pow(2, -4))).L2Norm();
        var hopt1 = 0.5 * Math.Pow(1e-5 / r1, 1.0 / 2.0);
        var hopt2 = 0.25 * Math.Pow(1e-5 / r2, 1.0 / 4.0);
        var R1 = FirstAlgorithm(SecOrderIt, 0, 5, y0, hopt1);
        var R2 = FirstAlgorithm(FourthOrderIt, 0, 5, y0, hopt2);
        var ytrue1 = BuildTrue(0, 5, hopt1);
        var ytrue2 = BuildTrue(0, 5, hopt2);
        var y1_ = ytrue1 - R1;
        var y2_ = ytrue2 - R2;
        var h1_ = ytrue1.x_vals.ToArray();
        var h2_ = ytrue2.x_vals.ToArray();
        var plt = new Plot();
        var scatter1 = plt.Add.Scatter(h1_, y1_);
        scatter1.LegendText = $"Second order RK, h_opt = {hopt1.ToString("F5")}";
        scatter1.MarkerStyle.Shape = ScottPlot.MarkerShape.FilledCircle;
        scatter1.MarkerSize = 3;

        var scatter2 = plt.Add.Scatter(h2_, y2_);
        scatter2.LegendText = $"\"The\" Runge-Kutta method, h_opt = {hopt2.ToString("F5")}";
        scatter2.MarkerStyle.Shape = ScottPlot.MarkerShape.FilledCircle;
        scatter2.MarkerSize = 3;

        plt.Title("||y_true - y|| при h_opt ~ x");
        plt.Axes.Bottom.Label.Text = "Координата x";
        plt.Axes.Left.Label.Text = "Евклидова норма погрешности";
        plt.ShowLegend(Alignment.UpperLeft);

        plt.Axes.AutoScale();

        plt.SavePng("../../../plots/second_plot.png", 850, 600);
    }

    void ThirdPlot()
    {
        var y0 = Vector<double>.Build.DenseOfArray(new double[] { 1, 1, 1, 1 });
        var res1 = SecondAlgorithm(0, 5, y0, MethodName.RK2, SecOrderIt);
        var res2 = SecondAlgorithm(0, 5, y0, MethodName.RK4, FourthOrderIt);
        var t = BuildTrueL(res1.x_vals);

        var multiplot = new Multiplot();
        multiplot.AddPlots(4);
        var x1_ = res1.x_vals.ToArray();
        var x2_ = res2.x_vals.ToArray();

        for (int i = 0; i < 4; i++)
        {
            var plot = multiplot.GetPlot(i);

            var C = plot.Add.Scatter(x1_, t[i].ToArray());
            C.LegendText = "True values";
            C.MarkerSize = 5;
            C.LineWidth = 5;
            C.Color = Colors.Green;

            var A = plot.Add.Scatter(x1_, res1[i].ToArray());
            A.LegendText = "RK2";
            A.MarkerSize = 2;
            A.LineWidth = 0;
            A.Color = Colors.Orange;

            var B = plot.Add.Scatter(x2_, res2[i].ToArray());
            B.LegendText = "RK4";
            B.MarkerSize = 2;
            B.LineWidth = 0;
            B.Color = Colors.Blue;

            plot.Axes.Bottom.Label.Text = "x";
            plot.Axes.Left.Label.Text = $"y{i + 1}";
            if (i == 0) { plot.ShowLegend(Alignment.UpperLeft); plot.Title("Построенные решения"); }
            else plot.HideLegend();
        }
        multiplot.Layout = new ScottPlot.MultiplotLayouts.Grid(rows: 2, columns: 2);
        multiplot.SavePng("../../../plots/third_plot.png", 1200, 800);

        var plt = new Plot();
        plt.Title("h ~ x");

        var D = plt.Add.Scatter(res1.h_x.ToArray(), res1.h.ToArray());
        D.MarkerShape = MarkerShape.FilledTriangleUp;
        D.MarkerSize = 5;
        D.LineWidth = 0;
        D.LegendText = "thrown h for RK2";
        D.Color = Colors.Red;

        var E = plt.Add.Scatter(res1.htrue_x.ToArray(), res1.htrue.ToArray());
        E.LegendText = "accepted h for RK2";
        E.MarkerSize = 2;
        E.LineWidth = 2;
        E.Color = Colors.Black;

        var F = plt.Add.Scatter(res2.h_x.ToArray(), res2.h.ToArray());
        F.MarkerShape = MarkerShape.Cross;
        F.MarkerSize = 5;
        F.LineWidth = 0;
        F.LegendText = "thrown h for RK4";
        F.Color = Colors.BlueViolet;

        var G = plt.Add.Scatter(res2.htrue_x.ToArray(), res2.htrue.ToArray());
        G.LegendText = "accepted h for RK4";
        G.MarkerSize = 2;
        G.LineWidth = 2;
        G.Color = Colors.YellowGreen;

        plt.Axes.Bottom.Label.Text = "x";
        plt.Axes.Left.Label.Text = "h";
        plt.ShowLegend(Alignment.MiddleLeft);
        plt.SavePng("../../../plots/fourth_plot.png", 850, 600);

        var u = BuildTrueL(res2.x_vals);

        var p = new Plot();
        p.Title("||y_true - y|| ~ x");
        var H = p.Add.Scatter(t.x_vals.ToArray(), t - res1);
        H.LegendText = "RK2";

        var I = p.Add.Scatter(u.x_vals.ToArray(), u - res2);
        I.LegendText = "RK4";

        p.Axes.Bottom.Label.Text = "x";
        p.Axes.Left.Label.Text = "Error L2 norm";
        p.ShowLegend(Alignment.MiddleLeft);
        p.SavePng("../../../plots/fifth_plot.png", 850, 600);
    }

    void FourthPlot()
    {
        var y0 = Vector<double>.Build.DenseOfArray(new double[] { 1, 1, 1, 1 });
        var count1 = Vector<double>.Build.Dense(5);
        var count2 = Vector<double>.Build.Dense(5);
        var x_vals = Vector<double>.Build.Dense(5);
        for (int i = 0; i < 5; i++)
        {
            double rtol = Math.Pow(10, i - 8);
            x_vals[i] = rtol;
            count1[i] = SecondAlgorithm(0, 5, y0, MethodName.RK2, SecOrderIt, rtol).count;
            count2[i] = SecondAlgorithm(0, 5, y0, MethodName.RK4, FourthOrderIt, rtol).count;
        }
        var x_ = x_vals.Map(Math.Log2).ToArray();
        var c1_ = count1.Map(Math.Log2).ToArray();
        var c2_ = count2.Map(Math.Log2).ToArray();

        var plt = new Plot();
        plt.Title("Число обращений к правой части от rtol");

        var sc1 = plt.Add.Scatter(x_, c1_);
        sc1.LegendText = "RK2";

        var sc2 = plt.Add.Scatter(x_, c2_);
        sc1.LegendText = "RK4";

        plt.Axes.Bottom.Label.Text = "Log2(rtol)";
        plt.Axes.Left.Label.Text = "Log2(count)";
        plt.SavePng("../../../plots/sixth_plot.png", 850, 600);
    }

    public void Plot()
    {
        FirstPlot();
        SecondPlot();
        ThirdPlot();
        FourthPlot();
    }

    delegate Vector<double> Method(double x0, Vector<double> y0, double h, out Vector<double> K1, Vector<double>? K = null);
    enum MethodName { RK2, RK4 }

    class Result1
    {
        List<double>[] arrays;
        public List<double> x_vals;
        public Result1()
        {
            arrays = new List<double>[4];
            for (int i = 0; i < 4; i++)
            {
                arrays[i] = new List<double>();
            }
            x_vals = new List<double>();
        }
        public Result1(double x0, Vector<double> y0) : this()
        {
            x_vals.Add(x0);
            for (int i = 0; i < 4; i++)
            {
                arrays[i].Add(y0[i]);
            }
        }
        public List<double> this[int i]
        {
            get => arrays[i];
            set => arrays[i] = value;
        }
        public void add(double x, Vector<double> y)
        {
            x_vals.Add(x);
            for (int i = 0; i < 4; i++)
            {
                arrays[i].Add(y[i]);
            }
        }
        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine("x_vals:");
            sb.AppendLine(string.Join(", ", x_vals.ToArray()));

            sb.AppendLine("y_vals:");
            for (int i = 0; i < arrays.Length; i++)
            {
                sb.AppendLine(string.Join(", ", this[i]));
            }
            return sb.ToString();
        }
        public Vector<double> ReturnLast() => Vector<double>.Build.DenseOfArray(new double[] { this[0][^1], this[1][^1], this[2][^1], this[3][^1] });
        public static double[] operator -(Result1 r1, Result1 r2)
        {
            int n = r1.x_vals.Count;
            var res = new double[n];
            for (int i = 0; i < 4; i++)
            {
                r1[i] = r1[i].Zip(r2[i], (a, b) => a - b).ToList();
            }
            for (int i = 0; i < n; i++)
            {
                res[i] = Vector<double>.Build.DenseOfArray(new double[] { r1[0][i], r1[1][i], r1[2][i], r1[3][i] }).L2Norm();
            }
            return res;
        }
    };

    class Result2 : Result1
    {
        public List<double> h_x;
        public List<double> htrue_x;
        public List<double> h;
        public List<double> htrue;
        public double count;
        public Result2() : base()
        {
            h = new List<double>();
            h_x = new List<double>();
            htrue = new List<double>();
            htrue_x = new List<double>();
            count = 0;
        }
        public Result2(double x0, Vector<double> y0) : base(x0, y0)
        {
            h = new List<double>();
            h_x = new List<double>();
            htrue = new List<double>();
            htrue_x = new List<double>();
            count = 0;
        }
        public void addh(double x, double y)
        {
            h_x.Add(x);
            h.Add(y);
        }
        public void addhtrue(double x, double y)
        {
            htrue_x.Add(x);
            htrue.Add(y);
        }
    }
}