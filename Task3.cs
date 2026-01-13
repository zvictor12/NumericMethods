public class Task3 : Task1
{
    const double A = 1.1;
    const double B = 2.5;
    const double RTOL = 1e-10;
    const double EXACT_VALUE1 = 14.27314090243354910687;
    const double EXACT_VALUE2 = 18.60294785731848208626;
    public class ResultTableRow
    {
        public int N { get; set; }
        public double h { get; set; }
        public double S { get; set; }
        public double? E { get; set; }
        public double? R { get; set; }
        public double? p { get; set; }
        public double? S_prime { get; set; }
        public double? E_prime { get; set; }
        public double? p_prime { get; set; }
        public double? p_r { get; set; }
    }

    double MidpointRule(int n, double a, double b)
    {
        double h = (a+b) / 2.0;
        return f(h) * (b - a);
    }

    double TrapezoidalRule(int n, double a, double b)
    {
        return 0.5*(f(a)+f(b))*(b-a);
    }

    double SimpsonRule(int n, double a, double b)
    {
        double h = (a+b) / 2;
        return (b-a) / 6 * (f(a)+4*f(h)+f(b));
    }

    List<ResultTableRow> ComputeIntegral(Methods methodName, Method method, double a, double b, int N = 1)
    {
        double EXACT_VALUE = (int)methodName < 3 ? EXACT_VALUE1 : EXACT_VALUE2;
        Console.WriteLine($"\n--- Method: {methodName.ToString()} ---");
        Console.WriteLine($"THEORETICAL VALUE: {EXACT_VALUE}");

        List<ResultTableRow> table = new List<ResultTableRow>();
        List<double> S_values = new List<double>();
        List<double> h_values = new List<double>();
        List<double> R_values = new List<double>();

        int iteration = 1;
        bool converged = false;

        Console.WriteLine($"{"N",6} {"h",12} {"S",20} {"E",12} {"R",12} {"p",8} {"S'",20} {"E'",12} {"p'",8} {"pr",8}");
        Console.WriteLine(new string('-', 130));

        while (!converged)
        {
            double S = OneIteration(N, a, b, method);
            double E = Math.Abs(EXACT_VALUE - S) / Math.Abs(EXACT_VALUE);

            S_values.Add(S);
            h_values.Add((b - a) / N);

            var row = new ResultTableRow
            {
                N = N,
                h = (b - a) / N,
                S = S,
                E = E
            };

            // Вычисляем дополнительные метрики (начиная с 3-й итерации)
            if (S_values.Count >= 3)
            {

                // Порядок сходимости по Эйткену
                double S_n = S_values[^1];
                double S_n1 = S_values[^2];
                double S_n2 = S_values[^3];
                row.p = Math.Log(Math.Abs((S_n - S_n1) / (S_n1 - S_n2))) / Math.Log(0.5);

                // Оценка погрешности по правилу Рунге
                row.R = Math.Abs(S_n - S_n1) / (Math.Pow(2, row.p.Value) - 1) / Math.Abs(S_n);
                R_values.Add(row.R.Value);

                // Уточнение по Ричардсону
                row.S_prime = S_n + (S_n - S_n1) / (Math.Pow(2, row.p.Value) - 1);
                row.E_prime = Math.Abs(EXACT_VALUE - row.S_prime.Value) / Math.Abs(EXACT_VALUE);

                // Порядок для уточненного значения
                row.p_prime = Math.Log(Math.Abs((row.S_prime.Value - S_n1) / (S_n1 - S_n2))) / Math.Log(0.5);
                

                // Наклон на графике (по последним 2 точкам)
                if (R_values.Count >= 2)
                {
                    row.p_r = (Math.Log10(R_values[^1]) - Math.Log10(R_values[^2])) /
                                   (Math.Log10(h_values[^1]) - Math.Log10(h_values[^2]));
                }
            }

            table.Add(row);

            // Вывод строки
            Console.WriteLine($"{row.N,6} {row.h,12:E3} {row.S,20:F15} {row.E,12:E4} " +
                             $"{(row.R.HasValue ? row.R.Value.ToString("E4") : "---"),12} " +
                             $"{(row.p.HasValue ? row.p.Value.ToString("F3") : "---"),8} " +
                             $"{(row.S_prime.HasValue ? row.S_prime.Value.ToString("F15") : "---"),20} " +
                             $"{(row.E_prime.HasValue ? row.E_prime.Value.ToString("E4") : "---"),12} " +
                             $"{(row.p_prime.HasValue ? row.p_prime.Value.ToString("F3") : "---"),8} " +
                             $"{(row.p_r.HasValue ? row.p_r.Value.ToString("F3") : "---"),8}");

            // Проверка сходимости
            if (row.R.HasValue && row.S_prime.HasValue && row.R.Value < RTOL * row.S_prime.Value)
            {
                converged = true;
            }

            N *= 2;
            iteration++;
        }
        return table;
    }

    private void PlotTask1Graph(Methods methodName, List<ResultTableRow> table, double p_theory)
    {
        var h_values = new List<double>();
        var R_values = new List<double>();
        for(int i = table.Count - 1;  i > 1; i--)
        {
            h_values.Add(table[i].h);
            R_values.Add(table[i].R.Value);
        }
        var plot = new Plot();
        plot.Title($"error ~ h: {methodName.ToString()}");
        plot.XLabel("log10(h)");
        plot.YLabel("log10(R)");

        var log_h = h_values.Select(x => Math.Log10(x)).ToArray();
        var log_R = R_values.Select(x => Math.Log10(x)).ToArray();

        var scatter = plot.Add.Scatter(log_h, log_R);
        scatter.LegendText = "Фактическая погрешность";
        scatter.LineWidth = 0;
        scatter.MarkerSize = 5;
        scatter.Color = Colors.Blue;

        // Опорная прямая с теоретическим наклоном
        double x1 = log_h[0];
        double x2 = log_h[^1];
        double y1 = log_R[0];
        double y2 = log_R[0] + p_theory * (x2 - x1);

        var theoryLine = plot.Add.Line(x1, y1, x2, y2);
        theoryLine.LegendText = $"Theoretical p = {p_theory.ToString("F2")}";
        theoryLine.LineWidth = 2;
        theoryLine.LineStyle.Pattern = LinePattern.Dashed;
        theoryLine.Color = Colors.Orange;

        plot.ShowLegend();
        plot.SavePng($"../../../intplots/{methodName.ToString()}.png", 800, 600);
    }

    private void PlotTask2Graph(Methods methodName, List<ResultTableRow> table, double[] referenceSlopes)
    {
        var h_values = new List<double>();
        var errors = new List<double>();
        for (int i = table.Count - 1; i > 1; i--)
        {
            h_values.Add(table[i].h);
            errors.Add(table[i].R.Value);
        }
        var plot = new Plot();
        plot.Title($"error ~ h: {methodName.ToString()}");
        plot.XLabel("log10(h)");
        plot.YLabel("log10(R)");

        var log_h = h_values.Select(x => Math.Log10(x)).ToArray();
        var log_err = errors.Select(x => Math.Log10(x)).ToArray();

        var scatter = plot.Add.Scatter(log_h, log_err);
        scatter.LegendText = methodName.ToString();
        scatter.LineWidth = 2;
        scatter.MarkerSize = 5;
        scatter.Color = Colors.Red;

        var colors = new[] { Colors.Blue, Colors.Green};
        for (int i = 0; i < referenceSlopes.Length; i++)
        {
            double slope = referenceSlopes[i];
            double x1 = log_h[0];
            double x2 = log_h[^1];
            double y1 = log_err[0];
            double y2 = y1 + slope * (x2 - x1);

            var line = plot.Add.Line(x1, y1, x2, y2);
            line.LegendText = $"p = {slope}";
            line.LineWidth = 1;
            line.LineStyle.Pattern = LinePattern.Dashed;
            line.Color = colors[i % colors.Length];
        }

        plot.ShowLegend();
        plot.SavePng($"../../../intplots/{methodName.ToString()}.png", 800, 600);
    }

    public void RunAllCalculations()
    {
        // Метод средних прямоугольников
        var mpTable = ComputeIntegral(Methods.Midpoint, MidpointRule, A, B);
        PlotTask1Graph(Methods.Midpoint, mpTable, 2.0);
        int Nopt = FindOptimalN(mpTable, 2.0);
        ComputeIntegral(Methods.Midpoint, MidpointRule, A, B, Nopt);

        Console.WriteLine("=".PadLeft(200, '='));
        Console.WriteLine("=".PadLeft(200, '='));

        // Метод трапеций
        var trapTable = ComputeIntegral(Methods.Trapezoid, TrapezoidalRule, A, B);
        PlotTask1Graph(Methods.Trapezoid, trapTable, 2.0);
        Nopt = FindOptimalN(trapTable, 2.0);
        ComputeIntegral(Methods.Trapezoid, TrapezoidalRule, A, B, Nopt);

        Console.WriteLine("=".PadLeft(200, '='));
        Console.WriteLine("=".PadLeft(200, '='));

        // Метод Симпсона
        var simpTable = ComputeIntegral(Methods.Simpson, SimpsonRule, A, B);
        PlotTask1Graph(Methods.Simpson, simpTable, 4.0);
        Nopt = FindOptimalN(simpTable, 4.0);
        ComputeIntegral(Methods.Simpson, SimpsonRule, A, B, Nopt);

        Console.WriteLine("=".PadLeft(200, '='));
        Console.WriteLine("=".PadLeft(200, '='));

        // Метод Ньютона - Котса
        var ncTable = ComputeIntegral(Methods.NewtonCotes, SmallNewtonCotes, A, B);
        PlotTask2Graph(Methods.NewtonCotes, ncTable, new double[] {3.0, 4.0});
        Nopt = FindOptimN(ncTable, Methods.NewtonCotes);
        ComputeIntegral(Methods.NewtonCotes, SmallNewtonCotes, A, B, Nopt);

        Console.WriteLine("=".PadLeft(200, '='));
        Console.WriteLine("=".PadLeft(200, '='));

        // Метод Гаусса
        var gTable = ComputeIntegral(Methods.Gauss, SmallGaussian, A, B);
        PlotTask1Graph(Methods.Gauss, gTable, 6.0);
        Nopt = FindOptimN(gTable, Methods.Gauss);
        ComputeIntegral(Methods.Gauss, SmallGaussian, A, B, Nopt);
    }

    private int FindOptimalN(List<ResultTableRow> table, double p_theory)
    {
        ResultTableRow? optimalRow = null;
        for (int i = 2; i < table.Count; i++)
        {
            if (table[i].p.HasValue)
            {
                double diff = Math.Abs(table[i].p.Value - p_theory);
                if (diff < 0.05 * p_theory)
                {
                    optimalRow = table[i];
                    break;
                }
            }
        }
        Console.WriteLine($"\nOPTIMAL N: {optimalRow?.N}, STEP: {optimalRow?.h}\n");
        return optimalRow is null ? 1 : optimalRow.N;
    }
    private int FindOptimN(List<ResultTableRow> table, Methods m)
    {
        ResultTableRow? optimalRow = null;
        for (int i = 2; i < table.Count; i++)
        {
            if (table[i].p.HasValue)
            {
                if(m == Methods.NewtonCotes)
                {
                    if(2.7 < table[i].p && table[i].p <= 4)
                    {
                        optimalRow = table[i];
                        break;
                    }
                }
                else
                {
                    if (5.7 < table[i].p && table[i].p < 6.7)
                    {
                        optimalRow = table[i];
                        break;
                    }
                }
            }
        }
        Console.WriteLine($"\nOPTIMAL N: {optimalRow?.N}, STEP: {optimalRow?.h}\n");
        return optimalRow is null ? 1 : optimalRow.N;
    }

    public new void print()
    {
        RunAllCalculations();
    }

    enum Methods { Midpoint, Trapezoid, Simpson, NewtonCotes, Gauss}
}