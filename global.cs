global using System;
global using System.Text;
global using System.Linq;
global using System.Collections.Generic;
global using MathNet;
global using MathNet.Numerics;
global using MathNet.Numerics.LinearAlgebra;
global using MathNet.Numerics.RootFinding;
global using MathNet.Numerics.Integration;
global using MathNet.Numerics.Differentiation;
global using MathNet.Numerics.Optimization;
global using AngouriMath;
global using AngouriMath.Extensions;
global using ScottPlot;

//public static class MultiDimNewtonRaphson
//{
//    // Решает систему f(x) = 0 методом Ньютона
//    // f: функция (вектор -> вектор)
//    // x0: начальное приближение (Vector<double>)
//    // tol: точность выхода
//    // maxIter: максимум итераций
//    // jacStep: шаг для численной производной
//    public static Vector<double> FindRoot(
//        Func<Vector<double>, Vector<double>> f,
//        Vector<double> x0,
//        double tol = 1e-7,
//        int maxIter = 100,
//        double jacStep = 1e-8)
//    {
//        var n = x0.Count;
//        var x = x0.Clone();

//        for (int iter = 0; iter < maxIter; iter++)
//        {
//            var fx = f(x);

//            // Проверка на выход
//            if (fx.L2Norm() < tol)
//                return x;

//            // Якобиан численно
//            var J = Matrix<double>.Build.Dense(n, n);
//            for (int j = 0; j < n; j++)
//            {
//                var xPlus = x.Clone();
//                xPlus[j] += jacStep;
//                var dF = (f(xPlus) - fx) / jacStep;
//                for (int i = 0; i < n; i++)
//                    J[i, j] = dF[i];
//            }

//            // Решаем J * dx = -f(x)
//            var dx = J.Solve(-fx);

//            x += dx;
//        }

//        throw new Exception("Newton-Raphson did not converge");
//    }
//}
