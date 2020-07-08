# MNK
using System;

namespace MNK
{

    public class Matrix
    {
    
        public double[,] Args { get; set; }
        public int Row { get; set; }
        public int Col { get; set; }
        public Matrix(double[] x)
        {
            Row = x.Length;
            Col = 1;
            Args = new double[Row, Col];
            for (int i = 0; i < Args.GetLength(0); i++)
                for (int j = 0; j < Args.GetLength(1); j++)
                    Args[i, j] = x[i];
        }
        
        public Matrix(double[,] x)
        {
            Row = x.GetLength(0);
            Col = x.GetLength(1);
            Args = new double[Row, Col];
            for (int i = 0; i < Args.GetLength(0); i++)
                for (int j = 0; j < Args.GetLength(1); j++)
                    Args[i, j] = x[i, j];
        }
        
        public Matrix(Matrix other)
        {
            this.Row = other.Row;
            this.Col = other.Col;
            Args = new double[Row, Col];
            for (int i = 0; i < Row; i++)
                for (int j = 0; j < Col; j++)
                    this.Args[i, j] = other.Args[i, j];
        }
        
        public Matrix Transposition()
        {
            double[,] t = new double[Col, Row];
            for (int i = 0; i < Row; i++)
                for (int j = 0; j < Col; j++)
                    t[j, i] = Args[i, j];
            return new Matrix(t);
        }
        
        public static Matrix operator *(Matrix m, double k)
        {
            Matrix ans = new Matrix(m);
            for (int i = 0; i < ans.Row; i++)
                for (int j = 0; j < ans.Col; j++)
                    ans.Args[i, j] = m.Args[i, j] * k;
            return ans;
        }
        
        public static Matrix operator *(Matrix m1, Matrix m2)
        {
            if (m1.Col != m2.Row) throw new ArgumentException("Multiplication of these two matrices can't be done!");
            double[,] ans = new double[m1.Row, m2.Col];
            for (int i = 0; i < m1.Row; i++)
            {
                for (int j = 0; j < m2.Col; j++)
                {
                    for (int k = 0; k < m2.Row; k++)
                    {
                        ans[i, j] += m1.Args[i, k] * m2.Args[k, j];
                    }
                }
            }
            return new Matrix(ans);
        }
        
        private Matrix getMinor(int row, int column)
        {
            if (Row != Col) throw new ArgumentException("Matrix should be square!");
            double[,] minor = new double[Row - 1, Col - 1];
            
            for (int i = 0; i < this.Row; i++)
            {
                for (int j = 0; j < this.Col; j++)
                {
                    if ((i != row) || (j != column))
                    {
                    
                        if (i > row && j < column) minor[i - 1, j] = this.Args[i, j];
                        if (i < row && j > column) minor[i, j - 1] = this.Args[i, j];
                        if (i > row && j > column) minor[i - 1, j - 1] = this.Args[i, j];
                        if (i < row && j < column) minor[i, j] = this.Args[i, j];
                    }
                }
            }
            return new Matrix(minor);
        }
        
        public static double Determ(Matrix m)
        {
        
            if (m.Row != m.Col) throw new ArgumentException("Matrix should be square!");
            double det = 0;
            
            int length = m.Row;
            if (length == 1) det = m.Args[0, 0];
            if (length == 2) det = m.Argsusing System;
           
            if (length > 2)
                for (int i = 0; i < m.Col; i++)
                    det += Math.Pow(-1, 0 + i) * m.Args[0, i] * Determ(m.getMinor(0, i));
            return det;
        }
        public Matrix MinorMatrix()
        {
            double[,] ans = new double[Row, Col];

            for (int i = 0; i < Row; i++)
                for (int j = 0; j < Col; j++)
                    ans[i, j] = Math.Pow(-1, i + j) * Determ(this.getMinor(i, j));
            return new Matrix(ans);
        }
        public Matrix InverseMatrix()
        {
            if (Math.Abs(Determ(this)) <= 0.000000001) throw new ArgumentException("Inverse matrix does not exist!");
            double k = 1 / Determ(this);
            Matrix minorMatrix = this.MinorMatrix();
            return minorMatrix * k;
        }
        public class LSM
        {
            // Массивы значений Х и У задаются как свойства
            public double[] X { get; set; }
            public double[] Y { get; set; }
            // Искомые коэффициенты полинома в данном случае, а в общем коэфф. при функциях
            private double[] coeff;
            public double[] Coeff { get { return coeff; } }
            // Конструктор класса. Примает 2 массива значений х и у
            public LSM(double[] x, double[] y)
            {
                if (x.Length != y.Length) throw new ArgumentException("X and Y arrays should be equal!");
                X = new double[x.Length];
                Y = new double[y.Length];
                for (int i = 0; i < x.Length; i++)
                {
                    X[i] = x[i];
                    Y[i] = y[i];
                }
            }
            // Метод Наименьших Квадратов
            // В качестве базисных функций используются степенные функции y = a0 * x^0 + a1 * x^1 + ... + am * x^m
            public void Polynomial(int m)
            {
                if (m <= 0) throw new ArgumentException("Порядок полинома должен быть больше 0");
                if (m >= X.Length) throw new ArgumentException("Порядок полинома должен быть на много меньше количества точек!");
                // массив для хранения значений базисных функций
                double[,] basic = new double[X.Length, m + 1];
                // заполнение массива для базисных функций
                for (int i = 0; i < basic.GetLength(0); i++)
                    for (int j = 0; j < basic.GetLength(1); j++)
                        basic[i, j] = Math.Pow(X[i], j);
                // Создание матрицы из массива значений базисных функций(МЗБФ)
                Matrix basicFuncMatr = new Matrix(basic);
                // Транспонирование МЗБФ
                Matrix transBasicFuncMatr = basicFuncMatr.Transposition();
                // Произведение транспонированного  МЗБФ на МЗБФ
                Matrix lambda = transBasicFuncMatr * basicFuncMatr;
                // Произведение транспонированого МЗБФ на следящую матрицу 
                Matrix beta = transBasicFuncMatr * new Matrix(Y);
                // Решение СЛАУ путем умножения обратной матрицы лямбда на бету
                Matrix a = lambda.InverseMatrix() * beta;
                // Присвоение значения полю класса 
                coeff = new double[a.Row];
                for (int i = 0; i < coeff.Length; i++)
                {
                    coeff[i] = a.Args[i, 0];
                }
            }
        }
        class Program
        {
            static void Main(string[] args)
            {
                //Исходные данные
                double[] x = new double[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
                double[] y = new double[] { 8, 8, 0, 5, 5, 2, 4, 3, 8, 9, 6 };
                // Создание экземляра класса LSM
                LSM myReg = new LSM(x, y);
                // Апроксимация заданных значений линейным полиномом
                myReg.Polynomial(1);
                // Вывод коэффициентов а0 и а1
                for (int i = 0; i < myReg.Coeff.Length; i++)
                {
                    Console.WriteLine(myReg.Coeff[i]);
                }
                Console.WriteLine();
                Console.ReadLine();
            }
        }
    }
}
