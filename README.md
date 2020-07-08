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
