using System;
using System.Collections.Generic;
using System.Text;

namespace GeneticsLab
{
    class PairWiseAlign
    {

        /// <summary>
        /// Align only 5000 characters in each sequence.
        /// </summary>
        private int MaxCharactersToAlign = 5001;
        private int indel = 5;
        private int match = -3;
        private int sub = 1;
        private Cell[,] table;
        private int[] row1;
        private int[] row2;

        /// <summary>
        /// this is the function you implement.
        /// </summary>
        /// <param name="sequenceA">the first sequence</param>
        /// <param name="sequenceB">the second sequence, may have length not equal to the length of the first seq.</param>
        /// <param name="resultTableSoFar">the table of alignment results that has been generated so far using pair-wise alignment</param>
        /// <param name="rowInTable">this particular alignment problem will occupy a cell in this row the result table.</param>
        /// <param name="columnInTable">this particular alignment will occupy a cell in this column of the result table.</param>
        /// <returns>the alignment score for sequenceA and sequenceB.  The calling function places the result in entry rowInTable,columnInTable
        /// of the ResultTable</returns>
        public int Align(GeneSequence sequenceA, GeneSequence sequenceB, ResultTable resultTableSoFar, int rowInTable, int columnInTable)
        {
            
            if (columnInTable < rowInTable)
            {
                int sizeA = 0;
                int sizeB = 0;

                string secA = '-' + sequenceA.Sequence;
                string secB = '-' + sequenceB.Sequence;

                if (secA.Length >= MaxCharactersToAlign && secB.Length >= MaxCharactersToAlign)
                {
                    sizeA = MaxCharactersToAlign;
                    row1 = new int[sizeA];
                    row2 = new int[sizeA];
                    BuildTable(secA.Substring(0, MaxCharactersToAlign), secB.Substring(0, MaxCharactersToAlign), sizeA);

                }
                else if (secA.Length < MaxCharactersToAlign && secB.Length >= MaxCharactersToAlign)
                {
                    sizeA = secA.Length;
                    sizeB = MaxCharactersToAlign;
                    row1 = new int[sizeA];
                    row2 = new int[sizeA];

                    BuildTable(secA, secB.Substring(0, MaxCharactersToAlign), sizeA);
                }
                else if (secA.Length >= MaxCharactersToAlign && secB.Length < MaxCharactersToAlign)
                {
                    sizeA = MaxCharactersToAlign;
                    sizeB = secB.Length;
                    row1 = new int[sizeA];
                    row2 = new int[sizeA];

                    BuildTable(secA.Substring(0, MaxCharactersToAlign), secB, sizeA);
                }
                else
                {
                    sizeA = secA.Length;
                    sizeB = secB.Length;
                    row1 = new int[sizeA];
                    row2 = new int[sizeA];

                    BuildTable(secA, secB, sizeA);
                }
                return row2[sizeA-1];
                /*
                int sizeA = 0;
                int sizeB = 0;

                string secA = '-' + sequenceA.Sequence;
                string secB = '-' + sequenceB.Sequence;

                if (secA.Length >= MaxCharactersToAlign && secB.Length >= MaxCharactersToAlign)
                {
                    sizeA = MaxCharactersToAlign;
                    sizeB = MaxCharactersToAlign;
                    table = new Cell[sizeA, sizeB];

                    Extract(secA.Substring(0, MaxCharactersToAlign), secB.Substring(0, MaxCharactersToAlign));

                }
                else if (secA.Length < MaxCharactersToAlign && secB.Length >= MaxCharactersToAlign)
                {
                    sizeA = secA.Length;
                    sizeB = MaxCharactersToAlign;
                    table = new Cell[sizeA, sizeB];

                    Extract(secA, secB.Substring(0, MaxCharactersToAlign));
                }
                else if (secA.Length >= MaxCharactersToAlign && secB.Length < MaxCharactersToAlign)
                {
                    sizeA = MaxCharactersToAlign;
                    sizeB = secB.Length;
                    table = new Cell[sizeA, sizeB];

                    Extract(secA.Substring(0, MaxCharactersToAlign), secB);
                }
                else
                {
                    sizeA = secA.Length;
                    sizeB = secB.Length;
                    table = new Cell[sizeA, sizeB];

                    Extract(secA, secB);
                }
                return table[sizeA - 1, sizeB - 1].leastCost;*/
            }
            else
            {
                return 0;
            }
        }

        public void BuildTable(string sequenceA, string sequenceB, int arrSize)
        {
            for (int j = 0; j < sequenceB.Length; j++)
            {
                for (int i = 0; i < sequenceA.Length; i++)
                {
                    if (j == 0)
                    {
                        row1[i] = Score(i, j, sequenceA[i], sequenceB[j]);
                    }
                    else
                    {
                        row2[i] = Score(i, j, sequenceA[i], sequenceB[j]);
                    }
                }
                if (j != 0)
                {
                    Array.Copy(row2, row1, arrSize);
                }
            }

            //PrintTable();
        }

        public void Extract(string sequenceA, string sequenceB)
        {
            for (int i = 0; i < sequenceA.Length; i++)
            {
                for (int j = 0; j < sequenceB.Length; j++)
                {
                    table[i, j] = ScoreExtract(i, j, sequenceA[i], sequenceB[j]);
                }

            }

            //PrintTable();
        }

        public int Score(int i, int j, char charA, char charB)
        {
            if (i == 0 && j == 0)
            {
                return 0;
            }
            else if (i == 0 && j > 0)
            {
                return row1[i] + indel;
            }
            else if (i > 0 && j == 0)
            {
                return row1[i - 1] + indel;
            }
            else if (i > 0 && j > 0)
            {
                int minScore = int.MaxValue;

                int first = row1[i] + indel;

                int second;
                if (charA == charB)
                {
                    second = row1[i - 1] + match;
                }
                else
                {
                    second = row1[i - 1] + sub;
                }

                int third = row2[i - 1] + indel;

                minScore = first;

                if (second < minScore)
                {
                    minScore = second;
                }
                if (third < minScore)
                {
                    minScore = third;
                }

                return minScore;
            }
            else
            {
                return 10000000;
            }
        }

        public Cell ScoreExtract(int i, int j, char charA, char charB)
        {
            if (i == 0 && j == 0)
            {
                return new Cell(0, null);
                //table[i, j] = new Cell(0, null);
            }
            else if (i == 0 && j > 0)
            {
                return new Cell(table[i, j - 1].leastCost + indel, table[i, j - 1]);
                //table[i, j] = new Cell(Score(i, j - 1), table[i, j - 1]);
            }
            else if (i > 0 && j == 0)
            {
                int test = table[i - 1, j].leastCost + indel;
                return new Cell(table[i - 1, j].leastCost + indel, table[i - 1, j]);
            }
            else if (i > 0 && j > 0)
            {
                int minScore = int.MaxValue;
                Cell minCell = null;

                int first = table[i - 1, j].leastCost + indel;

                int second;
                if (charA == charB)
                {
                    second = table[i - 1, j - 1].leastCost + match;
                }
                else
                {
                    second = table[i - 1, j - 1].leastCost + sub;
                }

                int third = table[i, j - 1].leastCost + indel;

                minScore = first;
                minCell = table[i - 1, j];

                if (second < minScore)
                {
                    minScore = second;
                    minCell = table[i - 1, j - 1];
                }
                if (third < minScore)
                {
                    minScore = third;
                    minCell = table[i, j - 1];
                }

                return new Cell(minScore, minCell);
            }
            else
            {
                return null;
            }
        }

        public void PrintTable()
        {
            int rowLength = table.GetLength(0);
            int colLength = table.GetLength(1);

            Console.WriteLine("\n------------------------------PRINTING TABLE -------------------- \n");
            for (int i = 0; i < rowLength; i++)
            {
                for (int j = 0; j < colLength; j++)
                {
                    Console.Write(string.Format("{0} ", table[i, j].leastCost));
                }
                Console.Write(Environment.NewLine + Environment.NewLine);
            }
            Console.WriteLine("\n------------------------------FINISHED PRINTING --------------------\n");

        }
    }
}
