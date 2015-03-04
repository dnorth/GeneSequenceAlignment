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
        private int MaxCharactersToExtract = 101;
        private int indel = 5;
        private int match = -3;
        private int sub = 1;
        private Cell[,] table;
        private int[] row1;
        private int[] row2;
        string first;
        string second;

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
            //Only get the top of the diagonal
            if (columnInTable < rowInTable)
            {
                int sizeA = 0;
                int sizeB = 0;

                //Add the arbitary value at the beginning to fill as 0
                string secA = '-' + sequenceA.Sequence;
                string secB = '-' + sequenceB.Sequence;

                //These are different cases to call BuildTable() based on whether or not sequenceA.Sequence's or sequenceB.Sequence's length are >= MaxCharactersToAlign 
                //We will create 2 arrays of size n, n being the y-value of the table (Meaning this is still in O(n) space because O(2n) is in O(n) ) 
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

                //In the end, we only care about the very last value in 
                return row2[sizeA-1];
            }
            else
            {
                return 0;
            }
        }

        //sequenceA - the sequence of the table in the vertical direction
        //sequenceB - the sequence of the table in the horizontal direction
        //arrSize - the size of the array, needed for Array.Copy()
        public void BuildTable(string sequenceA, string sequenceB, int arrSize)
        {
            //This double for loop will take O(n^2) time to run, but is the way I will loop through all posibilities in the table
            for (int j = 0; j < sequenceB.Length; j++)
            {
                for (int i = 0; i < sequenceA.Length; i++)
                {
                    //if we are working on the very first column, we will fill row1. This is not necessary but convenient for debugging purposes
                    if (j == 0)
                    {
                        row1[i] = Score(i, j, sequenceA[i], sequenceB[j]);
                    }
                    else
                    {
                        //Else we will fill row2. We're just filling the one column of the table at a time
                        row2[i] = Score(i, j, sequenceA[i], sequenceB[j]);
                    }
                }
                if (j != 0)
                {
                    //Once the column is filled, we will copy it over to the first column to reuse space. Then row2 will be written over with the next column of the table
                    Array.Copy(row2, row1, arrSize);
                }
            }

            //PrintTable();
        }

        //i- the y-value of the table
        //j- the x-value of the table
        //charA- the character at the x-position
        //charB- the character at the y-position
        public int Score(int i, int j, char charA, char charB)
        {
            //Every Possible Case of the recurrence relation

            //Case 1: 0 x-value and 0 y-value of the table
            if (i == 0 && j == 0)
            {
                return 0;
            }
            else if (i == 0 && j > 0)
            {
                //Case 2: if the y-value is 0 but the x-value is greater than 0, get the previous x-value with the same y-value and add an indel cost to it
                return row1[i] + indel;
            }
            else if (i > 0 && j == 0)
            {
                //Case 3: if the x-value is 0 but the y-value is greater than 0, get the same x-value (we have to be working on row1) with the previous y-value and add an indel cost to it
                return row1[i - 1] + indel;
            }
            else if (i > 0 && j > 0)
            {
                //Case 4: if the x-value is greater than 0 and the y-value is greater than 0, we get the minimum cost of three cases
                int minScore = int.MaxValue;

                //Case 4-1: get the previous x-value with the same y-value and add an indel cost to it
                int first = row1[i] + indel;

                //Case 4-2: get the previous x-value and the previous y-value and, if charA and charB match, add a match cost to it. If not, add a sub cost.
                int second;
                if (charA == charB)
                {
                    second = row1[i - 1] + match;
                }
                else
                {
                    second = row1[i - 1] + sub;
                }

                //Case 4-3: get the previous y-value and same x-value and add an indel cost to it
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

                //return the minimum of the 3
                return minScore;
            }
            else
            {
                //This should never happen... but just in case something goes terrible wrong... you'll notice with the number jump.
                return 10000000;
            }
        }

        //Used for testing to print the table 
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

         /// <summary>
        /// this is the function you implement.
        /// </summary>
        /// <param name="sequenceA">the first sequence</param>
        /// <param name="sequenceB">the second sequence, may have length not equal to the length of the first seq.</param>
        /// <param name="rowInTable">this particular alignment problem will occupy a cell in this row the result table.</param>
        /// <param name="columnInTable">this particular alignment will occupy a cell in this column of the result table.</param>
        /// <returns>the alignment score for sequenceA and sequenceB.  The calling function places the result in entry rowInTable,columnInTable
        /// of the ResultTable</returns>
        public string ExtractAlign(GeneSequence sequenceA, GeneSequence sequenceB, int rowInTable, int columnInTable)
        {

            //Only get the top of the diagonal
            if (columnInTable > rowInTable)
            {
                int sizeA = 0;
                int sizeB = 0;

                //Add the arbitary value at the beginning to fill as 0
                string secA = '-' + sequenceA.Sequence;
                string secB = '-' + sequenceB.Sequence;

                //These are different cases to call BuildTable() based on whether or not sequenceA.Sequence's or sequenceB.Sequence's length are >= MaxCharactersToAlign 
                //This time, because we're not worried about space, we will build a 2D array of sizeA and sizeB from the Cell class that I created
                if (secA.Length >= MaxCharactersToExtract && secB.Length >= MaxCharactersToExtract)
                {
                    sizeA = MaxCharactersToExtract;
                    sizeB = MaxCharactersToExtract;
                    table = new Cell[sizeA, sizeB];

                    Extract(secA.Substring(0, MaxCharactersToExtract), secB.Substring(0, MaxCharactersToExtract));

                }
                else if (secA.Length < MaxCharactersToExtract && secB.Length >= MaxCharactersToExtract)
                {
                    sizeA = secA.Length;
                    sizeB = MaxCharactersToExtract;
                    table = new Cell[sizeA, sizeB];

                    Extract(secA, secB.Substring(0, MaxCharactersToExtract));
                }
                else if (secA.Length >= MaxCharactersToExtract && secB.Length < MaxCharactersToExtract)
                {
                    sizeA = MaxCharactersToExtract;
                    sizeB = secB.Length;
                    table = new Cell[sizeA, sizeB];

                    Extract(secA.Substring(0, MaxCharactersToExtract), secB);
                }
                else
                {
                    sizeA = secA.Length;
                    sizeB = secB.Length;
                    table = new Cell[sizeA, sizeB];

                    Extract(secA, secB);
                }

                //When we have built the entire table we will get the final cell and trace back all its back-pointers, adding the character value given to it by the Cell class to the strings first and second
                //We will then concat the two strings so they show up nicely in the textbox and return the contatenation 
                convertToString(table[sizeA - 1, sizeB - 1]);
                first = first + "\r\n" + second;
                return first;
            }
            else
            {
               return "row: " +  rowInTable + " column: " +  columnInTable;
            }
        }

        //sequenceA - the sequence of the table in the vertical direction
        //sequenceB - the sequence of the table in the horizontal direction
        public void Extract(string sequenceA, string sequenceB)
        {
            //This double for loop will take O(n^2) time to run, but is the way I will loop through all posibilities in the table
            for (int i = 0; i < sequenceA.Length; i++)
            {
                for (int j = 0; j < sequenceB.Length; j++)
                {
                    //Save the Cell we get back from the function into the table
                    table[i, j] = ScoreExtract(i, j, sequenceA[i], sequenceB[j]);
                }

            }

            PrintTable();
        }

        //i- the x-value of the table
        //j- the y-value of the table
        //charA- the character at the x-position
        //charB- the character at the y-position
        public Cell ScoreExtract(int i, int j, char charA, char charB)
        {
            //Here we have all the cases of the original Scoring algorithm, we're just saving things in O(n^2) space because of the table. This algorithm is very similar the the one above
            //with a few things added to it because we are returning cells, and cells are returning their leastCost, backpointer, character for lining up the x-axis, character for lining up the y-axis
            if (i == 0 && j == 0)
            {
                return new Cell(0, null, '0', '0');
            }
            else if (i == 0 && j > 0)
            {
                return new Cell(table[i, j - 1].leastCost + indel, table[i, j - 1], '-', charB);
            }
            else if (i > 0 && j == 0)
            {
                int test = table[i - 1, j].leastCost + indel;
                return new Cell(table[i - 1, j].leastCost + indel, table[i - 1, j], charA, '-');
            }
            else if (i > 0 && j > 0)
            {
                int minScore = int.MaxValue;
                Cell minCell = null;
                char retA = '0';
                char retB = '0';

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
                retA = charA;
                retB = '-';


                if (second < minScore)
                {
                    minScore = second;
                    minCell = table[i - 1, j - 1];
                    retA = charA;
                    retB = charB;
                }
                if (third < minScore)
                {
                    minScore = third;
                    minCell = table[i, j - 1];
                    retA = '-';
                    retB = charB;
                }

                return new Cell(minScore, minCell, retA, retB);
            }
            else
            {
                return null;
            }
        }

        //finalCell - On the first call it should be the last cell in the table, on every recursive pass, it's the previous pointer of the original.
        public void convertToString(Cell finalCell)
        {
            //As long as the previous pointer isn't null, add the x-characters to a string and the y-characters to another string, recurse with the previous pointer
            if(finalCell.previous != null)
            {
                first = finalCell.aChar + first;
                second = finalCell.bChar + second;
                convertToString(finalCell.previous);
            }
        }
    }
}
