using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GeneticsLab
{
    //Class used for extraction. Gives us the least cost of the cell, the previous pointer, the character on the x-axis belonging to the alignment, and the character on the y-axis belonging to the alignment
    class Cell
    {
        public int leastCost { get; set; }
        public Cell previous { get; set; }
        public char aChar { get; set; }
        public char bChar { get; set; }

        public Cell(int leastCost, Cell previous, char aChar, char bChar)
        {
            this.leastCost = leastCost;
            this.previous = previous;
            this.aChar = aChar;
            this.bChar = bChar;
        }

    }
}
