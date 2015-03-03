using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace GeneticsLab
{
    class Cell
    {
        public int leastCost { get; set; }
        public Cell previous { get; set; }

        public Cell(int leastCost, Cell previous)
        {
            this.leastCost = leastCost;
            this.previous = previous;
        }

    }
}
