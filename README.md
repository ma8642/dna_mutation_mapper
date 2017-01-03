# dna_mutation_mapper
Software that discovers and displays order of mutation patterns in Breast Cancer fragments

This software was created over 10 weeks in the Schatz Lab at Cold Spring Harbor Laboratory.  It is a DNA analysis and visualization algorithm using Python and Graphviz software. In order to create this software, I created a new Graph class that implemented nodes with start and end “ports” which better represents DNA fragments. The program was able to reveal a clear and detailed representation of the mutation process of a convoluted area in the HER2+ breast cancer genome. Please see abstract.txt and ALFORD.pdf for further details on project scope and impact.

My mentor for that summer, Maria Nattestad, has recently published an open-use, sophisticated DNA mapping software, SplitThreader (http://splitthreader.com/), which uses the same graph structure that we came up with, and ideas inspired by our discussions about my project.

#To Run
1.  Download GraphViz (http://www.graphviz.org/), main.py, and graph_template.dot

2.  Open the terminal and navigate to the file that contains main.py and graph_template.dot

3.  Type into terminal:  python main.py filename

4.  A .png file with the filename specified will open revealing the DNA mutation mapping of the Her2 region.

![alt tag](https://github.com/ma8642/dna_mutation_mapper/blob/master/test.png)

You will also see a legend in the terminal describing the paths.

![alt tag](https://github.com/ma8642/dna_mutation_mapper/blob/master/legend.png)

This indicates that path 5, for example, can be attained by tracing along the green lines from chromosome fragment Y to chromosome fragment B17 to chromosome fragment Z, and that the minimum copy number along that path is 121.

