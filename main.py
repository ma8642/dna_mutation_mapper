#  Marley Alford code initially developed over 10 weeks, summer 2015

import string
import re
import random
import sys

# ''''''''''''''''''''EDGES''''''''''''''''''''''''
#  METHODS
#  str() -> prints Node$port [label=label, color=color, penwidth=pw, fontcolor=color]

class Edge(object):

    def __init__(self, node, w, color):  # inputs: parent node, edgeweight, color
        self.next = node  
        self.label = w
        # self.pos = node.pos  # this is the DNA position of
        if color == "":
            self.color = "black"
        else:
            self.color = color
        self.fontcolor = "red2"

    def str(self):
        if self.label > 400:
            return self.next.str(), "[label=\"%d\", color=\"%s\", penwidth=40, fontcolor=%s]" %(self.label, self.color, self.fontcolor)
        else:
            return self.next.str(), "[label=\"%d\", color=\"%s\", penwidth=%d, fontcolor=%s]" %(self.label, self.color, float(self.label)/10, self.fontcolor) # float(self.label)/10)


# ''''''''''''''''''''PORTS'''''''''''''''''''''''''
#  METHODS
#  give() ->  gives the children and edges of one port to another.  input:  other port to give attributes to.
#  add() -> connects one port to another by listing each in the other's children list, and making edges between them.
#  findEdge() -> given two ports, returns the edge between them.
#  getIndex() -> For graphing purposes.  The indices refer to ports start and end in a node.
#  minus() -> subtracts value sub from an edge between two ports.
#  str() -> prints Node$port

class Port(object):
    def __init__(self, name):
        self.n = name
        self.node = "NoParent"
        self.children = []  # child ports that connect to this port
        self.edges = []
        self.id = 0    # for graphing purposes
        self.pos = 0
        self.chrom = 0

    def give(self, other):
        self.id = other.id
        for child in self.children:
            if child not in other.children:
                other.children.append(child)
        for edge in self.edges:
            if edge not in other.edges:
                other.edges.append(edge)

    def add(self, other, w, color):
        self.children.append(other)    # puts the other into list
        other.children.append(self)  # puts self into children of other
        e = Edge(other, w, color)
        self.edges.append(e)   # adds e to list of edges
        e2 = Edge(self, w, color)
        other.edges.append(e2) # adds e to list of edges of other

    def findEdge(self, childport):
        for edge in self.edges:
            if edge.next.str() == childport.str():
                return edge

    def getIndex(self):           # for graphing purposes
        if self.n == "start":
            ind = "f1"
        else: # if self.n == "end"
            ind = "f2"
        return ind

    def minus(self, other, sub):
        edge = self.findEdge(other)
        same = other.findEdge(self)
        edge.label -= sub
        same.label -= sub

    def str(self):
        return "{}${}" .format(self.node, self.n)

# ''''''''''''''''''''NODES''''''''''''''''''''''''
#  METHODS
#  add() ->  gives the children and edges of one port to another.  input:  other port to give attributes to.
#  jump() -> connects one port to another by listing each in the other's children list, and making edges between them.
#  changeID() -> given two ports, returns the edge between them.

class Node(object):
    def __init__(self, name):
        self.n = name
        self.atts = {}  # node attribute, seq, chrom, pos
        self.ports = {}  # string: port object
        self.id = 0  # for graphing purposes

    def add(self, port):
        self.ports[port.n] = port
        port.node = self.n
        port.id = self.id
        if "chrom" in self.atts:
            port.chrom = self.atts["chrom"]
        if "DNApos" in self.atts:
            if port.n == "start":
                port.pos = self.atts["DNApos"][0]
            elif port.n == "end":
                port.pos = self.atts["DNApos"][1]

    def jump(self, port):    # jumps from current port to another in a node
        dest = ""
        if len(self.ports) > 1:
            if port.n == "end":
                dest = self.ports["start"]
            else:
                dest = self.ports["end"]
        else:
            dest = self.ports[port.n]
        return dest

    def changeID(self, ID):   # for graphing purposes
        self.id = ID
        for port in self.ports.values():
            port.id = ID

# ''''''''''''''''''''GRAPHS''''''''''''''''''''''''
#  METHODS
#  removeEdge() ->  gives the children and edges of one port to another.  input:  other port to give attributes to.
#  add() -> connects one port to another by listing each in the other's children list, and making edges between them.
#  findSink() -> given two ports, returns the edge between them.
#  copy() -> returns a new graph with all the attributes of the original
#  str() -> prints out all the edges in the graph.
#  
#  makeGraphBlueprint() -> writes a .dot file that can be made into a diagram of the graph.
#  readIn() -> takes a .dot file as input and turns it into a graph.
#  paths_from_to() -> recusively finds all possible paths from Source to Sink in a graph.  returns a list of all those paths.
#  findAllPaths() -> feeds parameters to paths_from_to() and the returns the list of all paths.
#  maxPathLen() -> finds the longest uninterrupted sequence in a given path.  Returns length of that path.,
#  allMaxPathLengths() -> Finds the longest uninterrupted sequence and min weighted edge of every path.  Returns two lists, one with paths, and one with the respective lengths of their longest uninterrupted sequence.
#  cap() -> Finds the minimum weighted edge of a graph.
#  sub() -> Subtracts the minimum weighted edge from all edges in a path.
#  findBestPaths() -> Finds thinnest point on chosen (longest) path (using cap()), and subtracts it (using sub()), printing out longest paths as it goes.
#  draw() -> Takes longest paths as input and then feeds this information to makeColorCode()
#  makeColorCode() -> makeColorCode() creates a graph and then calls makeGraphBlueprint()



class Graph(object):
    def __init__(self):
        self.Source = Port("Portal")
        self.Sink = self.Source
        self.nodes = {}  # creates dictionary of nodes in graph
        self.ports = {}  # creates dictionary of ports in graph

    def removeEdge(self, port1, port2):   # input the two nodes that the edge connects
        edge = port1.findEdge(port2)
        edge2 = port2.findEdge(port1)   # we have to find this copy of the edge too.
        port1.edges.remove(edge)
        port2.edges.remove(edge2)

    def add(self, node):
        if node.n not in self.nodes:
            self.nodes[node.n] = node  # append to dictionary
            if node.n == "Portal":  # If this is the first node we're inserting, then it's the Portal
                self.Source = node.ports["end"]  # this assigns source.  You'll have to assign sink in readin
            for port in node.ports.values():
                self.ports[port.str()] = port
        else:
            ourNode = self.nodes[node.n]
            for port in node.ports.values():  # put ports that weren't there already
                if port.n in ourNode.ports:           # if a particular port is already there
                    port.give(ourNode.ports[port.n])  # give attribute to existing node
                else:                               # if it's not already there
                    self.ports[port.str()] = port # put it into graph.ports
                    ourNode.ports[port.n] = port  # and into node.ports  (adding it to ourNode means that we are modifying the original node and not adding a new one.)

    def copy(self):
        newName = Graph()
        newName.Source = self.Source
        newName.Sink = self.Sink
        newName.nodes = self.nodes
        newName.ports = self.ports
        return newName

    def str(self):
        printed =[]  # keeps track of port combinations already printed so we don't print redundant things.
        print "Graph:"
        for key in self.nodes:
            if len(self.nodes[key].att) != 0:
                print key, self.nodes[key].att
        for port in self.ports.values():
            for edge in port.edges:
                p, ed = edge.str()  # extracting port "p" so we can compare it
                if port.str() + p not in printed and p + port.str() not in printed:  # checking if port combination to print hasn't already been printed
                    print "{} -> {} {}" .format(port.str(), p, ed)
                    printed.append(port.str() + p)

    def makeGraphBlueprint(self, filename):
        text_file = open("{}.dot".format(filename), "w" )
        text_file.write("graph structs\n")
        text_file.write("{label = \"DNA Fragments Mapping\"\n")
        text_file.write("\tgraph [rankdir = \"LR\" bgcolor = \"white\" style = \"filled\"];\n\tnode [shape = \"record\" style = \"filled\"];\n\tedge [label=\"Edge\" penwidth=12 fontcolor=\"red\"];\n")

        for node in self.nodes.values():
            if node.n != "Portal":
                text_file.write("\t\tstruct%d [label = \"%s |{<f1> start |<f2> end}\" shape = record fillcolor = \"goldenrod1\" pos = \"%s,%s!\"];\n" %(node.id, node.n, node.atts["x"], node.atts["y"])) 

        printed =[]  # keeps track of port combinations already printed so we don't print redundant things.
        for port in self.ports.values():
            # we want to not makeGraphBlueprint portal or any of the edges attached to it
            for edge in port.edges:
                p1, ed = edge.str()  # extracting port "p1" and separating it from edge "ed"
                port2 = edge.next
                if port.node == "Portal" or port2.node == "Portal":  # don't makeGraphBlueprint these
                        continue
                else:
                    if port.str() + port2.str() + ed not in printed and port2.str() + port.str() + ed not in printed:  # checking if port combination to print hasn't already been printed'''
                        text_file.write("\t struct{}:{} -- struct{}:{} {}; \n" .format(port.id, port.getIndex(), port2.id, port2.getIndex(), ed))
                        # print port.str(), edge.str(), "\t struct{}:{} -- struct{}:{} {}; \n" .format(port.id, empty, p2.id, nope, ed)   # DEBUG
                        printed.append(port.str() + port2.str() + ed)   
        text_file.write("}")
        text_file.close()

    def readIn(self, filename):
        f = open(filename, "r")   # open file
        nodesToAdd = []
        deets = {}
        DNApos = {}
        for line in f:
            if "node" in line:   # finding node attributes
                line = line.strip().split()
                name = line[1]
                x = line[4].replace(",","")
                y = line[7].replace("]","").replace(",","")
                DNAx = 0
                DNAy = 0
                chromosome = line[10].replace(",","")
                if len(line) > 13:     # DNApos gives the node positon in the chromosome
                    DNAx = line[13].replace(",","").replace("]","")
                    if len(line) > 16:
                        DNAy = line[16].replace(",","").replace("]","")
                deets[name] = x, y, DNAx, DNAy, chromosome   # Adds attributes to node attributes dictionary

            if "--" in line:    # look for lines with edge information
                line = line.strip().split()   # divide it up into ['start', '-', 'end [edge attribute]']
                xid = line[0]
                yid = line[2].replace("\\","")
                weight = 0
                color = ""
                if len(line) > 3:    # if line has labels and colors
                    for field in line:                  # add label attribute to edge
                        if 'label' in field:
                            regex = '.*\"(.*?)\".*'
                            matches = re.search(regex, field)
                            weight = int(matches.group(1))
                        if 'color' in field:            # add color attribute to edge
                            regex = '.*\"(.*?)\".*'
                            matches2 = re.search(regex, field)
                            color = str(matches2.group(1))
                        if 'penwidth' in field:         # add label attribute to edge (if 'label' isn't specified, this is backup)
                            regex = '.*\=(.*?)\].*'
                            matches3 = re.search(regex, field)
                            weight = int(matches3.group(1))
                        else:
                            continue
                node1 = xid.split("$")[0]
                node2 = yid.split("$")[0]
                port1 = xid.split("$")[1]
                port2 = yid.split("$")[1]

                port1 = Port(port1)
                node1 = Node(node1)
                x, y, DNAx, DNAy, chromosome = deets[node1.n]
                node1.atts["x"] = x
                node1.atts["y"] = y
                node1.atts["DNApos"] = DNAx, DNAy
                node1.atts["chrom"] = chromosome
                node1.add(port1)

                port2 = Port(port2)
                node2 = Node(node2)
                x, y, DNAx, DNAy, chromosome = deets[node2.n]
                node2.atts["x"] = x
                node2.atts["y"] = y
                node2.atts["DNApos"] = DNAx, DNAy
                node2.atts["chrom"] = chromosome
                node2.add(port2)
                port1.add(port2, weight, color)

                nodesToAdd.append(node1)                # and add node to graph
                nodesToAdd.append(node2)
        i = 0
        IDs = {}
        for node in nodesToAdd:    # assign ID to node and then add it to graph
            if node.n in IDs:
                node.changeID(IDs[node.n])   # if node id is already in IDs, change it to that
                self.add(node)
            else:
                IDs[node.n] = i   # if not, make an id for that node
                node.changeID(i)
                self.add(node)
            i += 1 
        # finally, assign Sink
        portal = self.nodes["Portal"]
        self.Sink = portal.ports["start"]       

    def paths_from_to(self, noReturn, source, dest, allpaths, path_so_far):
        if source.str() == dest.str():
            if path_so_far not in allpaths:
                allpaths.append(path_so_far)
        else:
            parentNode = self.nodes[source.node]
            if parentNode.n != "Portal":
                source = parentNode.jump(source)
            for port in source.children:
                if port.str() != noReturn:
                    self.paths_from_to(noReturn, port, dest, allpaths, path_so_far + [port.str()])

    def findAllPaths(self):
        source = self.Source
        origSource = source.str()
        dest = self.Sink
        allpaths = []
        path_so_far = [source.str()]
        self.paths_from_to(origSource, source, dest, allpaths, path_so_far)
        return allpaths

    def maxPathLen(self, path):
        currentChr = ""
        maxSeqLen = 0
        curSeqLen = 0
        for port in path:
            p = self.ports[port]  # next port we want to add
            node = self.nodes[p.node]
            start = int(node.atts["DNApos"][0])
            stop = int(node.atts["DNApos"][1])
            nlength = stop - start      # length of node that port is on
            if p.chrom == currentChr:
                curSeqLen += nlength
            else:
                currentChr = p.chrom
                curSeqLen = nlength

            if maxSeqLen < curSeqLen:
                maxSeqLen = curSeqLen
        return maxSeqLen

    def allMaxPathLengths(self):    # calls maxPathLen() 
        everyLength = []
        allPaths = self.findAllPaths()
        i = 0
        everyPath = [] # we have to make this because there might be some zero paths, and then allPaths and everyLength would disagree
        for path in allPaths:
            length = 0
            seq = []
            repeats = []
            maxPath = self.maxPathLen(path)
            minWeight = self.cap(path)  # minweight of the path
            if minWeight > 0:   # makes sure we aren't adding subtracted paths
                everyLength.append(maxPath)
                everyPath.append(path)
        return everyPath, everyLength

    def cap(self, path):  # finds min capacity of a path
            repeats = []  # counts repeated edges to take loops into account
            p1 = self.ports[path[0]]
            p2 = self.ports[path[1]]
            pathMin = p1.findEdge(p2).label  # smallest weighted edge on a particular path (we initialize it to be a random weight on the path)
            for port in range(len(path)-1):   # loop through edges of path and find smallest one
                parentp = path[port]   # find the ports
                childp = path[port+1]
                port1 = self.ports[parentp]
                port2 = self.ports[childp]
                repeats.append(port1.node + port2.node)  # put two copies of the adjacent nodes, so that future methods will find it no matter the order.
                repeats.append(port2.node + port1.node)

                parentNode = self.nodes[port1.node]
                if port1.node != "Portal":
                    port1 = parentNode.jump(port1)  # if port1 is on port start, jump to end (we want to make sure it is adjacent to next port).
                if port1.node != port2.node and repeats.count(port1.node+port2.node) > 1:
                    label = port1.findEdge(port2).label
                    occurences = repeats.count(port1.node+port2.node)
                    weight = label/occurences
                else:
                    weight = port1.findEdge(port2).label
                if weight < pathMin:
                    pathMin = weight
            return pathMin

    def sub(self, minWeight, path, status, i):
        subtracted = []
        for port in range(len(path)-1):
                parentp = path[port]   # find the ports in path
                childp = path[port+1]
                port1 = self.ports[parentp]
                port2 = self.ports[childp]

                parentNode = self.nodes[port1.node]
                if port1.node != "Portal":
                    port1 = parentNode.jump(port1)  # if port1 is on port 'start', we want to make sure it is adjacent to next node.
                if port1.findEdge(port2).label != 0:  # if edge is not weight 0
                    if port1 != port2:
                        port1.minus(port2, minWeight)  # subtract minweight using minus()
                    else:
                        port1.findEdge(port2).label -= minWeight   # this is for loops, otherwise it will subtract the minWeight from the same edge twice
                else:
                    status = "Done"

    def findBestPaths(self):  # finds thinnest point on chosen (longest) path (using cap()), and subtracts it (using sub()), printing out longest paths as it goes.
        status = "Continue"
        counter = 0
        smallest = 0
        bestPaths = []
        everyPath, everyLength = self.allMaxPathLengths()
        while status == "Continue":
            paths_so_far = []
            lens_so_far = []
            if len(everyLength) == 0:
                status = "Done"
            else:
                maxCap = max(everyLength)
                smallest = min(everyLength)
                ind = [i for i, j in enumerate(everyLength) if j == maxCap] # returns a list of the indices where the maxes are
                paths_so_far += [everyPath[ind[0]]]  # pick the first index
                lens_so_far += [everyLength[ind[0]]]
                length = everyLength[ind[0]]
                path = everyPath[ind[0]]
                minWeight = self.cap(path) # run cap() to find thickness of path
                bestPaths.append((length, path, minWeight))  # append longest path to the list of best paths
                self.sub(minWeight, path, status, ind) # subtract that path from the graph
                everyPath, everyLength = self.allMaxPathLengths() # find paths and lengths from modified graph
        return bestPaths


    def draw(self, filename, bestPaths):    # highlights findBestPaths() paths
        legend = {}
        colors = ["darkorchid", "deepskyblue", "springgreen3", "orchid1", "gold", "darkturquoise", "maroon1"]

        print "\n"
        print "Maximum Capacity Paths:\n"
        print "|    Path #    |        Color        |       Min Edge Capacity       |       Path        |\n"
        i = 0
        for path in bestPaths:  # pick maxCapacity paths
            pathStr = ""
            col = colors.pop(random.randrange(len(colors))) # pick a random color
            #print "Path:", path
            label = path[2]
            # print label
            # print label, path[1]
            mapping = path[1]
            for node in range (len(mapping)):
                if "Portal" not in mapping[node]:
                    if node < (len(mapping) - 2):
                        pathStr += mapping[node] + " -> "
                    else:
                        pathStr += mapping[node]

            
            print "|    " + str(i) + "    |        " + col + "       |       " + str(label) + "       |       " + pathStr + "    |\n"
            
            legend[col] = label, path[1]   # add data to legend
            i += 1
        self.makeColorCode(legend, filename) # pass data to makeColorCode

    def makeColorCode(self, legend, filename):
        newGraph = Graph()
        nodesToAdd = []
        for entry in legend:
            color = entry
            weight, path = legend[entry]
            for node in range(len(path)-1):
                node1 = path[node].split("$")[0]
                node2 = path[node+1].split("$")[0]
                port1 = path[node].split("$")[1]
                port2 = path[node+1].split("$")[1]

                if node1 != "Portal" and port1 == "start":
                    port1 = "end"
                elif node1 != "Portal" and port1 == "end":
                    port1 = "start"
                port1 = Port(port1)
                node1 = Node(node1)
                node1.add(port1)
                
                port2 = Port(port2)
                node2 = Node(node2)
                node2.add(port2)
                port1.add(port2, weight, color)

                if node1.n in self.nodes:
                    x = self.nodes[node1.n].atts["x"]
                    y = self.nodes[node1.n].atts["y"]
                    node1.atts["x"] = x
                    node1.atts["y"] = y

                if node2.n in self.nodes:
                    x = self.nodes[node2.n].atts["x"]
                    y = self.nodes[node2.n].atts["y"]
                    node2.atts["x"] = x
                    node2.atts["y"] = y
                    
                nodesToAdd.append(node1)                # and add node to graph
                nodesToAdd.append(node2)

        i = 0
        ids = {}
        for node in nodesToAdd:    # assign id to node and then add it to graph
            if node.n in ids:
                node.changeID(ids[node.n])
                newGraph.add(node)
            else:
                ids[node.n] = i
                node.changeID(i)
                newGraph.add(node)
            i += 1
        newGraph.makeGraphBlueprint(filename)


# ''''''''''''''''''''RUN SCRIPT''''''''''''''''''''''''

print "-------------------DNA_MUTATION_MAPPING-------------------"
filename = sys.argv[1]

template = "graph_template.dot"

# create graph from .dot data
graph = Graph()
graph.readIn(template)

# create graph
thePaths = graph.findBestPaths()
graph.draw(filename, thePaths)


# terminal commands to convert .dot file into .png file
import subprocess

text_to_run = '/usr/local/bin/neato -n2 -Tpng %s.dot > %s.png' % (filename, filename)
subprocess.call([text_to_run],shell=True)

# terminal commands to open .png
text_to_run = 'open *.png'
subprocess.call([text_to_run],shell=True)
