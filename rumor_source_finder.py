import matplotlib
matplotlib.use('Agg')
import random
import networkx
import matplotlib.pyplot as plt
from math import *
from tqdm import tqdm

class rumorFinder:
    # the function takes a filename as input. This file has network topology
    # e.g. node 1 3 4 means node 1 is connected to node 3 and node 4.
    # numNodes represents the number of nodes to be used in making the graph.
    # The function returns an adjacency list of the graph, which is a list of list. 
    def buildAdjacencyList(self, filename, numNodes):
        try:
            fd = open(filename, 'r')
        except:
            print("Error in opening file -", filename)
            exit()
        lines = fd.readlines()
        adjacency = [[] for _ in range(0, numNodes)]
        
        # reading each line from given data file
        # and creating adjacency list file's content are in form:
        # <source> <destination> <edge direction> <edge weight>
        for line in lines:
            lineSplit = line.split()
            
            # vertex name 1 in file is named 0 in adj list
            x, y = int(lineSplit[0])-1, int(lineSplit[1])-1
            if x >= numNodes:
                break
            elif y >= numNodes:
                continue
            else:
                adjacency[x].append(y)
                adjacency[y].append(x)
        fd.close()
        
        # cleaning graph by removing nodes with degree less than 3
        minDegree = 3
        deleteNodes = []

        for i in range(0, len(adjacency)):
            if len(adjacency[i]) < minDegree:
                deleteNodes.append(i)
                adjacency[i] = []

        for i in range(0, len(adjacency)):
            for delNode in deleteNodes:
                if delNode in adjacency[i]:
                    adjacency[i].remove(delNode)
        return adjacency
    
    # randomly selecting a source node from the graph to start rumor
    def selectSource(self,adjacency):
        src = random.randint(0,len(adjacency)-1)
        while len(adjacency[src]) == 0:
            src = random.randint(0,len(adjacency)-1)
        return src
    
    # This function simulates the rumor spreading in the graph
    # It takes graph and number of nodes to infect as input
    # then simulates the spreading of rumor from the source   
    # using SI model and returns the infected nodes and their 
    # corresponding nodes that infected them. 
    def rumorSpreading(self,src,adjacency,numNodesToInfect):
        infectionSource = [[] for _ in range(0,numNodesToInfect)]
        infected = [-1]*numNodesToInfect
        infected[0]=src

        i = 1
        while i < numNodesToInfect:
            susceptible = list()
            srcIdx = list()
            
            # identifying the susceptible nodes next to an
            # infected node
            for j in range(0,i):
                neighbour = adjacency[infected[j]]
                for k in range(0,len(neighbour)):
                    if neighbour[k] not in infected:
                        susceptible.append(neighbour[k])
                        srcIdx.append(j)

            # randomly selecting a susceptible node to infect
            randIdx = random.randint(0,len(susceptible)-1)
            infectedNew = susceptible[randIdx]
            src = srcIdx[randIdx]
            infected[i]=infectedNew

            # updating the source of infected node
            infectionSource[src].append(i)
            infectionSource[i]=[src]
            i += 1
        return infectionSource,infected

    # Function to calculate the size of a subtree and the cumulative product
    # of the size of the subtrees of all nodes in its subtree. This is 
    # calculated in an upward fashion. 
    def upwardMsgs(self, upmsg, infectionSource, parent, current):
        if len(infectionSource[current]) == 1:
            upmsg[parent][0] += 1
            upmsg[parent][1] *= upmsg[current][1]
        elif parent == current:
            for i in range(0, len(infectionSource[current])):
                # here all nodes in infectionSource[current] are considered as
                # children of current node
                upmsg = self.upwardMsgs(upmsg, infectionSource, current, infectionSource[current][i])
        else:
            for i in range(0, len(infectionSource[current])):
                childNode = infectionSource[current][i]
                if parent != childNode:
                    upmsg = self.upwardMsgs(upmsg, infectionSource, current, childNode)
                    upmsg[current][1] *= upmsg[childNode][1]
            upmsg[current][1] *= upmsg[current][0]
            upmsg[parent][0] += upmsg[current][0]
        return upmsg

    # Function to calculate the Rumor Centrality for each node using the message
    # values generated in upwardMsgs
    def downwardMsgs(self, msgDown, msgUp, infectionSource, parent, current):
        # number of nodes in the rumor graph
        N = len(infectionSource)

        # calculating Rumor Centrality when the node passed as parameter is not root
        if current != parent:
            msgDown[current] = msgDown[parent] + log(msgUp[current][0]) - log(N - msgUp[current][0])
            for i in range(0,len(infectionSource[current])):
                if infectionSource[current][i] != parent:
                    msgDown = self.downwardMsgs(msgDown, msgUp, infectionSource, current, infectionSource[current][i])

        # calculating Rumor Centrality when node passed as parameter is root
        else:
            msgDown[current] = log(factorial(N)) - log(N)
            for i in range(0,len(infectionSource[current])):
                msgDown[current] -= log(msgUp[infectionSource[current][i]][1])
            for j in range(0,len(infectionSource[current])):
                msgDown = self.downwardMsgs(msgDown, msgUp, infectionSource, current, infectionSource[current][j])
        return msgDown

    # This function finds the rumor source in the network by 
    # using rumor centrality. First the size of a subtree at a node
    # and the cumulative product of the size of the subtrees of all
    # nodes in its subtree is calculated by up messages. Then the rumor
    # centrality is calculated in downward fashion. The node with maximum
    # rumor centrality is then returned as the expected rumor centre.
    def findRumourCentre(self, infectionSource):
        root, rumorCenter = 0, -1
        # N is number of infected nodes
        N = len(infectionSource)
        upMsgs = [[1,1] for _ in range(0, N)]
        downMsgs = [1] * N                                                                                               
        upMsgs = self.upwardMsgs(upMsgs, infectionSource, root, root)
        downMsgs = self.downwardMsgs(downMsgs, upMsgs, infectionSource, root, root)
        rumorCenter = downMsgs.index(max(downMsgs))
        return rumorCenter
    
    # Function for plotting the graph
    def plotGraph(self, NMin, NMax, rumorError, NValues):
        plt.plot(NValues, [i*100 for i in rumorError], lw = 2.5, label = "Error in Detection")
        plt.axis([NMin, NMax, min(rumorError), 100])
        plt.xlabel("Number of Nodes Infected")
        plt.ylabel("Error in Detecting Source")
        plt.legend(loc = "lower right")
        plt.title("Error in Source Detection")
        plt.grid(which = "major")
        plt.savefig("graph.png")

def main():
    rfObj = rumorFinder()
    
    # minimum number of nodes to infect
    minNodes = 20

    # maximum number of nodes to infect
    maxNodes = 200

    # x-axis used for plotting graph
    xAxis = list()

    # list to store success rate of the rumor source prediction
    rumorRate = list()
    maxIterations = 50

    adjacency = rfObj.buildAdjacencyList("fb_out.txt", 5000)

    # taking readings of the rumor centrality by varying the number
    # number of nodes infected in the graph.
    for i in tqdm(range(minNodes, maxNodes+1, 10)):
        numOferrors = 0
        counter, j = 0, 0
        
        while j < maxIterations:
            while 1:
                try:
                    src = rfObj.selectSource(adjacency)
                    infectionSource, infected = rfObj.rumorSpreading(src, adjacency, i)
                    break
                except:
                    continue
            rumorCenter = infected[rfObj.findRumourCentre(infectionSource)]

            # checking if the expected rumor source is equal to
            # the actual rumor source.
            if src != rumorCenter:
                numOferrors += 1
            counter += 1
            j += 1
        
        # saving the error rate for results
        if counter != 0:
            rumorRate.append(numOferrors/counter)
            xAxis.append(i)
        
    # plotting results on graph
    rfObj.plotGraph(minNodes, maxNodes, rumorRate, xAxis)

if __name__ == "__main__":
    main()