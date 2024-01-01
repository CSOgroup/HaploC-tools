#https://codereview.stackexchange.com/questions/159381/mst-kruskals-algorithm-in-python
# Python program for Kruskal's algorithm to find Minimum Spanning Tree
import os,sys
#from collections import defaultdict

#Class to represent a graph
class Graph:
	def __init__(self,vertices):
        	self.V= vertices #No. of vertices
	        self.graph = [] # default dictionary to store graph
		self.edges = 0;

	# function to add an edge to graph
	def addEdge(self,u,v,w):
        	self.graph.append([u,v,w]); self.edges +=1;
	
	# function to add an edge to graph
	def addEdgeList(self,edgelist):
		for edge in edgelist:  self.graph.append(edge); self.edges +=1; 

    # A utility function to find set of an element i, (uses path compression technique)
	def find(self, parent, i):
		if parent[i] == i: return i
		return self.find(parent, parent[i])

    # A function that does union of two sets of x and y,  (uses union by rank)
	def union(self, parent, rank, x, y):
		xroot = self.find(parent, x);
		yroot = self.find(parent, y)

		# Attach smaller rank tree under root of high rank tree
		if rank[xroot] < rank[yroot]:   parent[xroot] = yroot
		elif rank[xroot] > rank[yroot]: parent[yroot] = xroot
		#If ranks are same, then make one as root and increment its rank by one
		else : parent[yroot] = xroot; rank[xroot] += 1

    # The main function to construct MST using Kruskal's algorithm
	def KruskalMST(self):
		result =[] #This will store the resultant MST
		i =0; e = 0 # An index variable

	        #Step 1:  Sort all the edges in increasing order of their # weight
		self.graph =  sorted(self.graph,key=lambda item: -1*item[2])

		parent = [] ; rank = []
	        # Create V subsets with single elements
		for node in xrange(self.V): parent.append(node); rank.append(0); 

	        # Number of edges to be taken is equal to V-1
		while e < self.V -1 and i < self.edges :
			# Step 2: Pick the smallest edge and increment the index
			# for next iteration
			u=  self.graph[i][0]; v = self.graph[i][1]; w = self.graph[i][2];
			x = self.find(parent, u)
			y = self.find(parent ,v)
			#print self.graph[i],x,y

			# If including this edge does't cause cycle, include and increment the index of result for next edge
			if x != y:
				e = e + 1  
				result.append(self.graph[i])
				self.union(parent, rank, x, y)          

			i = i + 1

		# print the contents of result[] to display the built MST
		#print "Following are the edges in the constructed MST"
		#for edge  in result:
			#print str(u) + " -- " + str(v) + " == " + str(weight)
			#if len(edge) ==4: print ("%d -- %d == %f %s" % (edge[0],edge[1],edge[2],edge[3]))
			#else: print ("%d -- %d == %f" % (edge[0],edge[1],edge[2]))
		return result;


def test_algo():
	g = Graph(15)
	g.addEdge(0, 1, 7)
	g.addEdge(0, 3, 3)
	g.addEdge(1, 2, 3)
	g.addEdge(1, 4, 2)
	g.addEdge(2, 5, 2)
	g.addEdge(3, 4, 3)
	g.addEdge(3, 6, 2)
	g.addEdge(4, 5, 5)
	g.addEdge(4, 7, 7)
	g.addEdge(5, 8, 3)
	g.addEdge(6, 7, 3)
	g.addEdge(6, 9, 1)
	g.addEdge(7, 8, 7)
	g.addEdge(7, 10, 3)
	g.addEdge(9, 12, 4)
	g.addEdge(9, 10, 1)
	g.addEdge(10, 11, 6)
	g.addEdge(10, 13, 7)
	g.addEdge(11, 12, 6)
	g.addEdge(11, 14, 2)
	g.addEdge(12, 13, 4)
	g.addEdge(13, 14, 5)
	g.KruskalMST()

#test_algo()
