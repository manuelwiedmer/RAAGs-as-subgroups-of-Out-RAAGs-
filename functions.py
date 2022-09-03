import networkx as nx
import itertools as it
import copy
import matplotlib.pyplot as plt
import random
import sys
import math
import pickle
from itertools import combinations
from itertools import permutations
import os

def star(vertex,Graph):
    """Returns the star of vertex in Graph.
    
    Graph : a networkx graph
    vertex : a vertex of Graph
    """
    
    st = set(Graph[vertex])|{vertex}
    return st

def conn_comps(Gamma,v):
    """Lists the connected components of Gamma\st(v).
    
    Gamma : a networkx graph
    v : a vertex of Gamma
    """
    
    Gamma_without_st = copy.deepcopy(Gamma)
    Gamma_without_st.remove_nodes_from(star(v,Gamma))
    return list(nx.connected_components(Gamma_without_st))

def is_Sil(Gamma,a,b):
    """Checks whether a and b form a SIL pair in the graph Gamma.
    
    Gamma : a networkx graph
    a : a vertex of Gamma
    b : a vertex of Gamma
    """

    for x in Gamma.neighbors(a):
        # first check whether they are adjacent
        if x==b:
            return False
    for comp1 in conn_comps(Gamma,a):
        for comp2 in conn_comps(Gamma,b):
            if comp1 == comp2:
                return True
    return False

def has_cycle(G):
    """Returns True if G has a cycle
    
    G : a networkx graph
    """
    if nx.cycle_basis(G):
        return(True)
    return(False)

def supp_graph(Gamma,vertex):
    """Returns the support graph of vertex.
    
    Gamma : a networkx graph
    vertex : a vertex of Gamma
    """

    suppgraph = nx.Graph()
    node = 0

    for comp in conn_comps(Gamma, vertex):
        suppgraph.add_node(node, vertices = comp)
        node = node + 1

    for n in suppgraph.nodes():
        for m in suppgraph.nodes():
            if n != m:
                for x in suppgraph.nodes[m]["vertices"]:
                    for comp in conn_comps(Gamma,x):
                        if suppgraph.nodes[n]["vertices"]==comp:
                            suppgraph.add_edge(n,m)
    
    return(suppgraph)

def supp_graph_number_of_vertices(Gamma,vertex):
    """Returns the number of vertices of the support graph of vertex.
    
    Gamma : a networkx graph
    vertex : a vertex of Gamma
    """    
    return(len(Gamma.nodes())-len(star(vertex,Gamma)))

def theta(Gamma):
    """Returns theta, the defining graph for the RAAG structure of PSO(A_Gamma) 
    
    Gamma : a networkx graph
    """

    theta = nx.Graph()
    node = 0
    for v in Gamma.nodes:
        sg = supp_graph(Gamma,v)
        if len(sg.edges) > 0:
            for supp_edge in sg.edges:
                supp_edge_start = supp_edge[0]
                supp_edge_end = supp_edge[1]
                theta.add_node(node,comp = [v,sg.nodes[supp_edge_start]["vertices"],sg.nodes[supp_edge_end]["vertices"]])
                node = node + 1
    for x in theta.nodes:
        for y in theta.nodes:
            added_edge = False
            if x !=y:
                if not is_Sil(Gamma,theta.nodes[x]["comp"][0],theta.nodes[y]["comp"][0]):
                    theta.add_edge(x,y)
                else:
                    if theta.nodes[x]["comp"][0] in theta.nodes[y]["comp"][1]:
                        non_dominant_comp1 = theta.nodes[y]["comp"][2]
                    elif theta.nodes[x]["comp"][0] in theta.nodes[y]["comp"][2]:
                        non_dominant_comp1 = theta.nodes[y]["comp"][1]
                    else:
                        theta.add_edge(x,y)
                        added_edge = True
                    if not added_edge:
                        if theta.nodes[y]["comp"][0] in theta.nodes[x]["comp"][1]:
                            non_dominant_comp2 = theta.nodes[x]["comp"][2]
                        elif theta.nodes[y]["comp"][0] in theta.nodes[x]["comp"][2]:
                            non_dominant_comp2 = theta.nodes[x]["comp"][1]
                        else:
                            theta.add_edge(x,y)
                            added_edge = True
                    if not added_edge:
                        #If the code arrives here, the pairs contain dominant components for one another.
                        if non_dominant_comp1 != non_dominant_comp2:
                            #other components do not coincide:
                            theta.add_edge(x,y)
                            
    #Adds also the vertices v_C^a
    for v in Gamma.nodes():
        i = nx.number_connected_components(supp_graph(Gamma,v))-1;
        for j in range(i):
            theta.add_node(node)
            for w in theta.nodes():
                if w != node:
                    theta.add_edge(w,node)
            node = node + 1
            
    return theta

def supp_graphs_forest(Gamma):
    """Check for a given graph Gamma, whether PSO(A_\Gamma) is a RAAG. Using a criterion of Day-Wade, this is equivalent to showing that all support graphs are forests. Returns True if all support graphs are forest, False otherwise.
    
    Gamma : a networkx graph
    """
    
    for v in Gamma.nodes():
        if has_cycle(supp_graph(Gamma,v))==True:
            return(False)
    return(True)

def theta_number_of_vertices(Gamma):
    """Only works if supp_graphs_forest(Gamma)==True, then it returns the number of edges of the corresponding theta graph
    
    Gamma : a networkx graph
    """
    
    n = 0
    for v in Gamma.nodes():
        A = supp_graph(Gamma,v)
        n += len(A.edges())
        n += max(nx.number_connected_components(A)-1,0)
    
    return(n)

def check_finite_index(G):
    """Checks whether PSO(A_Gamma) has finite index in Out(A_Gamma)
    
    Gamma : a networkx graph
    """
    for (u,v) in permutations(G.nodes(),2):
        if(u != v):
            A = star(u,G)
            A.remove(u)
            B = star(v,G)
            if(A.issubset(B)):
                return False
    return True

def progressBar(name, value, endvalue, bar_length = 25, width = 20):
    """Displays a progress bar with "name", tracking the progess of value/endvalue.
    """
    
    percent = float(value) / endvalue
    arrow = '-' * int(round(percent*bar_length) - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write("\r{0: <{1}} : [{2}]{3}%".format(name, width, arrow + spaces, int(round(percent*100))))
    sys.stdout.flush()
    if value == endvalue:
         sys.stdout.write('\n\n')
            
def allgraphs(i, n):
    """Creates the ith graph with n vertices (we use the binary representation of i to create a adjacency matrix)
    """
    G = nx.Graph()
    G.add_nodes_from(range(n))
    j = 0
    k = 1
    while(i != 0):
        if (i%2 != 0):
            G.add_edge(j,k)
            i -= 1
        i /= 2
        if (k != n-1):
            k += 1
        else:
            j += 1
            k = j+1
    return(G)

def progressBar_withnumber(name, value, endvalue, number, bar_length = 25, width = 20):
    """
    Displays a progress bar with "name", tracking the progess of value/endvalue.
    Also displayes the number of graphs found.
    """
    
    percent = float(value) / endvalue
    arrow = '-' * int(round(percent*bar_length) - 1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    sys.stdout.write("\r{0: <{1}} : [{2}]{3}%".format(name, width, arrow + spaces, int(round(percent*100)))+", found "+str(number)+" graphs")
    sys.stdout.flush()
    if value == endvalue:
         sys.stdout.write('\n\n')
            
def remove_duplicates(name_of_input_file):
    """Remove duplicate graphs from name_of_input_file
    """

    with open(name_of_input_file, "rb") as fp:
        graphs = pickle.load(fp)
    
    i_max = math.comb(len(graphs),2)
    i = 0
    for A in combinations(graphs,2):
        progressBar("Iterations",i,i_max)
        i += 1
        if(nx.is_isomorphic(A[0][0],A[1][0])):
            if A[1] in graphs:
                graphs.remove(A[1])

    progressBar("Iterations",i_max,i_max)
    
    if os.path.isfile(name_of_input_file):    
        os.remove(name_of_input_file) 
        
    with open(name_of_input_file, "wb") as fp:
        pickle.dump(graphs, fp)
        
def print_file (name_of_input_file, name_of_output_file):
    """Print the graphs in name_of_input_file to name_of_output_file
    """

    if os.path.isfile(name_of_output_file):
        os.remove(name_of_output_file)

    with open(name_of_input_file, "rb") as fp:
        graphs = pickle.load(fp)

    out=open(name_of_output_file,'a')


    if len(graphs)!=1:
        out.write(str(len(graphs))+" graphs have been found"+'\n'+'\n')
    else:
        out.write(str(len(graphs))+" graph has been found"+'\n'+'\n')

    for G in graphs:
        out.write("Number of vertices of Gamma: "+str(len(G[0].nodes))+'\n')
        out.write("Number of edges of Gamma: "+str(len(G[0].edges))+'\n')
        out.write("Nodes of Gamma: " + str(G[0].nodes)+'\n')
        out.write("Edges of Gamma: " + str(G[0].edges)+'\n')
        out.write("Number of vertices of Theta: "+str(len(G[1].nodes))+'\n')
        out.write("Number of edges of Theta: "+str(len(G[1].edges))+'\n')
        out.write("Nodes of Theta: " + str(G[1].nodes)+'\n')
        out.write("Edges of Theta: " + str(G[1].edges)+'\n')
        out.write('\n')
    out.close()
    
def extract_small_theta_graphs(name_of_input_file, name_of_output_file, n):
    """Extract of name_of_input_file all Theta graphs of size n and save them in name_of_output_file
    We save them with one of the smallest found Gamma graph (with respect to the number of vertices).
    """
    
    if os.path.isfile(name_of_output_file):
        os.remove(name_of_output_file)

    with open(name_of_input_file, "rb") as fp:
        graphs = pickle.load(fp)

    graphs_new = set()

    for A in graphs:
        if len(A[1].nodes()) == n:
            already_there = False
            graphs_to_remove = set()
            for B in graphs_new:
                if nx.is_isomorphic(A[1],B[1]):
                    if len(A[0]) < len(B[0]):
                        graphs_to_remove.add(B)
                    else:
                        already_there = True
            if already_there == False:
                graphs_new.add(A)
            for C in graphs_to_remove:
                if C in graphs_new:
                    graphs_new.remove(C)

    if os.path.isfile(name_of_output_file):
        os.remove(name_of_output_file)                

    with open(name_of_output_file, "wb") as fp:
        pickle.dump(graphs_new, fp)