#!/usr/bin/env python
from __future__ import division, absolute_import, print_function
import numpy as np
import networkx as nx
import random
import copy


__all__ = ['create_network', 'source_nodes', 'sink_nodes', 'plot_network']


def create_network(nnodes, nedges, maxdegree=None, maxoutdegree=None,
                   maxindegree=None, reverse=False, verbose=False):
    """
    This functions creates a directed acyclic graph with one source node, N nodes and M edges in total.
    Each node of the network is connected to the source node.
    Optionally the an upper degree of the nodes can be given to generate networks with are not too dense.
    The only exception from this upper bound is the source node.

    Definition
    ----------
    def create_network(nnodes, nedges, maxdegree=None)


    Input           Format                  Description
    -----           -----                   -----------
    nnodes          integer                 Number of nodes of the network. One node is the source not (no incoming edges, only outgoing).
                                            Multiple nodes might be sinks (no outgoing edges, only incoming).
    nedges          integer                 Number of (directed) edges of the network.
                                            Upper limit is ((maxdegree+1)*nnodes -  maxdegree - 1 ) / 2
                                               My calculations for that are:
                                               - source node has at most (nnodes-1) edges
                                               - all other (nnodes-1) nodes have at most maxdegree edges
                                               - all edges are now count twice because of directed property
                                            Lower limit is nnodes - 1
                                               Due to the fact that source node need to be connected at least to every other node
                                            TODO
                                               Redo this calculations when maxoutdegree and/or maxindegree are given
    maxdegree        integer                Upper degree of a node to otain not too dense networks. The only node is allowed to
                                            have a higher degree is the source node.
    maxoutdegree     integer                Upper limit of outgoing edges per node (except source node)
    maxindegree      integer                Upper limit of incoming edges per node (except source node)
    reverse          bool                   False: source node is the only source
                                            True:  after generatng the whole network (with one source node and maxdegree settings)
                                                   all edges are swapped such that one obtain a network with only one sink and multiple sources
                                                   Note: degree settings are made for one-source networks, i.e maxoutdegree -> maxindegree
                                                                                                               maxindegree  -> maxoutdegree
    verbose          bool                   If true print-outs, if false silent


    Description
    -----------
    (1) Adjacency matrix A[i,j]:
        - if A[i,j]=1 indicates that node i is connected with node j (i --> j)
        - if A[i,j]=0 indicates that connection i --> j does not exist
        A[i,j] != A[j,i] for directed graphs
    (2) Initialize A = 0
    (3) Fill only upper triangle of A (i < j) with nedges ones
        (this assures that graph is acyclic)
    (4) Diagonal needs to be zero such that there are no self-loops
    (5) Test if degree of all nodes (except first node = source node) have degree less eqaul max. degree given
    (6) Test if first node (source node) is connceted to all other nodes
        (assure that all nodes are connected and no subgraph exist)

    Restrictions
    ------------
    None

    Examples
    --------
    >>> random.seed('12345')
    >>> nnodes = 4
    >>> nedges = 5
    >>> G = create_network(nnodes, nedges, maxdegree=3)
    >>> nodes = G.nodes()
    >>> print('nodes: '+str(nodes))
    nodes: [0, 1, 2, 3]
    >>> edges = G.edges()
    >>> print('edges: '+str(edges))
    edges: [(0, 1), (0, 3), (1, 2), (1, 3), (2, 3)]
    >>> degrees = [ G.degree(ii) for ii in G.nodes() ]
    >>> print('degrees of nodes: '+str(degrees))
    degrees of nodes: [2, 3, 2, 3]

    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License.

    Copyright (c) 2016-2021 Juliane Mai

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.


    History
    -------
    Written,  Juliane Mai,    Oct 2016
    Modified, Matthias Cuntz, Sep 2021 - code refactoring
    """

    if maxdegree is None:
        maxdegree = nnodes - 1

    if maxoutdegree is None:
        maxoutdegree = nnodes - 1

    if maxindegree is None:
        maxindegree = nnodes - 1

    if (nedges < nnodes - 1):
        raise ValueError('Increase number of edges because no connected network is possible')

    if (nedges > int(((maxdegree+1)*nnodes - maxdegree - 1 ) / 2.0) ):
        raise ValueError('Decrease number of edges because at least one of the nodes will exceed upper limit of degrees')

    npos_edges = nnodes * (nnodes-1) / 2
    if (nedges > npos_edges):
        raise ValueError('Too many edges for given number of nodes! Upper limit is '+str(npos_edges))

    found_network = False
    while (not found_network):

        adjacency = np.zeros([nnodes, nnodes], dtype=int)
        cedge = 0
        while cedge < nedges:
            # integer random number of interval [0,npos_edges)
            rand = random.randrange(npos_edges)
            csum = np.cumsum(range(nnodes-1, 0, -1))

            # determine col and row for this number
            row = np.where(rand < csum)[0][0]
            if row == 0:
                sublist = np.arange(0, csum[row])
            else:
                sublist = np.arange(csum[row-1], csum[row])
            col = (row+1) + np.where(rand == sublist)[0][0]

            # check if this element is already set to "1",
            # if yes try another one
            if (adjacency[row, col] == 0):
                cedge += 1
                adjacency[row, col] = 1

        # calculate degrees of each node
        degrees     = [ np.sum(adjacency[:, ii]) + np.sum(adjacency[ii, :])
                        for ii in range(nnodes) ]
        degrees_in  = [ np.sum(adjacency[:, ii])
                        for ii in range(nnodes) ]
        degrees_out = [ np.sum(adjacency[ii, :])
                        for ii in range(nnodes) ]

        # print('degrees:     ',degrees)
        # print('degrees_in:  ',degrees_in)
        # print('degrees_out: ',degrees_out)

        # check if each node of network has degree less equal maxdegree
        degrees_satisified = ( all([ degrees[ii]     <= maxdegree
                                     for ii in range(1, nnodes) ]) and
                               all([ degrees_in[ii]  <= maxindegree
                                     for ii in range(1, nnodes) ]) and
                               all([ degrees_out[ii] <= maxoutdegree
                                     for ii in range(1, nnodes) ])
                               )

        if (degrees_satisified):
            if reverse:
                # swap all edges
                # only sink is last node
                adjacency = rotate(rotate(np.transpose(adjacency)))

            if verbose:
                print('degrees:     ', degrees)

            # create graph
            G = nx.MultiDiGraph()

            # add nodes to graph, IDs are 0...(n-1)
            G.add_nodes_from([ii for ii in range(nnodes)])
            nodes = G.nodes()
            if verbose:
                print('Added nodes: ', nodes)

            # add edges
            for ii in range(nnodes):
                for jj in range(ii+1, nnodes):
                    if adjacency[ii, jj] == 1:
                        G.add_edge(ii, jj)
            edges = G.edges()
            if verbose:
                print('Added edges: ', edges)

            if reverse:
                # check if all are connceted to sink node "nnodes-1"
                all_connected_to_0 = all([nx.has_path(G, ii, nnodes-1)
                                          for ii in range(1, nnodes)])
            else:
                # check if all are connceted to source node "0"
                all_connected_to_0 = all([nx.has_path(G, 0, ii)
                                          for ii in range(1, nnodes)])

            if all_connected_to_0:
                found_network = True
                # finished
            else:
                G.clear()
                if verbose:
                    print('New search required')
                    print('')
                # start new search
    return G


def source_nodes(G):
    """
    This functions returns all source nodes (no incoming edges) of a graph G.

    Definition
    ----------
    def source_nodes(G)


    Input           Format                  Description
    -----           -----                   -----------
    G               graph                   Networkx generated graph structure

    Description
    -----------
    None

    Restrictions
    ------------
    None

    Examples
    --------
    >>> random.seed('12345')
    >>> nnodes = 4
    >>> nedges = 5
    >>> G = create_network(nnodes, nedges, maxdegree=3)
    >>> sources = source_nodes(G)
    >>> print('source nodes: '+str(sources))
    source nodes: [0]

    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License.

    Copyright (c) 2016-2021 Juliane Mai

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.


    History
    -------
    Written,  Juliane Mai,    Oct 2016
    Modified, Matthias Cuntz, Sep 2021 - code refactoring
    """
    idx     = np.where( 0 == np.array([ G.in_degree(ii)
                                        for ii in G.nodes() ]))[0]
    sources = np.array(G.nodes())    # all nodes
    sources = sources[idx]           # only sources nodes
    return sources


def sink_nodes(G):
    """
    This functions returns all sink nodes (no outgoing edges) of a graph G.

    Definition
    ----------
    def sink_nodes(G)


    Input           Format                  Description
    -----           -----                   -----------
    G               graph                   Networkx generated graph structure

    Description
    -----------
    None

    Restrictions
    ------------
    None

    Examples
    --------
    >>> random.seed('12345')
    >>> nnodes = 4
    >>> nedges = 5
    >>> G = create_network(nnodes, nedges, maxdegree=3)
    >>> sinks = sink_nodes(G)
    >>> print('sink nodes: '+str(sinks))
    sink nodes: [3]

    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License.

    Copyright (c) 2016-2021 Juliane Mai

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.


    History
    -------
    Written,  Juliane Mai,    Oct 2016
    Modified, Matthias Cuntz, Sep 2021 - code refactoring
    """
    idx     = np.where( 0 == np.array([ G.out_degree(ii)
                                        for ii in G.nodes() ]))[0]
    sinks   = np.array(G.nodes())      # all nodes
    sinks   = sinks[idx]               # only sink nodes
    return sinks


def plot_network(G, fname_plot):
    """
    This functions plots a graph using either graphviz if available or standard
    matplotlib. Plots with graphviz are much more fancy since it reduces the
    number of intersections of plotted edges and draws round edges unstead of
    only straight lines.

    Definition
    ----------
    def plot_network(G, fname_plot)


    Input           Format                  Description
    -----           -----                   -----------
    G               graph                   Networkx generated graph structure
    fname_plot      string                  Filename of ouput pdf

    Description
    -----------
    None

    Restrictions
    ------------
    None

    Examples
    --------
    >>> random.seed('12345')
    >>> nnodes = 4
    >>> nedges = 5

    >>> G = create_network(nnodes, nedges, maxdegree=3)
    >>> fname_plot = 'test_dag.pdf'
    >>> plot_network(G, fname_plot)

    >>> import os
    >>> os.remove(fname_plot)

    License
    -------
    This file is part of the JAMS Python package, distributed under the MIT
    License.

    Copyright (c) 2016-2021 Juliane Mai

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.


    History
    -------
    Written,  Juliane Mai,    Oct 2016
    Modified, Matthias Cuntz, Sep 2021 - code refactoring
    """
    import matplotlib.pyplot as plt

    color   = [(0.237736, 0.340215, 0.575113),
               (0.404682, 0.548016, 0.244594),
               (0.786572, 0.753215, 0.297951)]
    sources = source_nodes(G)
    sinks   = sink_nodes(G)
    nodes = G.nodes()
    edges = G.edges()

    try:
        import pygraphviz as pgv

        A = pgv.AGraph(directed=True)
        A.add_edges_from(G.edges(),
                         color='grey', arrowhead='vee')

        alpha_fill  = 0.7
        alpha_bndr  = 0.0

        for ii in nodes:

            # determine color of node
            if ii in sinks:
                cc          = color[2]
                cc_hex_fill = "#%2x%2x%2x%2x" % (int(cc[0]*255),
                                                 int(cc[1]*255),
                                                 int(cc[2]*255),
                                                 int(alpha_fill*255))
                cc_hex_bndr = "#%2x%2x%2x%2x" % (int(cc[0]*255),
                                                 int(cc[1]*255),
                                                 int(cc[2]*255),
                                                 int(alpha_bndr*255))
            elif ii in sources:
                cc          = color[0]
                cc_hex_fill = "#%2x%2x%2x%2x" % (int(cc[0]*255),
                                                 int(cc[1]*255),
                                                 int(cc[2]*255),
                                                 int(alpha_fill*255))
                cc_hex_bndr = "#%2x%2x%2x%2x" % (int(cc[0]*255),
                                                 int(cc[1]*255),
                                                 int(cc[2]*255),
                                                 int(alpha_bndr*255))
            else:
                cc          = color[1]
                cc_hex_fill = "#%2x%2x%2x%2x" % (int(cc[0]*255),
                                                 int(cc[1]*255),
                                                 int(cc[2]*255),
                                                 int(alpha_fill*255))
                cc_hex_bndr = "#%2x%2x%2x%2x" % (int(cc[0]*255),
                                                 int(cc[1]*255),
                                                 int(cc[2]*255),
                                                 int(alpha_bndr*255))

            # plot node
            A.add_node(ii, color=cc_hex_bndr, fillcolor=cc_hex_fill,
                       style='filled')

        A.layout(prog='dot')
        A.draw(fname_plot)

    except ImportError:
        fig = plt.figure()

        pos   = nx.spring_layout(G)

        for ii in nodes:

            # determine color of node
            if ii in sinks:
                cc = color[2]
            elif ii in sources:
                cc = color[0]
            else:
                cc = color[1]

            # draw source node
            nx.draw_networkx_nodes(G, pos,
                                   nodelist=[ii],
                                   node_color=[cc],
                                   node_size=500,
                                   alpha=0.7)

        # draw edges
        nx.draw_networkx_edges(G, pos, alpha=0.7)

        # some math labels
        labels = {}
        for ii, inode in enumerate(nodes):
            labels[ii] = r'$'+str(inode)+'$'
        nx.draw_networkx_labels(G, pos, labels, font_size=16)

        plt.axis('off')
        fig.savefig(fname_plot)
        plt.close(fig)


def rotate(matrix):
    data = copy.deepcopy(matrix)
    for i in range(len(data)):
        for j in range(i+1, len(data[i])):
            temp = data[i][j]
            data[i][j] = data[j][i]
            data[j][i] = temp
    data = data[::-1]
    return data


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # nnodes    = 5
    # nedges    = 8
    # maxdegree = 3

    # # --------------------------
    # # create one network
    # # --------------------------
    # G = create_network(nnodes, nedges, maxdegree=maxdegree)
    # sink_nodes(G)
    # source_nodes(G)
    # plot_network(G, 'test_dag.pdf')

    # # --------------------------
    # # create 20 networks
    # # --------------------------
    # for inetwork in range(20):

    #     G = create_network(nnodes, nedges, maxdegree=maxdegree)

    #     # sources should be only "0" and should have at least one sink
    #     sources = source_nodes(G)
    #     sinks   = sink_nodes(G)

    #     # TODO add network to collection if not the same as an already existing: nx.is_isomorphic(G1, G2)

    #     # print some statistics
    #     print('number of nodes: ', G.number_of_nodes())
    #     print('number of edges: ', G.number_of_edges())
    #     print('source nodes:    ', sources)
    #     print('sink   nodes:    ', sinks)

    #     fname_plot = 'network_nodes_'+str(nnodes)+'_edges_'+str(nedges)+'_maxdeg_'+str(maxdegree)+'_'+str(inetwork+1)+'.pdf'
    #     plot_network(G,fname_plot)
