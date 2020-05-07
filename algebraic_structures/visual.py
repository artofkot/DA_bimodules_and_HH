# -*- coding: utf-8 -*- 
from .da_bimodule import  Bunch_of_arrows, DA_bimodule,cancel_pure_differential,da_arrow_to_str
from .da_bimodule import  da_randomly_cancel_until_possible
from .da_bimodule  import da_in_mod_gen, da_out_mod_gen, da_in_alg_tuple, da_out_alg_gen
import sys
reload(sys)
sys.setdefaultencoding('utf8')
import networkx as nx
import matplotlib.pyplot as plt
import pygraphviz as gv
from networkx.drawing.nx_agraph import write_dot

def arrow_to_label(arrow):
    if len(arrow)==4:
        return str(da_out_alg_gen(arrow)) + '-' + str(da_in_alg_tuple(arrow))
        

def draw_DA_bimodule(DA1):
    labels=[]
    graph=[]

    for generator1 in DA1.genset:
        for generator2 in DA1.genset:
            arrows=[arrow for arrow in DA1.da_arrows if (da_in_mod_gen(arrow)==generator1 and da_out_mod_gen(arrow)==generator2)]

            if len(arrows)!=0:
                x=''
                for ind, arrow in enumerate(arrows):
                    x= x + str(da_out_alg_gen(arrow)) + '_' + str(da_in_alg_tuple(arrow))
                    if ind+1!=len(arrows): x=x+' + '
                # print str(generator1) + '--' + x + '-->'+ str(generator2)
                labels.append(x)
                graph.append((da_in_mod_gen(arrow),da_out_mod_gen(arrow)))


    # labels=[arrow_to_label(arrow) for arrow in DA1.arrows]
    # graph=[(da_in_mod_gen(arrow),da_out_mod_gen(arrow)) for arrow in DA1.da_arrows]

    # draw_graph_using_gv(graph, nodes=DA1.genset, labels=labels)
    draw_graph(graph, nodes=DA1.genset, labels=labels)

def draw_D_bimodule(D):
    labels=[]
    graph=[]

    for generator1 in D.genset:
        for generator2 in D.genset:
            arrows=[arrow for arrow in D.left_d_arrows if (arrow[0]==generator1 and arrow[2]==generator2)]
            x=''
            for ind, arrow in enumerate(arrows):
                x= x + str(arrow[1])
                if ind+1!=len(arrows): x=x+' + '
            # print str(generator1) + '--' + x + '-->'+ str(generator2)
            labels.append(x)
            graph.append((arrow[0],arrow[2]))


    # labels=[arrow_to_label(arrow) for arrow in DA1.arrows]
    # graph=[(da_in_mod_gen(arrow),da_out_mod_gen(arrow)) for arrow in DA1.da_arrows]

    # draw_graph_using_gv(graph, nodes=DA1.genset, labels=labels)
    draw_graph(graph, nodes=D.genset, labels=labels)

def draw_chain_complex(C):
    graph=[(arrow[0],arrow[1]) for arrow in C.arrows]

    # draw_graph_using_gv(graph, nodes=DA1.genset, labels=labels)
    draw_graph(graph,nodes=C.genset)

# def draw_graph_using_gv(graph, nodes=None, labels=None, graph_layout='shell',
#                node_size=0, node_color='blue', node_alpha=0.8,
#                node_text_size=11,
#                edge_label_text_size=7,
#                edge_color='black', edge_alpha=1, edge_tickness=0.3,
#                edge_text_pos=0.4,
#                text_font='sans-serif'):
#     G = gv.AGraph(strict=False,directed=True)

#     for node in nodes:
#         G.add_node(node,label=node.name)

#     for i,edge in enumerate(graph):
#         G.add_edge(edge,label=labels[i],
#         arrowhead='open',
#         fontname='Courier',
#         fontsize='6')

#     G.layout(prog='circo') # use dot
#     G.draw('file.png')


def draw_graph(graph, nodes=None, labels=None, graph_layout='shell',
               node_size=0, node_color='blue', node_alpha=0.8,
               node_text_size=11,
               edge_label_text_size=7,
               edge_color='black', edge_alpha=1, edge_tickness=0.3,
               edge_text_pos=0.4,
               text_font='sans-serif'):

    # create networkx graph
    G=nx.MultiDiGraph()

    # add nodes
    if nodes!=None:
        G.add_nodes_from(nodes)

    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])

    # these are different layouts for the network you may try
    # shell seems to work best
    if graph_layout == 'spring':
        graph_pos=nx.spring_layout(G)
    elif graph_layout == 'spectral':
        graph_pos=nx.spectral_layout(G)
    elif graph_layout == 'random':
        graph_pos=nx.random_layout(G)
    else:
        graph_pos=nx.shell_layout(G)


    # draw graph
    nx.draw_networkx_nodes(G,graph_pos,node_size=node_size, 
                           alpha=node_alpha, node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,font_size=node_text_size,
                            font_family=text_font)

    if labels is None:
        labels = range(len(graph))

    edge_labels = dict(zip(graph, labels))
    nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels, 
                                 label_pos=edge_text_pos,font_size=edge_label_text_size,
                                 font_color='red')
    write_dot(G,'multi.dot')
    # A=nx.to_agraph(G)

    # show graph
    plt.show()

# graph = [(0, 1), (1, 5), (1, 7), (4, 5), (4, 8), (1, 6), (3, 7), (5, 9),
#          (2, 4), (0, 4), (2, 5), (3, 6), (8, 9)]

# # you may name your edge labels
# labels = map(chr, range(65, 65+len(graph)))
# print labels
# draw_graph(graph, labels)

# # if edge labels is not specified, numeric labels (0, 1, 2...) will be used
# # draw_graph(graph)