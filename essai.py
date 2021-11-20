import plotly.graph_objects as go
import networkx as nx
import pandas as pd
import numpy as np

genes = pd.read_csv('./BIOGRID-PROJECT-glioblastoma_project-GENES.projectindex.txt', sep="\\t")
interactions = pd.read_csv('./BIOGRID-PROJECT-glioblastoma_project-INTERACTIONS.tab3.txt', sep="\\t")

interactions = interactions.values
genes = genes.values

tenth = int(len(interactions)/10)
tenthInteractions = interactions[:tenth]

genesId = []

G = nx.random_geometric_graph(10, 0.125)
myG = nx.Graph()

for x in genes:
    myG.add_nodes_from([
        (str(x[0]), {'label' : x[3], 'interactions' : x[7]})
    ])
    genesId.append(x[0])

for x in tenthInteractions:
    if x[3] == x[4]:
        continue
    if x[3] not in genesId:
        continue
        myG.add_nodes_from([
            (str(x[3]), {'label' : '/'})
        ])
        genesId.append(x[3])
    if x[4] not in genesId:
        continue
        myG.add_nodes_from([
            (str(x[4]), {'label' : '/'})
        ])
        genesId.append(x[4])
    myG.add_edge(x[3], x[4])

print(G)
print(G.nodes)
print(myG.nodes)
for node in G.nodes():
    #print(myG.nodes)
    #print(G.nodes)
    print(node)
    break  


"""
edge_x = []
edge_y = []
for edge in myG.edges():
    x0, y0 = myG.nodes[edge[0]]['pos']
    x1, y1 = myG.nodes[edge[1]]['pos']
    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None)
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)

edge_trace = go.Scatter(
    x=edge_x, y=edge_y,
    line=dict(width=0.5, color='#888'),
    hoverinfo='none',
    mode='lines')

node_x = []
node_y = []
for node in myG.nodes():
    x, y = myG.nodes[node]['pos']
    node_x.append(x)
    node_y.append(y)

node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers',
    hoverinfo='text',
    marker=dict(
        showscale=True,
        # colorscale options
        #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
        #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
        #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
        colorscale='YlGnBu',
        reversescale=True,
        color=[],
        size=10,
        colorbar=dict(
            thickness=15,
            title='Node Connections',
            xanchor='left',
            titleside='right'
        ),
        line_width=2
    )
)

node_adjacencies = []
node_text = []
for node, adjacencies in enumerate(myG.adjacency()):
    node_adjacencies.append(len(adjacencies[1]))
    node_text.append('# of connections: '+str(len(adjacencies[1])))

node_trace.marker.color = node_adjacencies
node_trace.text = node_text

fig = go.Figure(data=[edge_trace, node_trace],
             layout=go.Layout(
                title='<br>Network graph made with Python',
                titlefont_size=16,
                showlegend=False,
                hovermode='closest',
                margin=dict(b=20,l=5,r=5,t=40),
                annotations=[ dict(
                    text="Python code: <a href='https://plotly.com/ipython-notebooks/network-graphs/'> https://plotly.com/ipython-notebooks/network-graphs/</a>",
                    showarrow=False,
                    xref="paper", yref="paper",
                    x=0.005, y=-0.002 ) ],
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )
fig.show()
"""