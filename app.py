from cytoscape import Cytoscape
import dash
import json
import time
from dash import dcc
from dash import html
import dash_cytoscape as cyto
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import numpy as np
import networkx as nx

genes = pd.read_csv('./BIOGRID-PROJECT-glioblastoma_project-GENES.projectindex.txt', sep="\\t",engine='python')
interactions = pd.read_csv('./BIOGRID-PROJECT-glioblastoma_project-INTERACTIONS.tab3.txt', sep="\\t",engine='python')


interactions = interactions[interactions['Entrez Gene Interactor A'] != '-']
interactions['Entrez Gene Interactor A'] = interactions['Entrez Gene Interactor A'].astype(int)



def networkGraph(nbr_edges):
    
    nodes = np.array(genes['ENTREZ GENE ID'])
    edges = np.array(interactions[['Entrez Gene Interactor A','Entrez Gene Interactor B']])
    edges = list(map(tuple,edges))

    df_attributes_e = interactions.drop(columns = ['Entrez Gene Interactor A','Entrez Gene Interactor B'])
    df_attributes_n = genes.drop(columns = ['ENTREZ GENE ID'])

    attributes_e = df_attributes_e.to_dict('records')
    attributes_n = df_attributes_n.to_dict('records')

    G = nx.Graph()
    dico_n = dict(zip(nodes,attributes_n))
    dico_e = dict(zip(edges,attributes_e))
    #print(edges)
    G.add_nodes_from(nodes)
    G.add_edges_from(edges[:nbr_edges])
    #print(G.edges)
    nx.set_node_attributes(G,dico_n)
    nx.set_edge_attributes(G,dico_e)
    return G
def  betweenness_centrality(G):
    return nx.betweenness_centrality(G)
def clustering_coefficient(G):
    return nx.average_clustering(G)
def minimum_spanning_tree(G):
    return nx.algorithms.minimum_spanning_edges(G)
def has_path(G,source,target):
    return nx.algorithms.has_path(G,source,target)
def shortest_path(G,source,target):
    return nx.algorithms.shortest_path(G,source,target)
def community(G):
    return nx.algorithms.community.greedy_modularity_communities(G)
G = networkGraph(500)
el = []

nodes = [
    {
        'data': {'id': str(id), 'label': '/','classes': 'black' },
        'position': {'x': 0, 'y': 0}
    }
    for id in G.nodes
]
edges = [
    {
        'data': {'source': str(source), 'target': str(target),'classes':'red'}, 
        'classes': 'red',
        'position': None
    }
    for source,target in G.edges
]
el = nodes+edges
app = dash.Dash(__name__)

genesId = []

genesOncogene = []
genesTumorSuppressor = []
genesCancerDriver = []

genesData = {}

styles = {
    'GRD' : {
        "width" : "0%",
        "height": "0vh",
        "margin": "0",
        "float": "left"
    },
    'CRC' : {
        "width" : "0%",
        "height": "0vh",
        "margin": "0",
        "float": "left"
    },
    'CCT' : {
        "width" : "0%",
        "height": "0vh",
        "margin": "0",
        "float": "left"
    },
    'BFT' : {
        "width" : "0%",
        "height": "0vh",
        "margin": "0",
        "float": "left"
    },
    'CSE' : {
        "width" : "100%",
        "height": "75vh",
        "margin": "0",
        "float": "left"
    }
}
"""
def loadData():
    for x in genes:
        genesId.append(x[0])
        el.append({
            'data': {'id': str(x[0]), 'label': str(x[3]), 'size': x[7]},
            'classes': 'known'
        })
        
        if (x[18] == "Oncogene"):
            genesOncogene.append(x)
        elif (x[18] == "Tumor Suppressor"):
            genesTumorSuppressor.append(x)
        elif (x[18] == "Cancer Driver"):
            genesCancerDriver.append(x)
        

    for x in interactions:
        if x[3] in genesId and x[4] in genesId:
            if x[3] != x[4]:
                if x[3] not in genesId:
                    el.append({
                        'data': {'id': str(x[3]), 'label': str(x[7])},
                        'classes': 'unknown'
                    })
                if x[4] not in genesId:
                    el.append({
                        'data': {'id': str(x[4]), 'label': str(x[8])},
                        'classes': 'unknown'
                    })
                el.append({
                    'data': {'source': str(x[3]), 'target': str(x[4])}, 
                    'classes': 'red',
                    'position': None
                })
"""

default_stylesheet=[
        # Class selectors
        {
            'selector': 'node',
            'style': {
                "width" : '7',
                "height" : '7'
            }
        },
        {
            'selector': 'edge',
            'style':{
                "line-color": "red",
                "width": 0.5
            }
        }
]

cyGrid = Cytoscape.newCyto('grid', el, "gridCytoscape")
cyConcentric = Cytoscape.newCyto('concentric', el, "concentricCytoscape")
cyCircle = Cytoscape.newCyto('circle', el, "circleCytoscape")
cyBreadthfirst = Cytoscape.newCyto('breadthfirst', el, "breadthfirstCytoscape")
cyCose = Cytoscape.newCyto('cose', el, "coseCytoscape")

global layoutActivated, addPressed, removePressed
layoutActivated = 1
addPressed = 0
removePressed = 0

def createBasicLayout():
    app.layout = html.Div([
        html.Div(className = "headerDiv",children=[
            html.H2("PROJECT : INFORMATION VISUALISATION"),
            dcc.Upload(
                id="upload-edges",
                multiple = False,
                children = html.Button('Import Interactions', id='import-edges', className="importButton")
            ),
            dcc.Upload(
                id="upload-nodes",
                multiple = False,
                children = html.Button('Import Genes', id='import-nodes', className="importButton")
            ),
        ]),

        html.Div(className = "leftSide", children = [
            html.Div(className = "layouts", children = [
                html.P("Layouts disponibles :"),
                
                dcc.Checklist(
                    id="checklistLayouts",
                    options=[
                        #available names = ('preset'), ('random'), 'grid', 'circle', 'concentric', 'breadthfirst', 'cose'
                        {'label': 'Grid', 'value': 'GRD'},
                        {'label': 'Circle', 'value': 'CRC'},
                        {'label': 'Concentric', 'value': 'CCT'},
                        {'label': 'Breadthfirst', 'value': 'BFT'},
                        {'label': 'Cose', 'value': 'CSE'}
                    ],
                    value = [],
                    labelStyle = {'display': 'block'}
                ),
                html.Button('Add layout', id="addLayout", n_clicks=0),
                html.Button('Remove layout', id="removeLayout", n_clicks=0)
            ]),
            html.Div(className = "metrics", children = [
                html.P("Metrics there."),
                dcc.Dropdown(
                    id='dropdown-metrics',
                    value='Choose a metric',
                    clearable=False,
                    options=[
                        {'label': name.capitalize(), 'value': name}
                        for name in ['Choose a metric','betweenness centrality','clustering coefficient','minimum spanning tree','shortest path','community']
                    ]
                ),
                html.Div(id='dd-output-container'),  
            ]),
        ]),
            
        html.Div(id="graphs", children = [
            dcc.Dropdown(
                id="dropdown1",
                className='dropdown',
                options=[
                    {'label': 'Grid', 'value': 'GRD'},
                    {'label': 'Circle', 'value': 'CRC'},
                    {'label': 'Concentric', 'value': 'CCT'},
                    {'label': 'Breadthfirst', 'value': 'BFT'},
                    {'label': 'Cose', 'value': 'CSE'}
                ],
                value='GRD'
            ),
            dcc.Dropdown(
                id="dropdown2",
                className='dropdown',
                options=[
                    {'label': 'Grid', 'value': 'GRD'},
                    {'label': 'Circle', 'value': 'CRC'},
                    {'label': 'Concentric', 'value': 'CCT'},
                    {'label': 'Breadthfirst', 'value': 'BFT'},
                    {'label': 'Cose', 'value': 'CSE'}
                ],
                value='GRD'
            ),
            dcc.Dropdown(
                id="dropdown3",
                className='dropdown',
                options=[
                    {'label': 'Grid', 'value': 'GRD'},
                    {'label': 'Circle', 'value': 'CRC'},
                    {'label': 'Concentric', 'value': 'CCT'},
                    {'label': 'Breadthfirst', 'value': 'BFT'},
                    {'label': 'Cose', 'value': 'CSE'}
                ],
                value='GRD'
            ),
            dcc.Dropdown(
                id="dropdown4",
                className='dropdown',
                options=[
                    {'label': 'Grid', 'value': 'GRD'},
                    {'label': 'Circle', 'value': 'CRC'},
                    {'label': 'Concentric', 'value': 'CCT'},
                    {'label': 'Breadthfirst', 'value': 'BFT'},
                    {'label': 'Cose', 'value': 'CSE'}
                ],
                value='GRD'
            ),
            dcc.Dropdown(
                id="dropdown5",
                className='dropdown',
                options=[
                    {'label': 'Grid', 'value': 'GRD'},
                    {'label': 'Circle', 'value': 'CRC'},
                    {'label': 'Concentric', 'value': 'CCT'},
                    {'label': 'Breadthfirst', 'value': 'BFT'},
                    {'label': 'Cose', 'value': 'CSE'}
                ],
                value='GRD'
            ),
            html.Div(id="graph1", children = [
                cyGrid
            ],
            style=styles['GRD']),
            html.Div(id="graph2", children = [
                cyCircle,
            ],
            style=styles['CRC']),
            html.Div(id="graph3", children = [
                cyConcentric,
            ],
            style=styles['CCT']),
            html.Div(id="graph4", children = [
                cyBreadthfirst,
            ],
            style=styles['BFT']),
            html.Div(id="graph5", children = [
                cyCose,
            ],
            style=styles['CSE'])
        ],
        ),
        html.Div(className = "filters", children = [
            html.P("Filters there.")
        ]),
        
        html.P(id='cytoscape-tapNodeData-json'),
        html.P(id="outputgraph"),
        html.P(id="outputnodes"),
        html.P(id="outputedges")
    ])

createBasicLayout()

@app.callback(  Output("graph1", 'style'),
                Output("graph2", 'style'),
                Output("graph3", 'style'),
                Output("graph4", 'style'),
                Output("graph5", 'style'),
                Output("dropdown1", 'style'),
                Output("dropdown2", 'style'),
                Output("dropdown3", 'style'),
                Output("dropdown4", 'style'),
                Output("dropdown5", 'style'),
                Input('addLayout', 'n_clicks'),
                Input('removeLayout', 'n_clicks'))
def update_output(add_n_clicks,remove_n_clicks):
    global layoutActivated, addPressed, removePressed
    if add_n_clicks != addPressed:
        layoutActivated = (layoutActivated + 1) if (layoutActivated < 5) else 5
        addPressed = add_n_clicks
    elif remove_n_clicks != removePressed:
        layoutActivated = (layoutActivated - 1) if (layoutActivated > 1) else 1
        removePressed = remove_n_clicks
    else:
        layoutActivated = 1
    width = str(99/layoutActivated) + "%"
    widthTab = {"width1":{"width":"0%", "display": "none"}, "width2":{"width":"0%", "display": "none"}, "width3":{"width":"0%", "display": "none"}, "width4":{"width":"0%", "display": "none"}, "width5":{"width":"0%", "display": "none"}}
    for i,wth in enumerate(widthTab):
        if(i < layoutActivated):
            widthTab[wth]["width"] = width
            widthTab[wth]["display"] = str("block")
        else:
            widthTab[wth]["width"] = "0%"
            widthTab[wth]["display"] = str("none")

    myReturn = {
            "width":widthTab["width1"]["width"], 
            "height":"65vh", 
            "float":"left", 
            "display" : widthTab["width1"]["display"]
        },{
            "width":widthTab["width2"]["width"], 
            "height":"65vh", 
            "float":"left", 
            "display" : widthTab["width2"]["display"]
        },{
            "width":widthTab["width3"]["width"], 
            "height":"65vh", 
            "float":"left", 
            "display" : widthTab["width3"]["display"]
        },{
            "width":widthTab["width4"]["width"], 
            "height":"65vh", 
            "float":"left", 
            "display" : widthTab["width4"]["display"]
        },{
            "width":widthTab["width5"]["width"], 
            "height":"65vh", 
            "float":"left", 
            "display" : widthTab["width5"]["display"]
        },{
            "width":widthTab["width1"]["width"], 
            "float":"left",
            "display" : widthTab["width1"]["display"]
        },{
            "width":widthTab["width2"]["width"], 
            "float":"left", 
            "display" : widthTab["width2"]["display"]
        },{
            "width":widthTab["width3"]["width"], 
            "float":"left", 
            "display" : widthTab["width3"]["display"]
        },{
            "width":widthTab["width4"]["width"], 
            "float":"left", 
            "display" : widthTab["width4"]["display"]
        },{
            "width":widthTab["width5"]["width"], 
            "float":"left", 
            "display" : widthTab["width5"]["display"]
        }
    return myReturn

@app.callback(  Output("gridCytoscape", 'layout'),
                Input('dropdown1', 'value'))
def update_output(value):
    layouts = {"GRD": "grid", "CRC": "circle", "CCT": "concentric", "BFT": "breadthfirst", "CSE": "cose"}
    return {
        'name': layouts[value],
        'animate': False
    }

@app.callback(Output('dd-output-container', 'children'),
              Output('gridCytoscape', 'stylesheet'), 
              Output('gridCytoscape', 'elements'),
              Input('dropdown-metrics', 'value'))
def update_metrics(metric):
    if (metric == 'betweenness centrality'):
        return 0,default_stylesheet,el
    if (metric == 'clustering coefficient'):
        return clustering_coefficient(G),default_stylesheet,el
    if (metric == 'minimum spanning tree'):
        edges = minimum_spanning_tree(G)
        edges = np.array(list(edges))[:,:2]
        new_edges = [
        {
            'data': {'source': str(source), 'target': str(target)}, 
            'classes': 'green',
            'position': None
        }
        for source,target in edges
        ]
        new_styles = [
        {
            
            'selector': '.green',
            'style':{
                "line-color": "green",
                "width": 0.5
            }
        }
        ]
        return "See minimum spanning tree on the figure",default_stylesheet+new_styles,el+new_edges
    if (metric == 'shortest path'):
        if (has_path(G, 5290, 308)):
            path = shortest_path(G,5290,308)
            edges = []
            length = len(path)
            for i in range (length-1):
                edges.append((path[i],path[i+1]))
            new_edges = [
            {
                'data': {'source': str(source), 'target': str(target)}, 
                'classes': 'green',
                'position': None
            }
            for source,target in edges
            ]
            new_styles = [
            {
                'selector': 'node',
                'style':{
                    'background-color': 'green',
                } 
            },   
            {               
                
                'selector': '.green',
                'style':{
                    "line-color": "green",
                    "width": 0.5
                }
            }
            ]
            return "The length of the path between source and target is "+str(length),default_stylesheet+new_styles,el+new_edges
        else:
            return "No path between source and target",default_stylesheet,el
    if (metric == 'community'):
        return 0,default_stylesheet,el
    else :
        return "Choose the proprety",default_stylesheet,el

"""
@app.callback(
                Output("graph1", 'style'),
                Output("graph2", 'style'),
                Output("graph3", 'style'),
                Output("graph4", 'style'),
                Output("graph5", 'style'),
                Input('checklistLayouts', 'value'))
def update_output(n_clicks):
    basicStyle = {
        'width': '0%', 
        'height': '100%',
        'border': '0px solid grey'
    }
    grdDisplay = "block"
    crcDisplay = "block"
    cctDisplay = "block"
    bftDisplay = "block"
    cseDisplay = "block"

    if "GRD" in n_clicks:
        cyGrid.style["width"] = str(99/len(n_clicks)) + "%"
    else:
        cyGrid.style["width"] = "0%"
        grdDisplay = "None"
    if "CRC" in n_clicks:
        cyCircle.style["width"] = str(99/len(n_clicks)) + "%"
    else:
        cyCircle.style["width"] = "0%"
        crcDisplay = "None"
    if "CCT" in n_clicks:
        cyConcentric.style["width"] = str(99/len(n_clicks)) + "%"
    else:
        cyConcentric.style["width"] = "0%"
        cctDisplay = "None"
    if "BFT" in n_clicks:
        cyBreadthfirst.style["width"] = str(99/len(n_clicks)) + "%"
    else:
        cyBreadthfirst.style["width"] = "0%"
        bftDisplay = "None"
    if "CSE" in n_clicks:
        cyCose.style["width"] = str(100/len(n_clicks)) + "%"
    else:
        cyCose.style["width"] = "0%"    
        cseDisplay = "None"
    grdWidth = {"width": cyGrid.style["width"], "height":"100%", "float":"left", "display":grdDisplay}
    crcWidth = {"width": cyCircle.style["width"], "height":"100%", "float":"left", "display":crcDisplay}
    cctWidth = {"width": cyConcentric.style["width"], "height":"100%", "float":"left", "display":cctDisplay}
    bftWidth = {"width": cyBreadthfirst.style["width"], "height":"100%", "float":"left", "display":bftDisplay}
    cseWidth = {"width": cyCose.style["width"], "height":"100%", "float":"left", "display":cseDisplay}

    print(grdWidth, crcWidth,cctWidth,bftWidth,cseWidth)
    return grdWidth, crcWidth, cctWidth,bftWidth,cseWidth
"""
"""
@app.callback(  Output('graph1', 'style'),
                Output('graph2', 'style'),
                Output('graph3', 'style'),
                Output('graph4', 'style'),
                Output('graph5', 'style'),
                Output('graph1', 'children'),
                Output('graph2', 'children'),
                Output('graph3', 'children'),
                Output('graph4', 'children'),
                Output('graph5', 'children'),
                Input('checklistLayouts', 'value'))
def update_output(n_clicks):
    newStyles = {
        'GRD' : {
            "width" : "0%",
            "height": "0vh",
            "margin": "0",
            "float": "left",
        },
        'CRC' : {
            "width" : "0%",
            "height": "0vh",
            "margin": "0",
            "float": "left",
        },
        'CCT' : {
            "width" : "0%",
            "height": "0vh",
            "margin": "0",
            "float": "left",
        },
        'BFT' : {
            "width" : "0%",
            "height": "0vh",
            "margin": "0",
            "float": "left",
        },
        'CSE' : {
            "width" : "0%",
            "height": "0vh",
            "margin": "0",
            "float": "left",
        }
    }

    if len(n_clicks) == 0:
        return styles["GRD"], styles["CRC"], styles["CCT"], styles["BFT"], styles["CSE"], None, None, None, None, None
    order = ["GRD", "CRC", "CCT", "BFT", "CSE"]
    n_clicks = sorted(n_clicks, key=lambda x:order.index(x))
    width = "0%"
    height = "0vh"
    lastWidth = "100%"
    lastHeight = "75vh"
    if len(n_clicks) == 2:
        width = "49%"
        height = "75vh"
        lastWidth = "49%"
        lastHeight = "75vh"
    if len(n_clicks) == 3:
        width = "49%"
        height = "37vh"
        lastWidth = "100%"
        lastHeight = "37vh"
    if len(n_clicks) == 4:
        width = "49%"
        height = "37vh"
        lastWidth = "49%"
        lastHeight = "37vh"
    if len(n_clicks) == 5:
        width = "50%"
        height = "24vh"
        lastWidth = "100%"
        lastHeight = "25vh"

    for graph in styles:
        if graph in n_clicks:
            if graph != n_clicks[len(n_clicks) - 1]:
                newStyles[graph]["width"] = width
                newStyles[graph]["height"] = height
            else:
                newStyles[graph]["width"] = lastWidth
                newStyles[graph]["height"] = lastHeight
        else:
            newStyles[graph]["width"] = "0%"
            newStyles[graph]["height"] = "0vh"  

    return newStyles["GRD"], newStyles["CRC"], newStyles["CCT"], newStyles["BFT"], newStyles["CSE"], cyGrid, cyCircle, cyConcentric, cyBreadthfirst, cyCose
"""

@app.callback(Output('outputnodes', 'children'),
              Input('upload-nodes', 'contents'),
              State('upload-nodes', 'filename'))
def update_output(list_of_contents, list_of_names):
    if list_of_contents is not None:
        print(pd.read_csv(list_of_names).values)

@app.callback(Output('outputedges', 'children'),
              Input('upload-edges', 'contents'),
              State('upload-edges', 'filename'))
def update_output(list_of_contents, list_of_names):
    if list_of_contents is not None:
        print(pd.read_csv(list_of_names).values)

"""
@app.callback(Output('cytoscape-tapNodeData-json', 'children'),
              Input('cytoscape-event-callbacks-1', 'tapNodeData'))
def displayTapNodeData(data):
    if json.dumps(data, indent=2) != "null":
        print(type(json.dumps(data, indent=2)))
        #return json.dumps(data, indent=2)
"""

#method_list = [method for method in dir(utils) if method.startswith('__') is False]
#print(method_list)


app.run_server(debug=True)