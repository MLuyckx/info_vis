from cytoscape import Cytoscape
import dash
import dash_bootstrap_components as dbc
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

#genes = pd.read_csv('./BIOGRID-PROJECT-glioblastoma_project-GENES.projectindex.txt', sep="\\t",engine='python')
#interactions = pd.read_csv('./BIOGRID-PROJECT-glioblastoma_project-INTERACTIONS.tab3.txt', sep="\\t",engine='python')



def networkGraph(nbr_edges, genes, interactions):
    newGenes = genes.copy()
    newGenes = newGenes.iloc[len(newGenes.index):,:]
    edges = np.array(interactions[['Entrez Gene Interactor A','Entrez Gene Interactor B']])
    for edge in edges[:nbr_edges]:
        toAdd1 = genes[genes["ENTREZ GENE ID"] == edge[0]]
        newGenes = newGenes.append(toAdd1)
    nodes = np.array(newGenes['ENTREZ GENE ID'])

    edges = list(map(tuple,edges))

    df_attributes_e = interactions.drop(columns = ['Entrez Gene Interactor A','Entrez Gene Interactor B'])
    df_attributes_n = newGenes.drop(columns = ['ENTREZ GENE ID'])

    attributes_e = df_attributes_e.to_dict('records')
    attributes_n = df_attributes_n.to_dict('records')

    G = nx.Graph()
    dico_n = dict(zip(nodes,attributes_n))
    dico_e = dict(zip(edges,attributes_e))

    G.add_nodes_from(nodes)
    G.add_edges_from(edges[:nbr_edges])

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

genes = pd.read_csv('./genes.csv', sep=";")
interactions = pd.read_csv('interactions.csv', sep=";")

interactions = interactions[interactions['Entrez Gene Interactor A'] != '-']
interactions['Entrez Gene Interactor A'] = interactions['Entrez Gene Interactor A'].astype(int)

global el,currentGraph, filterSelected
el = []
currentGraph = networkGraph(1000, genes, interactions)
filterSelected = None
def initialisation(G):
    global el
    nodes = [
        {
            'data': {'id': str(id), 'label': '/', 'size': G.degree[id]},
            'classes': 'black'
        }
        for id in G.nodes
    ]
    edges = [
        {
            'data': {'id':str(data['#BioGRID Interaction ID']), 'source':str(source), 'target':str(target)}, 
            'classes': 'red'
        }
        for source,target, data in G.edges.data()
    ]
    el = nodes+edges

initialisation(currentGraph)

global cyGrid,cyConcentric, cyCircle, cyBreadthfirst, cyCose
cyGrid = Cytoscape.newCyto('grid', el, "gridCytoscape")
cyConcentric = Cytoscape.newCyto('concentric', el, "concentricCytoscape")
cyCircle = Cytoscape.newCyto('circle', el, "circleCytoscape")
cyBreadthfirst = Cytoscape.newCyto('breadthfirst', el, "breadthfirstCytoscape")
cyCose = Cytoscape.newCyto('cose', el, "coseCytoscape")

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



app = dash.Dash(__name__)

genesId = []

genesOncogene = []
genesTumorSuppressor = []
genesCancerDriver = []

genesData = {}

default_stylesheet=[
    # Class selectors
    {
        'selector': 'node',
        'style': {
            "width" : "mapData(size, 0, 100, 3, 100)",
            "height" : "mapData(size, 0, 100, 3, 100)"
        }
    },
    {
        'selector': 'edge',
        'style':{
            "line-color": "red",
            "width":1
        }
    },
    {
        'selector': '.known',
        'style':{
            'background-color': "blue",
            "opacity": 0.9
        }
    },
    {
        'selector': '.unknown',
        'style':{
            "color": "grey",
            "opacity": 0.5
        }
    }
]

global layoutActivated, addPressed, removePressed
layoutActivated = 1
addPressed = 0
removePressed = 0

def createBasicLayout():
    app.layout = html.Div([
        html.Div(className = "headerDiv",children=[
            html.H2("PROJECT : INFORMATION VISUALISATION"),
            html.Button('Reload data', id='reloadData', className="reloadButton"),
            html.P(id='pid'),
            dcc.Upload(
                id="upload-edges",
                multiple = False,
                children = html.Button('Import Interactions', id='import-edges', className="importButton")
            ),
            dcc.Upload(
                id="upload-nodes",
                multiple = False,
                children = html.Button('Import Genes', id='import-nodes', className="importButton")
            )
        ]),

        html.Div(className = "leftSide", children = [
            html.Div(className = "layouts", children = [
                html.P("Affichage des layouts"),
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
                        for name in ['betweenness centrality','clustering coefficient','minimum spanning tree','shortest path','community']
                    ],
                    placeholder="Select a metric"
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
                placeholder="Select a layout"
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
                placeholder="Select a layout"
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
                placeholder="Select a layout"
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
                placeholder="Select a layout"
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
                placeholder="Select a layout"
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

        html.Div(className = "infoFilters", children = [
            html.P("Nothing there for the moment", id="infoFiltersText", className="infoFiltersText")
        ]),

        html.Div(className = "filters", children = [
            html.Button(
                html.Img(
                    src="./assets/star.png",
                    id="imageNeigColor"
                ),
                id='neigColor', 
                className="filterButton",
                title="Show neighbors"
            ),
            html.Button(
                html.Img(
                    src="./assets/classes.png",
                    id="imageClassColor"
                ),
                id='classColor', 
                className="filterButton",
                title="Show subcategories"
            ),
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
    if value==None:
        return None
    layouts = {"GRD": "grid", "CRC": "circle", "CCT": "concentric", "BFT": "breadthfirst", "CSE": "cose"}
    return {
        'name': layouts[value],
        'animate': False
    }

@app.callback(  Output("circleCytoscape", 'layout'),
                Input('dropdown2', 'value'))
def update_output(value):
    if value==None:
        return None
    layouts = {"GRD": "grid", "CRC": "circle", "CCT": "concentric", "BFT": "breadthfirst", "CSE": "cose"}
    return {
        'name': layouts[value],
        'animate': False
    }

@app.callback(  Output("concentricCytoscape", 'layout'),
                Input('dropdown3', 'value'))
def update_output(value):
    if value==None:
        return None
    layouts = {"GRD": "grid", "CRC": "circle", "CCT": "concentric", "BFT": "breadthfirst", "CSE": "cose"}
    return {
        'name': layouts[value],
        'animate': False
    }

@app.callback(  Output("breadthfirstCytoscape", 'layout'),
                Input('dropdown4', 'value'))
def update_output(value):
    if value==None:
        return None
    layouts = {"GRD": "grid", "CRC": "circle", "CCT": "concentric", "BFT": "breadthfirst", "CSE": "cose"}
    return {
        'name': layouts[value],
        'animate': False
    }

@app.callback(  Output("coseCytoscape", 'layout'),
                Input('dropdown5', 'value'))
def update_output(value):
    if value==None:
        return None
    layouts = {"GRD": "grid", "CRC": "circle", "CCT": "concentric", "BFT": "breadthfirst", "CSE": "cose"}
    return {
        'name': layouts[value],
        'animate': False
    }


@app.callback(Output('infoFiltersText', 'children'),
              Input('neigColor', 'n_clicks'),
              Input('classColor', 'n_clicks'))
def update_output(neigColor_clicks, classColor_clicks):
    global filterSelected
    ctx = dash.callback_context
    inputType = ctx.triggered[0]['prop_id'].split('.')[0]
    if inputType == 'neigColor':
        if filterSelected == 'neigColor':
            filterSelected = None
            return " "
        else:
            filterSelected = 'neigColor'
        return "Select a node on a graph."
    elif inputType == 'classColor':
        if filterSelected == 'classColor':
            filterSelected = None
            return " "
        else:
            filterSelected = 'classColor'
        return "Green : Oncogene |Yellow : Tumor Suppressor | Blue : Cancer Driver"
    else:
        return " "

@app.callback(Output('dd-output-container', 'children'),
              Output('gridCytoscape', 'stylesheet'), 
              Output('concentricCytoscape', 'stylesheet'),
              Output('circleCytoscape', 'stylesheet'),
              Output('breadthfirstCytoscape', 'stylesheet'),
              Output('coseCytoscape', 'stylesheet'),
              Output('gridCytoscape', 'elements'),
              Output('concentricCytoscape', 'elements'),
              Output('circleCytoscape', 'elements'),
              Output('breadthfirstCytoscape', 'elements'),
              Output('coseCytoscape', 'elements'),
              Input('dropdown-metrics', 'value'),
              Input('reloadData', 'n_clicks'),
              Input('gridCytoscape', 'tapNode'),
              Input('concentricCytoscape', 'tapNode'),
              Input('circleCytoscape', 'tapNode'),
              Input('breadthfirstCytoscape', 'tapNode'),
              Input('coseCytoscape', 'tapNode'),
              Input('classColor', 'n_clicks'))
def update_metrics(metric, n_clicks, node1, node2, node3, node4, node5, classColor):
    ctx = dash.callback_context
    inputType = ctx.triggered[0]['prop_id'].split('.')[0]
    global currentGraph
    ############################################################################
    
    # Metrics calculation

    ############################################################################
    if inputType == "dropdown-metrics":
        if (metric == 'betweenness centrality'):
            return 0,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el
        if (metric == 'clustering coefficient'):
            return clustering_coefficient(currentGraph),default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el
        if (metric == 'minimum spanning tree'):
            edges = minimum_spanning_tree(currentGraph)
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
                    "width": 1
                }
            }
            ]
            return "See minimum spanning tree on the figure", default_stylesheet+new_styles,default_stylesheet+new_styles,default_stylesheet+new_styles,default_stylesheet+new_styles,default_stylesheet+new_styles, el+new_edges, el+new_edges, el+new_edges, el+new_edges, el+new_edges

        if (metric == 'shortest path'):
            if (has_path(currentGraph, 5290, 308)):
                path = shortest_path(currentGraph,5290,308)
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
                        "width": 1
                    }
                }
                ]
                return "The length of the path between source and target is "+str(length), default_stylesheet+new_styles, default_stylesheet+new_styles,default_stylesheet+new_styles,default_stylesheet+new_styles,default_stylesheet+new_styles,el+new_edges, el+new_edges, el+new_edges, el+new_edges, el+new_edges
            else:
                return "No path between source and target", default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, el, el, el, el, el
                
        if (metric == 'community'):
            return 0, default_stylesheet, el, el, el, el, el
        else :
            return "Choose the proprety", default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, el, el, el, el, el
    ############################################################################
    
    # Reload Data

    ############################################################################
    elif inputType == "reloadData":
        global newNodesData, newEdgesData
        if len(newNodesData) and len(newEdgesData):
            newGraph = networkGraph(100, newNodesData, newEdgesData)
            initialisation(newGraph)
            currentGraph = newGraph
        return "Choose the proprety", default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, el, el, el, el, el
    ############################################################################
    
    # Filters processing

    ############################################################################
    elif inputType in ['gridCytoscape','concentricCytoscape','circleCytoscape','breadthfirstCytoscape','coseCytoscape']:
        if filterSelected == 'neigColor':
            newStyle = []
            checkNode = {
                'gridCytoscape': node1,
                'concentricCytoscape': node2,
                'circleCytoscape': node3,
                'breadthfirstCytoscape': node4,
                'coseCytoscape': node5
            }
            finalNode = checkNode[inputType]
            if finalNode:
                for edge in finalNode["edgesData"]:
                    if edge['source'] == finalNode['data']['id']:
                        newStyle.append({
                            "selector": 'node[id = "{}"]'.format(edge['target']),
                            "style": {
                                'background-color': "green",
                                'opacity': 0.9
                            }
                        })
                        newStyle.append({
                            "selector": 'edge[id= "{}"]'.format(edge['id']),
                            "style": {
                                "line-color": "green",
                                'opacity': 0.9,
                                'z-index': 5000
                            }
                        })
                    if edge['target'] == finalNode['data']['id']:
                        newStyle.append({
                            "selector": 'node[id = "{}"]'.format(edge['source']),
                            "style": {
                                'background-color': "green",
                                'opacity': 0.9
                            }
                        })
                        newStyle.append({
                            "selector": 'edge[id= "{}"]'.format(edge['id']),
                            "style": {
                                "line-color": "green",
                                'opacity': 0.9,
                                'z-index': 5000
                            }
                        })
                return "Choose the proprety", default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, el, el, el, el, el
            else:
                return "Choose the proprety", default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el
        else: 
           return "Choose the proprety",default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el 
    elif inputType == 'classColor':
        if filterSelected is None:
           return "Choose the proprety", default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, el, el, el, el, el 
        newStyle = []
        for node in currentGraph.nodes.data():
            if node[1] == {}:
                return "Choose the proprety", default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, el, el, el, el, el
            if node[1]['SUBCATEGORY VALUES'] == "Oncogene":
                newStyle.append({
                    "selector": 'node[id = "{}"]'.format(node[0]),
                    "style": {
                        'background-color': "green",
                        'opacity': 0.9
                    }
                })
            elif node[1]['SUBCATEGORY VALUES'] == "Tumor suppressor":
                newStyle.append({
                    "selector": 'node[id = "{}"]'.format(node[0]),
                    "style": {
                        'background-color': "yellow",
                        'opacity': 0.9
                    }
                })
            elif node[1]['SUBCATEGORY VALUES'] == "Cancer driver":
                newStyle.append({
                    "selector": 'node[id = "{}"]'.format(node[0]),
                    "style": {
                        'background-color': "blue",
                        'opacity': 0.9
                    }
                })
        return "Choose the proprety", default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, el, el, el, el, el
    else:
        return "Choose the proprety",default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el



global newNodesData, newEdgesData
newNodesData = []
newEdgesData = []

@app.callback(Output('import-nodes', 'children'),
              Input('upload-nodes', 'contents'),
              State('upload-nodes', 'filename'))
def update_output(list_of_contents, list_of_names):
    global newNodesData
    if list_of_contents is not None:
        newNodesData = pd.read_csv(list_of_names, sep=";")
        if list_of_names:
            return list_of_names
        else:
            return "Import Genes"
    return "Import Genes"

@app.callback(Output('import-edges', 'children'),
              Input('upload-edges', 'contents'),
              State('upload-edges', 'filename'))
def update_output(list_of_contents, list_of_names):
    global newEdgesData
    if list_of_contents is not None:
        newEdgesData = pd.read_csv(list_of_names, sep=";")
        newEdgesData = newEdgesData[newEdgesData['Entrez Gene Interactor A'] != '-']
        newEdgesData['Entrez Gene Interactor A'] = newEdgesData['Entrez Gene Interactor A'].astype(int)
        if list_of_names:
            return list_of_names
        else:
            return "Import Interactions"
    return "Import Interactions"

app.run_server(debug=True)