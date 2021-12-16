from cytoscape import Cytoscape
import dash
import dash_bootstrap_components as dbc
import json
import time
import random
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
def betweenness_centrality(G):
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

genesColumnsSublclasses = ["CATEGORY VALUES", "SUBCATEGORY VALUES"]
interactionsColumnsSubclasses = ["Experimental System", "Experimental System Type", "Author", "Throughput", "Modification", "Ontology Term Categories"]

global el,currentGraph, filterSelected
el = []
currentGraph = networkGraph(500, genes, interactions)
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
            #"width" : "mapData(size, 0, 100, 3, 100)",
            #"height" : "mapData(size, 0, 100, 3, 100)"
            "width" : 12,
            "height" : 12
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
                dcc.Input(id="input_node_1", type="text", placeholder="",  style={'display': 'none'}),
                dcc.Input(id="input_node_2", type="text", placeholder="", debounce=True,style={'display': 'none'}),
                html.Button(id='submit-button-state', n_clicks=0, children='Submit',hidden= True)
            ]),
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
            html.Div(className = "nodeColor", children = [
                html.Button(
                    html.Img(
                        src="./assets/classes.png",
                        id="imageClassColor"
                    ),
                    id='classColor', 
                    className="filterButton",
                    title="Show subcategories"
                ),
                dcc.Dropdown(
                    id="dropdownClasses",
                    className='filterDropdown',
                    options=[{'label': name, 'value' : name} for name in genesColumnsSublclasses],
                    optionHeight = 50,
                    placeholder="Select an attribute"
                )
            ]),
            html.Div(className = "edgeColor", children = [
                html.Button(
                    html.Img(
                        src="./assets/edges.png",
                        id="imageEdgeColor"
                    ),
                    id='edgeColor', 
                    className="filterButton",
                    title="Show edge subcategories"
                ),
                dcc.Dropdown(
                    id="dropdownEdges",
                    className='filterDropdown',
                    options=[{'label': name, 'value' : name} for name in interactionsColumnsSubclasses],
                    optionHeight = 50,
                    placeholder="Select an attribute"
                )
            ]),
        ]),
            
        html.Div(className = "infoFilters", children = [
            html.P("Nothing there for the moment", id="infoFiltersText", className="infoFiltersText")
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

        html.Div(id="dataInformation", children= [
            html.Table(id="dataInformationText")
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
    widthDropdown = str((99/layoutActivated)*2) + "%"
    widthTab = {"width1":{"width":"0%", "display": "none", "widthDropdown":"0%"}, "width2":{"width":"0%", "display": "none", "widthDropdown":"0%"}, "width3":{"width":"0%", "display": "none", "widthDropdown":"0%"}, "width4":{"width":"0%", "display": "none", "widthDropdown":"0%"}, "width5":{"width":"0%", "display": "none", "widthDropdown":"0%"}}
    for i,wth in enumerate(widthTab):
        if(i < layoutActivated):
            widthTab[wth]["width"] = width
            widthTab[wth]["display"] = str("block")
            widthTab[wth]["widthDropdown"] = widthDropdown
        else:
            widthTab[wth]["width"] = "0%"
            widthTab[wth]["display"] = str("none")
            widthTab[wth]["widthDropdown"] = "0%"


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
              Output('neigColor', 'children'),
              Output('classColor', 'children'),
              Output('edgeColor', 'children'),
              Input('neigColor', 'n_clicks'),
              Input('classColor', 'n_clicks'),
              Input('edgeColor', 'n_clicks'))
def update_output(neigColor_clicks, classColor_clicks, edgeColor_clicks):
    global filterSelected
    ctx = dash.callback_context
    inputType = ctx.triggered[0]['prop_id'].split('.')[0]
    if inputType == 'neigColor':
        if filterSelected == 'neigColor':
            filterSelected = None
            return " ", html.Img(
                src="./assets/star.png",
                id="imageNeigColor"
            ),html.Img(
                src="./assets/classes.png",
                id="imageClassColor"
            ),html.Img(
                src="./assets/edges.png",
                id="imageEdgeColor"
            )   
        else:
            filterSelected = 'neigColor'
        return "Select a node on a graph.", html.Img(
                src="./assets/starClicked.png",
                id="imageNeigColor"
            ),html.Img(
                src="./assets/classes.png",
                id="imageClassColor"
            ),html.Img(
                src="./assets/edges.png",
                id="imageEdgeColor"
            ) 
    elif inputType == 'classColor':
        if filterSelected == 'classColor':
            filterSelected = None
            return " ", html.Img(
                src="./assets/star.png",
                id="imageNeigColor"
            ),html.Img(
                src="./assets/classes.png",
                id="imageClassColor"
            ),html.Img(
                src="./assets/edges.png",
                id="imageEdgeColor"
            )  
        else:
            filterSelected = 'classColor'
        return " ", html.Img(
                src="./assets/star.png",
                id="imageNeigColor"
            ),html.Img(
                src="./assets/classesClicked.png",
                id="imageClassColor"
            ),html.Img(
                src="./assets/edges.png",
                id="imageEdgeColor"
            )
    elif inputType == 'edgeColor':
        if filterSelected == 'edgeColor':
            filterSelected = None
            return " ", html.Img(
                src="./assets/star.png",
                id="imageNeigColor"
            ),html.Img(
                src="./assets/classes.png",
                id="imageClassColor"
            ),html.Img(
                src="./assets/edges.png",
                id="imageEdgeColor"
            )
        else:
            filterSelected = 'edgeColor'
        return " ", html.Img(
                src="./assets/star.png",
                id="imageNeigColor"
            ),html.Img(
                src="./assets/classes.png",
                id="imageClassColor"
            ),html.Img(
                src="./assets/edgesClicked.png",
                id="imageEdgeColor"
            )
    else:
        return " ", html.Img(
            src="./assets/star.png",
            id="imageNeigColor"
        ),html.Img(
            src="./assets/classes.png",
            id="imageClassColor"
        ),html.Img(
            src="./assets/edges.png",
            id="imageEdgeColor"
        ) 

def update_shortest_path(G,input1,input2,style):
    if (has_path(G,int(input1), int(input2))):
        path = shortest_path(G,int(input1), int(input2))
        edges_list = []
        nodes_list = []
        length = len(path)
        newStyle = []
        for i in range (length-1):
            nodes_list.append(path[i])
            nodes_list.append(path[i+1])
            edges_list.append((path[i],path[i+1]))

        for edge in edges_list:
            for ed in G.edges.data():
                check1 = edge[0] == ed[0]
                check2 = edge[1] == ed[1]
                check3 = edge[0] == ed[1]
                check4 = edge[1] == ed[0]
                if (check1 and check2) or (check3 and check4):
                    newStyle.append({
                        "selector": '#{}'.format(ed[2]['#BioGRID Interaction ID']),
                        "style": {
                            "line-color": "blue",
                            'opacity': 0.9,
                            'z-index': 5000
                        }
                    })

        return "The length of the path between source " + str(input1) + " and target " + str(input2) + " is "+str(length-1),default_stylesheet+newStyle,default_stylesheet+newStyle,default_stylesheet+newStyle,default_stylesheet+newStyle,default_stylesheet+newStyle,el,el,el,el,el,False,style,style
    else:
        return "No path between source and target",default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el,False,style,style

global nodeShortestPath
nodeShortestPath = None

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
              Output('submit-button-state','hidden'),
              Output('input_node_1','style'),
              Output('input_node_2','style'),
              Input('dropdown-metrics', 'value'),
              Input('reloadData', 'n_clicks'),
              Input('gridCytoscape', 'tapNode'),
              Input('concentricCytoscape', 'tapNode'),
              Input('circleCytoscape', 'tapNode'),
              Input('breadthfirstCytoscape', 'tapNode'),
              Input('coseCytoscape', 'tapNode'),
              Input('gridCytoscape', 'tapEdge'),
              Input('concentricCytoscape', 'tapEdge'),
              Input('circleCytoscape', 'tapEdge'),
              Input('breadthfirstCytoscape', 'tapEdge'),
              Input('coseCytoscape', 'tapEdge'),
              Input('classColor', 'n_clicks'),
              Input('neigColor', 'n_clicks'),
              Input('edgeColor', 'n_clicks'),
              State('dropdownClasses', 'value'),
              State('dropdownEdges', 'value'),
              State('gridCytoscape', 'stylesheet'),
              State('concentricCytoscape', 'stylesheet'),
              State('circleCytoscape', 'stylesheet'),
              State('breadthfirstCytoscape', 'stylesheet'),
              State('coseCytoscape', 'stylesheet'),
              Input('submit-button-state','n_clicks'),
              State("input_node_1", "value"),
              State("input_node_2", "value"))
def update_metrics(metric, n_clicks, node1, node2, node3, node4, node5, edge1, edge2, edge3, edge4, edge5, classColor, neigColor, edgeColor, dropdownNodes, dropdownEdges, gridStylesheet, concentricStylesheet,circleStylesheet,breadthfirstStylesheet,coseStylesheet, submitClicks, input1, input2):
    style = {'marginRight':'10px'}
    ctx = dash.callback_context
    inputType = ctx.triggered[0]['prop_id'].split('.')[0]
    global currentGraph, nodeShortestPath
    ############################################################################
    
    # Metrics calculation

    ############################################################################

    if inputType == "dropdown-metrics":
        if (metric == 'betweenness centrality'):
            return "Select a node",default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el,False,style,{'display': 'none'}
        if (metric == 'clustering coefficient'):
            return "The average clustering coefficient of the graph is " + str(clustering_coefficient(currentGraph)) + ".",default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el,True,{'display': 'none'},{'display': 'none'}
        if (metric == 'minimum spanning tree'):
            edges = minimum_spanning_tree(currentGraph)
            edges = np.array(list(edges))[:,:2]
            newStyle = []
            for edge in edges:
                for ed in currentGraph.edges.data():
                    check1 = edge[0] == ed[0]
                    check2 = edge[1] == ed[1]
                    check3 = edge[0] == ed[1]
                    check4 = edge[1] == ed[0]
                    if (check1 and check2) or (check3 and check4):
                        newStyle.append({
                            "selector": '#{}'.format(ed[2]['#BioGRID Interaction ID']),
                            "style": {
                                "line-color": "blue",
                                'opacity': 0.9,
                                'z-index': 5000
                            }
                        })
            return "Minimum spanning tree is shown on the graph", default_stylesheet+newStyle,default_stylesheet+newStyle,default_stylesheet+newStyle,default_stylesheet+newStyle,default_stylesheet+newStyle, el, el, el, el, el,True,{'display': 'none'},{'display': 'none'}

        if (metric == 'shortest path'):
            return "Choose nodes",default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el,False,style,style
   
        if (metric == 'community'):
            comm = community(currentGraph)
            return "The graph is composed with " + str(len(comm)) + " communities.", default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, el, el, el, el, el,False,style,style
        else :
            return "Choose the proprety", default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, el, el, el, el, el, True,{'display': 'none'},{'display': 'none'}
    ############################################################################
    
    # Submit shortest path

    ############################################################################
    elif (ctx.triggered[0]['prop_id'].split('.')[0] == 'submit-button-state'):
        return update_shortest_path(currentGraph,input1,input2,style)
    ############################################################################
    
    # Reload Data

    ############################################################################
    elif inputType == "reloadData":
        global newNodesData, newEdgesData
        if len(newNodesData) and len(newEdgesData):
            newGraph = networkGraph(1000, newNodesData, newEdgesData)
            initialisation(newGraph)
            currentGraph = newGraph
        return "Choose the proprety", default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, el, el, el, el, el, True,{'display': 'none'},{'display': 'none'}
    ############################################################################
    
    # Filters processing

    ############################################################################
    elif inputType == 'neigColor':
        return "Choose the proprety", default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el,True,{'display': 'none'},{'display': 'none'}

    elif inputType in ['gridCytoscape','concentricCytoscape','circleCytoscape','breadthfirstCytoscape','coseCytoscape']:
        checkNode = {
            'gridCytoscape': node1,
            'concentricCytoscape': node2,
            'circleCytoscape': node3,
            'breadthfirstCytoscape': node4,
            'coseCytoscape': node5
        }
        if filterSelected == 'neigColor':
            newStyle = []
            finalNode = checkNode[inputType]
            if finalNode:
                for edge in finalNode["edgesData"]:
                    if edge['source'] == finalNode['data']['id']:
                        newStyle.append({
                            "selector": 'node[id = "{}"]'.format(edge['target']),
                            "style": {
                                'background-color': "blue",
                                'opacity': 0.9
                            }
                        })
                        newStyle.append({
                            "selector": 'edge[id= "{}"]'.format(edge['id']),
                            "style": {
                                "line-color": "blue",
                                'opacity': 0.9,
                                'z-index': 5000
                            }
                        })
                    if edge['target'] == finalNode['data']['id']:
                        newStyle.append({
                            "selector": 'node[id = "{}"]'.format(edge['source']),
                            "style": {
                                'background-color': "blue",
                                'opacity': 0.9
                            }
                        })
                        newStyle.append({
                            "selector": 'edge[id= "{}"]'.format(edge['id']),
                            "style": {
                                "line-color": "blue",
                                'opacity': 0.9,
                                'z-index': 5000
                            }
                        })
                return "Choose the proprety", default_stylesheet+newStyle,default_stylesheet+newStyle,default_stylesheet+newStyle,default_stylesheet+newStyle,default_stylesheet+newStyle, el, el, el, el, el,True,{'display': 'none'},{'display': 'none'}
            else:
                return "Choose the proprety",gridStylesheet,concentricStylesheet,circleStylesheet,breadthfirstStylesheet,coseStylesheet,el,el,el,el,el,True,{'display': 'none'},{'display': 'none'}
        else:
            finalNode = checkNode[inputType]
            if (metric == 'betweenness centrality'):
                dico = betweenness_centrality(currentGraph)
                node_id = finalNode['data']['id']
                beet_centr = dico[int(node_id)]
                return "The betweenness centrality of the node "+ str(node_id) + " is "+ str(beet_centr),default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el,False,style,{'display': 'none'}
            elif (metric == 'shortest path'):
                if finalNode:
                    if nodeShortestPath is None:
                        nodeShortestPath = finalNode
                        return "Select another node",default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el,False,style,style
                    else :
                        input1 = nodeShortestPath['data']['id']
                        input2 = finalNode['data']['id']
                        nodeShortestPath = None
                        return update_shortest_path(currentGraph,input1,input2,style)  
                else:
                   return "Select nodes",default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el,False,style,style 
            else: 
                return "HELLO",default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el,True,{'display': 'none'},{'display': 'none'} 
    elif inputType == 'classColor':
        if filterSelected is None:
           return "Choose the proprety", default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, default_stylesheet, el, el, el, el, el ,True,{'display': 'none'},{'display': 'none'} 
        newStyle = []
        categoryColor = {
            '-': "#808080",
            'Trancription':'blue',
            'Signal Transduction': 'magenta',
            'Cytoskeleton/Motility': 'green',
            'Chromatin Modification':'yellow',
            'DNA Damage Response': 'purple',
            'Other':'black',
            'Ubiquitin Proteasome System': 'orange',
            'Transporter': 'cyan',
            'Tumor suppressor': 'yellow',
            'Oncogene':'blue',
            'Cancer driver': 'green'
        }
        for node in currentGraph.nodes.data():
            if dropdownNodes == None:
                return "Choose the proprety", default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, default_stylesheet+newStyle, el, el, el, el, el,True,{'display': 'none'},{'display': 'none'} 
            if node[1] == {}:
                continue
            if node[1][dropdownNodes] in categoryColor:
                newStyle.append({
                    "selector": 'node[id = "{}"]'.format(node[0]),
                    "style": {
                        'background-color': categoryColor[node[1][dropdownNodes]],
                        'opacity': 0.9
                    }
                })
            else:
                categoryColor[node[1][dropdownNodes]] = ("#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]))
                newStyle.append({
                    "selector": 'node[id = "{}"]'.format(node[0]),
                    "style": {
                        'background-color': categoryColor[node[1][dropdownNodes]],
                        'opacity': 0.9
                    }
                })
        return "Choose the proprety",gridStylesheet+newStyle,concentricStylesheet+newStyle,circleStylesheet+newStyle,breadthfirstStylesheet+newStyle,coseStylesheet+newStyle, el, el, el, el, el,True,{'display': 'none'},{'display': 'none'} 
    elif inputType == 'edgeColor':
        if filterSelected is None:
           return "Choose the proprety",gridStylesheet,concentricStylesheet,circleStylesheet,breadthfirstStylesheet,coseStylesheet, el, el, el, el, el,True,{'display': 'none'},{'display': 'none'} 
        newStyle = []
        categoryColor = {
            '-': "#808080",
            'Trancription':'blue',
            'Signal Transduction': 'magenta',
            'Cytoskeleton/Motility': 'green',
            'Chromatin Modification':'yellow',
            'DNA Damage Response': 'purple',
            'Other':'black',
            'Ubiquitin Proteasome System': 'orange',
            'Transporter': 'cyan',
            'Tumor suppressor': 'yellow',
            'Oncogene':'blue',
            'Cancer driver': 'green'
        }
        for edge in currentGraph.edges.data():
            if dropdownEdges == None:
                return "Choose the proprety", gridStylesheet+newStyle,concentricStylesheet+newStyle,circleStylesheet+newStyle,breadthfirstStylesheet+newStyle,coseStylesheet+newStyle, el, el, el, el, el,True,{'display': 'none'},{'display': 'none'} 
            if edge[2] == {}:
                continue
            if edge[2][dropdownEdges] in categoryColor:
                newStyle.append({
                    "selector": '#{}'.format(edge[2]['#BioGRID Interaction ID']),
                    "style": {
                        'line-color': categoryColor[edge[2][dropdownEdges]]
                    }
                })
            else:
                categoryColor[edge[2][dropdownEdges]] = ("#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)]))
                newStyle.append({
                    "selector": '#{}'.format(edge[2]['#BioGRID Interaction ID']),
                    "style": {
                        'line-color': categoryColor[edge[2][dropdownEdges]]
                    }
                })
        return "Choose the proprety", gridStylesheet+newStyle,concentricStylesheet+newStyle,circleStylesheet+newStyle,breadthfirstStylesheet+newStyle,coseStylesheet+newStyle, el, el, el, el, el,True,{'display': 'none'},{'display': 'none'} 
          
    else:
        return "Choose the proprety",default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,default_stylesheet,el,el,el,el,el,True,{'display': 'none'},{'display': 'none'} 

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

@app.callback(Output('dataInformationText', 'children'),
              Input('gridCytoscape', 'tapNode'),
              Input('concentricCytoscape', 'tapNode'),
              Input('circleCytoscape', 'tapNode'),
              Input('breadthfirstCytoscape', 'tapNode'),
              Input('coseCytoscape', 'tapNode'),
              Input('gridCytoscape', 'tapEdge'),
              Input('concentricCytoscape', 'tapEdge'),
              Input('circleCytoscape', 'tapEdge'),
              Input('breadthfirstCytoscape', 'tapEdge'),
              Input('coseCytoscape', 'tapEdge'))
def displayTapNodeData(node1,node2,node3,node4,node5,edge1,edge2,edge3,edge4,edge5):
    ctx = dash.callback_context
    inputType = ctx.triggered[0]['prop_id']

    dictNode = {
        "gridCytoscape" : node1,
        "concentricCytoscape" : node2,
        "circleCytoscape" : node3,
        "breadthfirstCytoscape" : node4,
        "coseCytoscape" : node5
    }
    dictEdge = {
        "gridCytoscape" : edge1,
        "concentricCytoscape" : edge2,
        "circleCytoscape" : edge3,
        "breadthfirstCytoscape" : edge4,
        "coseCytoscape" : edge5
    }

    if inputType.split(".")[1] == "tapNode":
        for node in currentGraph.nodes.data():
            if int(node[0]) == int(dictNode[inputType.split(".")[0]]['data']['id']):
                if node[1] == {}:
                    return html.Table(
                        [html.Tr([html.Th("Attribute"), html.Th("Gene data")])] + 
                        [html.Tr([
                            html.Td(
                                "Id"
                            ),
                            html.Td(
                                dictNode[inputType.split(".")[0]]['data']['id']
                            )
                        ]),
                        html.Tr([
                            html.Td(
                                "Not available"
                            ),
                            html.Td(
                                "Not available"
                            )
                        ])]
                    )
                else:
                    return html.Table(
                        # Header
                        [html.Tr([html.Th("Attribute"), html.Th("Gene data")])] + 
                        # Body
                        [html.Tr([
                            html.Td(
                                i
                            ),
                            html.Td(
                                node[1][i]
                            )
                        ]) for i in node[1]]
                    )
    elif inputType.split(".")[1] == "tapEdge":
        for edge in currentGraph.edges.data():
            check1 = (int(edge[0]) == int(dictEdge[inputType.split(".")[0]]['data']['source']))
            check2 = (int(edge[1]) == int(dictEdge[inputType.split(".")[0]]['data']['target']))
            check3 = (int(edge[0]) == int(dictEdge[inputType.split(".")[0]]['data']['target']))
            check4 = (int(edge[1]) == int(dictEdge[inputType.split(".")[0]]['data']['source']))
            if (check1 and check2) or (check3 and check4):
                return html.Table(
                        # Header
                        [html.Tr([html.Th("Attribute"), html.Th("Gene data")])] + 
                        # Body
                        [html.Tr([
                            html.Td(
                                i
                            ),
                            html.Td(
                                edge[2][i]
                            )
                        ]) for i in edge[2]]
                    )
        return html.Table(
            [html.Tr([html.Th("Attribute"), html.Th("Information")])]
        )
    else:
        return html.Table(
            [html.Tr([html.Th("Attribute"), html.Th("Information")])]
        )

app.run_server(debug=True)