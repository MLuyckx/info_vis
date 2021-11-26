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

app = dash.Dash(__name__)

genes = pd.read_csv('./BIOGRID-PROJECT-glioblastoma_project-GENES.projectindex.txt', sep="\\t")
interactions = pd.read_csv('./BIOGRID-PROJECT-glioblastoma_project-INTERACTIONS.tab3.txt', sep="\\t")

interactions = interactions.values
genes = genes.values

tenth = int(len(interactions)/1)
tenthInteractions = interactions[:tenth]

el = []
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

def loadData():
    for x in genes:
        genesId.append(x[0])
        el.append({
            'data': {'id': str(x[0]), 'label': str(x[3]), 'size': x[7]},
            'classes': 'known'
        })
        """
        if (x[18] == "Oncogene"):
            genesOncogene.append(x)
        elif (x[18] == "Tumor Suppressor"):
            genesTumorSuppressor.append(x)
        elif (x[18] == "Cancer Driver"):
            genesCancerDriver.append(x)
        """

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

cyGrid = Cytoscape.newCyto('grid', el)
cyConcentric = Cytoscape.newCyto('concentric', el)
cyCircle = Cytoscape.newCyto('circle', el)
cyBreadthfirst = Cytoscape.newCyto('breadthfirst', el)
cyCose = Cytoscape.newCyto('cose', el)

global layoutActivated, addPressed, removePressed
layoutActivated = 1
addPressed = 0
removePressed = 0

def createBasicLayout():
    app.layout = html.Div([
        html.Div(className = "headerDiv",children=[
            html.H2("PROJECT : INFORMATION VISUALISATION"),
            html.P(id='dd-output-container'),
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
                    value = ["CSE"],
                    labelStyle = {'display': 'block'}
                ),
                html.Button('Add layout', id="addLayout", n_clicks=0),
                html.Button('Remove layout', id="removeLayout", n_clicks=0)
            ]),
            html.Div(className = "metrics", children = [
                html.P("Metrics there.")
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

loadData()
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

@app.callback(  Output("graph1", 'children'),
                Input('dropdown1', 'value'))
def update_output(value):
    layouts = {"GRD": cyGrid, "CRC": cyCircle, "CCT": cyConcentric, "BFT": cyBreadthfirst, "CSE": cyCose}
    return layouts[value]
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