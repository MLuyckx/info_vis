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

loadData()

cy = cyto.Cytoscape(
    id='cytoscape-event-callbacks-1',
    elements=el,
    layout={
        #available names = ('preset'), ('random'), 'grid', 'circle', 'concentric', 'breadthfirst', 'cose'
        'name': 'cose',
        'animate':False
    },
    style={
        'width': '100%', 
        'height': '100%',
        'border': '1px solid grey'
    },
    stylesheet=[
        # Class selectors
        {
            'selector': 'node',
            'style': {
                "width": "mapData(size, 0, 1000, 3, 40)",
                "height": "mapData(size, 0, 1000, 3, 40)",
                "opacity":0.9
            }
        },
        {
            'selector': 'edge',
            'style':{
                "line-color": "red",
                "width": 0.5
            }
        },
        {
            'selector': '.known',
            'style':{
                "color": "yellow",
            }
        },
        {
            'selector': '.unknown',
            'style':{
                "color": "blue",
            }
        }
    ]
)

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
    html.Div(className = "metrics", children = [
        html.P("Metrics there.")
    ]),
    html.Div(className = "graphs", children = [
        cy,
    ]),
    html.Div(className = "filters", children = [
        html.P("Filters there.")
    ]),
    
    html.P(id='cytoscape-tapNodeData-json'),
    html.P(id="outputnodes"),
    html.P(id="outputedges")
])

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
        #print(list_of_names)
        print(pd.read_csv(list_of_names).values)


#method_list = [method for method in dir(cy) if method.startswith('__') is False]
#print(method_list)
@app.callback(Output('cytoscape-tapNodeData-json', 'children'),
              Input('cytoscape-event-callbacks-1', 'tapNodeData'))
def displayTapNodeData(data):
    if json.dumps(data, indent=2) != "null":
        print(type(json.dumps(data, indent=2)))
        return json.dumps(data, indent=2)

app.run_server(debug=True)