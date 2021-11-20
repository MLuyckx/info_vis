import dash
import json
from dash import dcc
from dash import html
import dash_cytoscape as cyto
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
import numpy as np

genes = pd.read_csv('./BIOGRID-PROJECT-glioblastoma_project-GENES.projectindex.txt', sep="\\t")
interactions = pd.read_csv('./BIOGRID-PROJECT-glioblastoma_project-INTERACTIONS.tab3.txt', sep="\\t")

interactions = interactions.values
genes = genes.values

tenth = int(len(interactions)/50)
tenthInteractions = interactions[:tenth]

el = []
genesId = []
genesData = {}

"""for x in genes:
    #el.append({'data': {'id': str(x[0]), 'label': x[3], 'interactions': x[7]}})
    genesId.append(str(x[0]))
    genesData[str(x[0])] = {'label':x[3], 'interactions':x[7]}
i=0
"""
#print(genesData['108621']['label'])

for x in tenthInteractions:
    if x[3] == x[4]:
        continue
    if x[3] not in genesId:
        el.append({
            'data': {'id': str(x[3]), 'label': '/', 'interactions': 0},
            'position':{'x':0, 'y':0}
        })
        genesId.append(x[3])
    #else:
        #el.append({'data': {'id': str(x[3]), 'label': genesData[str(x[3])]['label']}})
    if x[4] not in genesId:
        el.append({
            'data': {'id': str(x[4]), 'label': '/', 'interactions': 0},
            'position':{'x':0, 'y':0}
        })
        genesId.append(x[4])
    el.append({
        'data': {'source': str(x[3]), 'target': str(x[4])}, 
        'classes': 'red',
        'position': None
    })

"""
for x in interactions:
    el.append({'data': {'source': str(x[3]), 'target': str(x[4])}})
"""

app = dash.Dash(__name__)

styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

cy = cyto.Cytoscape(
    id='cytoscape-event-callbacks-1',
    elements=el,
    layout={
        #available names = ('preset'), ('random'), 'grid', 'circle', 'concentric', 'breadthfirst', 'cose'
        'name': 'cose',
        'animate':False
    },
    style={'width': '100%', 'height': '500px'},
    stylesheet=[
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
)

app.layout = html.Div([
    html.P("Dash Cytoscape:"),
    cy,
    html.Pre(id='cytoscape-tapNodeData-json', style=styles['pre'])
])

method_list = [method for method in dir(cy) if method.startswith('__') is False]
print(method_list)

print(cy.available_properties)

@app.callback(Output('cytoscape-tapNodeData-json', 'children'),
              Input('cytoscape-event-callbacks-1', 'tapNodeData'))
def displayTapNodeData(data):
    return json.dumps(data, indent=2)

app.run_server(debug=True)