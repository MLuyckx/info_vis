import dash_cytoscape as cyto

class Cytoscape:
    def newCyto(type, elements, id):
        myCyto = cyto.Cytoscape(
            id=id,
            className='cytoscape-event-callbacks-1',
            elements=elements,
            layout={
                'name': type,
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
                        "width": "mapData(size, 0, 1000, 3, 30)",
                        "height": "mapData(size, 0, 1000, 3, 30)",
                        "opacity":0.9,
                        "color":"red"
                    }
                },
                {
                    'selector': 'edge',
                    'style':{
                        "line-color": "red",
                        "width": 0.5
                    }
                }
            ],
            responsive=True
        )
        return myCyto