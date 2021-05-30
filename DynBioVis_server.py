    #### IMPORTING
    
import networkx as nx
import numpy as np
import pandas as pd

from datetime import datetime

from sqlalchemy import create_engine
import sqlite3
import pickle, json

from tqdm import tqdm

from random import sample, seed

from bokeh.plotting import curdoc, figure, show, from_networkx
from bokeh.layouts import column, row
from bokeh.palettes import Spectral4, Turbo256
from bokeh.models.widgets import HTMLTemplateFormatter

from bokeh.models import (
    CustomJS, 
    Slider, 
    RangeSlider, 
    MultiLine, 
    Toggle,
    EdgesAndLinkedNodes, 
    NodesAndLinkedEdges,
    GraphRenderer, 
    Ellipse, 
    graphs, 
    LabelSet, 
    Circle, 
    Range1d,
    CheckboxGroup,
    DataTable,
    ColumnDataSource,
    TableColumn,
    TextInput,
    MultiChoice,
    Button,
    Div,
    Legend,
    HoverTool,
)

#seed(40)

    #### Opening databases and files

with open('DynBioVis_config.json', 'r') as f:
    tv_paths = json.load(f)

with open(
    tv_paths['networkx_network_path'], 'rb') as f:
    G_nx_test = pickle.load(f)
    
    
db_connection = sqlite3.connect(
    tv_paths['sql_db_path']
)
db_cursor = db_connection.cursor()


all_semTypes = set()
with open('SemanticTypes_2018AB.txt') as f:
    for line in f:
        all_semTypes.add(line.split('|')[-1].strip())

semType_palette_dict = dict(
    zip(
        list(all_semTypes), sample(Turbo256, len(all_semTypes))
    )
)


sql_query = \
    """
    SELECT source, records.pmid, timestamp, titles.title FROM records
    LEFT JOIN titles
    ON records.pmid = titles.pmid
    WHERE pair_text = "{} -> {}";
    """




    #### Supplimentary functions (Needed for callbacks)

def Get_nx_subgraph(
    nx_graph,
    date_start,
    date_end,
    nodename, # we allow partial nodename too
    radius,
    allowed_semTypes,
) -> nx.Graph:
    """Creates an ego_graph based on:
        * a given timeframe
        * central node
        * radius
    """
    
    timeframed_edgelist = [
        (u,v,d) for u,v,d in nx_graph.edges(data = True) if
            (d['date_added'] > date_start) and \
            (d['date_added'] <= date_end)
            #(nx_graph.nodes[u]['semType'].intersection(allowed_semTypes)) and \
            #(nx_graph.nodes[v]['semType'].intersection(allowed_semTypes))
    ]
    
    nx_graph_timeframed = nx.from_edgelist(
        timeframed_edgelist
    )
    
    target_nodes = []
    for single_nodename in nodename.split(';'):
        target_nodes += [n for n in nx_graph_timeframed.nodes if single_nodename in n.lower()]
    
    target_nodes_set = set(target_nodes)
    
    timeframed_edgelist_semTypes_filtered = [
        (u,v,d) for u,v,d in nx_graph_timeframed.edges(data = True) if
            (
                (nx_graph.nodes[u]['semType'].intersection(allowed_semTypes)) and \
                (nx_graph.nodes[v]['semType'].intersection(allowed_semTypes))
            ) 
            or \
            (
                (u in target_nodes_set) and \
                (nx_graph.nodes[v]['semType'].intersection(allowed_semTypes))
            )
            or \
            (
                (v in target_nodes_set) and \
                (nx_graph.nodes[u]['semType'].intersection(allowed_semTypes))
            )
    ]
    
    nx_graph_timeframed = nx.from_edgelist(
        timeframed_edgelist_semTypes_filtered
    )
    
    new_target_nodes = set(nx_graph_timeframed.nodes).intersection(target_nodes_set)
    
    #target_nodes = []
    #for single_nodename in nodename.split(';'):
    #    target_nodes += [n for n in nx_graph_timeframed.nodes if single_nodename in n.lower()]
    
    nx_graph_to_return = nx.Graph()
    for node in new_target_nodes:
        nx_graph_to_return.add_edges_from(
            nx.ego_graph(
                nx_graph_timeframed, 
                node, 
                radius
            ).edges()
        )
    
    for node in nx_graph_to_return.nodes:
        nx_graph_to_return.nodes[node]['semType'] = \
            list(nx_graph.nodes[node]['semType'])[0]
    
    return nx_graph_to_return, list(new_target_nodes)


def Get_edge_info(
    u,
    v,
) -> None:
    
    """Runs SQL query, which finds temporal information about a given edge
    Returns pandas df
    """
    
    uv_df = \
        pd.read_sql_query(
            sql_query.format(u, v),
            db_connection
        )
    
    vu_df = \
        pd.read_sql_query(
            sql_query.format(v, u),
        db_connection)
    
    uv_vu_df = pd.concat([uv_df, vu_df]).drop_duplicates()
    
    uv_vu_df['pair_text'] = f'{u} <-> {v}'
    
    return uv_vu_df.sort_values(by='timestamp')


def Get_all_edges_info() -> None:
    """Runs function Get_edge_info for all edges in the graph and writes it to disk."""
    
    df_list = []
    
    for (u,v) in tqdm(nx_graph.edges):
        df_list.append(Get_edge_info(u,v))

    all_edges_df = pd.concat(df_list).set_index('pair_text')
    
    now = datetime.now()
    dt_string = now.strftime("%m_%d_%Y-%H_%M")
    
    current_fname = ';'.join([cn.replace(' ', '_') for cn in centric_nodes]) + '_' + dt_string
    
    all_edges_df.to_csv(f'{current_fname}.csv')
    
    print(f'All edges info saved to: {current_fname}')
    
    return None
    
    
    
    
    #### Plotting part

tsl = '2018-01-01'
tsu = '2020-01-01'

nx_graph, centric_nodes = \
    Get_nx_subgraph(
        nx_graph=G_nx_test,
        date_start=tsl,
        date_end=tsu,
        nodename='lzheime',
        radius=2,
        allowed_semTypes=set(['Pathologic Function', 'Disease or Syndrome'])
    )




    #### Global variables

plot = figure(
    #plot_width=2000, 
    plot_height=800,
    x_range=Range1d(-1.1, 1.1), y_range=Range1d(-1.1, 1.1),
    tools='wheel_zoom, pan, tap',
    active_scroll='wheel_zoom',
    width_policy='max',
    #height_policy='max',
)

edge_history_CDS = ColumnDataSource(data=dict())

history_columns = [
    TableColumn(field=colname, title=colname, width=100) for colname in [
        'source', 'pair_text','pmid','timestamp','title'
    ]
]
history_columns[-1].width = 700

edge_history_datatable = DataTable(
    source=edge_history_CDS,
    columns=history_columns,
    width_policy='max',
    height=300,
)
    
allowed_semTypes = all_semTypes

hover_tooltips = [
    ("Semantic Type", "@semType"),
    ("Node Degree", "@node_degree"),
    ("Name", "@index")
]




    #### Updading part
    

plot.title.text = f"{centric_nodes} network (time period: from {tsl} to {tsu})"

nx_layout_dict = nx.spring_layout(nx_graph, seed=42)

## FAKE GLYPHS FOR LEGEND
current_colors_and_semTypes = []

for node in nx_graph.nodes:

    current_semType = nx_graph.nodes[node]['semType']
    current_color = semType_palette_dict[current_semType]

    nx_graph.nodes[node]['node_color'] = \
        semType_palette_dict[nx_graph.nodes[node]['semType']]

    current_colors_and_semTypes.append(
        (current_color, current_semType)
    )

    nx_graph.nodes[node]['node_size'] = \
        15 + 30*(nx_graph.degree[node]/30)
    
    nx_graph.nodes[node]['node_degree'] = \
        nx_graph.degree[node]

graph_renderer = from_networkx(nx_graph, nx_layout_dict, scale=1, center=(0, 0))
graph_renderer.node_renderer.data_source.data['x'] = np.array(list(nx_layout_dict.values()))[:,0]
graph_renderer.node_renderer.data_source.data['y'] = np.array(list(nx_layout_dict.values()))[:,1]

graph_renderer.edge_renderer.glyph = MultiLine(line_color="#CCCCCC", line_alpha=0.8, line_width=5)
graph_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
#graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)

graph_renderer.selection_policy = EdgesAndLinkedNodes()
#graph_renderer.inspection_policy = EdgesAndLinkedNodes()



## Adding "fake" legend data
current_colors_and_semTypes = np.array(list(set(current_colors_and_semTypes)))

legend_CDS = (
    dict(
        x=[0]*len(current_colors_and_semTypes),
        y=[0]*len(current_colors_and_semTypes),
        color=current_colors_and_semTypes[:,0],
        label=current_colors_and_semTypes[:,1],
    )
)
circle_renderer = plot.circle(
    x='x', y='y', radius=0.0, color='color', legend_field='label', source=legend_CDS)

labels = LabelSet(
    x='x', y='y',
    text='index',
    x_offset=5, y_offset=5,
    source = graph_renderer.node_renderer.data_source , 
    render_mode='canvas'
)

graph_renderer.node_renderer.glyph = Circle(size='node_size', fill_color='node_color')
#graph_renderer.edge_renderer.glyph = MultiLine(line_color="black", line_alpha=0.8, line_width=1)
plot.renderers.append(graph_renderer)
plot.renderers.append(labels)

labels.visible = False

node_hover_tool = HoverTool(
    tooltips=hover_tooltips
)
plot.add_tools(node_hover_tool)




    #### Data Table update

def CALLBACK_Click_on_edge(attr, old, new) -> None:
    """Callback function, activated when click on edge occurred."""
    #print(attr, old, new)
    
    if len(new) == 1:
        clicked_idx = new[0]
        u = graph_renderer.edge_renderer.data_source.data['start'][clicked_idx]
        v = graph_renderer.edge_renderer.data_source.data['end'][clicked_idx]
        print(f'clicked_on_edge: {u} <-> {v}')
        
        clicked_edge_hist_df = Get_edge_info(u,v)
        
        #print(clicked_edge_hist_df)
    
        ## Datatable updating part
        edge_history_CDS.data = clicked_edge_hist_df
        
    return None




    #### Controls column
    
checkbox = CheckboxGroup(labels=["Show labels"])
checkbox.js_on_click(CustomJS(args=dict(labels=labels, checkbox=checkbox), code="""labels.visible = 0 in checkbox.active;"""))
    
text_input = TextInput(value="", title="Term Search Box:")
text_input.js_on_change(
    "value", 
    CustomJS(
        code="""console.log('text_input: value=' + this.value, this.toString())"""
    )
)


text_input_from = TextInput(value="2018-01-01", title="Timeframe (lower bound):")
text_input_from.js_on_change(
    "value", 
    CustomJS(
        code="""console.log('text_input: value=' + this.value, this.toString())"""
    )
)


text_input_to = TextInput(value="2020-01-01", title="Timeframe (upper bound):")
text_input_to.js_on_change(
    "value", 
    CustomJS(
        code="""console.log('text_input: value=' + this.value, this.toString())"""
    )
)

LABELS = list(all_semTypes)
multi_choice = MultiChoice(
    value=[], 
    options=LABELS, 
    title='Allowed Semantic Types: ',
    #height=200,
    width = 285,
    #height_policy='fixed',
    #css_classes=['Scrollable'],
)
multi_choice.js_on_change("value", CustomJS(code="""
    console.log('multi_choice: value=' + this.value, this.toString())
"""))


slider_locality = Slider(
    start=1, 
    end=5, 
    value=1, 
    step=1, 
    title="Neighbors order (locality parameter)")
slider_locality.js_on_change("value", CustomJS(code="""
    console.log('slider: value=' + this.value, this.toString())
"""))


button_estimate = Button(label="Estimate Graph Size", button_type="primary")
string_estimate = """<b>Graph Size</b>:
    <p>Number of nodes: {} </p>
    <p>Number of edges: {} </p>"""

div_estimate = Div(
    text=string_estimate.format(
        nx_graph.number_of_nodes(),
        nx_graph.number_of_edges(),
    ),
)

button_search = Button(label="Reconstruct Graph", button_type="success")
#button_search.on_click(Search_handler)

button_save_info = Button(label="Save Temporal Edge Info", button_type="primary")




    #### Graph redraw callback
    
def Reconstruct_graph() -> None:
    """Rewrites global variable nx_graph. Does NOT redraw anything.
    """
    
    global nx_graph
    global centric_nodes
    
    global tsl
    global tsu
    
    search_term = text_input.value.lower()
    tsl = text_input_from.value
    tsu = text_input_to.value
    
    locality_radius = slider_locality.value
    
    print(f'Search term: {search_term}')
    
    if len(search_term) > 2: 
        if multi_choice.value:
            current_multi_choice = set(multi_choice.value)
        else:
            current_multi_choice = set(LABELS)

        nx_graph, centric_nodes = \
            Get_nx_subgraph(
                nx_graph=G_nx_test,
                date_start=tsl,
                date_end=tsu,
                nodename=search_term,
                radius=locality_radius,
                allowed_semTypes=current_multi_choice,
            )
    
    

def CALLBACK_Redraw_plot() -> None:
    """Updates graph when 'Reconstruct Graph' button is pressed.
    REDRAWS the dashboard.
    """
    
    #Reconstruct_graph()
    
    #global circle_renderer
    
    print('Edges:', nx_graph.number_of_edges())
    
    #plot.renderers.remove(graph_renderer)
    #plot.renderers.remove(labels)
    
    plot.title.text = f" {centric_nodes} network (time period: from {tsl} to {tsu})"
    #plot.legend.clear()
    
    nx_layout_dict_upd = nx.spring_layout(nx_graph, seed=42)
    
    ## Adding node attributes: color and size
    
    current_colors_and_semTypes = []
    
    for node in nx_graph.nodes:
        
        current_semType = nx_graph.nodes[node]['semType']
        current_color = semType_palette_dict[current_semType]
        
        nx_graph.nodes[node]['node_color'] = \
            semType_palette_dict[nx_graph.nodes[node]['semType']]
        
        current_colors_and_semTypes.append(
            (current_color, current_semType)
        )
        
        nx_graph.nodes[node]['node_size'] = \
        15 + 30*(nx_graph.degree[node]/30)
    
        nx_graph.nodes[node]['node_degree'] = \
            nx_graph.degree[node]
        
    ## Adding "fake" legend data
    current_colors_and_semTypes = np.array(list(set(current_colors_and_semTypes)))
    
    legend_CDS = (
        dict(
            x=[0]*len(current_colors_and_semTypes),
            y=[0]*len(current_colors_and_semTypes),
            color=current_colors_and_semTypes[:,0],
            label=current_colors_and_semTypes[:,1],
        )
    )
    
    #if circle_renderer in plot.renderers:
     #   plot.renderers.remove(circle_renderer)
    
    circle_renderer.data_source.data = legend_CDS
    
    #circle_renderer = \
     #   plot.circle(x='x', y='y', radius=0.0, color='color', legend_group='label', source=legend_CDS)
    
    print('legend_length:', len(plot.legend))
    
    graph_renderer_upd = from_networkx(
        nx_graph, 
        nx_layout_dict_upd, 
        scale=1, center=(0, 0))
    
    graph_renderer_upd.node_renderer.data_source.data['x'] = np.array(
        list(nx_layout_dict_upd.values()))[:,0]
    
    graph_renderer_upd.node_renderer.data_source.data['y'] = np.array(
        list(nx_layout_dict_upd.values()))[:,1]

    graph_renderer_upd.edge_renderer.glyph = MultiLine(line_color="#CCCCCC", line_alpha=0.8, line_width=5)
    graph_renderer_upd.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
    #graph_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)

    graph_renderer_upd.selection_policy = EdgesAndLinkedNodes()
    #graph_renderer.inspection_policy = EdgesAndLinkedNodes()

    labels_upd = LabelSet(
        x='x', y='y',
        text='index',
        x_offset=5, y_offset=5,
        source = graph_renderer_upd.node_renderer.data_source , 
        render_mode='canvas'
    )
    
    labels_upd.visible = False

    graph_renderer_upd.node_renderer.glyph = Circle(size='node_size', fill_color='node_color')
    #graph_renderer.edge_renderer.glyph = MultiLine(line_color="black", line_alpha=0.8, line_width=1)
    
    global graph_renderer
    global labels
    
    plot.renderers.remove(graph_renderer)
    plot.renderers.remove(labels)
    
    graph_renderer = graph_renderer_upd
    labels = labels_upd
    
    plot.renderers.append(graph_renderer)
    plot.renderers.append(labels)
    
    graph_renderer.edge_renderer.data_source.selected.on_change(
        "indices", 
        CALLBACK_Click_on_edge,
    )
    checkbox.js_on_click(CustomJS(args=dict(labels=labels, checkbox=checkbox), code="""labels.visible = 0 in checkbox.active;"""))
    checkbox.active = []
    
    plot.tools.pop()
    
    node_hover_tool = HoverTool(
        tooltips=hover_tooltips
    )
    plot.add_tools(node_hover_tool)

    return None

def CALLBACK_Estimate_graph_size() -> None:
    """Reconstructs graph and estimates its size. Does NOT update render.
    """
    Reconstruct_graph()
    string_current_graph_size = string_estimate.format(
        nx_graph.number_of_nodes(),
        nx_graph.number_of_edges()
    )
    div_estimate.text = string_current_graph_size
    print(string_current_graph_size)

    
    

    #### Connecting callbacks
    
graph_renderer.edge_renderer.data_source.selected.on_change(
    "indices", 
    CALLBACK_Click_on_edge,
)

button_search.on_click(
    CALLBACK_Redraw_plot,
)

button_estimate.on_click(
    CALLBACK_Estimate_graph_size,
)

button_save_info.on_click(
    Get_all_edges_info,
)




    #### Layout
    
left_column = column(
    text_input,
    checkbox,
    text_input_from,
    text_input_to,
    slider_locality,
    button_estimate,
    button_save_info,
    div_estimate,
    button_search,
    multi_choice,
)

right_column = column(
    plot, 
    edge_history_datatable,
    width_policy='max',
)

layout = row(left_column, right_column)




    #### Starting app

curdoc().add_root(layout)