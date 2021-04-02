#!/usr/bin/env python3


import numpy as np
import pandas as pd


import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots


import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

pio.templates.default = "plotly_white"

def select_sign(x):

    # PC change must be greater than 0
    if x.baseline_logAUC < x.log_auc:

        if x.pc_enrich < 0:
            return x.pc_enrich * -1

        else:
            return x.pc_enrich

    # PC Change must be less than 0
    else:
        if x.pc_enrich > 0:
            return x.pc_enrich * -1
        else:
            return x.pc_enrich

def calculate_enrichment(frame):
    frame['pc_enrich'] = (frame.log_auc - frame.baseline_logAUC) / np.abs(frame.baseline_logAUC)
    frame['pc_enrich'] = frame.apply(
        lambda x : select_sign(x),
        axis = 1
    )

    return frame

def calculate_speedup(frame):
    frame['speedup'] = frame.baseline_time / frame.time
    return frame

def generate_baseline(frame):
    baseline_frame = frame[frame.k == 1].\
        groupby(['receptor', 'match_type']).\
        apply(
            lambda x : pd.Series({
                "baseline_AUC" : x.auc.mean(),
                "baseline_logAUC" : x.log_auc.mean(),
                "baseline_time" : x.time.mean()
            })
        ).\
        reset_index()
    return baseline_frame

def make_box(subframe, x_val, y_val, rec, match_type, v = False, showlegend=False):
    d = {
        "c_match" : "#4A5385",
        "s_match" : "#942C2C"
    }
    match_name = {
        "c_match" : "Normal Match",
        "s_match" : "Scaled Match"
    }

    trace = go.Box(
        x = subframe[x_val],
        y = subframe[y_val],
        marker_color = d[match_type],
        name = match_name[match_type],
        visible = v,
        showlegend=showlegend,
        legendgroup = match_type
    )

    return trace

def load_timescores(fn):
    # Load in time/enrichment data
    time_scores = pd.read_csv(fn, sep="\t")

    # Define K from cluster_id
    time_scores['k'] = time_scores.cluster_id.apply(lambda x : int(x.split("_")[0][1:]))

    # build baseline statistics frame
    baseline_frame = generate_baseline(time_scores)

    # merge time/scores with baseline
    time_scores = time_scores.merge(baseline_frame)

    # calculate logAUC enrichment
    time_scores = calculate_enrichment(time_scores)

    # calculate speedup
    time_scores = calculate_speedup(time_scores)

    return time_scores

def load_msframe(fn):
    ms_frame = pd.read_csv(fn, sep="\t")
    ms_frame['ms_id'] = ms_frame.ms_id.apply(lambda x : "sph.{}".format(x))
    return ms_frame

def load_cooccurrence(fn):
    cooc_frame = pd.read_csv(fn, sep="\t")
    return cooc_frame

app = dash.Dash(__name__)


time_scores = load_timescores("../data/merged_time_and_enrichment.tab")
ms_frame = load_msframe("../data/sphere_usage.tab")
cooc_frame = load_cooccurrence("../data/co-occurrence.tab")

# unique receptors
receptors = time_scores.receptor.unique()

# cluster ids
cluster_ids = sorted(ms_frame.cluster_id.unique())



app.layout = html.Div([
    dcc.Tabs(id='tab-id', value='tab-1', children=[
        dcc.Tab(label='Global Timing and Enrichment Performance', value='tab-1'),
        dcc.Tab(label='Individual Timing and Enrichment Performance', value='tab-2'),
        dcc.Tab(label='Matching Sphere Usage', value='tab-3'),
    ]),
    html.Div(id='tab-content')
])

@app.callback(
    Output('tab-content', 'children'),
    Input('tab-id', 'value')
)
def render_content(tab):

    t1 = html.Div([
        html.Div([
            html.Div([
                dcc.RadioItems(
                    id = "Aggregate",
                    options = [{'label' : m, 'value' : m} for m in ["Aggregate", "All"]],
                    value = "All",
                    style={'float' : 'left', "width" : "100%", 'display' : 'block'}
                )
            ])
        ]),
        dcc.Graph(id = "Global_PercentChange", style = {'width' : '100%', 'display' : 'inline-block'}),
        dcc.Graph(id = "Global_Speedup", style = {'width' : '100%', 'display' : 'inline-block'})
    ])

    # Plot Performance Statistics
    t2 = html.Div([
        html.Div([
            html.Div([
                dcc.Dropdown(
                    id = "Receptor",
                    options = [{'label' : "Receptor : {}".format(r), 'value' : r} for r in receptors],
                    value = "AA2AR",
                    style={'float' : 'left', "width" : "51%", 'display' : 'block'}
                )
            ])
        ]),
        dcc.Graph(id = 'LogAUC', style = {'width' : '50%', 'float' : 'left', 'display' : 'inline-block'}),
        dcc.Graph(id = 'Enrichment', style = {'width' : '50%', 'display' : 'inline-block'}),
        dcc.Graph(id = 'Speedup', style = {'width' : '65%', 'display' : 'inline-block'}),
        dcc.Graph(id = 'Correlation', style = {'width' : '35%', 'display' : 'inline-block'})
    ])

    # Plot Single Run Statistics
    t3 = html.Div([
        html.Div([
            html.Div([
                dcc.Dropdown(
                    id = "Receptor",
                    options = [{'label' : "Receptor : {}".format(r), 'value' : r} for r in receptors],
                    value = "AA2AR",
                    style={'float' : 'left', "width" : "51%", 'display' : 'block'},
                    searchable=True
                ),
                dcc.Dropdown(
                    id = "Cluster ID",
                    options = [{'label' : "ID : {}".format(c), 'value' : c} for c in cluster_ids],
                    value = "k1_0",
                    style={'float' : 'left', "width" : "51%", 'display' : 'block'}
                ),
                dcc.RadioItems(
                    id = "MatchType",
                    options = [{'label' : m, 'value' : m} for m in ['c_match', 's_match', "mt"]],
                    value = "c_match",
                    style={'float' : 'left', "width" : "100%", 'display' : 'block'}
                )
            ])
        ]),
        dcc.Graph(
            id = 'Co-Occurrence',
            style = {
                'width' : '30%', 'float' : 'left',
                'height' : "35vh", 'column' : "1"
                }),
        dcc.Graph(
            id = 'SphereUsage',
            style = {
                'width' : '60%', 'float' : 'right',
                'height' : "60vh", 'column' : "2"
                }),
        dcc.Graph(
            id = 'Ligand-Decoy',
            style = {
                'width' : '40%', 'float' : 'left',
                'height' : "45vh", 'column' : "1"
                }),
        dcc.Graph(
            id = 'UsageBar',
            style = {
                'width' : '60%', 'float' : 'right',
                'height' : "20vh", 'column' : "1"
                })
    ])

    if tab == 'tab-1':
        return t1
    elif tab == 'tab-2':
        return t2
    else:
        return t3


###########################################
# Global Enrichment and Timing Statistics #
###########################################

agg_scores = time_scores.\
    groupby(['receptor','k', 'match_type']).\
    apply(
        lambda x : pd.Series({
            "AggEnrichMean" : x.pc_enrich.mean(),
            "AggEnrichMax" : x.pc_enrich.max(),
            "AggSpeedupMean" : x.speedup.mean(),
            "AggSpeedupMax" : x.speedup.max()
        })
    ).reset_index().\
    melt(id_vars = ['k', 'receptor', 'match_type'])

@app.callback(
    Output("Global_PercentChange", "figure"),
    Input("Aggregate", "value")
)
def Global_PercentChange(agg):
    if agg == "Aggregate":
        fig = px.box(
            agg_scores[agg_scores.variable.str.contains("Enrich")],
            x = 'k', y = 'value', points='all',
            color = 'match_type', hover_name = 'receptor',
            facet_col = 'variable'
        )
    else:
        fig = px.box(
            time_scores, x = 'k', y = 'pc_enrich',
            color = 'receptor', facet_row = 'match_type'
        )

    fig.update_xaxes(title = "Number of Subclusters (K-Means)")
    fig.update_yaxes(title = "%Change LogAUC")
    return fig

@app.callback(
    Output("Global_Speedup", "figure"),
    Input("Aggregate", "value")
)
def Global_PercentChange(agg):
    if agg == "Aggregate":
        fig = px.box(
            agg_scores[agg_scores.variable.str.contains("Speedup")],
            x = 'k', y = 'value', points='all',
            color = 'match_type', hover_name = 'receptor',
            facet_col = 'variable'
        )
    else:
        fig = px.box(
            time_scores, x = 'k', y = 'speedup',
            color = 'receptor', facet_row = 'match_type'
        )
    fig.update_xaxes(title = "Number of Subclusters (K-Means)")
    fig.update_yaxes(title = "Fold Speedup")
    return fig

#################################################
# Individidual Enrichment and Timing Statistics #
#################################################

@app.callback(
    Output("LogAUC", "figure"),
    Input("Receptor", "value")
)
def update_logAUC(rec):
    fig = px.box(
        time_scores[time_scores.receptor == rec],
        x = 'k', y = 'log_auc', points='outliers',
        color = "match_type", hover_name='cluster_id',
        color_discrete_sequence=['#7C80A3', '#B05B67']
    )
    fig.update_xaxes(title = "Number of Subclusters (K-Means)")
    fig.update_yaxes(title = "Adjusted LogAUC")
    return fig

@app.callback(
    Output("Enrichment", "figure"),
    Input("Receptor", "value")
)
def update_Enrichment(rec):
    fig = px.box(
        time_scores[time_scores.receptor == rec],
        x = 'k', y = 'pc_enrich',
        color = "match_type",
        color_discrete_sequence=['#7C80A3', '#B05B67']
    )
    fig.update_xaxes(title = "Number of Subclusters (K-Means)")
    fig.update_yaxes(title = "Percent Change in Adjusted LogAUC")
    return fig

@app.callback(
    Output("Speedup", "figure"),
    Input("Receptor", "value")
)
def update_Speedup(rec):
    fig = px.box(
        time_scores[time_scores.receptor == rec],
        x = 'k', y = 'speedup', points='outliers',
        color = "match_type", hover_name = 'cluster_id',
        color_discrete_sequence=['#7C80A3', '#B05B67']
    )
    fig.update_xaxes(title = "Number of Subclusters (K-Means)")
    fig.update_yaxes(title = "Fold Speedup")
    return fig

@app.callback(
    Output("Correlation", "figure"),
    Input("Receptor", "value")
)
def update_Correlation(rec):
    sub_frame = time_scores[time_scores.receptor == rec]
    plot_frame = sub_frame.\
        groupby(['receptor', 'match_type', 'cluster_id', 'k']).\
        apply(
            lambda x : pd.Series({
                'pc_enrich' : x.pc_enrich.mean(),
                'speedup' : x.speedup.mean()
            })
        ).reset_index()
    plot_frame['k'] = plot_frame['k'].astype(str)

    fig = px.scatter(
        plot_frame,
        x = 'speedup', y = 'pc_enrich',
        color = "k", hover_name = 'cluster_id',
        color_discrete_sequence = px.colors.sequential.Plasma_r,
        symbol = 'match_type'
    )

    fig.update_xaxes(title = "Fold Speedup")
    fig.update_yaxes(title = "Percent Change in Adjusted LogAUC")
    return fig



####################################
# Matching Sphere Usage Statistics #
####################################

@app.callback(
    Output("SphereUsage", "figure"),
    Input("Receptor", "value"),
    Input("MatchType", "value"),
    Input("Cluster ID", "value")
)
def update_sphere_usage(rec, mt, ci):
    sub_frame = ms_frame[
        (ms_frame.match_type == mt) &
        (ms_frame.receptor == rec) &
        (ms_frame.cluster_id == ci)
    ]

    fig = px.scatter_3d(
        sub_frame, x = 'x', y = 'y', z = 'z',
        color = 'Usage', text = 'ms_id',
        color_continuous_scale = "Magma"
    )

    fig.update_layout(title = "Matching Sphere Usage")
    return fig

@app.callback(
    Output("Ligand-Decoy", "figure"),
    Input("Receptor", "value"),
    Input("MatchType", "value"),
    Input("Cluster ID", "value")
)
def update_ligand_usage(rec, mt, ci):
    sub_frame = ms_frame[
        (ms_frame.match_type == mt) &
        (ms_frame.receptor == rec) &
        (ms_frame.cluster_id == ci)
    ]

    fig = px.scatter(
        sub_frame, x = 'Ligand_Usage', y = "Decoy_Usage",
        text = "ms_id"
    )
    fig.update_traces(textposition = 'top center', marker=dict(color = 'black'))
    fig.update_layout(title = "Matching Sphere Ligand/Decoy Usage")

    return fig

@app.callback(
    Output("UsageBar", "figure"),
    Input("Receptor", "value"),
    Input("MatchType", "value"),
    Input("Cluster ID", "value")
)
def update_ligand_usage(rec, mt, ci):
    sub_frame = ms_frame[
        (ms_frame.match_type == mt) &
        (ms_frame.receptor == rec) &
        (ms_frame.cluster_id == ci)
    ].sort_values("ms_id")
    total_usage = sub_frame.Usage.sum()
    sub_frame['fractional_usage'] = sub_frame.Usage / total_usage

    fig = px.bar(
        sub_frame, x = 'ms_id', y = "Usage",
        color = 'fractional_usage'
    )
    fig.update_layout(title = "Individual Sphere Usage")

    return fig

@app.callback(
    Output("Co-Occurrence", "figure"),
    Input("Receptor", "value"),
    Input("MatchType", "value"),
    Input("Cluster ID", "value")
)
def update_cooccurrence(rec, mt, ci):
    sub_frame = cooc_frame[
        (cooc_frame.match_type == mt) &
        (cooc_frame.receptor == rec) &
        (cooc_frame.cluster_id == ci)
    ]

    mat = sub_frame.iloc[:, 3:].values

    fig = go.Figure()
    trace = go.Heatmap(
        x = ["sph.{}".format(i) for i in range(mat.shape[0])],
        y = ["sph.{}".format(i) for i in range(mat.shape[0])],
        z = mat,
        colorscale = "Viridis"
    )
    fig.add_trace(trace)

    fig.update_layout(title = "Matching Sphere Co-Occurrence")
    return fig

if __name__ == '__main__':
    app.run_server(debug=False)
