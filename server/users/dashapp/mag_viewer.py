from dash import html, dcc, dash_table
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from django_plotly_dash import DjangoDash
from ..models import MAG
import json
import logging
# import dash_bio as dashbio  # Temporarily commented out
from Bio import SeqIO
from io import StringIO
import os
import traceback

logger = logging.getLogger(__name__)

app = DjangoDash('MagViewer', add_bootstrap_links=True)

def get_mag_data():
    """Helper function to get MAG data and log the results"""
    try:
        # Get all MAGs with detailed logging
        mags = MAG.objects.all()
        logger.info(f"Found {len(mags)} MAGs in database")
        
        # Debug each MAG's data
        for mag in mags:
            logger.info(f"\nMAG Details for {mag.name}:")
            logger.info(f"  Sample: {mag.sample_name}")
            logger.info(f"  Taxonomy: {mag.taxonomy}")
            logger.info(f"  Completeness: {mag.completeness}%")
            logger.info(f"  Contamination: {mag.contamination}%")
            logger.info(f"  Created at: {mag.created_at}")
            
            # Check file data
            logger.info("  File Status:")
            logger.info(f"    FASTA: {'Present' if mag.fasta_file and mag.fasta_file != 'empty' else 'Missing'}")
            if mag.fasta_file and mag.fasta_file != 'empty':
                logger.info(f"    FASTA size: {len(mag.fasta_file)} characters")
            
            logger.info(f"    GFF: {'Present' if mag.gff_file and mag.gff_file != 'empty' else 'Missing'}")
            if mag.gff_file and mag.gff_file != 'empty':
                logger.info(f"    GFF size: {len(mag.gff_file)} characters")
            
            logger.info(f"    FAI: {'Present' if mag.fai_file and mag.fai_file != 'empty' else 'Missing'}")
            
            # Check if all required files are present
            has_all_files = (
                mag.fasta_file and mag.fasta_file != 'empty' and
                mag.gff_file and mag.gff_file != 'empty'
            )
            logger.info(f"  Complete dataset: {'Yes' if has_all_files else 'No'}")
        
        # Return the data for the table
        mag_data = list(mags.values(
            'id', 'name', 'sample_name', 'taxonomy', 
            'completeness', 'contamination', 'strain_heterogeneity',
            'gene_count', 'trna_count', 'rrna_count', 'cds_count'
        ))
        return mag_data
        
    except Exception as e:
        logger.error(f"Error fetching MAGs: {str(e)}")
        logger.error(traceback.format_exc())
        return []

def get_fasta_sequence(mag_id):
    """Get the FASTA sequence for a MAG"""
    try:
        mag = MAG.objects.get(pk=mag_id)
        logger.info(f"\nProcessing FASTA for MAG {mag_id}:")
        logger.info(f"  Name: {mag.name}")
        logger.info(f"  Sample: {mag.sample_name}")
        logger.info(f"  Is Complete: {mag.is_complete}")
        
        # Debug the FASTA data
        if not mag.fasta_file or mag.fasta_file == 'empty':
            logger.error(f"No FASTA data found for MAG {mag_id}")
            return None
            
        logger.info(f"  FASTA size: {len(mag.fasta_file)} characters")
        
        # Parse the FASTA text with BioPython
        sequences = []
        try:
            fasta_io = StringIO(mag.fasta_file)
            for record in SeqIO.parse(fasta_io, 'fasta'):
                logger.info(f"  Found sequence: {record.id}")
                logger.info(f"    Length: {len(record.seq)} bp")
                logger.info(f"    Description: {record.description}")
                sequences.append({
                    'sequence': str(record.seq),
                    'id': record.id,
                    'description': record.description
                })
            fasta_io.close()
        except Exception as e:
            logger.error(f"Error parsing FASTA text: {str(e)}")
            logger.error(traceback.format_exc())
            return None
        
        if not sequences:
            logger.error("No sequences found in FASTA data")
        else:
            logger.info(f"Successfully parsed {len(sequences)} sequences")
            
        return sequences
    except Exception as e:
        logger.error(f"Error reading FASTA for MAG {mag_id}: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def parse_gff_data(gff_text):
    """Parse GFF text into a format suitable for the genome browser"""
    if not gff_text or gff_text == 'empty':
        logger.warning("No GFF data provided")
        return []
    
    try:
        features = []
        feature_counts = {'CDS': 0, 'tRNA': 0, 'rRNA': 0, 'other': 0}
        
        for line in gff_text.split('\n'):
            if line and not line.startswith('#'):
                parts = line.strip().split('\t')
                if len(parts) >= 8:
                    seqid, source, type_, start, end, score, strand, phase, *attrs = parts
                    
                    # Parse attributes
                    attr_dict = {}
                    if len(attrs) > 0:
                        for attr in attrs[0].split(';'):
                            if '=' in attr:
                                key, value = attr.split('=', 1)
                                attr_dict[key] = value
                    
                    feature = {
                        'ref': seqid,
                        'start': int(start),
                        'end': int(end),
                        'type': type_,
                        'strand': strand,
                        'name': attr_dict.get('ID', ''),
                        'description': attr_dict.get('product', ''),
                        'color': '#1f77b4'  # Default color
                    }
                    
                    # Color coding based on feature type
                    if type_ == 'CDS':
                        feature['color'] = '#2ca02c'  # Green
                        feature_counts['CDS'] += 1
                    elif type_ == 'tRNA':
                        feature['color'] = '#ff7f0e'  # Orange
                        feature_counts['tRNA'] += 1
                    elif type_ == 'rRNA':
                        feature['color'] = '#d62728'  # Red
                        feature_counts['rRNA'] += 1
                    else:
                        feature_counts['other'] += 1
                    
                    features.append(feature)
        
        # Log feature statistics
        logger.info("\nGFF Feature Statistics:")
        logger.info(f"  Total features: {len(features)}")
        logger.info(f"  CDS: {feature_counts['CDS']}")
        logger.info(f"  tRNA: {feature_counts['tRNA']}")
        logger.info(f"  rRNA: {feature_counts['rRNA']}")
        logger.info(f"  Other: {feature_counts['other']}")
        
        return features
    except Exception as e:
        logger.error(f"Error parsing GFF data: {str(e)}")
        logger.error(traceback.format_exc())
        return []

# Initialize the layout with data
initial_mags = get_mag_data()
initial_options = [
    {'label': f"{m['name']} (Completeness: {m['completeness']:.1f}%)" if m.get('completeness') is not None else m['name'], 
     'value': m['id']} 
    for m in initial_mags
]

# Create DataFrame for the table
df = pd.DataFrame(initial_mags)
if not df.empty:
    df = df.round({
        'completeness': 1,
        'contamination': 1,
        'strain_heterogeneity': 1
    })
    # Rename columns for better display
    df = df.rename(columns={
        'name': 'Name',
        'sample_name': 'Sample',
        'taxonomy': 'Taxonomy',
        'completeness': 'Completeness (%)',
        'contamination': 'Contamination (%)',
        'strain_heterogeneity': 'Strain Het. (%)',
        'gene_count': 'Genes',
        'trna_count': 'tRNA',
        'rrna_count': 'rRNA',
        'cds_count': 'CDS'
    })

app.layout = html.Div([
    dcc.Store(id='user-id', data=None),
    dbc.Container([
        dbc.Row([
            dbc.Col([
                html.H1("MAG Analysis Dashboard", className="text-center mb-4"),
            ])
        ]),
        dbc.Row([
            dbc.Col([
                html.H4("MAG Overview", className="mb-3"),
                dash_table.DataTable(
                    id='mag-table',
                    columns=[
                        {"name": i, "id": i} for i in df.columns if i != 'id'
                    ],
                    data=df.to_dict('records'),
                    sort_action="native",
                    sort_mode="multi",
                    filter_action="native",
                    style_table={'overflowX': 'auto'},
                    style_cell={
                        'textAlign': 'left',
                        'minWidth': '100px',
                        'maxWidth': '300px',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                    },
                    style_header={
                        'backgroundColor': 'rgb(230, 230, 230)',
                        'fontWeight': 'bold'
                    },
                    style_data_conditional=[
                        {
                            'if': {'row_index': 'odd'},
                            'backgroundColor': 'rgb(248, 248, 248)'
                        }
                    ],
                    tooltip_data=[
                        {
                            column: {'value': str(value), 'type': 'markdown'}
                            for column, value in row.items()
                        } for row in df.to_dict('records')
                    ],
                    tooltip_duration=None,
                    page_size=10
                ),
            ], width=12),
        ]),
        dbc.Row([
            dbc.Col([
                html.H4("Selected MAG Details", className="mt-4 mb-3"),
                dcc.Dropdown(
                    id='mag-selector',
                    options=initial_options,
                    value=initial_options[0]['value'] if initial_options else None,
                    placeholder='Select a MAG...',
                ),
            ])
        ]),
        dbc.Row([
            dbc.Col([
                dcc.Graph(id='quality-plot'),
            ], width=6),
            dbc.Col([
                dcc.Graph(id='annotation-stats'),
            ], width=6),
        ]),
        dbc.Row([
            dbc.Col([
                html.H4("Genome Browser", className="text-center mt-4"),
                dcc.Loading(
                    id="loading-genome",
                    type="default",
                    children=[
                        html.Div(id='genome-browser'),
                    ]
                ),
            ], width=12),
        ]),
    ])
])

@app.callback(
    Output('genome-browser', 'children'),
    Input('mag-selector', 'value')
)
def update_genome_browser(mag_id):
    if not mag_id:
        return html.Div("Please select a MAG")
    
    try:
        # Get the MAG data with detailed logging
        mag = MAG.objects.get(pk=mag_id)
        logger.info(f"\nLoading genome browser for MAG {mag_id}:")
        logger.info(f"  Name: {mag.name}")
        logger.info(f"  Sample: {mag.sample_name}")
        
        # Check FASTA data
        if not mag.fasta_file or mag.fasta_file == 'empty':
            logger.error(f"No FASTA data found for MAG {mag_id}")
            return html.Div([
                dbc.Alert([
                    html.H4("Missing FASTA Data", className="alert-heading"),
                    html.P("This MAG does not have FASTA sequence data. Please ensure the FASTA file was properly uploaded.")
                ], color="warning")
            ])
        
        # Check GFF data
        if not mag.gff_file or mag.gff_file == 'empty':
            logger.warning(f"No GFF annotation data found for MAG {mag_id}")
            features = []
        else:
            features = parse_gff_data(mag.gff_file)
            logger.info(f"  Found {len(features)} features in GFF file")
        
        # Get sequences
        sequences = get_fasta_sequence(mag_id)
        if not sequences:
            logger.error(f"Failed to parse FASTA data for MAG {mag_id}")
            return html.Div([
                dbc.Alert([
                    html.H4("FASTA Parsing Error", className="alert-heading"),
                    html.P("Could not parse the FASTA sequence data. The file may be corrupted.")
                ], color="danger")
            ])
        
        logger.info(f"  Found {len(sequences)} sequences in FASTA file")
        for i, seq in enumerate(sequences):
            logger.info(f"    Sequence {i+1}: {seq['id']} ({len(seq['sequence'])} bp)")
        
        # Create genome browser components
        genome_browsers = []
        for seq in sequences:
            sequence = seq['sequence']
            seq_id = seq['id']
            
            # Filter features for this sequence
            seq_features = [f for f in features if f['ref'] == seq_id]
            logger.info(f"    Features for {seq_id}: {len(seq_features)}")
            
            browser = html.Div([
                html.H5(f"Sequence: {seq_id}", className="mt-3"),
                html.P([
                    html.Strong("Length: "), f"{len(sequence):,} bp",
                    html.Br(),
                    html.Strong("Description: "), seq['description']
                ]),
                dbc.Card([
                    dbc.CardBody([
                        dashbio.GenomeBrowser(
                            id={'type': 'genome-browser', 'index': seq_id},
                            sequenceData=[{
                                'sequence': sequence,
                                'id': seq_id
                            }],
                            annotationData=seq_features,
                            width=1000,
                            height=400,
                            showNavigator=True,
                            showTrackLabel=True,
                            showGutter=True,
                            showTranslation=True,
                            compact=False,
                            scrollDelay=30,
                            showBases=True,
                            showRuler=True
                        ),
                        html.Div([
                            html.H6("Feature Types:", className="mt-3"),
                            html.Div([
                                html.Span("CDS", style={'color': '#2ca02c', 'marginRight': '15px'}),
                                html.Span("tRNA", style={'color': '#ff7f0e', 'marginRight': '15px'}),
                                html.Span("rRNA", style={'color': '#d62728'})
                            ])
                        ], className="mt-3") if features else html.Div([
                            dbc.Alert(
                                "No annotation data available. Features will not be displayed.",
                                color="info",
                                className="mt-3"
                            )
                        ])
                    ])
                ])
            ])
            genome_browsers.append(browser)
        
        return html.Div([
            dbc.Alert([
                html.H4("Genome Browser", className="alert-heading"),
                html.P([
                    "Navigate through the genome sequence and annotations below. ",
                    "Use mouse wheel to zoom, click and drag to pan. ",
                    "Hover over features to see details."
                ])
            ], color="info", className="mb-3"),
            html.Div(genome_browsers)
        ])
        
    except Exception as e:
        logger.error(f"Error in genome browser: {str(e)}")
        logger.error(traceback.format_exc())
        return html.Div([
            dbc.Alert([
                html.H4("Error Loading Genome Browser", className="alert-heading"),
                html.P(str(e)),
                html.Hr(),
                html.P("Check the server logs for more details.", className="mb-0")
            ], color="danger")
        ])

@app.callback(
    Output('mag-details', 'children'),
    Input('mag-selector', 'value')
)
def update_mag_details(mag_id):
    if not mag_id:
        return html.Div("Please select a MAG")
    
    try:
        mag = MAG.objects.get(pk=mag_id)
        return dbc.Card(
            dbc.CardBody([
                html.H4(f"MAG: {mag.name}", className="card-title"),
                html.H6(f"Sample: {mag.sample_name}", className="card-subtitle mb-2 text-muted"),
                dbc.Row([
                    dbc.Col([
                        html.P([
                            html.Strong("Quality Metrics:"),
                            html.Br(),
                            f"Completeness: {mag.completeness:.1f}% " if mag.completeness is not None else "Completeness: N/A",
                            html.Br(),
                            f"Contamination: {mag.contamination:.1f}% " if mag.contamination is not None else "Contamination: N/A",
                            html.Br(),
                            f"Strain Heterogeneity: {mag.strain_heterogeneity:.1f}% " if mag.strain_heterogeneity is not None else "Strain Heterogeneity: N/A",
                        ]),
                    ], width=6),
                    dbc.Col([
                        html.P([
                            html.Strong("Gene Statistics:"),
                            html.Br(),
                            f"Total Genes: {mag.gene_count}" if mag.gene_count is not None else "Total Genes: N/A",
                            html.Br(),
                            f"tRNA Genes: {mag.trna_count}" if mag.trna_count is not None else "tRNA Genes: N/A",
                            html.Br(),
                            f"rRNA Genes: {mag.rrna_count}" if mag.rrna_count is not None else "rRNA Genes: N/A",
                            html.Br(),
                            f"CDS: {mag.cds_count}" if mag.cds_count is not None else "CDS: N/A",
                        ]),
                    ], width=6),
                ]),
            ])
        )
    except MAG.DoesNotExist:
        return html.Div("MAG not found")
    except Exception as e:
        logger.error(f"Error in update_mag_details: {str(e)}")
        return html.Div(f"Error: {str(e)}")

@app.callback(
    Output('quality-plot', 'figure'),
    Input('mag-selector', 'value')
)
def update_quality_plot(mag_id):
    try:
        mags = get_mag_data()
        if not mags:
            return go.Figure()
        
        df = pd.DataFrame(mags)
        
        # Remove rows where completeness or contamination is None
        df = df.dropna(subset=['completeness', 'contamination'])
        
        if len(df) == 0:
            return go.Figure()

        fig = px.scatter(df, x='completeness', y='contamination',
                        hover_data=['name'],
                        title='MAG Quality Assessment')
        
        fig.update_layout(
            xaxis_title="Completeness (%)",
            yaxis_title="Contamination (%)",
            showlegend=True
        )
        
        # Add quality threshold lines
        fig.add_hline(y=5, line_dash="dash", line_color="red", annotation_text="5% Contamination")
        fig.add_vline(x=90, line_dash="dash", line_color="green", annotation_text="90% Completeness")
        
        # Highlight selected MAG if it exists
        if mag_id:
            selected_mag = next((m for m in mags if m['id'] == mag_id), None)
            if selected_mag and selected_mag.get('completeness') is not None and selected_mag.get('contamination') is not None:
                fig.add_trace(go.Scatter(
                    x=[selected_mag['completeness']],
                    y=[selected_mag['contamination']],
                    mode='markers',
                    marker=dict(size=15, color='red'),
                    name='Selected MAG'
                ))
        
        return fig
    except Exception as e:
        logger.error(f"Error in update_quality_plot: {str(e)}")
        return go.Figure()

@app.callback(
    Output('annotation-stats', 'figure'),
    Input('mag-selector', 'value')
)
def update_annotation_stats(mag_id):
    if not mag_id:
        return go.Figure()

    try:
        mag = MAG.objects.get(pk=mag_id)
        stats = {
            'Genes': mag.gene_count or 0,
            'tRNA': mag.trna_count or 0,
            'rRNA': mag.rrna_count or 0,
            'CDS': mag.cds_count or 0
        }
        
        fig = go.Figure(data=[
            go.Bar(
                x=list(stats.keys()),
                y=list(stats.values()),
                text=list(stats.values()),
                textposition='auto',
            )
        ])
        
        fig.update_layout(
            title=f"Gene Content Statistics for {mag.name}",
            yaxis_title="Count",
            showlegend=False
        )
        
        return fig
    except MAG.DoesNotExist:
        return go.Figure()
    except Exception as e:
        logger.error(f"Error in update_annotation_stats: {str(e)}")
        return go.Figure()

@app.callback(
    Output('taxonomy-info', 'children'),
    Input('mag-selector', 'value')
)
def update_taxonomy(mag_id):
    if not mag_id:
        return html.Div()
    
    try:
        mag = MAG.objects.get(pk=mag_id)
        return dbc.Card(
            dbc.CardBody([
                html.H5("Taxonomic Classification", className="card-title"),
                html.P(mag.taxonomy if mag.taxonomy else "No taxonomy information available"),
            ])
        )
    except MAG.DoesNotExist:
        return html.Div("MAG not found")
    except Exception as e:
        logger.error(f"Error in update_taxonomy: {str(e)}")
        return html.Div(f"Error: {str(e)}")
