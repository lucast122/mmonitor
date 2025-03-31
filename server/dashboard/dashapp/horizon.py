import dash
from dash import dcc, html
import dash_mantine_components as dmc
from django_plotly_dash import DjangoDash
from dash_iconify import DashIconify
from dash.dependencies import Input, Output, State
import pandas as pd
from datetime import datetime
from users.models import NanoporeRecord
from functools import lru_cache
from dash.exceptions import PreventUpdate

class Horizon:
    def __init__(self, user_id):
        self.app = DjangoDash('horizon')
        self.user_id = user_id
        self._cache = {}
        self.df = self._get_cached_data()
        # Get unique projects for dropdown
        self.projects = sorted(self.df['project_id'].unique())
        print(f"Initialized Horizon app with data shape: {self.df.shape}")
        self._init_layout()

    @lru_cache(maxsize=32)
    def _get_cached_data(self):
        if 'df' in self._cache:
            return self._cache['df']
            
        records = NanoporeRecord.objects.filter(user_id=self.user_id).select_related()
        df = pd.DataFrame.from_records(records.values())
        
        # Convert date to datetime if it's not already
        df['date'] = pd.to_datetime(df['date'])
        
        # Pre-compute common aggregations
        df['taxa_stats'] = df.groupby(['date', 'taxonomy'])['count'].transform(lambda x: {
            'sum': x.sum(),
            'mean': x.mean(),
            'std': x.std()
        })
        
        # Pre-compute date range
        df['date_range'] = {
            'min': df['date'].min(),
            'max': df['date'].max()
        }
        
        self._cache['df'] = df
        print(f"Fetched data with shape: {df.shape}")
        print(f"Columns: {df.columns}")
        return df

    @lru_cache(maxsize=32)
    def _get_cached_horizon_plot(self, date_range=None, selected_project=None, taxa_count=20):
        """Get cached horizon plot data with proper project filtering."""
        try:
            print(f"\n=== Debug: _get_cached_horizon_plot ===")
            print(f"Selected project: {selected_project}")
            print(f"Initial data shape: {self.df.shape}")
            
            df = self.df.copy()
            
            # Apply project filter before any calculations
            if selected_project and selected_project != 'ALL':
                print(f"Filtering for project: {selected_project}")
                print(f"Project IDs in data: {df['project_id'].unique()}")
                print(f"Project ID types in data: {df['project_id'].apply(type).unique()}")
                print(f"Selected project type: {type(selected_project)}")
                
                # Convert both to strings for comparison
                df = df[df['project_id'].astype(str) == str(selected_project)]
                print(f"Data shape after project filtering: {df.shape}")
                print(f"Remaining samples after filtering: {len(df['sample_id'].unique())}")
                print(f"Unique dates after filtering: {df['date'].nunique()}")

            if len(df) == 0:
                print("No data after filtering")
                return None

            # First calculate total counts per date to use for normalization
            total_counts = df.groupby('date')['count'].sum()
            print(f"Total counts shape: {total_counts.shape}")
            
            # Group by date and taxonomy, using count
            taxa_counts = df.groupby(['date', 'taxonomy'])['count'].sum().unstack(fill_value=0)
            print(f"Taxa counts shape: {taxa_counts.shape}")
            
            # Normalize counts by total sample size for each date
            taxa_abundances = taxa_counts.div(total_counts, axis=0) * 100  # Convert to percentage
            print(f"Taxa abundances shape: {taxa_abundances.shape}")
            
            # Calculate mean abundance for each taxon across all dates
            mean_abundances = taxa_abundances.mean()
            
            # Get top taxa by mean abundance
            top_taxa = mean_abundances.sort_values(ascending=False).head(taxa_count).index.tolist()
            print(f"Number of top taxa: {len(top_taxa)}")
            
            # Calculate differences from mean for each taxon (still in percentage terms)
            differences = taxa_abundances[top_taxa].subtract(mean_abundances[top_taxa], axis=1)
            print(f"Final differences shape: {differences.shape}")
            
            plot_data = differences.to_json(orient='split')
            return plot_data
            
        except Exception as e:
            print(f"Error in _get_cached_horizon_plot: {str(e)}")
            import traceback
            traceback.print_exc()
            return None

    def _get_horizon_plot_html(self, df=None, taxa_count=20):
        """Generate horizon plot HTML with proper project filtering."""
        if df is None:
            print("ERROR: No DataFrame provided")
            return "<p>No data available</p>"
        
        if len(df) == 0:
            print("ERROR: Empty DataFrame")
            return "<p>No data available for the selected filters</p>"
        
        print(f"\n=== Debug: _get_horizon_plot_html ===")
        print(f"Input data shape: {df.shape}")
        print(f"Sample taxa: {df['taxonomy'].head().tolist()}")
        
        # Calculate total counts per date for normalization
        total_counts = df.groupby('date')['count'].sum()
        print(f"Total counts shape: {total_counts.shape}")
        
        # Group by date and taxonomy, using count
        taxa_counts = df.groupby(['date', 'taxonomy'])['count'].sum().unstack(fill_value=0)
        print(f"Taxa counts shape: {taxa_counts.shape}")
        
        # Normalize counts by total sample size for each date
        taxa_abundances = taxa_counts.div(total_counts, axis=0) * 100  # Convert to percentage
        print(f"Taxa abundances shape: {taxa_abundances.shape}")
        
        # Calculate mean abundance for each taxon across all dates
        mean_abundances = taxa_abundances.mean()
        print(f"Mean abundances shape: {mean_abundances.shape}")
        
        # Get top taxa by mean abundance
        top_taxa = mean_abundances.sort_values(ascending=False).head(taxa_count).index.tolist()
        print(f"Number of top taxa: {len(top_taxa)}")
        print(f"Sample top taxa: {top_taxa[:5]}")
        
        # Calculate differences from mean for each taxon (still in percentage terms)
        differences = taxa_abundances[top_taxa].subtract(mean_abundances[top_taxa], axis=1)
        print(f"Final differences shape: {differences.shape}")
        
        plot_data = differences.to_json(orient='split')
        
        if plot_data is None:
            print("ERROR: Failed to generate plot data")
            return "<p>Error processing data</p>"
        
        html_content = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <script src="https://d3js.org/d3.v6.min.js"></script>
    <script src="https://unpkg.com/d3-horizon"></script>
    <script src="https://html2canvas.hertzen.com/dist/html2canvas.min.js"></script>
    <style>
    body { 
        font-family: Arial, sans-serif;
        margin: 0;
        padding: 0;
        overflow-x: hidden;
        width: 100%;
    }
    #plot-content {
        margin: 0;
        padding: 0;
        position: relative;
        width: 100%;
        overflow-x: hidden;
        background: white;
    }
    .horizon-rows {
        margin: 5px 0 30px 10px;
        padding: 0;
        overflow-x: hidden;
        width: calc(100% - 20px);
    }
    .horizon-row {
        display: flex;
        align-items: center;
        margin-bottom: 1px;
        padding: 0;
        height: 35px;
        position: relative;
        background: rgba(255,255,255,0.5);
        overflow-x: hidden;
        width: 100%;
    }
    .horizon-label {
        width: 320px;
        text-align: right;
        padding-right: 10px;
        font-size: 18px;
        font-style: italic;
        margin: 0;
        height: 35px;
        line-height: 35px;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
        position: absolute;
        left: 0;
        background: white;
        z-index: 2;
    }
    .horizon-chart {
        flex-grow: 1;
        margin-left: 300px;
        padding: 0;
        height: 35px;
        overflow-x: hidden;
        width: calc(100% - 300px);
    }
    .x-axis {
        margin: 5px 0 20px 330px;
        height: 30px;
        font-family: Arial, sans-serif;
        overflow-x: hidden;
        width: calc(100% - 330px); /* Slightly narrower than the chart */
    }
    .x-axis-label {
        text-align: center;
        margin-top: 10px;
        margin-left: 300px;
    }
    .y-axis-label {
        position: absolute;
        left: 20px;
        top: 50%;
        transform: rotate(-90deg) translateX(50%);
        transform-origin: left;
    }
    .plot-title {
        font-size: 24px;
        font-weight: bold;
        text-align: center;
        margin: 20px 0 30px 0;
        color: #333;
        padding-left: 300px;
    }
    .axis-label {
        font-size: 20px;
        font-weight: 500;
        color: #333;
    }
    .horizon-legend {
        position: relative;
        margin: 20px auto 30px auto;
        background: transparent;
        padding: 15px 20px;
        font-size: 20px;
        text-align: center;
        width: fit-content;
    }
    .legend-section {
        display: flex;
        justify-content: center;
        align-items: center;
        gap: 20px;
        margin-top: 10px;
    }
    .legend-column {
        display: flex;
        align-items: center;
        gap: 15px;
    }
    .legend-item {
        display: flex;
        align-items: center;
        gap: 4px;
    }
    .legend-color {
        width: 25px;
        height: 15px;
    }
    .legend-label {
        font-size: 20px;
        color: #333;
        white-space: nowrap;
    }
    .d3-horizon-tooltip {
        position: absolute;
        z-index: 10000 !important;
        background: rgba(0, 0, 0, 0.8);
        color: #fff;
        padding: 10px;
        border-radius: 5px;
        pointer-events: none;
        font-size: 14px;
        max-width: 200px;
        box-shadow: 0 0 10px rgba(0, 0, 0, 0.3);
    }
    </style>
</head>
<body>
    <div id="plot-content">
        <div class="plot-title">Species Abundance Changes Over Time</div>
        <div class="y-axis-label axis-label">Species</div>
        <div class="horizon-legend">
            <div style="font-weight: bold; margin-bottom: 10px; text-align: center;">Abundance Changes from Mean</div>
            <div class="legend-section">
                <div style="font-weight: 600; color: #666; margin-right: 10px;">Below Average</div>
                <div class="legend-column">
                    <div class="legend-item">
                        <div class="legend-color" style="background: rgb(180, 180, 230)"></div>
                        <span class="legend-label">0-10%</span>
                    </div>
                    <div class="legend-item">
                        <div class="legend-color" style="background: rgb(130, 130, 200)"></div>
                        <span class="legend-label">10-20%</span>
                    </div>
                    <div class="legend-item">
                        <div class="legend-color" style="background: rgb(80, 80, 170)"></div>
                        <span class="legend-label">20-30%</span>
                    </div>
                    <div class="legend-item">
                        <div class="legend-color" style="background: rgb(25, 25, 112)"></div>
                        <span class="legend-label">>30%</span>
                    </div>
                </div>
                <div style="font-weight: 600; color: #666; margin: 0 10px;">Above Average</div>
                <div class="legend-column">
                    <div class="legend-item">
                        <div class="legend-color" style="background: rgb(255, 200, 210)"></div>
                        <span class="legend-label">0-10%</span>
                    </div>
                    <div class="legend-item">
                        <div class="legend-color" style="background: rgb(255, 150, 160)"></div>
                        <span class="legend-label">10-20%</span>
                    </div>
                    <div class="legend-item">
                        <div class="legend-color" style="background: rgb(255, 90, 100)"></div>
                        <span class="legend-label">20-30%</span>
                    </div>
                    <div class="legend-item">
                        <div class="legend-color" style="background: rgb(255, 20, 60)"></div>
                        <span class="legend-label">>30%</span>
                    </div>
                </div>
            </div>
        </div>
        <div id="horizon-rows" class="horizon-rows"></div>
        <div id="x-axis" class="x-axis"></div>
        <div class="x-axis-label axis-label">Date</div>
    </div>
    <script>
    // Parse the data
    const data = """ + plot_data + """;
    const dates = data.index.map(d => new Date(d));
    const values = data.data;
    const taxa = data.columns;

    // Find max value for legend scaling
    const maxValue = Math.max(...values.flat().map(Math.abs));
    
    // Calculate the width based on container
    const container = document.getElementById('plot-content');
    const margin = { right: 60, left: 300 };  // Increased right margin
    const width = container.clientWidth - margin.left - margin.right;
    
    // Create tooltip div for hover information
    const tooltip = d3.select("body").append("div")
        .attr("class", "d3-horizon-tooltip")
        .style("opacity", 0)
        .style("position", "absolute")
        .style("z-index", 10000);
        
    // Track mouse position globally
    let mouseX = 0;
    let mouseY = 0;
    
    // Update mouse position on mousemove
    document.addEventListener('mousemove', function(event) {
        mouseX = event.pageX;
        mouseY = event.pageY;
    });
        
    // Create a container for each taxon
    taxa.forEach((taxon, i) => {
        // Create a container with label and chart
        const row = d3.select('#horizon-rows')
            .append('div')
            .attr('class', 'horizon-row')
            .style('position', 'relative');
            
        row.append('div')
            .attr('class', 'horizon-label')
            .text(taxon);
            
        const chartDiv = row.append('div')
            .attr('class', 'horizon-chart')
            .attr('id', 'horizon-' + i);
            
        // Get the data for this taxon in [timestamp, value] format
        const taxonData = dates.map((date, j) => [date.getTime(), values[j][i]]);

        try {
            // Create the chart using the Horizon constructor
            const chart = new Horizon(document.getElementById('horizon-' + i))
                .width(width)
                .height(35)
                .bands(5)
                .duration(2000)
                .mode('mirror')
                .positiveColors(['white', 'crimson'])
                .negativeColors(['white', 'midnightblue'])
                .tooltipContent(({ x, y }) => {
                    const date = new Date(x).toLocaleDateString();
                    const value = Math.abs(y).toFixed(2);
                    const direction = y > 0 ? 'above' : 'below';
                    
                    tooltip.transition()
                        .duration(200)
                        .style("opacity", .9);
                    
                    tooltip.html(`<b>${date}</b>: ${value}% ${direction} average`)
                        .style("left", (mouseX + 15) + "px")
                        .style("top", (mouseY - 40) + "px");
                        
                    return "";
                });

            // Set the data
            chart.data(taxonData);
            
            // Add mouseout event to hide tooltip
            document.getElementById('horizon-' + i).addEventListener('mouseout', function() {
                tooltip.transition()
                    .duration(500)
                    .style("opacity", 0);
            });
            
        } catch (error) {
            console.error("Error creating chart for taxon", taxon, ":", error);
        }
    });
    
    // Create x-axis
    const xScale = d3.scaleTime()
        .domain([dates[0], dates[dates.length - 1]])
        .range([0, width - 30]);  // Further reduce width to ensure alignment
    
    // Filter dates based on proximity in pixels
    const minLabelSpacing = 100;
    const filteredDates = [];
    
    if (dates.length > 0) {
        filteredDates.push(dates[0]);
        let lastIncludedDate = dates[0];
        
        for (let i = 1; i < dates.length; i++) {
            const currentDate = dates[i];
            const pixelDistance = Math.abs(xScale(currentDate) - xScale(lastIncludedDate));
            
            if (pixelDistance >= minLabelSpacing) {
                filteredDates.push(currentDate);
                lastIncludedDate = currentDate;
            }
        }
        
        if (dates.length > 1 && filteredDates[filteredDates.length - 1] !== dates[dates.length - 1]) {
            filteredDates.push(dates[dates.length - 1]);
        }
    }
        
    const xAxis = d3.axisBottom(xScale)
        .tickValues(filteredDates)
        .tickFormat(d3.timeFormat("%d %b"))
        .tickSize(-6);
        
    const xAxisSvg = d3.select("#x-axis")
        .append("svg")
        .attr("width", width - 30)  // Match the width from range
        .attr("height", 40)
        .append("g")
        .attr("transform", "translate(0,0)")
        .call(xAxis);
    
    xAxisSvg.select(".domain")
        .style("stroke", "#999")
        .style("stroke-width", "1px");
        
    xAxisSvg.selectAll(".tick line")
        .style("stroke", "#999")
        .style("stroke-width", "1px");
        
    xAxisSvg.selectAll(".tick text")
        .style("font-family", "Arial, sans-serif")
        .style("font-size", "16px")
        .style("fill", "#333")
        .style("text-anchor", "start")
        .attr("dx", "0.5em")
        .attr("dy", "0.3em")
        .attr("transform", "rotate(25)");

    // Handle window resize
    window.addEventListener('resize', function() {
        const newWidth = container.clientWidth - margin.left - margin.right;
        taxa.forEach((_, i) => {
            const chart = new Horizon(document.getElementById('horizon-' + i))
                .width(newWidth);
            chart.data(dates.map((date, j) => [date.getTime(), values[j][i]]));
        });
    });
    </script>
</body>
</html>
"""
        return html_content
    
    def _init_layout(self):
        """Initialize the layout of the horizon plot page."""
        self.app.layout = dmc.Container(
            dmc.Stack([
                dmc.Title("Horizon Plot", order=2, style={"marginBottom": "5px"}),
                
                # Add explanation for horizon plots
                dmc.Accordion(
                    children=[
                        dmc.AccordionItem(
                            [
                                dmc.AccordionControl("About Horizon Plots"),
                                dmc.AccordionPanel(
                                    dmc.Alert(
                                        title="About Horizon Plots",
                                        children=[
                                            "Horizon plots are a space-efficient way to visualize time series data:",
                                            html.Ul([
                                                html.Li("Each row represents a different taxon's abundance over time"),
                                                html.Li("Colors indicate the magnitude of values:"),
                                                html.Ul([
                                                    html.Li("Darker colors represent higher values"),
                                                    html.Li("Lighter colors represent lower values")
                                                ]),
                                                html.Li("Hover over any point to see exact values and timestamps"),
                                                html.Li("Use the controls below to filter by project and adjust the number of taxa shown")
                                            ], style={"marginBottom": "0", "marginTop": "0", "paddingTop": "0", 'maxWidth': '100%', 'width': '100%', 'minWidth': '100%'})
                                        ],
                                        color="blue",
                                        variant="light",
                                        style={"marginBottom": "5px", "padding": "5px 10px", 'maxWidth': '100%', 'width': '100%', 'minWidth': '100%'}
                                    )
                                )
                            ],
                            value="about"
                        )
                    ]
                ),

                # Controls group
                dmc.Group([
                    dmc.Select(
                        id='project-select-horizon',
                        label="Select Project",
                        data=[{'value': 'ALL', 'label': 'All Projects'}] + [
                            {'value': str(project_id), 'label': f'Project {project_id}'} 
                            for project_id in self.df['project_id'].unique() if pd.notna(project_id)
                        ],
                        value='ALL',
                        style={"width": "200px"}
                    ),
                    dmc.Select(
                        id='subproject-select-horizon',
                        label="Select Subproject",
                        data=[{'value': 'ALL', 'label': 'All Subprojects'}] + (
                            [{'value': str(subproject), 'label': f'Subproject {subproject}'} 
                             for subproject in self.df['subproject'].unique() if pd.notna(subproject)]
                            if 'subproject' in self.df.columns else []
                        ),
                        value='ALL',
                        style={"width": "200px"}
                    ),
                    dmc.NumberInput(
                        id='taxa-count',
                        label="Number of Taxa",
                        value=20,
                        min=1,
                        max=100,
                        style={"width": "150px"}
                    ),
                ], align="flex-end", style={"padding": "5px 0"}),
                
                # Update the plot container to ensure tooltips are visible
                dmc.Paper(
                    id='horizon-plot-container',
                    children=[
                        html.Iframe(
                            id='horizon-plot-iframe',
                            srcDoc=self._get_horizon_plot_html(df=self.df),
                            style={
                                'width': '100%',
                                'height': '100%',
                                'border': 'none',
                                'padding': '0',
                                'margin': '0',
                                'position': 'relative',
                                'z-index': 1
                            },
                            sandbox='allow-scripts allow-same-origin allow-popups allow-popups-to-escape-sandbox'
                        )
                    ],
                    shadow="sm",
                    p="0",
                    style={
                        "width": "95vw", 
                        "maxWidth": "90vw",
                        "marginLeft": "calc(-1 * var(--mantine-spacing-md))",
                        "marginRight": "calc(-1 * var(--mantine-spacing-md))",
                        "padding": "0",
                        "height": f"{50 + 50 * 20}px",
                        "position": "relative",
                        "z-index": 0,
                        "overflow": "hidden"
                    }
                ),
            ], spacing="xs"),
            fluid=True,
            p=0,
            style={'width': '100%', 'margin': '0', 'padding': '0', 'maxWidth': '100%', 'overflow': 'hidden'}
        )

        # Let's also update the HTML content to use full width
        html_style = """
    <style>
    body { 
        font-family: Arial, sans-serif;
        margin: 0;
        padding: 0;
        overflow-x: hidden;
        width: 100%;
    }
    #plot-content {
        margin: 0;
        padding: 0;
        position: relative;
        width: 100%;
        overflow-x: hidden;
        background: white;
    }
    .horizon-rows {
        margin: 5px 0 30px 10px;
        padding: 0;
        overflow-x: hidden;
        width: calc(100% - 20px);
    }
    """
        
        # Add callback for dynamic updates
        @self.app.callback(
            Output('horizon-plot-iframe', 'srcDoc'),
            [
                Input('project-select-horizon', 'value'),
                Input('subproject-select-horizon', 'value'),
                Input('taxa-count', 'value')
            ]
        )
        def update_plot(project, subproject, taxa_count):
            try:
                df = self.df.copy()
                
                # Filter by project if needed
                if project and project != 'ALL':
                    df = df[df['project_id'].astype(str) == str(project)]
                
                # Filter by subproject if needed
                if subproject and subproject != 'ALL' and 'subproject' in df.columns:
                    df = df[df['subproject'].astype(str) == str(subproject)]
                
                # Generate the plot HTML
                return self._get_horizon_plot_html(df=df, taxa_count=int(taxa_count) if taxa_count else 20)
            except Exception as e:
                import traceback
                traceback.print_exc()
                return f"<p>Error generating plot: {str(e)}</p>"
        
        
        