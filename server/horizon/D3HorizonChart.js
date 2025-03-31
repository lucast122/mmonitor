import React from 'react';
import PropTypes from 'prop-types';
import * as d3 from 'd3';

const dhc = require('d3-horizon-chart');
console.log('script loaded');

class D3HorizonChart extends React.Component {

    componentDidMount() {
        console.log('Component did mount');
        this.initD3Chart();
    }

    componentDidUpdate(prevProps) {
        console.log('Component did update');
        if (this.props.data !== prevProps.data) {
            this.updateD3Chart();
        }
    }

    initD3Chart() {
        console.log('Initializing D3 Chart with data:', this.props.data);

        // Assuming this.props.data is the data array for the horizon chart
        // generate some ransvm data
        var series = [];
        for (var i = 0, variance = 0; i < 1500; i++) {
            variance += (Math.random() - 0.5) / 10;
            series.push(Math.cos(i / 100) + variance);
        }
        //
        var horizonChart = dhc.horizonChart()
            .height(100)
            .title('Horizon, 4-band')
            .colors(['#313695', '#4575b4', '#74add1', '#abd9e9', '#fee090', '#fdae61', '#f46d43', '#d73027']);

        d3.select(this.el).selectAll('.horizon')
            .data([series])
            .enter().append('div')
            .attr('class', 'horizon')
            .each(function (d) {
                horizonChart.call(this, d);
            });


    }

    updateD3Chart() {
        console.log('Updating D3 Chart with new data:', this.props.data);

        // Add logic to update the chart with new data
    }

    render() {
        console.log('Rendering horizon chart in div:');
        return <div id="d3-horizon-chart-container" ref={el => this.el = el} />;
    }
}

D3HorizonChart.propTypes = {
    data: PropTypes.array.isRequired
};

export default D3HorizonChart;
