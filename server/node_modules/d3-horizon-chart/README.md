[![d3 horizon chart](http://kmandov.github.io/d3-horizon-chart/img/d3-horizon-charts-lead-01.png)](http://kmandov.github.io/d3-horizon-chart/)

# d3-horizon-chart

[![npm version](https://badge.fury.io/js/d3-horizon-chart.svg)](http://badge.fury.io/js/d3-horizon-chart)

`d3-horizon-chart` is a d3 plugin that draws horizon charts using an Html5 Canvas.
It provides an easy way to visualize large amounts of time series data.

The plugin is heavily inspired by [cubism.js](https://square.github.io/cubism/), but doesn't make assumptions about real-time metrics and back-end servers such as Cube/Graphite.
This makes it a good fit when you just want to plot your data as a horizon chart.

`d3-horizon-chart` follows the [latest plugin guidelines](https://bost.ocks.org/mike/d3-plugin/) using D3’s new 4.0 module pattern.


## Examples

- [Drawing a single horizon chart](http://bl.ocks.org/kmandov/a1abe4aa380fb8b4bd0b4c081a76ce13)
- [Loosing resolution when using area chart](http://bl.ocks.org/kmandov/5af65af3875c5c4afcdc0d675f60bb45)
- [Horizon charts preserve resolution at smaller heights](http://bl.ocks.org/kmandov/6e91165d4f32534ec4cab806b18b6684)
- [Visualizing stock data using horizon charts](http://bl.ocks.org/kmandov/dcb94f02e71e3560eb3255f2de3120e4)


Check out the [project page](http://kmandov.github.io/d3-horizon-chart/) for more examples.

## Installing

If you use NPM, `npm install d3-horizon-chart`. Otherwise, download the [latest release](https://github.com/kmandov/d3-horizon-chart/releases/latest).

## API Reference

<a name="horizon" href="#horizon">#</a> <b>horizonChart</b>()

Constructs a new horizonChart generator with the default settings.

<a name="horizon_height" href="#horizon_height">#</a> <i>horizon</i>.<b>height</b>([<i>height</i>])

If <i>height</i> is specified, sets the height of the chart to the specified number and returns this horizonChart generator. If <i>height</i> is not specified, returns the current height accessor, which defaults to:

```js
function height() {
  return 30;
}
```

