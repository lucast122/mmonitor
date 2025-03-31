// assets/download_svg.js
window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        download_svg: function(trigger) {
            const plots = document.querySelectorAll('.js-plotly-plot');
            plots.forEach(plot => {
                const svg = Plotly.Snapshot.toSVG(plot);
                // Process to trigger download for each SVG...
            });
            console.log("Downloading svg plots...")
            return ''; // You might want to return some status or indication of completion
        }
    }
});
