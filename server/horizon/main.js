// index.js or main.js
import React from 'react';
import ReactDOM from 'react-dom';
import D3HorizonChart from './D3HorizonChart'; // Adjust the import path as needed

document.addEventListener('DOMContentLoaded', () => {
  ReactDOM.render(
    <D3HorizonChart />,
    document.getElementById('d3-horizon-chart-container')
  );
});
