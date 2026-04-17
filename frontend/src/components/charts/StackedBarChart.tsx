import { useMemo } from 'react'
import Plot from 'react-plotly.js'
import type { TaxonomyAggregated } from '@/types'
import { getColorPalette } from '@/lib/utils'

interface StackedBarChartProps {
  data: TaxonomyAggregated[]
  plotType?: 'stacked_bar' | 'grouped_bar' | 'heatmap' | 'area'
}

export function StackedBarChart({ data, plotType = 'stacked_bar' }: StackedBarChartProps) {
  const plotData = useMemo(() => {
    if (!data || data.length === 0) return []

    // Get all unique taxa across all samples
    const allTaxa = new Set<string>()
    data.forEach((sample) => {
      sample.taxa.forEach((t) => allTaxa.add(t))
    })
    const taxaList = Array.from(allTaxa)

    // Get colors
    const colors = getColorPalette(taxaList.length)

    // Create traces for each taxon
    return taxaList.map((taxon, i) => {
      const y = data.map((sample) => {
        const idx = sample.taxa.indexOf(taxon)
        return idx >= 0 ? sample.relative_abundances[idx] : 0
      })

      return {
        name: taxon.length > 25 ? taxon.substring(0, 25) + '...' : taxon,
        x: data.map((s) => s.sample_name),
        y,
        type: plotType === 'area' ? 'scatter' : 'bar',
        fill: plotType === 'area' ? 'tonexty' : undefined,
        stackgroup: plotType === 'area' ? 'one' : undefined,
        marker: { color: colors[i] },
        hovertemplate: `<b>${taxon}</b><br>Sample: %{x}<br>Abundance: %{y:.2f}%<extra></extra>`,
      }
    })
  }, [data, plotType])

  const numTaxa = useMemo(() => {
    const allTaxa = new Set<string>()
    data.forEach((sample) => {
      sample.taxa.forEach((t) => allTaxa.add(t))
    })
    return allTaxa.size
  }, [data])

  const layout = useMemo(
    () => ({
      barmode: (plotType === 'stacked_bar' || plotType === 'area' ? 'stack' : 'group') as 'stack' | 'group',
      xaxis: {
        title: { text: 'Sample', standoff: 20 },
        tickangle: -45,
        automargin: true,
      },
      yaxis: {
        title: { text: 'Relative Abundance (%)' },
        range: [0, 100],
      },
      legend: {
        orientation: 'v' as const,
        y: 1,
        x: 1.02,
        xanchor: 'left' as const,
        yanchor: 'top' as const,
        font: { size: 10 },
        tracegroupgap: 2,
      },
      margin: { t: 30, b: 120, l: 60, r: 180 },
      height: Math.max(500, numTaxa * 15 + 200),
      autosize: true,
      hoverlabel: {
        bgcolor: 'white',
        bordercolor: '#ccc',
        font: { size: 12 },
      },
    }),
    [plotType, numTaxa]
  )

  if (plotType === 'heatmap') {
    // Build heatmap data
    const allTaxa = new Set<string>()
    data.forEach((sample) => {
      sample.taxa.forEach((t) => allTaxa.add(t))
    })
    const taxaList = Array.from(allTaxa)

    const z = taxaList.map((taxon) =>
      data.map((sample) => {
        const idx = sample.taxa.indexOf(taxon)
        return idx >= 0 ? sample.relative_abundances[idx] : 0
      })
    )

    // Truncate long taxon names for y-axis display
    const taxaLabels = taxaList.map((t) => (t.length > 30 ? t.substring(0, 30) + '...' : t))

    return (
      <Plot
        data={[
          {
            type: 'heatmap',
            z,
            x: data.map((s) => s.sample_name),
            y: taxaLabels,
            colorscale: 'Viridis',
            colorbar: { title: { text: 'Rel. Abundance (%)' } },
            hovertemplate: 'Sample: %{x}<br>Taxon: %{y}<br>Abundance: %{z:.2f}%<extra></extra>',
          },
        ]}
        layout={{
          xaxis: { title: { text: 'Sample' }, tickangle: -45, automargin: true },
          yaxis: { title: { text: 'Taxon' }, automargin: true },
          margin: { t: 30, b: 120, l: 200, r: 100 },
          height: Math.max(500, taxaList.length * 20 + 150),
          autosize: true,
        }}
        config={{ responsive: true, displayModeBar: true }}
        style={{ width: '100%' }}
      />
    )
  }

  return (
    <Plot
      data={plotData as Plotly.Data[]}
      layout={layout}
      config={{ responsive: true, displayModeBar: true }}
      style={{ width: '100%' }}
    />
  )
}
