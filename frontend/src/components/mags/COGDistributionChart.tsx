import { useMemo, useRef } from 'react'
import Plot from 'react-plotly.js'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { ChartExportButton } from '@/components/charts/ChartExportButton'
import type { COGDistributionItem } from '@/types'

interface COGDistributionChartProps {
  data: COGDistributionItem[]
  isLoading?: boolean
}

// Color scheme for COG categories
const COG_COLORS: Record<string, string> = {
  J: '#1f77b4', // Translation
  A: '#ff7f0e', // RNA processing
  K: '#2ca02c', // Transcription
  L: '#d62728', // Replication
  B: '#9467bd', // Chromatin
  D: '#8c564b', // Cell cycle
  V: '#e377c2', // Defense
  T: '#7f7f7f', // Signal transduction
  M: '#bcbd22', // Cell wall
  N: '#17becf', // Cell motility
  U: '#aec7e8', // Intracellular trafficking
  O: '#ffbb78', // Posttranslational modification
  X: '#98df8a', // Mobilome
  C: '#ff9896', // Energy production
  G: '#c5b0d5', // Carbohydrate metabolism
  E: '#c49c94', // Amino acid metabolism
  F: '#f7b6d2', // Nucleotide metabolism
  H: '#c7c7c7', // Coenzyme metabolism
  I: '#dbdb8d', // Lipid metabolism
  P: '#9edae5', // Inorganic ion transport
  Q: '#393b79', // Secondary metabolites
  R: '#637939', // General function
  S: '#8c6d31', // Function unknown
}

export function COGDistributionChart({ data, isLoading }: COGDistributionChartProps) {
  const cogChartRef = useRef<HTMLDivElement>(null)
  const chartData = useMemo(() => {
    if (!data || data.length === 0) return null

    const sorted = [...data].sort((a, b) => b.count - a.count)

    return {
      x: sorted.map((d) => d.count),
      y: sorted.map((d) => `${d.category}: ${d.description.slice(0, 30)}${d.description.length > 30 ? '...' : ''}`),
      type: 'bar' as const,
      orientation: 'h' as const,
      marker: {
        color: sorted.map((d) => COG_COLORS[d.category] || '#666666'),
      },
      hovertemplate: sorted.map(
        (d) => `<b>${d.category}</b>: ${d.description}<br>Count: ${d.count}<extra></extra>`
      ),
    }
  }, [data])

  if (isLoading) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>COG Category Distribution</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="flex h-64 items-center justify-center">
            <div className="animate-pulse text-muted-foreground">Loading...</div>
          </div>
        </CardContent>
      </Card>
    )
  }

  if (!data || data.length === 0) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>COG Category Distribution</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="flex h-64 items-center justify-center text-muted-foreground">
            No COG annotation data available
          </div>
        </CardContent>
      </Card>
    )
  }

  return (
    <Card>
      <CardHeader className="flex flex-row items-center justify-between">
        <CardTitle>COG Category Distribution</CardTitle>
        <ChartExportButton containerRef={cogChartRef} filename="cog-distribution" />
      </CardHeader>
      <CardContent>
        <div ref={cogChartRef}>
          <Plot
            data={chartData ? [chartData] : []}
            layout={{
              height: Math.max(400, data.length * 25),
              margin: { l: 250, r: 30, t: 10, b: 40 },
              xaxis: {
                title: { text: 'Gene Count' },
              },
              yaxis: {
                automargin: true,
              },
              paper_bgcolor: 'transparent',
              plot_bgcolor: 'transparent',
            }}
            config={{
              responsive: true,
              displayModeBar: false,
            }}
            style={{ width: '100%' }}
          />
        </div>
      </CardContent>
    </Card>
  )
}
