import { useEffect, useRef } from 'react'
import Horizon from 'd3-horizon'

interface HorizonData {
  dates: string[]
  taxa: string[]
  data: { date: string; taxon: string; value: number; percentage: number; mean: number }[]
  mean_abundances: Record<string, number>
}

interface HorizonChartProps {
  data: HorizonData
  bands?: number
  mode?: 'offset' | 'mirror'
}

export function HorizonChart({ data, bands = 4, mode = 'mirror' }: HorizonChartProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const chartsRef = useRef<Horizon[]>([])

  useEffect(() => {
    if (!containerRef.current || !data || data.taxa.length === 0 || data.dates.length === 0) {
      return
    }

    // Clear previous charts
    const container = containerRef.current
    container.innerHTML = ''
    chartsRef.current = []

    const { dates, taxa } = data

    // Create date index for x-axis positioning
    const dateIndex = new Map(dates.map((d, i) => [d, i]))

    // Group data by taxon
    const dataByTaxon = new Map<string, typeof data.data>()
    for (const d of data.data) {
      if (!dataByTaxon.has(d.taxon)) {
        dataByTaxon.set(d.taxon, [])
      }
      dataByTaxon.get(d.taxon)!.push(d)
    }

    // Find max absolute value for consistent scaling across all rows
    const maxAbsValue = Math.max(
      ...data.data.map((d) => Math.abs(d.value)),
      0.1 // Minimum to avoid division by zero
    )

    // Calculate dimensions
    const rowHeight = 30
    const labelWidth = 180
    const chartWidth = container.clientWidth - labelWidth - 40

    // Create wrapper with flex layout
    const wrapper = document.createElement('div')
    wrapper.style.display = 'flex'
    wrapper.style.flexDirection = 'column'
    wrapper.style.gap = '2px'
    wrapper.style.position = 'relative'
    container.appendChild(wrapper)

    // Add styles for tooltips to ensure they appear in foreground near cursor
    const styleEl = document.createElement('style')
    styleEl.textContent = `
      .float-tooltip-kap {
        z-index: 99999 !important;
        background: white !important;
        color: #333 !important;
        border: 1px solid #999 !important;
        border-radius: 6px !important;
        box-shadow: 0 4px 16px rgba(0,0,0,0.25) !important;
        padding: 10px 12px !important;
        font-size: 13px !important;
        pointer-events: none !important;
        max-width: 300px !important;
        line-height: 1.4 !important;
        overflow: visible !important;
      }
    `
    document.head.appendChild(styleEl)

    // Add legend at the top
    const legendDiv = document.createElement('div')
    legendDiv.style.display = 'flex'
    legendDiv.style.justifyContent = 'flex-end'
    legendDiv.style.alignItems = 'center'
    legendDiv.style.gap = '20px'
    legendDiv.style.padding = '8px 20px'
    legendDiv.style.fontSize = '12px'
    legendDiv.style.marginBottom = '10px'
    legendDiv.innerHTML = `
      <span style="font-weight: 600;">Deviation from Mean:</span>
      <span style="display: flex; align-items: center; gap: 4px;">
        <span style="display: inline-block; width: 40px; height: 14px; background: linear-gradient(to right, white, midnightblue); border: 1px solid #ccc;"></span>
        <span>Below mean</span>
      </span>
      <span style="display: flex; align-items: center; gap: 4px;">
        <span style="display: inline-block; width: 40px; height: 14px; background: linear-gradient(to right, white, crimson); border: 1px solid #ccc;"></span>
        <span>Above mean</span>
      </span>
    `
    wrapper.appendChild(legendDiv)

    // Create a row for each taxon
    taxa.forEach((taxon) => {
      const rowData = dataByTaxon.get(taxon) || []
      if (rowData.length === 0) return

      // Sort by date
      rowData.sort((a, b) => (dateIndex.get(a.date) || 0) - (dateIndex.get(b.date) || 0))

      // Create row container - overflow visible so tooltip isn't clipped
      const rowDiv = document.createElement('div')
      rowDiv.style.display = 'flex'
      rowDiv.style.alignItems = 'center'
      rowDiv.style.borderBottom = '1px solid #e5e7eb'
      rowDiv.style.overflow = 'visible'
      wrapper.appendChild(rowDiv)

      // Create label
      const labelDiv = document.createElement('div')
      labelDiv.style.width = `${labelWidth}px`
      labelDiv.style.flexShrink = '0'
      labelDiv.style.paddingRight = '10px'
      labelDiv.style.textAlign = 'right'
      labelDiv.style.fontSize = '11px'
      labelDiv.style.color = '#333'
      labelDiv.style.overflow = 'hidden'
      labelDiv.style.textOverflow = 'ellipsis'
      labelDiv.style.whiteSpace = 'nowrap'
      labelDiv.title = taxon
      labelDiv.textContent = taxon.length > 25 ? taxon.substring(0, 25) + '...' : taxon
      rowDiv.appendChild(labelDiv)

      // Create chart container - overflow visible so tooltip isn't clipped
      const chartDiv = document.createElement('div')
      chartDiv.style.width = `${chartWidth}px`
      chartDiv.style.height = `${rowHeight}px`
      chartDiv.style.position = 'relative'
      chartDiv.style.overflow = 'visible'
      rowDiv.appendChild(chartDiv)

      // Prepare data for d3-horizon: [x, y] pairs
      const horizonData: [number, number][] = rowData.map((d) => [
        dateIndex.get(d.date) || 0,
        d.value,
      ])

      // Create horizon chart
      const chart = new Horizon(chartDiv, { useCanvas: true })
        .data(horizonData)
        .width(chartWidth)
        .height(rowHeight)
        .bands(bands)
        .mode(mode)
        .yExtent(maxAbsValue)
        .positiveColors(['white', 'crimson'])
        .negativeColors(['white', 'midnightblue'])
        .tooltipContent(({ x, y }: { x: number; y: number }) => {
          const dateStr = dates[Math.round(x)] || ''
          const point = rowData.find((d) => d.date === dateStr)
          if (point) {
            return `
              <div style="padding: 4px;">
                <strong>${taxon}</strong><br/>
                <span style="color: #666;">Date:</span> ${point.date}<br/>
                <span style="color: #666;">Abundance:</span> ${point.percentage.toFixed(2)}%<br/>
                <span style="color: #666;">Mean:</span> ${point.mean.toFixed(2)}%<br/>
                <span style="color: ${point.value >= 0 ? 'crimson' : 'midnightblue'};">
                  ${point.value >= 0 ? '▲' : '▼'} ${Math.abs(point.value).toFixed(2)}% ${point.value >= 0 ? 'above' : 'below'} mean
                </span>
              </div>
            `
          }
          return `x: ${x}, y: ${y.toFixed(2)}`
        })

      chartsRef.current.push(chart)
    })

    // Add x-axis labels at the bottom
    const axisDiv = document.createElement('div')
    axisDiv.style.display = 'flex'
    axisDiv.style.marginLeft = `${labelWidth}px`
    axisDiv.style.marginTop = '10px'
    axisDiv.style.position = 'relative'
    axisDiv.style.height = '100px'  // Increased height for rotated labels
    axisDiv.style.width = `${chartWidth}px`
    wrapper.appendChild(axisDiv)

    // Show subset of date labels to avoid crowding
    const maxLabels = Math.min(dates.length, 10)
    const step = Math.max(1, Math.floor(dates.length / maxLabels))

    dates.forEach((date, i) => {
      if (i % step === 0 || i === dates.length - 1) {
        const labelSpan = document.createElement('span')
        labelSpan.style.position = 'absolute'
        labelSpan.style.left = `${(i / (dates.length - 1 || 1)) * chartWidth}px`
        labelSpan.style.top = '0'
        labelSpan.style.fontSize = '10px'
        labelSpan.style.color = '#666'
        labelSpan.style.transform = 'rotate(-45deg)'
        labelSpan.style.transformOrigin = 'top left'
        labelSpan.style.whiteSpace = 'nowrap'
        labelSpan.textContent = date
        axisDiv.appendChild(labelSpan)
      }
    })

    // Cleanup function
    return () => {
      chartsRef.current = []
      // Remove the style element we added
      if (styleEl.parentNode) {
        styleEl.parentNode.removeChild(styleEl)
      }
    }
  }, [data, bands, mode])

  return (
    <>
      {/* Global styles for d3-horizon tooltips */}
      <style>{`
        .float-tooltip-kap {
          z-index: 99999 !important;
          background: white !important;
          color: #333 !important;
          border: 1px solid #999 !important;
          border-radius: 6px !important;
          box-shadow: 0 4px 16px rgba(0,0,0,0.25) !important;
          padding: 10px 12px !important;
          font-size: 13px !important;
          line-height: 1.4 !important;
          max-width: 300px !important;
          pointer-events: none !important;
          overflow: visible !important;
        }
      `}</style>
      <div
        ref={containerRef}
        className="relative w-full"
        style={{ minHeight: data.taxa.length * 32 + 160, overflow: 'visible' }}
      />
    </>
  )
}
