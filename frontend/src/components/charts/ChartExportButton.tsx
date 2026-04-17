import { RefObject, useCallback } from 'react'
import { Button } from '@/components/ui/button'
import {
  DropdownMenu,
  DropdownMenuContent,
  DropdownMenuItem,
  DropdownMenuTrigger,
} from '@/components/ui/dropdown-menu'
import { Download } from 'lucide-react'

interface ChartExportButtonProps {
  containerRef: RefObject<HTMLDivElement>
  filename?: string
}

function triggerDownload(dataUrl: string, filename: string) {
  const link = document.createElement('a')
  link.href = dataUrl
  link.download = filename
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)
}

function triggerBlobDownload(blob: Blob, filename: string) {
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  document.body.appendChild(link)
  link.click()
  setTimeout(() => {
    document.body.removeChild(link)
    URL.revokeObjectURL(url)
  }, 100)
}

function getPlotlyElement(container: HTMLDivElement): HTMLElement | null {
  return container.querySelector('.js-plotly-plot')
}

async function getPlotly() {
  const mod = await import('plotly.js')
  // Dynamic import of CJS module: the Plotly object may be at .default or the module itself
  return (mod as any).default || mod
}

async function exportPlotlyAsPng(plotEl: HTMLElement, filename: string) {
  const Plotly = await getPlotly()
  const dataUrl = await Plotly.toImage(plotEl, {
    format: 'png',
    width: 1200,
    height: 800,
    scale: 2,
  })
  triggerDownload(dataUrl, filename)
}

async function exportPlotlyAsSvg(plotEl: HTMLElement, filename: string) {
  const Plotly = await getPlotly()
  const dataUrl = await Plotly.toImage(plotEl, {
    format: 'svg',
    width: 1200,
    height: 800,
  })
  const response = await fetch(dataUrl)
  const blob = await response.blob()
  triggerBlobDownload(blob, filename)
}

async function exportPlotlyAsPdf(plotEl: HTMLElement, filename: string) {
  const [Plotly, { jsPDF }] = await Promise.all([
    getPlotly(),
    import('jspdf'),
  ])
  const dataUrl = await Plotly.toImage(plotEl, {
    format: 'png',
    width: 1200,
    height: 800,
    scale: 2,
  })
  const pdf = new jsPDF({ orientation: 'landscape', unit: 'px', format: [1200, 800] })
  pdf.addImage(dataUrl, 'PNG', 0, 0, 1200, 800)
  pdf.save(filename)
}

async function exportCanvasAsPng(container: HTMLDivElement, filename: string) {
  const html2canvas = (await import('html2canvas')).default
  const canvas = await html2canvas(container, { backgroundColor: '#ffffff', scale: 2 })
  const dataUrl = canvas.toDataURL('image/png')
  triggerDownload(dataUrl, filename)
}

async function exportCanvasAsSvg(container: HTMLDivElement, filename: string) {
  const html2canvas = (await import('html2canvas')).default
  const canvas = await html2canvas(container, { backgroundColor: '#ffffff', scale: 2 })
  const dataUrl = canvas.toDataURL('image/png')
  const w = canvas.width
  const h = canvas.height
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}">
  <image href="${dataUrl}" width="${w}" height="${h}" />
</svg>`
  const blob = new Blob([svg], { type: 'image/svg+xml' })
  triggerBlobDownload(blob, filename)
}

async function exportCanvasAsPdf(container: HTMLDivElement, filename: string) {
  const [html2canvasMod, { jsPDF }] = await Promise.all([
    import('html2canvas'),
    import('jspdf'),
  ])
  const html2canvas = html2canvasMod.default
  const canvas = await html2canvas(container, { backgroundColor: '#ffffff', scale: 2 })
  const dataUrl = canvas.toDataURL('image/png')
  const w = canvas.width / 2
  const h = canvas.height / 2
  const pdf = new jsPDF({
    orientation: w > h ? 'landscape' : 'portrait',
    unit: 'px',
    format: [w, h],
  })
  pdf.addImage(dataUrl, 'PNG', 0, 0, w, h)
  pdf.save(filename)
}

export function ChartExportButton({
  containerRef,
  filename = 'chart',
}: ChartExportButtonProps) {
  const handleExport = useCallback(
    async (format: 'png' | 'svg' | 'pdf') => {
      const container = containerRef.current
      if (!container) return

      const plotEl = getPlotlyElement(container)
      const name = `${filename}.${format}`

      try {
        if (plotEl) {
          switch (format) {
            case 'png':
              await exportPlotlyAsPng(plotEl, name)
              break
            case 'svg':
              await exportPlotlyAsSvg(plotEl, name)
              break
            case 'pdf':
              await exportPlotlyAsPdf(plotEl, name)
              break
          }
        } else {
          switch (format) {
            case 'png':
              await exportCanvasAsPng(container, name)
              break
            case 'svg':
              await exportCanvasAsSvg(container, name)
              break
            case 'pdf':
              await exportCanvasAsPdf(container, name)
              break
          }
        }
      } catch (err) {
        console.error(`Failed to export chart as ${format}:`, err)
      }
    },
    [containerRef, filename]
  )

  return (
    <DropdownMenu>
      <DropdownMenuTrigger asChild>
        <Button variant="outline" size="icon" className="h-8 w-8">
          <Download className="h-4 w-4" />
        </Button>
      </DropdownMenuTrigger>
      <DropdownMenuContent align="end">
        <DropdownMenuItem onClick={() => handleExport('png')}>
          Export as PNG
        </DropdownMenuItem>
        <DropdownMenuItem onClick={() => handleExport('svg')}>
          Export as SVG
        </DropdownMenuItem>
        <DropdownMenuItem onClick={() => handleExport('pdf')}>
          Export as PDF
        </DropdownMenuItem>
      </DropdownMenuContent>
    </DropdownMenu>
  )
}
