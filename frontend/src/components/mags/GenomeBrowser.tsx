import { useEffect, useRef, useState } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'
import { AlertCircle, RefreshCw } from 'lucide-react'
import { magsApi } from '@/api/endpoints'
import type { MAGContig } from '@/types'

interface GenomeBrowserProps {
  magId: number
  magName: string
  contigs: MAGContig[]
}

// IGV.js is loaded dynamically
declare global {
  interface Window {
    igv: {
      createBrowser: (container: HTMLElement, options: object) => Promise<unknown>
      removeBrowser: (browser: unknown) => void
    }
  }
}

export function GenomeBrowser({ magId, magName, contigs }: GenomeBrowserProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const browserRef = useRef<unknown>(null)
  const [selectedContig, setSelectedContig] = useState<string>('')
  const [isLoading, setIsLoading] = useState(false)
  const [error, setError] = useState<string | null>(null)
  const [igvLoaded, setIgvLoaded] = useState(false)

  // Load IGV.js dynamically
  useEffect(() => {
    if (window.igv) {
      setIgvLoaded(true)
      return
    }

    const script = document.createElement('script')
    script.src = 'https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js'
    script.async = true
    script.onload = () => setIgvLoaded(true)
    script.onerror = () => setError('Failed to load IGV.js genome browser')
    document.head.appendChild(script)

    return () => {
      // Cleanup on unmount
      if (browserRef.current && window.igv) {
        try {
          window.igv.removeBrowser(browserRef.current)
        } catch {
          // Ignore errors during cleanup
        }
      }
    }
  }, [])

  // Initialize browser when IGV is loaded
  useEffect(() => {
    if (!igvLoaded || !containerRef.current || !magId) return

    const initBrowser = async () => {
      setIsLoading(true)
      setError(null)

      try {
        // Remove existing browser if any
        if (browserRef.current) {
          window.igv.removeBrowser(browserRef.current)
          browserRef.current = null
        }

        const fastaUrl = magsApi.getFastaUrl(magId)
        const gffUrl = magsApi.getGffUrl(magId)

        const options = {
          genome: {
            id: magName,
            name: magName,
            fastaURL: fastaUrl,
            indexURL: magsApi.getFastaIndexUrl(magId),
          },
          tracks: [
            {
              name: 'Genes',
              type: 'annotation',
              format: 'gff3',
              url: gffUrl,
              displayMode: 'EXPANDED',
              height: 200,
              color: '#3B82F6',
            },
          ],
          locus: selectedContig || (contigs.length > 0 ? contigs[0].contig_id : undefined),
        }

        browserRef.current = await window.igv.createBrowser(containerRef.current!, options)
      } catch (err) {
        console.error('Failed to initialize genome browser:', err)
        setError('Failed to initialize genome browser. The MAG files may not be available.')
      } finally {
        setIsLoading(false)
      }
    }

    initBrowser()
  }, [igvLoaded, magId, magName])

  // Navigate to contig when selected
  useEffect(() => {
    if (!browserRef.current || !selectedContig) return

    try {
      // IGV browser search/navigation
      const browser = browserRef.current as { search: (locus: string) => void }
      browser.search(selectedContig)
    } catch {
      // Ignore navigation errors
    }
  }, [selectedContig])

  const handleRefresh = () => {
    if (browserRef.current && window.igv) {
      window.igv.removeBrowser(browserRef.current)
      browserRef.current = null
    }
    setIgvLoaded(false)
    // Re-trigger the effect
    setTimeout(() => setIgvLoaded(true), 100)
  }

  if (error) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>Genome Browser</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="rounded-md border border-red-200 bg-red-50 p-4 text-red-800">
            <div className="flex items-center gap-2">
              <AlertCircle className="h-4 w-4" />
              <span>{error}</span>
            </div>
          </div>
          <Button onClick={handleRefresh} className="mt-4" variant="outline">
            <RefreshCw className="mr-2 h-4 w-4" />
            Retry
          </Button>
        </CardContent>
      </Card>
    )
  }

  return (
    <Card>
      <CardHeader>
        <div className="flex items-center justify-between">
          <CardTitle>Genome Browser</CardTitle>
          <div className="flex items-center gap-2">
            <Select value={selectedContig} onValueChange={setSelectedContig}>
              <SelectTrigger className="w-64">
                <SelectValue placeholder="Jump to contig..." />
              </SelectTrigger>
              <SelectContent>
                {contigs.map((contig) => (
                  <SelectItem key={contig.contig_id} value={contig.contig_id}>
                    {contig.contig_id} ({contig.length.toLocaleString()} bp)
                  </SelectItem>
                ))}
              </SelectContent>
            </Select>
            <Button variant="outline" size="icon" onClick={handleRefresh}>
              <RefreshCw className="h-4 w-4" />
            </Button>
          </div>
        </div>
      </CardHeader>
      <CardContent>
        {isLoading && (
          <div className="flex h-96 items-center justify-center">
            <div className="animate-pulse text-muted-foreground">Loading genome browser...</div>
          </div>
        )}
        <div
          ref={containerRef}
          className="min-h-[500px] rounded-lg border"
          style={{ display: isLoading ? 'none' : 'block' }}
        />
        {!isLoading && contigs.length === 0 && (
          <div className="flex h-96 items-center justify-center text-muted-foreground">
            No contig data available for this MAG
          </div>
        )}
      </CardContent>
    </Card>
  )
}
