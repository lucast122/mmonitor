import { useState, useRef, useEffect } from 'react'
import { useQuery, useMutation, useQueryClient } from '@tanstack/react-query'
import Plot from 'react-plotly.js'
import { correlationsApi, samplesApi } from '@/api/endpoints'
import { useFilterStore } from '@/stores/filterStore'
import { useAuthStore } from '@/stores/authStore'
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { ProjectFilter } from '@/components/filters/ProjectFilter'
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { ChartExportButton } from '@/components/charts/ChartExportButton'
import { Loader2, Upload, Download, Plus, Save, Trash2 } from 'lucide-react'
import { toast } from 'sonner'

interface MetadataRow {
  sample_id: string
  sample_name: string
  [key: string]: string | number | null
}

export function CorrelationsPage() {
  const { projectId } = useFilterStore()
  const { accessToken } = useAuthStore()
  const queryClient = useQueryClient()
  const fileInputRef = useRef<HTMLInputElement>(null)
  const correlationChartRef = useRef<HTMLDivElement>(null)

  const [method, setMethod] = useState<'pearson' | 'spearman' | 'kendall'>('pearson')
  const [pCutoff, setPCutoff] = useState(0.05)
  const [topNTaxa, setTopNTaxa] = useState(20)
  const [activeTab, setActiveTab] = useState('correlations')

  // Metadata editor state
  const [metadataRows, setMetadataRows] = useState<MetadataRow[]>([])
  const [metadataColumns, setMetadataColumns] = useState<string[]>([])
  const [newColumnName, setNewColumnName] = useState('')

  // Fetch samples for the project
  const { data: samplesData } = useQuery({
    queryKey: ['samples', projectId],
    queryFn: () => samplesApi.list(projectId || undefined),
    enabled: projectId !== null,
  })

  // Fetch existing metadata
  const { data: existingMetadata, isLoading: metadataLoading } = useQuery({
    queryKey: ['correlations-metadata', projectId],
    queryFn: () => correlationsApi.getMetadata(projectId || undefined),
    enabled: projectId !== null,
  })

  // Initialize metadata rows when data changes
  useEffect(() => {
    if (existingMetadata && existingMetadata.length > 0) {
      // Extract columns from metadata
      const cols = new Set<string>()
      existingMetadata.forEach((row) => {
        Object.keys(row).forEach((key) => {
          if (key !== 'sample_id' && key !== 'sample_name') {
            cols.add(key)
          }
        })
      })
      setMetadataColumns(Array.from(cols).sort())
      // Ensure each row has sample_id and sample_name
      const normalizedRows = existingMetadata.map((row) => ({
        sample_id: String(row.sample_id || ''),
        sample_name: String(row.sample_name || row.sample_id || ''),
        ...Object.fromEntries(
          Object.entries(row).filter(([key]) => key !== 'sample_id' && key !== 'sample_name')
        ),
      })) as MetadataRow[]
      setMetadataRows(normalizedRows)
    } else if (samplesData && samplesData.length > 0) {
      // Initialize with samples but no metadata columns
      const rows = samplesData.map((s) => ({
        sample_id: String(s.id),
        sample_name: s.name || String(s.id),
      }))
      setMetadataRows(rows)
      setMetadataColumns([])
    } else {
      // Reset when no data
      setMetadataRows([])
      setMetadataColumns([])
    }
  }, [existingMetadata, samplesData])

  // Fetch correlation results
  const {
    data: correlationData,
    isLoading: correlationLoading,
    error: correlationError,
  } = useQuery({
    queryKey: ['correlations-compute', projectId, method, pCutoff, topNTaxa],
    queryFn: () =>
      correlationsApi.compute({
        project_id: projectId || undefined,
        method,
        p_cutoff: pCutoff,
        top_n_taxa: topNTaxa,
      }),
    enabled: projectId !== null && activeTab === 'correlations',
    retry: false,
  })

  // Upload mutation
  const uploadMutation = useMutation({
    mutationFn: (file: File) => correlationsApi.uploadMetadata(projectId || undefined, file),
    onSuccess: (data) => {
      toast.success(`Metadata uploaded: ${data.created} created, ${data.updated} updated`)
      queryClient.invalidateQueries({ queryKey: ['correlations-metadata'] })
      queryClient.invalidateQueries({ queryKey: ['correlations-compute'] })
    },
    onError: (error: Error) => {
      toast.error(`Upload failed: ${error.message}`)
    },
  })

  // Save metadata mutation
  const saveMutation = useMutation({
    mutationFn: () => correlationsApi.saveMetadata(projectId || undefined, metadataRows),
    onSuccess: (data) => {
      toast.success(`Metadata saved: ${data.created} created, ${data.updated} updated`)
      queryClient.invalidateQueries({ queryKey: ['correlations-metadata'] })
      queryClient.invalidateQueries({ queryKey: ['correlations-compute'] })
    },
    onError: (error: Error) => {
      toast.error(`Save failed: ${error.message}`)
    },
  })

  const handleFileUpload = (event: React.ChangeEvent<HTMLInputElement>) => {
    const file = event.target.files?.[0]
    if (file) {
      uploadMutation.mutate(file)
    }
  }

  const handleAddColumn = () => {
    if (newColumnName && !metadataColumns.includes(newColumnName)) {
      setMetadataColumns([...metadataColumns, newColumnName])
      // Add the column to all rows with null value
      setMetadataRows(metadataRows.map((row) => ({ ...row, [newColumnName]: null })))
      setNewColumnName('')
    }
  }

  const handleRemoveColumn = (column: string) => {
    setMetadataColumns(metadataColumns.filter((c) => c !== column))
    setMetadataRows(
      metadataRows.map((row) => {
        const newRow = { ...row }
        delete newRow[column]
        return newRow
      })
    )
  }

  const handleCellChange = (rowIndex: number, column: string, value: string) => {
    const newRows = [...metadataRows]
    // Try to parse as number
    const numValue = parseFloat(value)
    newRows[rowIndex] = {
      ...newRows[rowIndex],
      [column]: isNaN(numValue) ? value : numValue,
    }
    setMetadataRows(newRows)
  }

  const handleExport = () => {
    const url = correlationsApi.getExportUrl({
      project_id: projectId || undefined,
      method,
      top_n_taxa: topNTaxa,
    })
    fetch(url, {
      headers: {
        Authorization: `Bearer ${accessToken}`,
      },
    })
      .then((response) => response.blob())
      .then((blob) => {
        const blobUrl = window.URL.createObjectURL(blob)
        const link = document.createElement('a')
        link.href = blobUrl
        link.download = `correlations_${method}.csv`
        link.click()
        window.URL.revokeObjectURL(blobUrl)
      })
  }

  // Build heatmap data
  const heatmapData =
    correlationData && correlationData.correlations
      ? [
          {
            type: 'heatmap' as const,
            z: correlationData.correlations,
            x: correlationData.metadata_fields,
            y: correlationData.taxa,
            colorscale: 'RdBu' as const,
            reversescale: true,
            zmin: -1,
            zmax: 1,
            colorbar: { title: { text: 'Correlation' } },
            hovertemplate:
              'Taxon: %{y}<br>Metadata: %{x}<br>Correlation: %{z:.3f}<extra></extra>',
          },
        ]
      : []

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-3xl font-bold">Correlation Analysis</h1>
        <p className="text-muted-foreground">
          Analyze correlations between taxonomic abundances and sample metadata
        </p>
      </div>

      {/* Filters */}
      <div className="flex flex-wrap gap-4">
        <div className="w-48">
          <label className="mb-1 block text-sm font-medium">Project</label>
          <ProjectFilter />
        </div>
        <div className="w-40">
          <label className="mb-1 block text-sm font-medium">Method</label>
          <Select value={method} onValueChange={(v) => setMethod(v as typeof method)}>
            <SelectTrigger>
              <SelectValue placeholder="Method" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="pearson">Pearson</SelectItem>
              <SelectItem value="spearman">Spearman</SelectItem>
              <SelectItem value="kendall">Kendall</SelectItem>
            </SelectContent>
          </Select>
        </div>
        <div className="w-32">
          <label className="mb-1 block text-sm font-medium">P-value Cutoff</label>
          <Input
            type="number"
            min={0}
            max={1}
            step={0.01}
            value={pCutoff}
            onChange={(e) => setPCutoff(parseFloat(e.target.value) || 0.05)}
          />
        </div>
        <div className="w-32">
          <label className="mb-1 block text-sm font-medium">Top N Taxa</label>
          <Input
            type="number"
            min={5}
            max={100}
            value={topNTaxa}
            onChange={(e) => setTopNTaxa(parseInt(e.target.value) || 20)}
          />
        </div>
      </div>

      {!projectId ? (
        <Card>
          <CardContent className="flex h-96 items-center justify-center text-muted-foreground">
            Select a project to view correlations
          </CardContent>
        </Card>
      ) : (
        <Tabs value={activeTab} onValueChange={setActiveTab} className="space-y-4">
          <TabsList>
            <TabsTrigger value="correlations">Correlation Heatmap</TabsTrigger>
            <TabsTrigger value="metadata">Metadata Editor</TabsTrigger>
          </TabsList>

          <TabsContent value="correlations" className="space-y-4">
            {/* Export button */}
            <div className="flex justify-end">
              <Button
                variant="outline"
                size="sm"
                onClick={handleExport}
                disabled={!correlationData || !correlationData.correlations}
              >
                <Download className="mr-2 h-4 w-4" />
                Export Correlations
              </Button>
            </div>

            <Card>
              <CardHeader className="flex flex-row items-center justify-between">
                <div>
                  <CardTitle>
                    Taxonomy-Metadata Correlations ({method.charAt(0).toUpperCase() + method.slice(1)})
                  </CardTitle>
                  {correlationData && (
                    <CardDescription>
                      Based on {correlationData.n_samples} samples with both taxonomy and metadata
                    </CardDescription>
                  )}
                </div>
                <ChartExportButton containerRef={correlationChartRef} filename="correlation-heatmap" />
              </CardHeader>
              <CardContent>
                <div ref={correlationChartRef}>
                  {correlationLoading ? (
                    <div className="flex h-96 items-center justify-center">
                      <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                    </div>
                  ) : correlationError ? (
                    <div className="flex h-96 flex-col items-center justify-center text-muted-foreground">
                      <p className="mb-4">
                        {(correlationError as Error).message || 'Failed to compute correlations'}
                      </p>
                      <p className="text-sm">
                        Make sure you have uploaded sample metadata in the "Metadata Editor" tab.
                      </p>
                    </div>
                  ) : heatmapData.length > 0 ? (
                    <Plot
                      data={heatmapData}
                      layout={{
                        xaxis: {
                          title: { text: 'Metadata Field' },
                          tickangle: -45,
                          automargin: true,
                        },
                        yaxis: {
                          title: { text: 'Taxon' },
                          automargin: true,
                        },
                        margin: { t: 30, b: 120, l: 200, r: 100 },
                        height: Math.max(500, (correlationData?.taxa.length || 10) * 25 + 150),
                        autosize: true,
                      }}
                      config={{ responsive: true, displayModeBar: true }}
                      style={{ width: '100%' }}
                    />
                  ) : (
                    <div className="flex h-96 flex-col items-center justify-center text-muted-foreground">
                      <p>No correlation data available.</p>
                      <p className="mt-2 text-sm">
                        Upload sample metadata in the "Metadata Editor" tab to compute correlations.
                      </p>
                    </div>
                  )}
                </div>
              </CardContent>
            </Card>
          </TabsContent>

          <TabsContent value="metadata" className="space-y-4">
            <Card>
              <CardHeader>
                <CardTitle>Sample Metadata</CardTitle>
                <CardDescription>
                  Upload a CSV file or edit metadata inline. The first column should be sample_id.
                  Only numeric metadata fields will be used for correlation analysis.
                </CardDescription>
              </CardHeader>
              <CardContent className="space-y-4">
                {/* Upload section */}
                <div className="flex items-center gap-4">
                  <input
                    type="file"
                    ref={fileInputRef}
                    accept=".csv"
                    onChange={handleFileUpload}
                    className="hidden"
                  />
                  <Button
                    variant="outline"
                    onClick={() => fileInputRef.current?.click()}
                    disabled={uploadMutation.isPending}
                  >
                    {uploadMutation.isPending ? (
                      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                    ) : (
                      <Upload className="mr-2 h-4 w-4" />
                    )}
                    Upload CSV
                  </Button>
                  <span className="text-sm text-muted-foreground">
                    CSV format: sample_id, field1, field2, ...
                  </span>
                </div>

                {/* Add column */}
                <div className="flex items-center gap-2">
                  <Input
                    placeholder="New column name"
                    value={newColumnName}
                    onChange={(e) => setNewColumnName(e.target.value)}
                    className="w-48"
                  />
                  <Button variant="outline" size="sm" onClick={handleAddColumn}>
                    <Plus className="mr-2 h-4 w-4" />
                    Add Column
                  </Button>
                </div>

                {/* Metadata table */}
                {metadataLoading ? (
                  <div className="flex h-48 items-center justify-center">
                    <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                  </div>
                ) : metadataRows.length > 0 ? (
                  <div className="max-h-96 overflow-auto rounded border">
                    <table className="w-full text-sm">
                      <thead className="sticky top-0 bg-muted">
                        <tr>
                          <th className="px-3 py-2 text-left font-medium">Sample</th>
                          {metadataColumns.map((col) => (
                            <th key={col} className="px-3 py-2 text-left font-medium">
                              <div className="flex items-center gap-2">
                                {col}
                                <button
                                  onClick={() => handleRemoveColumn(col)}
                                  className="text-muted-foreground hover:text-destructive"
                                >
                                  <Trash2 className="h-3 w-3" />
                                </button>
                              </div>
                            </th>
                          ))}
                        </tr>
                      </thead>
                      <tbody>
                        {metadataRows.map((row, rowIndex) => (
                          <tr key={row.sample_id} className="border-t">
                            <td className="px-3 py-2 font-medium">{row.sample_name}</td>
                            {metadataColumns.map((col) => (
                              <td key={col} className="px-3 py-2">
                                <Input
                                  type="text"
                                  value={row[col]?.toString() || ''}
                                  onChange={(e) => handleCellChange(rowIndex, col, e.target.value)}
                                  className="h-8 w-24"
                                  placeholder="—"
                                />
                              </td>
                            ))}
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                ) : (
                  <div className="flex h-48 items-center justify-center text-muted-foreground">
                    No samples found. Select a project with samples.
                  </div>
                )}

                {/* Save button */}
                {metadataRows.length > 0 && metadataColumns.length > 0 && (
                  <div className="flex justify-end">
                    <Button onClick={() => saveMutation.mutate()} disabled={saveMutation.isPending}>
                      {saveMutation.isPending ? (
                        <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                      ) : (
                        <Save className="mr-2 h-4 w-4" />
                      )}
                      Save Metadata
                    </Button>
                  </div>
                )}
              </CardContent>
            </Card>
          </TabsContent>
        </Tabs>
      )}
    </div>
  )
}
