import { useState, useRef } from 'react'
import { useQuery } from '@tanstack/react-query'
import Plot from 'react-plotly.js'
import { diversityApi, samplesApi } from '@/api/endpoints'
import { useFilterStore } from '@/stores/filterStore'
import { useAuthStore } from '@/stores/authStore'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { ProjectFilter } from '@/components/filters/ProjectFilter'
import { SampleFilter } from '@/components/filters/SampleFilter'
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { ChartExportButton } from '@/components/charts/ChartExportButton'
import { Loader2, Download } from 'lucide-react'

export function DiversityPage() {
  const { projectId, sampleIds } = useFilterStore()
  const { accessToken } = useAuthStore()
  const shannonChartRef = useRef<HTMLDivElement>(null)
  const observedSpeciesChartRef = useRef<HTMLDivElement>(null)
  const betaChartRef = useRef<HTMLDivElement>(null)
  const pcoaChartRef = useRef<HTMLDivElement>(null)
  const [metric, setMetric] = useState('bray_curtis')
  const [pcoaDimensions, setPcoaDimensions] = useState<'2d' | '3d'>('2d')
  const [selectedDate, setSelectedDate] = useState<string>('__all__')
  const [selectedSubproject, setSelectedSubproject] = useState<string>('__all__')

  // Fetch dates for the selected project
  const { data: datesData } = useQuery({
    queryKey: ['sample-dates', projectId],
    queryFn: () => samplesApi.getDates(projectId || undefined),
    enabled: projectId !== null,
  })

  // Fetch subprojects for the selected project
  const { data: subprojectsData } = useQuery({
    queryKey: ['sample-subprojects', projectId],
    queryFn: () => samplesApi.getSubprojects(projectId || undefined),
    enabled: projectId !== null,
  })

  // Convert special value to undefined for API calls
  const dateFilter = selectedDate === '__all__' ? undefined : selectedDate
  const subprojectFilter = selectedSubproject === '__all__' ? undefined : selectedSubproject

  const { data: alphaData, isLoading: alphaLoading } = useQuery({
    queryKey: ['diversity-alpha', projectId, sampleIds, dateFilter, subprojectFilter],
    queryFn: () =>
      diversityApi.getAlpha({
        project_id: projectId || undefined,
        sample_id: sampleIds.length > 0 ? sampleIds : undefined,
        date: dateFilter,
        subproject: subprojectFilter,
      }),
    enabled: projectId !== null,
  })

  const { data: betaData, isLoading: betaLoading } = useQuery({
    queryKey: ['diversity-beta', projectId, sampleIds, metric, dateFilter, subprojectFilter],
    queryFn: () =>
      diversityApi.getBeta({
        project_id: projectId || undefined,
        sample_id: sampleIds.length > 0 ? sampleIds : undefined,
        date: dateFilter,
        subproject: subprojectFilter,
        metric,
      }),
    enabled: projectId !== null,
  })

  const { data: pcoaData, isLoading: pcoaLoading } = useQuery({
    queryKey: ['diversity-pcoa', projectId, sampleIds, metric, dateFilter, subprojectFilter],
    queryFn: () =>
      diversityApi.getPCoA({
        project_id: projectId || undefined,
        sample_id: sampleIds.length > 0 ? sampleIds : undefined,
        date: dateFilter,
        subproject: subprojectFilter,
        metric,
        n_components: 3,
      }),
    enabled: projectId !== null,
  })

  const dates = datesData?.dates || []
  const subprojects = subprojectsData?.subprojects || []

  // Export handlers
  const handleExportAlpha = () => {
    const url = diversityApi.getExportAlphaUrl({
      project_id: projectId || undefined,
      sample_id: sampleIds.length > 0 ? sampleIds : undefined,
      date: dateFilter,
      subproject: subprojectFilter,
    })
    fetch(url, {
      headers: {
        'Authorization': `Bearer ${accessToken}`,
      },
    })
      .then(response => response.blob())
      .then(blob => {
        const blobUrl = window.URL.createObjectURL(blob)
        const link = document.createElement('a')
        link.href = blobUrl
        link.download = 'alpha_diversity.csv'
        link.click()
        window.URL.revokeObjectURL(blobUrl)
      })
  }

  const handleExportBeta = () => {
    const url = diversityApi.getExportBetaUrl({
      project_id: projectId || undefined,
      sample_id: sampleIds.length > 0 ? sampleIds : undefined,
      date: dateFilter,
      subproject: subprojectFilter,
      metric,
    })
    fetch(url, {
      headers: {
        'Authorization': `Bearer ${accessToken}`,
      },
    })
      .then(response => response.blob())
      .then(blob => {
        const blobUrl = window.URL.createObjectURL(blob)
        const link = document.createElement('a')
        link.href = blobUrl
        link.download = `beta_diversity_${metric}.csv`
        link.click()
        window.URL.revokeObjectURL(blobUrl)
      })
  }

  const handleExportPCoA = () => {
    const url = diversityApi.getExportPCoAUrl({
      project_id: projectId || undefined,
      sample_id: sampleIds.length > 0 ? sampleIds : undefined,
      date: dateFilter,
      subproject: subprojectFilter,
      metric,
      n_components: 3,
    })
    fetch(url, {
      headers: {
        'Authorization': `Bearer ${accessToken}`,
      },
    })
      .then(response => response.blob())
      .then(blob => {
        const blobUrl = window.URL.createObjectURL(blob)
        const link = document.createElement('a')
        link.href = blobUrl
        link.download = 'pcoa_coordinates.csv'
        link.click()
        window.URL.revokeObjectURL(blobUrl)
      })
  }

  // Alpha diversity bar chart data
  const alphaBarData = alphaData
    ? [
        {
          type: 'bar' as const,
          name: 'Shannon',
          x: alphaData.map((d) => d.sample_name),
          y: alphaData.map((d) => d.shannon),
          marker: { color: '#3b82f6' },
        },
        {
          type: 'bar' as const,
          name: 'Simpson',
          x: alphaData.map((d) => d.sample_name),
          y: alphaData.map((d) => d.simpson),
          marker: { color: '#10b981' },
        },
      ]
    : []

  const observedSpeciesData = alphaData
    ? [
        {
          type: 'bar' as const,
          name: 'Observed Species',
          x: alphaData.map((d) => d.sample_name),
          y: alphaData.map((d) => d.observed_species),
          marker: { color: '#8b5cf6' },
        },
      ]
    : []

  // Beta diversity heatmap data
  const betaHeatmapData = betaData
    ? [
        {
          type: 'heatmap' as const,
          z: betaData.distance_matrix,
          x: betaData.sample_names,
          y: betaData.sample_names,
          colorscale: 'RdBu' as const,
          reversescale: true,
          zmin: 0,
          zmax: 1,
          colorbar: { title: { text: 'Distance' } },
          hovertemplate:
            '%{y} vs %{x}<br>Distance: %{z:.3f}<extra></extra>',
        },
      ]
    : []

  // PCoA scatter plot data
  const pcoaScatterData2D =
    pcoaData && pcoaData.coordinates
      ? [
          {
            type: 'scatter' as const,
            mode: 'text+markers' as const,
            x: pcoaData.coordinates.map((c) => c[0]),
            y: pcoaData.coordinates.map((c) => c[1]),
            text: pcoaData.sample_names,
            textposition: 'top center' as const,
            marker: {
              size: 12,
              color: '#3b82f6',
              line: { width: 1, color: '#1e40af' },
            },
            hovertemplate:
              '%{text}<br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<extra></extra>',
          },
        ]
      : []

  const pcoaScatterData3D =
    pcoaData && pcoaData.coordinates && pcoaData.coordinates[0]?.length >= 3
      ? [
          {
            type: 'scatter3d' as const,
            mode: 'text+markers' as const,
            x: pcoaData.coordinates.map((c) => c[0]),
            y: pcoaData.coordinates.map((c) => c[1]),
            z: pcoaData.coordinates.map((c) => c[2] || 0),
            text: pcoaData.sample_names,
            marker: {
              size: 8,
              color: '#3b82f6',
              line: { width: 1, color: '#1e40af' },
            },
            hovertemplate:
              '%{text}<br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<br>PC3: %{z:.3f}<extra></extra>',
          },
        ]
      : []

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-3xl font-bold">Diversity Analysis</h1>
        <p className="text-muted-foreground">
          Explore alpha and beta diversity metrics for your samples
        </p>
      </div>

      {/* Filters */}
      <div className="flex flex-wrap gap-4">
        <div className="w-48">
          <label className="mb-1 block text-sm font-medium">Project</label>
          <ProjectFilter />
        </div>
        <div className="w-64">
          <label className="mb-1 block text-sm font-medium">Samples</label>
          <SampleFilter />
        </div>
        <div className="w-48">
          <label className="mb-1 block text-sm font-medium">Date</label>
          <Select value={selectedDate} onValueChange={setSelectedDate}>
            <SelectTrigger>
              <SelectValue placeholder="All dates" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="__all__">All dates</SelectItem>
              {dates.map((date) => (
                <SelectItem key={date} value={date}>
                  {date}
                </SelectItem>
              ))}
            </SelectContent>
          </Select>
        </div>
        <div className="w-48">
          <label className="mb-1 block text-sm font-medium">Subproject</label>
          <Select value={selectedSubproject} onValueChange={setSelectedSubproject}>
            <SelectTrigger>
              <SelectValue placeholder="All subprojects" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="__all__">All subprojects</SelectItem>
              {subprojects.map((sp) => (
                <SelectItem key={sp} value={sp}>
                  {sp}
                </SelectItem>
              ))}
            </SelectContent>
          </Select>
        </div>
        <div className="w-40">
          <label className="mb-1 block text-sm font-medium">Distance Metric</label>
          <Select value={metric} onValueChange={setMetric}>
            <SelectTrigger>
              <SelectValue placeholder="Metric" />
            </SelectTrigger>
            <SelectContent>
              <SelectItem value="bray_curtis">Bray-Curtis</SelectItem>
              <SelectItem value="jaccard">Jaccard</SelectItem>
            </SelectContent>
          </Select>
        </div>
      </div>

      {!projectId ? (
        <Card>
          <CardContent className="flex h-96 items-center justify-center text-muted-foreground">
            Select a project to view diversity data
          </CardContent>
        </Card>
      ) : (
        <Tabs defaultValue="alpha" className="space-y-4">
          <TabsList>
            <TabsTrigger value="alpha">Alpha Diversity</TabsTrigger>
            <TabsTrigger value="beta">Beta Diversity</TabsTrigger>
            <TabsTrigger value="pcoa">PCoA</TabsTrigger>
          </TabsList>

          <TabsContent value="alpha" className="space-y-4">
            {/* Export button */}
            <div className="flex justify-end">
              <Button
                variant="outline"
                size="sm"
                onClick={handleExportAlpha}
                disabled={!projectId || !alphaData || alphaData.length === 0}
              >
                <Download className="mr-2 h-4 w-4" />
                Export Alpha Diversity
              </Button>
            </div>

            <div className="grid gap-4 lg:grid-cols-2">
              <Card>
                <CardHeader className="flex flex-row items-center justify-between">
                  <CardTitle>Shannon & Simpson Diversity</CardTitle>
                  <ChartExportButton containerRef={shannonChartRef} filename="shannon-simpson-diversity" />
                </CardHeader>
                <CardContent>
                  <div ref={shannonChartRef}>
                    {alphaLoading ? (
                      <div className="flex h-80 items-center justify-center">
                        <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                      </div>
                    ) : alphaBarData.length > 0 ? (
                      <Plot
                        data={alphaBarData}
                        layout={{
                          barmode: 'group',
                          xaxis: { title: { text: 'Sample' }, tickangle: -45, automargin: true },
                          yaxis: { title: { text: 'Index Value' } },
                          margin: { t: 30, b: 100, l: 60, r: 30 },
                          height: 350,
                          autosize: true,
                          legend: { orientation: 'h' as const, y: 1.1 },
                        }}
                        config={{ responsive: true, displayModeBar: true }}
                        style={{ width: '100%' }}
                      />
                    ) : (
                      <div className="flex h-80 items-center justify-center text-muted-foreground">
                        No data available
                      </div>
                    )}
                  </div>
                </CardContent>
              </Card>

              <Card>
                <CardHeader className="flex flex-row items-center justify-between">
                  <CardTitle>Observed Species</CardTitle>
                  <ChartExportButton containerRef={observedSpeciesChartRef} filename="observed-species" />
                </CardHeader>
                <CardContent>
                  <div ref={observedSpeciesChartRef}>
                    {alphaLoading ? (
                      <div className="flex h-80 items-center justify-center">
                        <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                      </div>
                    ) : observedSpeciesData.length > 0 ? (
                      <Plot
                        data={observedSpeciesData}
                        layout={{
                          xaxis: { title: { text: 'Sample' }, tickangle: -45, automargin: true },
                          yaxis: { title: { text: 'Number of Species' } },
                          margin: { t: 30, b: 100, l: 60, r: 30 },
                          height: 350,
                          autosize: true,
                        }}
                        config={{ responsive: true, displayModeBar: true }}
                        style={{ width: '100%' }}
                      />
                    ) : (
                      <div className="flex h-80 items-center justify-center text-muted-foreground">
                        No data available
                      </div>
                    )}
                  </div>
                </CardContent>
              </Card>
            </div>

            {/* Alpha diversity table */}
            <Card>
              <CardHeader>
                <CardTitle>Alpha Diversity Metrics Table</CardTitle>
              </CardHeader>
              <CardContent>
                {alphaData && alphaData.length > 0 ? (
                  <div className="overflow-x-auto">
                    <table className="w-full text-sm">
                      <thead>
                        <tr className="border-b">
                          <th className="px-4 py-2 text-left">Sample</th>
                          <th className="px-4 py-2 text-right">Shannon</th>
                          <th className="px-4 py-2 text-right">Simpson</th>
                          <th className="px-4 py-2 text-right">Observed Species</th>
                          <th className="px-4 py-2 text-right">Pielou Evenness</th>
                        </tr>
                      </thead>
                      <tbody>
                        {alphaData.map((item) => (
                          <tr key={item.sample_id} className="border-b hover:bg-muted/50">
                            <td className="px-4 py-2">{item.sample_name}</td>
                            <td className="px-4 py-2 text-right">{item.shannon.toFixed(3)}</td>
                            <td className="px-4 py-2 text-right">{item.simpson.toFixed(3)}</td>
                            <td className="px-4 py-2 text-right">{item.observed_species}</td>
                            <td className="px-4 py-2 text-right">
                              {item.pielou_evenness?.toFixed(3) || 'N/A'}
                            </td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                ) : (
                  <div className="flex h-32 items-center justify-center text-muted-foreground">
                    No alpha diversity data available
                  </div>
                )}
              </CardContent>
            </Card>
          </TabsContent>

          <TabsContent value="beta">
            {/* Export button */}
            <div className="mb-4 flex justify-end">
              <Button
                variant="outline"
                size="sm"
                onClick={handleExportBeta}
                disabled={!projectId || !betaData || !betaData.distance_matrix}
              >
                <Download className="mr-2 h-4 w-4" />
                Export Distance Matrix
              </Button>
            </div>

            <Card>
              <CardHeader className="flex flex-row items-center justify-between">
                <CardTitle>
                  Beta Diversity Heatmap ({metric === 'bray_curtis' ? 'Bray-Curtis' : 'Jaccard'})
                </CardTitle>
                <ChartExportButton containerRef={betaChartRef} filename="beta-diversity-heatmap" />
              </CardHeader>
              <CardContent>
                <div ref={betaChartRef}>
                  {betaLoading ? (
                    <div className="flex h-96 items-center justify-center">
                      <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                    </div>
                  ) : betaHeatmapData.length > 0 ? (
                    <Plot
                      data={betaHeatmapData}
                      layout={{
                        xaxis: { title: { text: 'Sample' }, tickangle: -45, automargin: true },
                        yaxis: { title: { text: 'Sample' }, automargin: true },
                        margin: { t: 30, b: 120, l: 120, r: 80 },
                        height: 500,
                        autosize: true,
                      }}
                      config={{ responsive: true, displayModeBar: true }}
                      style={{ width: '100%' }}
                    />
                  ) : (
                    <div className="flex h-96 items-center justify-center text-muted-foreground">
                      No beta diversity data available (need at least 2 samples)
                    </div>
                  )}
                </div>
              </CardContent>
            </Card>
          </TabsContent>

          <TabsContent value="pcoa">
            {/* Export button */}
            <div className="mb-4 flex justify-end">
              <Button
                variant="outline"
                size="sm"
                onClick={handleExportPCoA}
                disabled={!projectId || !pcoaData || !pcoaData.coordinates}
              >
                <Download className="mr-2 h-4 w-4" />
                Export PCoA Coordinates
              </Button>
            </div>

            <Card>
              <CardHeader className="flex flex-row items-center justify-between">
                <CardTitle>
                  Principal Coordinates Analysis (PCoA)
                  {pcoaData?.explained_variance && (
                    <span className="ml-2 text-sm font-normal text-muted-foreground">
                      PC1: {pcoaData.explained_variance[0]?.toFixed(1)}%,
                      PC2: {pcoaData.explained_variance[1]?.toFixed(1)}%
                      {pcoaDimensions === '3d' && pcoaData.explained_variance[2] && (
                        <>, PC3: {pcoaData.explained_variance[2]?.toFixed(1)}%</>
                      )}
                    </span>
                  )}
                </CardTitle>
                <div className="flex items-center gap-2">
                  <Select
                    value={pcoaDimensions}
                    onValueChange={(v) => setPcoaDimensions(v as '2d' | '3d')}
                  >
                    <SelectTrigger className="w-24">
                      <SelectValue />
                    </SelectTrigger>
                    <SelectContent>
                      <SelectItem value="2d">2D</SelectItem>
                      <SelectItem value="3d">3D</SelectItem>
                    </SelectContent>
                  </Select>
                  <ChartExportButton containerRef={pcoaChartRef} filename="pcoa-plot" />
                </div>
              </CardHeader>
              <CardContent>
                <div ref={pcoaChartRef}>
                  {pcoaLoading ? (
                    <div className="flex h-96 items-center justify-center">
                      <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                    </div>
                  ) : pcoaDimensions === '2d' && pcoaScatterData2D.length > 0 ? (
                    <Plot
                      data={pcoaScatterData2D}
                      layout={{
                        xaxis: {
                          title: {
                            text: `PC1 (${pcoaData?.explained_variance?.[0]?.toFixed(1) || 0}%)`,
                          },
                        },
                        yaxis: {
                          title: {
                            text: `PC2 (${pcoaData?.explained_variance?.[1]?.toFixed(1) || 0}%)`,
                          },
                        },
                        margin: { t: 30, b: 60, l: 60, r: 30 },
                        height: 500,
                        autosize: true,
                        hovermode: 'closest' as const,
                      }}
                      config={{ responsive: true, displayModeBar: true }}
                      style={{ width: '100%' }}
                    />
                  ) : pcoaDimensions === '3d' && pcoaScatterData3D.length > 0 ? (
                    <Plot
                      data={pcoaScatterData3D}
                      layout={{
                        scene: {
                          xaxis: {
                            title: { text: `PC1 (${pcoaData?.explained_variance?.[0]?.toFixed(1) || 0}%)` },
                          },
                          yaxis: {
                            title: { text: `PC2 (${pcoaData?.explained_variance?.[1]?.toFixed(1) || 0}%)` },
                          },
                          zaxis: {
                            title: { text: `PC3 (${pcoaData?.explained_variance?.[2]?.toFixed(1) || 0}%)` },
                          },
                        },
                        margin: { t: 30, b: 30, l: 30, r: 30 },
                        height: 550,
                        autosize: true,
                      }}
                      config={{ responsive: true, displayModeBar: true }}
                      style={{ width: '100%' }}
                    />
                  ) : (
                    <div className="flex h-96 items-center justify-center text-muted-foreground">
                      No PCoA data available (need at least 3 samples)
                    </div>
                  )}
                </div>
              </CardContent>
            </Card>
          </TabsContent>
        </Tabs>
      )}
    </div>
  )
}
