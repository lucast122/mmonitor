import { useState, useRef } from 'react'
import { useQuery } from '@tanstack/react-query'
import Plot from 'react-plotly.js'
import { qcApi, samplesApi } from '@/api/endpoints'
import { useFilterStore } from '@/stores/filterStore'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
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
import { Loader2 } from 'lucide-react'
import { formatNumber } from '@/lib/utils'

export function QCPage() {
  const { projectId, sampleIds } = useFilterStore()
  const q20q30ChartRef = useRef<HTMLDivElement>(null)
  const readsChartRef = useRef<HTMLDivElement>(null)
  const qualityChartRef = useRef<HTMLDivElement>(null)
  const readLengthChartRef = useRef<HTMLDivElement>(null)
  const gcChartRef = useRef<HTMLDivElement>(null)
  const [activeTab, setActiveTab] = useState('overview')
  const [selectedDate, setSelectedDate] = useState<string>('__all__')
  const [selectedSubproject, setSelectedSubproject] = useState<string>('__all__')

  // Convert special value to undefined for API calls
  const dateFilter = selectedDate === '__all__' ? undefined : selectedDate
  const subprojectFilter = selectedSubproject === '__all__' ? undefined : selectedSubproject

  // Fetch dates and subprojects
  const { data: datesData } = useQuery({
    queryKey: ['sample-dates', projectId],
    queryFn: () => samplesApi.getDates(projectId || undefined),
    enabled: projectId !== null,
  })

  const { data: subprojectsData } = useQuery({
    queryKey: ['sample-subprojects', projectId],
    queryFn: () => samplesApi.getSubprojects(projectId || undefined),
    enabled: projectId !== null,
  })

  // Fetch QC list for table
  const { data: qcData, isLoading: qcLoading } = useQuery({
    queryKey: ['qc', projectId, dateFilter, subprojectFilter],
    queryFn: () => qcApi.list({
      project_id: projectId || undefined,
      date: dateFilter,
      subproject: subprojectFilter,
    }),
    enabled: projectId !== null,
  })

  // Fetch QC summary
  const { data: summary, isLoading: summaryLoading } = useQuery({
    queryKey: ['qc-summary', projectId, dateFilter, subprojectFilter],
    queryFn: () => qcApi.getSummary({
      project_id: projectId || undefined,
      date: dateFilter,
      subproject: subprojectFilter,
    }),
    enabled: projectId !== null,
  })

  // Fetch metrics list for plots
  const { data: metricsData, isLoading: metricsLoading } = useQuery({
    queryKey: ['qc-metrics-list', projectId, dateFilter, subprojectFilter],
    queryFn: () => qcApi.getMetricsList({
      project_id: projectId || undefined,
      date: dateFilter,
      subproject: subprojectFilter,
    }),
    enabled: projectId !== null,
  })

  const isLoading = qcLoading || summaryLoading || metricsLoading
  const qcList = qcData ? (Array.isArray(qcData) ? qcData : qcData?.results || []) : []
  const dates = datesData?.dates || []
  const subprojects = subprojectsData?.subprojects || []

  // Safe access to metrics data - ensure all arrays exist
  const safeMetricsData = metricsData &&
    metricsData.samples &&
    Array.isArray(metricsData.samples) &&
    Array.isArray(metricsData.mean_quality_scores) &&
    Array.isArray(metricsData.mean_read_lengths) &&
    Array.isArray(metricsData.mean_gc_contents) &&
    Array.isArray(metricsData.q20_scores) &&
    Array.isArray(metricsData.q30_scores) &&
    Array.isArray(metricsData.number_of_reads)
      ? metricsData
      : null

  // Filter by selected samples and deduplicate (keep latest entry per sample)
  const deduplicatedQcList = (() => {
    const latestBySample = new Map<string, typeof qcList[0]>()
    for (const qc of qcList) {
      // Later entries overwrite earlier ones, keeping the latest
      latestBySample.set(qc.sample_name, qc)
    }
    return Array.from(latestBySample.values())
  })()
  const filteredQcList = sampleIds.length > 0
    ? deduplicatedQcList.filter((qc) => sampleIds.includes(qc.sample_name))
    : deduplicatedQcList

  // Helper to safely filter paired arrays and deduplicate (keep latest entry per sample)
  const filterPairedArrays = (
    samples: string[],
    values: (number | null)[]
  ): { x: string[]; y: number[] } => {
    // Use a Map to keep only the latest value for each sample
    const latestValues = new Map<string, number>()
    for (let i = 0; i < Math.min(samples.length, values.length); i++) {
      if (values[i] !== null && values[i] !== undefined) {
        // Later entries overwrite earlier ones, keeping the latest
        latestValues.set(samples[i], values[i] as number)
      }
    }
    // Convert back to arrays
    const x: string[] = []
    const y: number[] = []
    for (const [sample, value] of latestValues) {
      x.push(sample)
      y.push(value)
    }
    return { x, y }
  }

  // Create bar chart data for quality scores (one value per sample)
  const qualityFiltered = safeMetricsData
    ? filterPairedArrays(safeMetricsData.samples, safeMetricsData.mean_quality_scores)
    : { x: [], y: [] }
  const qualityBarData = qualityFiltered.y.length > 0 ? [{
    type: 'bar' as const,
    x: qualityFiltered.x,
    y: qualityFiltered.y,
    name: 'Quality Score',
    marker: { color: '#3b82f6' },
    hovertemplate: 'Sample: %{x}<br>Quality: %{y:.2f}<extra></extra>',
  }] : []

  // Create bar chart data for read lengths (one value per sample)
  const readLengthFiltered = safeMetricsData
    ? filterPairedArrays(safeMetricsData.samples, safeMetricsData.mean_read_lengths)
    : { x: [], y: [] }
  const readLengthBarData = readLengthFiltered.y.length > 0 ? [{
    type: 'bar' as const,
    x: readLengthFiltered.x,
    y: readLengthFiltered.y,
    name: 'Read Length',
    marker: { color: '#10b981' },
    hovertemplate: 'Sample: %{x}<br>Length: %{y:.0f} bp<extra></extra>',
  }] : []

  // Create bar chart data for GC content (one value per sample)
  const gcFiltered = safeMetricsData
    ? filterPairedArrays(safeMetricsData.samples, safeMetricsData.mean_gc_contents)
    : { x: [], y: [] }
  const gcBarData = gcFiltered.y.length > 0 ? [{
    type: 'bar' as const,
    x: gcFiltered.x,
    y: gcFiltered.y,
    name: 'GC Content',
    marker: { color: '#ef4444' },
    hovertemplate: 'Sample: %{x}<br>GC: %{y:.2f}%<extra></extra>',
  }] : []

  // Create bar chart for Q20/Q30 scores (deduplicated)
  const q20Filtered = safeMetricsData
    ? filterPairedArrays(safeMetricsData.samples, safeMetricsData.q20_scores)
    : { x: [], y: [] }
  const q30Filtered = safeMetricsData
    ? filterPairedArrays(safeMetricsData.samples, safeMetricsData.q30_scores)
    : { x: [], y: [] }
  const qualityScoresBarData = q20Filtered.y.length > 0 ? [
    {
      type: 'bar' as const,
      name: 'Q20',
      x: q20Filtered.x,
      y: q20Filtered.y,
      marker: { color: '#3b82f6' },
      hovertemplate: 'Sample: %{x}<br>Q20: %{y:.2f}%<extra></extra>',
    },
    {
      type: 'bar' as const,
      name: 'Q30',
      x: q30Filtered.x,
      y: q30Filtered.y,
      marker: { color: '#10b981' },
      hovertemplate: 'Sample: %{x}<br>Q30: %{y:.2f}%<extra></extra>',
    },
  ] : []

  // Create bar chart for number of reads (deduplicated)
  const readsFiltered = safeMetricsData
    ? filterPairedArrays(safeMetricsData.samples, safeMetricsData.number_of_reads)
    : { x: [], y: [] }
  const readsBarData = readsFiltered.y.length > 0 ? [{
    type: 'bar' as const,
    name: 'Number of Reads',
    x: readsFiltered.x,
    y: readsFiltered.y,
    marker: { color: '#8b5cf6' },
    hovertemplate: 'Sample: %{x}<br>Reads: %{y:,.0f}<extra></extra>',
  }] : []

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-3xl font-bold">Quality Control</h1>
        <p className="text-muted-foreground">
          Analyze sequencing quality metrics across your samples
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
      </div>

      {!projectId ? (
        <Card>
          <CardContent className="flex h-96 items-center justify-center text-muted-foreground">
            Select a project to view QC data
          </CardContent>
        </Card>
      ) : isLoading ? (
        <Card>
          <CardContent className="flex h-96 items-center justify-center">
            <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
          </CardContent>
        </Card>
      ) : (
        <>
          {/* Summary Cards */}
          {summary && (
            <div className="grid gap-4 md:grid-cols-3 lg:grid-cols-6">
              <Card>
                <CardHeader className="pb-2">
                  <CardTitle className="text-sm font-medium text-muted-foreground">
                    Total Samples
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="text-2xl font-bold">{summary.total_samples}</div>
                </CardContent>
              </Card>
              <Card>
                <CardHeader className="pb-2">
                  <CardTitle className="text-sm font-medium text-muted-foreground">
                    Total Reads
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="text-2xl font-bold">
                    {formatNumber(summary.total_reads) || 'N/A'}
                  </div>
                </CardContent>
              </Card>
              <Card>
                <CardHeader className="pb-2">
                  <CardTitle className="text-sm font-medium text-muted-foreground">
                    Total Bases
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="text-2xl font-bold">
                    {formatNumber(summary.total_bases) || 'N/A'}
                  </div>
                </CardContent>
              </Card>
              <Card>
                <CardHeader className="pb-2">
                  <CardTitle className="text-sm font-medium text-muted-foreground">
                    Avg Quality
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="text-2xl font-bold">
                    {summary.avg_quality_score?.toFixed(1) || 'N/A'}
                  </div>
                </CardContent>
              </Card>
              <Card>
                <CardHeader className="pb-2">
                  <CardTitle className="text-sm font-medium text-muted-foreground">
                    Avg Read Length
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="text-2xl font-bold">
                    {summary.avg_read_length?.toFixed(0) || 'N/A'} bp
                  </div>
                </CardContent>
              </Card>
              <Card>
                <CardHeader className="pb-2">
                  <CardTitle className="text-sm font-medium text-muted-foreground">
                    Avg GC Content
                  </CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="text-2xl font-bold">
                    {summary.avg_gc_content?.toFixed(1) || 'N/A'}%
                  </div>
                </CardContent>
              </Card>
            </div>
          )}

          {/* Plots */}
          <Tabs value={activeTab} onValueChange={setActiveTab}>
            <TabsList>
              <TabsTrigger value="overview">Overview</TabsTrigger>
              <TabsTrigger value="quality">Quality Scores</TabsTrigger>
              <TabsTrigger value="lengths">Read Lengths</TabsTrigger>
              <TabsTrigger value="gc">GC Content</TabsTrigger>
              <TabsTrigger value="table">Data Table</TabsTrigger>
            </TabsList>

            <TabsContent value="overview" className="space-y-4">
              <div className="grid gap-4 lg:grid-cols-2">
                <Card>
                  <CardHeader className="flex flex-row items-center justify-between">
                    <CardTitle>Q20/Q30 Scores by Sample</CardTitle>
                    <ChartExportButton containerRef={q20q30ChartRef} filename="q20-q30-scores" />
                  </CardHeader>
                  <CardContent>
                    <div ref={q20q30ChartRef}>
                      {qualityScoresBarData.length > 0 ? (
                        <Plot
                          data={qualityScoresBarData}
                          layout={{
                            barmode: 'group',
                            xaxis: { title: { text: 'Sample' }, tickangle: -45 },
                            yaxis: { title: { text: 'Score (%)' }, range: [0, 100] },
                            margin: { t: 30, b: 100, l: 60, r: 30 },
                            height: 400,
                            autosize: true,
                            legend: { orientation: 'h' as const, y: 1.1 },
                          }}
                          config={{ responsive: true, displayModeBar: true }}
                          style={{ width: '100%' }}
                        />
                      ) : (
                        <div className="flex h-64 items-center justify-center text-muted-foreground">
                          No data available
                        </div>
                      )}
                    </div>
                  </CardContent>
                </Card>

                <Card>
                  <CardHeader className="flex flex-row items-center justify-between">
                    <CardTitle>Number of Reads by Sample</CardTitle>
                    <ChartExportButton containerRef={readsChartRef} filename="number-of-reads" />
                  </CardHeader>
                  <CardContent>
                    <div ref={readsChartRef}>
                      {readsBarData.length > 0 ? (
                        <Plot
                          data={readsBarData}
                          layout={{
                            xaxis: { title: { text: 'Sample' }, tickangle: -45 },
                            yaxis: { title: { text: 'Number of Reads' } },
                            margin: { t: 30, b: 100, l: 80, r: 30 },
                            height: 400,
                            autosize: true,
                          }}
                          config={{ responsive: true, displayModeBar: true }}
                          style={{ width: '100%' }}
                        />
                      ) : (
                        <div className="flex h-64 items-center justify-center text-muted-foreground">
                          No data available
                        </div>
                      )}
                    </div>
                  </CardContent>
                </Card>
              </div>
            </TabsContent>

            <TabsContent value="quality">
              <Card>
                <CardHeader className="flex flex-row items-center justify-between">
                  <CardTitle>Mean Quality Scores by Sample</CardTitle>
                  <ChartExportButton containerRef={qualityChartRef} filename="mean-quality-scores" />
                </CardHeader>
                <CardContent>
                  <div ref={qualityChartRef}>
                    {qualityBarData.length > 0 && qualityBarData[0].y.length > 0 ? (
                      <Plot
                        data={qualityBarData}
                        layout={{
                          xaxis: { title: { text: 'Sample' }, tickangle: -45 },
                          yaxis: { title: { text: 'Mean Quality Score' } },
                          margin: { t: 30, b: 100, l: 60, r: 30 },
                          height: 500,
                          autosize: true,
                          showlegend: false,
                        }}
                        config={{ responsive: true, displayModeBar: true }}
                        style={{ width: '100%' }}
                      />
                    ) : (
                      <div className="flex h-96 items-center justify-center text-muted-foreground">
                        No quality score data available
                      </div>
                    )}
                  </div>
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="lengths">
              <Card>
                <CardHeader className="flex flex-row items-center justify-between">
                  <CardTitle>Mean Read Lengths by Sample</CardTitle>
                  <ChartExportButton containerRef={readLengthChartRef} filename="mean-read-lengths" />
                </CardHeader>
                <CardContent>
                  <div ref={readLengthChartRef}>
                    {readLengthBarData.length > 0 && readLengthBarData[0].y.length > 0 ? (
                      <Plot
                        data={readLengthBarData}
                        layout={{
                          xaxis: { title: { text: 'Sample' }, tickangle: -45 },
                          yaxis: { title: { text: 'Mean Read Length (bp)' } },
                          margin: { t: 30, b: 100, l: 60, r: 30 },
                          height: 500,
                          autosize: true,
                          showlegend: false,
                        }}
                        config={{ responsive: true, displayModeBar: true }}
                        style={{ width: '100%' }}
                      />
                    ) : (
                      <div className="flex h-96 items-center justify-center text-muted-foreground">
                        No read length data available
                      </div>
                    )}
                  </div>
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="gc">
              <Card>
                <CardHeader className="flex flex-row items-center justify-between">
                  <CardTitle>Mean GC Content by Sample</CardTitle>
                  <ChartExportButton containerRef={gcChartRef} filename="mean-gc-content" />
                </CardHeader>
                <CardContent>
                  <div ref={gcChartRef}>
                    {gcBarData.length > 0 && gcBarData[0].y.length > 0 ? (
                      <Plot
                        data={gcBarData}
                        layout={{
                          xaxis: { title: { text: 'Sample' }, tickangle: -45 },
                          yaxis: { title: { text: 'Mean GC Content (%)' } },
                          margin: { t: 30, b: 100, l: 60, r: 30 },
                          height: 500,
                          autosize: true,
                          showlegend: false,
                        }}
                        config={{ responsive: true, displayModeBar: true }}
                        style={{ width: '100%' }}
                      />
                    ) : (
                      <div className="flex h-96 items-center justify-center text-muted-foreground">
                        No GC content data available
                      </div>
                    )}
                  </div>
                </CardContent>
              </Card>
            </TabsContent>

            <TabsContent value="table">
              <Card>
                <CardHeader>
                  <CardTitle>Sample Quality Metrics</CardTitle>
                </CardHeader>
                <CardContent>
                  {filteredQcList.length > 0 ? (
                    <div className="overflow-x-auto">
                      <table className="w-full text-sm">
                        <thead>
                          <tr className="border-b">
                            <th className="px-4 py-2 text-left">Sample</th>
                            <th className="px-4 py-2 text-right">Reads</th>
                            <th className="px-4 py-2 text-right">Mean Length</th>
                            <th className="px-4 py-2 text-right">Quality</th>
                            <th className="px-4 py-2 text-right">GC %</th>
                            <th className="px-4 py-2 text-right">Q20</th>
                            <th className="px-4 py-2 text-right">Q30</th>
                          </tr>
                        </thead>
                        <tbody>
                          {filteredQcList.map((qc) => (
                            <tr key={qc.id} className="border-b hover:bg-muted/50">
                              <td className="px-4 py-2">{qc.sample_name}</td>
                              <td className="px-4 py-2 text-right">
                                {formatNumber(qc.number_of_reads || 0)}
                              </td>
                              <td className="px-4 py-2 text-right">
                                {qc.mean_read_length?.toFixed(0) || 'N/A'}
                              </td>
                              <td className="px-4 py-2 text-right">
                                {qc.mean_quality_score?.toFixed(1) || 'N/A'}
                              </td>
                              <td className="px-4 py-2 text-right">
                                {qc.mean_gc_content?.toFixed(1) || 'N/A'}
                              </td>
                              <td className="px-4 py-2 text-right">
                                {qc.q20_score?.toFixed(1) || 'N/A'}
                              </td>
                              <td className="px-4 py-2 text-right">
                                {qc.q30_score?.toFixed(1) || 'N/A'}
                              </td>
                            </tr>
                          ))}
                        </tbody>
                      </table>
                    </div>
                  ) : (
                    <div className="flex h-48 items-center justify-center text-muted-foreground">
                      No QC data available
                    </div>
                  )}
                </CardContent>
              </Card>
            </TabsContent>
          </Tabs>
        </>
      )}
    </div>
  )
}
