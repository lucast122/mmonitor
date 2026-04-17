import { useState, useRef } from 'react'
import { useQuery } from '@tanstack/react-query'
import { taxonomyApi, horizonApi, samplesApi } from '@/api/endpoints'
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
import { Input } from '@/components/ui/input'
import { StackedBarChart } from '@/components/charts/StackedBarChart'
import { HorizonChart } from '@/components/charts/HorizonChart'
import { ChartExportButton } from '@/components/charts/ChartExportButton'
import { Loader2, Download } from 'lucide-react'

type PlotType = 'stacked_bar' | 'grouped_bar' | 'heatmap' | 'area'

export function TaxonomyPage() {
  const { projectId, sampleIds } = useFilterStore()
  const { accessToken } = useAuthStore()
  const taxonomyChartRef = useRef<HTMLDivElement>(null)
  const horizonChartRef = useRef<HTMLDivElement>(null)
  const [rank, setRank] = useState('tax_phylum')
  const [topN, setTopN] = useState(20)
  const [plotType, setPlotType] = useState<PlotType>('stacked_bar')
  const [selectedDate, setSelectedDate] = useState<string>('__all__')
  const [selectedSubproject, setSelectedSubproject] = useState<string>('__all__')
  const [horizonTaxaCount, setHorizonTaxaCount] = useState(20)
  const [horizonBands, setHorizonBands] = useState(4)
  const [horizonMode, setHorizonMode] = useState<'mirror' | 'offset'>('mirror')
  const [activeTab, setActiveTab] = useState('taxonomy')

  const { data: ranks } = useQuery({
    queryKey: ['taxonomy-ranks'],
    queryFn: taxonomyApi.getRanks,
  })

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

  const { data: taxonomyData, isLoading: taxonomyLoading } = useQuery({
    queryKey: ['taxonomy-aggregated', projectId, sampleIds, rank, topN, dateFilter, subprojectFilter],
    queryFn: () =>
      taxonomyApi.getAggregated({
        project_id: projectId || undefined,
        sample_id: sampleIds.length > 0 ? sampleIds : undefined,
        date: dateFilter,
        subproject: subprojectFilter,
        rank,
        top_n: topN,
        group_by: 'sample_id',
      }),
    enabled: projectId !== null,
  })

  // Fetch horizon data
  const { data: horizonData, isLoading: horizonLoading } = useQuery({
    queryKey: ['horizon-data', projectId, subprojectFilter, horizonTaxaCount],
    queryFn: () =>
      horizonApi.getData({
        project_id: projectId || undefined,
        subproject: subprojectFilter,
        taxa_count: horizonTaxaCount,
      }),
    enabled: projectId !== null && activeTab === 'horizon',
  })

  const handleDownloadRawCounts = () => {
    const url = taxonomyApi.getCountTableUrl({
      project_id: projectId || undefined,
      sample_id: sampleIds.length > 0 ? sampleIds : undefined,
      date: dateFilter,
      subproject: subprojectFilter,
      rank,
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
        link.download = 'raw_counts.csv'
        link.click()
        window.URL.revokeObjectURL(blobUrl)
      })
  }

  const handleDownloadRelativeAbundances = () => {
    const url = taxonomyApi.getNormalizedTableUrl({
      project_id: projectId || undefined,
      sample_id: sampleIds.length > 0 ? sampleIds : undefined,
      date: dateFilter,
      subproject: subprojectFilter,
      rank,
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
        link.download = 'relative_abundances.csv'
        link.click()
        window.URL.revokeObjectURL(blobUrl)
      })
  }

  const dates = datesData?.dates || []
  const subprojects = subprojectsData?.subprojects || []

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-3xl font-bold">Taxonomy Analysis</h1>
        <p className="text-muted-foreground">
          Visualize taxonomic composition and abundance changes across your samples
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
            Select a project to view taxonomy data
          </CardContent>
        </Card>
      ) : (
        <Tabs value={activeTab} onValueChange={setActiveTab}>
          <TabsList>
            <TabsTrigger value="taxonomy">Taxonomy Plot</TabsTrigger>
            <TabsTrigger value="horizon">Horizon Plot</TabsTrigger>
          </TabsList>

          <TabsContent value="taxonomy" className="space-y-4">
            {/* Taxonomy-specific controls */}
            <div className="flex flex-wrap items-end gap-4">
              <div className="w-48">
                <label className="mb-1 block text-sm font-medium">Taxonomic Rank</label>
                <Select value={rank} onValueChange={setRank}>
                  <SelectTrigger>
                    <SelectValue placeholder="Select rank" />
                  </SelectTrigger>
                  <SelectContent>
                    {ranks?.ranks.map((r) => (
                      <SelectItem key={r.value} value={r.value}>
                        {r.label}
                      </SelectItem>
                    ))}
                  </SelectContent>
                </Select>
              </div>
              <div className="w-32">
                <label className="mb-1 block text-sm font-medium">Top N Taxa</label>
                <Select value={topN.toString()} onValueChange={(v) => setTopN(parseInt(v))}>
                  <SelectTrigger>
                    <SelectValue placeholder="Top N" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="10">Top 10</SelectItem>
                    <SelectItem value="20">Top 20</SelectItem>
                    <SelectItem value="30">Top 30</SelectItem>
                    <SelectItem value="50">Top 50</SelectItem>
                  </SelectContent>
                </Select>
              </div>
              <div className="w-40">
                <label className="mb-1 block text-sm font-medium">Plot Type</label>
                <Select value={plotType} onValueChange={(v) => setPlotType(v as PlotType)}>
                  <SelectTrigger>
                    <SelectValue placeholder="Plot type" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="stacked_bar">Stacked Bar</SelectItem>
                    <SelectItem value="grouped_bar">Grouped Bar</SelectItem>
                    <SelectItem value="heatmap">Heatmap</SelectItem>
                    <SelectItem value="area">Area Chart</SelectItem>
                  </SelectContent>
                </Select>
              </div>
              <div className="flex gap-2">
                <Button
                  variant="outline"
                  size="sm"
                  onClick={handleDownloadRawCounts}
                  disabled={!projectId}
                >
                  <Download className="mr-2 h-4 w-4" />
                  Raw Counts
                </Button>
                <Button
                  variant="outline"
                  size="sm"
                  onClick={handleDownloadRelativeAbundances}
                  disabled={!projectId}
                >
                  <Download className="mr-2 h-4 w-4" />
                  Relative Abundances
                </Button>
              </div>
            </div>

            <Card>
              <CardHeader className="flex flex-row items-center justify-between">
                <CardTitle>Taxonomic Composition</CardTitle>
                <ChartExportButton containerRef={taxonomyChartRef} filename="taxonomic-composition" />
              </CardHeader>
              <CardContent>
                <div ref={taxonomyChartRef}>
                  {taxonomyLoading ? (
                    <div className="flex h-96 items-center justify-center">
                      <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                    </div>
                  ) : taxonomyData && Array.isArray(taxonomyData) && taxonomyData.length > 0 ? (
                    <StackedBarChart data={taxonomyData} plotType={plotType} />
                  ) : (
                    <div className="flex h-96 items-center justify-center text-muted-foreground">
                      No taxonomy data available
                    </div>
                  )}
                </div>
              </CardContent>
            </Card>
          </TabsContent>

          <TabsContent value="horizon" className="space-y-4">
            {/* Horizon-specific controls */}
            <div className="flex flex-wrap items-end gap-4">
              <div className="w-32">
                <label className="mb-1 block text-sm font-medium">Taxa Count</label>
                <Input
                  type="number"
                  min={5}
                  max={50}
                  value={horizonTaxaCount}
                  onChange={(e) => setHorizonTaxaCount(parseInt(e.target.value) || 20)}
                />
              </div>
              <div className="w-32">
                <label className="mb-1 block text-sm font-medium">Bands</label>
                <Select value={horizonBands.toString()} onValueChange={(v) => setHorizonBands(parseInt(v))}>
                  <SelectTrigger>
                    <SelectValue placeholder="Bands" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="2">2 bands</SelectItem>
                    <SelectItem value="3">3 bands</SelectItem>
                    <SelectItem value="4">4 bands</SelectItem>
                    <SelectItem value="5">5 bands</SelectItem>
                    <SelectItem value="6">6 bands</SelectItem>
                  </SelectContent>
                </Select>
              </div>
              <div className="w-32">
                <label className="mb-1 block text-sm font-medium">Mode</label>
                <Select value={horizonMode} onValueChange={(v) => setHorizonMode(v as 'mirror' | 'offset')}>
                  <SelectTrigger>
                    <SelectValue placeholder="Mode" />
                  </SelectTrigger>
                  <SelectContent>
                    <SelectItem value="mirror">Mirror</SelectItem>
                    <SelectItem value="offset">Offset</SelectItem>
                  </SelectContent>
                </Select>
              </div>
            </div>

            <Card>
              <CardHeader className="flex flex-row items-center justify-between">
                <CardTitle>Abundance Deviation from Mean Over Time</CardTitle>
                <ChartExportButton containerRef={horizonChartRef} filename="horizon-chart" />
              </CardHeader>
              <CardContent>
                <p className="mb-4 text-sm text-muted-foreground">
                  Each row shows a taxon's abundance deviation from its mean as a horizon chart.
                  <span className="ml-1" style={{ color: 'midnightblue' }}>Blue</span> indicates below-average abundance,
                  <span className="ml-1" style={{ color: 'crimson' }}>red</span> indicates above-average.
                  The chart uses {horizonBands} overlapping bands - darker colors represent larger deviations.
                </p>
                <div ref={horizonChartRef}>
                  {horizonLoading ? (
                    <div className="flex h-96 items-center justify-center">
                      <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                    </div>
                  ) : horizonData && horizonData.taxa.length > 0 ? (
                    <HorizonChart data={horizonData} bands={horizonBands} mode={horizonMode} />
                  ) : (
                    <div className="flex h-96 items-center justify-center text-muted-foreground">
                      No horizon data available. Need multiple dates/samples to show trends.
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
