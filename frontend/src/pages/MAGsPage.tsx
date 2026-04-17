import { useState } from 'react'
import { useQuery } from '@tanstack/react-query'
import { magsApi } from '@/api/endpoints'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { Label } from '@/components/ui/label'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { Loader2, Dna, BarChart3, Search, Eye } from 'lucide-react'
import {
  COGDistributionChart,
  FunctionalSummary,
  GeneSearchTable,
  GenomeBrowser,
} from '@/components/mags'
import type { MAG } from '@/types'

export function MAGsPage() {
  const [minCompleteness, setMinCompleteness] = useState(50)
  const [maxContamination, setMaxContamination] = useState(10)
  const [selectedMag, setSelectedMag] = useState<MAG | null>(null)
  const [activeTab, setActiveTab] = useState('overview')

  const { data: magsData, isLoading: magsLoading } = useQuery({
    queryKey: ['mags', minCompleteness, maxContamination],
    queryFn: () =>
      magsApi.list({
        min_completeness: minCompleteness,
        max_contamination: maxContamination,
      }),
  })

  const { data: magDetail, isLoading: detailLoading } = useQuery({
    queryKey: ['magDetail', selectedMag?.id],
    queryFn: () => magsApi.get(selectedMag!.id),
    enabled: !!selectedMag,
  })

  const { data: contigs, isLoading: contigsLoading } = useQuery({
    queryKey: ['magContigs', selectedMag?.id],
    queryFn: () => magsApi.getContigsSummary(selectedMag!.id),
    enabled: !!selectedMag,
  })

  const { data: cogDistribution, isLoading: cogLoading } = useQuery({
    queryKey: ['magCogDistribution', selectedMag?.id],
    queryFn: () => magsApi.getCogDistribution(selectedMag!.id),
    enabled: !!selectedMag && activeTab === 'functional',
  })

  const { data: pfamDomains } = useQuery({
    queryKey: ['magPfamDomains', selectedMag?.id],
    queryFn: () => magsApi.getPfamDomains(selectedMag!.id, 20),
    enabled: !!selectedMag && activeTab === 'functional',
  })

  const mags = magsData?.results || []

  const formatSize = (bytes: number | null) => {
    if (!bytes) return 'N/A'
    if (bytes >= 1000000) return `${(bytes / 1000000).toFixed(2)} Mb`
    if (bytes >= 1000) return `${(bytes / 1000).toFixed(2)} Kb`
    return `${bytes} bp`
  }

  return (
    <div className="space-y-6">
      <div>
        <h1 className="text-3xl font-bold">MAG Analysis</h1>
        <p className="text-muted-foreground">
          Explore metagenome-assembled genomes with functional annotations and genome browser
        </p>
      </div>

      <div className="flex flex-wrap gap-4">
        <div className="space-y-2">
          <Label htmlFor="minCompleteness">Min Completeness (%)</Label>
          <Input
            id="minCompleteness"
            type="number"
            min={0}
            max={100}
            value={minCompleteness}
            onChange={(e) => setMinCompleteness(Number(e.target.value))}
            className="w-32"
          />
        </div>
        <div className="space-y-2">
          <Label htmlFor="maxContamination">Max Contamination (%)</Label>
          <Input
            id="maxContamination"
            type="number"
            min={0}
            max={100}
            value={maxContamination}
            onChange={(e) => setMaxContamination(Number(e.target.value))}
            className="w-32"
          />
        </div>
      </div>

      <div className="grid gap-4 lg:grid-cols-3">
        {/* MAG List */}
        <Card className="lg:col-span-1">
          <CardHeader>
            <CardTitle>MAG List ({mags.length})</CardTitle>
          </CardHeader>
          <CardContent>
            {magsLoading ? (
              <div className="flex h-96 items-center justify-center">
                <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
              </div>
            ) : mags.length > 0 ? (
              <div className="max-h-[600px] space-y-2 overflow-y-auto">
                {mags.map((mag) => (
                  <div
                    key={mag.id}
                    className={`cursor-pointer rounded-lg border p-3 transition-colors hover:bg-muted ${
                      selectedMag?.id === mag.id ? 'border-primary bg-muted' : ''
                    }`}
                    onClick={() => setSelectedMag(mag)}
                  >
                    <div className="font-medium">{mag.name}</div>
                    <div className="text-sm text-muted-foreground">
                      Sample: {mag.sample_name}
                    </div>
                    <div className="mt-1 flex gap-4 text-xs">
                      <span
                        className={
                          (mag.completeness || 0) >= 90
                            ? 'text-green-600'
                            : (mag.completeness || 0) >= 50
                              ? 'text-yellow-600'
                              : 'text-red-600'
                        }
                      >
                        Comp: {mag.completeness?.toFixed(1) || 'N/A'}%
                      </span>
                      <span
                        className={
                          (mag.contamination || 0) <= 5
                            ? 'text-green-600'
                            : (mag.contamination || 0) <= 10
                              ? 'text-yellow-600'
                              : 'text-red-600'
                        }
                      >
                        Cont: {mag.contamination?.toFixed(1) || 'N/A'}%
                      </span>
                    </div>
                    {mag.genome_size && (
                      <div className="mt-1 text-xs text-muted-foreground">
                        Size: {formatSize(mag.genome_size)} | Contigs: {mag.n_contigs || 'N/A'}
                      </div>
                    )}
                  </div>
                ))}
              </div>
            ) : (
              <div className="flex h-48 items-center justify-center text-muted-foreground">
                No MAGs found matching criteria
              </div>
            )}
          </CardContent>
        </Card>

        {/* MAG Details with Tabs */}
        <div className="lg:col-span-2">
          {selectedMag ? (
            <Tabs value={activeTab} onValueChange={setActiveTab}>
              <TabsList className="mb-4 grid w-full grid-cols-4">
                <TabsTrigger value="overview" className="flex items-center gap-2">
                  <Eye className="h-4 w-4" />
                  Overview
                </TabsTrigger>
                <TabsTrigger value="functional" className="flex items-center gap-2">
                  <BarChart3 className="h-4 w-4" />
                  Functional
                </TabsTrigger>
                <TabsTrigger value="browser" className="flex items-center gap-2">
                  <Dna className="h-4 w-4" />
                  Browser
                </TabsTrigger>
                <TabsTrigger value="search" className="flex items-center gap-2">
                  <Search className="h-4 w-4" />
                  Gene Search
                </TabsTrigger>
              </TabsList>

              {/* Overview Tab */}
              <TabsContent value="overview">
                <Card>
                  <CardHeader>
                    <CardTitle>{selectedMag.name}</CardTitle>
                    <p className="text-sm text-muted-foreground">
                      Sample: {selectedMag.sample_name}
                    </p>
                  </CardHeader>
                  <CardContent>
                    {detailLoading ? (
                      <div className="flex h-48 items-center justify-center">
                        <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                      </div>
                    ) : (
                      <div className="space-y-6">
                        {/* Quality Metrics */}
                        <div>
                          <h4 className="mb-3 font-medium">Quality Metrics</h4>
                          <div className="grid grid-cols-2 gap-4 md:grid-cols-4">
                            <div className="rounded-lg border p-3">
                              <div className="text-sm text-muted-foreground">Completeness</div>
                              <div className="text-2xl font-bold">
                                {selectedMag.completeness?.toFixed(1) || 'N/A'}%
                              </div>
                            </div>
                            <div className="rounded-lg border p-3">
                              <div className="text-sm text-muted-foreground">Contamination</div>
                              <div className="text-2xl font-bold">
                                {selectedMag.contamination?.toFixed(1) || 'N/A'}%
                              </div>
                            </div>
                            <div className="rounded-lg border p-3">
                              <div className="text-sm text-muted-foreground">Strain Het.</div>
                              <div className="text-2xl font-bold">
                                {selectedMag.strain_heterogeneity?.toFixed(1) || 'N/A'}%
                              </div>
                            </div>
                            <div className="rounded-lg border p-3">
                              <div className="text-sm text-muted-foreground">GC Content</div>
                              <div className="text-2xl font-bold">
                                {selectedMag.gc_content?.toFixed(1) || 'N/A'}%
                              </div>
                            </div>
                          </div>
                        </div>

                        {/* Assembly Statistics */}
                        <div>
                          <h4 className="mb-3 font-medium">Assembly Statistics</h4>
                          <div className="grid grid-cols-2 gap-4 md:grid-cols-4">
                            <div className="rounded-lg border p-3">
                              <div className="text-sm text-muted-foreground">Genome Size</div>
                              <div className="text-xl font-bold">
                                {formatSize(selectedMag.genome_size)}
                              </div>
                            </div>
                            <div className="rounded-lg border p-3">
                              <div className="text-sm text-muted-foreground">Contigs</div>
                              <div className="text-xl font-bold">
                                {selectedMag.n_contigs?.toLocaleString() || 'N/A'}
                              </div>
                            </div>
                            <div className="rounded-lg border p-3">
                              <div className="text-sm text-muted-foreground">N50</div>
                              <div className="text-xl font-bold">
                                {formatSize(selectedMag.n50)}
                              </div>
                            </div>
                            <div className="rounded-lg border p-3">
                              <div className="text-sm text-muted-foreground">Gene Count</div>
                              <div className="text-xl font-bold">
                                {selectedMag.gene_count?.toLocaleString() || 'N/A'}
                              </div>
                            </div>
                          </div>
                        </div>

                        {/* Gene Counts */}
                        {magDetail && (
                          <div>
                            <h4 className="mb-3 font-medium">Gene Counts</h4>
                            <div className="grid grid-cols-3 gap-4">
                              <div className="rounded-lg border p-3">
                                <div className="text-sm text-muted-foreground">CDS</div>
                                <div className="text-xl font-bold">
                                  {magDetail.cds_count?.toLocaleString() || 'N/A'}
                                </div>
                              </div>
                              <div className="rounded-lg border p-3">
                                <div className="text-sm text-muted-foreground">tRNA</div>
                                <div className="text-xl font-bold">
                                  {magDetail.trna_count || 'N/A'}
                                </div>
                              </div>
                              <div className="rounded-lg border p-3">
                                <div className="text-sm text-muted-foreground">rRNA</div>
                                <div className="text-xl font-bold">
                                  {magDetail.rrna_count || 'N/A'}
                                </div>
                              </div>
                            </div>
                          </div>
                        )}

                        {/* Taxonomy */}
                        {selectedMag.taxonomy && Object.keys(selectedMag.taxonomy).length > 0 && (
                          <div>
                            <h4 className="mb-3 font-medium">Taxonomy</h4>
                            <div className="rounded-lg border p-3">
                              {Object.entries(selectedMag.taxonomy).map(([rank, value]) => (
                                <div key={rank} className="flex justify-between py-1">
                                  <span className="text-muted-foreground capitalize">{rank}:</span>
                                  <span className="font-medium">{value}</span>
                                </div>
                              ))}
                            </div>
                          </div>
                        )}
                      </div>
                    )}
                  </CardContent>
                </Card>
              </TabsContent>

              {/* Functional Tab */}
              <TabsContent value="functional" className="space-y-4">
                <FunctionalSummary
                  summary={magDetail?.functional_summary || null}
                  isLoading={detailLoading}
                />
                <COGDistributionChart
                  data={cogDistribution || []}
                  isLoading={cogLoading}
                />
                {pfamDomains && pfamDomains.length > 0 && (
                  <Card>
                    <CardHeader>
                      <CardTitle>Top PFAM Domains</CardTitle>
                    </CardHeader>
                    <CardContent>
                      <div className="space-y-2">
                        {pfamDomains.map((domain) => (
                          <div
                            key={domain.pfam}
                            className="flex items-center justify-between rounded border p-2"
                          >
                            <div>
                              <span className="font-mono text-sm">{domain.pfam}</span>
                              {domain.description && (
                                <span className="ml-2 text-sm text-muted-foreground">
                                  {domain.description}
                                </span>
                              )}
                            </div>
                            <span className="font-semibold">{domain.count}</span>
                          </div>
                        ))}
                      </div>
                    </CardContent>
                  </Card>
                )}
              </TabsContent>

              {/* Browser Tab */}
              <TabsContent value="browser">
                {contigsLoading ? (
                  <div className="flex h-96 items-center justify-center">
                    <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
                  </div>
                ) : (
                  <GenomeBrowser
                    magId={selectedMag.id}
                    magName={selectedMag.name}
                    contigs={contigs || []}
                  />
                )}
              </TabsContent>

              {/* Gene Search Tab */}
              <TabsContent value="search">
                <GeneSearchTable magId={selectedMag.id} />
              </TabsContent>
            </Tabs>
          ) : (
            <Card>
              <CardContent className="flex h-96 items-center justify-center text-muted-foreground">
                Select a MAG from the list to view details
              </CardContent>
            </Card>
          )}
        </div>
      </div>
    </div>
  )
}
