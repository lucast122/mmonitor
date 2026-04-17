import { useQuery } from '@tanstack/react-query'
import { projectsApi, qcApi } from '@/api/endpoints'
import { useFilterStore } from '@/stores/filterStore'
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from '@/components/ui/card'
import { ProjectFilter } from '@/components/filters/ProjectFilter'
import { Dna, BarChart3, Layers, Microscope, Loader2 } from 'lucide-react'
import { formatNumber } from '@/lib/utils'
import { Link } from 'react-router-dom'

export function DashboardPage() {
  const { projectId } = useFilterStore()

  const { data: projects, isLoading: projectsLoading } = useQuery({
    queryKey: ['projects'],
    queryFn: projectsApi.list,
  })

  const { data: qcSummary, isLoading: qcLoading } = useQuery({
    queryKey: ['qc-summary', projectId],
    queryFn: () => qcApi.getSummary({ project_id: projectId?.toString() }),
  })

  const stats = [
    {
      title: 'Total Samples',
      value: qcSummary?.total_samples || 0,
      icon: Dna,
      href: '/taxonomy',
    },
    {
      title: 'Total Reads',
      value: formatNumber(qcSummary?.total_reads || 0),
      icon: BarChart3,
      href: '/qc',
    },
    {
      title: 'Avg Quality Score',
      value: qcSummary?.avg_quality_score?.toFixed(1) || 'N/A',
      icon: Layers,
      href: '/qc',
    },
    {
      title: 'Projects',
      value: projects?.length || 0,
      icon: Microscope,
      href: '/taxonomy',
    },
  ]

  const isLoading = projectsLoading || qcLoading

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <div>
          <h1 className="text-3xl font-bold">Dashboard Overview</h1>
          <p className="text-muted-foreground">
            Welcome to MMonitor. View and analyze your metagenomic data.
          </p>
        </div>
        <div className="w-64">
          <ProjectFilter />
        </div>
      </div>

      {isLoading ? (
        <div className="flex items-center justify-center py-12">
          <Loader2 className="h-8 w-8 animate-spin text-muted-foreground" />
        </div>
      ) : (
        <>
          <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
            {stats.map((stat) => (
              <Link key={stat.title} to={stat.href}>
                <Card className="transition-colors hover:bg-muted/50">
                  <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
                    <CardTitle className="text-sm font-medium">{stat.title}</CardTitle>
                    <stat.icon className="h-4 w-4 text-muted-foreground" />
                  </CardHeader>
                  <CardContent>
                    <div className="text-2xl font-bold">{stat.value}</div>
                  </CardContent>
                </Card>
              </Link>
            ))}
          </div>

          <div className="grid gap-4 md:grid-cols-2">
            <Card>
              <CardHeader>
                <CardTitle>Quick Actions</CardTitle>
                <CardDescription>Navigate to different analysis modules</CardDescription>
              </CardHeader>
              <CardContent className="space-y-2">
                <Link
                  to="/taxonomy"
                  className="flex items-center gap-3 rounded-lg border p-3 transition-colors hover:bg-muted"
                >
                  <Dna className="h-5 w-5 text-primary" />
                  <div>
                    <div className="font-medium">Taxonomy Analysis</div>
                    <div className="text-sm text-muted-foreground">
                      View taxonomic composition charts
                    </div>
                  </div>
                </Link>
                <Link
                  to="/qc"
                  className="flex items-center gap-3 rounded-lg border p-3 transition-colors hover:bg-muted"
                >
                  <BarChart3 className="h-5 w-5 text-primary" />
                  <div>
                    <div className="font-medium">Quality Control</div>
                    <div className="text-sm text-muted-foreground">
                      Review sequencing quality metrics
                    </div>
                  </div>
                </Link>
                <Link
                  to="/diversity"
                  className="flex items-center gap-3 rounded-lg border p-3 transition-colors hover:bg-muted"
                >
                  <Layers className="h-5 w-5 text-primary" />
                  <div>
                    <div className="font-medium">Diversity Analysis</div>
                    <div className="text-sm text-muted-foreground">
                      Explore alpha and beta diversity
                    </div>
                  </div>
                </Link>
                <Link
                  to="/mags"
                  className="flex items-center gap-3 rounded-lg border p-3 transition-colors hover:bg-muted"
                >
                  <Microscope className="h-5 w-5 text-primary" />
                  <div>
                    <div className="font-medium">MAG Analysis</div>
                    <div className="text-sm text-muted-foreground">
                      Explore metagenome-assembled genomes
                    </div>
                  </div>
                </Link>
              </CardContent>
            </Card>

            <Card>
              <CardHeader>
                <CardTitle>Recent Projects</CardTitle>
                <CardDescription>Your available projects</CardDescription>
              </CardHeader>
              <CardContent>
                {projects && projects.length > 0 ? (
                  <div className="space-y-2">
                    {projects.slice(0, 5).map((project) => (
                      <div
                        key={project.id}
                        className="flex items-center justify-between rounded-lg border p-3"
                      >
                        <div>
                          <div className="font-medium">{project.name}</div>
                          <div className="text-sm text-muted-foreground">
                            {project.sample_count} samples
                          </div>
                        </div>
                      </div>
                    ))}
                  </div>
                ) : (
                  <p className="text-sm text-muted-foreground">
                    No projects found. Upload data to get started.
                  </p>
                )}
              </CardContent>
            </Card>
          </div>
        </>
      )}
    </div>
  )
}
