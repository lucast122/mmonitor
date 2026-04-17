import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import type { MAGFunctionalSummary } from '@/types'

interface FunctionalSummaryProps {
  summary: MAGFunctionalSummary | null
  isLoading?: boolean
}

export function FunctionalSummary({ summary, isLoading }: FunctionalSummaryProps) {
  if (isLoading) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>Functional Annotation Summary</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="flex h-32 items-center justify-center">
            <div className="animate-pulse text-muted-foreground">Loading...</div>
          </div>
        </CardContent>
      </Card>
    )
  }

  if (!summary) {
    return (
      <Card>
        <CardHeader>
          <CardTitle>Functional Annotation Summary</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="text-muted-foreground">No functional annotation data available</div>
        </CardContent>
      </Card>
    )
  }

  const stats = [
    {
      label: 'Genes with COG',
      value: summary.genes_with_cog,
      color: 'bg-blue-500',
    },
    {
      label: 'Genes with KEGG',
      value: summary.genes_with_kegg,
      color: 'bg-green-500',
    },
    {
      label: 'Genes with PFAM',
      value: summary.genes_with_pfam,
      color: 'bg-purple-500',
    },
    {
      label: 'Genes with GO',
      value: summary.genes_with_go,
      color: 'bg-orange-500',
    },
  ]

  const cogCategories = Object.keys(summary.cog_distribution || {}).length
  const keggKOs = Object.keys(summary.kegg_kos || {}).length
  const pfamDomains = Object.keys(summary.pfam_domains || {}).length

  return (
    <Card>
      <CardHeader>
        <CardTitle>Functional Annotation Summary</CardTitle>
      </CardHeader>
      <CardContent>
        <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
          {stats.map((stat) => (
            <div key={stat.label} className="rounded-lg border p-4">
              <div className="flex items-center gap-2">
                <div className={`h-3 w-3 rounded-full ${stat.color}`} />
                <span className="text-sm text-muted-foreground">{stat.label}</span>
              </div>
              <div className="mt-2 text-2xl font-bold">{stat.value.toLocaleString()}</div>
            </div>
          ))}
        </div>

        <div className="mt-6 grid gap-4 md:grid-cols-3">
          <div className="rounded-lg border p-4">
            <div className="text-sm text-muted-foreground">COG Categories</div>
            <div className="mt-1 text-xl font-semibold">{cogCategories}</div>
          </div>
          <div className="rounded-lg border p-4">
            <div className="text-sm text-muted-foreground">KEGG KOs</div>
            <div className="mt-1 text-xl font-semibold">{keggKOs}</div>
          </div>
          <div className="rounded-lg border p-4">
            <div className="text-sm text-muted-foreground">PFAM Domains</div>
            <div className="mt-1 text-xl font-semibold">{pfamDomains}</div>
          </div>
        </div>
      </CardContent>
    </Card>
  )
}
