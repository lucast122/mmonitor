import { useState } from 'react'
import { useQuery } from '@tanstack/react-query'
import { Search } from 'lucide-react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Input } from '@/components/ui/input'
import { Button } from '@/components/ui/button'
import { magsApi } from '@/api/endpoints'
import type { MAGAnnotation } from '@/types'

interface GeneSearchTableProps {
  magId: number
}

export function GeneSearchTable({ magId }: GeneSearchTableProps) {
  const [searchTerm, setSearchTerm] = useState('')
  const [cogFilter, setCogFilter] = useState('')
  const [pfamFilter, setPfamFilter] = useState('')
  const [keggFilter, setKeggFilter] = useState('')
  const [page, setPage] = useState(1)
  const pageSize = 20

  const { data, isLoading, refetch } = useQuery({
    queryKey: ['magGeneSearch', magId, searchTerm, cogFilter, pfamFilter, keggFilter, page],
    queryFn: () =>
      magsApi.geneSearch(magId, {
        search: searchTerm || undefined,
        cog: cogFilter || undefined,
        pfam: pfamFilter || undefined,
        kegg: keggFilter || undefined,
        page,
        page_size: pageSize,
      }),
    enabled: !!magId,
  })

  const handleSearch = () => {
    setPage(1)
    refetch()
  }

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      handleSearch()
    }
  }

  return (
    <Card>
      <CardHeader>
        <CardTitle>Gene Search</CardTitle>
      </CardHeader>
      <CardContent>
        <div className="mb-4 grid gap-2 md:grid-cols-5">
          <Input
            placeholder="Search product/gene ID..."
            value={searchTerm}
            onChange={(e) => setSearchTerm(e.target.value)}
            onKeyDown={handleKeyDown}
          />
          <Input
            placeholder="COG filter..."
            value={cogFilter}
            onChange={(e) => setCogFilter(e.target.value)}
            onKeyDown={handleKeyDown}
          />
          <Input
            placeholder="PFAM filter..."
            value={pfamFilter}
            onChange={(e) => setPfamFilter(e.target.value)}
            onKeyDown={handleKeyDown}
          />
          <Input
            placeholder="KEGG filter..."
            value={keggFilter}
            onChange={(e) => setKeggFilter(e.target.value)}
            onKeyDown={handleKeyDown}
          />
          <Button onClick={handleSearch}>
            <Search className="mr-2 h-4 w-4" />
            Search
          </Button>
        </div>

        {isLoading ? (
          <div className="flex h-64 items-center justify-center">
            <div className="animate-pulse text-muted-foreground">Loading...</div>
          </div>
        ) : (
          <>
            <div className="overflow-x-auto rounded-md border">
              <table className="w-full text-sm">
                <thead className="bg-muted/50">
                  <tr>
                    <th className="px-4 py-3 text-left font-medium">Gene ID</th>
                    <th className="px-4 py-3 text-left font-medium">Contig</th>
                    <th className="px-4 py-3 text-left font-medium">Position</th>
                    <th className="px-4 py-3 text-left font-medium">Strand</th>
                    <th className="px-4 py-3 text-left font-medium">Product</th>
                    <th className="px-4 py-3 text-left font-medium">COG</th>
                    <th className="px-4 py-3 text-left font-medium">PFAM</th>
                    <th className="px-4 py-3 text-left font-medium">KEGG</th>
                  </tr>
                </thead>
                <tbody className="divide-y">
                  {data?.results && data.results.length > 0 ? (
                    data.results.map((annotation: MAGAnnotation) => (
                      <tr key={annotation.id} className="hover:bg-muted/50">
                        <td className="px-4 py-2 font-mono text-sm">{annotation.gene_id}</td>
                        <td className="px-4 py-2 font-mono text-sm">{annotation.contig}</td>
                        <td className="px-4 py-2 text-sm">
                          {annotation.start.toLocaleString()}-{annotation.stop.toLocaleString()}
                        </td>
                        <td className="px-4 py-2">{annotation.strand}</td>
                        <td className="max-w-xs truncate px-4 py-2" title={annotation.product}>
                          {annotation.product}
                        </td>
                        <td className="px-4 py-2">{annotation.cog || '-'}</td>
                        <td className="px-4 py-2">{annotation.pfam || '-'}</td>
                        <td className="px-4 py-2">{annotation.kegg || '-'}</td>
                      </tr>
                    ))
                  ) : (
                    <tr>
                      <td colSpan={8} className="px-4 py-8 text-center text-muted-foreground">
                        No genes found
                      </td>
                    </tr>
                  )}
                </tbody>
              </table>
            </div>

            {data && data.count > pageSize && (
              <div className="mt-4 flex items-center justify-between">
                <div className="text-sm text-muted-foreground">
                  Showing {(page - 1) * pageSize + 1} to{' '}
                  {Math.min(page * pageSize, data.count)} of {data.count} genes
                </div>
                <div className="flex gap-2">
                  <Button
                    variant="outline"
                    size="sm"
                    disabled={page === 1}
                    onClick={() => setPage(page - 1)}
                  >
                    Previous
                  </Button>
                  <Button
                    variant="outline"
                    size="sm"
                    disabled={page * pageSize >= data.count}
                    onClick={() => setPage(page + 1)}
                  >
                    Next
                  </Button>
                </div>
              </div>
            )}
          </>
        )}
      </CardContent>
    </Card>
  )
}
