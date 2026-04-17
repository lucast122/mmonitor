import { useQuery } from '@tanstack/react-query'
import { samplesApi } from '@/api/endpoints'
import { useFilterStore } from '@/stores/filterStore'
import { Button } from '@/components/ui/button'
import {
  DropdownMenu,
  DropdownMenuCheckboxItem,
  DropdownMenuContent,
  DropdownMenuLabel,
  DropdownMenuSeparator,
  DropdownMenuTrigger,
} from '@/components/ui/dropdown-menu'
import { ChevronDown } from 'lucide-react'

export function SampleFilter() {
  const { projectId, sampleIds, toggleSampleId, setSampleIds } = useFilterStore()

  const { data: samples, isLoading } = useQuery({
    queryKey: ['samples', projectId],
    queryFn: () => samplesApi.list(projectId || undefined),
    enabled: projectId !== null,
  })

  const selectedCount = sampleIds.length

  return (
    <DropdownMenu>
      <DropdownMenuTrigger asChild>
        <Button variant="outline" className="w-full justify-between">
          {isLoading
            ? 'Loading...'
            : selectedCount > 0
              ? `${selectedCount} sample${selectedCount > 1 ? 's' : ''} selected`
              : 'Select samples'}
          <ChevronDown className="ml-2 h-4 w-4" />
        </Button>
      </DropdownMenuTrigger>
      <DropdownMenuContent className="w-64" align="start">
        <DropdownMenuLabel>Samples</DropdownMenuLabel>
        <DropdownMenuSeparator />
        {samples && samples.length > 0 && (
          <>
            <DropdownMenuCheckboxItem
              checked={sampleIds.length === samples.length}
              onCheckedChange={(checked) => {
                if (checked) {
                  setSampleIds(samples.map((s) => s.id))
                } else {
                  setSampleIds([])
                }
              }}
            >
              Select All
            </DropdownMenuCheckboxItem>
            <DropdownMenuSeparator />
          </>
        )}
        <div className="max-h-64 overflow-y-auto">
          {samples?.map((sample) => (
            <DropdownMenuCheckboxItem
              key={sample.id}
              checked={sampleIds.includes(sample.id)}
              onCheckedChange={() => toggleSampleId(sample.id)}
            >
              {sample.name}
            </DropdownMenuCheckboxItem>
          ))}
        </div>
        {(!samples || samples.length === 0) && !isLoading && (
          <div className="p-2 text-center text-sm text-muted-foreground">
            {projectId ? 'No samples found' : 'Select a project first'}
          </div>
        )}
      </DropdownMenuContent>
    </DropdownMenu>
  )
}
