import { useQuery } from '@tanstack/react-query'
import { projectsApi } from '@/api/endpoints'
import { useFilterStore } from '@/stores/filterStore'
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'

export function ProjectFilter() {
  const { projectId, setProjectId } = useFilterStore()

  const { data: projects, isLoading } = useQuery({
    queryKey: ['projects'],
    queryFn: projectsApi.list,
  })

  return (
    <Select
      value={projectId || ''}
      onValueChange={(value) => setProjectId(value || null)}
    >
      <SelectTrigger>
        <SelectValue placeholder={isLoading ? 'Loading...' : 'Select project'} />
      </SelectTrigger>
      <SelectContent>
        {projects?.map((project) => (
          <SelectItem key={project.id} value={project.id}>
            {project.name} ({project.sample_count} samples)
          </SelectItem>
        ))}
      </SelectContent>
    </Select>
  )
}
