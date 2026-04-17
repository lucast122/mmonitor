import { useQuery } from '@tanstack/react-query'
import { samplesApi } from '@/api/endpoints'
import { useFilterStore } from '@/stores/filterStore'
import { ProjectFilter } from './ProjectFilter'
import { SampleFilter } from './SampleFilter'
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select'

interface FilterBarProps {
  showDate?: boolean
  showSubproject?: boolean
  selectedDate?: string
  selectedSubproject?: string
  onDateChange?: (date: string) => void
  onSubprojectChange?: (subproject: string) => void
}

export function FilterBar({
  showDate = true,
  showSubproject = true,
  selectedDate = '__all__',
  selectedSubproject = '__all__',
  onDateChange,
  onSubprojectChange,
}: FilterBarProps) {
  const { projectId } = useFilterStore()

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

  const dates = datesData?.dates || []
  const subprojects = subprojectsData?.subprojects || []

  return (
    <div className="flex flex-wrap gap-4">
      <div className="w-48">
        <label className="mb-1 block text-sm font-medium">Project</label>
        <ProjectFilter />
      </div>
      <div className="w-64">
        <label className="mb-1 block text-sm font-medium">Samples</label>
        <SampleFilter />
      </div>
      {showDate && onDateChange && (
        <div className="w-48">
          <label className="mb-1 block text-sm font-medium">Date</label>
          <Select value={selectedDate} onValueChange={onDateChange}>
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
      )}
      {showSubproject && onSubprojectChange && (
        <div className="w-48">
          <label className="mb-1 block text-sm font-medium">Subproject</label>
          <Select value={selectedSubproject} onValueChange={onSubprojectChange}>
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
      )}
    </div>
  )
}
