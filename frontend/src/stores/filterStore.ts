import { create } from 'zustand'
import type { FilterState } from '@/types'

interface FilterStore extends FilterState {
  setProjectId: (projectId: string | null) => void
  setSampleIds: (sampleIds: string[]) => void
  toggleSampleId: (sampleId: string) => void
  setDateRange: (start: string | null, end: string | null) => void
  clearFilters: () => void
}

export const useFilterStore = create<FilterStore>((set) => ({
  projectId: null,
  sampleIds: [],
  dateRange: {
    start: null,
    end: null,
  },

  setProjectId: (projectId) =>
    set({
      projectId,
      sampleIds: [], // Clear samples when project changes
    }),

  setSampleIds: (sampleIds) => set({ sampleIds }),

  toggleSampleId: (sampleId) =>
    set((state) => ({
      sampleIds: state.sampleIds.includes(sampleId)
        ? state.sampleIds.filter((id) => id !== sampleId)
        : [...state.sampleIds, sampleId],
    })),

  setDateRange: (start, end) =>
    set({
      dateRange: { start, end },
    }),

  clearFilters: () =>
    set({
      projectId: null,
      sampleIds: [],
      dateRange: { start: null, end: null },
    }),
}))
