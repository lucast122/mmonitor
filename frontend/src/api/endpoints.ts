import apiClient from './client'
import type {
  LoginResponse,
  User,
  Project,
  Sample,
  TaxonomyAggregated,
  TaxonomicRank,
  QCMetrics,
  QCDistribution,
  QCSummary,
  AlphaDiversity,
  BetaDiversity,
  PCoAResult,
  RarefactionCurve,
  MAG,
  MAGDetail,
  MAGAnnotation,
  MAGContig,
  PaginatedResponse,
} from '@/types'

// Auth endpoints
export const authApi = {
  login: async (username: string, password: string): Promise<LoginResponse> => {
    const response = await apiClient.post('/auth/login/', { username, password })
    return response.data
  },

  register: async (data: {
    username: string
    email: string
    password: string
    password_confirm: string
  }): Promise<LoginResponse> => {
    const response = await apiClient.post('/auth/register/', data)
    return response.data
  },

  refresh: async (refreshToken: string): Promise<{ access: string }> => {
    const response = await apiClient.post('/auth/refresh/', { refresh: refreshToken })
    return response.data
  },

  getProfile: async (): Promise<User> => {
    const response = await apiClient.get('/auth/me/')
    return response.data
  },

  updateProfile: async (data: Partial<User>): Promise<User> => {
    const response = await apiClient.patch('/auth/me/', data)
    return response.data
  },
}

// Projects endpoints
export const projectsApi = {
  list: async (): Promise<Project[]> => {
    const response = await apiClient.get('/projects/')
    return response.data
  },

  get: async (id: string): Promise<Project> => {
    const response = await apiClient.get(`/projects/${id}/`)
    return response.data
  },
}

// Samples endpoints
export const samplesApi = {
  list: async (projectId?: string): Promise<Sample[]> => {
    const params = projectId ? { project_id: projectId } : {}
    const response = await apiClient.get('/samples/', { params })
    return response.data
  },

  get: async (id: string): Promise<Sample> => {
    const response = await apiClient.get(`/samples/${id}/`)
    return response.data
  },

  getDates: async (projectId?: string): Promise<{ dates: string[] }> => {
    const params = projectId ? { project_id: projectId } : {}
    const response = await apiClient.get('/samples/dates/', { params })
    return response.data
  },

  getSubprojects: async (projectId?: string): Promise<{ subprojects: string[] }> => {
    const params = projectId ? { project_id: projectId } : {}
    const response = await apiClient.get('/samples/subprojects/', { params })
    return response.data
  },

  listByFilters: async (params: {
    project_id?: string
    date?: string
    subproject?: string
  }): Promise<Sample[]> => {
    const response = await apiClient.get('/samples/by_filters/', { params })
    return response.data
  },
}

// Taxonomy endpoints
export const taxonomyApi = {
  list: async (params: {
    project_id?: string
    sample_id?: string[]
    page?: number
    page_size?: number
  }): Promise<PaginatedResponse<TaxonomyAggregated>> => {
    const response = await apiClient.get('/taxonomy/', { params })
    return response.data
  },

  getAggregated: async (params: {
    project_id?: string
    sample_id?: string[]
    date?: string
    subproject?: string
    rank?: string
    top_n?: number
    group_by?: string
  }): Promise<TaxonomyAggregated[] | TaxonomyAggregated> => {
    const response = await apiClient.get('/taxonomy/aggregated/', { params })
    return response.data
  },

  getRanks: async (): Promise<{ ranks: TaxonomicRank[] }> => {
    const response = await apiClient.get('/taxonomy/ranks/')
    return response.data
  },

  getCountTableUrl: (params: {
    project_id?: string
    sample_id?: string[]
    date?: string
    subproject?: string
    rank?: string
  }): string => {
    const searchParams = new URLSearchParams()
    if (params.project_id) searchParams.append('project_id', params.project_id)
    if (params.sample_id) params.sample_id.forEach(id => searchParams.append('sample_id', id))
    if (params.date) searchParams.append('date', params.date)
    if (params.subproject) searchParams.append('subproject', params.subproject)
    if (params.rank) searchParams.append('rank', params.rank)
    return `/api/v1/taxonomy/count_table/?${searchParams.toString()}`
  },

  getNormalizedTableUrl: (params: {
    project_id?: string
    sample_id?: string[]
    date?: string
    subproject?: string
    rank?: string
  }): string => {
    const searchParams = new URLSearchParams()
    if (params.project_id) searchParams.append('project_id', params.project_id)
    if (params.sample_id) params.sample_id.forEach(id => searchParams.append('sample_id', id))
    if (params.date) searchParams.append('date', params.date)
    if (params.subproject) searchParams.append('subproject', params.subproject)
    if (params.rank) searchParams.append('rank', params.rank)
    return `/api/v1/taxonomy/normalized_table/?${searchParams.toString()}`
  },
}

// QC endpoints
export const qcApi = {
  list: async (params: {
    project_id?: string
    sample_name?: string[]
    date?: string
    subproject?: string
    page?: number
    page_size?: number
  }): Promise<PaginatedResponse<QCMetrics> | QCMetrics[]> => {
    const response = await apiClient.get('/qc/', { params })
    return response.data
  },

  get: async (id: number): Promise<QCMetrics> => {
    const response = await apiClient.get(`/qc/${id}/`)
    return response.data
  },

  getDistributions: async (params: {
    project_id?: string
    sample_name?: string[]
    date?: string
    subproject?: string
  }): Promise<QCDistribution[]> => {
    const response = await apiClient.get('/qc/distributions/', { params })
    return response.data
  },

  getSummary: async (params: {
    project_id?: string
    date?: string
    subproject?: string
  }): Promise<QCSummary> => {
    const response = await apiClient.get('/qc/summary/', { params })
    return response.data
  },

  getViolinData: async (params: {
    project_id?: string
    date?: string
    subproject?: string
  }): Promise<{
    quality_scores: { sample: string; value: number }[]
    read_lengths: { sample: string; value: number }[]
    gc_contents: { sample: string; value: number }[]
  }> => {
    const response = await apiClient.get('/qc/violin_data/', { params })
    return response.data
  },

  getMetricsList: async (params: {
    project_id?: string
    date?: string
    subproject?: string
  }): Promise<{
    samples: string[]
    mean_quality_scores: (number | null)[]
    mean_read_lengths: (number | null)[]
    mean_gc_contents: (number | null)[]
    q20_scores: (number | null)[]
    q30_scores: (number | null)[]
    number_of_reads: (number | null)[]
  }> => {
    const response = await apiClient.get('/qc/metrics_list/', { params })
    return response.data
  },
}

// Diversity endpoints
export const diversityApi = {
  getAlpha: async (params: {
    project_id?: string
    sample_id?: string[]
    date?: string
    subproject?: string
  }): Promise<AlphaDiversity[]> => {
    const response = await apiClient.get('/diversity/alpha/', { params })
    return response.data
  },

  getBeta: async (params: {
    project_id?: string
    sample_id?: string[]
    date?: string
    subproject?: string
    metric?: string
  }): Promise<BetaDiversity> => {
    const response = await apiClient.get('/diversity/beta/', { params })
    return response.data
  },

  getPCoA: async (params: {
    project_id?: string
    sample_id?: string[]
    date?: string
    subproject?: string
    metric?: string
    n_components?: number
  }): Promise<PCoAResult> => {
    const response = await apiClient.get('/diversity/pcoa/', { params })
    return response.data
  },

  getRarefaction: async (params: {
    project_id?: string
    sample_id?: string[]
    date?: string
    subproject?: string
    n_steps?: number
    n_iterations?: number
  }): Promise<RarefactionCurve[]> => {
    const response = await apiClient.get('/diversity/rarefaction/', { params })
    return response.data
  },

  getExportAlphaUrl: (params: {
    project_id?: string
    sample_id?: string[]
    date?: string
    subproject?: string
  }): string => {
    const searchParams = new URLSearchParams()
    if (params.project_id) searchParams.append('project_id', params.project_id)
    if (params.sample_id) params.sample_id.forEach(id => searchParams.append('sample_id', id))
    if (params.date) searchParams.append('date', params.date)
    if (params.subproject) searchParams.append('subproject', params.subproject)
    return `/api/v1/diversity/export_alpha/?${searchParams.toString()}`
  },

  getExportBetaUrl: (params: {
    project_id?: string
    sample_id?: string[]
    date?: string
    subproject?: string
    metric?: string
  }): string => {
    const searchParams = new URLSearchParams()
    if (params.project_id) searchParams.append('project_id', params.project_id)
    if (params.sample_id) params.sample_id.forEach(id => searchParams.append('sample_id', id))
    if (params.date) searchParams.append('date', params.date)
    if (params.subproject) searchParams.append('subproject', params.subproject)
    if (params.metric) searchParams.append('metric', params.metric)
    return `/api/v1/diversity/export_beta/?${searchParams.toString()}`
  },

  getExportPCoAUrl: (params: {
    project_id?: string
    sample_id?: string[]
    date?: string
    subproject?: string
    metric?: string
    n_components?: number
  }): string => {
    const searchParams = new URLSearchParams()
    if (params.project_id) searchParams.append('project_id', params.project_id)
    if (params.sample_id) params.sample_id.forEach(id => searchParams.append('sample_id', id))
    if (params.date) searchParams.append('date', params.date)
    if (params.subproject) searchParams.append('subproject', params.subproject)
    if (params.metric) searchParams.append('metric', params.metric)
    if (params.n_components) searchParams.append('n_components', params.n_components.toString())
    return `/api/v1/diversity/export_pcoa/?${searchParams.toString()}`
  },
}

// Horizon endpoints (used within Taxonomy page)
export const horizonApi = {
  getData: async (params: {
    project_id?: string
    subproject?: string
    taxa_count?: number
  }): Promise<{
    dates: string[]
    taxa: string[]
    data: { date: string; taxon: string; value: number; percentage: number; mean: number }[]
    mean_abundances: Record<string, number>
  }> => {
    const response = await apiClient.get('/horizon/data/', { params })
    return response.data
  },
}

// MAGs endpoints
export const magsApi = {
  list: async (params: {
    sample_name?: string
    min_completeness?: number
    max_contamination?: number
    page?: number
    page_size?: number
  }): Promise<PaginatedResponse<MAG>> => {
    const response = await apiClient.get('/mags/', { params })
    return response.data
  },

  get: async (id: number): Promise<MAGDetail> => {
    const response = await apiClient.get(`/mags/${id}/`)
    return response.data
  },

  getContigs: async (id: number): Promise<MAGContig[]> => {
    const response = await apiClient.get(`/mags/${id}/contigs/`)
    return response.data
  },

  getAnnotations: async (
    id: number,
    params: {
      contig?: string
      search?: string
      page?: number
      page_size?: number
    }
  ): Promise<PaginatedResponse<MAGAnnotation>> => {
    const response = await apiClient.get(`/mags/${id}/annotations/`, { params })
    return response.data
  },

  getQualitySummary: async (sampleName?: string): Promise<
    { id: number; name: string; completeness: number; contamination: number; sample_name: string }[]
  > => {
    const params = sampleName ? { sample_name: sampleName } : {}
    const response = await apiClient.get('/mags/quality_summary/', { params })
    return response.data
  },

  getFunctionalSummary: async (id: number): Promise<{
    cog_distribution: Record<string, number>
    kegg_kos: Record<string, number>
    pfam_domains: Record<string, number>
    go_terms: Record<string, number>
    total_annotated: number
    genes_with_cog: number
    genes_with_kegg: number
    genes_with_pfam: number
  }> => {
    const response = await apiClient.get(`/mags/${id}/functional_summary/`)
    return response.data
  },

  getCogDistribution: async (id: number): Promise<
    { category: string; description: string; count: number }[]
  > => {
    const response = await apiClient.get(`/mags/${id}/cog_distribution/`)
    return response.data
  },

  getPfamDomains: async (id: number, topN?: number): Promise<
    { pfam: string; description: string; count: number }[]
  > => {
    const params = topN ? { top_n: topN } : {}
    const response = await apiClient.get(`/mags/${id}/pfam_domains/`, { params })
    return response.data
  },

  getContigsSummary: async (id: number): Promise<
    { id: number; contig_id: string; length: number; gc_content: number | null; gene_count: number }[]
  > => {
    const response = await apiClient.get(`/mags/${id}/contigs_summary/`)
    return response.data
  },

  geneSearch: async (
    id: number,
    params: {
      search?: string
      cog?: string
      pfam?: string
      kegg?: string
      contig?: string
      page?: number
      page_size?: number
    }
  ): Promise<PaginatedResponse<MAGAnnotation>> => {
    const response = await apiClient.get(`/mags/${id}/gene_search/`, { params })
    return response.data
  },

  getFastaUrl: (id: number): string => {
    return `/api/v1/mags/${id}/fasta/`
  },

  getGffUrl: (id: number): string => {
    return `/api/v1/mags/${id}/gff/`
  },

  getFastaIndexUrl: (id: number): string => {
    return `/api/v1/mags/${id}/fasta_index/`
  },
}

// Correlations endpoints
export const correlationsApi = {
  getMetadataFields: async (projectId?: string): Promise<{ fields: string[] }> => {
    const params = projectId ? { project_id: projectId } : {}
    const response = await apiClient.get('/correlations/metadata_fields/', { params })
    return response.data
  },

  getMetadata: async (projectId?: string): Promise<Array<Record<string, unknown>>> => {
    const params = projectId ? { project_id: projectId } : {}
    const response = await apiClient.get('/correlations/metadata/', { params })
    return response.data
  },

  uploadMetadata: async (projectId: string | undefined, file: File): Promise<{
    message: string
    created: number
    updated: number
  }> => {
    const formData = new FormData()
    formData.append('file', file)
    if (projectId) formData.append('project_id', projectId)
    const response = await apiClient.post('/correlations/upload_metadata/', formData, {
      headers: { 'Content-Type': 'multipart/form-data' },
    })
    return response.data
  },

  saveMetadata: async (projectId: string | undefined, metadata: Array<Record<string, unknown>>): Promise<{
    message: string
    created: number
    updated: number
  }> => {
    const response = await apiClient.post('/correlations/save_metadata/', {
      project_id: projectId,
      metadata,
    })
    return response.data
  },

  compute: async (params: {
    project_id?: string
    method?: 'pearson' | 'spearman' | 'kendall'
    p_cutoff?: number
    top_n_taxa?: number
  }): Promise<{
    taxa: string[]
    metadata_fields: string[]
    correlations: (number | null)[][]
    p_values: (number | null)[][]
    method: string
    p_cutoff: number
    n_samples: number
  }> => {
    const response = await apiClient.get('/correlations/compute/', { params })
    return response.data
  },

  getExportUrl: (params: {
    project_id?: string
    method?: string
    top_n_taxa?: number
  }): string => {
    const searchParams = new URLSearchParams()
    if (params.project_id) searchParams.append('project_id', params.project_id)
    if (params.method) searchParams.append('method', params.method)
    if (params.top_n_taxa) searchParams.append('top_n_taxa', params.top_n_taxa.toString())
    return `/api/v1/correlations/export/?${searchParams.toString()}`
  },
}
