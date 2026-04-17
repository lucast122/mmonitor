// User types
export interface User {
  id: number
  username: string
  email: string
  first_name: string
  last_name: string
  is_staff: boolean
  date_joined: string
}

export interface AuthTokens {
  access: string
  refresh: string
}

export interface LoginResponse {
  user: User
  tokens: AuthTokens
}

// Project types
export interface Project {
  id: string  // project_id can be string like 'R1', 'R2'
  name: string
  sample_count: number
  created_at?: string
}

// Sample types
export interface Sample {
  id: string  // sample_id can be string
  name: string
  project_id: string
  date?: string
  subproject?: string
}

// Taxonomy types
export interface TaxonomyRecord {
  id: number
  taxonomy: string
  abundance: number
  sample_id: string
  project_id: string
  subproject?: string
  date?: string
  tax_genus?: string
  tax_family?: string
  tax_order?: string
  tax_class?: string
  tax_phylum?: string
  tax_superkingdom?: string
}

export interface TaxonomyAggregated {
  sample_id: string
  sample_name: string
  taxa: string[]
  abundances: number[]
  relative_abundances: number[]
}

export interface TaxonomicRank {
  value: string
  label: string
}

// QC types
export interface QCMetrics {
  id: number
  sample_name: string
  project_id: string
  subproject_id?: string
  date?: string
  mean_read_length?: number
  median_read_length?: number
  mean_quality_score?: number
  mean_gc_content?: number
  number_of_reads?: number
  total_bases_sequenced?: number
  q20_score?: number
  q30_score?: number
}

export interface QCDistribution {
  sample_name: string
  read_lengths: number[]
  qualities: number[]
  gc_contents: number[]
}

export interface QCSummary {
  total_samples: number
  total_reads: number
  total_bases: number
  avg_quality_score?: number
  avg_read_length?: number
  avg_gc_content?: number
}

// Diversity types
export interface AlphaDiversity {
  sample_id: string
  sample_name: string
  shannon: number
  simpson: number
  observed_species: number
  chao1?: number
  pielou_evenness?: number
}

export interface BetaDiversity {
  sample_ids: string[]
  sample_names: string[]
  distance_matrix: number[][]
  metric: string
}

export interface PCoAResult {
  sample_ids: string[]
  sample_names: string[]
  coordinates: number[][]
  explained_variance: number[]
  metric: string
}

export interface RarefactionCurve {
  sample_id: string
  sample_name: string
  depths: number[]
  observed_species: number[]
  standard_deviation?: number[]
}

// MAG types
export interface MAG {
  id: number
  name: string
  taxonomy: Record<string, string> | null
  sample_name: string
  completeness: number | null
  contamination: number | null
  strain_heterogeneity: number | null
  genome_size: number | null
  n_contigs: number | null
  n50: number | null
  gc_content: number | null
  gene_count: number | null
  created_at: string
}

export interface MAGDetail extends MAG {
  trna_count: number | null
  rrna_count: number | null
  cds_count: number | null
  functional_summary: MAGFunctionalSummary | null
  annotation_summary: string | null
  quality_data: Record<string, unknown> | null
  annotations_count: number
  contigs_count: number
}

export interface MAGFunctionalSummary {
  cog_distribution: Record<string, number>
  kegg_kos: Record<string, number>
  pfam_domains: Record<string, number>
  genes_with_cog: number
  genes_with_kegg: number
  genes_with_pfam: number
  genes_with_go: number
}

export interface MAGAnnotation {
  id: number
  gene_id: string
  contig: string
  start: number
  stop: number
  strand: string
  product: string
  kegg?: string
  ec_number?: string
  ec_numbers: string[]
  cog?: string
  cog_category?: string
  pfam?: string
  pfam_description?: string
  go_terms?: string[]
  uniref: string[]
  uploaded_at: string
}

export interface MAGContig {
  id?: number
  contig_id: string
  length: number
  gc_content?: number | null
  gene_count?: number
  genes?: {
    gene_id: string
    start: number
    stop: number
    strand: string
    product: string
  }[]
}

export interface COGDistributionItem {
  category: string
  description: string
  count: number
}

export interface PFAMDomainItem {
  pfam: string
  description: string
  count: number
}

// Pagination types
export interface PaginatedResponse<T> {
  count: number
  next: string | null
  previous: string | null
  total_pages: number
  current_page: number
  results: T[]
}

// Filter types
export interface FilterState {
  projectId: string | null
  sampleIds: string[]
  dateRange: {
    start: string | null
    end: string | null
  }
}
