// Types matching the virome-qc passport JSON structure

export interface Passport {
  tool_version: string;
  profile: string;
  reads_input: number;
  reads_passed: number;
  survival_rate: number;
  unique_reads?: number;
  qc_survival_rate?: number;
  pairs_passed?: number;
  singletons?: number;
  pairs_merged?: number;
  modules: ModuleReport[];
  flags: QcFlag[];
  quality_tier: "PASS" | "WARN" | "FAIL";
  qa_stats?: AnalyticsSnapshot;
  provenance?: Provenance;
  ingestion?: Record<string, unknown>;
  parameters?: Record<string, unknown>;
  contamination_summary?: {
    biological_contamination_removed: number;
    biological_contamination_fraction: number;
    technical_artifacts_removed: number;
    technical_artifacts_fraction: number;
    total_removed: number;
    total_removed_fraction: number;
  };
}

/** Expected QC metric ranges from ViroForge calibration. */
export interface ExpectedRanges {
  survival?: [number, number];
  host_fraction?: [number, number];
  rrna_fraction?: [number, number];
  adapter_rate?: [number, number];
  duplication_rate?: [number, number];
  gc_content?: [number, number];
}

/** Extract expected ranges from passport parameters. */
export function getExpectedRanges(p: Passport): ExpectedRanges | null {
  const er = (p.parameters as Record<string, unknown> | undefined)
    ?.quality_assessment as Record<string, unknown> | undefined;
  const ranges = er?.expected_ranges as ExpectedRanges | undefined;
  if (!ranges) return null;
  // Fallback: if gc_content not in expected_ranges, try expected_gc_range
  if (!ranges.gc_content && er?.expected_gc_range) {
    const gc = er.expected_gc_range as [number, number];
    if (Array.isArray(gc) && gc.length === 2) {
      ranges.gc_content = gc;
    }
  }
  return ranges;
}

export interface ModuleReport {
  name: string;
  reads_processed: number;
  reads_removed: number;
  reads_modified: number;
  bases_removed: number;
  extra: Record<string, unknown>;
}

export interface QcFlag {
  code: string;
  message: string;
  severity: "PASS" | "WARN" | "FAIL";
}

export interface Provenance {
  timestamp: string;
  input_files: { path: string; size_bytes: number }[];
}

export interface AnalyticsSnapshot {
  summary: {
    reads_input: number;
    reads_passed: number;
    reads_failed: number;
    bases_input: number;
    bases_output: number;
    survival_rate: number;
  };
  per_position: {
    quality_before: PerPositionQuality;
    quality_after: PerPositionQuality;
    bases_before: PerPositionBases;
    bases_after: PerPositionBases;
  };
  distributions: {
    length_before: Histogram;
    length_after: Histogram;
    gc_content: Histogram;
    quality_scores: Histogram;
    insert_sizes: Histogram;
    trimmed_bases: Histogram;
  };
  adapters: { breakdown: { name: string; count: number }[] };
  duplication: {
    estimated_unique_sequences: number;
    estimated_duplication_rate: number;
    estimated_library_complexity: number;
  };
}

export interface PerPositionQuality {
  positions: {
    position: number;
    count: number;
    mean: number;
    q25: number;
    median: number;
    q75: number;
  }[];
}

export interface PerPositionBases {
  positions: {
    position: number;
    a: number;
    c: number;
    g: number;
    t: number;
    n: number;
    total: number;
  }[];
}

export interface Histogram {
  bin_edges: number[];
  counts: number[];
  underflow: number;
  overflow: number;
  total: number;
}
