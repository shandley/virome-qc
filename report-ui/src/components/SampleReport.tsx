import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { fmt, pct } from "@/lib/utils";
import type { Passport, Histogram, ErvAnalysis, ErvLocus } from "@/types";
import {
  AreaChart,
  Area,
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
} from "recharts";

// -- Helpers --

function fmtAxis(n: number): string {
  if (n >= 1_000_000) return (n / 1_000_000).toFixed(1) + "M";
  if (n >= 1_000) return (n / 1_000).toFixed(0) + "K";
  return n.toString();
}

function fmtBases(n: number): string {
  if (n >= 1e9) return (n / 1e9).toFixed(1) + " Gb";
  if (n >= 1e6) return (n / 1e6).toFixed(1) + " Mb";
  if (n >= 1e3) return (n / 1e3).toFixed(1) + " Kb";
  return n + " bp";
}

/** Compute weighted mean from a histogram */
function histMean(h: Histogram): number | null {
  if (!h || !h.counts || h.total === 0) return null;
  let sum = 0;
  for (let i = 0; i < h.counts.length; i++) {
    const binCenter =
      i + 1 < h.bin_edges.length
        ? (h.bin_edges[i] + h.bin_edges[i + 1]) / 2
        : h.bin_edges[i];
    sum += binCenter * h.counts[i];
  }
  return sum / h.total;
}

/** Compute weighted median from a histogram */
function histMedian(h: Histogram): number | null {
  if (!h || !h.counts || h.total === 0) return null;
  const half = h.total / 2;
  let cumulative = 0;
  for (let i = 0; i < h.counts.length; i++) {
    cumulative += h.counts[i];
    if (cumulative >= half) {
      return h.bin_edges[i];
    }
  }
  return h.bin_edges[h.bin_edges.length - 1];
}

const tooltipStyle = {
  backgroundColor: "var(--card)",
  border: "1px solid var(--border)",
  borderRadius: "var(--radius)",
  fontSize: 12,
};

// -- Stat Cards --

function StatCard({
  value,
  label,
  sub,
  color,
}: {
  value: string;
  label: string;
  sub?: string;
  color?: string;
}) {
  return (
    <Card>
      <CardContent className="p-4 text-center">
        <div
          className="text-2xl font-bold font-mono"
          style={color ? { color } : undefined}
        >
          {value}
        </div>
        <div className="text-[11px] font-semibold text-muted-foreground uppercase tracking-wider mt-1">
          {label}
        </div>
        {sub && (
          <div className="text-[10px] text-muted-foreground mt-1">{sub}</div>
        )}
      </CardContent>
    </Card>
  );
}

// -- Charts --

function HistogramChart({
  data,
  color,
  xLabel,
  yLabel = "Reads",
  formatX,
}: {
  data: Histogram;
  color: string;
  xLabel?: string;
  yLabel?: string;
  formatX?: (v: number) => string;
}) {
  const chartData = data.counts.map((count, i) => ({
    bin: formatX
      ? formatX(data.bin_edges[i])
      : data.bin_edges[i] < 2
        ? data.bin_edges[i].toFixed(2)
        : Math.round(data.bin_edges[i]).toString(),
    count,
  }));

  return (
    <ResponsiveContainer width="100%" height={220}>
      <BarChart data={chartData} margin={{ top: 5, right: 10, left: 10, bottom: 5 }}>
        <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
        <XAxis
          dataKey="bin"
          tick={{ fontSize: 10 }}
          className="fill-muted-foreground"
          interval="preserveStartEnd"
          label={xLabel ? { value: xLabel, position: "bottom", offset: -2, fontSize: 11 } : undefined}
        />
        <YAxis
          tick={{ fontSize: 10 }}
          className="fill-muted-foreground"
          tickFormatter={fmtAxis}
          label={yLabel ? { value: yLabel, angle: -90, position: "insideLeft", offset: 0, fontSize: 11 } : undefined}
        />
        <Tooltip
          contentStyle={tooltipStyle}
          formatter={(value: number) => [value.toLocaleString(), yLabel]}
        />
        <Bar dataKey="count" fill={color} radius={[2, 2, 0, 0]} opacity={0.85} />
      </BarChart>
    </ResponsiveContainer>
  );
}

function OverlayHistogramChart({
  before,
  after,
  xLabel,
  yLabel = "Reads",
  formatX,
}: {
  before: Histogram;
  after: Histogram;
  xLabel?: string;
  yLabel?: string;
  formatX?: (v: number) => string;
}) {
  const chartData = before.counts.map((count, i) => ({
    bin: formatX
      ? formatX(before.bin_edges[i])
      : before.bin_edges[i] < 2
        ? before.bin_edges[i].toFixed(2)
        : Math.round(before.bin_edges[i]).toString(),
    before: count,
    after: after.counts[i] || 0,
  }));

  return (
    <ResponsiveContainer width="100%" height={220}>
      <BarChart data={chartData} margin={{ top: 5, right: 10, left: 10, bottom: 5 }}>
        <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
        <XAxis
          dataKey="bin"
          tick={{ fontSize: 10 }}
          className="fill-muted-foreground"
          interval="preserveStartEnd"
          label={xLabel ? { value: xLabel, position: "bottom", offset: -2, fontSize: 11 } : undefined}
        />
        <YAxis
          tick={{ fontSize: 10 }}
          className="fill-muted-foreground"
          tickFormatter={fmtAxis}
          label={yLabel ? { value: yLabel, angle: -90, position: "insideLeft", offset: 0, fontSize: 11 } : undefined}
        />
        <Tooltip
          contentStyle={tooltipStyle}
          formatter={(value: number, name: string) => [value.toLocaleString(), name === "before" ? "Before QC" : "After QC"]}
        />
        <Legend iconSize={10} wrapperStyle={{ fontSize: 11 }} verticalAlign="top"
          formatter={(value: string) => value === "before" ? "Before QC" : "After QC"}
        />
        <Bar dataKey="before" fill="var(--chart-1)" opacity={0.5} radius={[2, 2, 0, 0]} />
        <Bar dataKey="after" fill="var(--chart-3)" opacity={0.85} radius={[2, 2, 0, 0]} />
      </BarChart>
    </ResponsiveContainer>
  );
}

function QualityProfileChart({
  positions,
  title,
}: {
  positions: { position: number; mean: number; median: number; q25: number; q75: number }[];
  title: string;
}) {
  if (!positions?.length) return null;

  const step = positions.length > 300 ? Math.ceil(positions.length / 300) : 1;
  const data = positions
    .filter((_, i) => i % step === 0)
    .map((p) => ({
      ...p,
      mean: Math.min(p.mean, 42),
      median: Math.min(p.median, 42),
      q25: Math.min(p.q25, 42),
      q75: Math.min(p.q75, 42),
      iqrBand: [Math.min(p.q25, 42), Math.min(p.q75, 42)] as [number, number],
    }));

  return (
    <Card>
      <CardHeader className="pb-2">
        <CardTitle className="text-sm">{title}</CardTitle>
      </CardHeader>
      <CardContent>
        <ResponsiveContainer width="100%" height={260}>
          <AreaChart data={data} margin={{ top: 5, right: 10, left: 10, bottom: 20 }}>
            <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
            <XAxis
              dataKey="position"
              tick={{ fontSize: 10 }}
              className="fill-muted-foreground"
              label={{ value: "Position", position: "bottom", offset: 2, fontSize: 11 }}
            />
            <YAxis
              domain={[0, 42]}
              allowDataOverflow={true}
              tick={{ fontSize: 10 }}
              className="fill-muted-foreground"
              label={{ value: "Q Score", angle: -90, position: "insideLeft", offset: 0, fontSize: 11 }}
            />
            <Tooltip
              contentStyle={tooltipStyle}
              formatter={(value: number | number[], name: string) => {
                if (Array.isArray(value)) return [`Q${value[0].toFixed(1)} - Q${value[1].toFixed(1)}`, "IQR"];
                return [`Q${value.toFixed(1)}`, name];
              }}
            />
            <Legend iconSize={10} wrapperStyle={{ fontSize: 11 }} verticalAlign="top" />
            <Area
              type="monotone"
              dataKey="iqrBand"
              fill="var(--chart-1)"
              fillOpacity={0.25}
              stroke="var(--chart-1)"
              strokeWidth={0.5}
              strokeOpacity={0.3}
              name="IQR (Q25-Q75)"
              legendType="rect"
            />
            <Area
              type="monotone"
              dataKey="median"
              fill="none"
              stroke="var(--chart-1)"
              strokeWidth={2}
              name="Median"
            />
            <Area
              type="monotone"
              dataKey="mean"
              fill="none"
              stroke="var(--chart-1)"
              strokeWidth={1}
              strokeDasharray="4 3"
              strokeOpacity={0.5}
              name="Mean"
            />
          </AreaChart>
        </ResponsiveContainer>
      </CardContent>
    </Card>
  );
}

function BaseCompositionChart({
  positions,
  title,
}: {
  positions: { position: number; a: number; c: number; g: number; t: number; n: number }[];
  title: string;
}) {
  if (!positions?.length) return null;

  const step = positions.length > 300 ? Math.ceil(positions.length / 300) : 1;
  const data = positions
    .filter((_, i) => i % step === 0)
    .map((p) => {
      const total = p.a + p.c + p.g + p.t + p.n;
      const norm = total > 0 ? 100 / total : 0;
      return {
        position: p.position,
        A: +(p.a * norm).toFixed(1),
        C: +(p.c * norm).toFixed(1),
        G: +(p.g * norm).toFixed(1),
        T: +(p.t * norm).toFixed(1),
        N: +(p.n * norm).toFixed(1),
      };
    });

  return (
    <Card>
      <CardHeader className="pb-2">
        <CardTitle className="text-sm">{title}</CardTitle>
      </CardHeader>
      <CardContent>
        <ResponsiveContainer width="100%" height={200}>
          <AreaChart data={data} margin={{ top: 5, right: 10, left: 10, bottom: 5 }}>
            <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
            <XAxis dataKey="position" tick={{ fontSize: 10 }} className="fill-muted-foreground" />
            <YAxis domain={[0, 100]} allowDataOverflow={true} tick={{ fontSize: 10 }} className="fill-muted-foreground" />
            <Tooltip
              contentStyle={tooltipStyle}
              formatter={(value: number) => [`${value}%`]}
            />
            <Legend iconSize={10} wrapperStyle={{ fontSize: 11 }} />
            <Area type="monotone" dataKey="A" stackId="1" fill="var(--chart-3)" stroke="none" />
            <Area type="monotone" dataKey="C" stackId="1" fill="var(--chart-2)" stroke="none" />
            <Area type="monotone" dataKey="G" stackId="1" fill="var(--chart-4)" stroke="none" />
            <Area type="monotone" dataKey="T" stackId="1" fill="var(--chart-5)" stroke="none" />
            <Area type="monotone" dataKey="N" stackId="1" fill="var(--muted-foreground)" fillOpacity={0.3} stroke="none" />
          </AreaChart>
        </ResponsiveContainer>
      </CardContent>
    </Card>
  );
}

// -- ERV Analysis --

function classificationColor(cls: string): string {
  switch (cls) {
    case "Endogenous": return "var(--chart-3)";
    case "Exogenous": return "var(--destructive)";
    default: return "var(--chart-4)";
  }
}

function ErvAnalysisCard({ erv, readsInput }: { erv: ErvAnalysis; readsInput: number }) {
  const { classifications: cls, loci } = erv;
  const totalClassified = cls.endogenous + cls.ambiguous + cls.exogenous;
  const ervFrac = erv.retroviral_reads_flagged / Math.max(readsInput, 1);

  // Sort loci: exogenous first, then ambiguous, then endogenous, by reads desc
  const sortedLoci = [...loci].sort((a, b) => {
    const order = { Exogenous: 0, Ambiguous: 1, Endogenous: 2 };
    const oa = order[a.classification] ?? 1;
    const ob = order[b.classification] ?? 1;
    return oa !== ob ? oa - ob : b.reads - a.reads;
  });

  // Only show top loci (up to 20)
  const displayLoci = sortedLoci.slice(0, 20);
  const hasMore = sortedLoci.length > 20;

  // Prepare data for classification donut-style summary
  const classData = [
    { label: "Endogenous", count: cls.endogenous, color: "var(--chart-3)" },
    { label: "Ambiguous", count: cls.ambiguous, color: "var(--chart-4)" },
    { label: "Exogenous", count: cls.exogenous, color: "var(--destructive)" },
  ];

  // CpG distribution for bar chart
  const cpgBins: number[] = new Array(10).fill(0);
  for (const l of loci) {
    const bin = Math.min(Math.floor(l.cpg_ratio * 10), 9);
    cpgBins[bin] += l.reads;
  }
  const cpgData = cpgBins.map((count, i) => ({
    bin: (i * 0.1).toFixed(1),
    count,
  }));

  return (
    <Card>
      <CardHeader className="pb-2">
        <CardTitle className="text-sm flex items-center gap-2">
          Retroviral Analysis
          {cls.exogenous > 0 && (
            <Badge variant="fail" className="text-[10px]">
              {cls.exogenous} exogenous
            </Badge>
          )}
        </CardTitle>
      </CardHeader>
      <CardContent className="space-y-4">
        {/* Summary stats */}
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
          <div className="text-center p-3 rounded-md bg-muted/50">
            <div className="text-lg font-bold font-mono">{fmt(erv.retroviral_reads_flagged)}</div>
            <div className="text-[10px] font-semibold text-muted-foreground uppercase">Retroviral Reads</div>
            <div className="text-[10px] text-muted-foreground">{pct(ervFrac)} of input</div>
          </div>
          <div className="text-center p-3 rounded-md bg-muted/50">
            <div className="text-lg font-bold font-mono">{erv.clusters_total}</div>
            <div className="text-[10px] font-semibold text-muted-foreground uppercase">Clusters</div>
            <div className="text-[10px] text-muted-foreground">{totalClassified} classified</div>
          </div>
          {classData.map(({ label, count, color }) => (
            <div key={label} className="text-center p-3 rounded-md bg-muted/50">
              <div className="text-lg font-bold font-mono" style={{ color }}>{count}</div>
              <div className="text-[10px] font-semibold text-muted-foreground uppercase">{label}</div>
              <div className="text-[10px] text-muted-foreground">
                {totalClassified > 0 ? pct(count / totalClassified) : "0%"} of clusters
              </div>
            </div>
          ))}
        </div>

        {/* Classification bar */}
        {totalClassified > 0 && (
          <div className="space-y-1">
            <div className="text-[11px] font-semibold text-muted-foreground uppercase tracking-wider">Classification Distribution</div>
            <div className="flex h-5 rounded-full overflow-hidden bg-muted">
              {classData.map(({ label, count, color }) => {
                const w = (count / totalClassified) * 100;
                if (w === 0) return null;
                return (
                  <div
                    key={label}
                    className="h-full flex items-center justify-center text-[9px] font-semibold text-white"
                    style={{ width: `${w}%`, backgroundColor: color, minWidth: w > 5 ? undefined : 4 }}
                    title={`${label}: ${count} (${w.toFixed(1)}%)`}
                  >
                    {w > 15 ? `${label} ${count}` : w > 8 ? `${count}` : ""}
                  </div>
                );
              })}
            </div>
          </div>
        )}

        {/* CpG distribution */}
        {loci.length > 0 && (
          <div>
            <div className="text-[11px] font-semibold text-muted-foreground uppercase tracking-wider mb-2">
              CpG Observed/Expected (reads by ratio)
            </div>
            <ResponsiveContainer width="100%" height={140}>
              <BarChart data={cpgData} margin={{ top: 5, right: 10, left: 10, bottom: 5 }}>
                <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
                <XAxis dataKey="bin" tick={{ fontSize: 10 }} className="fill-muted-foreground" />
                <YAxis tick={{ fontSize: 10 }} className="fill-muted-foreground" tickFormatter={fmtAxis} />
                <Tooltip contentStyle={tooltipStyle} formatter={(value: number) => [value.toLocaleString(), "Reads"]} />
                <Bar
                  dataKey="count"
                  radius={[2, 2, 0, 0]}
                  opacity={0.85}
                  fill="var(--chart-2)"
                />
              </BarChart>
            </ResponsiveContainer>
            <div className="text-[10px] text-muted-foreground mt-1">
              CpG &lt; 0.4 = endogenous (methylation-depleted). CpG &gt; 0.8 = exogenous (intact).
            </div>
          </div>
        )}

        {/* Loci table */}
        {displayLoci.length > 0 && (
          <div className="overflow-x-auto">
            <div className="text-[11px] font-semibold text-muted-foreground uppercase tracking-wider mb-2">
              Top Clusters
            </div>
            <table className="w-full text-xs">
              <thead>
                <tr className="border-b border-border">
                  <th className="text-left py-1.5 px-2 font-semibold text-muted-foreground">Match</th>
                  <th className="text-center py-1.5 px-2 font-semibold text-muted-foreground">Class</th>
                  <th className="text-right py-1.5 px-2 font-semibold text-muted-foreground">Reads</th>
                  <th className="text-right py-1.5 px-2 font-semibold text-muted-foreground">CpG O/E</th>
                  <th className="text-right py-1.5 px-2 font-semibold text-muted-foreground">Score</th>
                  <th className="text-right py-1.5 px-2 font-semibold text-muted-foreground">ORFs</th>
                </tr>
              </thead>
              <tbody>
                {displayLoci.map((l) => (
                  <tr key={l.cluster_id} className="border-b border-border/50 hover:bg-muted/30">
                    <td className="py-1.5 px-2 font-medium truncate max-w-[160px]" title={l.best_match}>
                      {l.best_match === "unknown" ? "—" : l.best_match}
                    </td>
                    <td className="py-1.5 px-2 text-center">
                      <span
                        className="inline-block px-1.5 py-0.5 rounded text-[10px] font-semibold"
                        style={{
                          color: classificationColor(l.classification),
                          backgroundColor: `color-mix(in srgb, ${classificationColor(l.classification)} 15%, transparent)`,
                        }}
                      >
                        {l.classification}
                      </span>
                    </td>
                    <td className="py-1.5 px-2 text-right font-mono">{l.reads}</td>
                    <td className="py-1.5 px-2 text-right font-mono">{l.cpg_ratio.toFixed(2)}</td>
                    <td className="py-1.5 px-2 text-right font-mono">{l.combined_score.toFixed(2)}</td>
                    <td className="py-1.5 px-2 text-right font-mono">{l.orf_intact}</td>
                  </tr>
                ))}
              </tbody>
            </table>
            {hasMore && (
              <div className="text-[10px] text-muted-foreground mt-1">
                Showing top 20 of {sortedLoci.length} clusters
              </div>
            )}
          </div>
        )}
      </CardContent>
    </Card>
  );
}

// -- Main Report --

export function SampleReport({ passport }: { passport: Passport }) {
  const p = passport;
  const ri = Math.max(p.reads_input, 1);
  const qa = p.qa_stats;

  const getModuleExtra = (name: string, key: string) => {
    const m = p.modules.find((m) => m.name === name);
    return (m?.extra?.[key] as number) || 0;
  };
  const getModuleRemoved = (name: string) =>
    p.modules.find((m) => m.name === name)?.reads_removed || 0;

  // Derived metrics from distributions
  const meanReadLen = qa?.distributions?.length_after
    ? histMean(qa.distributions.length_after)
    : null;
  const medianReadLen = qa?.distributions?.length_after
    ? histMedian(qa.distributions.length_after)
    : null;
  const meanQuality = qa?.distributions?.quality_scores
    ? histMean(qa.distributions.quality_scores)
    : null;
  const meanGC = qa?.distributions?.gc_content
    ? histMean(qa.distributions.gc_content)
    : null;

  const basesIn = qa?.summary?.bases_input || 0;
  const basesOut = qa?.summary?.bases_output || 0;
  const libraryComplexity = qa?.duplication?.estimated_library_complexity || 0;
  const hasPaired = (p.pairs_passed ?? 0) > 0;

  return (
    <div className="space-y-6">
      {/* Top-level summary: reads in/out, survival, bases */}
      <div className="grid grid-cols-2 sm:grid-cols-3 lg:grid-cols-6 gap-3">
        <StatCard
          value={fmt(p.reads_input)}
          label="Reads In"
          sub={basesIn > 0 ? fmtBases(basesIn) : `${p.reads_input.toLocaleString()} reads`}
        />
        <StatCard
          value={fmt(p.reads_passed)}
          label="Reads Out"
          sub={basesOut > 0 ? fmtBases(basesOut) : `${p.reads_passed.toLocaleString()} reads`}
        />
        <StatCard
          value={pct(p.survival_rate)}
          label="Survival"
          sub={`${fmt(p.reads_input - p.reads_passed)} removed`}
          color={p.quality_tier === "FAIL" ? "var(--destructive)" : p.quality_tier === "WARN" ? "var(--chart-4)" : "var(--chart-3)"}
        />
        {meanReadLen !== null && (
          <StatCard
            value={`${Math.round(meanReadLen)} bp`}
            label="Mean Length"
            sub={medianReadLen !== null ? `median ${Math.round(medianReadLen)} bp` : undefined}
          />
        )}
        {meanQuality !== null && (
          <StatCard
            value={`Q${meanQuality.toFixed(1)}`}
            label="Mean Quality"
            sub={meanQuality >= 30 ? "high quality" : meanQuality >= 20 ? "acceptable" : "low quality"}
            color={meanQuality >= 30 ? "var(--chart-3)" : meanQuality >= 20 ? "var(--chart-4)" : "var(--destructive)"}
          />
        )}
        {meanGC !== null && (
          <StatCard
            value={`${(meanGC * 100).toFixed(1)}%`}
            label="Mean GC"
          />
        )}
      </div>

      {/* Library metrics: contamination summary, complexity, duplication, paired-end */}
      <div className="grid grid-cols-2 sm:grid-cols-3 lg:grid-cols-4 gap-3">
        {p.contamination_summary && p.contamination_summary.biological_contamination_removed > 0 && (
          <StatCard
            value={pct(p.contamination_summary.biological_contamination_fraction)}
            label="Biological Contam."
            sub={`${fmt(p.contamination_summary.biological_contamination_removed)} reads (host + rRNA + contaminant)`}
            color={p.contamination_summary.biological_contamination_fraction > 0.10 ? "var(--destructive)" : p.contamination_summary.biological_contamination_fraction > 0.02 ? "var(--chart-4)" : "var(--chart-3)"}
          />
        )}
        {qa && (
          <StatCard
            value={pct(qa.duplication.estimated_duplication_rate)}
            label="Duplication"
            sub={`${fmt(qa.duplication.estimated_unique_sequences)} unique`}
          />
        )}
        {libraryComplexity > 0 && (
          <StatCard
            value={fmt(libraryComplexity)}
            label="Library Complexity"
            sub="estimated distinct molecules"
          />
        )}
        {hasPaired && (
          <StatCard
            value={fmt(p.pairs_passed!)}
            label="Pairs Output"
            sub={`${fmt(p.singletons || 0)} singletons (${pct((p.singletons || 0) / ri)})`}
          />
        )}
        {hasPaired && (p.pairs_merged || 0) > 0 && (
          <StatCard
            value={pct((p.pairs_merged || 0) / Math.max(p.pairs_passed!, 1))}
            label="Merged"
            sub={`${fmt(p.pairs_merged || 0)} overlapping pairs`}
          />
        )}
      </div>

      {/* Flags -- only shown when there are flags */}
      {p.flags.length > 0 && (
        <div className="space-y-2">
          {p.flags.map((f, i) => (
            <div
              key={i}
              className={`rounded-lg border p-3 text-sm ${
                f.severity === "FAIL"
                  ? "border-destructive bg-destructive/10"
                  : "border-chart-4 bg-chart-4/10"
              }`}
            >
              <span className="font-semibold">{f.code}</span> {f.message}
            </div>
          ))}
        </div>
      )}

      {/* Survival funnel */}
      <Card>
        <CardHeader className="pb-2">
          <CardTitle className="text-sm">Survival Funnel</CardTitle>
        </CardHeader>
        <CardContent>
          <div className="overflow-x-auto">
            <table className="w-full text-xs">
              <thead>
                <tr className="border-b border-border">
                  <th className="text-left py-2 px-3 font-semibold text-muted-foreground">Module</th>
                  <th className="text-right py-2 px-3 font-semibold text-muted-foreground">Processed</th>
                  <th className="text-right py-2 px-3 font-semibold text-muted-foreground">Removed</th>
                  <th className="text-right py-2 px-3 font-semibold text-muted-foreground">% of Input</th>
                  <th className="text-right py-2 px-3 font-semibold text-muted-foreground">Remaining</th>
                </tr>
              </thead>
              <tbody>
                {(() => {
                  let remaining = p.reads_input;
                  return p.modules.map((m, i) => {
                    const mpct = (m.reads_removed / ri) * 100;
                    remaining -= m.reads_removed;
                    return (
                      <tr key={i} className="border-b border-border/50 hover:bg-muted/30">
                        <td className="py-2 px-3 font-medium">{m.name}</td>
                        <td className="py-2 px-3 text-right font-mono">{fmt(m.reads_processed)}</td>
                        <td className="py-2 px-3 text-right font-mono">{fmt(m.reads_removed)}</td>
                        <td className="py-2 px-3 text-right font-mono">
                          <div className="flex items-center justify-end gap-2">
                            <div className="w-16 h-1.5 rounded-full bg-muted overflow-hidden">
                              <div
                                className="h-full rounded-full bg-primary"
                                style={{ width: `${Math.min(mpct * 2, 100)}%` }}
                              />
                            </div>
                            {mpct.toFixed(2)}%
                          </div>
                        </td>
                        <td className="py-2 px-3 text-right font-mono">{fmt(remaining)}</td>
                      </tr>
                    );
                  });
                })()}
              </tbody>
            </table>
          </div>
        </CardContent>
      </Card>

      {/* Adapters + Contaminants side by side */}
      <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
        {qa && qa.adapters.breakdown.length > 0 && (
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-sm">Adapters</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-2 gap-2">
                {qa.adapters.breakdown.map((a, i) => (
                  <div key={i} className="text-center p-3 rounded-md bg-muted/50">
                    <div className="text-lg font-bold font-mono text-primary">{fmt(a.count)}</div>
                    <div className="text-[10px] font-semibold text-muted-foreground uppercase">{a.name}</div>
                    <div className="text-[10px] text-muted-foreground">{pct(a.count / ri)}</div>
                  </div>
                ))}
              </div>
            </CardContent>
          </Card>
        )}

        {getModuleRemoved("contaminant") > 0 && (
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-sm">Contaminants</CardTitle>
            </CardHeader>
            <CardContent>
              <div className="grid grid-cols-3 gap-2">
                {[
                  { key: "rrna_removed", label: "rRNA", color: "text-chart-1",
                    subKeys: [
                      { key: "rrna_prokaryotic", label: "prokaryotic" },
                      { key: "rrna_eukaryotic", label: "eukaryotic" },
                    ],
                  },
                  { key: "phix_removed", label: "PhiX", color: "text-chart-2" },
                  { key: "vector_removed", label: "Vector", color: "text-chart-4" },
                ].map(({ key, label, color, subKeys }) => {
                  const val = getModuleExtra("contaminant", key);
                  if (!val) return null;
                  let subText = `${fmt(val)} reads`;
                  if (subKeys) {
                    const parts = subKeys
                      .map(({ key: sk, label: sl }) => {
                        const sv = getModuleExtra("contaminant", sk);
                        return sv > 0 ? `${fmt(sv)} ${sl}` : null;
                      })
                      .filter(Boolean);
                    if (parts.length > 0) subText = parts.join(", ");
                  }
                  return (
                    <div key={key} className="text-center p-3 rounded-md bg-muted/50">
                      <div className={`text-lg font-bold font-mono ${color}`}>{pct(val / ri)}</div>
                      <div className="text-[10px] font-semibold text-muted-foreground uppercase">{label}</div>
                      <div className="text-[10px] text-muted-foreground">{subText}</div>
                    </div>
                  );
                })}
              </div>
            </CardContent>
          </Card>
        )}
      </div>

      {/* Host + rRNA cards */}
      {(getModuleRemoved("host") > 0 || getModuleRemoved("rrna") > 0) && (
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
          {getModuleRemoved("host") > 0 && (() => {
            const hostRemoved = getModuleRemoved("host");
            const hostFrac = hostRemoved / ri;
            const ambiguous = getModuleExtra("host", "ambiguous_flagged");
            const hostModule = p.modules.find((m) => m.name === "host");
            const containmentDist = hostModule?.extra?.containment_distribution as
              | { bins: string[]; counts: number[]; zero_containment: number }
              | undefined;
            return (
              <Card>
                <CardHeader className="pb-2">
                  <CardTitle className="text-sm">Host Depletion</CardTitle>
                </CardHeader>
                <CardContent className="space-y-3">
                  <div className="grid grid-cols-2 gap-2">
                    <div className="text-center p-3 rounded-md bg-muted/50">
                      <div className={`text-lg font-bold font-mono ${hostFrac > 0.20 ? "text-destructive" : hostFrac > 0.01 ? "text-chart-4" : "text-chart-3"}`}>
                        {pct(hostFrac)}
                      </div>
                      <div className="text-[10px] font-semibold text-muted-foreground uppercase">Host Reads</div>
                      <div className="text-[10px] text-muted-foreground">{fmt(hostRemoved)} removed</div>
                    </div>
                    <div className="text-center p-3 rounded-md bg-muted/50">
                      <div className="text-lg font-bold font-mono">{fmt(ambiguous)}</div>
                      <div className="text-[10px] font-semibold text-muted-foreground uppercase">Ambiguous</div>
                      <div className="text-[10px] text-muted-foreground">{pct(ambiguous / ri)} flagged</div>
                    </div>
                  </div>
                  {containmentDist && containmentDist.counts.some((c: number) => c > 0) && (
                    <div>
                      <div className="text-[11px] font-semibold text-muted-foreground uppercase tracking-wider mb-2">
                        Containment Distribution (reads with &gt;0% host k-mers)
                      </div>
                      <ResponsiveContainer width="100%" height={120}>
                        <BarChart
                          data={containmentDist.bins.map((bin: string, i: number) => ({
                            bin,
                            count: containmentDist.counts[i],
                          }))}
                          margin={{ top: 5, right: 10, left: 10, bottom: 5 }}
                        >
                          <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
                          <XAxis dataKey="bin" tick={{ fontSize: 9 }} className="fill-muted-foreground" />
                          <YAxis tick={{ fontSize: 9 }} className="fill-muted-foreground" tickFormatter={fmtAxis} />
                          <Tooltip contentStyle={tooltipStyle} formatter={(value: number) => [value.toLocaleString(), "Reads"]} />
                          <Bar dataKey="count" fill="var(--chart-1)" radius={[2, 2, 0, 0]} opacity={0.85} />
                        </BarChart>
                      </ResponsiveContainer>
                      <div className="text-[10px] text-muted-foreground mt-1">
                        Green: kept. Yellow: ambiguous (20-50%). Red: removed as host (&gt;50%).
                        {containmentDist.zero_containment > 0 && ` ${fmt(containmentDist.zero_containment)} reads had 0% containment (not shown).`}
                      </div>
                    </div>
                  )}
                </CardContent>
              </Card>
            );
          })()}
          {getModuleRemoved("rrna") > 0 && (() => {
            const rrnaRemoved = getModuleRemoved("rrna");
            const rrnaFrac = rrnaRemoved / ri;
            return (
              <Card>
                <CardHeader className="pb-2">
                  <CardTitle className="text-sm">rRNA Screening</CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="text-center p-3 rounded-md bg-muted/50">
                    <div className={`text-lg font-bold font-mono ${rrnaFrac > 0.10 ? "text-destructive" : rrnaFrac > 0.01 ? "text-chart-4" : "text-chart-3"}`}>
                      {pct(rrnaFrac)}
                    </div>
                    <div className="text-[10px] font-semibold text-muted-foreground uppercase">rRNA Reads</div>
                    <div className="text-[10px] text-muted-foreground">{fmt(rrnaRemoved)} removed (SILVA)</div>
                  </div>
                </CardContent>
              </Card>
            );
          })()}
        </div>
      )}

      {/* ERV Analysis */}
      {p.erv_analysis && p.erv_analysis.retroviral_reads_flagged > 0 && (
        <ErvAnalysisCard erv={p.erv_analysis} readsInput={p.reads_input} />
      )}

      {/* Quality profiles */}
      {qa?.per_position && (
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
          <QualityProfileChart
            positions={qa.per_position.quality_before.positions}
            title="Quality Profile (Before QC)"
          />
          <QualityProfileChart
            positions={qa.per_position.quality_after.positions}
            title="Quality Profile (After QC)"
          />
        </div>
      )}

      {/* Base composition */}
      {qa?.per_position && (
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
          <BaseCompositionChart
            positions={qa.per_position.bases_before.positions}
            title="Base Composition (Before QC)"
          />
          <BaseCompositionChart
            positions={qa.per_position.bases_after.positions}
            title="Base Composition (After QC)"
          />
        </div>
      )}

      {/* Distributions */}
      {qa?.distributions && (
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-sm">Read Length Distribution</CardTitle>
            </CardHeader>
            <CardContent>
              {qa.distributions.length_after.total > 0 ? (
                <OverlayHistogramChart
                  before={qa.distributions.length_before}
                  after={qa.distributions.length_after}
                />
              ) : (
                <HistogramChart data={qa.distributions.length_before} color="var(--chart-1)" />
              )}
            </CardContent>
          </Card>
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-sm">GC Content</CardTitle>
            </CardHeader>
            <CardContent>
              <HistogramChart
                data={qa.distributions.gc_content}
                color="var(--chart-3)"
                formatX={(v) => `${(v * 100).toFixed(0)}%`}
              />
            </CardContent>
          </Card>
          <Card>
            <CardHeader className="pb-2">
              <CardTitle className="text-sm">Quality Scores</CardTitle>
            </CardHeader>
            <CardContent>
              <HistogramChart data={qa.distributions.quality_scores} color="var(--chart-2)" />
            </CardContent>
          </Card>
          {qa.distributions.insert_sizes.total > 0 && (
            <Card>
              <CardHeader className="pb-2">
                <CardTitle className="text-sm">Insert Sizes</CardTitle>
              </CardHeader>
              <CardContent>
                <HistogramChart data={qa.distributions.insert_sizes} color="var(--chart-4)" />
              </CardContent>
            </Card>
          )}
        </div>
      )}
    </div>
  );
}
