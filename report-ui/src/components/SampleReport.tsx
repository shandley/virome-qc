import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Badge } from "@/components/ui/badge";
import { fmt, pct } from "@/lib/utils";
import type { Passport, Histogram } from "@/types";
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
          label={yLabel ? { value: yLabel, angle: -90, position: "insideLeft", offset: 0, fontSize: 11 } : undefined}
        />
        <Tooltip
          contentStyle={{
            backgroundColor: "var(--card)",
            border: "1px solid var(--border)",
            borderRadius: "var(--radius)",
            fontSize: 12,
          }}
          formatter={(value: number) => [value.toLocaleString(), yLabel]}
        />
        <Bar dataKey="count" fill={color} radius={[2, 2, 0, 0]} opacity={0.85} />
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

  // Subsample if too many positions for smooth rendering
  const step = positions.length > 300 ? Math.ceil(positions.length / 300) : 1;
  const data = positions.filter((_, i) => i % step === 0);

  return (
    <Card>
      <CardHeader className="pb-2">
        <CardTitle className="text-sm">{title}</CardTitle>
      </CardHeader>
      <CardContent>
        <ResponsiveContainer width="100%" height={220}>
          <AreaChart data={data} margin={{ top: 5, right: 10, left: 10, bottom: 5 }}>
            <CartesianGrid strokeDasharray="3 3" className="stroke-border" />
            <XAxis
              dataKey="position"
              tick={{ fontSize: 10 }}
              className="fill-muted-foreground"
              label={{ value: "Position", position: "bottom", offset: -2, fontSize: 11 }}
            />
            <YAxis
              domain={[0, 42]}
              tick={{ fontSize: 10 }}
              className="fill-muted-foreground"
              label={{ value: "Q Score", angle: -90, position: "insideLeft", offset: 0, fontSize: 11 }}
            />
            <Tooltip
              contentStyle={{
                backgroundColor: "var(--card)",
                border: "1px solid var(--border)",
                borderRadius: "var(--radius)",
                fontSize: 12,
              }}
              formatter={(value: number, name: string) => [
                `Q${value.toFixed(1)}`,
                name,
              ]}
            />
            <Legend iconSize={10} wrapperStyle={{ fontSize: 11 }} />
            <Area
              type="monotone"
              dataKey="q75"
              stackId="band"
              fill="var(--chart-1)"
              fillOpacity={0.1}
              stroke="none"
              name="Q75"
            />
            <Area
              type="monotone"
              dataKey="q25"
              stackId="band"
              fill="var(--card)"
              fillOpacity={1}
              stroke="none"
              name="Q25"
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
    .map((p) => ({
      position: p.position,
      A: +(p.a * 100).toFixed(1),
      C: +(p.c * 100).toFixed(1),
      G: +(p.g * 100).toFixed(1),
      T: +(p.t * 100).toFixed(1),
      N: +(p.n * 100).toFixed(1),
    }));

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
            <YAxis domain={[0, 100]} tick={{ fontSize: 10 }} className="fill-muted-foreground" />
            <Tooltip
              contentStyle={{
                backgroundColor: "var(--card)",
                border: "1px solid var(--border)",
                borderRadius: "var(--radius)",
                fontSize: 12,
              }}
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

export function SampleReport({ passport }: { passport: Passport }) {
  const p = passport;
  const ri = Math.max(p.reads_input, 1);
  const qa = p.qa_stats;

  // Extract module metrics
  const getModuleExtra = (name: string, key: string) => {
    const m = p.modules.find((m) => m.name === name);
    return (m?.extra?.[key] as number) || 0;
  };
  const getModuleRemoved = (name: string) =>
    p.modules.find((m) => m.name === name)?.reads_removed || 0;

  const dedupRemoved = getModuleRemoved("dedup");
  const dedupRate = dedupRemoved / ri;

  return (
    <div className="space-y-6">
      {/* Summary cards */}
      <div className="grid grid-cols-2 sm:grid-cols-3 lg:grid-cols-5 gap-3">
        <StatCard value={fmt(p.reads_input)} label="Input" sub={`${p.reads_input.toLocaleString()} reads`} />
        <StatCard value={fmt(p.reads_passed)} label="Passed" sub={`${p.reads_passed.toLocaleString()} reads`} />
        <StatCard
          value={pct(p.survival_rate)}
          label="Survival"
          sub={`${fmt(p.reads_input - p.reads_passed)} removed`}
          color={p.quality_tier === "FAIL" ? "var(--destructive)" : p.quality_tier === "WARN" ? "var(--chart-4)" : "var(--chart-3)"}
        />
        {qa && (
          <StatCard
            value={pct(qa.duplication.estimated_duplication_rate)}
            label="Duplication"
            sub={`${fmt(qa.duplication.estimated_unique_sequences)} unique`}
          />
        )}
        {(p.pairs_passed ?? 0) > 0 && (
          <StatCard
            value={fmt(p.pairs_passed!)}
            label="Pairs"
            sub={`${fmt(p.pairs_merged || 0)} merged`}
          />
        )}
      </div>

      {/* Flags */}
      {p.flags.length > 0 ? (
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
      ) : (
        <div className="rounded-lg border border-chart-3/30 bg-chart-3/5 p-3 text-sm text-center text-muted-foreground">
          No quality flags raised
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

      {/* Adapters + Contaminants */}
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
                  { key: "rrna_removed", label: "rRNA", color: "text-chart-1" },
                  { key: "phix_removed", label: "PhiX", color: "text-chart-2" },
                  { key: "vector_removed", label: "Vector", color: "text-chart-4" },
                ].map(({ key, label, color }) => {
                  const val = getModuleExtra("contaminant", key);
                  if (!val) return null;
                  return (
                    <div key={key} className="text-center p-3 rounded-md bg-muted/50">
                      <div className={`text-lg font-bold font-mono ${color}`}>{pct(val / ri)}</div>
                      <div className="text-[10px] font-semibold text-muted-foreground uppercase">{label}</div>
                      <div className="text-[10px] text-muted-foreground">{fmt(val)} reads</div>
                    </div>
                  );
                })}
              </div>
            </CardContent>
          </Card>
        )}
      </div>

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
              <CardTitle className="text-sm">Read Length</CardTitle>
            </CardHeader>
            <CardContent>
              <HistogramChart data={qa.distributions.length_before} color="var(--chart-1)" />
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
