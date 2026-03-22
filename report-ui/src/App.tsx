import { Moon, Sun } from "lucide-react";
import { useState } from "react";
import { Badge } from "@/components/ui/badge";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { SampleReport } from "@/components/SampleReport";
import type { Passport } from "./types";
import { fmt, pct } from "@/lib/utils";

function ThemeToggle() {
  const [dark, setDark] = useState(
    document.documentElement.classList.contains("dark")
  );
  const toggle = () => {
    document.documentElement.classList.toggle("dark");
    const isDark = document.documentElement.classList.contains("dark");
    localStorage.setItem("virome-qc-theme", isDark ? "dark" : "light");
    setDark(isDark);
  };
  return (
    <button
      onClick={toggle}
      className="inline-flex items-center gap-2 rounded-md border border-border bg-secondary px-3 py-1.5 text-xs font-medium text-secondary-foreground hover:bg-muted transition-colors"
    >
      {dark ? <Sun className="h-3.5 w-3.5" /> : <Moon className="h-3.5 w-3.5" />}
      {dark ? "Light" : "Dark"}
    </button>
  );
}

function tierVariant(tier: string): "pass" | "warn" | "fail" {
  return tier.toLowerCase() as "pass" | "warn" | "fail";
}

export function App({
  data,
  isBatch,
}: {
  data: Passport | Passport[] | null;
  isBatch: boolean;
}) {
  if (!data) {
    return (
      <div className="flex items-center justify-center h-screen">
        <p className="text-muted-foreground">No passport data found.</p>
      </div>
    );
  }

  if (isBatch) {
    return <BatchView passports={data as Passport[]} />;
  }

  return <SingleView passport={data as Passport} />;
}

function SingleView({ passport }: { passport: Passport }) {
  const sampleName =
    passport.provenance?.input_files?.[0]?.path
      ?.split("/")
      .pop()
      ?.replace(/\.fastq\.gz$|\.fq\.gz$|\.fastq$/, "") || "Sample";

  return (
    <div className="min-h-screen bg-background">
      <div className="max-w-6xl mx-auto px-4 py-6 space-y-6">
        {/* Header */}
        <div className="flex items-start justify-between">
          <div>
            <h1 className="text-xl font-bold tracking-tight">{sampleName}</h1>
            <p className="text-sm text-muted-foreground">
              {passport.profile} -- virome-qc v{passport.tool_version}
            </p>
          </div>
          <div className="flex items-center gap-3">
            <ThemeToggle />
            <Badge variant={tierVariant(passport.quality_tier)}>
              {passport.quality_tier}
            </Badge>
          </div>
        </div>

        <SampleReport passport={passport} />
      </div>
    </div>
  );
}

function BatchView({ passports }: { passports: Passport[] }) {
  return (
    <div className="min-h-screen bg-background">
      <div className="max-w-6xl mx-auto px-4 py-6 space-y-6">
        <div className="flex items-start justify-between">
          <div>
            <h1 className="text-xl font-bold tracking-tight">
              Batch Report
            </h1>
            <p className="text-sm text-muted-foreground">
              {passports.length} samples
            </p>
          </div>
          <ThemeToggle />
        </div>

        {/* Summary table */}
        <Card>
          <CardHeader>
            <CardTitle className="text-sm">Samples</CardTitle>
          </CardHeader>
          <CardContent>
            <div className="overflow-x-auto">
              <table className="w-full text-xs">
                <thead>
                  <tr className="border-b border-border">
                    <th className="text-left py-2 px-3 font-semibold text-muted-foreground">Sample</th>
                    <th className="text-right py-2 px-3 font-semibold text-muted-foreground">Input</th>
                    <th className="text-right py-2 px-3 font-semibold text-muted-foreground">Passed</th>
                    <th className="text-right py-2 px-3 font-semibold text-muted-foreground">Survival</th>
                    <th className="text-right py-2 px-3 font-semibold text-muted-foreground">Dedup</th>
                    <th className="text-right py-2 px-3 font-semibold text-muted-foreground">rRNA</th>
                    <th className="text-center py-2 px-3 font-semibold text-muted-foreground">Tier</th>
                  </tr>
                </thead>
                <tbody>
                  {passports.map((p, i) => {
                    const name =
                      p.provenance?.input_files?.[0]?.path?.split("/").pop() ||
                      `Sample ${i + 1}`;
                    const dedupRate =
                      p.modules.find((m) => m.name === "dedup")?.reads_removed ||
                      0;
                    const rrnaRate =
                      (p.modules
                        .find((m) => m.name === "contaminant")
                        ?.extra?.rrna_removed as number) || 0;
                    return (
                      <tr key={i} className="border-b border-border/50 hover:bg-muted/50">
                        <td className="py-2 px-3 font-medium">{name}</td>
                        <td className="py-2 px-3 text-right font-mono">{fmt(p.reads_input)}</td>
                        <td className="py-2 px-3 text-right font-mono">{fmt(p.reads_passed)}</td>
                        <td className="py-2 px-3 text-right font-mono">{pct(p.survival_rate)}</td>
                        <td className="py-2 px-3 text-right font-mono">
                          {pct(dedupRate / Math.max(p.reads_input, 1))}
                        </td>
                        <td className="py-2 px-3 text-right font-mono">
                          {pct(rrnaRate / Math.max(p.reads_input, 1))}
                        </td>
                        <td className="py-2 px-3 text-center">
                          <Badge variant={tierVariant(p.quality_tier)} className="text-[10px]">
                            {p.quality_tier}
                          </Badge>
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
          </CardContent>
        </Card>

        {/* Individual sample reports */}
        {passports.map((p, i) => {
          const name =
            p.provenance?.input_files?.[0]?.path?.split("/").pop() ||
            `Sample ${i + 1}`;
          return (
            <div key={i}>
              <h2 className="text-base font-semibold mt-8 mb-3 flex items-center gap-2">
                {name}
                <Badge variant={tierVariant(p.quality_tier)} className="text-[10px]">
                  {p.quality_tier}
                </Badge>
              </h2>
              <SampleReport passport={p} />
            </div>
          );
        })}
      </div>
    </div>
  );
}
