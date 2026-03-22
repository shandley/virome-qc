//! Multi-sample batch report generator
//!
//! Reads passport.json files from multiple samples and generates a single
//! interactive HTML dashboard for study-level QC comparison.

use crate::report::Passport;
use anyhow::Result;
use std::path::Path;

/// Summary of one sample for the batch table
#[derive(Debug, serde::Serialize)]
struct SampleSummary {
    name: String,
    reads_input: u64,
    reads_passed: u64,
    survival_rate: f64,
    adapter_rate: f64,
    dedup_rate: f64,
    rrna_rate: f64,
    host_rate: f64,
    quality_tier: String,
    mean_gc: f64,
}

/// Generate a batch report from a directory of passport files
pub fn generate_batch_report(input_dir: &Path, output_path: &Path) -> Result<()> {
    // Collect all passport files
    let mut passports: Vec<(String, Passport)> = Vec::new();

    for entry in std::fs::read_dir(input_dir)? {
        let entry = entry?;
        let path = entry.path();

        // Look for passport.json in subdirectories or directly
        let passport_path = if path.is_dir() {
            path.join("passport.json")
        } else if path.extension().map(|e| e == "json").unwrap_or(false) {
            path.clone()
        } else {
            continue;
        };

        if !passport_path.exists() {
            continue;
        }

        let contents = std::fs::read_to_string(&passport_path)?;
        if let Ok(passport) = serde_json::from_str::<Passport>(&contents) {
            let name = path
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string();
            passports.push((name, passport));
        }
    }

    if passports.is_empty() {
        anyhow::bail!("No passport.json files found in {}", input_dir.display());
    }

    passports.sort_by(|a, b| a.0.cmp(&b.0));

    // Extract per-sample summaries
    let summaries: Vec<SampleSummary> = passports
        .iter()
        .map(|(name, p)| {
            let ri = p.reads_input.max(1) as f64;

            let mut adapter_rate = 0.0;
            let mut dedup_rate = 0.0;
            let mut rrna_rate = 0.0;
            let mut host_rate = 0.0;

            for m in &p.modules {
                match m.name.as_str() {
                    "adapter" => {
                        adapter_rate = m
                            .extra
                            .get("adapters_found_3prime")
                            .and_then(|v| v.as_u64())
                            .unwrap_or(0) as f64
                            / ri;
                    }
                    "dedup" => {
                        dedup_rate = m.reads_removed as f64 / ri;
                    }
                    "contaminant" => {
                        rrna_rate = m
                            .extra
                            .get("rrna_removed")
                            .and_then(|v| v.as_u64())
                            .unwrap_or(0) as f64
                            / ri;
                    }
                    "host" => {
                        host_rate = m.reads_removed as f64 / ri;
                    }
                    _ => {}
                }
            }

            let mean_gc = p
                .qa_stats
                .as_ref()
                .map(|qa| {
                    let gc = &qa.distributions.gc_content;
                    if gc.total == 0 {
                        return 0.0;
                    }
                    let mut sum = 0.0;
                    let mut count = 0u64;
                    for (i, &c) in gc.counts.iter().enumerate() {
                        if i < gc.bin_edges.len() {
                            sum += gc.bin_edges[i] * c as f64;
                            count += c;
                        }
                    }
                    if count > 0 {
                        sum / count as f64
                    } else {
                        0.0
                    }
                })
                .unwrap_or(0.0);

            SampleSummary {
                name: name.clone(),
                reads_input: p.reads_input,
                reads_passed: p.reads_passed,
                survival_rate: p.survival_rate,
                adapter_rate,
                dedup_rate,
                rrna_rate,
                host_rate,
                quality_tier: format!("{:?}", p.quality_tier).to_uppercase(),
                mean_gc,
            }
        })
        .collect();

    let summaries_json = serde_json::to_string(&summaries)?;
    let n = summaries.len();

    let html = format!(
        r##"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>virome-qc Batch Report ({n} samples)</title>
<link href="https://fonts.googleapis.com/css2?family=Montserrat:wght@400;500;600;700&family=Fira+Code:wght@400;500&display=swap" rel="stylesheet">
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.7/dist/chart.umd.min.js"></script>
<style>
:root {{
  --bg: oklch(0.9578 0.0058 264.5321); --fg: oklch(0.4355 0.0430 279.3250);
  --card: oklch(1.0000 0 0); --border: oklch(0.8083 0.0174 271.1982);
  --muted: oklch(0.9060 0.0117 264.5071); --muted-fg: oklch(0.5471 0.0343 279.0837);
  --chart-1: oklch(0.5547 0.2503 297.0156); --chart-2: oklch(0.6820 0.1448 235.3822);
  --chart-3: oklch(0.6250 0.1772 140.4448); --chart-4: oklch(0.6920 0.2041 42.4293);
  --chart-5: oklch(0.7141 0.1045 33.0967);
  --pass: oklch(0.6250 0.1772 140.4448); --warn: oklch(0.6920 0.2041 42.4293);
  --fail: oklch(0.5505 0.2155 19.8095);
  --radius: 0.5rem; --shadow: 0 1px 3px 0 rgb(0 0 0 / 0.06);
  --font-sans: 'Montserrat', system-ui, sans-serif; --font-mono: 'Fira Code', monospace;
}}
.dark {{
  --bg: oklch(0.2155 0.0254 284.0647); --fg: oklch(0.8787 0.0426 272.2767);
  --card: oklch(0.2429 0.0304 283.9110); --border: oklch(0.3240 0.0319 281.9784);
  --muted: oklch(0.2973 0.0294 276.2144); --muted-fg: oklch(0.7510 0.0396 273.9320);
  --chart-1: oklch(0.7871 0.1187 304.7693); --chart-2: oklch(0.8467 0.0833 210.2545);
  --chart-3: oklch(0.8577 0.1092 142.7153); --chart-4: oklch(0.8237 0.1015 52.6294);
  --pass: oklch(0.8577 0.1092 142.7153); --warn: oklch(0.8237 0.1015 52.6294);
  --fail: oklch(0.7556 0.1297 2.7642); --shadow: 0 1px 3px 0 rgb(0 0 0 / 0.2);
}}
*,*::before,*::after {{ margin:0;padding:0;box-sizing:border-box; }}
body {{ font-family:var(--font-sans); background:var(--bg); color:var(--fg); font-size:13px; transition:background .2s,color .2s; }}
.container {{ max-width:1300px; margin:0 auto; padding:24px; }}
h1 {{ font-size:22px; font-weight:700; letter-spacing:-0.02em; }}
h2 {{ font-size:15px; font-weight:600; margin:28px 0 12px; padding-bottom:6px; border-bottom:1px solid var(--border); }}
.header {{ display:flex; justify-content:space-between; align-items:center; margin-bottom:24px; }}
.header-meta {{ color:var(--muted-fg); font-size:13px; margin-top:4px; }}
.theme-toggle {{ background:var(--muted); border:1px solid var(--border); border-radius:var(--radius); padding:6px 12px; cursor:pointer; font-family:var(--font-sans); font-size:12px; font-weight:500; color:var(--muted-fg); display:inline-flex; align-items:center; gap:6px; transition:all .15s; }}
.theme-toggle:hover {{ background:var(--border); color:var(--fg); }}
table {{ width:100%; border-collapse:collapse; background:var(--card); border:1px solid var(--border); border-radius:var(--radius); overflow:hidden; box-shadow:var(--shadow); font-size:12px; }}
th {{ background:var(--muted); text-align:left; padding:8px 10px; font-size:10px; font-weight:600; color:var(--muted-fg); text-transform:uppercase; letter-spacing:0.05em; cursor:pointer; user-select:none; }}
th:hover {{ background:var(--border); }}
td {{ padding:6px 10px; border-top:1px solid var(--border); }}
td.num {{ text-align:right; font-family:var(--font-mono); font-size:11px; }}
.tier {{ display:inline-block; padding:2px 8px; border-radius:10px; font-size:10px; font-weight:600; color:white; }}
.tier-PASS {{ background:var(--pass); }} .tier-WARN {{ background:var(--warn); }} .tier-FAIL {{ background:var(--fail); }}
.chart-panel {{ background:var(--card); border:1px solid var(--border); border-radius:var(--radius); padding:20px; box-shadow:var(--shadow); margin-bottom:16px; transition:background .2s,border-color .2s; }}
.chart-panel h3 {{ font-size:13px; font-weight:600; margin-bottom:12px; }}
.chart-row {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; }}
@media(max-width:768px) {{ .chart-row {{ grid-template-columns:1fr; }} }}
canvas {{ max-height:260px; }}
footer {{ margin-top:32px; padding-top:12px; border-top:1px solid var(--border); color:var(--muted-fg); font-size:11px; text-align:center; }}
</style>
</head>
<body>
<div class="container">
<div class="header">
  <div><h1>virome-qc Batch Report</h1><div class="header-meta">{n} samples</div></div>
  <button class="theme-toggle" onclick="toggleTheme()" id="theme-btn"><span id="theme-label">Dark</span></button>
</div>

<h2>Sample Overview</h2>
<table id="sample-table">
<thead>
<tr>
  <th onclick="sortTable(0)">Sample</th>
  <th onclick="sortTable(1)">Input</th>
  <th onclick="sortTable(2)">Passed</th>
  <th onclick="sortTable(3)">Survival</th>
  <th onclick="sortTable(4)">Adapter</th>
  <th onclick="sortTable(5)">Dedup</th>
  <th onclick="sortTable(6)">rRNA</th>
  <th onclick="sortTable(7)">Host</th>
  <th onclick="sortTable(8)">GC</th>
  <th onclick="sortTable(9)">Tier</th>
</tr>
</thead>
<tbody id="table-body"></tbody>
</table>

<h2>Cross-Sample Distributions</h2>
<div class="chart-row">
  <div class="chart-panel"><h3>Survival Rate</h3><canvas id="ch-survival"></canvas></div>
  <div class="chart-panel"><h3>Duplication Rate</h3><canvas id="ch-dedup"></canvas></div>
  <div class="chart-panel"><h3>rRNA Fraction</h3><canvas id="ch-rrna"></canvas></div>
</div>
<div class="chart-row">
  <div class="chart-panel"><h3>Adapter Rate</h3><canvas id="ch-adapter"></canvas></div>
  <div class="chart-panel"><h3>Host Fraction</h3><canvas id="ch-host"></canvas></div>
  <div class="chart-panel"><h3>GC Content</h3><canvas id="ch-gc"></canvas></div>
</div>

<footer>Generated by virome-qc</footer>
</div>

<script>
const samples = {summaries_json};
const charts = [];

function fmt(n) {{ return n >= 1e6 ? (n/1e6).toFixed(1)+'M' : n >= 1e3 ? (n/1e3).toFixed(1)+'K' : n.toString(); }}
function pct(v) {{ return (v*100).toFixed(1)+'%'; }}
function css(v) {{ return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); }}

function populateTable() {{
  const body = document.getElementById('table-body');
  body.innerHTML = '';
  samples.forEach(s => {{
    const tr = document.createElement('tr');
    tr.innerHTML = `<td>${{s.name}}</td><td class="num">${{fmt(s.reads_input)}}</td><td class="num">${{fmt(s.reads_passed)}}</td><td class="num">${{pct(s.survival_rate)}}</td><td class="num">${{pct(s.adapter_rate)}}</td><td class="num">${{pct(s.dedup_rate)}}</td><td class="num">${{pct(s.rrna_rate)}}</td><td class="num">${{pct(s.host_rate)}}</td><td class="num">${{(s.mean_gc*100).toFixed(1)}}%</td><td><span class="tier tier-${{s.quality_tier}}">${{s.quality_tier}}</span></td>`;
    body.appendChild(tr);
  }});
}}

let sortCol = -1, sortAsc = true;
function sortTable(col) {{
  if (sortCol === col) sortAsc = !sortAsc;
  else {{ sortCol = col; sortAsc = true; }}
  const keys = ['name','reads_input','reads_passed','survival_rate','adapter_rate','dedup_rate','rrna_rate','host_rate','mean_gc','quality_tier'];
  samples.sort((a,b) => {{
    let va = a[keys[col]], vb = b[keys[col]];
    if (typeof va === 'string') return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
    return sortAsc ? va - vb : vb - va;
  }});
  populateTable();
}}

function makeBarChart(id, values, color, label) {{
  const el = document.getElementById(id);
  if (!el) return;
  const c1 = css(color);
  const labels = samples.map(s => s.name);
  const ch = new Chart(el, {{
    type: 'bar',
    data: {{ labels, datasets: [{{ data: values, backgroundColor: c1+'cc', borderColor: c1, borderWidth: 1, borderRadius: 2 }}] }},
    options: {{
      responsive: true, maintainAspectRatio: false, animation: {{ duration: 200 }},
      plugins: {{ legend: {{ display: false }}, tooltip: {{
        backgroundColor: css('--card'), titleColor: css('--fg'), bodyColor: css('--muted-fg'),
        borderColor: css('--border'), borderWidth: 1, padding: 8, cornerRadius: 6,
        callbacks: {{ label: ctx => label + ': ' + (ctx.parsed.y*100).toFixed(2) + '%' }}
      }} }},
      scales: {{
        x: {{ display: samples.length <= 30, ticks: {{ color: css('--muted-fg'), font: {{ size: 9 }}, maxRotation: 45 }}, grid: {{ display: false }} }},
        y: {{ ticks: {{ color: css('--muted-fg'), callback: v => (v*100).toFixed(0)+'%' }}, grid: {{ color: css('--border') }} }}
      }}
    }}
  }});
  charts.push(ch);
}}

function renderCharts() {{
  charts.forEach(c => c.destroy()); charts.length = 0;
  makeBarChart('ch-survival', samples.map(s => s.survival_rate), '--chart-1', 'Survival');
  makeBarChart('ch-dedup', samples.map(s => s.dedup_rate), '--chart-2', 'Dedup');
  makeBarChart('ch-rrna', samples.map(s => s.rrna_rate), '--chart-3', 'rRNA');
  makeBarChart('ch-adapter', samples.map(s => s.adapter_rate), '--chart-4', 'Adapter');
  makeBarChart('ch-host', samples.map(s => s.host_rate), '--chart-5', 'Host');
  makeBarChart('ch-gc', samples.map(s => s.mean_gc), '--chart-1', 'GC');
}}

function toggleTheme() {{
  document.documentElement.classList.toggle('dark');
  const isDark = document.documentElement.classList.contains('dark');
  localStorage.setItem('virome-qc-theme', isDark ? 'dark' : 'light');
  document.getElementById('theme-label').textContent = isDark ? 'Light' : 'Dark';
  renderCharts();
}}

(function() {{
  const saved = localStorage.getItem('virome-qc-theme');
  const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
  if (saved === 'dark' || (!saved && prefersDark)) {{
    document.documentElement.classList.add('dark');
    document.getElementById('theme-label').textContent = 'Light';
  }}
  Chart.defaults.font.family = "'Montserrat'";
  populateTable();
  renderCharts();
}})();
</script>
</body>
</html>"##,
        n = n,
        summaries_json = summaries_json,
    );

    std::fs::write(output_path, html)?;
    Ok(())
}
