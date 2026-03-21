//! HTML report generator -- self-contained single-file QC dashboard
//!
//! Generates an interactive HTML report using Chart.js for interactive
//! charts with tooltips, hover, and legend toggling. Catppuccin theme
//! with dark mode support. Single file with CDN fallback.

use crate::report::Passport;
use anyhow::Result;
use std::path::Path;

/// Generate a self-contained HTML report from a passport
pub fn generate_html_report(passport: &Passport, output_path: &Path) -> Result<()> {
    let passport_json = serde_json::to_string(passport)?;
    let html = build_html(&passport_json, passport);
    std::fs::write(output_path, html)?;
    Ok(())
}

fn build_html(passport_json: &str, passport: &Passport) -> String {
    let survival_pct = passport.survival_rate * 100.0;
    let tier_class = match passport.quality_tier {
        crate::report::passport::QualityTier::Pass => "pass",
        crate::report::passport::QualityTier::Warn => "warn",
        crate::report::passport::QualityTier::Fail => "fail",
    };

    // Module funnel rows
    let mut funnel_rows = String::new();
    let mut reads_remaining = passport.reads_input;
    for m in &passport.modules {
        let pct = if passport.reads_input > 0 {
            m.reads_removed as f64 / passport.reads_input as f64 * 100.0
        } else { 0.0 };
        reads_remaining = reads_remaining.saturating_sub(m.reads_removed);
        let bar_w = (pct * 2.0).min(100.0);
        funnel_rows.push_str(&format!(
            r#"<tr><td class="font-medium">{}</td><td class="num">{}</td><td class="num">{}</td><td class="num"><div class="bar-cell"><span class="mini-bar" style="width:{:.0}%"></span>{:.2}%</div></td><td class="num">{}</td></tr>"#,
            m.name, fmt(m.reads_processed), fmt(m.reads_removed), bar_w, pct, fmt(reads_remaining),
        ));
    }

    // Contaminant cards
    let mut contaminant_cards = String::new();
    for m in &passport.modules {
        if m.name == "contaminant" {
            if let Some(rrna) = m.extra.get("rrna_removed").and_then(|v| v.as_u64()) {
                let prok = m.extra.get("rrna_prokaryotic").and_then(|v| v.as_u64()).unwrap_or(0);
                let euk = m.extra.get("rrna_eukaryotic").and_then(|v| v.as_u64()).unwrap_or(0);
                contaminant_cards.push_str(&format!(
                    r#"<div class="card"><div class="card-value" style="color:var(--chart-1)">{}</div><div class="card-label">rRNA</div><div class="card-sub">{} prokaryotic / {} eukaryotic</div></div>"#,
                    fmt(rrna), fmt(prok), fmt(euk)
                ));
            }
            if let Some(phix) = m.extra.get("phix_removed").and_then(|v| v.as_u64()) {
                contaminant_cards.push_str(&format!(
                    r#"<div class="card"><div class="card-value" style="color:var(--chart-2)">{}</div><div class="card-label">PhiX</div></div>"#, fmt(phix)
                ));
            }
            if let Some(vec) = m.extra.get("vector_removed").and_then(|v| v.as_u64()) {
                contaminant_cards.push_str(&format!(
                    r#"<div class="card"><div class="card-value" style="color:var(--chart-4)">{}</div><div class="card-label">Vector</div></div>"#, fmt(vec)
                ));
            }
        }
    }

    let flags_html = if passport.flags.is_empty() {
        r#"<div class="card" style="text-align:center"><span class="badge pass">No quality flags raised</span></div>"#.to_string()
    } else {
        passport.flags.iter().map(|f| {
            let cls = format!("{:?}", f.severity).to_lowercase();
            format!(r#"<div class="alert {cls}"><strong>{}</strong> {}</div>"#, f.code, f.message)
        }).collect::<Vec<_>>().join("\n")
    };

    let paired_html = if passport.pairs_passed > 0 {
        format!(
            r#"<div class="card"><div class="card-value">{}</div><div class="card-label">Pairs</div><div class="card-sub">{} pairs</div></div>
            <div class="card"><div class="card-value">{}</div><div class="card-label">Singletons</div></div>
            <div class="card"><div class="card-value">{}</div><div class="card-label">Merged</div></div>"#,
            fmt(passport.pairs_passed), passport.pairs_passed, fmt(passport.singletons), fmt(passport.pairs_merged),
        )
    } else { String::new() };

    let dup_html = passport.qa_stats.as_ref().map(|qa| {
        format!(
            r#"<div class="card"><div class="card-value">{:.1}%</div><div class="card-label">Duplication</div><div class="card-sub">{} unique</div></div>"#,
            qa.duplication.estimated_duplication_rate * 100.0, fmt(qa.duplication.estimated_unique_sequences),
        )
    }).unwrap_or_default();

    // Adapter breakdown cards
    let mut adapter_cards = String::new();
    if let Some(ref qa) = passport.qa_stats {
        for a in &qa.adapters.breakdown {
            let pct = if passport.reads_input > 0 { a.count as f64 / passport.reads_input as f64 * 100.0 } else { 0.0 };
            adapter_cards.push_str(&format!(
                r#"<div class="card"><div class="card-value" style="color:var(--chart-1)">{}</div><div class="card-label">{}</div><div class="card-sub">{:.2}% of reads</div></div>"#,
                fmt(a.count), a.name, pct
            ));
        }
    }
    if adapter_cards.is_empty() {
        adapter_cards = r#"<div class="card" style="text-align:center"><div class="card-sub">No adapters detected</div></div>"#.to_string();
    }

    // Contaminant section -- hide if not enabled
    let contaminant_section = if contaminant_cards.is_empty() {
        String::new()
    } else {
        format!(r#"<h2>Contaminants</h2><div class="grid">{contaminant_cards}</div>"#)
    };

    // Insert size visibility is handled in JS based on data

    format!(r##"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>virome-qc Report</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Montserrat:wght@400;500;600;700&family=Fira+Code:wght@400;500&display=swap" rel="stylesheet">
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.7/dist/chart.umd.min.js"></script>
<style>
:root {{
  --bg: oklch(0.9578 0.0058 264.5321); --fg: oklch(0.4355 0.0430 279.3250);
  --card: oklch(1.0000 0 0); --card-fg: oklch(0.4355 0.0430 279.3250);
  --primary: oklch(0.5547 0.2503 297.0156); --primary-fg: oklch(1.0000 0 0);
  --muted: oklch(0.9060 0.0117 264.5071); --muted-fg: oklch(0.5471 0.0343 279.0837);
  --border: oklch(0.8083 0.0174 271.1982);
  --chart-1: oklch(0.5547 0.2503 297.0156); --chart-2: oklch(0.6820 0.1448 235.3822);
  --chart-3: oklch(0.6250 0.1772 140.4448); --chart-4: oklch(0.6920 0.2041 42.4293);
  --chart-5: oklch(0.7141 0.1045 33.0967);
  --pass: oklch(0.6250 0.1772 140.4448); --warn: oklch(0.6920 0.2041 42.4293);
  --fail: oklch(0.5505 0.2155 19.8095);
  --radius: 0.5rem;
  --shadow: 0 1px 3px 0 rgb(0 0 0 / 0.06), 0 1px 2px -1px rgb(0 0 0 / 0.06);
  --font-sans: 'Montserrat', system-ui, sans-serif;
  --font-mono: 'Fira Code', ui-monospace, monospace;
}}
.dark {{
  --bg: oklch(0.2155 0.0254 284.0647); --fg: oklch(0.8787 0.0426 272.2767);
  --card: oklch(0.2429 0.0304 283.9110); --card-fg: oklch(0.8787 0.0426 272.2767);
  --primary: oklch(0.7871 0.1187 304.7693); --muted: oklch(0.2973 0.0294 276.2144);
  --muted-fg: oklch(0.7510 0.0396 273.9320); --border: oklch(0.3240 0.0319 281.9784);
  --chart-1: oklch(0.7871 0.1187 304.7693); --chart-2: oklch(0.8467 0.0833 210.2545);
  --chart-3: oklch(0.8577 0.1092 142.7153); --chart-4: oklch(0.8237 0.1015 52.6294);
  --chart-5: oklch(0.9226 0.0238 30.4919);
  --pass: oklch(0.8577 0.1092 142.7153); --warn: oklch(0.8237 0.1015 52.6294);
  --fail: oklch(0.7556 0.1297 2.7642);
  --shadow: 0 1px 3px 0 rgb(0 0 0 / 0.2), 0 1px 2px -1px rgb(0 0 0 / 0.2);
}}
*,*::before,*::after {{ margin:0; padding:0; box-sizing:border-box; }}
body {{ font-family:var(--font-sans); background:var(--bg); color:var(--fg); line-height:1.6; font-size:14px; transition:background .2s,color .2s; }}
.container {{ max-width:1100px; margin:0 auto; padding:32px 24px; }}
.header {{ display:flex; justify-content:space-between; align-items:flex-start; margin-bottom:32px; }}
.header h1 {{ font-size:22px; font-weight:700; letter-spacing:-0.02em; }}
.header-meta {{ color:var(--muted-fg); font-size:13px; margin-top:4px; }}
.header-actions {{ display:flex; align-items:center; gap:12px; }}
h2 {{ font-size:16px; font-weight:600; margin:36px 0 16px; padding-bottom:8px; border-bottom:1px solid var(--border); }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr)); gap:12px; margin-bottom:20px; }}
.card {{ background:var(--card); border:1px solid var(--border); border-radius:var(--radius); padding:16px; box-shadow:var(--shadow); text-align:center; transition:background .2s,border-color .2s; }}
.card-value {{ font-size:26px; font-weight:700; font-family:var(--font-mono); color:var(--fg); }}
.card-label {{ font-size:11px; color:var(--muted-fg); margin-top:2px; font-weight:600; text-transform:uppercase; letter-spacing:0.05em; }}
.card-sub {{ font-size:11px; color:var(--muted-fg); margin-top:4px; }}
.badge {{ display:inline-flex; padding:5px 14px; border-radius:999px; font-weight:600; font-size:13px; color:white; }}
.badge.pass {{ background:var(--pass); }} .badge.warn {{ background:var(--warn); }} .badge.fail {{ background:var(--fail); }}
.alert {{ padding:12px 16px; margin-bottom:8px; border-radius:var(--radius); border:1px solid; font-size:13px; }}
.alert.warn {{ background:oklch(0.95 0.04 85); border-color:var(--warn); }}
.alert.fail {{ background:oklch(0.95 0.03 20); border-color:var(--fail); }}
.dark .alert.warn {{ background:oklch(0.30 0.04 85); }} .dark .alert.fail {{ background:oklch(0.30 0.03 20); }}
table {{ width:100%; border-collapse:collapse; background:var(--card); border:1px solid var(--border); border-radius:var(--radius); overflow:hidden; box-shadow:var(--shadow); margin-bottom:16px; }}
th {{ background:var(--muted); text-align:left; padding:10px 14px; font-size:11px; font-weight:600; color:var(--muted-fg); text-transform:uppercase; letter-spacing:0.05em; }}
td {{ padding:8px 14px; border-top:1px solid var(--border); font-size:13px; }}
td.num {{ text-align:right; font-family:var(--font-mono); font-size:12px; }}
.font-medium {{ font-weight:500; }}
.bar-cell {{ display:flex; align-items:center; justify-content:flex-end; gap:8px; }}
.mini-bar {{ display:inline-block; height:6px; background:var(--primary); border-radius:3px; opacity:0.5; }}
.chart-panel {{ background:var(--card); border:1px solid var(--border); border-radius:var(--radius); padding:20px; box-shadow:var(--shadow); margin-bottom:16px; transition:background .2s,border-color .2s; }}
.chart-panel h3 {{ font-size:13px; font-weight:600; margin-bottom:12px; }}
.chart-row {{ display:grid; grid-template-columns:1fr 1fr; gap:16px; }}
@media(max-width:768px) {{ .chart-row {{ grid-template-columns:1fr; }} }}
canvas {{ max-height:260px; }}
.theme-toggle {{ background:var(--muted); border:1px solid var(--border); border-radius:var(--radius); padding:6px 12px; cursor:pointer; font-family:var(--font-sans); font-size:12px; font-weight:500; color:var(--muted-fg); display:inline-flex; align-items:center; gap:6px; transition:all .15s; }}
.theme-toggle:hover {{ background:var(--border); color:var(--fg); }}
footer {{ margin-top:48px; padding-top:16px; border-top:1px solid var(--border); color:var(--muted-fg); font-size:12px; text-align:center; }}
</style>
</head>
<body>
<div class="container">

<div class="header">
  <div>
    <h1>virome-qc</h1>
    <div class="header-meta">{profile} &middot; v{version}</div>
  </div>
  <div class="header-actions">
    <button class="theme-toggle" onclick="toggleTheme()" id="theme-btn"><span id="theme-label">Dark</span></button>
    <span class="badge {tier_class}">{tier}</span>
  </div>
</div>

<div class="grid">
  <div class="card"><div class="card-value">{reads_input}</div><div class="card-label">Input</div><div class="card-sub">{reads_input_exact} reads</div></div>
  <div class="card"><div class="card-value">{reads_passed}</div><div class="card-label">Passed</div><div class="card-sub">{reads_passed_exact} reads</div></div>
  <div class="card"><div class="card-value" style="color:var(--{tier_class})">{survival_pct:.1}%</div><div class="card-label">Survival</div><div class="card-sub">{reads_removed} removed</div></div>
  {dup_html}
  {paired_html}
</div>

{flags_html}

<h2>Survival Funnel</h2>
<table>
<tr><th>Module</th><th>Processed</th><th>Removed</th><th>% of Input</th><th>Remaining</th></tr>
{funnel_rows}
</table>

<h2>Adapters</h2>
<div class="grid">{adapter_cards}</div>

{contaminant_section}

<h2>Quality Profiles</h2>
<div class="chart-row">
  <div class="chart-panel"><h3>Before QC</h3><canvas id="chart-qbefore"></canvas></div>
  <div class="chart-panel"><h3>After QC</h3><canvas id="chart-qafter"></canvas></div>
</div>

<h2>Base Composition</h2>
<div class="chart-row">
  <div class="chart-panel"><h3>Before QC</h3><canvas id="chart-bbefore"></canvas></div>
  <div class="chart-panel"><h3>After QC</h3><canvas id="chart-bafter"></canvas></div>
</div>

<h2>Distributions</h2>
<div class="chart-row">
  <div class="chart-panel"><h3>Read Length (before/after)</h3><canvas id="chart-length"></canvas></div>
  <div class="chart-panel"><h3>GC Content</h3><canvas id="chart-gc"></canvas></div>
</div>
<div class="chart-row">
  <div class="chart-panel"><h3>Quality Scores</h3><canvas id="chart-qscores"></canvas></div>
  <div class="chart-panel" id="panel-insert" style="display:none"><h3>Insert Size</h3><canvas id="chart-insert"></canvas></div>
  <div class="chart-panel" id="panel-trimmed"><h3>Bases Trimmed</h3><canvas id="chart-trimmed"></canvas></div>
</div>

<footer>Generated by virome-qc v{version}</footer>
</div>

<script type="application/json" id="passport-data">{passport_json}</script>
<script>
const data = JSON.parse(document.getElementById('passport-data').textContent);
const qa = data.qa_stats || {{}};
const charts = [];

function css(v) {{ return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); }}
function colors() {{
  return {{
    c1: css('--chart-1'), c2: css('--chart-2'), c3: css('--chart-3'),
    c4: css('--chart-4'), c5: css('--chart-5'),
    fg: css('--fg'), muted: css('--muted-fg'), border: css('--border'), card: css('--card'),
  }};
}}

function chartDefaults(cl) {{
  return {{
    responsive: true,
    maintainAspectRatio: false,
    animation: {{ duration: 300 }},
    plugins: {{
      legend: {{ display: false }},
      tooltip: {{
        backgroundColor: cl.card,
        titleColor: cl.fg,
        bodyColor: cl.muted,
        borderColor: cl.border,
        borderWidth: 1,
        padding: 10,
        titleFont: {{ family: "'Montserrat'", weight: '600', size: 12 }},
        bodyFont: {{ family: "'Fira Code'", size: 11 }},
        cornerRadius: 6,
      }},
    }},
    scales: {{
      x: {{ ticks: {{ color: cl.muted, font: {{ family: "'Fira Code'", size: 10 }} }}, grid: {{ color: cl.border, drawBorder: false }} }},
      y: {{ ticks: {{ color: cl.muted, font: {{ family: "'Fira Code'", size: 10 }} }}, grid: {{ color: cl.border, drawBorder: false }} }},
    }},
  }};
}}

function makeQualityChart(canvasId, positions) {{
  if (!positions || !positions.length) return;
  const cl = colors();
  const labels = positions.map(p => p.position);
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;
  const chart = new Chart(ctx, {{
    type: 'line',
    data: {{
      labels,
      datasets: [
        {{ label: 'Q75', data: positions.map(p => p.q75), borderWidth: 0, backgroundColor: cl.c1 + '20', fill: '+1', pointRadius: 0 }},
        {{ label: 'Q25', data: positions.map(p => p.q25), borderWidth: 0, backgroundColor: 'transparent', fill: false, pointRadius: 0 }},
        {{ label: 'Median', data: positions.map(p => p.median || p.mean), borderColor: cl.c1, borderWidth: 2, fill: false, pointRadius: 0, tension: 0.2 }},
        {{ label: 'Mean', data: positions.map(p => p.mean), borderColor: cl.c1 + '60', borderWidth: 1, borderDash: [4,3], fill: false, pointRadius: 0, tension: 0.2 }},
      ],
    }},
    options: {{
      ...chartDefaults(cl),
      plugins: {{
        ...chartDefaults(cl).plugins,
        legend: {{ display: true, labels: {{ color: cl.muted, font: {{ family: "'Montserrat'", size: 11 }}, boxWidth: 14, padding: 12 }} }},
        tooltip: {{
          ...chartDefaults(cl).plugins.tooltip,
          callbacks: {{ label: (ctx) => `${{ctx.dataset.label}}: Q${{ctx.parsed.y.toFixed(1)}}` }},
        }},
      }},
      scales: {{
        ...chartDefaults(cl).scales,
        y: {{ ...chartDefaults(cl).scales.y, min: 0, max: 42, title: {{ display: true, text: 'Phred Score', color: cl.muted, font: {{ family: "'Montserrat'", size: 11 }} }} }},
        x: {{ ...chartDefaults(cl).scales.x, title: {{ display: true, text: 'Position', color: cl.muted, font: {{ family: "'Montserrat'", size: 11 }} }} }},
      }},
      interaction: {{ mode: 'index', intersect: false }},
    }},
  }});
  charts.push(chart);
}}

function makeBasesChart(canvasId, positions) {{
  if (!positions || !positions.length) return;
  const cl = colors();
  const labels = positions.map(p => p.position);
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;
  const chart = new Chart(ctx, {{
    type: 'bar',
    data: {{
      labels,
      datasets: [
        {{ label: 'A', data: positions.map(p => (p.a||0)*100), backgroundColor: cl.c3 }},
        {{ label: 'C', data: positions.map(p => (p.c||0)*100), backgroundColor: cl.c2 }},
        {{ label: 'G', data: positions.map(p => (p.g||0)*100), backgroundColor: cl.c4 }},
        {{ label: 'T', data: positions.map(p => (p.t||0)*100), backgroundColor: cl.c5 }},
        {{ label: 'N', data: positions.map(p => (p.n||0)*100), backgroundColor: cl.muted + '40' }},
      ],
    }},
    options: {{
      ...chartDefaults(cl),
      plugins: {{
        ...chartDefaults(cl).plugins,
        legend: {{ display: true, labels: {{ color: cl.muted, font: {{ family: "'Montserrat'", size: 11 }}, boxWidth: 12, padding: 10 }} }},
        tooltip: {{ ...chartDefaults(cl).plugins.tooltip, callbacks: {{ label: (ctx) => `${{ctx.dataset.label}}: ${{ctx.parsed.y.toFixed(1)}}%` }} }},
      }},
      scales: {{
        x: {{ ...chartDefaults(cl).scales.x, stacked: true, ticks: {{ ...chartDefaults(cl).scales.x.ticks, maxTicksLimit: 15 }} }},
        y: {{ ...chartDefaults(cl).scales.y, stacked: true, max: 100, title: {{ display: true, text: '%', color: cl.muted }} }},
      }},
    }},
  }});
  charts.push(chart);
}}

function makeHistChart(canvasId, dist, color) {{
  if (!dist || !dist.counts) return;
  const cl = colors();
  const edges = dist.bin_edges || [];
  const labels = dist.counts.map((_, i) => {{
    const v = edges[i];
    return v !== undefined ? (v < 2 ? v.toFixed(2) : Math.round(v).toString()) : i.toString();
  }});
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;
  const chart = new Chart(ctx, {{
    type: 'bar',
    data: {{
      labels,
      datasets: [{{ data: dist.counts, backgroundColor: color + 'cc', borderColor: color, borderWidth: 1, borderRadius: 2 }}],
    }},
    options: {{
      ...chartDefaults(cl),
      scales: {{
        x: {{ ...chartDefaults(cl).scales.x, ticks: {{ ...chartDefaults(cl).scales.x.ticks, maxTicksLimit: 12 }} }},
        y: {{ ...chartDefaults(cl).scales.y, title: {{ display: true, text: 'Reads', color: cl.muted }} }},
      }},
      plugins: {{
        ...chartDefaults(cl).plugins,
        tooltip: {{
          ...chartDefaults(cl).plugins.tooltip,
          callbacks: {{ title: (items) => items[0] ? `Bin: ${{items[0].label}}` : '', label: (ctx) => ctx.parsed.y.toLocaleString() + ' reads' }},
        }},
      }},
    }},
  }});
  charts.push(chart);
}}

function makeOverlayHistChart(canvasId, dist1, dist2, color1, color2, label1, label2) {{
  if (!dist1 || !dist1.counts) return;
  const cl = colors();
  const edges = dist1.bin_edges || [];
  const labels = dist1.counts.map((_, i) => {{
    const v = edges[i];
    return v !== undefined ? (v < 2 ? v.toFixed(2) : Math.round(v).toString()) : i.toString();
  }});
  const ctx = document.getElementById(canvasId);
  if (!ctx) return;
  const datasets = [{{ label: label1, data: dist1.counts, backgroundColor: color1 + '80', borderColor: color1, borderWidth: 1, borderRadius: 2 }}];
  if (dist2 && dist2.counts) {{
    datasets.push({{ label: label2, data: dist2.counts, backgroundColor: color2 + '80', borderColor: color2, borderWidth: 1, borderRadius: 2 }});
  }}
  const chart = new Chart(ctx, {{
    type: 'bar',
    data: {{ labels, datasets }},
    options: {{
      ...chartDefaults(cl),
      plugins: {{
        ...chartDefaults(cl).plugins,
        legend: {{ display: true, labels: {{ color: cl.muted, font: {{ family: "'Montserrat'", size: 11 }}, boxWidth: 12, padding: 10 }} }},
        tooltip: {{ ...chartDefaults(cl).plugins.tooltip, callbacks: {{ label: (ctx) => `${{ctx.dataset.label}}: ${{ctx.parsed.y.toLocaleString()}}` }} }},
      }},
      scales: {{
        x: {{ ...chartDefaults(cl).scales.x, ticks: {{ ...chartDefaults(cl).scales.x.ticks, maxTicksLimit: 12 }} }},
        y: {{ ...chartDefaults(cl).scales.y, title: {{ display: true, text: 'Reads', color: cl.muted }} }},
      }},
    }},
  }});
  charts.push(chart);
}}

function renderAll() {{
  charts.forEach(c => c.destroy());
  charts.length = 0;
  const cl = colors();
  if (qa.per_position) {{
    makeQualityChart('chart-qbefore', (qa.per_position.quality_before||{{}}).positions);
    makeQualityChart('chart-qafter', (qa.per_position.quality_after||{{}}).positions);
    makeBasesChart('chart-bbefore', (qa.per_position.bases_before||{{}}).positions);
    makeBasesChart('chart-bafter', (qa.per_position.bases_after||{{}}).positions);
  }}
  if (qa.distributions) {{
    const d = qa.distributions;
    // Overlay before/after read length
    if (d.length_before && d.length_after) {{
      makeOverlayHistChart('chart-length', d.length_before, d.length_after, cl.c1, cl.c3, 'Before', 'After');
    }} else if (d.length_before) {{
      makeHistChart('chart-length', d.length_before, cl.c1);
    }}
    makeHistChart('chart-gc', d.gc_content, cl.c3);
    makeHistChart('chart-qscores', d.quality_scores, cl.c2);
    // Insert size chart (show panel if data exists)
    if (d.insert_sizes && d.insert_sizes.total > 0) {{
      document.getElementById('panel-insert').style.display = '';
      makeHistChart('chart-insert', d.insert_sizes, cl.c2);
    }}
    makeHistChart('chart-trimmed', d.trimmed_bases, cl.c4);
  }}
}}

// Theme
function toggleTheme() {{
  document.documentElement.classList.toggle('dark');
  const isDark = document.documentElement.classList.contains('dark');
  localStorage.setItem('virome-qc-theme', isDark ? 'dark' : 'light');
  document.getElementById('theme-label').textContent = isDark ? 'Light' : 'Dark';
  renderAll();
}}
(function() {{
  const saved = localStorage.getItem('virome-qc-theme');
  const prefersDark = window.matchMedia('(prefers-color-scheme: dark)').matches;
  if (saved === 'dark' || (!saved && prefersDark)) {{
    document.documentElement.classList.add('dark');
    document.getElementById('theme-label').textContent = 'Light';
  }}
  Chart.defaults.font.family = "'Montserrat'";
  renderAll();
}})();
</script>
</body>
</html>"##,
        profile = passport.profile,
        version = passport.tool_version,
        tier_class = tier_class,
        tier = format!("{:?}", passport.quality_tier).to_uppercase(),
        reads_input = fmt(passport.reads_input),
        reads_input_exact = passport.reads_input,
        reads_passed = fmt(passport.reads_passed),
        reads_passed_exact = passport.reads_passed,
        reads_removed = fmt(passport.reads_input.saturating_sub(passport.reads_passed)),
        survival_pct = survival_pct,
        dup_html = dup_html,
        paired_html = paired_html,
        flags_html = flags_html,
        funnel_rows = funnel_rows,
        adapter_cards = adapter_cards,
        contaminant_section = contaminant_section,
        passport_json = passport_json,
    )
}

fn fmt(n: u64) -> String {
    if n >= 1_000_000 { format!("{:.1}M", n as f64 / 1_000_000.0) }
    else if n >= 1_000 { format!("{:.1}K", n as f64 / 1_000.0) }
    else { format!("{n}") }
}
