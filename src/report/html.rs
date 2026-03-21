//! HTML report generator -- self-contained single-file QC dashboard
//!
//! Generates an interactive HTML report from passport data using a
//! Catppuccin-inspired theme with shadcn-style components. Single file
//! with embedded data, inline CSS, and vanilla JS SVG charts.
//! Google Fonts loaded for Montserrat + Fira Code (degrades gracefully offline).

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

    // Build module funnel rows
    let mut funnel_rows = String::new();
    let mut reads_remaining = passport.reads_input;
    for m in &passport.modules {
        let pct = if passport.reads_input > 0 {
            m.reads_removed as f64 / passport.reads_input as f64 * 100.0
        } else {
            0.0
        };
        reads_remaining = reads_remaining.saturating_sub(m.reads_removed);
        let bar_width = (pct * 2.0).min(100.0); // scale for visual
        funnel_rows.push_str(&format!(
            r#"<tr>
              <td class="font-medium">{name}</td>
              <td class="num">{processed}</td>
              <td class="num">{removed}</td>
              <td class="num"><div class="bar-cell"><span class="mini-bar" style="width: {bar_width:.0}%"></span>{pct:.2}%</div></td>
              <td class="num">{remaining}</td>
            </tr>"#,
            name = m.name,
            processed = fmt(m.reads_processed),
            removed = fmt(m.reads_removed),
            bar_width = bar_width,
            pct = pct,
            remaining = fmt(reads_remaining),
        ));
    }

    // Contaminant detail cards
    let mut contaminant_cards = String::new();
    for m in &passport.modules {
        if m.name == "contaminant" {
            if let Some(rrna) = m.extra.get("rrna_removed").and_then(|v| v.as_u64()) {
                let prok = m.extra.get("rrna_prokaryotic").and_then(|v| v.as_u64()).unwrap_or(0);
                let euk = m.extra.get("rrna_eukaryotic").and_then(|v| v.as_u64()).unwrap_or(0);
                contaminant_cards.push_str(&format!(
                    r#"<div class="card"><div class="card-value chart-1">{}</div><div class="card-label">rRNA</div><div class="card-sub">{} prokaryotic / {} eukaryotic</div></div>"#,
                    fmt(rrna), fmt(prok), fmt(euk)
                ));
            }
            if let Some(phix) = m.extra.get("phix_removed").and_then(|v| v.as_u64()) {
                contaminant_cards.push_str(&format!(
                    r#"<div class="card"><div class="card-value chart-2">{}</div><div class="card-label">PhiX</div></div>"#,
                    fmt(phix)
                ));
            }
            if let Some(vec) = m.extra.get("vector_removed").and_then(|v| v.as_u64()) {
                contaminant_cards.push_str(&format!(
                    r#"<div class="card"><div class="card-value chart-4">{}</div><div class="card-label">Vector</div></div>"#,
                    fmt(vec)
                ));
            }
        }
    }

    // Flags
    let flags_html = if passport.flags.is_empty() {
        r#"<div class="card" style="text-align:center"><span class="badge pass">No quality flags raised</span></div>"#.to_string()
    } else {
        passport.flags.iter().map(|f| {
            let cls = format!("{:?}", f.severity).to_lowercase();
            format!(r#"<div class="alert {cls}"><strong>{}</strong> {}</div>"#, f.code, f.message)
        }).collect::<Vec<_>>().join("\n")
    };

    // Paired-end cards
    let paired_html = if passport.pairs_passed > 0 {
        format!(
            r#"<div class="card"><div class="card-value">{}</div><div class="card-label">Pairs passed</div></div>
            <div class="card"><div class="card-value">{}</div><div class="card-label">Singletons</div></div>
            <div class="card"><div class="card-value">{}</div><div class="card-label">Merged</div></div>"#,
            fmt(passport.pairs_passed), fmt(passport.singletons), fmt(passport.pairs_merged),
        )
    } else {
        String::new()
    };

    // Duplication info
    let dup_html = passport.qa_stats.as_ref().map(|qa| {
        let dup = &qa.duplication;
        format!(
            r#"<div class="card"><div class="card-value">{:.1}%</div><div class="card-label">Est. duplication</div><div class="card-sub">{} unique sequences</div></div>"#,
            dup.estimated_duplication_rate * 100.0,
            fmt(dup.estimated_unique_sequences),
        )
    }).unwrap_or_default();

    format!(r##"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>virome-qc Report</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Montserrat:wght@400;500;600;700&family=Fira+Code:wght@400;500&display=swap" rel="stylesheet">
<style>
:root {{
  --bg: oklch(0.9578 0.0058 264.5321);
  --fg: oklch(0.4355 0.0430 279.3250);
  --card: oklch(1.0000 0 0);
  --card-fg: oklch(0.4355 0.0430 279.3250);
  --primary: oklch(0.5547 0.2503 297.0156);
  --primary-fg: oklch(1.0000 0 0);
  --secondary: oklch(0.8575 0.0145 268.4756);
  --muted: oklch(0.9060 0.0117 264.5071);
  --muted-fg: oklch(0.5471 0.0343 279.0837);
  --accent: oklch(0.6820 0.1448 235.3822);
  --accent-fg: oklch(1.0000 0 0);
  --destructive: oklch(0.5505 0.2155 19.8095);
  --border: oklch(0.8083 0.0174 271.1982);
  --chart-1: oklch(0.5547 0.2503 297.0156);
  --chart-2: oklch(0.6820 0.1448 235.3822);
  --chart-3: oklch(0.6250 0.1772 140.4448);
  --chart-4: oklch(0.6920 0.2041 42.4293);
  --chart-5: oklch(0.7141 0.1045 33.0967);
  --pass: oklch(0.6250 0.1772 140.4448);
  --warn: oklch(0.6920 0.2041 42.4293);
  --fail: oklch(0.5505 0.2155 19.8095);
  --radius: 0.5rem;
  --shadow: 0px 4px 6px 0px hsl(240 30% 25% / 0.08), 0px 1px 2px -1px hsl(240 30% 25% / 0.08);
  --font-sans: 'Montserrat', system-ui, sans-serif;
  --font-mono: 'Fira Code', ui-monospace, monospace;
}}
*, *::before, *::after {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: var(--font-sans); background: var(--bg); color: var(--fg); line-height: 1.6; font-size: 14px; }}
.container {{ max-width: 1100px; margin: 0 auto; padding: 32px 24px; }}

/* Header */
.header {{ display: flex; justify-content: space-between; align-items: flex-start; margin-bottom: 32px; }}
.header h1 {{ font-size: 22px; font-weight: 700; letter-spacing: -0.02em; color: var(--fg); }}
.header-meta {{ color: var(--muted-fg); font-size: 13px; margin-top: 4px; }}
.badge {{ display: inline-flex; align-items: center; padding: 6px 16px; border-radius: 999px; font-weight: 600; font-size: 13px; letter-spacing: 0.02em; }}
.badge.pass {{ background: var(--pass); color: white; }}
.badge.warn {{ background: var(--warn); color: white; }}
.badge.fail {{ background: var(--fail); color: white; }}

/* Section headings */
h2 {{ font-size: 16px; font-weight: 600; margin: 36px 0 16px 0; padding-bottom: 8px; border-bottom: 1px solid var(--border); letter-spacing: -0.01em; }}

/* Cards */
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 12px; margin-bottom: 20px; }}
.card {{ background: var(--card); border: 1px solid var(--border); border-radius: var(--radius); padding: 16px; box-shadow: var(--shadow); text-align: center; }}
.card-value {{ font-size: 26px; font-weight: 700; font-family: var(--font-mono); color: var(--fg); }}
.card-value.chart-1 {{ color: var(--chart-1); }}
.card-value.chart-2 {{ color: var(--chart-2); }}
.card-value.chart-3 {{ color: var(--chart-3); }}
.card-value.chart-4 {{ color: var(--chart-4); }}
.card-value.primary {{ color: var(--primary); }}
.card-label {{ font-size: 12px; color: var(--muted-fg); margin-top: 2px; font-weight: 500; text-transform: uppercase; letter-spacing: 0.04em; }}
.card-sub {{ font-size: 11px; color: var(--muted-fg); margin-top: 4px; }}

/* Alerts / flags */
.alert {{ padding: 12px 16px; margin-bottom: 8px; border-radius: var(--radius); border: 1px solid; font-size: 13px; }}
.alert.warn {{ background: oklch(0.95 0.04 85); border-color: var(--warn); }}
.alert.fail {{ background: oklch(0.95 0.03 20); border-color: var(--fail); }}

/* Table */
table {{ width: 100%; border-collapse: collapse; background: var(--card); border: 1px solid var(--border); border-radius: var(--radius); overflow: hidden; box-shadow: var(--shadow); margin-bottom: 16px; }}
th {{ background: var(--muted); text-align: left; padding: 10px 14px; font-size: 12px; font-weight: 600; color: var(--muted-fg); text-transform: uppercase; letter-spacing: 0.04em; }}
td {{ padding: 9px 14px; border-top: 1px solid var(--border); font-size: 13px; }}
td.num {{ text-align: right; font-family: var(--font-mono); font-size: 12px; }}
.font-medium {{ font-weight: 500; }}
.bar-cell {{ display: flex; align-items: center; justify-content: flex-end; gap: 8px; }}
.mini-bar {{ display: inline-block; height: 8px; background: var(--primary); border-radius: 4px; opacity: 0.6; min-width: 0; }}

/* Charts */
.chart-panel {{ background: var(--card); border: 1px solid var(--border); border-radius: var(--radius); padding: 20px; box-shadow: var(--shadow); margin-bottom: 16px; }}
.chart-row {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
@media (max-width: 768px) {{ .chart-row {{ grid-template-columns: 1fr; }} }}
svg {{ width: 100%; height: auto; }}
.axis-label {{ font-size: 10px; fill: var(--muted-fg); font-family: var(--font-mono); }}
.axis-line {{ stroke: var(--border); stroke-width: 0.5; }}

/* Footer */
footer {{ margin-top: 48px; padding-top: 16px; border-top: 1px solid var(--border); color: var(--muted-fg); font-size: 12px; text-align: center; }}
</style>
</head>
<body>
<div class="container">

<div class="header">
  <div>
    <h1>virome-qc</h1>
    <div class="header-meta">{profile} &middot; v{version}</div>
  </div>
  <span class="badge {tier_class}">{tier}</span>
</div>

<div class="grid">
  <div class="card"><div class="card-value">{reads_input}</div><div class="card-label">Input reads</div></div>
  <div class="card"><div class="card-value">{reads_passed}</div><div class="card-label">Passed</div></div>
  <div class="card"><div class="card-value primary" style="color: var(--{tier_class})">{survival_pct:.1}%</div><div class="card-label">Survival</div></div>
  {dup_html}
  {paired_html}
</div>

{flags_html}

<h2>Survival Funnel</h2>
<table>
<tr><th>Module</th><th>Processed</th><th>Removed</th><th>% of Input</th><th>Remaining</th></tr>
{funnel_rows}
</table>

<h2>Contaminants</h2>
<div class="grid">
{contaminant_cards}
</div>

<h2>Quality Profiles</h2>
<div class="chart-row">
  <div class="chart-panel"><h3 style="font-size:13px;font-weight:600;margin-bottom:12px">Per-Position Quality (before QC)</h3><div id="chart-quality-before"></div></div>
  <div class="chart-panel"><h3 style="font-size:13px;font-weight:600;margin-bottom:12px">Per-Position Quality (after QC)</h3><div id="chart-quality-after"></div></div>
</div>

<h2>Base Composition</h2>
<div class="chart-row">
  <div class="chart-panel"><h3 style="font-size:13px;font-weight:600;margin-bottom:12px">Before QC</h3><div id="chart-bases-before"></div></div>
  <div class="chart-panel"><h3 style="font-size:13px;font-weight:600;margin-bottom:12px">After QC</h3><div id="chart-bases-after"></div></div>
</div>

<h2>Distributions</h2>
<div class="chart-row">
  <div class="chart-panel"><h3 style="font-size:13px;font-weight:600;margin-bottom:12px">Read Length</h3><div id="chart-length"></div></div>
  <div class="chart-panel"><h3 style="font-size:13px;font-weight:600;margin-bottom:12px">GC Content</h3><div id="chart-gc"></div></div>
</div>
<div class="chart-row">
  <div class="chart-panel"><h3 style="font-size:13px;font-weight:600;margin-bottom:12px">Quality Scores</h3><div id="chart-qscores"></div></div>
  <div class="chart-panel"><h3 style="font-size:13px;font-weight:600;margin-bottom:12px">Bases Trimmed</h3><div id="chart-trimmed"></div></div>
</div>

<footer>Generated by virome-qc v{version}</footer>
</div>

<script type="application/json" id="passport-data">
{passport_json}
</script>

<script>
const data = JSON.parse(document.getElementById('passport-data').textContent);
const qa = data.qa_stats || {{}};
const cs = getComputedStyle(document.documentElement);
const C = {{
  chart1: cs.getPropertyValue('--chart-1').trim() || '#7c3aed',
  chart2: cs.getPropertyValue('--chart-2').trim() || '#0ea5e9',
  chart3: cs.getPropertyValue('--chart-3').trim() || '#22c55e',
  chart4: cs.getPropertyValue('--chart-4').trim() || '#f59e0b',
  chart5: cs.getPropertyValue('--chart-5').trim() || '#f97316',
  border: cs.getPropertyValue('--border').trim() || '#e2e8f0',
  muted: cs.getPropertyValue('--muted-fg').trim() || '#64748b',
}};
const BASES = {{ a: C.chart3, c: C.chart2, g: C.chart4, t: C.chart5, n: C.muted }};

function svg(w, h) {{
  const s = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  s.setAttribute('viewBox', `0 0 ${{w}} ${{h}}`);
  return s;
}}
function el(tag, attrs) {{
  const e = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  return e;
}}

function drawHist(id, bins, counts, color) {{
  const box = document.getElementById(id);
  if (!box || !bins || !counts) return;
  const W=560, H=180, P={{l:45,r:15,t:8,b:30}};
  const s = svg(W, H);
  const mx = Math.max(...counts, 1);
  const n = counts.length, bw = (W-P.l-P.r)/n;
  // grid lines
  for (let i=0; i<=4; i++) {{
    const y = P.t + (H-P.t-P.b) * i / 4;
    s.appendChild(el('line', {{x1:P.l, x2:W-P.r, y1:y, y2:y, class:'axis-line'}}));
  }}
  for (let i=0; i<n; i++) {{
    const h = (counts[i]/mx)*(H-P.t-P.b);
    s.appendChild(el('rect', {{x:P.l+i*bw, y:H-P.b-h, width:Math.max(bw-1,1), height:h, fill:color, rx:1, opacity:0.85}}));
  }}
  const step = Math.max(1, Math.floor(n/8));
  for (let i=0; i<bins.length-1; i+=step) {{
    const t = el('text', {{x:P.l+i*bw+bw/2, y:H-4, class:'axis-label', 'text-anchor':'middle'}});
    t.textContent = typeof bins[i]==='number' ? (bins[i]<2 ? bins[i].toFixed(2) : Math.round(bins[i])) : bins[i];
    s.appendChild(t);
  }}
  box.appendChild(s);
}}

function drawQuality(id, positions) {{
  const box = document.getElementById(id);
  if (!box || !positions || !positions.length) return;
  const W=560, H=200, P={{l:45,r:15,t:12,b:30}};
  const s = svg(W, H);
  const n = positions.length;
  const xs = i => P.l + i*(W-P.l-P.r)/n;
  const ys = v => P.t + (1-v/42)*(H-P.t-P.b);
  // grid
  for (const q of [0,10,20,30,40]) {{
    const y = ys(q);
    s.appendChild(el('line', {{x1:P.l, x2:W-P.r, y1:y, y2:y, class:'axis-line'}}));
    const t = el('text', {{x:P.l-6, y:y+4, class:'axis-label', 'text-anchor':'end'}});
    t.textContent = 'Q'+q; s.appendChild(t);
  }}
  // Q25-Q75 band
  if (positions[0].q25 !== undefined) {{
    let d = '';
    for (let i=0;i<n;i++) d += (i?'L':'M')+xs(i).toFixed(1)+','+ys(positions[i].q75).toFixed(1);
    for (let i=n-1;i>=0;i--) d += 'L'+xs(i).toFixed(1)+','+ys(positions[i].q25).toFixed(1);
    s.appendChild(el('path', {{d:d+'Z', fill:C.chart1, opacity:0.15}}));
  }}
  // median line
  let md = '';
  for (let i=0;i<n;i++) md += (i?'L':'M')+xs(i).toFixed(1)+','+ys(positions[i].median||positions[i].mean).toFixed(1);
  s.appendChild(el('path', {{d:md, fill:'none', stroke:C.chart1, 'stroke-width':1.5}}));
  // mean line
  let ml = '';
  for (let i=0;i<n;i++) ml += (i?'L':'M')+xs(i).toFixed(1)+','+ys(positions[i].mean).toFixed(1);
  s.appendChild(el('path', {{d:ml, fill:'none', stroke:C.chart1, 'stroke-width':1, 'stroke-dasharray':'3,3', opacity:0.6}}));
  // x-axis
  const step = Math.max(1, Math.floor(n/10));
  for (let i=0;i<n;i+=step) {{
    const t = el('text', {{x:xs(i), y:H-4, class:'axis-label', 'text-anchor':'middle'}});
    t.textContent = positions[i].position; s.appendChild(t);
  }}
  box.appendChild(s);
}}

function drawBases(id, positions) {{
  const box = document.getElementById(id);
  if (!box || !positions || !positions.length) return;
  const W=560, H=160, P={{l:45,r:15,t:8,b:32}};
  const s = svg(W, H);
  const n = positions.length, bw = (W-P.l-P.r)/n, ch = H-P.t-P.b;
  for (let i=0;i<n;i++) {{
    const p = positions[i]; let y = P.t;
    for (const [b, c] of [['a',BASES.a],['c',BASES.c],['g',BASES.g],['t',BASES.t],['n',BASES.n]]) {{
      const h = (p[b]||0)*ch;
      if (h>0.3) s.appendChild(el('rect', {{x:P.l+i*bw, y:y, width:Math.max(bw,1), height:h, fill:c}}));
      y += h;
    }}
  }}
  // legend
  let lx = P.l;
  for (const [label, c] of [['A',BASES.a],['C',BASES.c],['G',BASES.g],['T',BASES.t],['N',BASES.n]]) {{
    s.appendChild(el('rect', {{x:lx, y:H-14, width:10, height:10, fill:c, rx:2}}));
    const t = el('text', {{x:lx+14, y:H-5, class:'axis-label'}});
    t.textContent = label; s.appendChild(t); lx += 36;
  }}
  box.appendChild(s);
}}

// Render
if (qa.per_position) {{
  drawQuality('chart-quality-before', (qa.per_position.quality_before||{{}}).positions);
  drawQuality('chart-quality-after', (qa.per_position.quality_after||{{}}).positions);
  drawBases('chart-bases-before', (qa.per_position.bases_before||{{}}).positions);
  drawBases('chart-bases-after', (qa.per_position.bases_after||{{}}).positions);
}}
if (qa.distributions) {{
  const d = qa.distributions;
  if (d.length_before) drawHist('chart-length', d.length_before.bin_edges, d.length_before.counts, C.chart1);
  if (d.gc_content) drawHist('chart-gc', d.gc_content.bin_edges, d.gc_content.counts, C.chart3);
  if (d.quality_scores) drawHist('chart-qscores', d.quality_scores.bin_edges, d.quality_scores.counts, C.chart2);
  if (d.trimmed_bases) drawHist('chart-trimmed', d.trimmed_bases.bin_edges, d.trimmed_bases.counts, C.chart4);
}}
</script>
</body>
</html>"##,
        profile = passport.profile,
        version = passport.tool_version,
        tier_class = tier_class,
        tier = format!("{:?}", passport.quality_tier).to_uppercase(),
        reads_input = fmt(passport.reads_input),
        reads_passed = fmt(passport.reads_passed),
        survival_pct = survival_pct,
        dup_html = dup_html,
        paired_html = paired_html,
        flags_html = flags_html,
        funnel_rows = funnel_rows,
        contaminant_cards = contaminant_cards,
        passport_json = passport_json,
    )
}

fn fmt(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.1}K", n as f64 / 1_000.0)
    } else {
        format!("{n}")
    }
}
