//! HTML report generator -- self-contained single-file QC dashboard
//!
//! Generates an interactive HTML report from passport data. The report
//! is a single file with embedded data and inline SVG charts. No external
//! dependencies (no CDN, no JavaScript libraries).

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
    let tier_color = match passport.quality_tier {
        crate::report::passport::QualityTier::Pass => "#22c55e",
        crate::report::passport::QualityTier::Warn => "#eab308",
        crate::report::passport::QualityTier::Fail => "#ef4444",
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
        funnel_rows.push_str(&format!(
            "<tr><td>{}</td><td class='num'>{}</td><td class='num'>{}</td><td class='num'>{:.2}%</td><td class='num'>{}</td></tr>\n",
            m.name,
            format_number(m.reads_processed),
            format_number(m.reads_removed),
            pct,
            format_number(reads_remaining),
        ));
    }

    // Build contaminant detail
    let mut contaminant_detail = String::new();
    for m in &passport.modules {
        if m.name == "contaminant" {
            if let Some(rrna) = m.extra.get("rrna_removed").and_then(|v| v.as_u64()) {
                let prok = m.extra.get("rrna_prokaryotic").and_then(|v| v.as_u64()).unwrap_or(0);
                let euk = m.extra.get("rrna_eukaryotic").and_then(|v| v.as_u64()).unwrap_or(0);
                contaminant_detail.push_str(&format!(
                    "<div class='stat-card'><div class='stat-value'>{}</div><div class='stat-label'>rRNA removed<br><small>{} prokaryotic, {} eukaryotic</small></div></div>",
                    format_number(rrna), format_number(prok), format_number(euk)
                ));
            }
            if let Some(phix) = m.extra.get("phix_removed").and_then(|v| v.as_u64()) {
                contaminant_detail.push_str(&format!(
                    "<div class='stat-card'><div class='stat-value'>{}</div><div class='stat-label'>PhiX removed</div></div>",
                    format_number(phix)
                ));
            }
            if let Some(vec) = m.extra.get("vector_removed").and_then(|v| v.as_u64()) {
                contaminant_detail.push_str(&format!(
                    "<div class='stat-card'><div class='stat-value'>{}</div><div class='stat-label'>Vector removed</div></div>",
                    format_number(vec)
                ));
            }
        }
    }

    // Build flags section
    let flags_html = if passport.flags.is_empty() {
        "<p class='no-flags'>No quality flags raised.</p>".to_string()
    } else {
        passport.flags.iter().map(|f| {
            let severity_class = format!("{:?}", f.severity).to_lowercase();
            format!("<div class='flag {severity_class}'><strong>{}</strong>: {}</div>", f.code, f.message)
        }).collect::<Vec<_>>().join("\n")
    };

    // Paired-end section
    let paired_html = if passport.pairs_passed > 0 {
        format!(
            "<div class='stat-card'><div class='stat-value'>{}</div><div class='stat-label'>Pairs passed</div></div>\
             <div class='stat-card'><div class='stat-value'>{}</div><div class='stat-label'>Singletons</div></div>\
             <div class='stat-card'><div class='stat-value'>{}</div><div class='stat-label'>Pairs merged</div></div>",
            format_number(passport.pairs_passed),
            format_number(passport.singletons),
            format_number(passport.pairs_merged),
        )
    } else {
        String::new()
    };

    format!(r##"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>virome-qc Report</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: #f8fafc; color: #1e293b; line-height: 1.5; }}
.container {{ max-width: 1200px; margin: 0 auto; padding: 24px; }}
h1 {{ font-size: 24px; font-weight: 600; margin-bottom: 4px; }}
h2 {{ font-size: 18px; font-weight: 600; margin: 32px 0 16px 0; padding-bottom: 8px; border-bottom: 2px solid #e2e8f0; }}
h3 {{ font-size: 15px; font-weight: 600; margin: 20px 0 8px 0; }}
.header {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 24px; }}
.header-meta {{ color: #64748b; font-size: 14px; }}
.tier-badge {{ display: inline-block; padding: 4px 16px; border-radius: 20px; font-weight: 700; font-size: 14px; color: white; }}
.stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 16px; margin-bottom: 24px; }}
.stat-card {{ background: white; border-radius: 8px; padding: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); text-align: center; }}
.stat-value {{ font-size: 28px; font-weight: 700; color: #0f172a; }}
.stat-label {{ font-size: 13px; color: #64748b; margin-top: 4px; }}
table {{ width: 100%; border-collapse: collapse; background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.1); margin-bottom: 16px; }}
th {{ background: #f1f5f9; text-align: left; padding: 10px 12px; font-size: 13px; font-weight: 600; color: #475569; }}
td {{ padding: 8px 12px; border-top: 1px solid #f1f5f9; font-size: 14px; }}
td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
.flag {{ padding: 10px 14px; margin-bottom: 8px; border-radius: 6px; font-size: 14px; }}
.flag.warn {{ background: #fef9c3; border-left: 4px solid #eab308; }}
.flag.fail {{ background: #fee2e2; border-left: 4px solid #ef4444; }}
.no-flags {{ color: #22c55e; font-weight: 500; }}
.chart-container {{ background: white; border-radius: 8px; padding: 16px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); margin-bottom: 16px; }}
svg {{ width: 100%; height: auto; }}
.bar {{ fill: #3b82f6; }}
.bar-before {{ fill: #94a3b8; }}
.bar-after {{ fill: #3b82f6; }}
.axis-label {{ font-size: 11px; fill: #64748b; }}
.chart-title {{ font-size: 14px; font-weight: 600; fill: #1e293b; }}
footer {{ margin-top: 40px; padding-top: 16px; border-top: 1px solid #e2e8f0; color: #94a3b8; font-size: 12px; text-align: center; }}
</style>
</head>
<body>
<div class="container">

<div class="header">
  <div>
    <h1>virome-qc Quality Report</h1>
    <div class="header-meta">Profile: {profile} | virome-qc v{version}</div>
  </div>
  <div>
    <span class="tier-badge" style="background: {tier_color}">{tier}</span>
  </div>
</div>

<div class="stats-grid">
  <div class="stat-card"><div class="stat-value">{reads_input}</div><div class="stat-label">Input reads</div></div>
  <div class="stat-card"><div class="stat-value">{reads_passed}</div><div class="stat-label">Reads passed</div></div>
  <div class="stat-card"><div class="stat-value" style="color: {tier_color}">{survival_pct:.1}%</div><div class="stat-label">Survival rate</div></div>
  {paired_html}
</div>

<h2>Quality Flags</h2>
{flags_html}

<h2>Survival Funnel</h2>
<table>
<tr><th>Module</th><th>Processed</th><th>Removed</th><th>% of Input</th><th>Remaining</th></tr>
{funnel_rows}
</table>

<h2>Contaminant Breakdown</h2>
<div class="stats-grid">
{contaminant_detail}
</div>

<h2>Per-Position Quality</h2>
<div class="chart-container" id="chart-quality"></div>

<h2>Per-Position Base Composition</h2>
<div class="chart-container" id="chart-bases"></div>

<h2>Read Length Distribution</h2>
<div class="chart-container" id="chart-length"></div>

<h2>GC Content Distribution</h2>
<div class="chart-container" id="chart-gc"></div>

<footer>Generated by virome-qc v{version}</footer>

</div>

<script type="application/json" id="passport-data">
{passport_json}
</script>

<script>
const data = JSON.parse(document.getElementById('passport-data').textContent);
const qa = data.qa_stats || {{}};

function makeSVG(width, height) {{
  const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', `0 0 ${{width}} ${{height}}`);
  svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
  return svg;
}}

function drawBarChart(containerId, bins, counts, color, xLabel) {{
  const el = document.getElementById(containerId);
  if (!el || !bins || !counts) return;
  const W = 700, H = 220, pad = {{l:50,r:20,t:10,b:40}};
  const svg = makeSVG(W, H);
  const maxVal = Math.max(...counts, 1);
  const n = counts.length;
  const bw = (W - pad.l - pad.r) / n;
  for (let i = 0; i < n; i++) {{
    const h = (counts[i] / maxVal) * (H - pad.t - pad.b);
    const x = pad.l + i * bw;
    const y = H - pad.b - h;
    const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    rect.setAttribute('x', x); rect.setAttribute('y', y);
    rect.setAttribute('width', Math.max(bw - 1, 1)); rect.setAttribute('height', h);
    rect.setAttribute('fill', color);
    svg.appendChild(rect);
  }}
  // X-axis labels (every 5th bin)
  for (let i = 0; i < bins.length; i += Math.max(1, Math.floor(n/10))) {{
    const txt = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    txt.setAttribute('x', pad.l + i * bw); txt.setAttribute('y', H - 5);
    txt.setAttribute('class', 'axis-label'); txt.setAttribute('text-anchor', 'middle');
    txt.textContent = bins[i] !== undefined ? (typeof bins[i] === 'number' ? bins[i].toFixed(bins[i] < 1 ? 2 : 0) : bins[i]) : '';
    svg.appendChild(txt);
  }}
  el.appendChild(svg);
}}

function drawLineChart(containerId, positions, yKey, color, yLabel) {{
  const el = document.getElementById(containerId);
  if (!el || !positions || positions.length === 0) return;
  const W = 700, H = 250, pad = {{l:50,r:20,t:20,b:40}};
  const svg = makeSVG(W, H);
  const n = positions.length;
  const xScale = (W - pad.l - pad.r) / n;
  const maxY = Math.max(...positions.map(p => p[yKey] || 0), 1);
  const minY = 0;
  const yRange = maxY - minY || 1;

  // Draw line
  let pathD = '';
  for (let i = 0; i < n; i++) {{
    const x = pad.l + i * xScale;
    const y = pad.t + (1 - (positions[i][yKey] - minY) / yRange) * (H - pad.t - pad.b);
    pathD += (i === 0 ? 'M' : 'L') + x.toFixed(1) + ',' + y.toFixed(1);
  }}
  const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
  path.setAttribute('d', pathD);
  path.setAttribute('fill', 'none'); path.setAttribute('stroke', color);
  path.setAttribute('stroke-width', '1.5');
  svg.appendChild(path);

  // Q25/Q75 band if available
  if (positions[0].q25 !== undefined) {{
    let bandD = '';
    for (let i = 0; i < n; i++) {{
      const x = pad.l + i * xScale;
      const y = pad.t + (1 - (positions[i].q75 - minY) / yRange) * (H - pad.t - pad.b);
      bandD += (i === 0 ? 'M' : 'L') + x.toFixed(1) + ',' + y.toFixed(1);
    }}
    for (let i = n - 1; i >= 0; i--) {{
      const x = pad.l + i * xScale;
      const y = pad.t + (1 - (positions[i].q25 - minY) / yRange) * (H - pad.t - pad.b);
      bandD += 'L' + x.toFixed(1) + ',' + y.toFixed(1);
    }}
    bandD += 'Z';
    const band = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    band.setAttribute('d', bandD);
    band.setAttribute('fill', color); band.setAttribute('opacity', '0.15');
    svg.appendChild(band);
  }}

  // Axes
  for (let v = 0; v <= 40; v += 10) {{
    const y = pad.t + (1 - v / yRange) * (H - pad.t - pad.b);
    const txt = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    txt.setAttribute('x', pad.l - 8); txt.setAttribute('y', y + 4);
    txt.setAttribute('class', 'axis-label'); txt.setAttribute('text-anchor', 'end');
    txt.textContent = 'Q' + v;
    svg.appendChild(txt);
  }}
  for (let i = 0; i < n; i += Math.max(1, Math.floor(n / 10))) {{
    const txt = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    txt.setAttribute('x', pad.l + i * xScale); txt.setAttribute('y', H - 5);
    txt.setAttribute('class', 'axis-label'); txt.setAttribute('text-anchor', 'middle');
    txt.textContent = positions[i].position;
    svg.appendChild(txt);
  }}
  el.appendChild(svg);
}}

function drawStackedBases(containerId, positions) {{
  const el = document.getElementById(containerId);
  if (!el || !positions || positions.length === 0) return;
  const W = 700, H = 200, pad = {{l:50,r:20,t:10,b:40}};
  const svg = makeSVG(W, H);
  const n = positions.length;
  const bw = (W - pad.l - pad.r) / n;
  const colors = {{a: '#22c55e', c: '#3b82f6', g: '#f59e0b', t: '#ef4444', n: '#94a3b8'}};
  const chartH = H - pad.t - pad.b;

  for (let i = 0; i < n; i++) {{
    const p = positions[i];
    const x = pad.l + i * bw;
    let y = pad.t;
    for (const [base, color] of [['a', colors.a], ['c', colors.c], ['g', colors.g], ['t', colors.t], ['n', colors.n]]) {{
      const h = (p[base] || 0) * chartH;
      if (h > 0.5) {{
        const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        rect.setAttribute('x', x); rect.setAttribute('y', y);
        rect.setAttribute('width', Math.max(bw, 1)); rect.setAttribute('height', h);
        rect.setAttribute('fill', color);
        svg.appendChild(rect);
      }}
      y += h;
    }}
  }}
  // Legend
  let lx = pad.l;
  for (const [label, color] of [['A', colors.a], ['C', colors.c], ['G', colors.g], ['T', colors.t], ['N', colors.n]]) {{
    const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    rect.setAttribute('x', lx); rect.setAttribute('y', H - 15);
    rect.setAttribute('width', 12); rect.setAttribute('height', 12);
    rect.setAttribute('fill', color);
    svg.appendChild(rect);
    const txt = document.createElementNS('http://www.w3.org/2000/svg', 'text');
    txt.setAttribute('x', lx + 16); txt.setAttribute('y', H - 5);
    txt.setAttribute('class', 'axis-label');
    txt.textContent = label;
    svg.appendChild(txt);
    lx += 40;
  }}
  el.appendChild(svg);
}}

// Render charts
if (qa.per_position) {{
  const qb = qa.per_position.quality_before;
  if (qb && qb.positions) drawLineChart('chart-quality', qb.positions, 'mean', '#3b82f6');
  const bb = qa.per_position.bases_before;
  if (bb && bb.positions) drawStackedBases('chart-bases', bb.positions);
}}
if (qa.distributions) {{
  const ld = qa.distributions.length_before;
  if (ld) drawBarChart('chart-length', ld.bin_edges, ld.counts, '#3b82f6');
  const gc = qa.distributions.gc_content;
  if (gc) drawBarChart('chart-gc', gc.bin_edges, gc.counts, '#22c55e');
}}
</script>
</body>
</html>"##,
        profile = passport.profile,
        version = passport.tool_version,
        tier_color = tier_color,
        tier = format!("{:?}", passport.quality_tier).to_uppercase(),
        reads_input = format_number(passport.reads_input),
        reads_passed = format_number(passport.reads_passed),
        survival_pct = survival_pct,
        paired_html = paired_html,
        flags_html = flags_html,
        funnel_rows = funnel_rows,
        contaminant_detail = contaminant_detail,
        passport_json = passport_json,
    )
}

fn format_number(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{:.1}K", n as f64 / 1_000.0)
    } else {
        format!("{n}")
    }
}
