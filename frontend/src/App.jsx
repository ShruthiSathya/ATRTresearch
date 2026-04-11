/**
 * App.jsx — ATRT Drug Repurposing Frontend v2.0
 *
 * Field name contract (matches _normalise_candidate() in save_results.py):
 *   candidate.name                    — drug name
 *   candidate.mechanism               — mechanism of action
 *   candidate.score                   — composite score (0–1)
 *   candidate.tissue_expression_score — GSE70678 tissue score (0–1)
 *   candidate.depmap_score            — DepMap CRISPR essentiality (0–1)
 *   candidate.escape_bypass_score     — resistance bypass score (0–1)
 *   candidate.ppi_score               — STRING-DB PPI proximity (0–1)
 *   candidate.bbb_penetrance          — "HIGH" | "MODERATE" | "LOW" | "UNKNOWN"
 *   candidate.ezh2_boosted            — bool: EZH2 inhibitor (×1.40 boost applied)
 *   candidate.aurka_boosted           — bool: AURKA inhibitor
 *   candidate.ic50_validated          — bool
 *   candidate.ic50_um                 — float µM or null
 *   candidate.ic50_cell_line          — string e.g. "BT16"
 *   candidate.has_generic             — bool
 *   candidate.cmap_score              — CMap reversal score (0–1, 0.50 = neutral prior)
 *
 * API endpoint contract (matches main.py /analyze):
 *   POST /analyze → { success, candidates[], filtered_count, filtered_drugs[], stats }
 *   POST /generate_ai_analysis → { success, analysis }
 */

import { useState, useEffect, useRef } from 'react';

// Vite proxies /api → http://localhost:8000 (see vite.config.js)
const API_BASE = '/api';

/* ─── Field accessors matching backend schema ─────────────────────────────── */
const getDrugName = (c) => c.name || c.drug_name || 'Unknown';
const getScore    = (c) => c.score ?? c.composite_score ?? 0;
const getBBB      = (c) => c.bbb_penetrance || c.bbb_status || 'UNKNOWN';

const BBB_COLOR = {
  HIGH:     { bg: '#e4f0e8', color: '#2a7a4b', label: 'BBB HIGH' },
  MODERATE: { bg: '#f0ece4', color: '#7a6e2a', label: 'BBB MOD'  },
  LOW:      { bg: '#fce8e4', color: '#c8502a', label: 'BBB LOW'   },
  UNKNOWN:  { bg: '#f0f0f0', color: '#9c9894', label: 'BBB ?'     },
};

const fmtPct = (n) => n != null ? `${Math.round((n ?? 0) * 100)}%` : '—';
const fmtScore = (n) => n != null ? n.toFixed(3) : '—';
const fmtIC50  = (n) => n != null ? `${n} µM` : '—';

/* ─── Styles ──────────────────────────────────────────────────────────────── */
const STYLES = `
@import url('https://fonts.googleapis.com/css2?family=Fraunces:ital,wght@0,300;0,600;0,900;1,300;1,600&family=JetBrains+Mono:wght@300;400;500&family=Instrument+Sans:wght@300;400;500&display=swap');

*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }

:root {
  --ink:       #0d0c0b;
  --ink-soft:  #46423e;
  --ink-faint: #9c9894;
  --paper:     #f6f3ee;
  --card:      #ffffff;
  --accent:    #bf4b2a;
  --accent-lt: #f2ebe6;
  --green:     #28744a;
  --green-lt:  #e5f0ea;
  --amber:     #7a6520;
  --amber-lt:  #f2ece0;
  --rule:      #ddd9d3;
  --radius:    2px;
  --ff-display: 'Fraunces', Georgia, serif;
  --ff-body:    'Instrument Sans', sans-serif;
  --ff-mono:    'JetBrains Mono', monospace;
  --shadow:     0 1px 3px rgba(13,12,11,.06), 0 6px 24px rgba(13,12,11,.08);
  --shadow-lg:  0 2px 8px rgba(13,12,11,.10), 0 16px 48px rgba(13,12,11,.12);
}

html { scroll-behavior: smooth; }

body {
  background: var(--paper);
  color: var(--ink);
  font-family: var(--ff-body);
  font-weight: 300;
  line-height: 1.6;
  min-height: 100vh;
}

/* ── Layout ── */
.app { max-width: 820px; margin: 0 auto; padding: 72px 24px 120px; }

/* ── Header ── */
.header { margin-bottom: 60px; }

.eyebrow {
  font-family: var(--ff-mono);
  font-size: 10px;
  letter-spacing: 0.20em;
  text-transform: uppercase;
  color: var(--ink-faint);
  display: flex;
  align-items: center;
  gap: 8px;
  margin-bottom: 20px;
}
.status-orb {
  width: 7px; height: 7px;
  border-radius: 50%;
  flex-shrink: 0;
  transition: background 0.4s;
}

.headline {
  font-family: var(--ff-display);
  font-size: clamp(44px, 8vw, 72px);
  font-weight: 900;
  line-height: 1.03;
  color: var(--ink);
  margin-bottom: 16px;
}
.headline em { color: var(--accent); font-style: italic; font-weight: 300; }

.sub {
  font-size: 15px;
  color: var(--ink-soft);
  max-width: 520px;
  font-weight: 300;
  line-height: 1.7;
}

/* ── Search ── */
.search-section { margin-bottom: 12px; }
.search-label {
  display: block;
  font-family: var(--ff-mono);
  font-size: 10px;
  letter-spacing: 0.14em;
  text-transform: uppercase;
  color: var(--ink-faint);
  margin-bottom: 10px;
}
.search-row { display: flex; gap: 10px; }
.search-input {
  flex: 1;
  background: var(--card);
  border: 1.5px solid var(--rule);
  border-radius: var(--radius);
  padding: 15px 20px;
  font-family: var(--ff-body);
  font-size: 16px;
  font-weight: 300;
  color: var(--ink);
  outline: none;
  transition: border-color 0.2s;
}
.search-input::placeholder { color: var(--ink-faint); }
.search-input:focus { border-color: var(--accent); }
.search-btn {
  background: var(--ink);
  color: var(--paper);
  border: none;
  border-radius: var(--radius);
  padding: 15px 32px;
  font-family: var(--ff-mono);
  font-size: 11px;
  font-weight: 500;
  letter-spacing: 0.10em;
  text-transform: uppercase;
  cursor: pointer;
  white-space: nowrap;
  transition: background 0.18s, transform 0.1s;
}
.search-btn:hover:not(:disabled) { background: var(--accent); }
.search-btn:active:not(:disabled) { transform: scale(0.97); }
.search-btn:disabled { opacity: 0.4; cursor: not-allowed; }

/* ── Filters ── */
.filters {
  display: flex;
  gap: 24px;
  flex-wrap: wrap;
  padding: 16px 0 0;
  border-top: 1px solid var(--rule);
  margin-bottom: 52px;
}
.filter-group { display: flex; flex-direction: column; gap: 7px; }
.filter-label {
  font-family: var(--ff-mono);
  font-size: 9px;
  letter-spacing: 0.15em;
  text-transform: uppercase;
  color: var(--ink-faint);
}
.pill-row { display: flex; gap: 5px; flex-wrap: wrap; }
.pill {
  font-family: var(--ff-mono);
  font-size: 10px;
  padding: 5px 11px;
  border: 1px solid var(--rule);
  border-radius: 20px;
  background: var(--card);
  color: var(--ink-soft);
  cursor: pointer;
  letter-spacing: 0.04em;
  transition: all 0.15s;
  user-select: none;
}
.pill:hover { border-color: var(--accent); color: var(--accent); }
.pill-on { background: var(--accent-lt); border-color: var(--accent); color: var(--accent); }
.filter-select {
  background: var(--card);
  border: 1px solid var(--rule);
  border-radius: var(--radius);
  padding: 5px 10px;
  font-family: var(--ff-mono);
  font-size: 10px;
  color: var(--ink-soft);
  cursor: pointer;
  outline: none;
}
.filter-select:focus { border-color: var(--accent); }

/* ── Error ── */
.error-bar {
  margin-bottom: 28px;
  padding: 12px 16px;
  background: #fef2f2;
  border-left: 3px solid #dc2626;
  font-family: var(--ff-mono);
  font-size: 12px;
  color: #dc2626;
  border-radius: 0 var(--radius) var(--radius) 0;
}

/* ── Results header ── */
.results-header {
  border-top: 2px solid var(--ink);
  padding-top: 28px;
  margin-bottom: 32px;
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  gap: 12px;
  flex-wrap: wrap;
}
.results-title {
  font-family: var(--ff-display);
  font-size: 30px;
  font-weight: 600;
  color: var(--ink);
  text-transform: capitalize;
}
.results-meta {
  font-family: var(--ff-mono);
  font-size: 10px;
  color: var(--ink-faint);
  letter-spacing: 0.10em;
}

/* ── Stats bar ── */
.stats-bar {
  display: flex;
  gap: 0;
  border: 1.5px solid var(--rule);
  border-radius: var(--radius);
  overflow: hidden;
  margin-bottom: 28px;
  background: var(--card);
}
.stat-cell {
  flex: 1;
  padding: 12px 16px;
  border-right: 1px solid var(--rule);
  display: flex;
  flex-direction: column;
  gap: 3px;
}
.stat-cell:last-child { border-right: none; }
.stat-cell-label {
  font-family: var(--ff-mono);
  font-size: 9px;
  letter-spacing: 0.12em;
  text-transform: uppercase;
  color: var(--ink-faint);
}
.stat-cell-value {
  font-family: var(--ff-mono);
  font-size: 18px;
  font-weight: 500;
  color: var(--ink);
}

/* ── Drug list ── */
.drug-list { display: flex; flex-direction: column; gap: 3px; }

/* ── Drug row ── */
.drug-row {
  background: var(--card);
  border: 1.5px solid var(--rule);
  border-radius: var(--radius);
  overflow: hidden;
  transition: box-shadow 0.18s, border-color 0.18s;
  cursor: pointer;
}
.drug-row:hover { box-shadow: var(--shadow); border-color: #ccc8c2; }

.drug-main {
  display: flex;
  align-items: flex-start;
  padding: 18px 22px;
  gap: 14px;
}
.drug-rank {
  font-family: var(--ff-mono);
  font-size: 10px;
  color: var(--ink-faint);
  min-width: 24px;
  padding-top: 4px;
  flex-shrink: 0;
}
.drug-info { flex: 1; min-width: 0; }
.drug-name {
  font-family: var(--ff-display);
  font-size: 20px;
  font-weight: 600;
  color: var(--ink);
  margin-bottom: 3px;
}
.drug-mech {
  font-size: 12px;
  color: var(--ink-faint);
  font-weight: 300;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  max-width: 400px;
}
.drug-tags {
  display: flex;
  gap: 5px;
  flex-shrink: 0;
  flex-wrap: wrap;
  padding-top: 3px;
}
.tag {
  font-family: var(--ff-mono);
  font-size: 9px;
  letter-spacing: 0.07em;
  padding: 3px 8px;
  border-radius: 2px;
  font-weight: 500;
  white-space: nowrap;
}
.tag-ezh2   { background: #fce8e4; color: #bf4b2a; }
.tag-aurka  { background: #ebe4f5; color: #5a3a9a; }
.tag-generic { background: var(--green-lt); color: var(--green); }

.score-cluster {
  display: flex;
  gap: 18px;
  flex-shrink: 0;
  padding-top: 2px;
}
.mini-meter { display: flex; flex-direction: column; gap: 5px; min-width: 72px; }
.mini-meter-labels { display: flex; justify-content: space-between; }
.mini-meter-name {
  font-family: var(--ff-mono);
  font-size: 8px;
  letter-spacing: 0.12em;
  text-transform: uppercase;
  color: var(--ink-faint);
}
.mini-meter-val {
  font-family: var(--ff-mono);
  font-size: 12px;
  font-weight: 500;
  color: var(--ink);
}
.bar-track { width: 72px; height: 2px; background: var(--rule); border-radius: 2px; overflow: hidden; }
.bar-fill  { height: 100%; border-radius: 2px; transition: width 0.5s cubic-bezier(.16,1,.3,1); }
.bar-composite { background: var(--ink); }
.bar-depmap    { background: var(--accent); }

.chevron {
  font-size: 10px;
  color: var(--ink-faint);
  flex-shrink: 0;
  margin-left: 4px;
  padding-top: 6px;
  transition: transform 0.2s;
  align-self: flex-start;
}
.chevron.open { transform: rotate(180deg); }

/* ── Expanded detail ── */
.drug-detail {
  border-top: 1px solid var(--rule);
  padding: 22px 22px 22px 60px;
  background: #faf9f7;
  display: flex;
  flex-direction: column;
  gap: 18px;
  animation: slideIn 0.15s ease;
}

@keyframes slideIn {
  from { opacity: 0; transform: translateY(-5px); }
  to   { opacity: 1; transform: translateY(0); }
}

.detail-scores-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(110px, 1fr));
  gap: 10px;
}
.score-cell {
  background: var(--card);
  border: 1px solid var(--rule);
  border-radius: var(--radius);
  padding: 10px 12px;
}
.score-cell-label {
  font-family: var(--ff-mono);
  font-size: 8px;
  letter-spacing: 0.14em;
  text-transform: uppercase;
  color: var(--ink-faint);
  margin-bottom: 4px;
}
.score-cell-value {
  font-family: var(--ff-mono);
  font-size: 20px;
  font-weight: 500;
  color: var(--ink);
}
.score-cell-sub {
  font-family: var(--ff-mono);
  font-size: 9px;
  color: var(--ink-faint);
  margin-top: 2px;
}

.ic50-block {
  background: var(--card);
  border: 1px solid var(--rule);
  border-radius: var(--radius);
  padding: 10px 14px;
  display: flex;
  gap: 20px;
  flex-wrap: wrap;
  align-items: center;
}
.ic50-label {
  font-family: var(--ff-mono);
  font-size: 9px;
  letter-spacing: 0.12em;
  text-transform: uppercase;
  color: var(--ink-faint);
  margin-bottom: 3px;
}
.ic50-value { font-family: var(--ff-mono); font-size: 15px; font-weight: 500; color: var(--green); }
.ic50-source { font-family: var(--ff-mono); font-size: 9px; color: var(--ink-faint); }

.rationale-block {
  border-left: 2px solid #b0ccb8;
  padding-left: 14px;
}
.rationale-label {
  font-family: var(--ff-mono);
  font-size: 9px;
  letter-spacing: 0.14em;
  text-transform: uppercase;
  color: var(--green);
  margin-bottom: 6px;
}
.rationale-text { font-size: 13px; color: var(--ink-soft); line-height: 1.7; font-weight: 300; }
.rationale-spinner { display: flex; align-items: center; gap: 8px; }
.r-dot {
  width: 4px; height: 4px;
  background: var(--green);
  border-radius: 50%;
  animation: rbounce 1.2s ease-in-out infinite;
  opacity: 0.3;
}
.r-dot:nth-child(2) { animation-delay: 0.2s; }
.r-dot:nth-child(3) { animation-delay: 0.4s; }
@keyframes rbounce {
  0%, 80%, 100% { opacity: 0.2; transform: scale(0.8); }
  40% { opacity: 0.9; transform: scale(1); }
}
.rationale-spinner-text {
  font-family: var(--ff-mono);
  font-size: 10px;
  color: var(--ink-faint);
  letter-spacing: 0.1em;
}

/* ── Filtered drugs ── */
.filtered-block {
  margin-top: 36px;
  padding-top: 20px;
  border-top: 1px solid var(--rule);
}
.filtered-summary {
  font-family: var(--ff-mono);
  font-size: 10px;
  color: var(--ink-faint);
  letter-spacing: 0.12em;
  cursor: pointer;
  text-transform: uppercase;
  list-style: none;
}
.filtered-table { width: 100%; border-collapse: collapse; margin-top: 12px; }
.filtered-table th {
  font-family: var(--ff-mono);
  font-size: 9px;
  letter-spacing: 0.12em;
  text-transform: uppercase;
  color: var(--ink-faint);
  text-align: left;
  padding: 6px 8px;
  border-bottom: 1px solid var(--rule);
}
.filtered-table td {
  font-family: var(--ff-mono);
  font-size: 11px;
  color: var(--ink-soft);
  padding: 7px 8px;
  border-bottom: 1px solid var(--rule);
}

/* ── Loading ── */
.loading {
  padding: 52px 0;
  display: flex;
  align-items: center;
  gap: 14px;
}
.dots { display: flex; gap: 5px; }
.dot {
  width: 5px; height: 5px;
  background: var(--ink);
  border-radius: 50%;
  animation: rbounce 1.2s ease-in-out infinite;
}
.dot:nth-child(2) { animation-delay: 0.2s; }
.dot:nth-child(3) { animation-delay: 0.4s; }
.loading-text {
  font-family: var(--ff-mono);
  font-size: 11px;
  letter-spacing: 0.14em;
  text-transform: uppercase;
  color: var(--ink-faint);
}

/* ── Footer ── */
.footer {
  margin-top: 80px;
  padding-top: 24px;
  border-top: 1px solid var(--rule);
  font-family: var(--ff-mono);
  font-size: 9px;
  color: var(--ink-faint);
  letter-spacing: 0.10em;
  line-height: 2;
}

/* ── Responsive ── */
@media (max-width: 620px) {
  .score-cluster { display: none; }
  .drug-main { padding: 14px 16px; gap: 10px; }
  .drug-detail { padding-left: 16px; }
  .search-row { flex-direction: column; }
}
`;

/* ─── AI rationale fetch ───────────────────────────────────────────────────── */
async function fetchRationale(candidate, disease) {
  const res = await fetch(`${API_BASE}/generate_ai_analysis`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      target: getDrugName(candidate),
      disease,
      mechanism: candidate.mechanism || '',
    }),
  });
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  const data = await res.json();
  if (!data.success) throw new Error(data.error || 'AI analysis failed');
  return data.analysis || 'Rationale unavailable.';
}

/* ─── DrugRow component ───────────────────────────────────────────────────── */
function DrugRow({ candidate, rank, disease }) {
  const [open, setOpen]                   = useState(false);
  const [rationale, setRationale]         = useState(null);
  const [rationaleLoading, setRLoading]   = useState(false);
  const rationaleRef                      = useRef(null);

  const name    = getDrugName(candidate);
  const score   = getScore(candidate);
  const bbb     = getBBB(candidate);
  const bbbStyle = BBB_COLOR[bbb] || BBB_COLOR.UNKNOWN;

  const handleToggle = async () => {
    const next = !open;
    setOpen(next);
    if (next && !rationale && !rationaleLoading) {
      setRLoading(true);
      try {
        const text = await fetchRationale(candidate, disease);
        setRationale(text);
      } catch {
        setRationale('Could not generate rationale. Ensure the backend AI endpoint is configured.');
      } finally {
        setRLoading(false);
      }
    }
  };

  return (
    <div className="drug-row">
      <div className="drug-main" onClick={handleToggle}>
        <span className="drug-rank">{String(rank).padStart(2, '0')}</span>

        <div className="drug-info">
          <div className="drug-name">{name}</div>
          {candidate.mechanism && (
            <div className="drug-mech">{candidate.mechanism}</div>
          )}
        </div>

        {/* Tags */}
        <div className="drug-tags">
          <span className="tag" style={{ background: bbbStyle.bg, color: bbbStyle.color }}>
            {bbbStyle.label}
          </span>
          {candidate.ezh2_boosted && (
            <span className="tag tag-ezh2">EZH2↑ ×1.40</span>
          )}
          {candidate.aurka_boosted && (
            <span className="tag tag-aurka">AURKA↑</span>
          )}
          {candidate.has_generic && (
            <span className="tag tag-generic">GENERIC</span>
          )}
        </div>

        {/* Mini score meters */}
        <div className="score-cluster">
          <div className="mini-meter">
            <div className="mini-meter-labels">
              <span className="mini-meter-name">Composite</span>
              <span className="mini-meter-val">{fmtPct(score)}</span>
            </div>
            <div className="bar-track">
              <div className="bar-fill bar-composite" style={{ width: fmtPct(score) }} />
            </div>
          </div>
          <div className="mini-meter">
            <div className="mini-meter-labels">
              <span className="mini-meter-name">DepMap</span>
              <span className="mini-meter-val">{fmtPct(candidate.depmap_score)}</span>
            </div>
            <div className="bar-track">
              <div className="bar-fill bar-depmap" style={{ width: fmtPct(candidate.depmap_score) }} />
            </div>
          </div>
        </div>

        <span className={`chevron ${open ? 'open' : ''}`}>▾</span>
      </div>

      {open && (
        <div className="drug-detail" ref={rationaleRef}>
          {/* Score breakdown grid */}
          <div className="detail-scores-grid">
            {[
              {
                label: 'Composite',
                value: fmtScore(score),
                sub: 'tissue×0.40 + depmap×0.35 + escape×0.20 + ppi×0.05',
              },
              {
                label: 'Tissue / GSE70678',
                value: fmtScore(candidate.tissue_expression_score),
                sub: '49 ATRT samples (Torchia 2015)',
              },
              {
                label: 'DepMap CRISPR',
                value: fmtScore(candidate.depmap_score),
                sub: 'BT16/BT37/G401/A204 (Behan 2019)',
              },
              {
                label: 'Escape Bypass',
                value: fmtScore(candidate.escape_bypass_score),
                sub: 'SMARCB1-loss resistance nodes',
              },
              {
                label: 'PPI Network',
                value: fmtScore(candidate.ppi_score),
                sub: 'STRING-DB v11.5 proximity',
              },
              {
                label: 'CMap L1000',
                value: fmtScore(candidate.cmap_score),
                sub: candidate.cmap_score != null && Math.abs((candidate.cmap_score ?? 0.5) - 0.5) < 0.01
                  ? 'Neutral prior (not in L1000)'
                  : 'SMARCB1-loss reversal score',
              },
            ].map(({ label, value, sub }) => (
              <div className="score-cell" key={label}>
                <div className="score-cell-label">{label}</div>
                <div className="score-cell-value">{value}</div>
                <div className="score-cell-sub">{sub}</div>
              </div>
            ))}
          </div>

          {/* IC50 validation (only if verified primary literature data exists) */}
          {candidate.ic50_validated && candidate.ic50_um != null && (
            <div className="ic50-block">
              <div>
                <div className="ic50-label">Published IC₅₀ (ATRT cell lines)</div>
                <div className="ic50-value">{fmtIC50(candidate.ic50_um)}</div>
              </div>
              {candidate.ic50_cell_line && (
                <div>
                  <div className="ic50-label">Cell line</div>
                  <div className="ic50-source">{candidate.ic50_cell_line}</div>
                </div>
              )}
              {candidate.ic50_source && (
                <div style={{ flex: 1 }}>
                  <div className="ic50-label">Source</div>
                  <div className="ic50-source">{candidate.ic50_source}</div>
                </div>
              )}
            </div>
          )}

          {/* AI rationale */}
          <div className="rationale-block">
            <div className="rationale-label">AI Biological Rationale (SMARCB1-null context)</div>
            {rationaleLoading ? (
              <div className="rationale-spinner">
                <div className="r-dot"/><div className="r-dot"/><div className="r-dot"/>
                <span className="rationale-spinner-text">Generating rationale...</span>
              </div>
            ) : (
              <p className="rationale-text">
                {rationale || 'Click to expand — rationale will load automatically.'}
              </p>
            )}
          </div>
        </div>
      )}
    </div>
  );
}

/* ─── Stats bar ───────────────────────────────────────────────────────────── */
function StatsBar({ stats }) {
  if (!stats) return null;
  const smarcb1Pct = stats.total_samples
    ? Math.round((stats.smarcb1_loss_count / stats.total_samples) * 100)
    : null;

  const cells = [
    { label: 'Screened', value: stats.n_screened ?? stats.n_drugs_screened ?? '—' },
    {
      label: 'SMARCB1 loss',
      value: smarcb1Pct != null
        ? `${smarcb1Pct}%`
        : `${stats.smarcb1_loss_count ?? 0}`,
      sub: smarcb1Pct != null ? `of ${stats.total_samples} samples` : '~95% expected',
    },
    { label: 'EZH2 boosted', value: stats.n_ezh2_boosted ?? '—' },
    { label: 'DepMap source', value: stats.depmap_source?.includes('CSV') ? 'Live 24Q4' : 'Fallback' },
  ].filter(c => c.value !== '—' && c.value !== undefined);

  return (
    <div className="stats-bar">
      {cells.slice(0, 4).map(({ label, value, sub }) => (
        <div className="stat-cell" key={label}>
          <span className="stat-cell-label">{label}</span>
          <span className="stat-cell-value">{value}</span>
          {sub && <span style={{ fontFamily: 'var(--ff-mono)', fontSize: '9px', color: 'var(--ink-faint)' }}>{sub}</span>}
        </div>
      ))}
    </div>
  );
}

/* ─── Main App ────────────────────────────────────────────────────────────── */
const SUBGROUPS = ['Pan', 'TYR', 'SHH', 'MYC'];
const LOCATIONS = [
  { value: 'unknown_location', label: 'Location unknown' },
  { value: 'infratentorial',   label: 'Infratentorial (~50%)' },
  { value: 'supratentorial',   label: 'Supratentorial (~35%)' },
];

export default function App() {
  const [query,    setQuery]    = useState('atrt');
  const [subgroup, setSubgroup] = useState('Pan');
  const [location, setLocation] = useState('unknown_location');
  const [loading,  setLoading]  = useState(false);
  const [results,  setResults]  = useState(null);
  const [error,    setError]    = useState(null);
  const [backendOk, setBackend] = useState(null);

  useEffect(() => {
    fetch(`${API_BASE}/`)
      .then(r => r.ok ? r.json() : Promise.reject(r.status))
      .then(d => setBackend(d.pipeline_ready ?? true))
      .catch(() => setBackend(false));
  }, []);

  const handleSubmit = async (e) => {
    e.preventDefault();
    const q = query.trim();
    if (!q) return;
    setLoading(true); setError(null); setResults(null);

    try {
      const res = await fetch(`${API_BASE}/analyze`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          disease_name: q.toLowerCase(),
          min_score:   0.20,
          max_results: 15,
          subgroup: subgroup === 'Pan' ? undefined : subgroup,
          location,
        }),
      });
      if (!res.ok) throw new Error(`HTTP ${res.status} from /analyze`);
      const data = await res.json();
      if (!data.success) throw new Error(data.error || 'Pipeline returned an error.');
      setResults(data);
    } catch (err) {
      setError(
        err.message.includes('fetch')
          ? `Cannot reach backend at ${API_BASE}. Is uvicorn running on :8000?`
          : err.message
      );
    } finally {
      setLoading(false);
    }
  };

  const candidates = results?.candidates ?? [];

  return (
    <>
      <style>{STYLES}</style>
      <div className="app">

        {/* Header */}
        <header className="header">
          <p className="eyebrow">
            <span
              className="status-orb"
              style={{
                background:
                  backendOk === null ? '#ccc' :
                  backendOk           ? '#28744a' : '#bf4b2a',
              }}
            />
            ATRT Drug Repurposing Engine · GSE70678 · DepMap 24Q4
          </p>
          <h1 className="headline">
            Multi-omic<br /><em>drug repurposing</em>
          </h1>
          <p className="sub">
            SMARCB1-deficient ATRT pipeline. Ranks candidates by tissue expression (GSE70678),
            CRISPR essentiality (DepMap, rhabdoid lines), resistance bypass, and PPI proximity.
            EZH2 inhibitors carry a synthetic-lethality boost ×1.40 (Knutson 2013 PNAS).
          </p>
        </header>

        {/* Search */}
        <div className="search-section">
          <label className="search-label" htmlFor="disease-q">Disease or indication</label>
          <form className="search-row" onSubmit={handleSubmit}>
            <input
              id="disease-q"
              className="search-input"
              type="text"
              value={query}
              onChange={e => setQuery(e.target.value)}
              placeholder="e.g. atrt, rhabdoid tumor, SMARCB1-deficient"
              autoComplete="off"
              spellCheck="false"
            />
            <button
              type="submit"
              className="search-btn"
              disabled={loading || !query.trim()}
            >
              {loading ? 'Running…' : 'Analyse →'}
            </button>
          </form>
        </div>

        {/* Filters */}
        <div className="filters">
          <div className="filter-group">
            <span className="filter-label">Molecular subgroup</span>
            <div className="pill-row">
              {SUBGROUPS.map(s => (
                <button
                  key={s}
                  type="button"
                  className={`pill ${subgroup === s ? 'pill-on' : ''}`}
                  onClick={() => setSubgroup(s)}
                >
                  {s === 'Pan' ? 'Pan-ATRT' : `ATRT-${s}`}
                </button>
              ))}
            </div>
          </div>

          <div className="filter-group">
            <span className="filter-label">Tumour location</span>
            <select
              className="filter-select"
              value={location}
              onChange={e => setLocation(e.target.value)}
            >
              {LOCATIONS.map(({ value, label }) => (
                <option key={value} value={value}>{label}</option>
              ))}
            </select>
          </div>
        </div>

        {error && <div className="error-bar">⚠ {error}</div>}

        {loading && (
          <div className="loading">
            <div className="dots">
              <div className="dot"/><div className="dot"/><div className="dot"/>
            </div>
            <span className="loading-text">Scoring compounds · DepMap · Tissue · Escape bypass</span>
          </div>
        )}

        {results && !loading && (
          <>
            <div className="results-header">
              <h2 className="results-title">{query}</h2>
              <span className="results-meta">
                {candidates.length} candidates
                {results.filtered_count > 0 && ` · ${results.filtered_count} safety-filtered`}
                {subgroup !== 'Pan' && ` · subgroup: ATRT-${subgroup}`}
              </span>
            </div>

            <StatsBar stats={results.stats} />

            <div className="drug-list">
              {candidates.length === 0 ? (
                <p style={{ fontFamily: 'var(--ff-mono)', fontSize: '12px', color: 'var(--ink-faint)', padding: '20px 0' }}>
                  No candidates returned above the minimum score threshold. Try lowering the score or checking backend logs.
                </p>
              ) : (
                candidates.map((c, i) => (
                  <DrugRow
                    key={getDrugName(c) + i}
                    candidate={c}
                    rank={i + 1}
                    disease={query}
                  />
                ))
              )}
            </div>

            {/* Filtered drugs */}
            {results.filtered_drugs?.length > 0 && (
              <div className="filtered-block">
                <details>
                  <summary className="filtered-summary">
                    {results.filtered_drugs.length} drug(s) removed by safety filter ▾
                  </summary>
                  <table className="filtered-table">
                    <thead>
                      <tr><th>Drug</th><th>Reason</th><th>Severity</th></tr>
                    </thead>
                    <tbody>
                      {results.filtered_drugs.map((d, i) => (
                        <tr key={i}>
                          <td>{d.drug_name}</td>
                          <td>{d.reason}</td>
                          <td>{d.severity}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </details>
              </div>
            )}
          </>
        )}

        <footer className="footer">
          DATA · GSE70678 (Torchia 2015, Cancer Cell PMID 26609405) ·
          DepMap 24Q4 (BT16/BT37/G401/A204) · STRING-DB v11.5 ·
          CMap L1000 (Subramanian 2017, Cell PMID 29195078)
          <br />
          SCORING · composite = tissue×0.40 + depmap×0.35 + escape×0.20 + ppi×0.05 ·
          EZH2 inhibitors ×1.40 (Knutson 2013, PNAS PMID 23620515) ·
          Confidence reported as [conservative, optimistic] (PBTC-047 precedent, Monje 2023)
        </footer>
      </div>
    </>
  );
}