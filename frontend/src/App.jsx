import { useState, useEffect } from 'react';

const API_BASE = './api';

/* ─── helpers ─────────────────────────────────────────────── */
const pct = (n) => `${Math.round((n ?? 0) * 100)}%`;

const noveltyScore = (c) => {
  const base = 1 - (c.score ?? 0);
  const bbbBoost = c.bbb_status === 'HIGH' ? 0.08 : 0;
  return Math.max(0, Math.min(1, base + bbbBoost));
};

const successRate = (c) => c.score ?? c.composite_score ?? 0;

/* ─── styles ──────────────────────────────────────────────── */
const styles = `
  @import url('https://fonts.googleapis.com/css2?family=DM+Serif+Display:ital@0;1&family=DM+Mono:wght@300;400;500&family=DM+Sans:wght@300;400;500&display=swap');

  *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }

  :root {
    --ink: #0f0e0d;
    --ink-soft: #4a4743;
    --ink-faint: #9c9894;
    --paper: #f7f5f2;
    --paper-warm: #ede9e2;
    --paper-card: #ffffff;
    --accent: #c8502a;
    --accent-light: #f0e8e4;
    --success: #2a7a4b;
    --success-light: #e4f0e8;
    --mid: #7a6e2a;
    --mid-light: #f0ece4;
    --rule: #e0dbd4;
    --shadow-hover: 0 2px 8px rgba(15,14,13,0.1), 0 8px 32px rgba(15,14,13,0.09);
    --radius: 3px;
    --font-display: 'DM Serif Display', Georgia, serif;
    --font-body: 'DM Sans', sans-serif;
    --font-mono: 'DM Mono', monospace;
  }

  body {
    background: var(--paper);
    color: var(--ink);
    font-family: var(--font-body);
    font-weight: 300;
    min-height: 100vh;
  }

  .app { max-width: 780px; margin: 0 auto; padding: 64px 24px 120px; }

  .header { margin-bottom: 56px; }

  .header-eyebrow {
    font-family: var(--font-mono);
    font-size: 11px;
    font-weight: 400;
    letter-spacing: 0.18em;
    text-transform: uppercase;
    color: var(--ink-faint);
    margin-bottom: 16px;
    display: flex;
    align-items: center;
    gap: 6px;
  }

  .status-dot {
    display: inline-block;
    width: 6px; height: 6px;
    border-radius: 50%;
    flex-shrink: 0;
  }

  .header-title {
    font-family: var(--font-display);
    font-size: clamp(40px, 7vw, 62px);
    line-height: 1.05;
    color: var(--ink);
    margin-bottom: 14px;
  }

  .header-title em { color: var(--accent); font-style: italic; }

  .header-sub {
    font-size: 15px;
    color: var(--ink-soft);
    font-weight: 300;
    line-height: 1.6;
    max-width: 480px;
  }

  .search-block { margin-bottom: 24px; }

  .search-label {
    display: block;
    font-family: var(--font-mono);
    font-size: 11px;
    letter-spacing: 0.15em;
    text-transform: uppercase;
    color: var(--ink-faint);
    margin-bottom: 10px;
  }

  .search-row { display: flex; gap: 12px; align-items: stretch; }

  .search-input {
    flex: 1;
    background: var(--paper-card);
    border: 1.5px solid var(--rule);
    border-radius: var(--radius);
    padding: 14px 18px;
    font-family: var(--font-body);
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
    padding: 14px 28px;
    font-family: var(--font-mono);
    font-size: 12px;
    font-weight: 500;
    letter-spacing: 0.1em;
    text-transform: uppercase;
    cursor: pointer;
    white-space: nowrap;
    transition: background 0.2s, transform 0.1s;
  }

  .search-btn:hover:not(:disabled) { background: var(--accent); }
  .search-btn:active:not(:disabled) { transform: scale(0.98); }
  .search-btn:disabled { opacity: 0.45; cursor: not-allowed; }

  .filters {
    display: flex;
    gap: 20px;
    flex-wrap: wrap;
    margin-bottom: 48px;
    padding-top: 14px;
    border-top: 1px solid var(--rule);
  }

  .filter-group { display: flex; flex-direction: column; gap: 6px; }

  .filter-label {
    font-family: var(--font-mono);
    font-size: 9px;
    letter-spacing: 0.14em;
    text-transform: uppercase;
    color: var(--ink-faint);
  }

  .pill-row { display: flex; gap: 4px; flex-wrap: wrap; }

  .pill {
    font-family: var(--font-mono);
    font-size: 10px;
    padding: 4px 10px;
    border: 1px solid var(--rule);
    border-radius: 20px;
    background: var(--paper-card);
    color: var(--ink-soft);
    cursor: pointer;
    transition: all 0.15s;
    letter-spacing: 0.04em;
  }

  .pill:hover { border-color: var(--accent); color: var(--accent); }
  .pill-active { background: var(--accent-light); border-color: var(--accent); color: var(--accent); }

  .filter-select {
    background: var(--paper-card);
    border: 1px solid var(--rule);
    border-radius: var(--radius);
    padding: 5px 10px;
    font-family: var(--font-mono);
    font-size: 11px;
    color: var(--ink-soft);
    outline: none;
    cursor: pointer;
  }

  .filter-select:focus { border-color: var(--accent); }

  .error-bar {
    margin-bottom: 32px;
    padding: 12px 16px;
    background: #fef2f2;
    border-left: 3px solid #dc2626;
    font-family: var(--font-mono);
    font-size: 13px;
    color: #dc2626;
    border-radius: 0 var(--radius) var(--radius) 0;
  }

  .results-header {
    border-top: 1.5px solid var(--ink);
    padding-top: 28px;
    margin-bottom: 36px;
    display: flex;
    align-items: baseline;
    justify-content: space-between;
    gap: 16px;
    flex-wrap: wrap;
  }

  .results-title { font-family: var(--font-display); font-size: 28px; color: var(--ink); text-transform: capitalize; }
  .results-meta { font-family: var(--font-mono); font-size: 11px; color: var(--ink-faint); letter-spacing: 0.1em; }

  .regimen-list { display: flex; flex-direction: column; gap: 2px; }

  .regimen-row {
    background: var(--paper-card);
    border: 1.5px solid var(--rule);
    border-radius: var(--radius);
    overflow: hidden;
    transition: box-shadow 0.2s, border-color 0.2s;
    cursor: pointer;
  }

  .regimen-row:hover { box-shadow: var(--shadow-hover); border-color: #ccc8c2; }

  .regimen-main {
    display: flex;
    align-items: center;
    padding: 18px 22px;
    gap: 16px;
  }

  .regimen-rank {
    font-family: var(--font-mono);
    font-size: 11px;
    color: var(--ink-faint);
    min-width: 28px;
    flex-shrink: 0;
    align-self: flex-start;
    padding-top: 3px;
  }

  .regimen-name-block { flex: 1; min-width: 0; }

  .regimen-name {
    font-family: var(--font-display);
    font-size: 18px;
    color: var(--ink);
    margin-bottom: 2px;
    white-space: normal;
    word-break: break-word;
  }

  .regimen-mechanism {
    font-size: 12px;
    color: var(--ink-faint);
    font-weight: 300;
    margin-top: 4px;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }

  .regimen-tags { display: flex; gap: 6px; flex-shrink: 0; align-self: flex-start; padding-top: 2px; flex-wrap: wrap; }

  .tag {
    font-family: var(--font-mono);
    font-size: 10px;
    letter-spacing: 0.08em;
    padding: 3px 8px;
    border-radius: 2px;
    font-weight: 500;
    white-space: nowrap;
  }

  .tag-high   { background: var(--success-light); color: var(--success); }
  .tag-medium { background: var(--mid-light); color: var(--mid); }
  .tag-low    { background: #fce8e4; color: var(--accent); }

  .regimen-meters { display: flex; gap: 20px; flex-shrink: 0; align-self: flex-start; padding-top: 2px; }

  .meter { display: flex; flex-direction: column; align-items: flex-end; gap: 5px; min-width: 80px; }

  .meter-labels { display: flex; justify-content: space-between; width: 100%; }

  .meter-title {
    font-family: var(--font-mono);
    font-size: 9px;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    color: var(--ink-faint);
  }

  .meter-value { font-family: var(--font-mono); font-size: 13px; font-weight: 500; color: var(--ink); }

  .bar-track { width: 80px; height: 3px; background: var(--rule); border-radius: 2px; overflow: hidden; }

  .bar-fill { height: 100%; border-radius: 2px; transition: width 0.6s cubic-bezier(0.16,1,0.3,1); }

  .bar-success { background: var(--success); }
  .bar-novelty { background: var(--accent); }

  .chevron { color: var(--ink-faint); font-size: 11px; flex-shrink: 0; transition: transform 0.2s; margin-left: 4px; align-self: flex-start; padding-top: 5px; }
  .chevron.open { transform: rotate(180deg); }

  .regimen-detail {
    border-top: 1px solid var(--rule);
    padding: 20px 22px 20px 66px;
    background: #faf9f7;
    display: flex;
    flex-direction: column;
    gap: 14px;
    animation: slideDown 0.18s ease;
  }

  @keyframes slideDown {
    from { opacity: 0; transform: translateY(-6px); }
    to   { opacity: 1; transform: translateY(0); }
  }

  .detail-row { display: flex; gap: 32px; flex-wrap: wrap; }

  .detail-stat { display: flex; flex-direction: column; gap: 2px; }

  .detail-stat-label {
    font-family: var(--font-mono);
    font-size: 9px;
    letter-spacing: 0.14em;
    text-transform: uppercase;
    color: var(--ink-faint);
  }

  .detail-stat-value { font-family: var(--font-mono); font-size: 15px; font-weight: 500; color: var(--ink); }

  .detail-one-liner {
    font-size: 13px;
    color: var(--ink-soft);
    font-weight: 300;
    line-height: 1.6;
    border-left: 2px solid var(--accent);
    padding-left: 12px;
  }

  .rationale-block { border-left: 2px solid #b8d4b8; padding-left: 12px; }

  .rationale-label {
    font-family: var(--font-mono);
    font-size: 9px;
    letter-spacing: 0.14em;
    text-transform: uppercase;
    color: var(--success);
    margin-bottom: 5px;
  }

  .rationale-text { font-size: 13px; color: var(--ink-soft); font-weight: 300; line-height: 1.7; }

  .rationale-loading { display: flex; align-items: center; gap: 8px; }
  .rationale-dots { display: flex; gap: 4px; }

  .r-dot {
    width: 4px; height: 4px;
    background: var(--success);
    border-radius: 50%;
    animation: pulse 1.2s ease-in-out infinite;
    opacity: 0.5;
  }
  .r-dot:nth-child(2) { animation-delay: 0.2s; }
  .r-dot:nth-child(3) { animation-delay: 0.4s; }

  @keyframes pulse {
    0%, 80%, 100% { opacity: 0.2; transform: scale(0.85); }
    40% { opacity: 0.8; transform: scale(1); }
  }

  .rationale-loading-text {
    font-family: var(--font-mono);
    font-size: 10px;
    color: var(--ink-faint);
    letter-spacing: 0.1em;
  }

  .filtered-section { margin-top: 32px; border-top: 1px solid var(--rule); padding-top: 16px; }

  .filtered-section summary {
    font-family: var(--font-mono);
    font-size: 10px;
    color: var(--ink-faint);
    letter-spacing: 0.1em;
    cursor: pointer;
    text-transform: uppercase;
  }

  .filtered-table { width: 100%; border-collapse: collapse; margin-top: 12px; }

  .filtered-table th {
    font-family: var(--font-mono);
    font-size: 9px;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    color: var(--ink-faint);
    text-align: left;
    padding: 6px 8px;
    border-bottom: 1px solid var(--rule);
  }

  .filtered-table td {
    font-family: var(--font-mono);
    font-size: 11px;
    color: var(--ink-soft);
    padding: 6px 8px;
    border-bottom: 1px solid var(--rule);
  }

  .loading-state { padding: 48px 0; display: flex; align-items: center; gap: 16px; }
  .loading-dots { display: flex; gap: 6px; }

  .dot {
    width: 5px; height: 5px;
    background: var(--ink);
    border-radius: 50%;
    animation: pulse 1.2s ease-in-out infinite;
  }
  .dot:nth-child(2) { animation-delay: 0.2s; }
  .dot:nth-child(3) { animation-delay: 0.4s; }

  .loading-text {
    font-family: var(--font-mono);
    font-size: 12px;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    color: var(--ink-faint);
  }

  .footer {
    margin-top: 80px;
    padding-top: 24px;
    border-top: 1px solid var(--rule);
    font-family: var(--font-mono);
    font-size: 10px;
    color: var(--ink-faint);
    letter-spacing: 0.1em;
    line-height: 1.8;
  }

  @media (max-width: 600px) {
    .regimen-meters { display: none; }
    .regimen-main { padding: 14px 16px; gap: 10px; }
    .regimen-detail { padding-left: 16px; }
    .search-row { flex-direction: column; }
    .filters { gap: 12px; }
  }
`;

/* ─── AI rationale fetcher ────────────────────────────────── */
async function fetchRationale(candidate, disease) {
  const res = await fetch(`${API_BASE}/generate_ai_analysis`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      target: candidate.drug_name || candidate.name,
      disease,
      mechanism: candidate.mechanism || candidate.drug_class || 'Not specified',
    }),
  });
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  const data = await res.json();
  return data.analysis ?? 'Rationale unavailable.';
}

/* ─── CandidateRow ────────────────────────────────────────── */
function CandidateRow({ candidate, rank, disease }) {
  const [open, setOpen] = useState(false);
  const [rationale, setRationale] = useState(null);
  const [rationaleLoading, setRationaleLoading] = useState(false);

  const name = candidate.drug_name || candidate.name || 'Unknown';
  const sr = successRate(candidate);
  const novelty = noveltyScore(candidate);
  const mechanism = candidate.mechanism || candidate.drug_class || '';
  const bbb = candidate.bbb_status || candidate.bbb;
  const bbbTag = bbb
    ? { HIGH: 'tag-high', MODERATE: 'tag-medium', LOW: 'tag-low' }[bbb]
    : null;

  const handleToggle = async () => {
    const next = !open;
    setOpen(next);
    if (next && !rationale && !rationaleLoading) {
      setRationaleLoading(true);
      try {
        const text = await fetchRationale(candidate, disease);
        setRationale(text);
      } catch {
        setRationale('Could not generate rationale at this time.');
      } finally {
        setRationaleLoading(false);
      }
    }
  };

  return (
    <div className="regimen-row">
      <div className="regimen-main" onClick={handleToggle}>
        <span className="regimen-rank">{String(rank).padStart(2, '0')}</span>

        <div className="regimen-name-block">
          <div className="regimen-name">{name}</div>
          {mechanism && <div className="regimen-mechanism">{mechanism}</div>}
        </div>

        <div className="regimen-tags">
          {bbbTag && <span className={`tag ${bbbTag}`}>BBB {bbb}</span>}
        </div>

        <div className="regimen-meters">
          <div className="meter">
            <div className="meter-labels">
              <span className="meter-title">Score</span>
              <span className="meter-value">{pct(sr)}</span>
            </div>
            <div className="bar-track">
              <div className="bar-fill bar-success" style={{ width: pct(sr) }} />
            </div>
          </div>
          <div className="meter">
            <div className="meter-labels">
              <span className="meter-title">Novelty</span>
              <span className="meter-value">{pct(novelty)}</span>
            </div>
            <div className="bar-track">
              <div className="bar-fill bar-novelty" style={{ width: pct(novelty) }} />
            </div>
          </div>
        </div>

        <span className={`chevron ${open ? 'open' : ''}`}>▾</span>
      </div>

      {open && (
        <div className="regimen-detail">
          <div className="detail-row">
            {candidate.score != null && (
              <div className="detail-stat">
                <span className="detail-stat-label">Composite Score</span>
                <span className="detail-stat-value">{candidate.score.toFixed(3)}</span>
              </div>
            )}
            {candidate.depmap_score != null && (
              <div className="detail-stat">
                <span className="detail-stat-label">DepMap</span>
                <span className="detail-stat-value">{candidate.depmap_score.toFixed(3)}</span>
              </div>
            )}
            {candidate.tissue_score != null && (
              <div className="detail-stat">
                <span className="detail-stat-label">Tissue</span>
                <span className="detail-stat-value">{candidate.tissue_score.toFixed(3)}</span>
              </div>
            )}
            {candidate.escape_score != null && (
              <div className="detail-stat">
                <span className="detail-stat-label">Escape Bypass</span>
                <span className="detail-stat-value">{candidate.escape_score.toFixed(3)}</span>
              </div>
            )}
            {candidate.ppi_score != null && (
              <div className="detail-stat">
                <span className="detail-stat-label">PPI</span>
                <span className="detail-stat-value">{candidate.ppi_score.toFixed(3)}</span>
              </div>
            )}
          </div>

          {candidate.indication && (
            <p className="detail-one-liner">{candidate.indication}</p>
          )}

          <div className="rationale-block">
            <div className="rationale-label">AI Biological Rationale</div>
            {rationaleLoading ? (
              <div className="rationale-loading">
                <div className="rationale-dots">
                  <div className="r-dot" /><div className="r-dot" /><div className="r-dot" />
                </div>
                <span className="rationale-loading-text">Generating rationale...</span>
              </div>
            ) : (
              <p className="rationale-text">{rationale}</p>
            )}
          </div>
        </div>
      )}
    </div>
  );
}

/* ─── App ─────────────────────────────────────────────────── */
const SUBGROUPS = ['Pan', 'TYR', 'SHH', 'MYC'];
const LOCATIONS = ['unknown_location', 'infratentorial', 'supratentorial', 'spinal'];

export default function App() {
  const [query, setQuery] = useState('');
  const [subgroup, setSubgroup] = useState('Pan');
  const [location, setLocation] = useState('unknown_location');
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState(null);
  const [error, setError] = useState(null);
  const [backendOk, setBackendOk] = useState(null);

  useEffect(() => {
    fetch(`${API_BASE}/`)
      .then(r => r.json())
      .then(d => setBackendOk(d.pipeline_ready ?? true))
      .catch(() => setBackendOk(false));
  }, []);

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!query.trim()) return;
    setLoading(true);
    setError(null);
    setResults(null);

    try {
      const res = await fetch(`${API_BASE}/analyze`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          disease_name: query.trim().toLowerCase(),
          min_score: 0.2,
          max_results: 15,
          subgroup: subgroup === 'Pan' ? undefined : subgroup,
          location,
        }),
      });
      const data = await res.json();
      if (!data.success) setError(data.error || 'Pipeline returned an error.');
      else setResults(data);
    } catch {
      setError(`Cannot reach backend at ${API_BASE}. Is uvicorn running?`);
    } finally {
      setLoading(false);
    }
  };

  const candidates = results?.candidates ?? [];

  return (
    <>
      <style>{styles}</style>
      <div className="app">

        <header className="header">
          <p className="header-eyebrow">
            <span
              className="status-dot"
              style={{
                background: backendOk === null ? '#ccc' : backendOk ? '#2a7a4b' : '#c8502a',
              }}
            />
            ATRT Drug Repurposing Engine · v1.0
          </p>
          <h1 className="header-title">
            Drug repurposing,<br /><em>reimagined.</em>
          </h1>
          <p className="header-sub">
            Enter a disease to surface ranked drug candidates scored by multi-omic pipeline —
            tissue expression, CRISPR essentiality, escape bypass, and PPI proximity.
          </p>
        </header>

        <div className="search-block">
          <label className="search-label" htmlFor="disease-input">Target disease</label>
          <form className="search-row" onSubmit={handleSubmit}>
            <input
              id="disease-input"
              type="text"
              className="search-input"
              value={query}
              onChange={e => setQuery(e.target.value)}
              placeholder="e.g. atrt, rhabdoid tumor, medulloblastoma"
              autoComplete="off"
              spellCheck="false"
            />
            <button type="submit" className="search-btn" disabled={loading || !query.trim()}>
              {loading ? 'Analysing...' : 'Analyse →'}
            </button>
          </form>
        </div>

        <div className="filters">
          <div className="filter-group">
            <span className="filter-label">Subgroup</span>
            <div className="pill-row">
              {SUBGROUPS.map(s => (
                <button
                  key={s}
                  className={`pill ${subgroup === s ? 'pill-active' : ''}`}
                  onClick={() => setSubgroup(s)}
                  type="button"
                >
                  {s === 'Pan' ? 'Pan-ATRT' : s}
                </button>
              ))}
            </div>
          </div>

          <div className="filter-group">
            <span className="filter-label">Location</span>
            <select
              className="filter-select"
              value={location}
              onChange={e => setLocation(e.target.value)}
            >
              {LOCATIONS.map(l => (
                <option key={l} value={l}>{l.replace(/_/g, ' ')}</option>
              ))}
            </select>
          </div>
        </div>

        {error && <div className="error-bar">Error: {error}</div>}

        {loading && (
          <div className="loading-state">
            <div className="loading-dots">
              <div className="dot" /><div className="dot" /><div className="dot" />
            </div>
            <span className="loading-text">Screening compounds · Scoring · Filtering</span>
          </div>
        )}

        {results && !loading && (
          <>
            <div className="results-header">
              <h2 className="results-title">{query}</h2>
              <span className="results-meta">
                {candidates.length} candidates · ranked by composite score
                {results.filtered_count > 0 && ` · ${results.filtered_count} safety-filtered`}
              </span>
            </div>

            <div className="regimen-list">
              {candidates.map((c, i) => (
                <CandidateRow
                  key={c.drug_name || i}
                  candidate={c}
                  rank={i + 1}
                  disease={query}
                />
              ))}
            </div>

            {results.filtered_drugs?.length > 0 && (
              <div className="filtered-section">
                <details>
                  <summary>{results.filtered_drugs.length} drugs removed by safety filter</summary>
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
          Data sources: GSE70678 · DepMap 24Q4 · STRING-DB · CMap L1000 · OpenTargets · ClinicalTrials.gov
          <br />
          Scoring: tissue × 0.40 + DepMap × 0.35 + escape × 0.20 + PPI × 0.05 · EZH2 ×1.40 (Knutson 2013)
        </footer>
      </div>
    </>
  );
}
