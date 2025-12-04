-- Synchronism Research Memory Schema
-- Created: December 3, 2025
-- Purpose: Accumulated research knowledge for Synchronism project
-- Pattern: Modeled after SAGE memory with SNARC salience scoring

PRAGMA foreign_keys = ON;

-- ============================================================================
-- Core Findings Table
-- ============================================================================

CREATE TABLE IF NOT EXISTS findings (
    id TEXT PRIMARY KEY,

    -- Temporal
    created TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    updated TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    session_id TEXT,                    -- e.g., "78", "79"

    -- Hierarchical Domain (fractal contextualization)
    domain TEXT NOT NULL,               -- cosmology, galaxy_physics, coherence_math, implementation, validation
    subdomain TEXT,                     -- e.g., under galaxy_physics: rotation_curves, tully_fisher, scaling_relations

    -- Content
    title TEXT NOT NULL,
    summary TEXT NOT NULL,
    details TEXT,                       -- Full explanation, derivation, etc.

    -- Type & Status
    finding_type TEXT NOT NULL CHECK(finding_type IN (
        'derivation',       -- Mathematical derivation from first principles
        'validation',       -- Empirical test result
        'prediction',       -- Testable prediction
        'connection',       -- Link to other theory (MOND, etc.)
        'methodology',      -- Research method or standard
        'failure',          -- What didn't work (lessons)
        'insight',          -- Conceptual breakthrough
        'question'          -- Open question for future work
    )),

    status TEXT DEFAULT 'active' CHECK(status IN (
        'active',           -- Current understanding
        'superseded',       -- Replaced by newer finding
        'deprecated',       -- No longer believed correct
        'speculative'       -- Hypothesis, not validated
    )),

    validation_level TEXT DEFAULT 'theoretical' CHECK(validation_level IN (
        'proven',           -- Mathematically proven
        'validated',        -- Tested on data (e.g., SPARC)
        'theoretical',      -- Derived but not tested
        'speculative',      -- Proposed, needs work
        'falsified'         -- Tested and failed
    )),

    -- SNARC Salience Scores (0.0 - 1.0)
    surprise REAL DEFAULT 0.5 CHECK(surprise BETWEEN 0.0 AND 1.0),
    novelty REAL DEFAULT 0.5 CHECK(novelty BETWEEN 0.0 AND 1.0),
    arousal REAL DEFAULT 0.5 CHECK(arousal BETWEEN 0.0 AND 1.0),
    reward REAL DEFAULT 0.5 CHECK(reward BETWEEN 0.0 AND 1.0),
    conflict REAL DEFAULT 0.0 CHECK(conflict BETWEEN 0.0 AND 1.0),

    -- Composite salience (weighted sum - matches SAGE weights)
    salience REAL GENERATED ALWAYS AS (
        0.25 * surprise +
        0.20 * novelty +
        0.30 * arousal +
        0.10 * reward +
        0.15 * conflict
    ) STORED,

    -- Decay & Access
    decay_rate REAL DEFAULT 0.01,
    last_accessed TIMESTAMP,
    access_count INTEGER DEFAULT 0,

    -- References
    source_files TEXT,                  -- JSON array of related files
    related_findings TEXT,              -- JSON array of finding IDs
    tags TEXT                           -- JSON array of keywords
);

-- ============================================================================
-- Parameters Table (Synchronism-specific)
-- ============================================================================

CREATE TABLE IF NOT EXISTS parameters (
    name TEXT PRIMARY KEY,
    symbol TEXT,                        -- e.g., "γ", "A", "B"

    -- Values
    value_derived REAL,
    value_empirical REAL,
    value_current REAL,                 -- What we're using now

    -- Derivation status
    derivation_status TEXT CHECK(derivation_status IN (
        'first_principles',             -- Fully derived
        'semi_empirical',               -- Partially derived
        'empirical',                    -- Fitted from data
        'unknown'                       -- Not yet understood
    )),
    derivation_session TEXT,            -- Session that derived it
    derivation_summary TEXT,            -- How it was derived

    -- Physical meaning
    physical_meaning TEXT,
    units TEXT,

    -- Metadata
    created TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- Predictions Table (Testable claims)
-- ============================================================================

CREATE TABLE IF NOT EXISTS predictions (
    id TEXT PRIMARY KEY,

    -- Content
    title TEXT NOT NULL,
    description TEXT NOT NULL,
    quantitative_claim TEXT,            -- e.g., "0.28 dex offset"

    -- Testability
    test_method TEXT,                   -- How to test
    required_data TEXT,                 -- What data needed
    discriminates_from TEXT,            -- What theory it distinguishes from (MOND, ΛCDM)

    -- Status
    status TEXT DEFAULT 'untested' CHECK(status IN (
        'untested',
        'partially_tested',
        'confirmed',
        'falsified',
        'inconclusive'
    )),
    test_result TEXT,
    test_session TEXT,

    -- Priority (for research focus)
    priority TEXT DEFAULT 'medium' CHECK(priority IN ('critical', 'high', 'medium', 'low')),

    -- SNARC (for pruning)
    salience REAL DEFAULT 0.5,

    -- Metadata
    source_session TEXT,
    created TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- ============================================================================
-- Validation Results (SPARC tests, etc.)
-- ============================================================================

CREATE TABLE IF NOT EXISTS validations (
    id TEXT PRIMARY KEY,

    -- What was tested
    test_name TEXT NOT NULL,
    dataset TEXT,                       -- e.g., "SPARC", "Santos-Santos"

    -- Parameters used
    parameters_json TEXT,               -- JSON of parameter values

    -- Results
    success_rate REAL,
    median_chi_squared REAL,
    sample_size INTEGER,

    -- Comparison
    comparison_model TEXT,              -- What we compared against
    comparison_result TEXT,             -- Better/worse/same

    -- Metadata
    session_id TEXT,
    created TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    notes TEXT
);

-- ============================================================================
-- Cross-References (connections between findings)
-- ============================================================================

CREATE TABLE IF NOT EXISTS connections (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    finding_a TEXT NOT NULL,
    finding_b TEXT NOT NULL,
    relationship TEXT,                  -- "supports", "contradicts", "extends", "supersedes"
    notes TEXT,
    created TIMESTAMP DEFAULT CURRENT_TIMESTAMP,

    FOREIGN KEY (finding_a) REFERENCES findings(id),
    FOREIGN KEY (finding_b) REFERENCES findings(id)
);

-- ============================================================================
-- Indexes for fast retrieval
-- ============================================================================

CREATE INDEX IF NOT EXISTS idx_findings_domain ON findings(domain);
CREATE INDEX IF NOT EXISTS idx_findings_type ON findings(finding_type);
CREATE INDEX IF NOT EXISTS idx_findings_salience ON findings(salience DESC);
CREATE INDEX IF NOT EXISTS idx_findings_session ON findings(session_id);
CREATE INDEX IF NOT EXISTS idx_predictions_status ON predictions(status);
CREATE INDEX IF NOT EXISTS idx_predictions_priority ON predictions(priority);

-- ============================================================================
-- Full-text search on findings
-- ============================================================================

CREATE VIRTUAL TABLE IF NOT EXISTS findings_fts USING fts5(
    title,
    summary,
    details,
    tags,
    content='findings',
    content_rowid='rowid'
);

-- Triggers to keep FTS in sync
CREATE TRIGGER IF NOT EXISTS findings_ai AFTER INSERT ON findings BEGIN
    INSERT INTO findings_fts(rowid, title, summary, details, tags)
    VALUES (NEW.rowid, NEW.title, NEW.summary, NEW.details, NEW.tags);
END;

CREATE TRIGGER IF NOT EXISTS findings_ad AFTER DELETE ON findings BEGIN
    INSERT INTO findings_fts(findings_fts, rowid, title, summary, details, tags)
    VALUES('delete', OLD.rowid, OLD.title, OLD.summary, OLD.details, OLD.tags);
END;

CREATE TRIGGER IF NOT EXISTS findings_au AFTER UPDATE ON findings BEGIN
    INSERT INTO findings_fts(findings_fts, rowid, title, summary, details, tags)
    VALUES('delete', OLD.rowid, OLD.title, OLD.summary, OLD.details, OLD.tags);
    INSERT INTO findings_fts(rowid, title, summary, details, tags)
    VALUES (NEW.rowid, NEW.title, NEW.summary, NEW.details, NEW.tags);
END;
