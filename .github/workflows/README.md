# Synchronism GitHub Actions Workflows

## Workflows Overview

### 1. `update_synchronism.yml` - Governance Auto-Update
- **Trigger**: Daily at 12:00 UTC or manual
- **Purpose**: Runs the AI governance system to propose and review whitepaper improvements
- **Actions**:
  - Runs `scripts/governance/main.py update`
  - Uses AI participants (Claude, GPT) to review sections
  - Commits approved changes automatically
- **Secrets Required**:
  - `SYNC_GH_PAT` - GitHub Personal Access Token
  - `CLAUDE_API` - Anthropic API key
  - `GPT_API` - OpenAI API key

### 2. `build_whitepaper.yml` - Whitepaper Builder
- **Trigger**: 
  - Push to main branch (when whitepaper files change)
  - Pull requests affecting whitepaper
  - Manual trigger
- **Purpose**: Build all whitepaper formats and deploy to GitHub Pages
- **Actions**:
  - Builds Markdown version
  - Builds PDF version (with pandoc/LaTeX)
  - Builds Web version (HTML + assets)
  - Uploads artifacts for each format
  - Deploys to GitHub Pages (`docs/` directory)
  - Generates build summary
- **Features**:
  - CI-safe (skips git pull operations)
  - Continues on PDF build issues (LaTeX can be finicky)
  - Creates detailed build summary
  - Preserves all artifacts

## Local Testing

To test the build process locally:

```bash
cd whitepaper

# Test individual builds
bash make-md.sh
bash make-pdf.sh
bash make-web-clean.sh

# For CI environment simulation
CI=true bash make-md.sh
```

## Secrets Configuration

Required repository secrets:
1. `SYNC_GH_PAT` - Personal Access Token with repo permissions
2. `CLAUDE_API` - Anthropic API key for Claude
3. `GPT_API` - OpenAI API key for GPT

## Artifacts

Build artifacts are available for download from the Actions tab:
- `synchronism-markdown` - Complete markdown file
- `synchronism-pdf` - PDF document
- `synchronism-web` - Web version with all assets

## GitHub Pages

The whitepaper is automatically deployed to:
https://dp-web4.github.io/Synchronism/whitepaper/

## Troubleshooting

### PDF Build Issues
LaTeX dependencies can sometimes fail. The workflow continues despite PDF issues to ensure other formats are built.

### Git Conflicts
The build scripts include conflict detection but in CI, we skip git operations entirely since Actions handles checkout.

### Missing Artifacts
Check the build logs in the Actions tab. Each build step logs its status.