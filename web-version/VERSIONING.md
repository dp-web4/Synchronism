# Synchronism Web Version - Auto-Versioning System

## Overview

The Synchronism web version uses an automatic versioning system that combines a manually-controlled main version with dynamic timestamps.

## Version Format

```
{main-version}.yy.mm.dd.hh:mm
```

Example: `V0.25.06.27.03:45`

- **main-version**: Stored in `main-version.json`, manually updated for major releases
- **yy**: Two-digit year
- **mm**: Two-digit month  
- **dd**: Two-digit day
- **hh**: Two-digit hour (24-hour format)
- **mm**: Two-digit minutes

## Configuration

### main-version.json
```json
{
  "mainVersion": "V0",
  "lastUpdated": "2025-06-27T03:45:00Z"
}
```

### Manual Updates
To update the main version:
1. Edit `main-version.json`
2. Change `mainVersion` to the new version (e.g., "V1", "V2")
3. Update `lastUpdated` to current ISO timestamp

### Automatic Updates
The timestamp portion (`.yy.mm.dd.hh:mm`) updates automatically whenever:
- The page loads
- Content is modified dynamically
- `synchronismNav.refreshVersion()` is called

## Implementation Details

### JavaScript API
```javascript
// Access the navigation instance
synchronismNav.refreshVersion();  // Updates the displayed version

// The version is automatically updated on page load
// and can be manually refreshed after content changes
```

### HTML Structure
The version is displayed in elements with class `document-version`:
```html
<div class="version-badge document-version">V0.24.09.28.11.00</div>
```

## Usage Notes

1. **Page Load**: Version updates automatically when the page loads
2. **Content Changes**: Call `synchronismNav.refreshVersion()` after making dynamic content changes
3. **Main Version**: Update `main-version.json` for major releases or milestones
4. **Fallback**: If `main-version.json` cannot be loaded, defaults to "V0"

## Benefits

- **Traceability**: Each viewing shows exactly when it was accessed
- **Version Control**: Main version tracks major releases
- **Automatic Timestamps**: No manual timestamp updates needed
- **Dynamic Display**: Version reflects the current state at viewing time