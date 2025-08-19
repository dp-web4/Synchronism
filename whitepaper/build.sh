#!/bin/bash

# build.sh - Main build script for Synchronism whitepaper
# Usage: 
#   ./build.sh          - Build all formats
#   ./build.sh md       - Build markdown only
#   ./build.sh pdf      - Build PDF only
#   ./build.sh web      - Build web version only
#   ./build.sh clean    - Clean build directory
#   ./build.sh rebuild  - Clean and rebuild everything

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}✓${NC} $1"
}

print_error() {
    echo -e "${RED}✗${NC} $1"
}

print_info() {
    echo -e "${BLUE}ℹ${NC} $1"
}

print_header() {
    echo ""
    echo -e "${YELLOW}═══════════════════════════════════════════════════${NC}"
    echo -e "${YELLOW}   $1${NC}"
    echo -e "${YELLOW}═══════════════════════════════════════════════════${NC}"
    echo ""
}

# Function to clean build directory
clean_build() {
    print_header "Cleaning build directory"
    if [ -d "build" ]; then
        rm -rf build/*
        print_status "Cleaned build directory"
    else
        print_info "Build directory doesn't exist"
    fi
}

# Function to build markdown
build_markdown() {
    print_header "Building Markdown"
    if bash make-md.sh; then
        print_status "Markdown build successful"
    else
        print_error "Markdown build failed"
        return 1
    fi
}

# Function to build PDF
build_pdf() {
    print_header "Building PDF"
    if bash make-pdf.sh; then
        print_status "PDF build successful"
    else
        print_error "PDF build failed"
        return 1
    fi
}

# Function to build web
build_web() {
    print_header "Building Web Version"
    if bash make-web-clean.sh; then
        print_status "Web build successful"
    else
        print_error "Web build failed"
        return 1
    fi
}

# Function to build all formats
build_all() {
    print_header "Building All Formats"
    
    build_markdown
    build_pdf
    build_web
    
    print_header "Build Complete"
    
    # Show summary
    echo "📦 Build Artifacts:"
    echo ""
    
    if [ -f "build/Synchronism_Whitepaper_Complete.md" ]; then
        size=$(du -h "build/Synchronism_Whitepaper_Complete.md" | cut -f1)
        echo "  📝 Markdown: build/Synchronism_Whitepaper_Complete.md ($size)"
    fi
    
    if [ -f "build/Synchronism_Whitepaper.pdf" ]; then
        size=$(du -h "build/Synchronism_Whitepaper.pdf" | cut -f1)
        echo "  📕 PDF:      build/Synchronism_Whitepaper.pdf ($size)"
    fi
    
    if [ -d "build/web-clean" ]; then
        count=$(find build/web-clean -name "*.html" | wc -l)
        echo "  🌐 Web:      build/web-clean/index.html ($count HTML files)"
    fi
    
    echo ""
    echo "📤 GitHub Pages Location:"
    echo "  • ../docs/whitepaper/"
    echo "    - Markdown: Synchronism_Whitepaper_Complete.md"
    echo "    - PDF:      Synchronism_Whitepaper.pdf"
    echo "    - Web:      index.html"
}

# Main script logic
case "$1" in
    clean)
        clean_build
        ;;
    md|markdown)
        build_markdown
        ;;
    pdf)
        build_pdf
        ;;
    web)
        build_web
        ;;
    rebuild)
        clean_build
        build_all
        ;;
    "")
        build_all
        ;;
    *)
        echo "Usage: $0 [clean|md|pdf|web|rebuild]"
        echo ""
        echo "Options:"
        echo "  clean    - Clean build directory"
        echo "  md       - Build markdown only"
        echo "  pdf      - Build PDF only"
        echo "  web      - Build web version only"
        echo "  rebuild  - Clean and rebuild everything"
        echo "  (none)   - Build all formats"
        exit 1
        ;;
esac

echo ""
print_info "Build timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""