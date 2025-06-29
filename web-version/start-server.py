#!/usr/bin/env python3
"""
Simple HTTP server for Synchronism web interface
Run from the web-version directory: python start-server.py
"""

import http.server
import socketserver
import os
import sys

# Change to the directory containing this script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

PORT = 8000

class MyHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Add CORS headers to allow local file access
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        super().end_headers()

Handler = MyHTTPRequestHandler

try:
    with socketserver.TCPServer(("", PORT), Handler) as httpd:
        print(f"✅ Synchronism Web Server starting on port {PORT}")
        print(f"🌐 Open: http://localhost:{PORT}")
        print(f"📁 Serving from: {os.getcwd()}")
        print(f"🛑 Press Ctrl+C to stop")
        httpd.serve_forever()
except KeyboardInterrupt:
    print("\n🛑 Server stopped")
except OSError as e:
    if e.errno == 98:  # Address already in use
        print(f"❌ Port {PORT} is already in use. Try a different port or stop the existing server.")
    else:
        print(f"❌ Error starting server: {e}")
    sys.exit(1)