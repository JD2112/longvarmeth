#!/usr/bin/env python3
import http.server
import socketserver
import sys

class JBrowseHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Enable CORS
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'X-Requested-With, Content-Type, Range')
        # Allow Range requests for JBrowse to fetch chunks of BAM/VCF files
        self.send_header('Accept-Ranges', 'bytes')
        # Disable caching to force browser to apply new headers immediately
        self.send_header('Cache-Control', 'no-cache, no-store, must-revalidate')
        self.send_header('Pragma', 'no-cache')
        self.send_header('Expires', '0')
        super().end_headers()

    def do_OPTIONS(self):
        # Handle preflight CORS requests
        self.send_response(200)
        self.end_headers()

    def guess_type(self, path):
        # Prevent Python from adding Content-Encoding: gzip to .gz/tbi files,
        # which causes the browser to decompress them automatically, breaking JBrowse's bgzf reader.
        if path.endswith('.gz') or path.endswith('.tbi') or path.endswith('.bgz') or path.endswith('.csi'):
            return 'application/octet-stream'
        return super().guess_type(path)

class ThreadingHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True

def main():
    port = 8000
    if len(sys.argv) > 1:
        try:
            port = int(sys.argv[1])
        except ValueError:
            print(f"Invalid port: {sys.argv[1]}. Defaulting to 8000.", file=sys.stderr)
            
    handler = JBrowseHTTPRequestHandler
    
    # Allow port reuse
    ThreadingHTTPServer.allow_reuse_address = True
    with ThreadingHTTPServer(("", port), handler) as httpd:
        print(f"Serving JBrowse files on port {port}...")
        print(f"CORS is enabled, cache-control is disabled, and the server is MULTI-THREADED.")
        print(f"Open your JBrowse instance and configure tracks using: http://localhost:{port}/...")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nServer stopped.")

if __name__ == '__main__':
    main()
