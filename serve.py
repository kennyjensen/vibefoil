"""Serve Vibefoil locally with the same custom 404 page as GitHub Pages."""

import argparse
from functools import partial
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path


class PagesRequestHandler(SimpleHTTPRequestHandler):
    def send_error(self, code, message=None, explain=None):
        if code != HTTPStatus.NOT_FOUND:
            return super().send_error(code, message, explain)
        try:
            body = (Path(self.directory) / "404.html").read_bytes()
        except OSError:
            return super().send_error(code, message, explain)

        self.send_response(code, message)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        if self.command != "HEAD":
            self.wfile.write(body)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8080)
    parser.add_argument("--bind", default="127.0.0.1")
    args = parser.parse_args()
    handler = partial(PagesRequestHandler, directory=str(Path(__file__).resolve().parent))
    with ThreadingHTTPServer((args.bind, args.port), handler) as server:
        print(f"Serving Vibefoil at http://{args.bind}:{server.server_port}", flush=True)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass


if __name__ == "__main__":
    main()
