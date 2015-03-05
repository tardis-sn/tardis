import os

from tornado.ioloop import IOLoop
from tornado.web import Application, RequestHandler, StaticFileHandler
from tornado.httpserver import HTTPServer

TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), "templates")
STATIC_PATH = os.path.join(os.path.dirname(__file__), "static")

class MainHandler(RequestHandler):
    def get(self):
        self.render("index.html")

def make_app():
    return Application([
            (r"/", MainHandler),
            (r'/static/(.*)', StaticFileHandler, {'path': STATIC_PATH})
        ],
        template_path = TEMPLATE_PATH, debug=True)

def main():
    app = make_app()
    http_server = HTTPServer(app)
    port = int(os.environ.get("PORT", 8888))
    http_server.listen(port)
    IOLoop.instance().start()
    

if __name__ == "__main__":
    main()