import os

from tornado.ioloop import IOLoop
from tornado.web import Application, RequestHandler
from tornado.httpserver import HTTPServer

class IndexHandler(RequestHandler):
    def get(self):
        self.render("index.html")

def main():
    settings = {
        "template_path":os.path.join(os.path.dirname(__file__), "templates"),
        "static_path": os.path.join(os.path.dirname(__file__), "static"),
        "debug": True
    }
    application = Application([
        (r"/", IndexHandler)
    ], **settings)
    http_server = HTTPServer(application)
    port = int(os.environ.get("PORT", 3000))
    http_server.listen(port)
    IOLoop.instance().start()

if __name__ == "__main__":
    main()