import os
import json
import yaml

from tornado.ioloop import IOLoop
from tornado.web import Application, RequestHandler, StaticFileHandler
from tornado.httpserver import HTTPServer

TEMPLATE_PATH = os.path.join(os.path.dirname(__file__), "templates")
STATIC_PATH = os.path.join(os.path.dirname(__file__), "static")

class MainHandler(RequestHandler):
    def get(self):
        self.render("index.html")

class YamlHandler(RequestHandler):
    def post(self):
        config_file = file('input_data.yml', 'w')
        config_file.write('tardis_config_version: v1.0\n')
        data = json.loads(self.request.body)
        yml = yaml.safe_dump(data, config_file)
        config_file.close()
        self.set_header("Content-Type", "application/json")
        self.set_status(201)

def make_app():
    return Application([
            (r"/", MainHandler),
            (r"/yaml", YamlHandler),
            (r'/static/(.*)', StaticFileHandler, {'path': STATIC_PATH})
        ],
        template_path = TEMPLATE_PATH)

def main():
    app = make_app()
    http_server = HTTPServer(app)
    port = int(os.environ.get("PORT", 5000))
    http_server.listen(port)
    print "server listen"
    IOLoop.instance().start()
    

if __name__ == "__main__":
    main()