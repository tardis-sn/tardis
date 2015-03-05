import tornado.ioloop
import tornado.web
import tornado.template
import os
from forms.forms import TardisForm
from jinja2 import Environment, FileSystemLoader

templates = tornado.template.Loader("templates")
TEMPLATE_FILE = "index.html"
templateLoader = FileSystemLoader( searchpath="templates/" )
templateEnv = Environment(loader=templateLoader,auto_reload=True)
template = templateEnv.get_template(TEMPLATE_FILE)

class MainHandler(tornado.web.RequestHandler):
    def get(self):
        form = TardisForm()
        # loader = tornado.template.Loader("templates")
        # self.write(loader.load("index.html").generate(form=form))
        self.write(template.render(form=form))
        # self.render("index.html",form=form)
    def post(self):
        form = TardisForm(self.request.arguments)
        details = ''
        if form.validate():
            for f in self.request.arguments:
                details += self.get_argument(f, default=None, strip=False)
            self.write(details)
        else:
            self.set_status(400)
            self.write(form.errors)

application = tornado.web.Application([
    (r"/", MainHandler),
    (r'/(.*)', tornado.web.StaticFileHandler, {'path': "static"})
    ],
    template_path="templates",
    debug=True)


if __name__ == "__main__":
    application.listen(int(os.environ.get("PORT", 5000)))
    tornado.ioloop.IOLoop.instance().start()
