from TardisForm.views import FormHandler
import os
import tornado.web

urls = [
    (r"/form", FormHandler),
    (r'/static/(.*)', tornado.web.StaticFileHandler, {'path': os.path.join(os.path.dirname(__file__), 'static')})
]