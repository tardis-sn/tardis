__author__ = 'Saurabh'
import os
import tornado.web
import tornado.template
from util import ordered_load, get_label_from_property

class FormHandler(tornado.web.RequestHandler):
    def get(self):
        yaml_config = ordered_load(open(os.path.join(os.path.dirname(__file__), 'tardis_config_file.yml')).read())
        self.render(
            'config_form.html',
            yaml_config=yaml_config,
            yaml_config_name='root',
            render_form=self.render_form,
            get_label_from_property=get_label_from_property
        )

    def render_form(self, yaml_config, yaml_config_name):
        return self.render_string(
            'config_form_generator.html',
            yaml_config=yaml_config,
            yaml_config_name=yaml_config_name,
            render_form=self.render_form,
            get_label_from_property=get_label_from_property
        )
