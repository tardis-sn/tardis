__author__ = 'Saurabh'
import os

TEMPLATE_DIR= os.path.join(os.path.dirname(__file__), "template")

STATIC_FILES_DIR = os.path.join(os.path.dirname(__file__), "static")

settings = {
    "static_path": STATIC_FILES_DIR,
    "template_path": TEMPLATE_DIR,
    "cookie_secret": "__TODO:_GENERATE_YOUR_OWN_RANDOM_VALUE_HERE__",
    "login_url": "/login",
    "xsrf_cookies": True,
    "debug": True
}