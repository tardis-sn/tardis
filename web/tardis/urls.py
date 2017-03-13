from django.conf.urls import patterns, include, url

urlpatterns = patterns('',
    url(r'^$', 'tardis.views.home', name='home'),
)
