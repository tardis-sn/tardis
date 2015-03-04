from flask.ext.wtf import Form
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired

class SupernovaForm(Form):
    luminosity_requested = StringField('luminosity requested', validators=[DataRequired()])
    time_explosion = StringField('time explosion', validators=[DataRequired()])
    distance = StringField('distance')
    luminosity_wavelength_start = StringField(' luminosity wavelength start')
    luminosity_wavelength_end = StringField('luminosity wavelength end')
    submit = SubmitField('Submit')
