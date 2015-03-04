import os
import StringIO

from flask import Flask, render_template, send_file
import yaml

from forms import SupernovaForm

app = Flask(__name__)

# WARNING:In production remove the secret key from here and set it to an
# environment variable
app.secret_key = 'secret-key-here'


@app.route('/', methods=['GET', 'POST'])
def index():
    supernova_form = SupernovaForm(csrf_enabled=False)
    if supernova_form.validate_on_submit():

        d = {
            'tardis_config_version': 'v1.0',
            'supernova': {
                'luminosity_requested': supernova_form.luminosity_requested.data,
                'time_explosion': supernova_form.time_explosion.data,
                'distance': supernova_form.distance.data,
                'luminosity_wavelength_start': supernova_form.luminosity_wavelength_start.data,
                'luminosity_wavelength_end': supernova_form.luminosity_wavelength_end.data,
            }
        }

        f = StringIO.StringIO()
        f.write(yaml.safe_dump(d, default_flow_style=False))
        f.seek(0)

        return send_file(f, attachment_filename='test.yml', as_attachment=True)

    else:
        return render_template('supernova.html', form=supernova_form)

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
