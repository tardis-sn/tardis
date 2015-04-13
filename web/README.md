Web Interface for TARDIS
========================

This is a basic schema to generate a web-form for TARDIS.

A live version is available [here](https://whispering-atoll-6188.herokuapp.com/).

Working/Not-Working
-------------------

###Working:
- Form generation from the yaml config file
- YAML data generation from the form

###Not Working:
- Input Validation
- Using default values if data not present
- Start, Stop, Num fields in quantity_ranged_sampled values.

Deploying
---------

Run `python app.py` and open [http://127.0.0.1:5000](http://127.0.0.1:5000) in your browser.

Dependencies
------------
- Torndado
- Jinja
- WTForms

