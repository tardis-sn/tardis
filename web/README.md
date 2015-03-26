Tardis Web Interface
=====================

> This is a mock web intrface for genrating the yaml config files which are used to run simulations in tardis. For now the input elements were generated manually, I'm working on generating them using a schema file. The interface is live at http://tardis-web-mockup.herokuapp.com/

Local Installation
--------------
```sh
git clone https://github.com/drreddy/tardis.git
cd tardis
git checkout web-interface
pip install -r requirements.txt
python server.py
```

the application will be live at http://127.0.0.1:3000

### Work done so far
---
1. Created a full config file generator by using input elements based on the tardis_config.yml file
2. Deployed a test tardis server live on heroku at http://tardis-server.herokuapp.com/

### To Do (Future work)
---
1. Implement automatic input elements genration using a schema file and add client side validation.
2. Create a web service with tardis running in the back-end for running simulations and showing results on the client side, based on the config file genrated on client side
3. Create a deploy to  heroku feature from this repo.

A little bit of code explanation is done [here](https://github.com/drreddy/tardis/tree/web-interface/web/codeExplanation.md)
