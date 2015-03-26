# Code Documentation

> A small explanation is provided on how the code works

1. A global js object called jsonConf is present which contains all the user inputs.
2. The user inputs are accessed using on change events on the input elements in jquery
3. Finally to create the yaml file the js object is used and is converted to a yaml file using blob web api and yaml.js

The structure of the js object is maintained using the following data-* tags present in the input element:

1. data-parentTags: Used to maintain the hierarchy according to the dafault config file
2. data-field: Used to represent the field name of the current input element in the final yaml file
3. data-type: Used to represent the type of the field (i.e int,float,quantity,string)

Finally to tackle the container fields :
1. Initially the containers and not shown
2. Once a user chooses a container the required container is shown using data-targetCont tag in jquery
3. Finally in the js object the field associated with the containers are recreated every time a different container is selected