:: only really needed for win and py27,
:: we can drop this line when we drop py27
del /f /q NEWS.txt
%PYTHON% -m pip install . --no-deps -vv
