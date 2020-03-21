import pandas as pd
import hiplot
import pylab
import matplotlib.pyplot as plot

def conv_plus_plot(filename):
  read_file = pd.read_csv(filename) #input the file path which maybe a text document, hdf5 file, or any other format accepted by pandas
  read_file.to_csv(filename, index=None) #input the target path for csv read_file
  
  with open(filename) as f: #Enter the csv file name
      a = [{k: v for k,v in row.items()} for row in csv.DictReader(f,skipinitialspace=True)]
  #execute on ipython --pylab
  hip.Experiment.from_iterable(a).display(force_full_width=True)
  pylab.show()
  plt.show()
