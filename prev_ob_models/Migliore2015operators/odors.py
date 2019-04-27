import fileinput

class odor():
  def __init__(self, index, name, glom_weights):
    self.index = index
    self.name = name
    self.glom_weights = glom_weights

odors = {} # by name

for line in fileinput.input('input-odors.txt'):
  data = line.split('\t')
  odors.update({data[0]: odor(fileinput.lineno(), data[0], [float(i) for i in data[1:]])})

if __name__ == '__main__':
  for name in odors:
    print name, odors[name].index
  from mayavi.mlab import barchart, show
  barchart([odors[name].glom_weights for name in odors])
  show()
