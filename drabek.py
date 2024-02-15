import sys
import re
import random
from math import sqrt

###############################################################################

def mmin(x, y):
  if x is None:
    return y
  else:
    return min(x, y)

def mmax(x, y):
  if x is None:
    return y
  else:
    return max(x, y)

###############################################################################

def least_squares_linear_fit(points):
  sum_x = sum_y = sum_x2 = sum_xy = 0.0
  sum_1 = float(len(points))
  for x, y in points:
    sum_x += x
    sum_y += y
    sum_x2 += x * x
    sum_xy += x * y

  determinant = sum_x * sum_x - sum_1 * sum_x2
  slope = (sum_y * sum_x - sum_1 * sum_xy) / determinant
  intercept = (sum_x * sum_xy - sum_y * sum_x2) / determinant

  return slope, intercept

def pearson_r(x_y):
  mean_x = mean_y = 0.0
  for x, y in x_y:
    mean_x += x
    mean_y += y
  mean_x /= float(len(x_y))
  mean_y /= float(len(x_y))

  s1 = s2 = s3 = 0.0
  for x, y in x_y:
    dx = x - mean_x
    dy = y - mean_y
    s1 += dx * dy
    s2 += dx * dx
    s3 += dy * dy

  return s1 / sqrt(s2 * s3)

###############################################################################

def all(iterable):
  for element in iterable:
    if not element:
      return False
  return True

def any(iterable):
  for element in iterable:
    if element:
      return True
  return False

###############################################################################

class ID_gen:
  def __init__(self, prefix='', first_id=0):
    self.prefix = prefix
    self.id = first_id - 1

  def __call__(self):
    self.id += 1
    return '%s%s' % (self.prefix, self.id)

###############################################################################

def incr(d, x, c=1):
  try:
    d[x] += c
  except:
    d[x] = c

################################################################################
# argmax and family

def index_min(l):
  """returns the index of the smallest element in l
     (the first if there are ties)"""
  if len(l) == 0:
    raise ValueError(l)
  best_i = 0
  best_val = l[0]
  for i in range(1, len(l)):
    if l[i] < best_val:
      best_val = l[i]
      best_i = i
  return best_i

def index_max(l):
  """returns the index of the largest element in l
     (the first if there are ties)"""
  if len(l) == 0:
    raise ValueError(l)
  best_i = 0
  best_val = l[0]
  for i in range(1, len(l)):
    if l[i] > best_val:
      best_val = l[i]
      best_i = i
  return best_i

def index_min_random_tie(l):
  """returns the index of the smallest element in l
     (a randomly selected one if there are ties)"""
  if len(l) == 0:
    raise ValueError(l)
  best_i = 0
  best_val = l[0]
  num_same = 1
  for i in range(1, len(l)):
    if l[i] < best_val:
      best_val = l[i]
      best_i = i
      num_same = 1
    elif l[i] == best_val:
      num_same += 1
      if random.randrange(num_same) == 0:
        best_i = i
  return best_i

def index_max_random_tie(l):
  """returns the index of the largest element in l
     (a randomly selected one if there are ties)"""
  if len(l) == 0:
    raise ValueError(l)
  best_i = 0
  best_val = l[0]
  num_same = 1
  for i in range(1, len(l)):
    if l[i] > best_val:
      best_val = l[i]
      best_i = i
      num_same = 1
    elif l[i] == best_val:
      num_same += 1
      if random.randrange(num_same) == 0:
        best_i = i
  return best_i

def index_min_list(d):
  """returns a list of the indices of the (tied) largest elements in l"""
  assert 0

def index_max_list(d):
  """returns a list of the indices of the (tied) largest elements in l"""
  return index_min_list(map(lambda x: -x, l))

def key_min(d):
  """returns the key whose value is the smallest in d
     (the first if there are ties)"""
  if len(d) == 0:
    raise ValueError(d)
  best_key = None
  best_val = None
  first_round = 1
  for key, val in d.iteritems():
    if first_round or val < best_val:
      first_round = 0
      best_key = key
      best_val = val
  return best_key

def key_min_random_tie(d):
  """returns the key whose value is the smallest in d
  (a randomly selected one if there are ties)"""
  if len(d) == 0:
    raise ValueError(d)
  best_key = None
  best_val = None
  num_same = None
  for key, val in d.iteritems():
    if num_same == None or val < best_val:
      best_key = key
      best_val = val
      num_same = 1
    elif val == best_val:
      num_same += 1
      if random.randrange(num_same) == 0:
        best_key = key
  return best_key

def key_min_list(d):
  """returns a list of the indices of the (tied) largest elements in l"""
  if len(d) == 0:
    raise ValueError(d)
  best_keys = None
  best_val = None
  for key, val in d.iteritems():
    if best_keys == None or val < best_val:
      best_keys = [key]
      best_val = val
    elif val == best_val:
      best_keys.append(key)
  return best_keys

def key_max(d):
  """returns the key whose value is the largest in d
     (the first if there are ties)"""
  if len(d) == 0:
    raise ValueError(d)
  best_key = None
  best_val = None
  for key, val in d.iteritems():
    if best_val == None or val > best_val:
      best_key = key
      best_val = val
  return best_key

def key_max_random_tie(d):
  """returns the index of the smallest element in l
     (a randomly selected one if there are ties)"""
  if len(d) == 0:
    raise ValueError(d)
  best_key = None
  best_val = None
  num_same = None
  for key, val in d.iteritems():
    if num_same == None or val > best_val:
      best_key = key
      best_val = val
      num_same = 1
    elif val == best_val:
      num_same += 1
      if random.randrange(num_same) == 0:
        best_key = key
  return best_key

def key_max_list(d):
  """returns a list of the indices of the (tied) largest elements in l"""
  if len(d) == 0:
    raise ValueError(d)
  best_keys = None
  best_val = None
  for key, val in d.iteritems():
    if best_keys == None or val > best_val:
      best_keys = [key]
      best_val = val
    elif val == best_val:
      best_keys.append(key)
  return best_keys

################################################################################

def read_fasta(file=sys.stdin):
  id_line_pat = re.compile('^>\\s*(\\S+).*')
  id = None
  sequence_parts = []
  for line in file:
    match = id_line_pat.match(line)
    if match:
      if id is not None or sequence_parts:
        yield id, ''.join(sequence_parts)
      id, = match.groups()
      sequence_parts = []
    else:
      sequence_parts.append(line.strip())
  yield id, ''.join(sequence_parts)

def read_fasta_fancy(file=sys.stdin):
  id_line_pat = re.compile('^>\\s*(\\S+)(.*)')
  id = None
  rest_of_head_line = None
  sequence_parts = []
  for line in file:
    match = id_line_pat.match(line)
    if match:
      if id is not None or sequence_parts:
        yield id, rest_of_head_line, ''.join(sequence_parts)
      id, rest_of_head_line = match.groups()
      sequence_parts = []
    else:
      sequence_parts.append(line.strip())
  yield id, rest_of_head_line, ''.join(sequence_parts)

def read_fasta_agnostic(file=sys.stdin):
  id = None
  sequence_parts = []
  for line in file:
    if line.startswith('>'):
      if id is not None or sequence_parts:
        yield id, ''.join(sequence_parts)
      id = line[1:-1]
      sequence_parts = []
    else:
      sequence_parts.append(line.strip())
  yield id, ''.join(sequence_parts)

def write_fasta_entry(label, sequence, file=sys.stdout, max_line_length=60):
  print('>' + label, file=file)
  position = 0
  while position < len(sequence):
    print(sequence[position:position + max_line_length], file=file)
    position += max_line_length

################################################################################

def mean(xs):
  return sum(xs) / float(len(xs))

def variance(xs):
  average = mean(xs)
  return mean([(x - average) ** 2 for x in xs])

def standard_deviation(xs):
  return sqrt(variance(xs))

def median(xs):
  xs = sorted(xs)
  halfway = len(xs) / 2.0
  if halfway % 1 == 0:
    return (xs[int(halfway) - 1] + xs[int(halfway)]) / 2.0
  else:
    return xs[int(halfway)]

################################################################################

def print_tabbed(*fields):
  print('\t'.join(map(str, fields)))

################################################################################

def read_phylip_dist(file):
  result = {}
  num_objects = None
  object_names = []
  for line in file:
    fields = line.split()
    if num_objects is None:
      num_objects, = map(int, fields)
      continue
    new_object_name = fields[0]
    distances = map(float, fields[1:])
    assert len(object_names) == len(distances), (len(object_names), len(distances))
    for old_object_name, distance in zip(object_names, distances):
      result[new_object_name, old_object_name] = distance
      result[old_object_name, new_object_name] = distance
    object_names.append(new_object_name)
  return result

################################################################################

def get_quantile(data, fraction=0.5):
  data = list(sorted(data))
  position = fraction * (len(data) - 1)

  if position % 1 == 0.0:
    return data[int(position)]
  else:
    return 0.5 * (data[int(position)] + data[int(position) + 1])

def find_approximate_quantiles(data, num_bins=2):
  target_bin_size = len(data) / float(num_bins)
  target_bin_ratio = 1.0 / num_bins
  x = 1 / float(len(data))

  data.append(data[-1] + 1)

  bin_num = 1
  prev_datum = None
  encompassed = 0
  for i, datum in enumerate(data):
    if prev_datum is None:
      bin_lower = datum
    else:
      if datum != prev_datum and i >= target_bin_size * bin_num:
        #print '\t'.join(map(str, (bin_lower, prev_datum, prev_datum - bin_lower, i - encompassed)))
        yield bin_lower, prev_datum
        bin_lower = datum
        encompassed = i
        bin_num += 1

    prev_datum = datum

################################################################################

def get_reverse_complement(seq, base_complement={}):
  if not base_complement:
    for base_1, base_2 in (
        ('a', 't'),
        ('c', 'g'),
        ('w', 'w'),
        ('s', 's'),
        ('m', 'k'),
        ('r', 'y'),
        ('b', 'v'),
        ('d', 'h'),
        ('n', 'n'),
        ('-', '-'),
        ('.', '.')):
      base_complement[base_1] = base_2
      base_complement[base_2] = base_1
      base_complement[base_1.upper()] = base_2.upper()
      base_complement[base_2.upper()] = base_1.upper()
  return ''.join(base_complement[base] for base in reversed(seq))

def translate_codon(codon, translation_table={}, translation_table_file_name='/home/edrabek/misc/standard-dna-code'):
  if not translation_table:
    translation_table_file = open(translation_table_file_name)
    for line in translation_table_file:
      codon, amino_acid = line[:-1].split('\t')
      translation_table[codon] = amino_acid

  return translation_table[codon]

