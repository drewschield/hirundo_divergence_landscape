import sys
from numpy import array
from numpy import median
from numpy import mean
from numpy import log2

auto = []
za_norm = []

for line in open(sys.argv[1], 'r'):
	acov = float(line.split()[3])
	auto.append(acov)

auto_med = median(array(auto))

for line in open(sys.argv[2], 'r'):
	zcov = line.split()[3]
	if float(zcov) > 0.00:
		za = log2(float(zcov)/float(auto_med))
		za_norm.append(za)

za_mean = mean(array(za_norm))

sample = str(sys.argv[1])
sample = sample.split('/')[1]
sample = sample.split('.')[0]

if za_mean < -0.2:
	message = 'is likely female'
	print za_mean, sample, message
else:
	message = 'is likely male'
	print za_mean, sample, message
