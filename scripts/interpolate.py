
def loadLogspaceMatrix(fname):
	with open(fname) as f:
		mat = [[math.exp(float(x)) for x in line.split()] for line in f]
	return mat

def loadMatrix(fname):
	with open(fname) as f:
		mat = [[float(x) for x in line.split()] for line in f]
	return mat

def sum(a,b):
	return [a[i]+b[i] for i in range(len(a))]
	
def prod(x,a):
	return [x*a[i] for i in range(len(a))]
	
def interpolate(rho,a,b):
	return sum(prod((1-rho), a),prod(rho,b))
	
if __name__ == "__main__":
	import sys
	import math
	
	if len(sys.argv) < 4:
		print 'USAGE: [global parameter file name] [local parameter file name] [rho] [output file name] *[logspace(true/false)]'
		sys.exit(0)
	
	global_param_file = sys.argv[1]		
	local_param_file = sys.argv[2]
	rho = float(sys.argv[3])
	output_file = sys.argv[4]
	logspace = sys.argv[5] in ('true')
	
	if logspace:
		gp = loadLogspaceMatrix(global_param_file) 	#global parameter in log space
		lp = loadLogspaceMatrix(local_param_file)	#local parameter in log space
	else:
		print 'loading not log space matrix!!!!!!'
		gp = loadMatrix(global_param_file) 	#global parameter
		lp = loadMatrix(local_param_file)	#local parameter
	
	out = open(output_file, 'w')
	for i in range(len(gp)):
		r = interpolate(rho, gp[i], lp[i])
		if logspace:
			r = map(lambda x: math.log(x), r)
		out.write('\t'.join(map(lambda f: '{:.4f}'.format(f),r)) + '\n')
		
