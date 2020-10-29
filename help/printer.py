import sys

def printnow(s, newline = True):
	s = str(s)
	if newline:
		s += '\n'
	sys.stdout.write(s)
	sys.stdout.flush()

def printerr(s, newline = True):
	s = str(s)
	if newline:
		s += '\n'
	sys.stderr.write(s)
	sys.stderr.flush()
