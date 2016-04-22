# Convert string to lower case up to the first occurrence of a separator
def lowerToSep(string, separator='='):
	line=string.partition(separator)
	string=str(line[0]).lower()+str(line[1])+str(line[2])
	return string