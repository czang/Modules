import bisect

def is_list_sorted(List):
        """
        Check if sorted in ascending order.
        input is a list of values.
        output: sorted =1 or 0
        """
        sorted = 1;
        for index in range(0, len(List)-1):
                if List[index] > List[index+1]:
                        sorted = 0;
        return sorted;

def detect_dip_from_list(List, t1 = 0.1):
	d1 = derivative(List)
	d2 = derivative(d1)
	index_result = []
	for i in range(2, len(List)-2):
		i1 = i - 1
		i2 = i1 - 1
		m = float(List[i] + List[i1] + List[i2]) / 3
		if d1[i1] <= t1 * m and d2[i2] > 0:
			index_result.append(i)
	return index_result

def detect_peak_from_list(List, t1 = 0.1):
	d1 = derivative(List)
	d2 = derivative(d1)
	index_result = []
	for i in range(2, len(List)-2):
		i1 = i - 1
		i2 = i1 - 1
		m = float(List[i] + List[i1] + List[i2]) / 3
		if d1[i1] <= t1 * m and d2[i2] < 0:
			index_result.append(i)
	return index_result

def detect_dips(List, dip_index_list, peak_index_list, depth=0.5):
	

def get_sublist_index(List, start, end):
	if is_list_sorted(List) == 0:
		List.sort()
	i = bisect_left(List, start)
	j = bisect_right(List, end)
	return i, j

def wiglists2datalist(List1,List2, step=10):
	assert len(List1) == len(List2)
	result = []
	result.append(List2[0])
	flag = List1[0]
	i = 1
	while i < len(List1):
		if List1[i] == flag + step: 
			result.append(List2[i])
			flag = List1[i]
			i += 1
		else:
			while List1[i] >= flag + step:
				result.append(0)
				flag += step
			i += 1
	assert len(result) * (step-1) == List1[-1] - List[1]
	return result

def derivative(List):
	result = []
	for i in range(1,len(List)-1):
		result.append(float(List[i+1] - List[i-1])/2)
	assert len(result) + 2 == len(List)
	return result


def find_dips(wig_x, wig_y, island_start_list, island_end_list, valley_depth=0.5):
	

def peak_5(List):
	assert len(List) == 5
	if List[2] > List[1] and List[2] > List[3]:
		if List[1] > List[0] and List[3] > List[4]:
			return 1
		else:
			return -1
	else:
		return -1

def dip_5(List):
	assert len(List) == 5
	if List[2] < List[1] and List[2] < List[3]:
		if List[1] < List[0] and List[3] < List[4]:
			return 1
		else:
			return -1
	else:
		return -1
