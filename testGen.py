
def create_test_case(sample_len, virus_len):
    file = open("test1.in", 'w') 
    file.write("%d\n" % (1))

    s = ""
    for i in range(sample_len - 1):
        if (i < sample_len // 2):
            s += 'A'
        else:
            s += 'A'
    
    
    v = ""
    for i in range(virus_len):
        v += 'A'

    file.write("%s %s" % (s, v))

def create_test_phread(sample_len):
    file = open("test1.in", 'w') 
    file.write("%d\n" % (1))

    s = ""
    for i in range(sample_len):
        s += "B"
    
    file.write("%s %d %d" % (s, sample_len // 2, 0))

#create_test_case(200001, 100000)
create_test_phread(200000)