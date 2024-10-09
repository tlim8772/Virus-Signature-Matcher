
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

create_test_case(200001, 100000)