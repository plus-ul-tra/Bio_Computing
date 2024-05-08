msa =list()

s1 = "012345678901234567890123456789012345678901234567890123456789"
s2 = "------------------------------------------------------------"
s3 = "dickdickdickdickdickdickdickdickdickdickdickdickdickdickdick"
s4 = "gyumgyumgyumgyumgyumgyumgyumgyumgyumgyumgyumgyumgyumgyumgyum"
st = "************-----------------*******************************"
msa.append(s1)
msa.append(s2)
msa.append(s3)
msa.append(s4)

def output_data(filename, msa, star):
    output_file =open(filename, 'w')
    s =""
    
    if len(msa[0]) > 50 :
        msa.append(star)
        chunk_num = len(msa[0]) // 50
        # 180 의 경우 몫 3
        for i in range(0,chunk_num):
            for line in msa:
                # 한줄씩 가져와서
                s+=line[i*50:(i+1)*50] +'\n'



    output_file.write(s)
    output_file.close()

    return 0

output = 'testoutput.txt'
output_data(output,msa,st)