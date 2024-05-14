import numpy as np
import time
import sys
import os

def get_input_data(filename):
    input_file = open(filename,'r')
    comment_list = list()
    gene_list = list()
    current_sequence = ''
    
    first_char = input_file.read(1)
    input_file.seek(0)

    if first_char == '>':
    # \n를 고려하여 다음 comment, 즉 '>' 발견 전까지 모두 dna sequence로 설정.
        for line in input_file:
            line = line.strip()

            if line.startswith('>'):
                comment_list.append(line)

                if current_sequence:
                    gene_list.append(current_sequence)

                current_sequence = ''

            else:
                    current_sequence += line
    elif first_char =="":
        gene_list.append('nodata')
        code = -2
        return code, gene_list
    
    else:
        # FASTA format이 아닌 경우(#########했음#########)
        gene_list.append('nofasta')
        code = -1
        input_file.close()

        return code, gene_list

    if current_sequence:
      gene_list.append(current_sequence)
    else :
       # 아무 데이터도 없는 경우
       gene_list.append('nodata')
       code = -1
       return code, gene_list
       


    input_file.close()    

    # 예외코드 -3, -4는 해당 함수에서 선별
    exception_code, final_data = exception_check(gene_list)

    # input data에 등장하는 모든 comment와 sequence list, return
    return exception_code, final_data




def exception_check(data_list):
    # 여기까지 넘어 왔다. fasta format을 만족하는 데이터 들이 있긴하다
    # 동작 : list 개수, upper 로 변환,
    new_data_list = list()
    code = 0
    if len(data_list) == 1:
        # data가 하나 뿐이다 최소 2개는 있어야함
        code = -3
        return code, data_list
    
    for data in data_list:
        new_data = data.upper()
        if any(base not in 'ATGC' for base in new_data):
            code = -4

            return code, data_list
        else:
            new_data_list.append(new_data)
    
        
    # 아무 예외 없는 정상 데이터의 경우 code = 0
    return code, new_data_list

# 정상 알고리즘 출력
def result_output(filename, data):
    output_file = open(filename, 'w')

    output_file.close()
    return 0



def exception_output(filename, code):
    output_file = open(filename, 'w')
    # -1 : no fasta format
    # -2 : no data in file
    # -3 : single sequence
    # -4 : no dna sequence
    if code == -1:
        s = 'No correct format'
    elif code == -2:
        s = 'No DNA sequence'
    elif code == -3:    
        s = 'Need more sequences'
    elif code == -4:
        s = 'No DNA sequence'

    output_file.write(s)
    output_file.close()    

    return 0


def pattern_finding(sequences):



    return 0



def main():
    input_filename ='input.txt'
    output_filename ='output.txt'
    #input_filename = sys.argv[1]    
    #output_filename = sys.argv[2]
    #pattern_length = int.(sys.argv[3])
    if not os.path.exists(input_filename):
        print("No input file")
        return 0
    
    ex_code, dna_sequence = get_input_data(input_filename)
    
    print(ex_code)
    if ex_code != 0:
        exception_output(output_filename,ex_code)
        return 0
    #start_time= time.time()
    result = pattern_finding(dna_sequence)
    #end_time = time.time()
    #elapsed_time = (end_time - start_time) * 1000000
    #print(f"Elapsed Time: {elapsed_time:f} microseconds")
    result_output(output_filename, result)
    
    return 0




if __name__ == '__main__':
    main()