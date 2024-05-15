"""
2019253022 홍한울
"""

import sys
import re
import time


def get_input_data(filename):

    input_file = open(filename,'r')
    comment_list = list()
    gene_list = list()
    current_sequence = ''
    
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


    input_file.close()    

    # input data에 등장하는 모든 comment와 sequence list, return
    return comment_list, gene_list

    

def finding_repeat_segment(dna_sequences):
    # dna_sequences = sequence list
    longest_len = 0
    longest_segment =''
    sequence = dna_sequences[0]

    # R.E
   
    pattern = re.compile(r'(.+)(?:\1)+')
    if dna_sequences == 'No DNA sequence':
        longest_segment ='No DNA sequence'
        
    else:
            # DNA 서열이 아닌 문자가 sequence에 포함되는 경우
            if any(base not in "ATCGatcg" for base in sequence):
                longest_segment = 'No correct format'

            segment = pattern.findall(sequence)
            if segment:
                #segment 중  len가 가장 큰 값
                max_segment = max(segment, key =len)

            if longest_len < len(max_segment):
                longest_len = len(max_segment)
                longest_segment = max_segment
                
    return longest_segment

def output_data(filename, segment):
    file = open(filename,'w')

    # longest repeated segment가 2이상인 경우만 출력

    # 빈 파일의 경우 (DNA없는 경우)
    if 'No DNA sequence' == segment:
        s = segment + ' (empty file)'
        file.write(s)
    elif 'No correct format' == segment:
        s = segment
        file.write(s)
    else:
        if len(segment) > 1:
            s = 'longest repeated segment : ' + segment +'\n' +'len : ' + str(len(segment))
            file.write(s)

        # 반복되는 항목이 없는 경우
        elif len(segment) == 1 :
            s ='No repeats'    
            file.write(s)
    
    file.close()
    return 0 




def main():
    #input_filename = 'input.txt'
    #output_filename = 'output.txt' 
    input_filename = sys.argv[1]    
    output_filename = sys.argv[2]

    comment, dna_sequence = get_input_data(input_filename)

    start_time= time.time()
    longest_segment = finding_repeat_segment(dna_sequence)
    end_time = time. time()

    elapsed_time = end_time - start_time
    output_data(output_filename, longest_segment)
    print(f"Elapsed Time: {elapsed_time:f} microseconds")


    return 0

if __name__ == '__main__':
    main()