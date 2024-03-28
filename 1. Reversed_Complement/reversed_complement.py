"""
2019253022
홍한울
"""
import sys 

def get_input_data(filename):
    input_file = open(filename, 'r')
    comment_list = list()
    # gene
    gene_list = list()
    
    for line in input_file:

     if '>' in line:
        comment_list.append(line.strip())
     
     else :
        gene_list.append(line.strip())

    input_file.close()    
    return comment_list,gene_list
     

def reverse_complement(sequence):
    reversed_complement=list()
    # dna 문자열 List

    # dictionary key, value로 변환할 dna sequence 맵핑
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 'T', 't': 'A', 'c': 'G', 'g': 'C'}
    
    # 전체 sequence에서 한 sequence만 
    for line in sequence:
       # 빈문자열 초기화
       reversed_sequence = ''

       # 문자열을 역방향으로 순차적으로 하나의 dna를 상보적 관계 매칭
       for dna in reversed(line): 
             
             # if 염기서열이 아닌 것이 있는 경우
             if dna not in complement:
                 reversed_sequence = 'No DNA sequence'
                 break
         
             reversed_sequence += complement[dna]
          # dictionary의 key값으로 상보관계의 value 입력
     
       reversed_complement.append(reversed_sequence)
        
    return reversed_complement

def output_to_file(filename, comment,result):

    file =open(filename,'w')
    for i, com in enumerate(comment):
       
       #최종 출력 : > 이후의 comment와 해당 comment의 DNA 서열의 reversed complement 출력

       s = com[1:] + " : " + result[i] + "\n"
       file.write(s)

    file.close()
    return 0

def main():
    #input_filename = 'input.txt'
    input_filename = sys.argv[1]
    #output_filename = 'output.txt'    
    output_filename = sys.argv[2]

    comment, dna_sequence = get_input_data(input_filename)
    # algorithm 동작
    result = reverse_complement(dna_sequence)

    output_to_file(output_filename,comment,result)

if __name__ == '__main__':
    main()