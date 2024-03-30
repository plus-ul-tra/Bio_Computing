import sys
import re

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

    
    if current_sequence:
      gene_list.append(current_sequence)
    else :
       # 아무 데이터도 없는 경우
       gene_list = 'no DNA sequence'


    input_file.close()    

    # input data에 등장하는 모든 comment와 sequence list, return
    return comment_list, gene_list


def finding_repeat_segment(dna_sequences):
    # dna_sequences = sequence list
    longest_len = 0
    longest_segment =''

    # R.E
    pattern = re.compile(r'(.+?)\1+')

    for sequence in dna_sequences:
        segment = pattern.findall(sequence)
        if segment:
            #segment 중  len가 가장 큰 값
            max_segment = max(segment, key =len)

        if longest_len < len(max_segment):
            longest_len = len(max_segment)
            longest_segment =(max_segment)
            

    


    return longest_segment

def output_data(filename, segment):
    file = open(filename,'w')
    if segment:
        s = 'longest repeated segment : ' + segment +'\n' +'len : ' + str(len(segment))
        file.write(s)
    else :
        s ='no repeat'    
        file.write(s)
    file.close()
    return 0 


def main(): 
    input_filename = 'test.txt'
    #input_filename = 'input.txt'
    output_filename = 'output.txt'
    #input_filename = sys.argv[1]    
    #output_filename = sys.argv[2]
    comment, dna_sequence = get_input_data(input_filename)
    
    longest_segment = finding_repeat_segment(dna_sequence)
    output_data(output_filename, longest_segment)
    return 0

if __name__ == '__main__':
    main()