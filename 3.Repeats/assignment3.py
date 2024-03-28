import sys

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


def finding_repeat_segment(dna_sequence):


    return 0




def main(): 
    input_filename = 'input.txt'
    output_filename = 'output.txt'
    #input_filename = sys.argv[1]    
    #output_filename = sys.argv[2]
    comment, dna_sequence = get_input_data(input_filename)
    
    finding_repeat_segment(dna_sequence)

    return 0

if __name__ == '__main__':
    main()