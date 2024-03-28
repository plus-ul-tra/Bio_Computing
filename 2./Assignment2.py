import sys

def get_input_data(filename):
    input_file = open(filename , 'r')
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


def counting_pattern(comment, sequence):
    max_count = 0
    highest_sequence = ''
    highest_comment = ''
    max_rate = 0
    dna_set = ('A', 'G', 'T', 'C', 'a', 'g', 't', 'c')

    for i, line in enumerate(sequence):
        g_count = 0
        c_count = 0
        f_count = 0
        for dna in line:
            
            # sequence 에 DNA가 아닌 값이 하나라도 있는 경우 바로 비정상 값 return
            if dna not in dna_set:
                highest_comment = 'no DNA sequence'
                highest_sequence = 'no DNA sequence'
                max_rate = 0
                break
            
            f_count += 1
            if dna == 'G' or dna == 'g':
                g_count += 1
            elif dna == 'C' or dna == 'c':
                c_count += 1
        else:
            # 정상적인 데이터의 경우 count를 통한 rate 반환
            if max_count < g_count + c_count:
                max_count = g_count + c_count
                max_rate = (g_count + c_count) / f_count
                highest_comment = comment[i]
                highest_sequence = line

    return highest_comment, highest_sequence, max_rate


def output_to_file(filename, comment, result, rate):
   file =open(filename,'w')
   
   # 비정상 데이터의 경우 최종 출력 
   if result =='no DNA sequence' :
       file.write(result)
   else:
    # 정상 데이터의 경우
    s = comment[1:] + ', G, C rate : ' + str(rate)
    file.write(s)

   file.close()
   return 0


def main():
    #input_filename = 'input.txt'
    #output_filename = 'output.txt'
    input_filename = sys.argv[1]    
    output_filename = sys.argv[2]

    comment, dna_sequence = get_input_data(input_filename)
    # algorithm 동작
    highest_comment, result, rate = counting_pattern(comment, dna_sequence)

    output_to_file(output_filename,highest_comment,result, rate)

if __name__ == '__main__':
    main()