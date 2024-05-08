import numpy as np
import time
import sys
import os

def check_error(data):
    modified_data =list()

    if len(data) == 1:
        if data == "nofasta":
            return -2
        elif data == "no protein":
            return -1
        else :
            # 정상이나 시퀸스가 하나인 경우
            return -3
    
    
    for line in data:
        new_line=line.upper()
        modified_data.append(new_line)
        if any(base not in "CSTPAGNDEQHRKMILVFYW" for base in new_line):
            return -1

    return modified_data

def get_blosum62_matrix(filename):
    blosum62 = open(filename, 'r')
    amino_list = list()
    score_matrix = list()
    lines = blosum62.readlines()[1:]
    for line in lines:
        row=list(line[1:].split())
        score_matrix.append(row)

    return score_matrix



def get_input_data(filename):
    input_file = open(filename,'r')
    comment_list = list()
    protein_list = list()
    current_sequence = ''
    
    first_char = input_file.read(1)
    input_file.seek(0)

    if first_char =='>':
    # \n를 고려하여 다음 comment, 즉 '>' 발견 전까지 모두 dna sequence로 설정.
        for line in input_file:
            line = line.strip()

            if line.startswith('>'):
                comment_list.append(line)

                if current_sequence:
                    protein_list.append(current_sequence)

                current_sequence = ''

            else:
                    current_sequence += line
    else:
        # FASTA format이 아닌 경우
        protein_list.append('nofasta') 
        return protein_list

    if current_sequence:
      protein_list.append(current_sequence)
    else :
       # 아무 데이터도 없는 경우
       protein_list = 'no protein'
       clear_data = -4
       # -4 리턴
       return clear_data
    
    # 에러 코드 -1 ~ -3 확인
    clear_data = check_error(protein_list)


    input_file.close()    

    # input data에 등장하는 모든 comment와 sequence list, return
    return clear_data


# pairwise alignment 쌍에 대한 blosum62를 이용하여 scoring
def get_blosum_score(amino1, amino2,blosum62_score_matrix):
    amino_acid = "CSTPAGNDEQHRKMILVFYW"
    index1 = amino_acid.find(amino1)
    index2 = amino_acid.find(amino2)
    score = int(blosum62_score_matrix[index1][index2])
    
    return score


def pairwise_sequence(seq1, seq2, _blosum62_score_matrix, gap_penalty = -5):
    n = len(seq1)
    m = len(seq2)

    # score_matrix 초기화
    score_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # back trace 시 작성할 score_matrix 생성
    traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # 초기 gap penalty 설정
    for i in range(1, n + 1):
        score_matrix[i][0] = i * gap_penalty
        traceback_matrix[i][0] = 1
    for j in range(1, m + 1):
        score_matrix[0][j] = j * gap_penalty
        traceback_matrix[0][j] = 2

    # 각 자리간의 경우로 score_matrix 채우기
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # match 된 경우에서 뒤에 더해지는 score를 blosum62의 table을 매김.
            # 이부분만 수정 -> seq1[i-1], seq2[j-1]을 받아서 blosum62의 점수 반환하는 함수사용
            value = get_blosum_score(seq1[i - 1],seq2[j - 1],_blosum62_score_matrix)
            # (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            match = score_matrix[i - 1][j - 1] + value
            # delete, insert는 모두 gap (옆, 밑 이동을 의미)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
            if score_matrix[i][j] == match:
                # Diagonal
                traceback_matrix[i][j] = 0  

            elif score_matrix[i][j] == delete:
                # Up
                traceback_matrix[i][j] = 1  
            else:
                # Left
                traceback_matrix[i][j] = 2  

    # 최적경로 back trace
    aligned_seq1 = ''
    aligned_seq2 = ''
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback_matrix[i][j] == 0:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and traceback_matrix[i][j] == 1:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1
    return score_matrix[n][m], aligned_seq1, aligned_seq2

def global_sequence_alignment(seq1, seq2, match_score=2, mismatch_score=-2, gap_penalty=-1):
    n = len(seq1)
    m = len(seq2)

    # matrix 초기화
    score_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # back trace 시 작성할 matrix 생성
    traceback_matrix = np.zeros((n + 1, m + 1), dtype=int)

    # 초기 gap penalty 설정
    for i in range(1, n + 1):
        score_matrix[i][0] = i * gap_penalty
        traceback_matrix[i][0] = 1
    for j in range(1, m + 1):
        score_matrix[0][j] = j * gap_penalty
        traceback_matrix[0][j] = 2

    # 각 자리간의 경우로 matrix 채우기
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
            if score_matrix[i][j] == match:
                # Diagonal
                traceback_matrix[i][j] = 0  

            elif score_matrix[i][j] == delete:
                # Up
                traceback_matrix[i][j] = 1  
            else:
                # Left
                traceback_matrix[i][j] = 2  

    # 최적경로 back trace
    aligned_seq1 = ''
    aligned_seq2 = ''
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback_matrix[i][j] == 0:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and traceback_matrix[i][j] == 1:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2

def select_center_sequence(sequence_list, score_matrix):

    max_strength = 0
    for s1 in sequence_list:
        strength = 0
        for s2 in sequence_list:
            if s1 != s2:
                #본인을 제외한 나머지 sequence와 pairwise 구함.
                score, align1,align2 = pairwise_sequence(s1,s2, score_matrix) 
                strength += score
        
        if strength > max_strength :
            max_strength = strength
            center_sequence = s1


    return center_sequence

# 최종 MSA list 생성.
def multiple_sequence_alignment(protein_list, blosum62_score_matrix):
    center_list = list()
    non_center_list = list()
    msa =list()
    center = select_center_sequence(protein_list, blosum62_score_matrix)
    for s1 in protein_list:
        if center != s1:
           # center가 아닌 다른 녀석들과의 aligned 진행
           score, aligned_center, aligned_non_center = pairwise_sequence(center, s1, blosum62_score_matrix)
           # center_list의 첫번째 value = 처음의 psa 0 : 1, 2, 3
           center_list.append(aligned_center)
           # non_center_list center가 아닌 나머지 alignment 0, 1, 2, 3
           non_center_list.append(aligned_non_center)

    msa.append(center_list[0])
    for s1 in non_center_list:
        # x(0): center , y(1) ,z(2) ,w(3), k(4)
        msa.append(s1)
        
    for unmodified_center in center_list[1:] :
        # 둘의 결과는 같음 -> 같은 문자열에서 갭의 차이만 있었기 때문에
        new_center1, new_center2 = global_sequence_alignment(msa[0],unmodified_center)
        # 기준변경
        msa[0] = new_center1
        for index, s1 in enumerate(msa[1:],start =1):
            same_center, modified = global_sequence_alignment(msa[0],s1)
            msa[index] = modified
            msa[0] = same_center
        # 새로운 center가 등장할 때 마다. 나머지를 update, 신규 append
            # 기준을 바꾸는게 우선, 그 후 update,
    len_str = len(msa[0])
    for i in range(1, len(msa)):  # 첫번째 요소를 제외한 나머지 요소들에 대해
        if len(msa[i]) < len_str:  # len보다 작으면
             msa[i] += '-' * (len_str - len(msa[i]))  # 문자열 뒤에 '-' 추가
    

    return msa

# 같은 column * 문자열 생성
def searching_same_char(MSA_list):

    strings = MSA_list
    
    length = len(strings[0])
    star_list=list()
    
    for i in range(length):
        
        characters = [s[i] for s in strings]
        
        if all(char == characters[0] and char !='-' for char in characters):
            star_list.append("*")   
        else:
            star_list.append(" ")

    result = ''.join(star_list)

    return result

def output_data(filename, msa, star):
    output_file =open(filename, 'w')
    s =""
    if len(msa[0]) > 50 :
            msa.append(star)
            chunk_num = len(msa[0]) // 5
        # 180 의 경우 몫 3
            for i in range(0,chunk_num):
                for line in msa:
                # 한줄씩 가져와서
                    s+=line[i*50:(i+1)*50] +'\n'
    else:
        for line in msa:
            s += line +'\n'
        s += star

    output_file.write(s)
    output_file.close()

    return 0

def output_error(filename, error_code):
    output = open(filename, 'w')

    if error_code == -1:
        s = 'there is no protein acid'
    if error_code == -2:
        s = 'No protein sequence'
    if error_code == -3:
        s = "Need more sequence"
    if error_code == -4:
        s= "No FASTA format"

    output.write(s)
    output.close()
    return 0


def main():
    #input_filename ='input.txt'
    #input_filename ='er1.txt'
    #output_filename ='output.txt'
    blosum62_file = 'blosum62.txt'
    input_filename = sys.argv[1]    
    output_filename = sys.argv[2]

    #blosum62 생성
    score_matrix = get_blosum62_matrix(blosum62_file)
    protein_list = get_input_data(input_filename)
    
    if protein_list == -1:
        output_error(output_filename, protein_list)
        # -1의 경우, 프로틴이 아닌 다른 문자
        return 0
    elif protein_list == -2:
        output_error(output_filename, protein_list)
        # 데이터 없는 경우
        return 0
    elif protein_list == -3:
        output_error(output_filename, protein_list)
    elif protein_list == -4 :
        output_error(output_filename, protein_list)

    msa = multiple_sequence_alignment(protein_list, score_matrix)
    star_list = searching_same_char(msa)
    output_data(output_filename,msa, star_list)
    

    return 0

if __name__ == '__main__':
    main()