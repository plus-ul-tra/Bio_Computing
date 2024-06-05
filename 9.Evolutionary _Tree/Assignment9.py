from collections import defaultdict
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
       # 2의 제곱수가 아닌 경우
    if (len(data_list) & (len(data_list) - 1)) != 0:
        code = -5
        return code, data_list
    
    first_length = len(data_list[0]) 

    for data in data_list:
        new_data = data.upper()

        if len(data) != first_length:
            code = -6
            return code, data_list
        
        if any(base not in 'ATGC' for base in new_data):
            code = -4

            return code, data_list
        else:
            new_data_list.append(new_data)
    
        
    # 아무 예외 없는 정상 데이터의 경우 code = 0
    return code, new_data_list


def fitch_algorithm(sequences):
    """Fitch 알고리즘을 적용하여 내부 노드에 대한 최적의 서열을 찾습니다."""
    n = len(sequences)
    m = len(sequences[0])
    tree = defaultdict(list)
    
    # 노드 초기화
    nodes = {i: set(seq) for i, seq in enumerate(sequences)}
    parents = {i: None for i in range(2 * n - 1)}
    
    # 초기 트리 구성 (리프 노드만)
    for i in range(n):
        tree[i] = []
    
    current_node = n
    while len(nodes) > 1:
        # 두 노드를 선택
        (a, aset), (b, bset) = sorted(nodes.items())[:2]
        del nodes[a]
        del nodes[b]
        
        # 새로운 노드 생성
        intersection = aset & bset
        if intersection:
            nodes[current_node] = intersection
        else:
            nodes[current_node] = aset | bset
        
        tree[current_node] = [a, b]
        parents[a] = current_node
        parents[b] = current_node
        current_node += 1

    root = current_node - 1
    
    # 최적 서열 결정
    def backtrack(node):
        if node < n:
            return list(sequences[node])
        left, right = tree[node]
        left_seq = backtrack(left)
        right_seq = backtrack(right)
        optimal_seq = []
        for l, r in zip(left_seq, right_seq):
            if l == r:
                optimal_seq.append(l)
            else:
                optimal_seq.append(min(l, r))  # 임의의 선택
        return optimal_seq
    
    optimal_sequences = []
    for node in range(2 * n - 1):
        if node < n:
            optimal_sequences.append(sequences[node])
        else:
            optimal_sequences.append(''.join(backtrack(node)))

    return optimal_sequences


def calculate_parsimony(sequences, tree):
    score = 0
    for node, children in tree.items():
        if children:
            left, right = children
            for i in range(len(sequences[0])):
                if sequences[node][i] != sequences[left][i]:
                    score += 1
                if sequences[node][i] != sequences[right][i]:
                    score += 1
    return score

def _tree(sequences, optimal_sequences):
    tree = defaultdict(list)
    for i in range(len(optimal_sequences) - 1):
        tree[i] = []

    n = len(sequences)
    current_node = n
    while current_node < 2 * n - 1:
        tree[current_node] = [current_node - n, current_node - n + 1]
        current_node += 1

    return tree

def result_output(filename, parsimony_score, optimal_sequences):
    output_file = open(filename, 'w')
    s = 'parsimony_score : ' + str(parsimony_score) + '\nsequences : ' + ", ".join(map(str,optimal_sequences))
    print(s)
    output_file.write(s)
    output_file.close()

    return 0

def exception_output(filename, code):
    output_file = open(filename, 'w')
    # -1 : no fasta format
    # -2 : no data in file
    # -3 : single sequence
    # -4 : no dna sequence
    # -5 : no 2^n sequences
    # -6 : Incorrect length of sequences
    if code == -1:
        s = 'No correct format'
    elif code == -2:
        s = 'No DNA sequence'
    elif code == -3:    
        s = 'Need more sequences'
    elif code == -4:
        s = 'No DNA sequence'
    elif code == -5:
        s = 'Incorrect number of sequences'
    elif code == -6:
        s = 'Incorrect length of sequences'
    output_file.write(s)
    output_file.close()
    return 0

def main():
    #input_filename ='input.txt'
    #output_filename ='output.txt'
    input_filename = sys.argv[1]    
    output_filename = sys.argv[2]

    if not os.path.exists(input_filename):
        print("No input file")
        return 0
    

    ex_code, sequences = get_input_data(input_filename)

    # 예외 판별 및 출력
    if ex_code != 0:
        exception_output(output_filename, ex_code)
        return 0
    
    start_time= time.time()
    optimal_sequences = fitch_algorithm(sequences)
    evo_tree = _tree(sequences, optimal_sequences)
    parsimony_score = calculate_parsimony(optimal_sequences, evo_tree)
    end_time = time.time()
    elapsed_time = (end_time - start_time) * 1000000
    print(elapsed_time)
    #정상 출력
    result_output(output_filename, parsimony_score, optimal_sequences)

    return 0


if __name__ == '__main__':
    main()