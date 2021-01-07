#include <iostream>
#include <utility>
#include <cassert>

// Let us define as a C-almost-ascending subsequence (Caas) such a subsequence
// of a given sequence, that each of its elements is greater than all the previous
// elements diminished by C. The task here is to find the longest common Caas
// of two sequences, or, if its length is proven to be more than 20, answer
// with 20. Due to the assumption that no sequence longer than 20 has to be memorized,
// this can be done in O(mn) time, where n, m are the lengths of the sequences.

// My solution is a combination of the classical LAS algorithm and the dynamic
// algorithm for finding the longest common subsequence of two sequences.

using sequence = std::pair<int, int*>;

// How many cells are occupied in a given sequence
int occupied_in(const sequence& seq){
    return std::get<0>(seq);
}

int* local_array(const sequence& seq){
    return std::get<1>(seq);
}

int MAX_LENGTH(){
    static int res = 20;
    return res;
}

int* create_1D_array(){
    int max = MAX_LENGTH();
    int* res = (int*) malloc(max*sizeof(int));
    for(int i = 0; i < max; i++) res[i] = 0;
    return res;
}

sequence create_new_sequence(){
    sequence res;
    std::get<0>(res) = 0;
    std::get<1>(res) = create_1D_array();
    return res;
}

sequence create_from_array(int* arr, int size){
    sequence res;
    std::get<0>(res) = size;
    std::get<1>(res) = arr;
    return res;
}

int max(int a, int b){
    return a > b ? a : b;
}

int min(int a, int b){
    return a < b ? a : b;
}

void print_seq(sequence seq){
    int size = occupied_in(seq);
    int* arr = local_array(seq);
    for(int i = 0; i < size; i++) std::cout << arr[i] << " ";
    std::cout << "\n";
}

sequence merge_sequences(sequence seqx, sequence seqy){

    int* res_arr = (int*) malloc(MAX_LENGTH()*sizeof(int));
    int lenx = occupied_in(seqx);
    int lenny = occupied_in(seqy);
    int* arrx = local_array(seqx);
    int* arry = local_array(seqy);
    int a = min(lenx, lenny);
    int k;
    for(int i = 0; i < a; i++){
        k = min(arrx[i], arry[i]);
        res_arr[i] = k;
    }
    if(lenx < lenny){
        for(int i = a; i < lenny; i++) res_arr[i] = arry[i];
    }
    if(lenny < lenx){
        for(int i = a; i < lenx; i++) res_arr[i] = arrx[i];
    }
    sequence res = create_from_array(res_arr, max(lenx, lenny));
    return res;
}

bool within_downward_tolerance(int max, int n, int c){
    return n + c >= max && n <= max;
}

sequence copy_and_insert(sequence seq, int n, int c, bool* foundMAX){
    int size = occupied_in(seq);
    int* prev_arr = local_array(seq);
    int* res_arr = (int*) malloc (MAX_LENGTH() * sizeof(int));
    int res_size = size;
    if(size == 0){
        res_arr[0] = n;
        res_size = 1;
        return create_from_array(res_arr, res_size);
    }
    if(n >= prev_arr[size-1]-c){
        ++res_size;
        if(res_size >= 20){
            *foundMAX = true;
            return create_from_array(prev_arr, res_size);
        }
    }
    int i = size;
    while(i > 0){
        if(n >= prev_arr[i-1]){
            if(i == size || prev_arr[i] > n){
                res_arr[i] = n;
            }
            i--;
            break;
        }
        if(within_downward_tolerance(prev_arr[i-1], n, c)){
            res_arr[i] = prev_arr[i-1];
            i--;
         }
        else{
            if(i < size) res_arr[i] = prev_arr[i];
            i--;
        }
    }
    while(i > 0){
        res_arr[i] = prev_arr[i];
        i--;
    }
    if(n <= prev_arr[0]) res_arr[0] = n;
    else res_arr[0] = prev_arr[0];
    sequence res = create_from_array(res_arr, res_size);
    return res;
}

sequence* create_1D_seq_array(int m){
    sequence* res = (sequence*) malloc(m*sizeof(sequence));
    return res;
}

sequence** create_2D_seq_array(int n, int m){
    sequence** res = (sequence**) malloc(n*sizeof(sequence*));
    for(int i = 0; i < n; i++) res[i] = create_1D_seq_array(m);
    return res;
}

// Maximal common almost-ascending sequence
int max_caas(int n, int m, int c, int* A, int* B){
    sequence** seqs = create_2D_seq_array(n+1, m+1);
    for(int i = 0; i <= n; i++) seqs[i][0] = create_new_sequence();
    for(int j = 1; j <= m; j++) seqs[0][j] = create_new_sequence();
    bool foundMAX = false;
    for(int i = 1; i <= n; i++){
        for(int j = 1; j <= m; j++){
            if(A[i-1] == B[j-1]){
                seqs[i][j] = copy_and_insert(seqs[i-1][j-1], A[i-1], c, &foundMAX);
                if(foundMAX){
                    return MAX_LENGTH();
                }
            }
            else{
                seqs[i][j] = merge_sequences(seqs[i-1][j], seqs[i][j-1]);
            }
        }
    }
    return occupied_in(seqs[n][m]);
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    int n;
    int m;
    int c;
    std::cin >> n;
    std::cin >> m;
    std::cin >> c;
    int* A = (int*) malloc(n*sizeof(int));
    int* B = (int*) malloc(m * sizeof(int));
    for(int i = 0; i < n; i++) std::cin >> A[i];
    for(int i = 0; i < m; i++) std::cin >> B[i];
    std::cout <<  max_caas(n, m, c, A, B);
    return 0;
}
