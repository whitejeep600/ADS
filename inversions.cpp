#include <iostream>
#include <tuple>
#include <cassert>

// For a permutation a of numbers 1, 2... n,
// let us define as k-inversion such a sequence of indexes
// x1, x2, ... xk that x1 < x2 < ... < xk, but a(x1) > a(x2) > ... > a(xk) - for
// example, 2-inversions are 'ordinary' inversions of a permutation. The task is to 
// find the number of k-inversions in a given permutation.

// To do this efficiently, I used segment trees. 

using seg_tree = std::tuple< long long,  long long, long long*,  long long*>;

// for convenience, nodes are numbered 1-n (!)
 long long right_son( long long* values,  long long node_no,  long long total_nodes){
     long long res = 2*node_no+1;
    if(res > total_nodes) return -1;
    if(values[res-1] == -1) return -1;
    return res;
}

 long long left_son( long long* values,  long long node_no,  long long total_nodes){
     long long res = 2*node_no;
    if(res > total_nodes) return -1;
    if(values[res-1] == -1) return -1;
    return res;
}

// assuming that a given node does have a father
 long long father( long long node_no){
    return node_no/2;
}

 long long node_no_of_ith_leaf(seg_tree tree,  long long k){
    static  long long beginning = std::get<0>(tree) - std::get<1>(tree);
    return beginning + k;
}

 long long ceil( long long a,  long long b){
    if(a%b == 0) return a / b;
    else return a/b + 1;
}

 long long aux_nodes_for_n_leaves( long long n){
     long long res = 1;
    do{
        res *= 2;
    }
    while(res < n);

    return res-1+n;
}

 long long nodes_for_n_leaves( long long n){
    static  long long res = aux_nodes_for_n_leaves(n);
    return res;
}

 long long local_value(seg_tree &tree,  long long node_no){
    return std::get<2>(tree)[node_no-1];
}

 long long right_bound_for_node(seg_tree& tree,  long long node_no){
    return std::get<3>(tree)[node_no-1];
}

int max(long long a, long long b){
     return a > b ? a : b;
 }

 long long fill_right_bound_for_node( long long* values,  long long node_id,  long long* nodes,  long long total_nodes){
    if(values[node_id-1] == -1) return -1;
    if(nodes[node_id-1] != -1) return nodes[node_id-1];
    if(right_son(values, node_id, total_nodes) != -1) {
         long long res = fill_right_bound_for_node(values, right_son(values, node_id, total_nodes), nodes, total_nodes);
        fill_right_bound_for_node(values, left_son(values, node_id, total_nodes)
                , nodes, total_nodes);
        nodes[node_id - 1] = res;
        return res;
    }
    else if(left_son(values, node_id, total_nodes) != -1){
         long long res = fill_right_bound_for_node(values, left_son(values, node_id, total_nodes), nodes, total_nodes);
        nodes[node_id - 1] = res;
        return res;
    }
    else{
        assert(false);
        return 0;
    }
}

void fill_right_bounds( long long* values,  long long* nodes,  long long total_nodes,  long long leaves){
    for( long long i = 0; i < leaves; i++) nodes[i + total_nodes - leaves ] = i+1;
    fill_right_bound_for_node(values, 1, nodes, total_nodes);
}

void fill_parentage( long long* values,  long long start){
    while(start != 0){
        if(values[start-1] == 0) return;
        values[start-1] = 0;
        start /= 2;
    }
    values[start] = 0;
}

seg_tree init_tree_n_leaves( long long n){
     long long total_nodes = nodes_for_n_leaves(n);
    seg_tree result;
    std::get<0>(result) = total_nodes;
    std::get<1>(result) = n;

     long long* values = ( long long*) malloc(total_nodes*sizeof( long long));
    for( long long i = 0; i < total_nodes; i++) values[i] = -1;
    for( long long i = 0; i < n; i++){
        fill_parentage(values, i+total_nodes-n+1);
    }
    std::get<2>(result) = values;

     long long* right_bounds = ( long long*) malloc(total_nodes*sizeof( long long));
    for( long long i = 0; i < total_nodes; i++) right_bounds[i] = -1;
    fill_right_bounds(values, right_bounds, total_nodes, n);
    std::get<3>(result) = right_bounds;

    return result;
}

void clear_tree(seg_tree& tree){
    for( long long i = 0; i < std::get<0>(tree); i++){
        if(std::get<2>(tree)[i] != -1) std::get<2>(tree)[i] = 0;
    }
}

 long long left_sons_right_bound(seg_tree& tree,  long long node_no){
    return right_bound_for_node(tree, left_son(std::get<2>(tree), node_no, std::get<0>(tree)));
}

 long long max_of_segment_for_node(seg_tree& tree,  long long left,  long long right,  long long node_no,  long long local_left,  long long local_right){
    if(left > right) return 0;
    static  long long total_nodes = std::get<0>(tree);
    if(node_no == -1 || node_no > total_nodes) return 0;
    if(local_left > right || local_right < left) return 0;
    if(left <= local_left && right >= local_right){
        return local_value(tree, node_no);
    }
    else {
         long long left_son_no = left_son(std::get<2>(tree), node_no, total_nodes);
         long long lsrb = right_bound_for_node(tree, left_son_no);
         long long right_son_no = right_son(std::get<2>(tree), node_no, total_nodes);
         long long from_left = max_of_segment_for_node(tree, left, right, left_son_no, local_left, lsrb) % 1000000000;
         long long from_right = sum_of_segment_for_node(tree, left, right, right_son_no, lsrb+1, local_right) % 1000000000;
        return from_left + from_right;
    }
}

 long long max_of_segment(seg_tree& tree,  long long left,  long long right){
    static  long long total_nodes = std::get<0>(tree);
    static  long long leaves = std::get<1>(tree);
    if(total_nodes == 1){
        if(left <= 1 && right >= 1) return local_value(tree, 1);
    }
     long long lsrb = left_sons_right_bound(tree, 1);
     long long from_left = max_of_segment_for_node(tree, left, right, 2, 1, lsrb);
     long long from_right = max_of_segment_for_node(tree, left, right, 3, lsrb+1, leaves);
    return max(from_left, from_right);
}

void insert_n_with_weight_k(seg_tree& tree,  long long n,  long long k){
     long long i = node_no_of_ith_leaf(tree, n);
    while(i != 1){
        std::get<2>(tree)[i-1] += k;
        i = father(i);
    }
    std::get<2>(tree)[0] += k;
}

 long long sum_of_vector( long long* vector,  long long n){
     long long res = 0;
    for( long long i = 0; i < n; i++){
        res += vector[i];
    }
    return res;
}

// Changes the vector of i-inversions of a given permutation into
// the vector of (i+1)-inversions.
void increase_k_for_inversion(seg_tree tree,  long long* vector,  long long n,  long long* a){
     long long aux;
    for( long long i = 0; i < n; i++){
        aux = sum_of_segment(tree, a[i]+1, n);
        insert_n_with_weight_k(tree, a[i], vector[i]);
        vector[i] = aux;
    }
}

 long long k_inversions( long long k, seg_tree tree,  long long* vector,  long long n,  long long* a){
    for( long long i = 0; i < k-1; i++){
        if(i != 0) clear_tree(tree);
        increase_k_for_inversion(tree, vector, n, a);
    }
    return sum_of_vector(vector, n);
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    long long n;
    std::cin >> n;
    long long k;
    std::cin >> k;

    long long a[n];
    for( long long i = 0; i < n; i++) std::cin >> a[i];

    long long inv_vector[n]; // The inversion vector of the permutation
    for( long long i = 0; i < n; i++) inv_vector[i] = 1;
    seg_tree tree = init_tree_n_leaves(n);

    return 0;
}