#include <iostream>
#include <tuple>
#include <cassert>
#include <climits>

// For an n x m rectangle of unitary cucumber fields 1x1, each of which has an
// integer fertility coefficient, we need to answer a certain number of queries
// about subrectangles included in that rectangle, asking what the difference
// between the most and the least fertile field is in a given rectangle.

// In order to do this efficiently, I have used 2D segment trees (i. e. segment
// trees of segment trees).

// total size, number of leaves, min_values, max_values, right bound
using seg_tree = std::tuple< int,  int, int*, int*, int*>;

// total size, number of leaves, seg_trees of values, right bound
using super_seg_tree = std::tuple<int, int, seg_tree**, int*>;


// for convenience, nodes are numbered 1-n (!)
int right_son( int* values,  int node_no,  int total_nodes){
    int res = 2*node_no+1;
    if(res > total_nodes) return -1;
    if(values[res-1] == -1) return -1;
    return res;
}

int sup_right_son(seg_tree** trees,  int node_no,  int total_nodes){
    int res = 2*node_no+1;
    if(res > total_nodes) return -1;
    if(trees[res-1] == nullptr) return -1;
    return res;
}

int fill_sup_right_son(int node_no,  int total_nodes){
    int res = 2*node_no+1;
    if(res > total_nodes) return -1;
    return res;
}

int fill_sup_left_son( int node_no,  int total_nodes){
    int res = 2*node_no;
    if(res > total_nodes) return -1;
    return res;
}

int left_son( int* values,  int node_no,  int total_nodes){
    int res = 2*node_no;
    if(res > total_nodes) return -1;
    if(values[res-1] == -1) return -1;
    return res;
}

int sup_left_son( seg_tree** trees,  int node_no,  int total_nodes){
    int res = 2*node_no;
    if(res > total_nodes) return -1;
    if(trees[res-1] == nullptr) return -1;
    return res;
}

// assuming that a given node does have a father
int father( int node_no){
    return node_no/2;
}

int node_no_of_ith_leaf(seg_tree& tree,  int k){
    static int beginning = std::get<0>(tree) - std::get<1>(tree);
    return beginning + k;
}

int sup_node_no_of_ith_leaf(super_seg_tree tree,  int k){
    static int beginning = std::get<0>(tree) - std::get<1>(tree);
    return beginning + k;
}

int aux_nodes_for_n_leaves( int n){
    int res = 1;
    do{
        res *= 2;
    }
    while(res < n);

    return res-1+n;
}

int nodes_for_n_leaves( int n){
    static int res = aux_nodes_for_n_leaves(n);
    return res;
}

int sup_nodes_for_n_leaves( int n){
    static int res = aux_nodes_for_n_leaves(n);
    return res;
}

int local_min_value(const seg_tree &tree,  int node_no){
    return std::get<2>(tree)[node_no-1];
}

int local_max_value(const seg_tree &tree,  int node_no){
    return std::get<3>(tree)[node_no-1];
}

int min(int a, int b){
    return a < b ? a : b;
}

int max(int a, int b){
    return a > b ? a : b;
}

int right_bound_for_node(const seg_tree& tree,  int node_no){
    return std::get<4>(tree)[node_no-1];
}

int sup_right_bound_for_node(const super_seg_tree& tree,  int node_no){
    return std::get<3>(tree)[node_no-1];
}

int fill_right_bound_for_node(int* values,  int node_id,  int* nodes,  int total_nodes){
    if(values[node_id-1] == -1) return -1;
    if(nodes[node_id-1] != -1){
        int res = nodes[node_id-1];
        return res;
    }
    if(right_son(values, node_id, total_nodes) != -1) {
        int res = fill_right_bound_for_node(values, right_son(values, node_id, total_nodes), nodes, total_nodes);
        fill_right_bound_for_node(values, left_son(values, node_id, total_nodes)
                , nodes, total_nodes);
        nodes[node_id - 1] = res;
        return res;
    }
    else if(left_son(values, node_id, total_nodes) != -1){
        int res = fill_right_bound_for_node(values, left_son(values, node_id, total_nodes), nodes, total_nodes);
        nodes[node_id - 1] = res;
        return res;
    }
    else{
        assert(false);
        return 0;
    }
}

int sup_fill_right_bound_for_node(seg_tree** trees,  int node_id,  int* nodes,  int total_nodes){
    if(node_id > total_nodes){
        return -1;
    }
    if(nodes[node_id-1] != -1){
        int res = nodes[node_id-1];
        return res;
    }
    int rson = fill_sup_right_son(node_id, total_nodes);
    if(rson != -1) {
        int res = sup_fill_right_bound_for_node(trees, rson, nodes, total_nodes);
        if(res != -1) {
            sup_fill_right_bound_for_node(trees, fill_sup_left_son(node_id, total_nodes), nodes, total_nodes);
            nodes[node_id - 1] = res;
            return res;
        }
    }
    int lson = fill_sup_left_son(node_id, total_nodes);
    if(lson != -1){
        int res = sup_fill_right_bound_for_node(trees, fill_sup_left_son( node_id, total_nodes), nodes, total_nodes);
        nodes[node_id - 1] = res;
        return res;
    }
    else{
        nodes[node_id-1] = -1;
        return -1;
    }


}

void fill_right_bounds( int* values,  int* nodes,  int total_nodes,  int leaves){
    for( int i = 0; i < leaves; i++) nodes[i + total_nodes - leaves ] = i+1;
    fill_right_bound_for_node(values, 1, nodes, total_nodes);
}

void sup_fill_right_bounds(seg_tree** trees,  int* nodes,  int total_nodes,  int leaves){
    for( int i = 0; i < leaves; i++){
        nodes[i + total_nodes - leaves ] = i+1;
    }
    sup_fill_right_bound_for_node(trees, 1, nodes, total_nodes);
}

void fill_parentage(int* min_values, int* max_values,  int start){
    while(start != 0){
        if(max_values[start-1] == 0) return;
        min_values[start-1] = INT_MAX;
        max_values[start-1] = 0;
        start /= 2;
    }
    max_values[start] = 0;
}

seg_tree* init_tree_n_leaves( int n){
    int total_nodes = nodes_for_n_leaves(n);
    seg_tree* result = (seg_tree*) malloc(sizeof(seg_tree));
    std::get<0>(*result) = total_nodes;
    std::get<1>(*result) = n;

    int* min_values = ( int*) malloc(total_nodes*sizeof( int));
    int* max_values = ( int*) malloc(total_nodes*sizeof( int));
    for( int i = 0; i < total_nodes; i++){
        min_values[i] = INT_MAX;
        max_values[i] = -1;
    }
    for( int i = 0; i < n; i++){
        fill_parentage(min_values, max_values, i+total_nodes-n+1);
    }
    std::get<2>(*result) = min_values;
    std::get<3>(*result) = max_values;

    int* right_bounds = ( int*) malloc(total_nodes*sizeof( int));
    for( int i = 0; i < total_nodes; i++) right_bounds[i] = -1;
    fill_right_bounds(max_values, right_bounds, total_nodes, n);
    std::get<4>(*result) = right_bounds;

    return result;
}


int left_sons_right_bound(const seg_tree& tree,  int node_no){
    return right_bound_for_node(tree, left_son(std::get<2>(tree), node_no, std::get<0>(tree)));
}

int sup_left_sons_right_bound(const super_seg_tree& tree,  int node_no){
    return sup_right_bound_for_node(tree, sup_left_son(std::get<2>(tree), node_no, std::get<0>(tree)));
}

int max_of_segment_for_node(const seg_tree& tree,  int left,  int right,  int node_no,  int local_left,  int local_right){
    if(left > right) return 0;
    static int total_nodes = std::get<0>(tree);
    if(node_no == -1 || node_no > total_nodes) return 0;
    if(local_left > right || local_right < left) return 0;
    if(left <= local_left && right >= local_right){
        return local_max_value(tree, node_no);
    }
    else {
        int left_son_no = left_son(std::get<2>(tree), node_no, total_nodes);
        int lsrb = right_bound_for_node(tree, left_son_no);
        int right_son_no = right_son(std::get<2>(tree), node_no, total_nodes);
        int from_left = max_of_segment_for_node(tree, left, right, left_son_no, local_left, lsrb);
        int from_right = max_of_segment_for_node(tree, left, right, right_son_no, lsrb+1, local_right);
        return max(from_left, from_right);
    }
}

int min_of_segment_for_node(const seg_tree& tree,  int left,  int right,  int node_no,  int local_left,  int local_right){
    if(left > right) return 0;
    static int total_nodes = std::get<0>(tree);
    if(node_no == -1 || node_no > total_nodes) return INT_MAX;
    if(local_left > right || local_right < left) return INT_MAX;
    if(left <= local_left && right >= local_right){
        return local_min_value(tree, node_no);
    }
    else {
        int left_son_no = left_son(std::get<2>(tree), node_no, total_nodes);
        int lsrb = right_bound_for_node(tree, left_son_no);
        int right_son_no = right_son(std::get<2>(tree), node_no, total_nodes);
        int from_left = min_of_segment_for_node(tree, left, right, left_son_no, local_left, lsrb);
        int from_right = min_of_segment_for_node(tree, left, right, right_son_no, lsrb+1, local_right);
        return min(from_left, from_right);
    }
}


int min_of_segment(const seg_tree& tree,  int left,  int right){
    static int total_nodes = std::get<0>(tree);
    static int leaves = std::get<1>(tree);
    if(total_nodes == 1){
        if(left <= 1 && right >= leaves) return local_min_value(tree, 1);
    }
    int lsrb = left_sons_right_bound(tree, 1);
    int from_left = min_of_segment_for_node(tree, left, right, 2, 1, lsrb);
    int from_right = min_of_segment_for_node(tree, left, right, 3, lsrb+1, leaves);
    int res = min(from_right, from_left);
    return res;
}

int sup_min_of_segment_for_node(const super_seg_tree& tree, int super_left, int super_right,
                                      int sub_left,  int sub_right,
                                      int node_no, int local_left, int local_right){
    static int total_nodes = std::get<0>(tree);
    if(node_no == -1 || node_no > total_nodes) return INT_MAX;
    if(local_left > super_right || local_right < super_left) return INT_MAX;
    if(super_left <= local_left && super_right >= local_right){
        return min_of_segment(*(std::get<2>(tree)[node_no-1]), sub_left, sub_right);
    }
    else {
        int left_son_no = sup_left_son(std::get<2>(tree), node_no, total_nodes);
        int lsrb = sup_right_bound_for_node(tree, left_son_no);
        int right_son_no = sup_right_son(std::get<2>(tree), node_no, total_nodes);
        int from_left = sup_min_of_segment_for_node(tree, super_left, super_right, sub_left, sub_right, left_son_no, local_left, lsrb);
        int from_right = sup_min_of_segment_for_node(tree, super_left, super_right, sub_left,  sub_right, right_son_no, lsrb+1, local_right);
        return min(from_left, from_right);
    }
}

int sup_min_of_segment(super_seg_tree& tree, int super_left, int super_right,
                             int sub_left, int sub_right){
    static int total_nodes = std::get<0>(tree);
    static int leaves = std::get<1>(tree);
    if(total_nodes == 1){
        if(super_left <= 1 && super_right >= leaves) return min_of_segment(*(std::get<2>(tree)[1]), sub_left, sub_right);
    }
    int lsrb = sup_left_sons_right_bound(tree, 1);
    int from_left = sup_min_of_segment_for_node(tree, super_left, super_right, sub_left, sub_right, 2, 1, lsrb);
    int from_right = sup_min_of_segment_for_node(tree, super_left, super_right, sub_left, sub_right, 3, lsrb+1, leaves);
    int res = min(from_right, from_left);
    return res;
}

// for subtrees
int max_of_segment(const seg_tree& tree, int left, int right){
    static int total_nodes = std::get<0>(tree);
    static int leaves = std::get<1>(tree);
    if(total_nodes == 1){
        if(left <= 1 && right >= leaves) return local_min_value(tree, 1);
    }
    int lsrb = left_sons_right_bound(tree, 1);
    int from_left = max_of_segment_for_node(tree, left, right, 2, 1, lsrb);
    int from_right = max_of_segment_for_node(tree, left, right, 3, lsrb+1, leaves);
    int res = max(from_right, from_left);
    return res;
}
// for supertrees
int sup_max_of_segment_for_node(super_seg_tree& tree, int super_left, int super_right,
                                      int sub_left,  int sub_right,
                                      int node_no, int local_left, int local_right){
    static int total_nodes = std::get<0>(tree);
    if(node_no == -1 || node_no > total_nodes) return 0;
    if(local_left > super_right || local_right < super_left) return 0;
    if(super_left <= local_left && super_right >= local_right){
        return max_of_segment(*(std::get<2>(tree)[node_no-1]), sub_left, sub_right);
    }
    else {
        int left_son_no = sup_left_son(std::get<2>(tree), node_no, total_nodes);
        int lsrb = sup_right_bound_for_node(tree, left_son_no);
        int right_son_no = sup_right_son(std::get<2>(tree), node_no, total_nodes);
        int from_left = sup_max_of_segment_for_node(tree, super_left, super_right, sub_left, sub_right, left_son_no, local_left, lsrb);
        int from_right = sup_max_of_segment_for_node(tree, super_left, super_right, sub_left,  sub_right, right_son_no, lsrb+1, local_right);
        return max(from_left, from_right);
    }
}

int sup_max_of_segment(super_seg_tree& tree, int super_left, int super_right,
                             int sub_left, int sub_right){
    static int total_nodes = std::get<0>(tree);
    static int leaves = std::get<1>(tree);
    if(total_nodes == 1){
        if(super_left <= 1 && super_right >= leaves) return max_of_segment(*(std::get<2>(tree)[1]), sub_left, sub_right);
    }
    int lsrb = sup_left_sons_right_bound(tree, 1);
    int from_left = sup_max_of_segment_for_node(tree, super_left, super_right, sub_left, sub_right, 2, 1, lsrb);
    int from_right = sup_max_of_segment_for_node(tree, super_left, super_right, sub_left, sub_right, 3, lsrb+1, leaves);
    int res = max(from_right, from_left);
    return res;
}

void insert_n_with_weight_k(seg_tree& tree,  int n,  int k){
    int i = node_no_of_ith_leaf(tree, n);
    bool found_bigger = false;
    bool found_smaller = false;
    while(i != 1 && (!found_bigger || !found_smaller)){
        if(k < local_min_value(tree, i)){
            std::get<2>(tree)[i-1] = k;
        }
        else found_smaller = true;
        if(k > local_max_value(tree, i)){
            std::get<3>(tree)[i-1] = k;
        }
        else found_bigger = true;
        i = father(i);
    }
    if(k < local_min_value(tree, 1)) {
        std::get<2>(tree)[0] = k;
    }
    if(k > local_max_value(tree, 1)) {
        std::get<3>(tree)[0] = k;
    }

}

seg_tree* copy_seg_tree(seg_tree* tree){
    seg_tree* result = (seg_tree*) malloc(sizeof(seg_tree));
    std::get<0>(*result) = std::get<0>(*tree);
    std::get<1>(*result) = std::get<1>(*tree);
    int* min_values = (int*) malloc(std::get<0>(*tree)*sizeof(int));
    int* max_values = (int*) malloc(std::get<0>(*tree)*sizeof(int));
    int* right_bounds = (int*) malloc(std::get<0>(*tree)*sizeof(int));
    for(int i = 0; i < std::get<0>(*tree); i++){
        min_values[i] = std::get<2>(*tree)[i];
        max_values[i] = std::get<3>(*tree)[i];
        right_bounds[i] = std::get<4>(*tree)[i];
    }
    std::get<2>(*result) = min_values;
    std::get<3>(*result) = max_values;
    std::get<4>(*result) = right_bounds;
    return result;
}

seg_tree* seg_tree_merge(seg_tree* birch, seg_tree* poplar){
    if(poplar == nullptr)
        return birch == nullptr ? nullptr : copy_seg_tree(birch);
    assert(birch != nullptr);
    seg_tree* result = (seg_tree*) malloc(sizeof(seg_tree));
    std::get<0>(*result) = std::get<0>(*birch);
    std::get<1>(*result) = std::get<1>(*birch);
    int* min_values = (int*) malloc(std::get<0>(*birch)*sizeof(int));
    int* max_values = (int*) malloc(std::get<0>(*birch)*sizeof(int));
    int* right_bounds = (int*) malloc(std::get<0>(*birch)*sizeof(int));
    for(int i = 0; i < std::get<0>(*birch); i++){
        min_values[i] = min(std::get<2>(*birch)[i], std::get<2>(*poplar)[i]);
        max_values[i] = max(std::get<3>(*birch)[i], std::get<3>(*poplar)[i]);
        right_bounds[i] = std::get<4>(*birch)[i];
    }
    std::get<2>(*result) = min_values;
    std::get<3>(*result) = max_values;
    std::get<4>(*result) = right_bounds;
    return result;
}

super_seg_tree construct_super_tree(int total_leaves, seg_tree** sub_trees){
    super_seg_tree res;
    int total_nodes = sup_nodes_for_n_leaves(total_leaves);
    std::get<0>(res) = total_nodes;
    std::get<1>(res) = total_leaves;
    int* right_bound = (int*) malloc(total_nodes* sizeof(int));
    for( int i = 0; i < total_nodes; i++) right_bound[i] = -1;
    int first_leaf = sup_node_no_of_ith_leaf(res, 1);
    seg_tree** trees = (seg_tree**) malloc(total_nodes * sizeof(seg_tree*));
    for(int i = 0; i < total_nodes; i++) trees[i] = nullptr;
    for(int i = 0; i < total_leaves; i++) trees[first_leaf-1+i] = sub_trees[i];
    sup_fill_right_bounds(trees, right_bound, total_nodes, total_leaves);
    int i;
    int j = 1;
    while(j < total_nodes) j *= 2;
    j /= 4;
    while(j > 0){
        i = j;
        while(i <= 2*j-1) {
            trees[i - 1] = seg_tree_merge(2*i -1 >= total_nodes || trees[i * 2 - 1] == nullptr? nullptr :trees[i * 2 - 1],
                                          2*i >= total_nodes || trees[i * 2] == nullptr? nullptr : trees[i * 2]);
            i++;
        }
        j /= 2;
    }
    std::get<2>(res) = trees;
    std::get<3>(res) = right_bound;
    return res;
}

void process_query(super_seg_tree& sup_tree, int sub_left, int super_left, int sub_right, int super_right){
    int a = sup_max_of_segment(sup_tree, super_left, super_right, sub_left,sub_right );
    int b = sup_min_of_segment(sup_tree, super_left, super_right, sub_left,sub_right );
    std::cout << a-b << "\n";
}

void process_queries(super_seg_tree& sup_tree, int k) {
    int x1;
    int y1;
    int x2;
    int y2;
    for (int i = 0; i < k; i++) {
        std::cin >> x1;
        std::cin >> y1;
        std::cin >> x2;
        std::cin >> y2;
        process_query(sup_tree, y1 + 1, x1 + 1, y2 + 1, x2 + 1);
    }
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    int n;
    std::cin >> n;
    int m;
    std::cin >> m;
    int k;
    std::cin >> k;
    int a;
    seg_tree** basic_trees = (seg_tree**) malloc(n * sizeof(seg_tree*));
    seg_tree* new_tree;
    for(int i = 0; i < n; i++){
        new_tree = init_tree_n_leaves(m);
        for(int j = 1; j <= m; j++){
            std::cin >> a;
            insert_n_with_weight_k(*new_tree, j, a);
        }
        basic_trees[i] = new_tree;
    }
    super_seg_tree sup_tree = construct_super_tree(n, basic_trees);
    process_queries(sup_tree, k);
    return 0;
}